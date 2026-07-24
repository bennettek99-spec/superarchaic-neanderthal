"""
significance.py — is the Neanderthal signal at Denisovan CDDRs above chance?

Turns the descriptive z-tails into inference for the success-criterion question:
does the confound-controlled Neanderthal depth at Denisovan-specific CDDRs exceed
what a spatially-matched null produces (Model A), or not?

Two tests, both on the chr16-22 window tables (default 50 kb):

1. CIRCULAR-SHIFT PERMUTATION NULL (headline). Statistic S = mean, over Denisovan-
   specific CDDR windows, of the mean Neanderthal residual-excess depth (each
   Neanderthal's div_*_afr residualized on div_den_afr, z-scored — i.e. depth NOT
   explained by the shared per-window rate). Null: independently circular-shift the
   Neanderthal-excess track within each chromosome (preserving its autocorrelation)
   and recompute S over the fixed CDDR positions, N times. Empirical p = P(S_null >= S_obs).
   A large p => Neanderthals are NOT unusually deep at Denisovan CDDRs => no Model A.

2. LEAVE-ONE-CHROMOSOME-OUT JACKKNIFE on the four focal deep-divergence rates and on
   S, to show the conclusion is not driven by any single chromosome.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"
OUT = ROOT / "results" / "cddr"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vcflib import robust_z


def _resid_z(y, x):
    ok = np.isfinite(y) & np.isfinite(x)
    b1, b0 = np.polyfit(x[ok], y[ok], 1)
    out = np.full_like(y, np.nan, dtype=float)
    out[ok] = robust_z((y - (b0 + b1 * x))[ok])
    return out


def load(chroms, win):
    fr = []
    for c in chroms:
        p = TABLES / f"windows.chr{c}.win{win}.tsv"
        if p.exists():
            fr.append(pd.read_csv(p, sep="\t"))
    df = pd.concat(fr, ignore_index=True)
    df = df[(df["callable_frac"] >= 0.30) & np.isfinite(df["div_den_afr"])].copy()
    df["chrom"] = df["chrom"].astype(int)
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    nea = df[["div_alt_afr", "div_vin_afr", "div_chag_afr"]].mean(axis=1)
    df["den_excess"] = df["div_den_afr"] - nea
    df["z_den_excess"] = robust_z(df["den_excess"].to_numpy())
    # confound-controlled Neanderthal depth = residual of each nea_afr on den_afr, averaged
    exc = np.column_stack([_resid_z(df[f"div_{s}_afr"].to_numpy(), df["div_den_afr"].to_numpy())
                           for s in ("alt", "vin", "chag")])
    df["nea_exc"] = np.nanmean(exc, axis=1)
    df["artifact"] = ((df["mean_rm"].fillna(0) > 0.6) | (df["mean_map20"].fillna(1) < 0.7) |
                      (df["callable_frac"] < 0.40))
    return df


def _S(nea_exc, cddr_mask):
    v = nea_exc[cddr_mask]
    v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else np.nan


def permutation(df, z_cddr=3.0, n_perm=2000, seed=1, exclude_artifact=True):
    rng = np.random.default_rng(seed)
    cddr = (df["z_den_excess"].to_numpy() >= z_cddr)
    if exclude_artifact:
        cddr &= ~df["artifact"].to_numpy()
    nea = df["nea_exc"].to_numpy()
    chrom = df["chrom"].to_numpy()
    S_obs = _S(nea, cddr)
    # per-chromosome index blocks for independent circular shifts
    blocks = [np.where(chrom == c)[0] for c in np.unique(chrom)]
    null = np.empty(n_perm)
    for i in range(n_perm):
        shifted = nea.copy()
        for idx in blocks:
            shifted[idx] = np.roll(nea[idx], rng.integers(1, len(idx)))
        null[i] = _S(shifted, cddr)
    p = (np.sum(null >= S_obs) + 1) / (n_perm + 1)
    return S_obs, null, p, int(cddr.sum())


def jackknife(df, z_cddr=3.0):
    rows = []
    chroms = sorted(df["chrom"].unique())
    for drop in [None] + chroms:
        d = df if drop is None else df[df["chrom"] != drop]
        u = np.isfinite(d["div_den_afr"])
        rates = {g: float(np.mean(robust_z(d.loc[u, f"div_{g}_afr"].to_numpy()) >= z_cddr))
                 for g in ("den", "alt", "vin", "chag")}
        cddr = (d["z_den_excess"].to_numpy() >= z_cddr) & ~d["artifact"].to_numpy()
        rows.append(dict(drop=("none" if drop is None else f"chr{drop}"),
                         **{f"rate_{g}": rates[g] for g in rates},
                         S=_S(d["nea_exc"].to_numpy(), cddr)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=["16", "17", "18", "19", "20", "21", "22"])
    ap.add_argument("--win", type=int, default=50000)
    ap.add_argument("--z-cddr", type=float, default=3.0)
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args()
    df = load(args.chroms, args.win)
    print(f"windows: {len(df)} usable ({args.win} bp) across chr {','.join(args.chroms)}")

    S_obs, null, p, n_cddr = permutation(df, z_cddr=args.z_cddr, n_perm=args.n_perm)
    print("\n== Circular-shift permutation null (Model-A co-location test) ==")
    print(f"  Denisovan-specific CDDRs (non-artifact): {n_cddr}")
    print(f"  S_obs (mean Neanderthal residual-excess z at CDDRs): {S_obs:+.3f}")
    print(f"  null mean {null.mean():+.3f}  sd {null.std():.3f}  95%[{np.percentile(null,2.5):+.3f},{np.percentile(null,97.5):+.3f}]")
    print(f"  empirical p (S_obs >= null): {p:.3f}")
    verdict = ("NO excess -> no Model-A signal" if p > 0.05
               else "excess above chance -> candidate Model-A signal (examine further)")
    print(f"  => {verdict}")

    jk = jackknife(df, z_cddr=args.z_cddr)
    print("\n== Leave-one-chromosome-out jackknife ==")
    with pd.option_context("display.float_format", lambda v: f"{v:.4f}"):
        print(jk.to_string(index=False))
    full = jk.iloc[0]
    print(f"\n  focal-rate ordering (full set): den={full.rate_den:.4f} "
          f"alt={full.rate_alt:.4f} vin={full.rate_vin:.4f} chag={full.rate_chag:.4f}")
    print(f"  Denisovan is {'NOT ' if full.rate_den <= max(full.rate_alt, full.rate_vin, full.rate_chag) else ''}"
          f"the most deep-enriched archaic "
          f"({'no Denisovan-specific excess' if full.rate_den <= full.rate_alt else 'Denisovan-elevated'}).")
    OUT.mkdir(parents=True, exist_ok=True)
    jk.to_csv(OUT / f"jackknife.win{args.win}.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame({"null_S": null}).to_csv(OUT / f"permnull.win{args.win}.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT/f'jackknife.win{args.win}.tsv'} and permnull.win{args.win}.tsv")


if __name__ == "__main__":
    main()
