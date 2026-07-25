"""
significance.py — is the Neanderthal signal at Denisovan CDDRs above chance?

CORRECTED: the original version of this test was circular
--------------------------------------------------------
It selected candidate windows on

    z_den_excess >= 3,   den_excess = div_den_afr - mean(div_alt_afr, div_vin_afr, div_chag_afr)

and then reported, as the test statistic, the mean Neanderthal residual depth at those
same windows. Selecting windows for having HIGH Denisovan and LOW Neanderthal
divergence and then measuring how low Neanderthal divergence is there cannot come out
any other way: the reported S_obs = -3.35 against a null s.d. of 0.14 is roughly 25
sigma, which is a machine-precision statement about the selection rule, not about
Neanderthals. The old p = 1.000 was not evidence against Model A.

Two independent things are fixed here:

* SELECTION. Candidate windows are now flagged on Denisovan depth ALONE, normalized by
  the external human-vs-ancestral substitution rate (R_den = div_den_afr / sub_rate,
  src/ratemap.py). Nothing about the Neanderthals enters the selection.
* ORTHOGONALITY. The Neanderthal statistic is the residual of R_nea on R_den taken in
  QUANTILE BINS of R_den rather than by a single straight line. A linear residual is
  orthogonal to its regressor only on average; in a heavy upper tail a curved
  relationship leaves a systematic residual exactly where the candidates are. Binned
  residualization removes that by construction.

Note what this test can and cannot do even when correct: it is still a CO-LOCATION
test, and simulation puts its power under a true Model A near 6%, because superarchaic
ancestry sits at different loci in each lineage after ~420 ky. It is retained as a
Model-B probe and as a reference point. The tests with real Model-A power live in
src/modelA.py.

Two tests, both on the window tables (default 50 kb):

1. CIRCULAR-SHIFT PERMUTATION NULL. Statistic S = mean binned-residual Neanderthal
   depth over candidate windows. Null: independently circular-shift the Neanderthal
   track within each chromosome (preserving its autocorrelation) and recompute S over
   the fixed candidate positions. Empirical p = P(S_null >= S_obs).

2. LEAVE-ONE-CHROMOSOME-OUT JACKKNIFE on the four focal deep-divergence rates and on
   S, to show no single chromosome drives the conclusion.
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


def _binned_resid_z(y, x, nbins=None, per_bin=50):
    """Residual of y on x taken WITHIN quantile bins of x, then z-scored.

    A straight-line residual is orthogonal to x only on average. Any curvature in the
    y-on-x relationship leaves a systematic residual in the upper tail of x -- which is
    precisely where the candidate windows are selected, so a linear residual would
    reintroduce the selection artefact this module exists to remove. Subtracting a bin
    mean makes the residual mean zero inside every stratum of x.

    The bin COUNT matters more than it looks. With equal-count bins and a long-tailed
    x, the topmost bin spans an enormous x-range, and the candidate windows all sit at
    its upper edge -- so a coarse binning leaves most of the bias in place exactly
    where it does damage. Measured on curved synthetic data, the mean residual over the
    top 1% of x was +6.6 with 20 bins, +3.0 with 50, and 0.00 with 100+. Bins are
    therefore sized to ~`per_bin` points (default 50), not fixed at some round number.
    """
    y = np.asarray(y, float); x = np.asarray(x, float)
    out = np.full(len(y), np.nan)
    ok = np.isfinite(y) & np.isfinite(x)
    if nbins is None:
        nbins = int(np.clip(ok.sum() // per_bin, 20, 400))
    if ok.sum() < 10 * nbins:
        return _resid_z(y, x)
    edges = np.quantile(x[ok], np.linspace(0, 1, nbins + 1))
    edges[0] -= 1e-9; edges[-1] += 1e-9
    which = np.digitize(x, edges) - 1
    resid = np.full(len(y), np.nan)
    for b in range(nbins):
        m = ok & (which == b)
        if m.sum() >= 5:
            resid[m] = y[m] - np.mean(y[m])
    good = np.isfinite(resid)
    out[good] = robust_z(resid[good])
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

    # rate-normalized depth: the Neanderthal-INDEPENDENT selection variable
    if "sub_rate" in df.columns and np.isfinite(df["sub_rate"]).any():
        sr = df["sub_rate"].to_numpy(float)
        sr = np.where(np.isfinite(sr) & (sr > 0), sr, np.nan)
    else:
        sr = np.ones(len(df))
    for s in ("den", "alt", "vin", "chag"):
        df[f"R_{s}"] = df[f"div_{s}_afr"].to_numpy(float) / sr
    df["R_nea"] = df[["R_alt", "R_vin", "R_chag"]].mean(axis=1)
    df["z_R_den"] = robust_z(df["R_den"].to_numpy())

    # LEGACY selection variable, retained only to reproduce and expose the old bias
    nea = df[["div_alt_afr", "div_vin_afr", "div_chag_afr"]].mean(axis=1)
    df["den_excess"] = df["div_den_afr"] - nea
    df["z_den_excess"] = robust_z(df["den_excess"].to_numpy())

    # Neanderthal statistic, orthogonal to the selection variable by construction
    df["nea_exc"] = _binned_resid_z(df["R_nea"].to_numpy(), df["R_den"].to_numpy())
    # legacy Neanderthal statistic (linear residual on raw Denisovan divergence)
    exc = np.column_stack([_resid_z(df[f"div_{s}_afr"].to_numpy(), df["div_den_afr"].to_numpy())
                           for s in ("alt", "vin", "chag")])
    df["nea_exc_legacy"] = np.nanmean(exc, axis=1)

    df["artifact"] = ((df["mean_rm"].fillna(0) > 0.6) | (df["mean_map20"].fillna(1) < 0.7) |
                      (df["callable_frac"] < 0.40))
    return df


def _S(nea_exc, cddr_mask):
    v = nea_exc[cddr_mask]
    v = v[np.isfinite(v)]
    return float(v.mean()) if len(v) else np.nan


def permutation(df, z_cddr=3.0, n_perm=2000, seed=1, exclude_artifact=True,
                mode="unbiased"):
    """Circular-shift null for Neanderthal depth at Denisovan candidate windows.

    mode='unbiased' (default) selects candidates on rate-normalized DENISOVAN depth
    alone and measures a binned residual, so selection and statistic are orthogonal.
    mode='legacy' reproduces the original circular test -- candidates chosen on
    Denisovan-minus-Neanderthal excess, statistic measured on Neanderthal depth -- and
    is provided only to demonstrate the artefact.
    """
    rng = np.random.default_rng(seed)
    if mode == "legacy":
        cddr = (df["z_den_excess"].to_numpy() >= z_cddr)
        nea = df["nea_exc_legacy"].to_numpy()
    else:
        cddr = (df["z_R_den"].to_numpy() >= z_cddr)
        nea = df["nea_exc"].to_numpy()
    if exclude_artifact:
        cddr &= ~df["artifact"].to_numpy()
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


def wholebin_control(df, per_bin=50, ks=(1, 2, 3, 5)):
    """S computed over the top k WHOLE residual bins instead of above a z threshold.

    The binned residual has mean zero inside each bin by construction, so a candidate
    set made of complete bins carries no partial-bin bias. Any excess that survives here
    is not an artefact of where the threshold happened to fall. Measured on chr13-22,
    the z>=3 selection gives S = +0.415 while whole-bin selections of comparable size
    give +0.03 to +0.10 -- i.e. most of the apparent excess is the cut, not the data.
    """
    x = df["R_den"].to_numpy(float)
    y = df["R_nea"].to_numpy(float)
    r = _binned_resid_z(y, x, per_bin=per_bin)
    order = np.argsort(-x)
    n = len(x)
    nb = int(np.clip(n // per_bin, 20, 400))
    per = max(n // nb, 1)
    out = {}
    for k in ks:
        idx = order[:k * per]
        out[f"top {k} whole bin(s) ({len(idx)} win)"] = float(np.nanmean(r[idx]))
    return out


def jackknife(df, z_cddr=3.0):
    rows = []
    chroms = sorted(df["chrom"].unique())
    for drop in [None] + chroms:
        d = df if drop is None else df[df["chrom"] != drop]
        u = np.isfinite(d["div_den_afr"])
        rates = {g: float(np.mean(robust_z(d.loc[u, f"div_{g}_afr"].to_numpy()) >= z_cddr))
                 for g in ("den", "alt", "vin", "chag")}
        # select on rate-normalized Denisovan depth, matching permutation(mode='unbiased').
        # Selecting here on z_den_excess while reporting the corrected nea_exc statistic
        # would mix the biased selection with the fixed statistic.
        cddr = (d["z_R_den"].to_numpy() >= z_cddr) & ~d["artifact"].to_numpy()
        rows.append(dict(drop=("none" if drop is None else f"chr{drop}"),
                         **{f"rate_{g}": rates[g] for g in rates},
                         S=_S(d["nea_exc"].to_numpy(), cddr)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+",
                    default=[str(c) for c in range(13, 23)])
    ap.add_argument("--win", type=int, default=50000)
    ap.add_argument("--z-cddr", type=float, default=3.0)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--show-legacy", action="store_true",
                    help="also run the original biased test, to expose the artefact")
    args = ap.parse_args()
    df = load(args.chroms, args.win)
    print(f"windows: {len(df)} usable ({args.win} bp) across chr {','.join(args.chroms)}")

    S_obs, null, p, n_cddr = permutation(df, z_cddr=args.z_cddr, n_perm=args.n_perm,
                                         mode="unbiased")
    print("\n== Circular-shift permutation null, de-biased (co-location test) ==")
    print(f"  candidates selected on rate-normalized Denisovan depth alone: {n_cddr}")
    print(f"  S_obs (mean binned-residual Neanderthal depth at candidates): {S_obs:+.3f}")
    print(f"  null mean {null.mean():+.3f}  sd {null.std():.3f}  95%[{np.percentile(null,2.5):+.3f},{np.percentile(null,97.5):+.3f}]")
    print(f"  empirical p (S_obs >= null): {p:.3f}")

    wb = wholebin_control(df)
    print("\n  WHOLE-BIN CONTROL (the number that matters):")
    for k, v in wb.items():
        print(f"    {k:34s} S = {v:+.3f}")
    print("  A z-threshold cuts THROUGH a residual bin, so the candidates are the upper")
    print("  part of the bins they fall in and inherit a positive residual. Defining the")
    print("  same candidate count at whole-bin boundaries removes that, and S collapses")
    print("  by roughly an order of magnitude. The de-biased S_obs above is therefore")
    print("  still mostly selection, not signal -- much smaller than the legacy artefact,")
    print("  but not evidence.")
    print("\n  This test cannot be rescued by better conditioning. Condition fully on")
    print("  Denisovan depth and it is vacuous by construction; condition less and it is")
    print("  dominated by the ~0.73 Denisovan-Neanderthal correlation that shared ancestry")
    print("  produces under EVERY model, including the null. Its permutation null is the")
    print("  wrong null. Treat the output as descriptive only; the powered, externally")
    print("  calibrated tests are src/modelA.py (T1) with sims/power_calibration.py.")

    if args.show_legacy:
        Sl, nl, pl, ncl = permutation(df, z_cddr=args.z_cddr, n_perm=args.n_perm,
                                      mode="legacy")
        print("\n== LEGACY (biased) version, for comparison only ==")
        print(f"  candidates: {ncl}   S_obs {Sl:+.3f}   null {nl.mean():+.3f} +/- {nl.std():.3f}"
              f"   p={pl:.3f}")
        print(f"  S_obs sits {abs(Sl-nl.mean())/max(nl.std(),1e-9):.0f} null s.d. from the null "
              f"mean. That is the selection rule, not biology: candidates were chosen for")
        print("  having LOW Neanderthal divergence and the statistic then measures it.")

    jk = jackknife(df, z_cddr=args.z_cddr)
    print("\n== Leave-one-chromosome-out jackknife ==")
    with pd.option_context("display.float_format", lambda v: f"{v:.4f}"):
        print(jk.to_string(index=False))
    full = jk.iloc[0]
    print(f"\n  focal-rate ordering (full set): den={full.rate_den:.4f} "
          f"alt={full.rate_alt:.4f} vin={full.rate_vin:.4f} chag={full.rate_chag:.4f}")
    nea_max = max(full.rate_alt, full.rate_vin, full.rate_chag)
    if full.rate_den > nea_max:
        print("  Denisovan has the highest focal deep-window rate (Denisovan-elevated).")
    else:
        print("  Denisovan does NOT have the highest focal deep-window rate "
              "-> no Denisovan-specific excess.")
    print("  NB: equal rates across all four archaics is what BOTH the null and Model A")
    print("      predict (Model A lifts every lineage together), so this ordering")
    print("      discriminates Model B only -- it is not a test of Model A.")
    OUT.mkdir(parents=True, exist_ok=True)
    jk.to_csv(OUT / f"jackknife.win{args.win}.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame({"null_S": null}).to_csv(OUT / f"permnull.win{args.win}.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT/f'jackknife.win{args.win}.tsv'} and permnull.win{args.win}.tsv")


if __name__ == "__main__":
    main()
