"""
neanderthal_compare.py — classify CDDRs by Neanderthal sharing (Models A vs B).

Consumes the window tables from scan_windows.py and, for every Candidate Deep
Divergence Region (flagged on Denisovan-vs-African divergence), asks whether the
three Neanderthals (Altai, Vindija, Chagyrskaya) SHARE the deep-divergence signal.

Two per-Neanderthal signals per window, z-scored against the genome-wide (usable-
window) distribution:
  * nea_afr_z : Neanderthal-vs-African divergence (its own coalescence depth vs moderns)
  * den_nea_z : Denisovan-vs-that-Neanderthal divergence

Confound control. Denisovan-vs-African and Neanderthal-vs-African divergence are
baseline-correlated (shared ancestry + per-window mutation/recombination-rate
variation), so a raw "Neanderthal is also deep" test is biased. We therefore also
compute the residual of nea_afr after regressing it on den_afr genome-wide; a CDDR
only counts as Neanderthal-shared if the Neanderthal is deep BEYOND that prediction
(residual excess), which is the signature of genuinely shared deep ancestry rather
than a shared rate artifact.

Interpretation (never "superarchaic"; always candidate/consistent-with):
  Model A (Neandersovan): at CDDRs, Neanderthals show EXCESS depth vs modern AND are
    NOT unusually diverged from Denisovan there (they carry the same deep haplotype).
  Model B (Denisovan-only): at CDDRs, Neanderthals are NOT deep vs modern AND ARE
    unusually diverged from Denisovan (Denisovan alone carries the deep haplotype).

Categories per CDDR:
  1  shared by Denisovan + ALL 3 Neanderthals   (Model-A-consistent)
  2  shared by Denisovan + SOME Neanderthals
  3  unique to Denisovan                         (Model-B-consistent)
  4  inconclusive (Neanderthal callability/data insufficient)
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
sys.path.insert(0, str(ROOT.parent / "archaic-introgression"))
try:
    from archaic.windows import robust_z
except Exception:                                        # pragma: no cover
    def robust_z(x):
        x = np.asarray(x, float); m = np.nanmedian(x)
        s = 1.4826 * np.nanmedian(np.abs(x - m))
        return (x - m) / s if s else np.zeros_like(x)

NEAS = {"alt": "Altai", "vin": "Vindija", "chag": "Chagyrskaya"}


def _resid_excess(y, x):
    """Residual of y after OLS on x (both 1-D, NaN-aware), z-scored. Excess = high y for its x."""
    y = np.asarray(y, float); x = np.asarray(x, float)
    ok = np.isfinite(y) & np.isfinite(x)
    out = np.full(len(y), np.nan)
    if ok.sum() < 10:
        return out
    b1, b0 = np.polyfit(x[ok], y[ok], 1)
    resid = y - (b0 + b1 * x)
    out[ok] = robust_z(resid[ok])
    return out


def classify(win: int, chroms, z_hi=2.0, z_lo=1.0, min_callable_frac=0.30):
    frames = []
    for c in chroms:
        p = TABLES / f"windows.chr{c}.win{win}.tsv"
        if p.exists():
            frames.append(pd.read_csv(p, sep="\t"))
    if not frames:
        raise FileNotFoundError(f"no window tables for win{win} (run scan_windows.py)")
    df = pd.concat(frames, ignore_index=True)

    usable = (df["callable_frac_den"] >= min_callable_frac) & np.isfinite(df["div_den_afr"])
    base = df[usable].copy()

    # genome-wide z baselines (computed on usable windows, mapped back)
    for s in NEAS:
        df[f"z_{s}_afr"] = np.nan
        df.loc[usable, f"z_{s}_afr"] = robust_z(base[f"div_{s}_afr"].to_numpy())
        df[f"z_den_{s}"] = np.nan
        df.loc[usable, f"z_den_{s}"] = robust_z(base[f"div_den_{s}"].to_numpy())
        # confound-controlled excess of Neanderthal depth given Denisovan depth
        df[f"exc_{s}_afr"] = np.nan
        df.loc[usable, f"exc_{s}_afr"] = _resid_excess(
            base[f"div_{s}_afr"].to_numpy(), base["div_den_afr"].to_numpy())

    def nea_callable_frac(row, s):
        cb = row.get(f"callable_{s}", np.nan)
        return cb / win if np.isfinite(cb) else np.nan

    cddr = df[df["cddr"]].copy()
    cats, shared_counts, notes = [], [], []
    for _, r in cddr.iterrows():
        shared = 0
        inconclusive = 0
        per = []
        for s in NEAS:
            cf = nea_callable_frac(r, s)
            if not np.isfinite(cf) or cf < min_callable_frac or not np.isfinite(r[f"z_{s}_afr"]):
                inconclusive += 1
                per.append(f"{s}:NA")
                continue
            # shared-deep: excess depth vs modern AND not unusually diverged from Denisovan
            is_shared = (r[f"exc_{s}_afr"] >= z_hi) and (r[f"z_den_{s}"] <= z_lo)
            is_denonly = (r[f"z_den_{s}"] >= z_hi) and (r[f"exc_{s}_afr"] < z_lo)
            if is_shared:
                shared += 1
                per.append(f"{s}:shared")
            elif is_denonly:
                per.append(f"{s}:denonly")
            else:
                per.append(f"{s}:ambig")
        if inconclusive == len(NEAS):
            cat = 4
        elif shared == len(NEAS) - inconclusive and shared >= 1:
            cat = 1 if inconclusive == 0 else 2
        elif shared >= 1:
            cat = 2
        elif shared == 0 and inconclusive < len(NEAS):
            cat = 3
        else:
            cat = 4
        cats.append(cat)
        shared_counts.append(shared)
        notes.append(";".join(per))
    cddr["category"] = cats
    cddr["n_nea_shared"] = shared_counts
    cddr["nea_detail"] = notes

    OUT.mkdir(parents=True, exist_ok=True)
    cols = (["chrom", "start", "end", "z_den_afr", "div_den_afr", "category",
             "n_nea_shared", "nea_detail"]
            + [f"exc_{s}_afr" for s in NEAS] + [f"z_den_{s}" for s in NEAS]
            + ["mean_map20", "mean_rm", "mean_ur", "callable_frac_den"])
    cols = [c for c in cols if c in cddr.columns]
    cddr = cddr.sort_values("z_den_afr", ascending=False)
    outp = OUT / f"neanderthal_classified.win{win}.tsv"
    cddr[cols].to_csv(outp, sep="\t", index=False, float_format="%.6g")

    n = len(cddr)
    counts = {k: int((cddr["category"] == k).sum()) for k in (1, 2, 3, 4)}
    conc = n - counts[4]
    print(f"\nwin{win}: {n} CDDRs across chr {','.join(map(str,chroms))}")
    print(f"  Category 1 (Den + ALL Neanderthals, Model-A-consistent): {counts[1]}")
    print(f"  Category 2 (Den + SOME Neanderthals):                    {counts[2]}")
    print(f"  Category 3 (Denisovan-only, Model-B-consistent):         {counts[3]}")
    print(f"  Category 4 (inconclusive):                               {counts[4]}")
    if conc:
        print(f"  --> of {conc} conclusive CDDRs, "
              f"{100*(counts[1]+counts[2])/conc:.0f}% show ANY Neanderthal sharing, "
              f"{100*counts[1]/conc:.0f}% shared by all three.")
    print(f"  wrote {outp}")
    return cddr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wins", nargs="+", type=int, default=[20000, 50000, 100000])
    ap.add_argument("--chroms", nargs="+", default=["21", "22"])
    ap.add_argument("--z-hi", type=float, default=2.0)
    ap.add_argument("--z-lo", type=float, default=1.0)
    args = ap.parse_args()
    for w in args.wins:
        try:
            classify(w, args.chroms, z_hi=args.z_hi, z_lo=args.z_lo)
        except FileNotFoundError as e:
            print(e)


if __name__ == "__main__":
    main()
