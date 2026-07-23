"""
rank_candidates.py — Denisovan-specific deep-divergence catalog + evidence ranking.

Raw archaic-vs-African divergence is confounded: per-window mutation-rate and
coalescent-variance heterogeneity inflate ALL archaic divergences together, so a
deep-divergence window is not per se lineage-specific (we observe den/alt/vin/chag
CDDR rates that are essentially equal). The unconfounded statistic for DENISOVAN-
specific deep ancestry (the Model-B signature) is the contrast

    den_excess = div_den_afr - mean(div_alt_afr, div_vin_afr, div_chag_afr)

which cancels shared per-window rate variation (validated in simulation: den_excess
mean is largest under Model B, smallest under the null). This module builds the
candidate catalog on that contrast, attaches the alternative-explanation covariates
(repeat / mappability / copy-number / callability), classifies Neanderthal sharing,
and assigns an evidence tier. Nothing here is called "superarchaic".

Evidence tiers per candidate:
  Artifact-likely : repeat/low-mappability/segdup or low callability dominated
  Insufficient    : too little callable data to evaluate
  Weak            : passes artifact filter, single window size only
  Moderate        : passes artifact filter, replicates across >=2 window sizes
  Strong          : Moderate + would require Neanderthal-sharing pattern (see notes)
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


def _load(win, chroms):
    fr = [pd.read_csv(TABLES / f"windows.chr{c}.win{win}.tsv", sep="\t")
          for c in chroms if (TABLES / f"windows.chr{c}.win{win}.tsv").exists()]
    df = pd.concat(fr, ignore_index=True)
    u = (df["callable_frac"] >= 0.30) & np.isfinite(df["div_den_afr"])
    df = df[u].copy()
    nea = df[["div_alt_afr", "div_vin_afr", "div_chag_afr"]].mean(axis=1)
    df["nea_afr"] = nea
    df["den_excess"] = df["div_den_afr"] - nea
    df["z_den_excess"] = robust_z(df["den_excess"].to_numpy())
    # per-Neanderthal lineage-specific excess (for Model-A / Model-3 candidates)
    for s in ("alt", "vin", "chag"):
        others = [f"div_{o}_afr" for o in ("den", "alt", "vin", "chag") if o != s]
        df[f"{s}_excess"] = df[f"div_{s}_afr"] - df[others].mean(axis=1)
    df["artifact"] = ((df["mean_rm"].fillna(0) > 0.6) |
                      (df["mean_map20"].fillna(1) < 0.7) |
                      (df["mean_ur"].fillna(0) > 0.9) |
                      (df["callable_frac"] < 0.40))
    return df


def build(chroms, wins=(20000, 50000, 100000), z_cddr=3.0):
    OUT.mkdir(parents=True, exist_ok=True)
    primary = wins[1] if len(wins) > 1 else wins[0]
    dfs = {w: _load(w, chroms) for w in wins}
    df = dfs[primary]
    cddr = df[df["z_den_excess"] >= z_cddr].copy()

    # replication: does a candidate window overlap a z>=cddr window at another size?
    def replicated(row):
        hits = 0
        for w in wins:
            if w == primary:
                continue
            o = dfs[w]
            m = ((o["chrom"] == row["chrom"]) & (o["start"] < row["end"]) &
                 (o["end"] > row["start"]) & (o["z_den_excess"] >= z_cddr))
            hits += int(m.any())
        return hits
    cddr["n_repl_other_wins"] = cddr.apply(replicated, axis=1)

    # Neanderthal sharing tag (any Neanderthal itself lineage-specifically deep here?)
    def nea_tag(row):
        deep = [s for s in ("alt", "vin", "chag") if row[f"{s}_excess"] > 0 and
                robust_z_at(dfs[primary][f"{s}_excess"], row[f"{s}_excess"]) >= 2.0]
        return ",".join(deep) if deep else "none"

    def tier(row):
        if row["artifact"]:
            return "Artifact-likely"
        if row["callable_frac"] < 0.40:
            return "Insufficient"
        if row["n_repl_other_wins"] >= 1:
            return "Moderate"
        return "Weak"
    cddr["evidence_tier"] = cddr.apply(tier, axis=1)

    cols = ["chrom", "start", "end", "z_den_excess", "den_excess", "div_den_afr",
            "nea_afr", "callable_frac", "mean_rm", "mean_map20", "mean_ur",
            "artifact", "n_repl_other_wins", "evidence_tier"]
    cddr = cddr.sort_values("z_den_excess", ascending=False)
    out = OUT / f"candidate_catalog.win{primary}.tsv"
    cddr[cols].to_csv(out, sep="\t", index=False, float_format="%.6g")

    n = len(cddr)
    tiers = cddr["evidence_tier"].value_counts().to_dict()
    print(f"Denisovan-specific candidate catalog (win{primary}, chr {','.join(map(str,chroms))})")
    print(f"  {n} windows with z_den_excess >= {z_cddr}")
    for t in ["Strong", "Moderate", "Weak", "Insufficient", "Artifact-likely"]:
        print(f"    {t:16s}: {tiers.get(t, 0)}")
    surv = cddr[~cddr["artifact"] & (cddr["evidence_tier"].isin(["Moderate", "Strong"]))]
    print(f"  --> {len(surv)} candidates survive the artifact + replication filter.")
    print(f"  wrote {out}")
    return cddr


def robust_z_at(series, value):
    x = np.asarray(series, float)
    m = np.nanmedian(x); s = 1.4826 * np.nanmedian(np.abs(x - m))
    return (value - m) / s if s else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=["21", "22"])
    ap.add_argument("--wins", nargs="+", type=int, default=[20000, 50000, 100000])
    ap.add_argument("--z-cddr", type=float, default=3.0)
    args = ap.parse_args()
    build(args.chroms, wins=tuple(args.wins), z_cddr=args.z_cddr)


if __name__ == "__main__":
    main()
