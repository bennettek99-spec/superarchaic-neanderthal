"""
scan_windows.py — windowed divergence scan for Candidate Deep Divergence Regions.

Builds, per chromosome and per window size (20/50/100 kb), a table of pairwise
sequence divergence between the Denisovan, the three Neanderthals, and modern
African/European panels, plus callability and quality/annotation aggregates. The
primary coalescence-depth proxy is Denisovan-vs-African divergence (per the DEEP
2026 result that superarchaic ancestry elevates local archaic-modern divergence);
windows in its extreme upper tail — surviving callability and quality filters — are
flagged as CDDRs. Classifying whether Neanderthals share them (Model A) vs not
(Model B) happens downstream (neanderthal_compare.py).

Method (efficient all-sites handling):
  * DENOMINATORS (pairwise-callable bp per window) come from the BED masks.
  * NUMERATORS (differences) are summed over the union of variant sites only;
    invariant hom-reference sites contribute zero and are filled from the masks.
All divergence algebra is vectorized over sites (bases encoded as 0..3).

Outputs:
  results/tables/windows.chr{C}.win{W}.tsv      (every window's metrics)
  results/cddr/cddr.chr{C}.win{W}.tsv           (flagged CDDRs)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vcflib as V

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
MASKS = ROOT / "data" / "masks"
TABLES = ROOT / "results" / "tables"
CDDRDIR = ROOT / "results" / "cddr"

# window-scan reuses the main pipeline's robust z / empirical p if importable
try:
    sys.path.insert(0, str(ROOT.parent / "archaic-introgression"))
    from archaic.windows import robust_z, empirical_p
except Exception:                                    # pragma: no cover - fallback
    def robust_z(x):
        x = np.asarray(x, float); med = np.nanmedian(x)
        mad = np.nanmedian(np.abs(x - med)); s = 1.4826 * mad
        if not np.isfinite(s) or s == 0:
            s = np.nanstd(x)
        return (x - med) / s if s else np.zeros_like(x)

    def empirical_p(x):
        x = np.asarray(x, float); n = np.isfinite(x).sum()
        if n == 0:
            return np.full_like(x, np.nan)
        r = np.argsort(np.argsort(x))
        return np.clip(2.0 * np.minimum(r + 1, n - r) / n, 1.0 / n, 1.0)

ARCHAICS = ["Denisova", "AltaiNea", "Vindija33.19", "Chagyrskaya"]
SHORT = {"Denisova": "den", "AltaiNea": "alt", "Vindija33.19": "vin", "Chagyrskaya": "chag"}
# mask files follow the EVA FilterBed subdir names (Altai, not AltaiNea)
MASK_STEM = {"Denisova": "Denisova", "AltaiNea": "Altai",
             "Vindija33.19": "Vindija33.19", "Chagyrskaya": "Chagyrskaya"}
_B2C = {"A": 0, "C": 1, "G": 2, "T": 3}
MISS = -2


def _code(series) -> np.ndarray:
    """Map a pandas Series of base strings to int8 codes (A/C/G/T=0..3, else MISS)."""
    return series.map(lambda b: _B2C.get(b, MISS)).astype("int64").to_numpy()


def load_chromosome(chrom: str):
    """Outer-merge archaic variant caches + modern freqs into one union table."""
    dfs = {}
    for g in ARCHAICS:
        p = CACHE / f"{g}.chr{chrom}.variants.tsv.gz"
        if not p.exists():
            raise FileNotFoundError(f"missing cache {p} (run extract_variants.py)")
        s = SHORT[g]
        d = pd.read_csv(p, sep="\t", dtype=str)
        keep = {"pos": "pos", "ref": "ref", "a1": f"{s}_a1", "a2": f"{s}_a2"}
        if g == "Denisova":                          # keep the annotation columns once
            keep.update({"map20": "map20", "rm": "rm", "ur": "ur",
                         "cpg": "cpg", "bsc": "bsc", "canc": "canc"})
        d = d[list(keep)].rename(columns=keep)
        d["pos"] = d["pos"].astype(np.int64)
        dfs[g] = d
    mod = pd.read_csv(CACHE / f"1000G.chr{chrom}.freq.tsv.gz", sep="\t", dtype=str)
    mod["pos"] = mod["pos"].astype(np.int64)
    mod = mod.rename(columns={"ref": "mref", "alt": "malt"})

    u = dfs["Denisova"]
    for g in ARCHAICS[1:]:
        u = u.merge(dfs[g], on="pos", how="outer", suffixes=("", "_dup"))
        if "ref_dup" in u.columns:                   # coalesce the shared hg19 ref
            u["ref"] = u["ref"].fillna(u["ref_dup"]); u = u.drop(columns=["ref_dup"])
    u = u.merge(mod, on="pos", how="outer")
    u["ref"] = u["ref"].fillna(u["mref"])
    u = u.sort_values("pos").reset_index(drop=True)
    return u


def _fill_genome(u, s, mask, refc, pos0):
    """Return (a0,a1) code arrays for genome `s`: cached alleles, else hom-ref if
    callable, else MISS. Non-SNP (indel/N) alleles collapse to MISS."""
    a1 = _code(u[f"{s}_a1"].fillna("."))
    a2 = _code(u[f"{s}_a2"].fillna("."))
    has_call = (a1 != MISS) & (a2 != MISS)
    callable_here = V.mask_membership(mask, pos0)
    homref = callable_here & ~has_call & (refc != MISS)
    a1 = np.where(homref, refc, a1)
    a2 = np.where(homref, refc, a2)
    # anything still without two valid SNP alleles -> missing for divergence
    bad = (a1 == MISS) | (a2 == MISS)
    a1 = np.where(bad, MISS, a1); a2 = np.where(bad, MISS, a2)
    return a1, a2


def _dxy_arch(x0, x1, y0, y1):
    """Per-site dxy between two diploid archaic calls; NaN where either missing."""
    miss = (x0 == MISS) | (y0 == MISS)
    mm = ((x0 != y0).astype(np.float64) + (x0 != y1) + (x1 != y0) + (x1 != y1)) * 0.25
    mm[miss] = np.nan
    return mm


def _dxy_mod(x0, x1, refc, altc, p):
    """Per-site dxy between a diploid archaic call and a modern panel (freq p of altc)."""
    miss = (x0 == MISS)

    def one(a):
        out = np.where(a == altc, 1.0 - p, np.where(a == refc, p, 1.0))
        return out
    mm = 0.5 * (one(x0) + one(x1))
    mm[miss] = np.nan
    return mm


def _win_sum(pos0, values, win, nbins):
    """Sum `values` (NaN-aware) into disjoint windows of width `win`; also count."""
    wid = (pos0 // win).astype(np.int64)
    finite = np.isfinite(values)
    s = np.bincount(wid[finite], weights=values[finite], minlength=nbins)
    n = np.bincount(wid[finite], minlength=nbins)
    return s[:nbins], n[:nbins]


def scan(chrom: str, win_sizes=(20000, 50000, 100000), min_callable_frac=0.30):
    u = load_chromosome(chrom)
    pos0 = (u["pos"].to_numpy() - 1).astype(np.int64)
    refc = _code(u["ref"].fillna("."))

    masks = {g: V.read_mask(MASKS / f"{MASK_STEM[g]}.chr{chrom}.mask.bed.gz", chrom) for g in ARCHAICS}
    # COMMON mask = sites callable in ALL four archaics. Cross-genome divergences are
    # computed only here so that den/alt/vin/chag comparisons use identical sites and
    # denominators (the 2014 Altai/Denisovan and 2016-18 Vindija/Chagyrskaya callability
    # masks differ; per-genome denominators would make divergences non-comparable).
    common = masks[ARCHAICS[0]]
    for g in ARCHAICS[1:]:
        common = V.intervals_intersect(common, masks[g])
    in_common = V.mask_membership(common, pos0)
    genc = {}
    for g in ARCHAICS:
        s = SHORT[g]
        genc[s] = _fill_genome(u, s, masks[g], refc, pos0)

    # modern panels
    malt = _code(u["malt"].fillna("."))
    p_afr = pd.to_numeric(u["afr"], errors="coerce").fillna(0.0).to_numpy()
    p_eur = pd.to_numeric(u["eur"], errors="coerce").fillna(0.0).to_numpy()

    # per-site divergences
    site = {}
    pairs_arch = [("den", "alt"), ("den", "vin"), ("den", "chag"),
                  ("alt", "vin"), ("alt", "chag"), ("vin", "chag")]
    for a, b in pairs_arch:
        site[f"div_{a}_{b}"] = _dxy_arch(*genc[a], *genc[b])
    for s in ("den", "alt", "vin", "chag"):
        site[f"div_{s}_afr"] = _dxy_mod(genc[s][0], genc[s][1], refc, malt, p_afr)
        site[f"div_{s}_eur"] = _dxy_mod(genc[s][0], genc[s][1], refc, malt, p_eur)
    # restrict every divergence numerator to the common callable mask
    for k in list(site):
        site[k] = np.where(in_common, site[k], np.nan)

    # annotation (from Denisovan-carried columns; representative within-window)
    map20 = pd.to_numeric(u["map20"], errors="coerce").to_numpy()
    rm = pd.to_numeric(u["rm"], errors="coerce").to_numpy()
    ur = pd.to_numeric(u["ur"], errors="coerce").to_numpy()
    bsc = pd.to_numeric(u["bsc"], errors="coerce").to_numpy()

    TABLES.mkdir(parents=True, exist_ok=True)
    CDDRDIR.mkdir(parents=True, exist_ok=True)
    chrom_len = int(u["pos"].max())
    outs = {}
    for win in win_sizes:
        nbins = chrom_len // win + 1
        rows = {"chrom": [str(chrom)] * nbins,
                "start": (np.arange(nbins) * win).astype(np.int64),
                "end": (np.arange(nbins) * win + win).astype(np.int64)}
        # numerators + counts per window
        for key, vals in site.items():
            s, n = _win_sum(pos0, vals, win, nbins)
            rows[f"_num_{key}"] = s
            rows[f"_n_{key}"] = n
        # single common-mask denominator for every divergence (comparable across genomes)
        common_bp = np.array(
            [V.callable_bp_in_window(common, i * win, i * win + win) for i in range(nbins)],
            dtype=np.int64)
        rows["callable_common"] = common_bp
        # per-genome callable bp kept for reporting / missingness diagnostics
        for g in ARCHAICS:
            rows[f"callable_{SHORT[g]}"] = np.array(
                [V.callable_bp_in_window(masks[g], i * win, i * win + win) for i in range(nbins)],
                dtype=np.int64)

        df = pd.DataFrame(rows)
        denom = np.where(common_bp > 0, common_bp, np.nan)
        for a, b in pairs_arch:
            df[f"div_{a}_{b}"] = df[f"_num_div_{a}_{b}"] / denom
        for s in ("den", "alt", "vin", "chag"):
            for pop in ("afr", "eur"):
                df[f"div_{s}_{pop}"] = df[f"_num_div_{s}_{pop}"] / denom
        # annotation aggregates
        for name, arr in (("map20", map20), ("rm", rm), ("ur", ur), ("bsc", bsc)):
            sN, nN = _win_sum(pos0, arr, win, nbins)
            df[f"mean_{name}"] = sN / np.where(nN > 0, nN, np.nan)
        df["n_sites"] = df[[c for c in df.columns if c.startswith("_n_div_den_afr")]].sum(axis=1)
        df["callable_frac"] = common_bp / win
        df["callable_frac_den"] = df["callable_common"] / win  # back-compat alias

        df = df[[c for c in df.columns if not c.startswith("_")]]
        # CDDR flagging on Denisovan-vs-African divergence
        usable = (df["callable_frac_den"] >= min_callable_frac) & np.isfinite(df["div_den_afr"])
        df["z_den_afr"] = np.nan
        df.loc[usable, "z_den_afr"] = robust_z(df.loc[usable, "div_den_afr"].to_numpy())
        df["p_den_afr"] = np.nan
        df.loc[usable, "p_den_afr"] = empirical_p(df.loc[usable, "div_den_afr"].to_numpy())
        # Mappability/repeat control is already applied at the SITE level by the
        # mq25/mapab100 callability masks (the common mask). The window-mean Map20/RM
        # (computed over the repeat-enriched variant sites) is NOT used as a hard gate
        # -- it is kept as a reported covariate so the evidence ranking can label a
        # candidate "artifact-likely". A CDDR is a usable window in the deep tail; each
        # candidate is then examined against alternative explanations downstream.
        df["cddr"] = usable & (df["z_den_afr"] >= 3.0)
        df["low_quality"] = (df["mean_map20"].fillna(1) < 0.5) | (df["mean_rm"].fillna(0) > 0.8)

        tp = TABLES / f"windows.chr{chrom}.win{win}.tsv"
        df.to_csv(tp, sep="\t", index=False, float_format="%.6g")
        cd = df[df["cddr"]].sort_values("z_den_afr", ascending=False)
        cd.to_csv(CDDRDIR / f"cddr.chr{chrom}.win{win}.tsv", sep="\t", index=False, float_format="%.6g")
        print(f"chr{chrom} win{win}: {int(usable.sum())} usable windows, "
              f"{int(df['cddr'].sum())} CDDRs  (median div_den_afr="
              f"{np.nanmedian(df.loc[usable,'div_den_afr']):.5f})")
        outs[win] = df
    return outs


def _g(short):
    for k, v in SHORT.items():
        if v == short:
            return k
    raise KeyError(short)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=["21", "22"])
    ap.add_argument("--wins", nargs="+", type=int, default=[20000, 50000, 100000])
    ap.add_argument("--min-callable-frac", type=float, default=0.30)
    args = ap.parse_args()
    for c in args.chroms:
        scan(c, win_sizes=tuple(args.wins), min_callable_frac=args.min_callable_frac)


if __name__ == "__main__":
    main()
