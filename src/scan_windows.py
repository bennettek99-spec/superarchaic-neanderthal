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

POLARIZED SITE PATTERNS (added for the Model-A test battery)
------------------------------------------------------------
Divergence alone cannot separate Model A from the null, because superarchaic ancestry
entering the Neanderthal+Denisovan ancestor raises Denisovan AND Neanderthal depth
together -- any Den-vs-Nea contrast cancels it by construction. So each site is also
POLARIZED against the ancestral allele (CAnc from the Altai/Denisovan INFO, falling
back to the 1000G EPO AA field) and tallied into lineage-specific derived-allele
patterns, restricted to sites where all four archaics are callable and the derived
allele is absent from Africans:

  pat_nea_all   derived in ALL THREE Neanderthals, ancestral in Denisovan
                -> the Model-A signature: deep ancestry retained on the Neanderthal
                   side of the 420 ka split but drifted out of Denisovans
  pat_den_only  derived in Denisovan alone            -> the Model-B signature
  pat_all_arch  derived in all four                   -> ordinary Nea-Den shared ancestry
  pat_{alt,vin,chag}_only  single-Neanderthal private -> drift/error baseline

These are counts of DEEP-BRANCH mutations, and their per-window CLUSTERING (not their
mean, which is confounded with the split time) is what distinguishes introgressed
haplotype blocks from incomplete lineage sorting. Compositions such as
pat_nea_all / n_der_any are rate-free: local mutation rate scales numerator and
denominator alike.

Outputs:
  results/tables/windows.chr{C}.win{W}.tsv      (every window's metrics)
  results/cddr/cddr.chr{C}.win{W}.tsv           (flagged CDDRs)
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vcflib as V
from vcflib import robust_z, empirical_p

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
MASKS = ROOT / "data" / "masks"
TABLES = ROOT / "results" / "tables"
CDDRDIR = ROOT / "results" / "cddr"
RATEDIR = ROOT / "results" / "ratemap"

ARCHAICS = ["Denisova", "AltaiNea", "Vindija33.19", "Chagyrskaya"]
SHORT = {"Denisova": "den", "AltaiNea": "alt", "Vindija33.19": "vin", "Chagyrskaya": "chag"}
# mask files follow the EVA FilterBed subdir names (Altai, not AltaiNea)
MASK_STEM = {"Denisova": "Denisova", "AltaiNea": "Altai",
             "Vindija33.19": "Vindija33.19", "Chagyrskaya": "Chagyrskaya"}
MISS = int(V.MISS)
_code = V.code_bases                                 # vectorized base -> int8 code


def _num(series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float32)


def load_chromosome(chrom: str) -> dict:
    """Align the four archaic caches + the 1000G frequency cache onto one position axis.

    Returns a dict of parallel numpy arrays rather than a merged DataFrame. Four
    successive pandas outer-merges on object-dtype columns were the memory and time
    bottleneck on the larger chromosomes; a single sorted union plus searchsorted
    scatter is O(n log n) with int64 keys and int8 payloads.
    """
    raw, pos_sets = {}, []
    for g in ARCHAICS:
        p = CACHE / f"{g}.chr{chrom}.variants.tsv.gz"
        if not p.exists():
            raise FileNotFoundError(f"missing cache {p} (run extract_variants.py)")
        cols = ["pos", "ref", "a1", "a2"]
        if g in ("Denisova", "AltaiNea"):             # only these carry EPO annotation
            cols += ["canc", "map20", "rm", "ur", "cpg", "bsc"]
        d = pd.read_csv(p, sep="\t", usecols=cols, dtype=str)
        d["pos"] = d["pos"].astype(np.int64)
        raw[g] = d
        pos_sets.append(d["pos"].to_numpy())
    mod = pd.read_csv(CACHE / f"1000G.chr{chrom}.freq.tsv.gz", sep="\t", dtype=str)
    mod["pos"] = mod["pos"].astype(np.int64)
    pos_sets.append(mod["pos"].to_numpy())

    pos = np.unique(np.concatenate(pos_sets))
    n = len(pos)
    out = {"pos": pos, "pos0": pos - 1}

    def scatter(dst, src_pos, values):
        dst[np.searchsorted(pos, src_pos)] = values
        return dst

    refc = np.full(n, MISS, dtype=np.int8)
    ancc = np.full(n, MISS, dtype=np.int8)
    for g in ARCHAICS:
        d = raw[g]
        gp = d["pos"].to_numpy()
        for tag, col in (("a1", "a1"), ("a2", "a2")):
            arr = np.full(n, MISS, dtype=np.int8)
            out[f"{SHORT[g]}_{tag}"] = scatter(arr, gp, _code(d[col].fillna(".")))
        # hg19 REF and the chimp-human ancestral base are position-level facts: take
        # them from whichever genome supplies them, first valid wins.
        rc = _code(d["ref"].fillna("."))
        need = (refc[np.searchsorted(pos, gp)] == MISS) & (rc != MISS)
        idx = np.searchsorted(pos, gp)
        refc[idx[need]] = rc[need]
        if "canc" in d.columns:
            ac = _code(d["canc"].fillna("."))
            need = (ancc[idx] == MISS) & (ac != MISS)
            ancc[idx[need]] = ac[need]
        if g == "Denisova":
            for name in ("map20", "rm", "ur", "cpg", "bsc"):
                arr = np.full(n, np.nan, dtype=np.float32)
                out[name] = scatter(arr, gp, _num(d[name]))

    mp = mod["pos"].to_numpy()
    midx = np.searchsorted(pos, mp)
    mrefc = _code(mod["ref"].fillna("."))
    need = (refc[midx] == MISS) & (mrefc != MISS)
    refc[midx[need]] = mrefc[need]
    aac = _code(mod["aa"].fillna(".")) if "aa" in mod.columns else np.full(len(mp), MISS, np.int8)
    need = (ancc[midx] == MISS) & (aac != MISS)
    ancc[midx[need]] = aac[need]

    malt = np.full(n, MISS, dtype=np.int8)
    scatter(malt, mp, _code(mod["alt"].fillna(".")))
    in_1000g = np.zeros(n, dtype=bool)
    in_1000g[midx] = True
    p_afr = np.zeros(n, dtype=np.float32)
    p_eur = np.zeros(n, dtype=np.float32)
    scatter(p_afr, mp, np.nan_to_num(_num(mod["afr"])))
    scatter(p_eur, mp, np.nan_to_num(_num(mod["eur"])))

    out.update(refc=refc, ancc=ancc, malt=malt, p_afr=p_afr, p_eur=p_eur,
               in_1000g=in_1000g)
    for name in ("map20", "rm", "ur", "cpg", "bsc"):
        out.setdefault(name, np.full(n, np.nan, dtype=np.float32))
    return out


def _fill_genome(u, s, mask, refc, pos0):
    """Return (a0,a1) code arrays for genome `s`: cached alleles, else hom-ref if
    callable, else MISS. Non-SNP (indel/N) alleles collapse to MISS."""
    a1 = u[f"{s}_a1"].copy()
    a2 = u[f"{s}_a2"].copy()
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


def polarized_patterns(u, genc, in_common):
    """Per-site polarized derived-allele patterns across the four archaics.

    Returns a dict of boolean per-site arrays. A site qualifies only when it is in the
    common callable mask, has a confident ancestral base, all four archaics are called,
    and the African derived frequency is resolvable. `afr0` means the derived allele is
    entirely absent from the African panel, which is what makes a pattern lineage-
    specific rather than shared standing variation.
    """
    ancc, refc, malt = u["ancc"], u["refc"], u["malt"]
    anc_ok = ancc != MISS

    # African derived frequency. Sites absent from the 1000G biallelic-SNP cache are
    # monomorphic-reference in moderns (the convention extract_modern.py documents),
    # so their derived frequency is 0 or 1 depending on whether REF is the ancestral.
    q = np.full(len(ancc), np.nan, dtype=np.float32)
    mono = ~u["in_1000g"] & anc_ok & (refc != MISS)
    q[mono] = np.where(refc[mono] == ancc[mono], 0.0, 1.0)
    poly = u["in_1000g"] & anc_ok & (refc != MISS) & (malt != MISS)
    alt_is_anc = poly & (malt == ancc)
    ref_is_anc = poly & (refc == ancc)
    q[alt_is_anc] = 1.0 - u["p_afr"][alt_is_anc]
    q[ref_is_anc] = u["p_afr"][ref_is_anc]
    # triallelic w.r.t. the ancestral base (neither REF nor ALT is ancestral): unusable
    q[poly & ~alt_is_anc & ~ref_is_anc] = np.nan

    der, called = {}, {}
    for s in ("den", "alt", "vin", "chag"):
        a1, a2 = genc[s]
        called[s] = (a1 != MISS) & (a2 != MISS)
        der[s] = called[s] & anc_ok & ((a1 != ancc) | (a2 != ancc))

    base = in_common & anc_ok & np.isfinite(q)
    for s in called:
        base = base & called[s]
    afr0 = base & (q == 0)

    D, A, Vv, C = der["den"], der["alt"], der["vin"], der["chag"]
    NEA_ALL = A & Vv & C
    return {
        "n_pol": base,                                   # polarizable, all-4-callable
        "n_der_any": base & (D | A | Vv | C | (q > 0)),  # rate-proportional denominator
        "pat_den_only": afr0 & D & ~A & ~Vv & ~C,        # Model-B signature
        "pat_nea_all": afr0 & ~D & NEA_ALL,              # Model-A signature
        "pat_nea_any": afr0 & ~D & (A | Vv | C),
        "pat_all_arch": afr0 & D & NEA_ALL,              # ordinary Nea-Den shared ancestry
        "pat_alt_only": afr0 & ~D & A & ~Vv & ~C,        # single-lineage drift/error
        "pat_vin_only": afr0 & ~D & ~A & Vv & ~C,
        "pat_chag_only": afr0 & ~D & ~A & ~Vv & C,
    }


def _win_sum(pos0, values, win, nbins):
    """Sum `values` (NaN-aware) into disjoint windows of width `win`; also count."""
    wid = (pos0 // win).astype(np.int64)
    finite = np.isfinite(values)
    s = np.bincount(wid[finite], weights=values[finite], minlength=nbins)
    n = np.bincount(wid[finite], minlength=nbins)
    return s[:nbins], n[:nbins]


def _win_count(pos0, flag, win, nbins):
    """Count True entries of a boolean per-site array into fixed-width windows."""
    wid = (pos0[flag] // win).astype(np.int64)
    return np.bincount(wid, minlength=nbins)[:nbins].astype(np.int64)


def scan(chrom: str, win_sizes=(20000, 50000, 100000), min_callable_frac=0.30):
    u = load_chromosome(chrom)
    pos0 = u["pos0"]
    refc = u["refc"]

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
    malt, p_afr, p_eur = u["malt"], u["p_afr"], u["p_eur"]

    # polarized derived-allele patterns (Model-A battery; see the module docstring)
    patterns = polarized_patterns(u, genc, in_common)

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
    map20, rm, ur, bsc = u["map20"], u["rm"], u["ur"], u["bsc"]

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
        # polarized site-pattern counts per window
        for key, flag in patterns.items():
            rows[key] = _win_count(pos0, flag, win, nbins)
        # single common-mask denominator for every divergence (comparable across genomes)
        common_bp = V.callable_bp_windows(common, win, nbins)
        rows["callable_common"] = common_bp
        # per-genome callable bp kept for reporting / missingness diagnostics
        for g in ARCHAICS:
            rows[f"callable_{SHORT[g]}"] = V.callable_bp_windows(masks[g], win, nbins)

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

        # ---- external mutation-rate denominator (src/ratemap.py), if built ---------
        # div_x_afr / sub_rate removes local mutation-rate variation WITHOUT
        # contrasting the archaics against each other, so unlike den_excess it leaves
        # Model A's shared Neanderthal+Denisovan signal intact. Absent a rate map the
        # columns are NaN and downstream falls back to the contrast statistics.
        rp = RATEDIR / f"rate.chr{chrom}.win{win}.tsv"
        if rp.exists():
            rate = pd.read_csv(rp, sep="\t", usecols=["start", "sub_rate", "repeat_frac",
                                                      "anc_valid_bp", "hum_sub"])
            df = df.merge(rate, on="start", how="left")
        else:
            for c in ("sub_rate", "repeat_frac", "anc_valid_bp", "hum_sub"):
                df[c] = np.nan
        sr = df["sub_rate"].to_numpy(dtype=float)
        sr = np.where(np.isfinite(sr) & (sr > 0), sr, np.nan)
        for s in ("den", "alt", "vin", "chag"):
            df[f"R_{s}_afr"] = df[f"div_{s}_afr"] / sr      # rate-normalized depth
        df["R_nea_afr"] = df[["R_alt_afr", "R_vin_afr", "R_chag_afr"]].mean(axis=1)
        # rate-free polarized compositions (local mutation rate cancels in the ratio)
        nder = df["n_der_any"].replace(0, np.nan)
        for k in ("pat_nea_all", "pat_den_only", "pat_all_arch",
                  "pat_alt_only", "pat_vin_only", "pat_chag_only"):
            df[f"f_{k[4:]}"] = df[k] / nder

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


def _scan_one(job):
    """Worker entry point: chromosomes are fully independent, so they parallelize."""
    chrom, wins, mcf = job
    try:
        scan(chrom, win_sizes=wins, min_callable_frac=mcf)
        return (chrom, None)
    except Exception as e:
        return (chrom, f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=["21", "22"])
    ap.add_argument("--wins", nargs="+", type=int, default=[20000, 50000, 100000])
    ap.add_argument("--min-callable-frac", type=float, default=0.30)
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel chromosomes; each worker peaks near 2-3 GB, so on a "
                         "16 GB laptop 4 is a sensible ceiling")
    args = ap.parse_args()
    jobs = [(c, tuple(args.wins), args.min_callable_frac) for c in args.chroms]
    if args.jobs <= 1 or len(jobs) == 1:
        for j in jobs:
            _scan_one(j)
        return
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for chrom, err in ex.map(_scan_one, jobs):
            if err:
                print(f"chr{chrom}: FAILED — {err}", file=sys.stderr)


if __name__ == "__main__":
    main()
