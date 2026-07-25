"""
modelA.py — locus-free tests with real power for Model A (Neandersovan).

WHY THE PILOT DESIGN COULD NOT ANSWER THE QUESTION
---------------------------------------------------
The pilot asked "do Neanderthals show excess depth AT Denisovan-flagged windows?".
Simulation put that test's power at ~6% under a TRUE Model A, and the reason is
structural, not statistical bad luck: superarchaic ancestry entering the common
Neanderthal+Denisovan ancestor is retained at DIFFERENT loci in the two lineages
after ~420 ky of drift and recombination. Co-location is simply the wrong
expectation, so adding chromosomes to that test buys almost nothing.

Worse, the headline statistic

    den_excess = div_den_afr - mean(div_alt_afr, div_vin_afr, div_chag_afr)

cancels Model A exactly -- under Model A both terms rise together. And selecting
candidate windows on den_excess and then measuring mean Neanderthal residual depth at
those windows selects "high Denisovan, LOW Neanderthal" and then reports that
Neanderthals are low: a guaranteed negative result independent of any biology.

WHAT REPLACES IT
----------------
Three tests that do not require the two lineages to share loci, ordered by how much
they discriminate. All are computed on the same window-table schema, so they run
unchanged on real chromosomes and on simulated ones (sims/superarchaic_sim.py),
which is what turns them into calibrated power statements.

T1  CLADE CLUSTERING  (the headline; discriminates Model A from the null)
    Polarized pattern `pat_nea_all` = derived in ALL THREE Neanderthals, ancestral in
    Denisovan, absent in Africans. Under Model A those are deep-branch mutations sitting
    on introgressed HAPLOTYPE BLOCKS, so they are spatially clustered and overdispersed
    across windows. Under incomplete lineage sorting they are scattered. The raw amount
    of clustering is not specific -- rate variation, local Ne and mapping artifacts all
    cluster -- so the statistic is a CONTRAST of clustering against patterns that share
    those confounds but carry no Model-A signal:
        pat_alt_only  single-Neanderthal private   -> pure drift/error control
        pat_den_only  Denisovan private            -> Model-B signal
        pat_all_arch  all four archaics            -> ordinary Nea-Den shared ancestry
    A Model-A genome has pat_nea_all clustered like pat_den_only and unlike
    pat_alt_only. Confounds move all four together and cancel in the contrast.

T2  TAIL SYMMETRY  (discriminates Model B from everything else)
    Rate-normalized depth R_x = div_x_afr / sub_rate, with sub_rate the external
    human-vs-ancestral substitution density (src/ratemap.py) -- a denominator that
    removes mutation-rate heterogeneity WITHOUT contrasting the archaics against each
    other, so it leaves Model A intact. Tail mass of R_x at a common threshold is
    compared across the four archaics with a block jackknife. Model B lifts Denisovan
    alone; Model A lifts all four equally and is therefore NOT separable from the null
    by this test -- which is exactly why T1 exists and why the pilot's "all four rates
    are equal" observation was never evidence against Model A.

T3  CROSS-LINEAGE COVARIANCE  (discriminates Model A from Model 3, separate pulses)
    The aggregate, unthresholded version of the pilot's co-location test: the weighted
    covariance of the Denisovan and Neanderthal deep-residual tracks over ALL windows.
    Thresholding at z>=3 discarded the ~94% of shared regions the classifier missed;
    a covariance keeps every window's contribution. Its null is NOT permutation --
    Neanderthal and Denisovan depth are correlated (~0.7) by shared ancestry under
    every model -- so it is interpreted against the simulated M0/M1/M2/M3 baselines.

Uncertainty everywhere is a Busing-Meijer-van der Leeden weighted block jackknife over
contiguous ~5 Mb blocks, which respects linkage; nominal per-window independence would
badly understate the error.

Usage:
  python src/modelA.py --chroms 13 14 ... 22 --win 50000
  python src/modelA.py --tables results/sims/windows_M1_neandersovan.tsv   # sim input
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"
OUT = ROOT / "results" / "modelA"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vcflib import robust_z  # noqa: E402

# Regions excluded from the primary scan and reported separately. The MHC carries
# trans-species balancing selection: genuinely ancient haplotypes that are NOT archaic
# introgression, and the strongest deep-divergence signal in the human genome. Leaving
# it in would let one locus drive every statistic here.
SPECIAL_REGIONS = {"MHC": ("6", 28_000_000, 34_000_000)}

PATTERNS = ["pat_nea_all", "pat_den_only", "pat_all_arch",
            "pat_alt_only", "pat_vin_only", "pat_chag_only"]
LINEAGES = ["den", "alt", "vin", "chag"]


# --------------------------------------------------------------------- loading

def load(chroms=None, win=50000, tables=None, min_callable=0.30, min_der=30,
         exclude_special=True):
    """Load window tables (real or simulated) and attach the derived statistics."""
    if tables:
        frames = [pd.read_csv(p, sep="\t") for p in tables]
    else:
        frames = []
        for c in chroms:
            p = TABLES / f"windows.chr{c}.win{win}.tsv"
            if p.exists():
                frames.append(pd.read_csv(p, sep="\t"))
            else:
                print(f"  (no table for chr{c} win{win})", file=sys.stderr)
    if not frames:
        raise FileNotFoundError("no window tables found")
    df = pd.concat(frames, ignore_index=True)
    df["chrom"] = df["chrom"].astype(str)

    excluded = pd.DataFrame()
    if exclude_special:
        hit = np.zeros(len(df), bool)
        for name, (c, lo, hi) in SPECIAL_REGIONS.items():
            m = (df["chrom"] == str(c)) & (df["end"] > lo) & (df["start"] < hi)
            if m.any():
                df.loc[m, "special_region"] = name
            hit |= m.to_numpy()
        excluded = df[hit].copy()
        df = df[~hit].copy()

    keep = (df["callable_frac"] >= min_callable) & np.isfinite(df["div_den_afr"])
    if "n_der_any" in df.columns:
        keep &= df["n_der_any"] >= min_der
    df = df[keep].sort_values(["chrom", "start"]).reset_index(drop=True)

    # rate-normalized depth; falls back to raw divergence when no rate map exists
    if "sub_rate" in df.columns and np.isfinite(df["sub_rate"]).any():
        sr = df["sub_rate"].to_numpy(float)
        sr = np.where(np.isfinite(sr) & (sr > 0), sr, np.nan)
        df["_rate_used"] = "sub_rate"
    else:
        sr = np.ones(len(df))
        df["_rate_used"] = "none"
    for x in LINEAGES:
        df[f"R_{x}"] = df[f"div_{x}_afr"].to_numpy(float) / sr
    df["R_nea"] = df[["R_alt", "R_vin", "R_chag"]].mean(axis=1)
    return df, excluded


# ------------------------------------------------------- block jackknife (SEs)

def make_blocks(df, block_bp=5_000_000):
    """Contiguous ~5 Mb blocks within chromosomes; linkage lives inside a block."""
    key = df["chrom"].astype(str) + ":" + (df["start"] // block_bp).astype(str)
    return pd.factorize(key)[0]


def block_jackknife(stat_fn, df, blocks, weights=None):
    """Weighted delete-one-block jackknife (Busing, Meijer & van der Leeden 1999).

    stat_fn(sub_df) -> float. Returns (theta_hat, theta_jack, se). Weighting by each
    block's share of the data is what makes unequal block sizes safe.
    """
    theta_hat = stat_fn(df)
    ublocks = np.unique(blocks)
    g = len(ublocks)
    if g < 3 or not np.isfinite(theta_hat):
        return theta_hat, np.nan, np.nan
    w = np.ones(len(df)) if weights is None else np.asarray(weights, float)
    n = w.sum()
    thetas, hs = [], []
    for b in ublocks:
        sel = blocks != b
        m_j = w[~sel].sum()
        if m_j <= 0 or m_j >= n:
            continue
        t = stat_fn(df[sel])
        if not np.isfinite(t):
            continue
        thetas.append(t)
        hs.append(n / m_j)
    if len(thetas) < 3:
        return theta_hat, np.nan, np.nan
    thetas = np.array(thetas)
    hs = np.array(hs)
    theta_J = g * theta_hat - np.sum((1 - 1 / hs) * thetas)
    pseudo = hs * theta_hat - (hs - 1) * thetas
    var = np.mean((pseudo - theta_J) ** 2 / (hs - 1))
    return theta_hat, theta_J, float(np.sqrt(max(var, 0)))


# ------------------------------------------------------- clustering statistics

def dispersion(counts, offset):
    """Quasi-Poisson dispersion phi of `counts` given an exposure `offset`.

    phi = mean over windows of (obs - exp)^2 / exp, with exp = offset * (sum obs /
    sum offset). phi == 1 is Poisson scatter (what independent ILS mutations give);
    phi >> 1 means events pile into some windows, as a haplotype block of deep
    ancestry produces.

    CAUTION: phi is NOT comparable between patterns with different abundances. Under
    the usual clustered model Var = mu + alpha*mu^2, so phi = 1 + alpha*mu -- a pattern
    with twice the count shows twice the excess phi at identical relative clustering.
    Use relative_dispersion() for cross-pattern contrasts.
    """
    c = np.asarray(counts, float)
    o = np.asarray(offset, float)
    ok = np.isfinite(c) & np.isfinite(o) & (o > 0)
    c, o = c[ok], o[ok]
    if len(c) < 10 or c.sum() <= 0:
        return np.nan
    exp = o * (c.sum() / o.sum())
    return float(np.mean((c - exp) ** 2 / exp))


def relative_dispersion(counts, offset):
    """Mean-scale-free clustering coefficient alpha, from Var = mu + alpha * mu^2.

    Moment estimator alpha = sum[(obs-exp)^2 - exp] / sum[exp^2]. This is the
    negative-binomial relative dispersion: alpha == 0 is pure Poisson scatter and
    alpha is invariant to how abundant the pattern is, so alpha(pat_nea_all) and
    alpha(pat_alt_only) are directly comparable even though the two differ severalfold
    in count. This is the statistic the Model-A contrast is built on.
    """
    c = np.asarray(counts, float)
    o = np.asarray(offset, float)
    ok = np.isfinite(c) & np.isfinite(o) & (o > 0)
    c, o = c[ok], o[ok]
    if len(c) < 10 or c.sum() <= 0:
        return np.nan
    exp = o * (c.sum() / o.sum())
    num = np.sum((c - exp) ** 2 - exp)
    den = np.sum(exp ** 2)
    return float(num / den) if den > 0 else np.nan


def latent_autocorr(df, count_col, offset_col="n_der_any", lag=1, win=None):
    """Lag-`lag` autocorrelation of the LATENT intensity behind a count series.

    Write x_i = lambda_i + eps_i with eps_i independent Poisson sampling noise. Then
        Cov(x_i, x_i+lag) = Cov(lambda_i, lambda_i+lag)      (noise cancels at lag > 0)
        Var(x) - E[x]     = Var(lambda)
    so Cov(x_i, x_i+lag) / (Var(x) - E[x]) estimates the autocorrelation of lambda
    itself. Raw autocorrelation of counts is attenuated by Poisson noise in proportion
    to how RARE the pattern is, which would make an abundant pattern look more
    clustered than a rare one at identical biology -- exactly the artefact that makes
    the alpha values non-comparable across patterns of different abundance. This
    statistic is free of that attenuation and is the robust member of the T1 battery.
    """
    c = df[count_col].to_numpy(float)
    o = df[offset_col].to_numpy(float)
    ok = np.isfinite(c) & np.isfinite(o) & (o > 0)
    if ok.sum() < 50:
        return np.nan
    # normalize away the exposure so only clustering of the RATE remains
    x = c[ok] / (o[ok] / np.mean(o[ok]))
    sub = df[ok].copy()
    sub["_x"] = x
    num = _lag_cov(sub, "_x", lag=lag, win=win)
    den = np.var(x) - np.mean(c[ok])
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0:
        return np.nan
    return float(num / den)


def _lag_cov(df, col, lag=1, win=None):
    """Lag-`lag` autocovariance over genuinely adjacent windows within chromosomes."""
    x = df[col].to_numpy(float)
    chrom = df["chrom"].to_numpy()
    start = df["start"].to_numpy()
    if win is None:
        d = np.diff(np.unique(start))
        win = int(np.min(d)) if len(d) else 1
    a, b = [], []
    for c in np.unique(chrom):
        m = chrom == c
        s, v = start[m], x[m]
        order = np.argsort(s)
        s, v = s[order], v[order]
        adj = (s[lag:] - s[:-lag]) == win * lag
        va, vb = v[:-lag][adj], v[lag:][adj]
        good = np.isfinite(va) & np.isfinite(vb)
        a.append(va[good]); b.append(vb[good])
    a, b = np.concatenate(a), np.concatenate(b)
    if len(a) < 20:
        return np.nan
    mu = np.mean(np.concatenate([a, b]))
    return float(np.mean((a - mu) * (b - mu)))


def spatial_autocorr(df, col, lag=1, win=None):
    """Lag-`lag` spatial autocorrelation of a per-window composition within chromosomes.

    Only genuinely adjacent windows contribute (gaps from unusable windows are not
    silently bridged). Introgressed blocks span several windows and autocorrelate;
    ILS does not.
    """
    x = df[col].to_numpy(float)
    chrom = df["chrom"].to_numpy()
    start = df["start"].to_numpy()
    if win is None:
        d = np.diff(np.unique(start))
        win = int(np.min(d)) if len(d) else 1
    a, b = [], []
    for c in np.unique(chrom):
        m = chrom == c
        s, v = start[m], x[m]
        order = np.argsort(s)
        s, v = s[order], v[order]
        adj = (s[lag:] - s[:-lag]) == win * lag
        va, vb = v[:-lag][adj], v[lag:][adj]
        good = np.isfinite(va) & np.isfinite(vb)
        a.append(va[good])
        b.append(vb[good])
    a, b = np.concatenate(a), np.concatenate(b)
    if len(a) < 20 or a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def T1_clade_clustering(df, blocks, lag=1, win=50000):
    """Clustering of each polarized pattern, with block-jackknife SEs.

    The Model-A read is a CONTRAST: pat_nea_all clustered like pat_den_only and unlike
    the single-Neanderthal drift control pat_alt_only. Shared confounds (rate, Ne,
    mappability) inflate every row together and so cancel between rows.
    """
    rows = []
    for pat in PATTERNS:
        if pat not in df.columns:
            continue
        phi, _, phi_se = block_jackknife(
            lambda d, p=pat: dispersion(d[p], d["n_der_any"]), df, blocks)
        alpha, _, alpha_se = block_jackknife(
            lambda d, p=pat: relative_dispersion(d[p], d["n_der_any"]), df, blocks)
        comp = f"f_{pat[4:]}"
        if comp not in df.columns:
            df[comp] = df[pat] / df["n_der_any"].replace(0, np.nan)
        ac, _, ac_se = block_jackknife(
            lambda d, c=comp: spatial_autocorr(d, c, lag=lag, win=win), df, blocks)
        lac, _, lac_se = block_jackknife(
            lambda d, p=pat: latent_autocorr(d, p, lag=lag, win=win), df, blocks)
        rows.append(dict(pattern=pat, total=int(df[pat].sum()),
                         per_window=float(df[pat].mean()),
                         alpha=alpha, alpha_se=alpha_se,
                         latent_ac=lac, latent_ac_se=lac_se,
                         dispersion=phi, dispersion_se=phi_se,
                         autocorr=ac, autocorr_se=ac_se))
    return pd.DataFrame(rows)


def contrast_alpha(df, blocks, pat_a, pat_b):
    """alpha(pat_a) - alpha(pat_b), jackknifed as ONE statistic.

    Differencing two separately-jackknifed estimates would ignore their covariance
    (both respond to the same blocks); recomputing the difference inside each delete-
    one-block replicate propagates it correctly and gives a much tighter, honest SE.
    """
    def stat(d):
        a = relative_dispersion(d[pat_a], d["n_der_any"])
        b = relative_dispersion(d[pat_b], d["n_der_any"])
        return a - b
    est, _, se = block_jackknife(stat, df, blocks)
    return dict(a=pat_a, b=pat_b, diff=est, se=se,
                z=est / se if se and np.isfinite(se) and se > 0 else np.nan)


# ------------------------------------------------------------- tail statistics

def tail_mass(df, col, thr):
    v = df[col].to_numpy(float)
    v = v[np.isfinite(v)]
    return float(np.mean(v >= thr)) if len(v) else np.nan


def T2_tail_symmetry(df, blocks, q=0.99):
    """Deep-tail mass per lineage on the rate-normalized depth R_x.

    The threshold is a quantile of the POOLED four-lineage distribution, so the
    comparison is symmetric by construction. Model B lifts Denisovan alone. Model A
    lifts all four together and is NOT distinguishable from the null here -- this test
    exists to rule Model B in or out, not to test Model A.
    """
    pooled = np.concatenate([df[f"R_{x}"].to_numpy(float) for x in LINEAGES])
    pooled = pooled[np.isfinite(pooled)]
    thr = float(np.quantile(pooled, q))
    rows = []
    for x in LINEAGES:
        est, _, se = block_jackknife(
            lambda d, c=f"R_{x}": tail_mass(d, c, thr), df, blocks)
        rows.append(dict(lineage=x, tail_mass=est, se=se, threshold=thr))
    out = pd.DataFrame(rows)
    nea = out[out.lineage != "den"]
    d = out[out.lineage == "den"].iloc[0]
    # Denisovan minus the Neanderthal mean, with a jackknife SE on the contrast itself
    def contrast(sub):
        p = np.concatenate([sub[f"R_{x}"].to_numpy(float) for x in LINEAGES])
        t = float(np.quantile(p[np.isfinite(p)], q))
        return tail_mass(sub, "R_den", t) - np.mean(
            [tail_mass(sub, f"R_{x}", t) for x in ("alt", "vin", "chag")])
    c_est, _, c_se = block_jackknife(contrast, df, blocks)
    return out, dict(contrast=c_est, se=c_se,
                     z=c_est / c_se if c_se and np.isfinite(c_se) and c_se > 0 else np.nan,
                     den=d.tail_mass, nea_mean=float(nea.tail_mass.mean()))


def T3_cross_lineage_cov(df, blocks):
    """Aggregate (unthresholded) covariance of the Denisovan and Neanderthal depth tracks.

    Reported with a block-jackknife SE. Interpreted against the simulated M0/M1/M2/M3
    baselines, not against a permutation null: the two tracks are strongly correlated
    under EVERY model because the lineages share ancestry, so zero is not the null.
    """
    def stat(d):
        a = d["R_den"].to_numpy(float)
        b = d["R_nea"].to_numpy(float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 20:
            return np.nan
        return float(np.corrcoef(a[ok], b[ok])[0, 1])
    est, _, se = block_jackknife(stat, df, blocks)
    return dict(corr_den_nea=est, se=se)


def legacy_colocation(df, z_cddr=3.0):
    """The pilot's thresholded co-location test, kept only to quantify what it misses.

    Reported alongside the new statistics so the ~6%-power gap is visible rather than
    implicit. Candidate windows are flagged on Denisovan depth ALONE (not on a
    Denisovan-minus-Neanderthal contrast), which removes the selection bias that made
    the pilot's S_obs mechanically negative.
    """
    z_den = robust_z(df["R_den"].to_numpy(float))
    cddr = z_den >= z_cddr
    z_nea = robust_z(df["R_nea"].to_numpy(float))
    if cddr.sum() == 0:
        return dict(n_cddr=0, frac_shared=np.nan, mean_z_nea_at_cddr=np.nan)
    return dict(n_cddr=int(cddr.sum()),
                frac_shared=float(np.mean(z_nea[cddr] >= 2.0)),
                mean_z_nea_at_cddr=float(np.nanmean(z_nea[cddr])))


# --------------------------------------------------------------------- driver

def run(df, excluded=None, win=50000, label="observed", q=0.99, block_bp=5_000_000):
    blocks = make_blocks(df, block_bp)
    n_blocks = len(np.unique(blocks))
    print(f"\n{'='*74}\n{label}: {len(df):,} usable windows, {n_blocks} jackknife blocks "
          f"({block_bp/1e6:.0f} Mb), rate denominator = {df['_rate_used'].iloc[0]}\n{'='*74}")

    t1 = T1_clade_clustering(df, blocks, win=win)
    print("\n-- T1  clade clustering (Model A vs null) " + "-" * 32)
    print("   pattern         total  per-win   alpha (rel.disp)     latent_ac(lag1)   phi")
    for _, r in t1.iterrows():
        print(f"   {r['pattern']:14s} {r['total']:6d}  {r['per_window']:7.2f}  "
              f"{r['alpha']:7.3f} +/- {r['alpha_se']:<6.3f} "
              f"{r['latent_ac']:+7.3f} +/- {r['latent_ac_se']:<6.3f} {r['dispersion']:6.2f}")
    got = {r["pattern"]: r for _, r in t1.iterrows()}
    contrasts = []
    for ctl in ("pat_alt_only", "pat_all_arch", "pat_den_only"):
        if "pat_nea_all" in got and ctl in got:
            contrasts.append(contrast_alpha(df, blocks, "pat_nea_all", ctl))
    if contrasts:
        print("\n   CONTRASTS on alpha (jackknifed as a single statistic):")
        for c in contrasts:
            print(f"     alpha(pat_nea_all) - alpha({c['b']:<13s}) = {c['diff']:+.3f} "
                  f"+/- {c['se']:.3f}   z = {c['z']:+.2f}")
        print("   vs pat_alt_only : drift/error control -> Model A predicts > 0")
        print("   vs pat_all_arch : same depth/genome count on the SHARED archaic branch;")
        print("                     note Model A also feeds this pattern, so it is a")
        print("                     conservative (signal-contaminated) control")
        print("   vs pat_den_only : Model B's own signature -> Model A predicts ~0,")
        print("                     Model B predicts < 0")

    t2, t2c = T2_tail_symmetry(df, blocks, q=q)
    print(f"\n-- T2  tail symmetry at the pooled q={q} threshold " + "-" * 22)
    for _, r in t2.iterrows():
        print(f"   R_{r['lineage']:<5s} tail mass {r['tail_mass']:.4f} +/- {r['se']:.4f}")
    print(f"   Denisovan - Neanderthal-mean = {t2c['contrast']:+.4f} "
          f"+/- {t2c['se']:.4f}  (z = {t2c['z']:+.2f})")
    print("   z >> 0 favours Model B; z ~ 0 is equally consistent with Model A and the null.")

    t3 = T3_cross_lineage_cov(df, blocks)
    print(f"\n-- T3  cross-lineage covariance " + "-" * 42)
    print(f"   corr(R_den, R_nea) = {t3['corr_den_nea']:+.4f} +/- {t3['se']:.4f}")
    print("   Compare against the simulated M0/M1/M2/M3 values (sims/superarchaic_sim.py).")

    leg = legacy_colocation(df)
    print(f"\n-- legacy co-location test (for comparison only) " + "-" * 25)
    print(f"   {leg['n_cddr']} CDDRs on R_den; {100*leg['frac_shared']:.0f}% show Neanderthal "
          f"co-elevation (simulation puts this test's power under a true Model A near 6%)")

    OUT.mkdir(parents=True, exist_ok=True)
    tag = label.replace(" ", "_").replace("/", "_")
    t1.to_csv(OUT / f"T1_clustering.{tag}.win{win}.tsv", sep="\t", index=False,
              float_format="%.6g")
    t2.to_csv(OUT / f"T2_tailmass.{tag}.win{win}.tsv", sep="\t", index=False,
              float_format="%.6g")
    summary = dict(label=label, n_windows=len(df), n_blocks=n_blocks,
                   rate_denominator=df["_rate_used"].iloc[0],
                   **{f"alpha_{p}": got[p]["alpha"] for p in got},
                   **{f"alpha_se_{p}": got[p]["alpha_se"] for p in got},
                   **{f"ac_{p}": got[p]["autocorr"] for p in got},
                   # per-window abundances: sims/power_calibration.py checks these
                   # against the simulated ones, since power numbers are only
                   # trustworthy if the simulator reproduces the observed spectrum
                   **{p: got[p]["per_window"] for p in got},
                   **{f"contrast_{c['b']}": c["diff"] for c in contrasts},
                   **{f"contrast_z_{c['b']}": c["z"] for c in contrasts},
                   tail_contrast=t2c["contrast"], tail_contrast_se=t2c["se"],
                   tail_contrast_z=t2c["z"], corr_den_nea=t3["corr_den_nea"],
                   corr_den_nea_se=t3["se"], legacy_n_cddr=leg["n_cddr"],
                   legacy_frac_shared=leg["frac_shared"])
    pd.DataFrame([summary]).to_csv(OUT / f"summary.{tag}.win{win}.tsv", sep="\t",
                                   index=False, float_format="%.6g")

    if excluded is not None and len(excluded):
        print(f"\n-- excluded special regions (reported, not tested) " + "-" * 23)
        for name, grp in excluded.groupby("special_region"):
            zr = robust_z(np.concatenate([df[f"R_{x}"].to_numpy(float) for x in LINEAGES]))
            print(f"   {name}: {len(grp)} windows, median div_den_afr="
                  f"{grp['div_den_afr'].median():.5f} vs genome "
                  f"{df['div_den_afr'].median():.5f} "
                  f"({grp['div_den_afr'].median()/df['div_den_afr'].median():.2f}x)")
        excluded.to_csv(OUT / f"special_regions.win{win}.tsv", sep="\t", index=False,
                        float_format="%.6g")
    print(f"\nwrote {OUT}/*.{tag}.win{win}.tsv")
    return summary


def sensitivity(df, win=50000, block_bp=5_000_000):
    """Is the T1 contrast an artefact of data quality, or of the control chosen?

    Two things must be checked before the contrast means anything:

    1. QUALITY. Mappability, repeat content and uneven per-genome callability all
       cluster spatially and could manufacture clustering in pat_nea_all. If the
       contrast survives progressively harsher quality filters, they are not the cause.
    2. CHOICE OF CONTROL. pat_alt_only, pat_vin_only and pat_chag_only are supposed to
       be exchangeable single-Neanderthal drift controls, but Altai split from the
       Vindija/Chagyrskaya ancestor ~150 ka while those two split ~95 ka, so Altai has
       a far longer private branch (~26k private sites against ~4k and ~2k here). alpha
       is poorly determined for the rare ones, and the contrast changes sign depending
       on which is used. Reporting all three is mandatory, not optional.
    """
    rows = []

    def add(name, d):
        if len(d) < 500:
            return
        b = make_blocks(d, block_bp)
        c = contrast_alpha(d, b, "pat_nea_all", "pat_alt_only")
        rows.append(dict(subset=name, n=len(d),
                         alpha_nea_all=relative_dispersion(d["pat_nea_all"], d["n_der_any"]),
                         alpha_alt_only=relative_dispersion(d["pat_alt_only"], d["n_der_any"]),
                         contrast=c["diff"], se=c["se"], z=c["z"]))

    add("all windows", df)
    for t in (0.4, 0.5, 0.6):
        add(f"callable_frac>={t}", df[df["callable_frac"] >= t])
    cf = df[[f"callable_{x}" for x in LINEAGES]].to_numpy(float)
    even = cf.min(axis=1) / np.maximum(cf.max(axis=1), 1)
    add("balanced callability (>=0.9)", df[even >= 0.9])
    if "mean_map20" in df.columns:
        add("non-artifact windows", df[(df["mean_map20"].fillna(1) >= 0.7) &
                                       (df["mean_rm"].fillna(0) <= 0.6)])
    if "repeat_frac" in df.columns:
        add("repeat_frac<0.5", df[df["repeat_frac"].fillna(0) < 0.5])
    sens = pd.DataFrame(rows)

    blocks = make_blocks(df, block_bp)
    ctl = pd.DataFrame([contrast_alpha(df, blocks, "pat_nea_all", c)
                        for c in ("pat_alt_only", "pat_vin_only", "pat_chag_only")])
    print("\n-- sensitivity of the T1 contrast " + "-" * 40)
    print(f"   {'subset':32s} {'n':>6s} {'contrast':>9s} {'z':>7s}")
    for _, r in sens.iterrows():
        print(f"   {r['subset']:32s} {r['n']:6d} {r['contrast']:+9.3f} {r['z']:+7.2f}")
    print("\n   control choice (single-Neanderthal drift patterns are NOT exchangeable):")
    for _, r in ctl.iterrows():
        print(f"     vs {r['b']:15s} {r['diff']:+8.3f} +/- {r['se']:.3f}   z = {r['z']:+6.2f}")
    OUT.mkdir(parents=True, exist_ok=True)
    sens.to_csv(OUT / f"T1_sensitivity.win{win}.tsv", sep="\t", index=False,
                float_format="%.6g")
    ctl.to_csv(OUT / f"T1_control_choice.win{win}.tsv", sep="\t", index=False,
               float_format="%.6g")
    return sens, ctl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+",
                    default=[str(c) for c in range(13, 23)])
    ap.add_argument("--win", type=int, default=50000)
    ap.add_argument("--tables", nargs="+", default=None,
                    help="explicit window-table paths (e.g. simulated ones)")
    ap.add_argument("--label", default=None)
    ap.add_argument("--q", type=float, default=0.99)
    ap.add_argument("--block-mb", type=float, default=5.0)
    ap.add_argument("--keep-special", action="store_true",
                    help="do NOT exclude the MHC (default: excluded, reported separately)")
    ap.add_argument("--sensitivity", action="store_true",
                    help="also run the T1 quality / control-choice sensitivity analysis")
    args = ap.parse_args()
    df, excl = load(args.chroms, args.win, tables=args.tables,
                    exclude_special=not args.keep_special)
    label = args.label or ("sim" if args.tables else "chr" + ",".join(args.chroms))
    run(df, excl, win=args.win, label=label, q=args.q,
        block_bp=int(args.block_mb * 1e6))
    if args.sensitivity:
        sensitivity(df, win=args.win, block_bp=int(args.block_mb * 1e6))


if __name__ == "__main__":
    main()
