"""
superarchaic_sim.py — msprime simulations of Models 0-4 for CDDR calibration.

Simulates the archaic demographic history with a deeply diverged ("superarchaic")
ghost lineage under five models, samples the Denisovan, three Neanderthals (as
ancient tips), and modern African/European panels, then computes the SAME windowed
divergence statistics used on real data (scan_windows). Purpose: set the CDDR
false-positive rate (Model 0), quantify power (Models 1-3), check the ancient-
structure confound (Model 4), and — the crux — measure whether Denisovan and
Neanderthal divergence CO-ELEVATE at the same windows (Model A) or not (Model B).

Models (superarchaic gene flow, backward-in-time mass migrations SOURCE->SUPER):
  M0_none              no superarchaic gene flow (null)
  M1_neandersovan      pulse into the Nea+Den ancestor (NSOV) before they split  -> Model A
  M2_denisovan_only    pulse into the Denisovan lineage after the split          -> Model B
  M3_separate          independent pulses into Denisovan AND into Neanderthals   -> both, not shared-by-descent
  M4_ancient_structure continuous deep migration (structure), no discrete pulse  -> confound

Times in years (converted to generations at 29 yr/gen); demographic values are
literature-guided (Prufer 2014/2017; Rogers 2020; Mafessoni 2020) and deliberately
standard, not fitted.
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import msprime

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "sims"
RATEDIR = ROOT / "results" / "ratemap"
sys.path.insert(0, str(ROOT / "src"))
from vcflib import robust_z

GEN = 29.0
def _g(years):
    return years / GEN

# split times (years)
T_CHIMP = 6_500_000
T_SUPER = 2_000_000      # superarchaic separates from the modern+archaic ancestor
T_HN    = 620_000        # modern (HUM) vs Neandersovan (NSOV)
T_NDSPLIT = 420_000      # Neanderthal (NEA) vs Denisovan (DEN)
T_NEA1  = 150_000        # Altai vs (Vindija,Chagyrskaya) ancestor
T_NEA2  = 95_000         # Vindija vs Chagyrskaya
T_AFREUR = 65_000
T_NEA_INTO_EUR = 50_000
# ancient sampling times (years before present)
S_ALTAI, S_CHAG, S_VIN, S_DEN = 115_000, 80_000, 55_000, 72_000
# effective sizes
NE = dict(CHIMP=10000, SUPER=15000, HUM=12000, AFR=15000, EUR=8000,
          NSOV=3000, NEA=2500, DEN=2500, ALTAI=1500, VIN=1500, CHAG=1500,
          VC=2000, HN=15000, ANC0=15000, ROOT=15000)


def build(model: str, f_super: float = 0.04, f_den: float = 0.03,
          f_nea: float = 0.03, m_struct: float = 2e-5,
          t_hn: float = None, ne_nsov: float = None) -> msprime.Demography:
    """Build the demography. `t_hn` and `ne_nsov` override the module defaults so the
    Neandersovan branch can be fitted (see fit_demography.py) rather than assumed."""
    t_hn = T_HN if t_hn is None else t_hn
    d = msprime.Demography()
    for name, size in NE.items():
        d.add_population(name=name,
                         initial_size=(ne_nsov if name == "NSOV" and ne_nsov else size))
    # Neanderthal->EUR pulse (realism; ~2%)
    d.add_mass_migration(time=_g(T_NEA_INTO_EUR), source="EUR", dest="VIN", proportion=0.02)

    # superarchaic gene flow per model (SOURCE lineages -> SUPER, backward in time)
    if model == "M1_neandersovan":
        # midpoint of the Neandersovan branch. This was hardcoded at 500 ka, which is
        # silently invalid once t_hn drops below it -- NSOV does not exist before t_hn,
        # so a fitted branch length would have placed the pulse outside its own branch.
        d.add_mass_migration(time=_g(0.5 * (T_NDSPLIT + t_hn)), source="NSOV",
                             dest="SUPER", proportion=f_super)
    elif model == "M2_denisovan_only":
        d.add_mass_migration(time=_g(300_000), source="DEN", dest="SUPER", proportion=f_super)
    elif model == "M3_separate":
        d.add_mass_migration(time=_g(300_000), source="DEN", dest="SUPER", proportion=f_den)
        d.add_mass_migration(time=_g(250_000), source="NEA", dest="SUPER", proportion=f_nea)
    elif model == "M4_ancient_structure":
        # continuous deep exchange between SUPER and the NSOV branch (structure, no pulse)
        d.add_migration_rate_change(time=_g(T_NDSPLIT + 10_000), rate=m_struct,
                                    source="NSOV", dest="SUPER")
        d.add_migration_rate_change(time=_g(t_hn), rate=0, source="NSOV", dest="SUPER")
    elif model != "M0_none":
        raise ValueError(model)

    # topology (younger -> older); sort_events() handles ordering.
    # ((VIN,CHAG) at T_NEA2 -> VC), (VC,ALTAI at T_NEA1 -> NEA), (NEA,DEN -> NSOV).
    d.add_population_split(time=_g(T_NEA2), derived=["VIN", "CHAG"], ancestral="VC")
    d.add_population_split(time=_g(T_NEA1), derived=["VC", "ALTAI"], ancestral="NEA")
    d.add_population_split(time=_g(T_NDSPLIT), derived=["NEA", "DEN"], ancestral="NSOV")
    d.add_population_split(time=_g(t_hn), derived=["HUM", "NSOV"], ancestral="HN")
    d.add_population_split(time=_g(T_AFREUR), derived=["AFR", "EUR"], ancestral="HUM")
    d.add_population_split(time=_g(T_SUPER), derived=["SUPER", "HN"], ancestral="ANC0")
    d.add_population_split(time=_g(T_CHIMP), derived=["ANC0", "CHIMP"], ancestral="ROOT")
    d.sort_events()
    return d


# The African panel size is NOT cosmetic. Every polarized pattern is conditioned on the
# derived allele being ABSENT from Africans, and that filter is far easier to pass
# against 20 individuals than against the 661 AFR individuals behind 1000G's AFR_AF.
# With a 20-individual panel the simulated pat_all_arch came out ~8x the observed rate
# purely from panel size. Matching 1000G makes the simulated and real patterns
# comparable, which is the whole point of running both through one statistics module.
N_AFR_1000G, N_EUR_1000G = 661, 503

SAMPLES = [("DEN", 1, S_DEN), ("ALTAI", 1, S_ALTAI), ("VIN", 1, S_VIN),
           ("CHAG", 1, S_CHAG), ("AFR", N_AFR_1000G, 0), ("EUR", 20, 0), ("CHIMP", 1, 0)]


def simulate(model, seq_len=20_000_000, mu=1.25e-8, rho=1e-8, seed=1, f_super=0.04):
    demo = build(model, f_super=f_super, f_den=f_super, f_nea=f_super)
    samples = [msprime.SampleSet(n, population=p, time=_g(t)) for p, n, t in SAMPLES]
    ts = msprime.sim_ancestry(samples=samples, demography=demo, sequence_length=seq_len,
                              recombination_rate=rho, random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed + 1)
    return mts


# ===========================================================================
# Realistic-heterogeneity simulation emitting the REAL window-table schema.
#
# The pilot's simulated null produced a z>=3 window rate of ~0.7% against ~3.5%
# observed, so it could only calibrate qualitatively: it assumed a constant mutation
# rate and constant recombination, and real genomes have neither. Here the per-window
# mutation-rate multipliers are RESAMPLED FROM THE MEASURED human-vs-ancestral
# substitution map (results/ratemap/), so the simulated overdispersion is the observed
# overdispersion rather than a guess. Recombination gets lognormal variation, which
# controls how long a coalescent block stays intact and therefore how much the
# clustering statistics scatter under the null.
#
# The output columns match results/tables/windows.chr*.win*.tsv exactly, so
# src/modelA.py runs on simulated and real data through the same code path -- which is
# what converts the test battery into calibrated power and false-positive rates.
# ===========================================================================

EMPIRICAL_SUB_RATE = 0.0068          # median human-vs-ancestral substitutions/bp


def empirical_rate_multipliers(n_windows, win, rng, fallback_cv=0.25):
    """Per-window mutation-rate multipliers resampled from the measured rate map.

    Falls back to a lognormal with `fallback_cv` when no rate map has been built yet,
    so the simulator still runs on a fresh checkout.
    """
    files = sorted(RATEDIR.glob(f"rate.chr*.win{win}.tsv"))
    vals = []
    for f in files:
        d = pd.read_csv(f, sep="\t", usecols=["sub_rate", "anc_valid_bp"])
        d = d[(d["anc_valid_bp"] > 0.3 * win) & np.isfinite(d["sub_rate"])]
        vals.append(d["sub_rate"].to_numpy())
    if vals:
        v = np.concatenate(vals)
        v = v[np.isfinite(v) & (v > 0)]
        m = np.median(v)
        mult = rng.choice(v / m, size=n_windows, replace=True)
        return mult, "empirical"
    sigma = np.sqrt(np.log(1 + fallback_cv ** 2))
    return rng.lognormal(-0.5 * sigma ** 2, sigma, n_windows), "lognormal-fallback"


def _rate_map(mult, win, base, seq_len):
    pos = np.arange(len(mult) + 1, dtype=float) * win
    pos[-1] = max(seq_len, pos[-2] + 1)
    return msprime.RateMap(position=pos, rate=base * mult)


def simulate_hetero(model, seq_len=30_000_000, win=50_000, mu=1.25e-8, rho=1e-8,
                    seed=1, f_super=0.06, rate_hetero=True, recomb_cv=0.6,
                    t_hn=None, ne_nsov=None):
    """Ancestry + mutations under heterogeneous rates. Returns (mts, rate_multipliers)."""
    rng = np.random.default_rng(seed)
    nwin = int(np.ceil(seq_len / win))
    if rate_hetero:
        mult, _src = empirical_rate_multipliers(nwin, win, rng)
        sigma = np.sqrt(np.log(1 + recomb_cv ** 2))
        rmult = rng.lognormal(-0.5 * sigma ** 2, sigma, nwin)
    else:
        mult = np.ones(nwin)
        rmult = np.ones(nwin)
    demo = build(model, f_super=f_super, f_den=f_super, f_nea=f_super,
                 t_hn=t_hn, ne_nsov=ne_nsov)
    samples = [msprime.SampleSet(n, population=p, time=_g(t)) for p, n, t in SAMPLES]
    ts = msprime.sim_ancestry(
        samples=samples, demography=demo, sequence_length=seq_len,
        recombination_rate=_rate_map(rmult, win, rho, seq_len),
        random_seed=seed, discrete_genome=True)
    mts = msprime.sim_mutations(ts, rate=_rate_map(mult, win, mu, seq_len),
                                random_seed=seed + 100000)
    return mts, mult


def sim_window_table(mts, mult, win, seq_len, seed, chrom="S"):
    """Compute the SAME per-window statistics the real scan produces.

    msprime states are already polarized (0 ancestral, 1 derived), so the polarized
    site patterns are exact here -- which is what makes the simulated values a valid
    calibration for the real ones, where polarization comes from the EPO ancestral call.
    """
    idx = _hap_index(mts)
    pos = np.array([v.site.position for v in mts.variants()])
    G = (mts.genotype_matrix() > 0)
    nbins = int(np.ceil(seq_len / win))
    wid = np.clip((pos // win).astype(int), 0, nbins - 1)

    def dosage(pop):                       # derived alleles carried (0..2) per site
        return G[:, idx[pop]].sum(axis=1)

    def freq(pop):
        return G[:, idx[pop]].mean(axis=1)

    p_afr = freq("AFR")
    p_eur = freq("EUR")
    arch = {"den": "DEN", "alt": "ALTAI", "vin": "VIN", "chag": "CHAG"}
    dos = {k: dosage(v) for k, v in arch.items()}
    nhap = {k: len(idx[v]) for k, v in arch.items()}

    site = {}
    for a, b in [("den", "alt"), ("den", "vin"), ("den", "chag"),
                 ("alt", "vin"), ("alt", "chag"), ("vin", "chag")]:
        pa, pb = dos[a] / nhap[a], dos[b] / nhap[b]
        site[f"div_{a}_{b}"] = pa * (1 - pb) + pb * (1 - pa)
    for s in arch:
        ps = dos[s] / nhap[s]
        site[f"div_{s}_afr"] = ps * (1 - p_afr) + (1 - ps) * p_afr
        site[f"div_{s}_eur"] = ps * (1 - p_eur) + (1 - ps) * p_eur

    der = {s: dos[s] > 0 for s in arch}
    afr0 = p_afr == 0
    NEA_ALL = der["alt"] & der["vin"] & der["chag"]
    pat = {
        "n_pol": np.ones(len(pos), bool),
        "n_der_any": der["den"] | NEA_ALL | der["alt"] | der["vin"] | der["chag"] | (p_afr > 0),
        "pat_den_only": afr0 & der["den"] & ~der["alt"] & ~der["vin"] & ~der["chag"],
        "pat_nea_all": afr0 & ~der["den"] & NEA_ALL,
        "pat_nea_any": afr0 & ~der["den"] & (der["alt"] | der["vin"] | der["chag"]),
        "pat_all_arch": afr0 & der["den"] & NEA_ALL,
        "pat_alt_only": afr0 & ~der["den"] & der["alt"] & ~der["vin"] & ~der["chag"],
        "pat_vin_only": afr0 & ~der["den"] & ~der["alt"] & der["vin"] & ~der["chag"],
        "pat_chag_only": afr0 & ~der["den"] & ~der["alt"] & ~der["vin"] & der["chag"],
    }

    out = {"chrom": [str(chrom)] * nbins,
           "start": np.arange(nbins) * win,
           "end": np.arange(nbins) * win + win}
    for k, v in site.items():
        out[k] = np.bincount(wid, weights=v, minlength=nbins)[:nbins] / win
    for k, v in pat.items():
        out[k] = np.bincount(wid[v], minlength=nbins)[:nbins]

    df = pd.DataFrame(out)
    # Simulated rate denominator, measured the same noisy way as the real one: a
    # Poisson count of human-lineage substitutions in the window, not the true rate.
    rng = np.random.default_rng(seed + 999)
    lam = mult[:nbins] * EMPIRICAL_SUB_RATE * win
    hum_sub = rng.poisson(lam)
    df["anc_valid_bp"] = win
    df["hum_sub"] = hum_sub
    df["sub_rate"] = hum_sub / win
    df["true_rate_mult"] = mult[:nbins]
    df["callable_common"] = win
    df["callable_frac"] = 1.0
    df["callable_frac_den"] = 1.0
    for s in arch:
        df[f"callable_{s}"] = win
    for c in ("mean_map20", "mean_rm", "mean_ur", "mean_bsc"):
        df[c] = np.nan
    df["repeat_frac"] = 0.0
    return df


def _one_rep(job):
    model, seed, seq_len, win, f_super, rate_hetero = job
    mts, mult = simulate_hetero(model, seq_len=seq_len, win=win, seed=seed,
                                f_super=f_super, rate_hetero=rate_hetero)
    return sim_window_table(mts, mult, win, seq_len, seed, chrom=f"{seed}")


def battery(models=("M0_none", "M1_neandersovan", "M2_denisovan_only",
                    "M3_separate", "M4_ancient_structure"),
            seeds=(1, 2, 3), seq_len=30_000_000, win=50_000, f_super=0.06,
            rate_hetero=True, jobs=1, tag=""):
    """Run each model through the REAL statistics battery (src/modelA.py).

    Each seed is an independent 'chromosome', so the block jackknife and the
    per-chromosome circular shifts behave as they do on real data.
    """
    sys.path.insert(0, str(ROOT / "src"))
    import modelA

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for m in models:
        jobs_list = [(m, sd, seq_len, win, f_super, rate_hetero) for sd in seeds]
        if jobs > 1 and len(jobs_list) > 1:
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                dfs = list(ex.map(_one_rep, jobs_list))
        else:
            dfs = [_one_rep(j) for j in jobs_list]
        df = pd.concat(dfs, ignore_index=True)
        p = OUT / f"windows_{m}{tag}.tsv"
        df.to_csv(p, sep="\t", index=False, float_format="%.6g")
        sub, _ = modelA.load(tables=[p], win=win, exclude_special=False)
        summ = modelA.run(sub, None, win=win, label=f"{m}{tag}")
        summ["model"] = m
        summ["f_super"] = f_super
        rows.append(summ)
    s = pd.DataFrame(rows)
    s.to_csv(OUT / f"battery_summary{tag}.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"\nwrote {OUT/f'battery_summary{tag}.tsv'}")
    return s


def _hap_index(mts):
    """Map population name -> array of haplotype (sample node) column indices."""
    pop_id = {p.metadata.get("name", p.id) if p.metadata else p.id: p.id
              for p in mts.populations()}
    name_by_id = {v: k for k, v in pop_id.items()}
    idx = {}
    for j, node in enumerate(mts.samples()):
        nm = name_by_id[mts.node(node).population]
        idx.setdefault(nm, []).append(j)
    return {k: np.array(v) for k, v in idx.items()}


def window_divergence(mts, win=50_000):
    idx = _hap_index(mts)
    pos = np.array([v.site.position for v in mts.variants()])
    G = mts.genotype_matrix().astype(np.int8)            # (n_sites, n_hap), 0=anc 1=der
    seq_len = int(mts.sequence_length)

    def freq(pop):
        return G[:, idx[pop]].mean(axis=1)

    def dxy_pop_pop(a, b):
        # mean over hap pairs of allele mismatch = pa(1-pb)+pb(1-pa)
        pa, pb = freq(a), freq(b)
        return pa * (1 - pb) + pb * (1 - pa)

    site = {
        "den_afr": dxy_pop_pop("DEN", "AFR"),
        "alt_afr": dxy_pop_pop("ALTAI", "AFR"),
        "vin_afr": dxy_pop_pop("VIN", "AFR"),
        "chag_afr": dxy_pop_pop("CHAG", "AFR"),
        "den_alt": dxy_pop_pop("DEN", "ALTAI"),
        "den_vin": dxy_pop_pop("DEN", "VIN"),
        "den_chag": dxy_pop_pop("DEN", "CHAG"),
        "alt_vin": dxy_pop_pop("ALTAI", "VIN"),
    }
    nbins = seq_len // win
    wid = (pos // win).astype(int)
    wid[wid >= nbins] = nbins - 1
    rows = {"start": np.arange(nbins) * win}
    for k, v in site.items():
        s = np.bincount(wid, weights=v, minlength=nbins)[:nbins]
        rows[k] = s / win
    df = pd.DataFrame(rows)
    df["nea_afr"] = df[["alt_afr", "vin_afr", "chag_afr"]].mean(axis=1)
    df["den_nea"] = df[["den_alt", "den_vin", "den_chag"]].mean(axis=1)
    return df


def analyze(models=("M0_none", "M1_neandersovan", "M2_denisovan_only",
                    "M3_separate", "M4_ancient_structure"),
            seq_len=20_000_000, win=50_000, z_thr=3.0, seed=1):
    OUT.mkdir(parents=True, exist_ok=True)
    summ = []
    for m in models:
        mts = simulate(m, seq_len=seq_len, seed=seed)
        df = window_divergence(mts, win=win)
        df["z_den"] = robust_z(df["den_afr"].to_numpy())
        df["z_nea"] = robust_z(df["nea_afr"].to_numpy())
        df["cddr"] = df["z_den"] >= z_thr
        # Model A vs B discriminator: do CDDR windows also elevate Neanderthal divergence?
        co = df.loc[df["cddr"], "z_nea"]
        corr = float(np.corrcoef(df["z_den"], df["z_nea"])[0, 1]) if len(df) > 2 else np.nan
        summ.append(dict(
            model=m, n_win=len(df),
            cddr_rate=float(df["cddr"].mean()),
            median_den_afr=float(df["den_afr"].median()),
            median_nea_afr=float(df["nea_afr"].median()),
            mean_z_nea_at_cddr=float(co.mean()) if len(co) else np.nan,
            frac_cddr_with_nea_elev=float((co >= 1.0).mean()) if len(co) else np.nan,
            corr_zden_znea=corr,
        ))
        df.to_csv(OUT / f"windows_{m}.tsv", sep="\t", index=False, float_format="%.6g")
        print(f"{m:22s} cddr_rate={summ[-1]['cddr_rate']:.4f}  "
              f"z_nea@cddr={summ[-1]['mean_z_nea_at_cddr']:.2f}  "
              f"corr(zden,znea)={corr:.2f}")
    s = pd.DataFrame(summ)
    s.to_csv(OUT / "model_summary.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"\nwrote {OUT/'model_summary.tsv'}")
    return s


def _resid_z(y, x):
    from numpy import polyfit
    ok = np.isfinite(y) & np.isfinite(x)
    b1, b0 = polyfit(x[ok], y[ok], 1)
    return robust_z(y - (b0 + b1 * x))


def discriminate(models=("M0_none", "M1_neandersovan", "M2_denisovan_only",
                         "M3_separate", "M4_ancient_structure"),
                 seeds=(1, 2, 3), seq_len=30_000_000, win=50_000,
                 f_super=0.06, z_cddr=3.0, z_hi=2.0, z_lo=1.0):
    """Does the CONFOUND-CONTROLLED discriminator separate Model A from Model B?

    Pools windows across seeds. Flags CDDRs on Denisovan-vs-African divergence, then
    classifies each CDDR by the residual excess of Neanderthal-vs-African depth given
    Denisovan depth (shared-ancestry signature, Model A) vs Denisovan-vs-Neanderthal
    excess (Denisovan-only signature, Model B). Reports the mix per model + the M0 FPR.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for m in models:
        dfs = []
        for sd in seeds:
            mts = simulate(m, seq_len=seq_len, seed=sd, f_super=f_super)
            dfs.append(window_divergence(mts, win=win))
        df = pd.concat(dfs, ignore_index=True)
        df["z_den"] = robust_z(df["den_afr"].to_numpy())
        df["exc_nea"] = _resid_z(df["nea_afr"].to_numpy(), df["den_afr"].to_numpy())
        df["z_den_nea"] = robust_z(df["den_nea"].to_numpy())
        cddr = df[df["z_den"] >= z_cddr]
        shared = ((cddr["exc_nea"] >= z_hi) & (cddr["z_den_nea"] <= z_lo)).sum()
        denonly = ((cddr["z_den_nea"] >= z_hi) & (cddr["exc_nea"] < z_lo)).sum()
        n = len(cddr)
        rows.append(dict(model=m, n_win=len(df), n_cddr=n,
                         cddr_rate=n / len(df),
                         pct_shared=100 * shared / n if n else np.nan,
                         pct_denonly=100 * denonly / n if n else np.nan))
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "discriminator_summary.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"\nCONFOUND-CONTROLLED DISCRIMINATOR (f_super={f_super}, seeds={list(seeds)}, "
          f"{seq_len/1e6:.0f}Mb x{len(seeds)})")
    print(f"{'model':22s} {'n_cddr':>7s} {'cddr_rate':>9s} {'%shared':>8s} {'%den-only':>9s}")
    for r in rows:
        print(f"{r['model']:22s} {r['n_cddr']:7d} {r['cddr_rate']:9.4f} "
              f"{r['pct_shared']:8.0f} {r['pct_denonly']:9.0f}")
    print("\nExpectation: M1(Model A) high %shared; M2(Model B) high %den-only; "
          "M0 low CDDR rate (= FPR).")
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["analyze", "discriminate"], default="analyze")
    ap.add_argument("--seq-len", type=int, default=20_000_000)
    ap.add_argument("--win", type=int, default=50_000)
    ap.add_argument("--z-thr", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3])
    ap.add_argument("--f-super", type=float, default=0.06)
    ap.add_argument("--models", nargs="+", default=None)
    args = ap.parse_args()
    if args.mode == "discriminate":
        kw = dict(seeds=tuple(args.seeds), seq_len=args.seq_len, win=args.win,
                  f_super=args.f_super)
        if args.models:
            kw["models"] = tuple(args.models)
        discriminate(**kw)
    else:
        kw = dict(seq_len=args.seq_len, win=args.win, z_thr=args.z_thr, seed=args.seed)
        if args.models:
            kw["models"] = tuple(args.models)
        analyze(**kw)


if __name__ == "__main__":
    main()
