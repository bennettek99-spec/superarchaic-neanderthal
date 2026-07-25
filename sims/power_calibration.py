"""
power_calibration.py — replicate-level false-positive rate and power for the Model-A battery.

The single-pass battery in superarchaic_sim.py gives ONE value of each test statistic
per model, with a block-jackknife SE computed inside that one dataset. That is enough
to see the direction of an effect and nothing more: it cannot say how often the null
produces a value that large, which is the only thing that licenses a claim.

This module simulates many INDEPENDENT datasets per model, each sized to the real
analysis, and records the test statistics for each. From those distributions:

  threshold  = 95th percentile of the statistic under M0 (the null)
  FPR        = P(statistic > threshold | M0)          -- 0.05 by construction, a check
  power(Mk)  = P(statistic > threshold | Mk)          -- what the pilot could not state
  p_obs(Mk)  = where the OBSERVED value falls in model Mk's distribution

That last column is the one that matters: it says which models can and cannot produce
what the real chromosomes show, rather than which model the point estimate happens to
sit nearest.

Two limitations this run does NOT remove, both of which belong in any write-up:
  * M1 (Model A) and M4 (ancient structure) produce very similar clustering. Separating
    them needs an explicit structure model fitted to the data, not more replicates.
  * The simulated polarized-pattern ABUNDANCES do not yet match the observed ones
    (pat_all_arch is several-fold too high), which means the demography -- chiefly the
    Neandersovan branch length and Ne -- is not yet fitted. --report-abundance prints
    the comparison so the gap is visible; treat the power numbers as calibrated in
    shape, not in absolute scale, until it is closed.

Usage (overnight):
  python sims/power_calibration.py --reps 20 --jobs 6
  python sims/power_calibration.py --reps 20 --f-super 0.02 0.04 0.06 0.10   # power curve
"""
from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "sims"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "sims"))

MODELS = ("M0_none", "M1_neandersovan", "M2_denisovan_only",
          "M3_separate", "M4_ancient_structure")

# Statistics carried through from src/modelA.py. The first is the headline.
STATS = ["contrast_pat_alt_only", "corr_den_nea", "alpha_pat_nea_all",
         "alpha_pat_den_only", "tail_contrast_z", "legacy_frac_shared"]


def _one_replicate(job):
    """Simulate one whole dataset and return its battery summary as a dict."""
    model, rep, f_super, n_chunks, chunk_bp, win = job
    import superarchaic_sim as S
    import modelA
    dfs = []
    for k in range(n_chunks):
        seed = 1_000_000 * (rep + 1) + 1000 * k + hash(model) % 997
        mts, mult = S.simulate_hetero(model, seq_len=chunk_bp, win=win, seed=seed,
                                      f_super=f_super)
        d = S.sim_window_table(mts, mult, win, chunk_bp, seed, chrom=f"{k}")
        dfs.append(d)
    df = pd.concat(dfs, ignore_index=True)
    # modelA.load() reads from disk; build the identical frame in memory instead so a
    # replicate never touches the filesystem (hundreds of these run in parallel)
    sub = _prepare_inmemory(df, modelA)
    blocks = modelA.make_blocks(sub, 5_000_000)
    t1 = modelA.T1_clade_clustering(sub, blocks, win=win)
    got = {r["pattern"]: r for _, r in t1.iterrows()}
    c = modelA.contrast_alpha(sub, blocks, "pat_nea_all", "pat_alt_only")
    _, t2c = modelA.T2_tail_symmetry(sub, blocks)
    t3 = modelA.T3_cross_lineage_cov(sub, blocks)
    leg = modelA.legacy_colocation(sub)
    return dict(model=model, rep=rep, f_super=f_super, n_windows=len(sub),
                contrast_pat_alt_only=c["diff"], contrast_z=c["z"],
                corr_den_nea=t3["corr_den_nea"],
                alpha_pat_nea_all=got["pat_nea_all"]["alpha"],
                alpha_pat_alt_only=got["pat_alt_only"]["alpha"],
                alpha_pat_den_only=got["pat_den_only"]["alpha"],
                alpha_pat_all_arch=got["pat_all_arch"]["alpha"],
                mean_pat_nea_all=got["pat_nea_all"]["per_window"],
                mean_pat_all_arch=got["pat_all_arch"]["per_window"],
                mean_pat_den_only=got["pat_den_only"]["per_window"],
                tail_contrast_z=t2c["z"], legacy_frac_shared=leg["frac_shared"])


def _prepare_inmemory(df, modelA):
    """Attach modelA's derived columns to an in-memory simulated window table."""
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    sr = df["sub_rate"].to_numpy(float)
    sr = np.where(np.isfinite(sr) & (sr > 0), sr, np.nan)
    for x in modelA.LINEAGES:
        df[f"R_{x}"] = df[f"div_{x}_afr"].to_numpy(float) / sr
    df["R_nea"] = df[["R_alt", "R_vin", "R_chag"]].mean(axis=1)
    df["_rate_used"] = "sub_rate"
    return df[df["n_der_any"] >= 30].sort_values(["chrom", "start"]).reset_index(drop=True)


def run(models=MODELS, reps=20, f_supers=(0.06,), n_chunks=15, chunk_bp=20_000_000,
        win=50_000, jobs=4, observed=None):
    OUT.mkdir(parents=True, exist_ok=True)
    jobs_list = [(m, r, f, n_chunks, chunk_bp, win)
                 for m in models for f in f_supers for r in range(reps)
                 if not (m == "M0_none" and f != f_supers[0])]   # M0 ignores f_super
    print(f"{len(jobs_list)} replicate datasets "
          f"({n_chunks}x{chunk_bp/1e6:.0f} Mb = {n_chunks*chunk_bp/1e6:.0f} Mb each), "
          f"{jobs} workers")
    t0 = time.time()
    rows = []
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for i, r in enumerate(ex.map(_one_replicate, jobs_list), 1):
                rows.append(r)
                if i % 5 == 0:
                    el = time.time() - t0
                    print(f"  {i}/{len(jobs_list)} done ({el/60:.1f} min, "
                          f"eta {el/i*(len(jobs_list)-i)/60:.1f} min)", flush=True)
    else:
        for j in jobs_list:
            rows.append(_one_replicate(j))
    reps_df = pd.DataFrame(rows)
    reps_df.to_csv(OUT / "power_replicates.tsv", sep="\t", index=False, float_format="%.6g")

    print("\n" + "=" * 78)
    print("REPLICATE-LEVEL CALIBRATION")
    print("=" * 78)
    summ = []
    for stat in STATS:
        if stat not in reps_df.columns:
            continue
        null = reps_df[(reps_df.model == "M0_none")][stat].dropna()
        if len(null) < 5:
            continue
        thr = float(np.quantile(null, 0.95))
        print(f"\n{stat}   (M0 95th percentile = {thr:+.4f})")
        obs = observed.get(stat) if observed else None
        if obs is not None:
            print(f"   observed = {obs:+.4f}")
        print(f"   {'model':22s} {'mean':>9s} {'sd':>8s} {'power':>7s}  {'P(sim>=obs)':>12s}")
        for m in models:
            v = reps_df[(reps_df.model == m)][stat].dropna()
            if not len(v):
                continue
            power = float(np.mean(v > thr))
            pobs = float(np.mean(v >= obs)) if obs is not None else np.nan
            tag = "  <- FPR" if m == "M0_none" else ""
            print(f"   {m:22s} {v.mean():+9.4f} {v.std():8.4f} {power:7.2f}  "
                  f"{pobs:12.2f}{tag}")
            summ.append(dict(statistic=stat, model=m, mean=v.mean(), sd=v.std(),
                             threshold=thr, power=power, observed=obs, p_obs=pobs,
                             n_reps=len(v)))
    s = pd.DataFrame(summ)
    s.to_csv(OUT / "power_summary.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"\nwrote {OUT/'power_replicates.tsv'} and power_summary.tsv "
          f"({(time.time()-t0)/60:.1f} min)")
    return reps_df, s


def report_abundance(reps_df, observed_table=None):
    """Compare simulated and observed polarized-pattern abundances.

    A power number is only trustworthy if the simulator reproduces the data it is
    calibrating. Large abundance gaps mean the demography needs fitting first.
    """
    cols = ["mean_pat_nea_all", "mean_pat_all_arch", "mean_pat_den_only"]
    print("\n-- simulated vs observed pattern abundance (per 50 kb window) --")
    sim = reps_df[reps_df.model == "M0_none"][cols].mean()
    for c in cols:
        line = f"   {c:22s} sim {sim[c]:8.2f}"
        if observed_table is not None:
            o = observed_table.get(c.replace("mean_", ""))
            if o is not None:
                line += f"   observed {o:8.2f}   ratio {sim[c]/o:5.2f}x"
        print(line)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--models", nargs="+", default=list(MODELS))
    ap.add_argument("--f-super", nargs="+", type=float, default=[0.06])
    ap.add_argument("--n-chunks", type=int, default=15,
                    help="independent 'chromosomes' per replicate dataset")
    ap.add_argument("--chunk-mb", type=float, default=20.0)
    ap.add_argument("--win", type=int, default=50000)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--observed", default=None,
                    help="results/modelA/summary.*.tsv to locate the real value in each "
                         "model's simulated distribution")
    args = ap.parse_args()
    obs = None
    obs_tab = None
    if args.observed and Path(args.observed).exists():
        o = pd.read_csv(args.observed, sep="\t").iloc[0].to_dict()
        obs = {k: o[k] for k in STATS if k in o}
        obs_tab = o
    reps_df, _ = run(models=tuple(args.models), reps=args.reps,
                     f_supers=tuple(args.f_super), n_chunks=args.n_chunks,
                     chunk_bp=int(args.chunk_mb * 1e6), win=args.win,
                     jobs=args.jobs, observed=obs)
    report_abundance(reps_df, obs_tab)


if __name__ == "__main__":
    main()
