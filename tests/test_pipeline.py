"""Unit tests for the superarchaic pipeline core (no downloaded data required).

Run: python -m pytest tests/ -q     (or: python tests/test_pipeline.py)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(SRC))
import vcflib as V           # noqa: E402
import scan_windows as S     # noqa: E402

C = {"A": 0, "C": 1, "G": 2, "T": 3}


# ---- genotype parsing ------------------------------------------------------

def test_sample_call_homref_and_het():
    r = V.Rec("21", 1, "G", [], "AC=0", "GT:DP:GQ", "0/0:10:40")
    al, dp, gq = V.sample_call(r, min_dp=4, min_gq=30)
    assert al == ("G", "G") and dp == 10 and gq == 40.0
    r = V.Rec("21", 2, "G", ["A"], "AC=1", "GT:DP:GQ", "0/1:12:55")
    assert V.sample_call(r)[0] == ("G", "A")


def test_sample_call_filters_and_missing():
    low = V.Rec("21", 3, "G", ["A"], ".", "GT:DP:GQ", "1/1:2:9")
    assert V.sample_call(low, min_dp=4)[0] is None
    miss = V.Rec("21", 4, "G", [], ".", "GT:DP:GQ", "./.:0:0")
    assert V.sample_call(miss)[0] is None


def test_multiallelic_resolves_by_gt_index():
    r = V.Rec("21", 5, "G", ["A", "T"], ".", "GT:DP:GQ", "1/2:20:50")
    assert V.sample_call(r)[0] == ("A", "T")


# ---- divergence algebra (scalar + vectorized agree) ------------------------

def test_pairwise_and_freq_scalar():
    assert V.pairwise_mismatch(("G", "G"), ("A", "A")) == 1.0
    assert V.pairwise_mismatch(("G", "A"), ("G", "A")) == 0.5
    assert abs(V.mismatch_to_freq(("A", "A"), "G", "A", 0.25) - 0.75) < 1e-12
    assert abs(V.mismatch_to_freq(("G", "G"), "G", "A", 0.25) - 0.25) < 1e-12


def test_vectorized_dxy_matches_scalar():
    for x0, x1, y0, y1, exp in [("G", "G", "A", "A", 1.0), ("G", "G", "G", "G", 0.0),
                                ("G", "A", "G", "A", 0.5), ("A", "C", "C", "A", 0.5)]:
        v = S._dxy_arch(np.array([C[x0]]), np.array([C[x1]]),
                        np.array([C[y0]]), np.array([C[y1]]))[0]
        assert abs(v - exp) < 1e-9
    assert np.isnan(S._dxy_arch(np.array([S.MISS]), np.array([S.MISS]),
                                np.array([0]), np.array([0]))[0])


def test_vectorized_dxy_mod():
    v = S._dxy_mod(np.array([C["A"]]), np.array([C["A"]]),
                   np.array([C["G"]]), np.array([C["A"]]), np.array([0.25]))[0]
    assert abs(v - 0.75) < 1e-9
    # monomorphic modern (p=0, alt==ref): a differing archaic mismatches fully
    v = S._dxy_mod(np.array([C["T"]]), np.array([C["T"]]),
                   np.array([C["G"]]), np.array([C["G"]]), np.array([0.0]))[0]
    assert abs(v - 1.0) < 1e-9


# ---- interval / mask math --------------------------------------------------

def test_interval_math():
    a = np.array([[0, 100], [200, 300]], dtype=np.int64)
    b = np.array([[50, 250]], dtype=np.int64)
    assert V.intervals_intersect(a, b).tolist() == [[50, 100], [200, 250]]
    assert V.callable_bp_in_window(a, 0, 1000) == 200
    assert V.callable_bp_in_window(a, 50, 250) == 100
    mem = V.mask_membership(a, [0, 99, 100, 200, 300])
    assert mem.tolist() == [True, True, False, True, False]


def test_win_sum_nan_aware():
    pos0 = np.array([0, 10, 50000, 50001, 100000])
    vals = np.array([1.0, np.nan, 2.0, 3.0, 4.0])
    s, n = S._win_sum(pos0, vals, 50000, 3)
    assert s.tolist() == [1.0, 5.0, 4.0] and n.tolist() == [1, 2, 1]


# ---- vectorized replacements must match the naive implementations -----------

def test_merge_intervals():
    a = np.array([[0, 100], [50, 60], [200, 300], [300, 310], [400, 401]])
    assert V.merge_intervals(a).tolist() == [[0, 100], [200, 310], [400, 401]]
    assert len(V.merge_intervals(np.zeros((0, 2), dtype=np.int64))) == 0


def test_callable_bp_windows_matches_brute_force():
    """The O(n) cumsum path replaced a per-window Python loop; they must agree."""
    rng = np.random.default_rng(0)
    for _ in range(50):
        st = np.sort(rng.integers(0, 1000, 30))
        mask = V.merge_intervals(np.column_stack([st, st + rng.integers(1, 40, 30)]))
        flat = np.zeros(2000, bool)
        for s, e in mask:
            flat[s:e] = True
        win, nb = 137, 12
        got = V.callable_bp_windows(mask, win, nb)
        exp = [int(flat[i * win:(i + 1) * win].sum()) for i in range(nb)]
        assert got.tolist() == exp
        # the scalar entry point must agree with the vectorized one
        assert V.callable_bp_in_window(mask, 0, win) == exp[0]


def test_code_bases_rejects_multichar_alleles():
    """An indel allele like 'AT' must NOT be silently truncated to 'A'."""
    got = V.code_bases(np.array(["A", "C", "G", "T", ".", "AT", "N", "a"], dtype=object))
    assert got.tolist() == [0, 1, 2, 3, -2, -2, -2, 0]


# ---- Model-A statistics ----------------------------------------------------

def test_relative_dispersion_is_mean_scale_free():
    """alpha must not change when a pattern is simply more abundant.

    This is the property phi lacks (phi = 1 + alpha*mu), and the reason cross-pattern
    contrasts are built on alpha rather than on phi.
    """
    import modelA as M
    rng = np.random.default_rng(1)
    n = 4000
    lam = rng.gamma(shape=4.0, scale=1 / 4.0, size=n)      # Var(lambda)/E^2 = 0.25
    off = np.full(n, 100.0)
    a_low = M.relative_dispersion(rng.poisson(2.0 * lam), off)
    a_high = M.relative_dispersion(rng.poisson(20.0 * lam), off)
    assert abs(a_low - 0.25) < 0.08, a_low
    assert abs(a_high - 0.25) < 0.05, a_high
    # phi, by contrast, scales with abundance -- documenting why it is not used
    p_low = M.dispersion(rng.poisson(2.0 * lam), off)
    p_high = M.dispersion(rng.poisson(20.0 * lam), off)
    assert p_high > 2 * p_low


def test_relative_dispersion_zero_for_poisson():
    import modelA as M
    rng = np.random.default_rng(2)
    counts = rng.poisson(5.0, 5000)
    assert abs(M.relative_dispersion(counts, np.full(5000, 100.0))) < 0.02


def test_block_jackknife_recovers_known_se():
    """On independent blocks the block jackknife SE must match the analytic SE."""
    import modelA as M
    rng = np.random.default_rng(3)
    n = 3000
    df = pd.DataFrame({"chrom": "1", "start": np.arange(n) * 50000,
                       "v": rng.normal(0, 1, n)})
    blocks = M.make_blocks(df, block_bp=5_000_000)
    est, _, se = M.block_jackknife(lambda d: float(d["v"].mean()), df, blocks)
    analytic = 1.0 / np.sqrt(n)
    assert abs(est - df["v"].mean()) < 1e-12
    assert 0.5 * analytic < se < 2.0 * analytic, (se, analytic)


def test_binned_residual_is_orthogonal_in_the_tail():
    """The fix for the selection bias: residuals must have ~zero mean in x's tail.

    A linear residual leaves a systematic offset there when y-on-x is curved, which is
    exactly what made the original permutation test report a 25-sigma artefact.
    """
    import significance as sig
    rng = np.random.default_rng(4)
    x = rng.gamma(2.0, 1.0, 8000)
    y = 0.5 * x + 0.3 * x ** 2 + rng.normal(0, 1, 8000)     # deliberately curved
    tail = x >= np.quantile(x, 0.99)
    lin = sig._resid_z(y, x)
    binned = sig._binned_resid_z(y, x)          # adaptive ~50 points per bin
    assert abs(np.nanmean(binned[tail])) < abs(np.nanmean(lin[tail]))
    assert abs(np.nanmean(binned[tail])) < 0.5, np.nanmean(binned[tail])
    # and the coarse binning that motivated the adaptive default must FAIL this,
    # so the regression is caught if someone re-pins nbins to a small constant
    coarse = sig._binned_resid_z(y, x, nbins=20)
    assert abs(np.nanmean(coarse[tail])) > 1.0


def test_polarized_patterns_are_mutually_exclusive():
    """A site can satisfy at most one lineage-specific pattern."""
    rng = np.random.default_rng(5)
    n = 500
    anc = rng.integers(0, 4, n).astype(np.int8)
    u = {"ancc": anc, "refc": anc.copy(),
         "malt": np.full(n, S.MISS, np.int8),
         "p_afr": np.zeros(n, np.float32),
         "in_1000g": np.zeros(n, bool)}
    genc = {}
    for s in ("den", "alt", "vin", "chag"):
        a = anc.copy()
        flip = rng.random(n) < 0.3
        a[flip] = (a[flip] + 1) % 4
        genc[s] = (a, a.copy())
    pats = S.polarized_patterns(u, genc, np.ones(n, bool))
    excl = ["pat_den_only", "pat_nea_all", "pat_alt_only", "pat_vin_only",
            "pat_chag_only", "pat_all_arch"]
    stacked = np.column_stack([pats[k] for k in excl]).sum(axis=1)
    assert stacked.max() <= 1, "patterns overlap"
    # every flagged site must be polarizable and inside the callable base set
    for k in excl:
        assert not np.any(pats[k] & ~pats["n_pol"])


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"OK  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
