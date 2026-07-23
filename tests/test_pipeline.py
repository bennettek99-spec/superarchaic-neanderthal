"""Unit tests for the superarchaic pipeline core (no downloaded data required).

Run: python -m pytest tests/ -q     (or: python tests/test_pipeline.py)
"""
import sys
from pathlib import Path

import numpy as np

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


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"OK  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
