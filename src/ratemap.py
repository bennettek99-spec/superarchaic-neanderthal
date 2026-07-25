"""
ratemap.py — per-window LOCAL MUTATION RATE from human-vs-ancestral substitutions.

Why this module exists (the central statistical fix)
----------------------------------------------------
Per-window mutation-rate heterogeneity inflates EVERY archaic-vs-modern divergence at
a locus together. The pilot handled that by contrasting genomes against each other:

    den_excess = div_den_afr - mean(div_alt_afr, div_vin_afr, div_chag_afr)

That cancels the shared rate -- but it also cancels Model A exactly, because under
Model A (superarchaic gene flow into the Neanderthal+Denisovan ancestor) BOTH terms
rise together. Model A was structurally invisible to the pilot's headline statistic,
which is the real reason the co-location test looked so underpowered.

The fix is an EXTERNAL rate denominator that does not involve the archaics at all:
the density of substitutions on the human lineage since the human-chimp ancestor,
i.e. positions where the hg19 reference differs from the Ensembl EPO 6-primate
ancestral sequence. That accumulates over ~6 Myr, so at 50 kb it is ~300 events per
window (~6% Poisson noise) versus ~72 events for div_den_afr -- a near-noiseless
local rate estimate. Dividing archaic depth by it removes mutation-rate variation
while leaving Model A's signal intact.

Only high-confidence ancestral calls count: the EPO FASTA writes a base UPPERCASE
when all outgroups agree and lowercase when the call is weaker, so case is
meaningful and is preserved. hg19's own lowercase is RepeatMasker soft-masking and
is case-folded for the substitution count -- but it is also tallied separately,
which yields a free genome-wide per-window repeat fraction (the exact
RepeatMasker track docs/FINDINGS.md listed as missing).

Everything is restricted to the COMMON callable mask (intersection of the four
archaic FilterBeds) so the denominator covers the same sites as the numerators.

Laptop-friendly acquisition, matching src/stagger.sh: the hg19 chromosome FASTA is
downloaded (~50-90 MB gz), reduced to per-window counts, then DELETED. The EPO
ancestral FASTAs are extracted once from the tarball already in data/ancestral/.

Usage:
  python src/ratemap.py --chroms 6 13 14 ...        # download -> count -> delete
  python src/ratemap.py --chroms 6 --keep-fasta     # keep the hg19 FASTA
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vcflib as V

ROOT = Path(__file__).resolve().parent.parent
ANC = ROOT / "data" / "ancestral"
MASKS = ROOT / "data" / "masks"
REFDIR = ROOT / "data" / "hg19"
OUT = ROOT / "results" / "ratemap"

ANC_TARBALL = ANC / "homo_sapiens_ancestor_GRCh37_e71.tar.bz2"
ANC_DIR = ANC / "homo_sapiens_ancestor_GRCh37_e71"
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr{c}.fa.gz"

MASK_STEM = ["Denisova", "Altai", "Vindija33.19", "Chagyrskaya"]

# Ancestral codes: case IS confidence, so uppercase-only maps to a base.
_ANC_LUT = np.full(256, V.MISS, dtype=np.int8)
for _b, _c in zip(b"ACGT", range(4)):
    _ANC_LUT[_b] = _c


def _read_fasta_bytes(path, want_chrom=None) -> np.ndarray:
    """Read one chromosome's sequence from a (possibly gz/bz2) FASTA as raw uint8.

    Case is PRESERVED -- callers decide whether it means confidence (EPO ancestral)
    or soft-masking (hg19). If want_chrom is given, only that record is returned.
    """
    buf = bytearray()
    grabbing = want_chrom is None
    want = None if want_chrom is None else str(want_chrom).replace("chr", "")
    with V.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if buf:
                    break
                if want is not None:
                    toks = line[1:].strip().replace(":", " ").split()
                    grabbing = want in toks or f"chr{want}" in toks
                continue
            if grabbing:
                buf.extend(line.rstrip("\n").encode("ascii", "replace"))
    return np.frombuffer(bytes(buf), dtype=np.uint8)


def ancestral_path(chrom: str) -> Path:
    """Path to the EPO ancestral FASTA for one chromosome, extracting it on demand.

    bzip2 is not seekable, so pulling one member still decompresses the whole 800 MB
    archive. Extract EVERY chromosome on the first miss rather than paying that cost
    once per chromosome.
    """
    hits = list(ANC.glob(f"**/homo_sapiens_ancestor_{chrom}.fa"))
    if hits:
        return hits[0]
    if not ANC_TARBALL.exists():
        raise FileNotFoundError(
            f"{ANC_TARBALL} not found. Fetch the Ensembl EPO ancestral set for GRCh37 "
            f"(release 71) into data/ancestral/ first.")
    print("  extracting the ancestral tarball (one-off, all chromosomes, ~10 min)...",
          flush=True)
    with tarfile.open(ANC_TARBALL, "r:bz2") as tf:
        tf.extractall(ANC, filter="data")
    hits = list(ANC.glob(f"**/homo_sapiens_ancestor_{chrom}.fa"))
    if not hits:
        raise FileNotFoundError(f"chr{chrom} not present inside {ANC_TARBALL.name}")
    return hits[0]


def fetch_hg19(chrom: str) -> Path:
    """Download the hg19 chromosome FASTA from UCSC (resumable) if not present."""
    REFDIR.mkdir(parents=True, exist_ok=True)
    out = REFDIR / f"chr{chrom}.fa.gz"
    if out.exists() and out.stat().st_size > 1_000_000:
        return out
    url = UCSC.format(c=chrom)
    print(f"  downloading {url}", flush=True)
    subprocess.run(["curl", "-fL", "-C", "-", "--retry", "8", "--retry-delay", "5",
                    "--retry-all-errors", "--connect-timeout", "30",
                    "--no-progress-meter", "-o", str(out), url], check=True)
    return out


def common_mask(chrom: str) -> np.ndarray:
    m = V.read_mask(MASKS / f"{MASK_STEM[0]}.chr{chrom}.mask.bed.gz", chrom)
    for stem in MASK_STEM[1:]:
        m = V.intervals_intersect(
            m, V.read_mask(MASKS / f"{stem}.chr{chrom}.mask.bed.gz", chrom))
    return m


def mask_to_flags(mask: np.ndarray, length: int) -> np.ndarray:
    """Expand disjoint intervals to a per-base boolean array of `length`."""
    d = np.zeros(length + 1, dtype=np.int8)
    if len(mask):
        s = np.clip(mask[:, 0], 0, length)
        e = np.clip(mask[:, 1], 0, length)
        np.add.at(d, s, 1)
        np.add.at(d, e, -1)
    return np.cumsum(d)[:length].astype(bool)


def window_sum(x: np.ndarray, win: int, nbins: int) -> np.ndarray:
    """Sum a long boolean/numeric array into `nbins` fixed-width windows.

    Reduces the full-width rows in place rather than building a chromosome-length
    int array, which matters when x has ~171 million entries (chr6).
    """
    out = np.zeros(nbins, dtype=np.int64)
    full = min(len(x) // win, nbins)
    if full:
        out[:full] = x[:full * win].reshape(full, win).sum(axis=1, dtype=np.int64)
    if full < nbins and full * win < len(x):
        out[full] = x[full * win:].sum(dtype=np.int64)
    return out


def build(chrom: str, wins=(20000, 50000, 100000), keep_fasta=False) -> dict:
    t0 = time.time()
    print(f"chr{chrom}: building rate map")
    anc_raw = _read_fasta_bytes(ancestral_path(chrom))
    ref_path = fetch_hg19(chrom)
    ref_raw = _read_fasta_bytes(ref_path, want_chrom=chrom)
    L = min(len(anc_raw), len(ref_raw))
    if len(anc_raw) != len(ref_raw):
        print(f"  note: ancestral {len(anc_raw):,} bp vs hg19 {len(ref_raw):,} bp; "
              f"using the common prefix ({L:,} bp)")
    anc_raw, ref_raw = anc_raw[:L], ref_raw[:L]

    anc = _ANC_LUT[anc_raw]                      # uppercase-only == high confidence
    ref = V._BASE_LUT[ref_raw]                   # case-folded (hg19 case == RepeatMasker)
    repeat = (ref_raw >= ord("a")) & (ref_raw <= ord("z"))
    del anc_raw, ref_raw

    inmask = mask_to_flags(common_mask(chrom), L)
    valid = inmask & (anc >= 0) & (ref >= 0)
    sub = valid & (anc != ref)
    del anc, ref

    out = {}
    OUT.mkdir(parents=True, exist_ok=True)
    for win in wins:
        nbins = L // win + 1
        n_valid = window_sum(valid, win, nbins)
        n_sub = window_sum(sub, win, nbins)
        n_rep = window_sum(repeat, win, nbins)
        n_mask = window_sum(inmask, win, nbins)
        df = pd.DataFrame({
            "chrom": str(chrom),
            "start": (np.arange(nbins) * win).astype(np.int64),
            "end": (np.arange(nbins) * win + win).astype(np.int64),
            "anc_valid_bp": n_valid,        # callable AND confidently polarizable
            "hum_sub": n_sub,               # hg19 != ancestral among those
            "repeat_bp": n_rep,             # RepeatMasker soft-masked (free, genome-wide)
            "mask_bp": n_mask,
        })
        # local substitution rate per bp; the Model-A-preserving rate denominator
        df["sub_rate"] = df["hum_sub"] / df["anc_valid_bp"].replace(0, np.nan)
        df["repeat_frac"] = df["repeat_bp"] / win
        df["anc_cov"] = df["anc_valid_bp"] / df["mask_bp"].replace(0, np.nan)
        p = OUT / f"rate.chr{chrom}.win{win}.tsv"
        df.to_csv(p, sep="\t", index=False, float_format="%.6g")
        out[win] = df
        good = df["anc_valid_bp"] > 0
        print(f"  win{win}: {int(good.sum())} windows with ancestral coverage, "
              f"median sub_rate={np.nanmedian(df.loc[good,'sub_rate']):.5f}, "
              f"median anc_cov={np.nanmedian(df.loc[good,'anc_cov']):.2f} -> {p.name}")
    if not keep_fasta:
        ref_path.unlink(missing_ok=True)
        print(f"  deleted {ref_path.name} (rate map retained)")
    print(f"  chr{chrom} done in {time.time()-t0:.0f}s")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", required=True)
    ap.add_argument("--wins", nargs="+", type=int, default=[20000, 50000, 100000])
    ap.add_argument("--keep-fasta", action="store_true",
                    help="keep the downloaded hg19 FASTA (default: delete after counting)")
    args = ap.parse_args()
    for c in args.chroms:
        try:
            build(c, wins=tuple(args.wins), keep_fasta=args.keep_fasta)
        except Exception as e:                       # keep going across chromosomes
            print(f"chr{c}: FAILED — {type(e).__name__}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
