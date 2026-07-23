"""
vcflib.py — dependency-light streaming readers for the high-coverage archaic VCFs.

No bcftools / cyvcf2 / pysam in this environment, so everything here is pure Python
(gzip + numpy). Design notes specific to the Max Planck EVA archaic VCFs:

* Altai Neanderthal / Denisovan ("*.mod.vcf.gz", Prufer 2014) are ALL-SITES VCFs.
  FILTER is "." — do NOT trust it. Callability comes from the per-sample GT/DP/GQ and,
  preferably, the external FilterBed callability masks. Their INFO carries rich
  Ensembl-EPO annotation at aligned sites: CAnc (chimp-human ancestor base = the
  polarizing ancestral allele), Map20 (Duke 20-mer mappability), RM (repeat-masked
  flag), UR (copy-number/segdup control region), SysErr/SysErrHCB (systematic error),
  CpG, mSC/pSC/GRP/bSC (conservation / background-selection scores).
* Vindija33.19 ("chrN_mq25_mapab100.vcf.gz") and Chagyrskaya ("chrN.noRB.vcf.gz")
  are already mappability/MQ filtered; schema is checked at read time, not assumed.
* 1000 Genomes phase3 INFO already exposes AFR_AF / EUR_AF / EAS_AF etc., so modern
  allele frequencies are read straight from INFO — no need to genotype 2504 samples.

Efficiency strategy the rest of the pipeline relies on:
  - callable-site DENOMINATORS per window come from BED masks (interval arithmetic),
  - divergence NUMERATORS are accumulated by scanning VARIANT sites only.
This avoids parsing tens of millions of invariant lines per genome.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import Iterator

import numpy as np

# ----------------------------------------------------------------------------- IO


def open_text(path):
    """Open a possibly-gzipped text file for line iteration."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return open(path, "r", encoding="ascii", errors="replace")


def parse_info(info: str) -> dict:
    """Parse a VCF INFO string into {key: value}; bare flags map to True."""
    d = {}
    if info == "." or not info:
        return d
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        elif field:
            d[field] = True
    return d


# ------------------------------------------------------------------- record model


@dataclass(slots=True)
class Rec:
    chrom: str
    pos: int
    ref: str
    alt: list      # list of ALT allele strings ([] if ".")
    info: str      # raw INFO string (parse lazily via parse_info)
    fmt: str       # FORMAT string
    sample: str    # the (single) sample field for archaic VCFs

    def alleles(self):
        """Return the allele table: index 0 = REF, 1.. = ALT alleles."""
        return [self.ref] + self.alt


def iter_records(path, keep_only_variant: bool = False) -> Iterator[Rec]:
    """Stream a single-sample VCF. If keep_only_variant, skip lines with ALT '.'.

    Yields Rec for the FIRST sample column only (archaic VCFs are single-sample).
    """
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            alt_field = f[4]
            if keep_only_variant and (alt_field == "." or alt_field == ""):
                continue
            alt = [] if alt_field == "." else alt_field.split(",")
            yield Rec(f[0], int(f[1]), f[3], alt, f[7], f[8], f[9])


# --------------------------------------------------------------- genotype parsing


def _fmt_index(fmt: str, key: str) -> int:
    keys = fmt.split(":")
    return keys.index(key) if key in keys else -1


def sample_call(rec: Rec, min_dp: int = 0, min_gq: float = 0.0):
    """Return (allele_bases, dp, gq) for a single-sample record, or (None, dp, gq).

    allele_bases is a tuple of 1 or 2 base strings drawn from REF/ALT by the GT
    indices (e.g. ('G','A') for a 0/1 het where REF=G, ALT=A). Returns None alleles
    when GT is missing, or when DP/GQ are below thresholds (treated as not-callable).
    Phasing '|' vs '/' is ignored. For an all-sites record with ALT '.' and GT 0/0,
    alleles == (REF, REF).
    """
    parts = rec.sample.split(":")
    gi = _fmt_index(rec.fmt, "GT")
    if gi < 0 or gi >= len(parts):
        return None, None, None
    gt = parts[gi].replace("|", "/")
    di = _fmt_index(rec.fmt, "DP")
    qi = _fmt_index(rec.fmt, "GQ")
    dp = _to_int(parts[di]) if 0 <= di < len(parts) else None
    gq = _to_float(parts[qi]) if 0 <= qi < len(parts) else None
    if gt in ("./.", ".", "./") or gt == "":
        return None, dp, gq
    if dp is not None and dp < min_dp:
        return None, dp, gq
    if gq is not None and gq < min_gq:
        return None, dp, gq
    alleles = rec.alleles()
    try:
        idx = [int(a) for a in gt.split("/") if a != "."]
    except ValueError:
        return None, dp, gq
    if not idx:
        return None, dp, gq
    try:
        bases = tuple(alleles[i] for i in idx)
    except IndexError:
        return None, dp, gq
    return bases, dp, gq


def _to_int(s):
    try:
        return int(s)
    except (ValueError, TypeError):
        return None


def _to_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


# ------------------------------------------------------------- divergence algebra


def alt_dosage(alleles, ref: str, alt: str):
    """Copies (0/1/2) of a specified `alt` base in a genotype `alleles` tuple."""
    if alleles is None:
        return None
    return sum(1 for a in alleles if a == alt)


def pairwise_mismatch(a_alleles, b_alleles) -> float:
    """Expected mismatch between one random allele from each diploid call (dxy site).

    a_alleles, b_alleles are tuples of base strings (len 1 or 2). Returns the mean
    over all allele pairs of 1[base_a != base_b] in [0,1], or None if either missing.
    """
    if a_alleles is None or b_alleles is None:
        return None
    n = 0
    mm = 0
    for x in a_alleles:
        for y in b_alleles:
            n += 1
            if x != y:
                mm += 1
    return mm / n if n else None


def mismatch_to_freq(a_alleles, ref: str, alt: str, p_alt: float):
    """Expected mismatch between a diploid call and a random allele from a panel.

    The panel carries base `alt` at frequency p_alt and base `ref` at 1-p_alt.
    Returns dxy(archaic, panel) at this site in [0,1], or None if the call is missing.
    """
    if a_alleles is None:
        return None
    tot = 0.0
    for x in a_alleles:
        # P(panel allele != x)
        if x == alt:
            tot += (1.0 - p_alt)
        elif x == ref:
            tot += p_alt
        else:  # x is a third allele: mismatches both ref and alt
            tot += 1.0
    return tot / len(a_alleles)


# --------------------------------------------------------------------- BED masks


def read_mask(path, chrom: str = None) -> np.ndarray:
    """Read a (gzipped) BED of callable intervals -> (N,2) int array of [start,end).

    BED is 0-based half-open. If `chrom` is given, keep only that chromosome
    (accepts '21' or 'chr21'). Files are typically single-chromosome already.
    """
    rows = []
    want = None
    if chrom is not None:
        want = {chrom, f"chr{chrom}", chrom.replace("chr", "")}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            if want is not None and f[0] not in want:
                continue
            rows.append((int(f[1]), int(f[2])))
    if not rows:
        return np.zeros((0, 2), dtype=np.int64)
    a = np.array(rows, dtype=np.int64)
    return a[np.argsort(a[:, 0])]


def intervals_intersect(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Intersect two sorted (N,2) half-open interval sets -> merged intersection."""
    if len(a) == 0 or len(b) == 0:
        return np.zeros((0, 2), dtype=np.int64)
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i, 0], b[j, 0])
        hi = min(a[i, 1], b[j, 1])
        if lo < hi:
            out.append((lo, hi))
        if a[i, 1] < b[j, 1]:
            i += 1
        else:
            j += 1
    return np.array(out, dtype=np.int64) if out else np.zeros((0, 2), dtype=np.int64)


def mask_membership(mask: np.ndarray, pos0) -> np.ndarray:
    """Vectorized: is each 0-based position inside any half-open mask interval?

    mask is (N,2) sorted by start. pos0 is an int array of 0-based coordinates.
    Returns a boolean array. (Convert 1-based VCF POS with pos0 = POS - 1.)
    """
    pos0 = np.asarray(pos0, dtype=np.int64)
    if len(mask) == 0:
        return np.zeros(len(pos0), dtype=bool)
    starts = mask[:, 0]
    ends = mask[:, 1]
    idx = np.searchsorted(starts, pos0, side="right") - 1
    ok = idx >= 0
    out = np.zeros(len(pos0), dtype=bool)
    j = idx[ok]
    out[ok] = (pos0[ok] >= starts[j]) & (pos0[ok] < ends[j])
    return out


def callable_bp_in_window(mask: np.ndarray, wstart: int, wend: int) -> int:
    """Total callable base pairs of `mask` (sorted intervals) within [wstart,wend)."""
    if len(mask) == 0:
        return 0
    lo = np.searchsorted(mask[:, 1], wstart, side="right")
    tot = 0
    for k in range(lo, len(mask)):
        s, e = int(mask[k, 0]), int(mask[k, 1])
        if s >= wend:
            break
        tot += min(e, wend) - max(s, wstart)
    return int(tot)


# --------------------------------------------------------------- ancestral FASTA


def read_ancestral_fasta(path, chrom: str) -> str:
    """Read one chromosome's ancestral sequence from an (optionally gz) FASTA.

    Ensembl EPO ancestral FASTAs use headers like '>ANCESTOR_for_chromosome:...:21:...'.
    Returns the sequence string (0-based indexable; uppercase). '-'/'.'/'N' mark
    positions with no confident ancestral call.
    """
    seq = []
    want = str(chrom).replace("chr", "")
    grabbing = False
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if grabbing:
                    break
                header = line[1:].strip()
                # match ':21:' or ' 21 ' or endswith the chromosome token
                toks = header.replace(":", " ").split()
                grabbing = want in toks
                continue
            if grabbing:
                seq.append(line.strip())
    return "".join(seq).upper()


# --------------------------------------------------------------------- self-test


def _selftest():
    # genotype parsing
    r_hom = Rec("21", 100, "G", [], "AC=0", "GT:DP:GQ", "0/0:10:40")
    al, dp, gq = sample_call(r_hom, min_dp=4, min_gq=30)
    assert al == ("G", "G") and dp == 10 and gq == 40.0, al

    r_het = Rec("21", 200, "G", ["A"], "AC=1", "GT:DP:GQ", "0/1:12:55")
    al, _, _ = sample_call(r_het)
    assert al == ("G", "A"), al
    assert alt_dosage(al, "G", "A") == 1

    r_low = Rec("21", 300, "G", ["A"], ".", "GT:DP:GQ", "1/1:2:9")
    assert sample_call(r_low, min_dp=4)[0] is None            # filtered by DP
    r_miss = Rec("21", 400, "G", [], ".", "GT:DP:GQ", "./.:0:0")
    assert sample_call(r_miss)[0] is None

    # divergence algebra
    assert pairwise_mismatch(("G", "G"), ("A", "A")) == 1.0
    assert pairwise_mismatch(("G", "G"), ("G", "G")) == 0.0
    assert pairwise_mismatch(("G", "A"), ("G", "A")) == 0.5
    # archaic hom-ALT vs a panel that is mostly REF -> high mismatch (0.75)
    assert abs(mismatch_to_freq(("A", "A"), "G", "A", 0.25) - 0.75) < 1e-12
    # archaic hom-REF vs panel with ALT freq 0.25 -> mismatch == 0.25
    assert abs(mismatch_to_freq(("G", "G"), "G", "A", 0.25) - 0.25) < 1e-12
    assert abs(mismatch_to_freq(("G", "A"), "G", "A", 0.5) - 0.5) < 1e-12

    # interval math
    a = np.array([[0, 100], [200, 300]], dtype=np.int64)
    b = np.array([[50, 250]], dtype=np.int64)
    inter = intervals_intersect(a, b)
    assert inter.tolist() == [[50, 100], [200, 250]], inter.tolist()
    assert callable_bp_in_window(a, 0, 1000) == 200
    assert callable_bp_in_window(a, 50, 250) == 100  # [50,100)+[200,250)
    mem = mask_membership(a, [0, 99, 100, 150, 200, 299, 300])
    assert mem.tolist() == [True, True, False, False, True, True, False], mem.tolist()

    # INFO parsing
    info = parse_info("AC=0;RM;Map20=0.75;CAnc=G")
    assert info["RM"] is True and info["Map20"] == "0.75" and info["CAnc"] == "G"

    print("vcflib self-test: OK")


if __name__ == "__main__":
    _selftest()
