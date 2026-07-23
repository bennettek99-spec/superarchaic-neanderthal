"""
extract_modern.py — compact the 1000 Genomes phase3 VCF to modern allele frequencies.

We never genotype the 2504 samples: 1000G phase3 INFO already exposes per-super-
population ALT frequencies (AFR_AF/EUR_AF/EAS_AF/AMR_AF/SAS_AF) and the EPO ancestral
allele (AA). We keep biallelic SNPs only and store, per chromosome:

  pos  ref  alt  afr  eur  eas  aa

Modern callability is treated as genome-wide (phase3 is high-coverage); a position
absent here is taken as monomorphic-reference in moderns (ALT freq 0) by scan_windows.

File resolution:
  chr22 is already present locally in the workspace ../vcf; chr21 downloads to
  data/modern. Override with --vcf.
"""
from __future__ import annotations

import argparse
import gzip
import time
from pathlib import Path

import vcflib as V

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
WORKSPACE_VCF = ROOT.parent / "vcf"     # ALL.chrN.phase3...vcf.gz live here (chr 1,2,8,12,19,22)
MODERN_DIR = ROOT / "data" / "modern"


def default_vcf(chrom: str) -> Path:
    local = MODERN_DIR / f"ALL.chr{chrom}.phase3.vcf.gz"
    if local.exists():
        return local
    hits = list(WORKSPACE_VCF.glob(f"ALL.chr{chrom}.phase3*vcf.gz"))
    return hits[0] if hits else local


def _aa_base(aa: str) -> str:
    # 1000G AA field looks like "G|||" or ".|||"; first token, uppercased.
    if not aa or aa == ".":
        return "."
    b = aa.split("|", 1)[0].upper()
    return b if b in ("A", "C", "G", "T") else "."


def extract(chrom: str, vcf: Path = None, force: bool = False) -> Path:
    vcf = Path(vcf) if vcf else default_vcf(chrom)
    CACHE.mkdir(parents=True, exist_ok=True)
    out = CACHE / f"1000G.chr{chrom}.freq.tsv.gz"
    if out.exists() and not force:
        print(f"  SKIP {out.name} (exists)")
        return out
    if not vcf.exists():
        print(f"  MISSING modern VCF for chr{chrom}: {vcf}")
        return out
    t0 = time.time()
    n_in = n_out = 0
    tmp = out.with_suffix(".tmp.gz")
    with gzip.open(tmp, "wt", encoding="ascii") as w:
        w.write("pos\tref\talt\tafr\teur\teas\taa\n")
        with V.open_text(vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                n_in += 1
                ref, alt = f[3], f[4]
                if len(ref) != 1 or len(alt) != 1 or alt == ".":
                    continue  # biallelic SNPs only
                i = V.parse_info(f[7])
                if i.get("VT") not in (None, "SNP") and "," in i.get("VT", ""):
                    continue
                w.write(f"{f[1]}\t{ref}\t{alt}\t{i.get('AFR_AF','.')}\t"
                        f"{i.get('EUR_AF','.')}\t{i.get('EAS_AF','.')}\t{_aa_base(i.get('AA','.'))}\n")
                n_out += 1
                if n_in % 1_000_000 == 0:
                    print(f"    chr{chrom}: {n_in:,} lines ({time.time()-t0:.0f}s)", flush=True)
    tmp.replace(out)
    print(f"  chr{chrom}: {n_out:,} biallelic SNPs of {n_in:,} lines in "
          f"{time.time()-t0:.0f}s -> {out.name}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=["21", "22"])
    ap.add_argument("--vcf", default=None, help="explicit VCF path (single chrom)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    for c in args.chroms:
        extract(c, vcf=args.vcf, force=args.force)


if __name__ == "__main__":
    main()
