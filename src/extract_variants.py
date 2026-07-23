"""
extract_variants.py — one-pass compaction of each archaic all-sites VCF.

Streams every archaic VCF once, keeps only sites carrying an ALT allele, and writes
a compact gzipped TSV of *resolved* genotype bases plus a few QC/annotation columns.
Downstream (scan_windows.py) then works from these small caches + the BED masks,
never re-touching the multi-GB VCFs.

Why this is safe for divergence:
  At a site that is variant in genome X but absent from genome Y's cache, Y is either
  homozygous-reference (if the position is in Y's callability mask) or not-callable.
  The window scan resolves that from the masks, so caching variant sites only loses
  no information about pairwise differences.

Output columns (per genome, per chromosome):
  pos ref alt a1 a2 dp gq canc map20 rm ur cpg bsc afr eur asn
where a1/a2 are the called allele BASES ('.' if missing/filtered). The remaining
columns come from the Ensembl-EPO INFO annotation carried by the Altai/Denisovan
'.mod' VCFs (canc=chimp-human ancestor base, map20=Duke mappability, rm=repeat-mask
flag, ur=copy-number/segdup flag, cpg=CpG flag, bsc=B-score, afr/eur/asn=1000G
super-population ALT freqs). They are '.' for Vindija/Chagyrskaya, which lack them;
because these are position-level genomic annotations they only need one source.

Usage:
  python src/extract_variants.py --chroms 21 22
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import vcflib as V

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
CACHE = ROOT / "data" / "cache"

# genome -> (filename template using {c}, min_dp, min_gq). Thresholds are gentle;
# the masks do the heavy callability lifting. Chagyrskaya/Vindija GQ may be absent
# -> sample_call simply ignores a missing FORMAT key.
GENOMES = {
    "Denisova":      ("Denisova.{c}.vcf.gz",      5, 30),
    "AltaiNea":      ("AltaiNea.{c}.vcf.gz",       5, 30),
    "Vindija33.19":  ("Vindija33.19.{c}.vcf.gz",   3, 0),
    "Chagyrskaya":   ("Chagyrskaya.{c}.vcf.gz",    3, 0),
}


def extract(genome: str, chrom: str, force: bool = False) -> Path:
    tmpl, min_dp, min_gq = GENOMES[genome]
    vcf = RAW / tmpl.format(c=chrom)
    CACHE.mkdir(parents=True, exist_ok=True)
    out = CACHE / f"{genome}.chr{chrom}.variants.tsv.gz"
    if out.exists() and not force:
        print(f"  SKIP {out.name} (exists)")
        return out
    if not vcf.exists():
        print(f"  MISSING {vcf.name} — not downloaded yet", file=sys.stderr)
        return out
    t0 = time.time()
    n_in = n_out = 0
    tmp = out.with_suffix(".tmp.gz")
    cols = ["pos", "ref", "alt", "a1", "a2", "dp", "gq",
            "canc", "map20", "rm", "ur", "cpg", "bsc", "afr", "eur", "asn"]
    with gzip.open(tmp, "wt", encoding="ascii") as w:
        w.write("\t".join(cols) + "\n")
        for rec in V.iter_records(vcf, keep_only_variant=True):
            n_in += 1
            alleles, dp, gq = V.sample_call(rec, min_dp=min_dp, min_gq=min_gq)
            if alleles is None:
                a1 = a2 = "."
            elif len(alleles) == 1:
                a1 = a2 = alleles[0]
            else:
                a1, a2 = alleles[0], alleles[1]
            i = V.parse_info(rec.info)
            row = [
                str(rec.pos), rec.ref, ",".join(rec.alt), a1, a2,
                "" if dp is None else str(dp), "" if gq is None else str(gq),
                i.get("CAnc", "."), i.get("Map20", "."),
                "1" if i.get("RM") else "0", "1" if i.get("UR") else "0",
                "1" if i.get("CpG") else "0", i.get("bSC", "."),
                i.get("AFR_AF", "."), i.get("EUR_AF", "."), i.get("ASN_AF", "."),
            ]
            w.write("\t".join(row) + "\n")
            n_out += 1
            if n_in % 500_000 == 0:
                print(f"    {genome} chr{chrom}: {n_in:,} variant lines "
                      f"({time.time()-t0:.0f}s)", flush=True)
    tmp.replace(out)
    print(f"  {genome} chr{chrom}: {n_out:,} sites kept of {n_in:,} variant lines "
          f"in {time.time()-t0:.0f}s -> {out.name}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=["21", "22"])
    ap.add_argument("--genomes", nargs="+", default=list(GENOMES))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    for c in args.chroms:
        print(f"chr{c}:")
        for g in args.genomes:
            extract(g, c, force=args.force)


if __name__ == "__main__":
    main()
