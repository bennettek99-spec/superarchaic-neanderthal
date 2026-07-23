# Dataset Inventory (Deliverable 2)

All coordinates are **hg19 / GRCh37 / 1000g**, consistent across every source (verified by
inspection). Large files are git-ignored; this document is the manifest.

## High-coverage archaic genomes — Max Planck EVA (`https://cdna.eva.mpg.de/neandertal/`)

| Genome | Role | Coverage | Reference | Path (EVA) | chr21+22 pilot |
|---|---|---|---|---|---|
| Altai Denisovan (`DenisovaPinky`) | Denisovan | ~30× | Meyer 2012; Prüfer 2014 | `altai/Denisovan/DenisovaPinky.hg19_1000g.{c}.mod.vcf.gz` | 0.77 + 0.74 GB |
| Altai Neanderthal | Neanderthal | ~52× | Prüfer 2014 | `altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.{c}.mod.vcf.gz` | 0.90 + 0.86 GB |
| Vindija 33.19 | Neanderthal | ~30× | Prüfer 2017 | `Vindija/VCF/Vindija33.19/chr{c}_mq25_mapab100.vcf.gz` | 0.62 + 0.58 GB |
| Chagyrskaya 8 | Neanderthal | ~27× | Mafessoni 2020 | `Chagyrskaya/VCF/chr{c}.noRB.vcf.gz` | 0.58 + 0.53 GB |

These are **all-sites VCFs** (invariant positions included). The Altai/Denisovan 2014 `.mod`
files carry rich Ensembl-EPO INFO annotation at aligned sites — `CAnc` (chimp–human ancestor =
polarizing ancestral allele), `Map20` (Duke 20-mer mappability), `RM` (repeat-masked), `UR`
(copy-number/"unique region" flag — present on most sites, **not** an exclusion), `CpG`,
`bSC`/`mSC`/`pSC`/`GRP` (background-selection/conservation), and embedded 1000G AFR/EUR/ASN
frequencies. Vindija/Chagyrskaya callsets lack this INFO (position-level annotation is taken
from one source only).

Genome-wide sizes (for planning a scaled run): Altai ≈ 70 GB, Denisovan ≈ 63 GB,
Vindija ≈ 48 GB, Chagyrskaya ≈ 48 GB → **≈ 230–260 GB** for all four.

## Callability masks (FilterBed, mq25 / mapab100)

Per-genome BED of **callable** intervals (0-based). `Vindija/FilterBed/{Altai,Denisova,
Vindija33.19}/chr{c}_mask.bed.gz` and `Chagyrskaya/FilterBed/chr{c}_mask.bed.gz`. Altai chr21
callable ≈ 21.3 Mb (~46% of assembled chr21). The pipeline intersects all four into a **common
mask**; every cross-genome divergence is computed only there, so archaic comparisons use
identical sites (essential — the 2014 vs 2016-18 masks differ).

## Modern humans — 1000 Genomes phase 3

`ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` (EBI). chr22
already present locally in the workspace `../vcf`; chr21 fetched to `data/modern/`. Super-
population ALT frequencies (`AFR_AF`, `EUR_AF`, `EAS_AF`) and the EPO ancestral allele (`AA`)
are read from INFO — the 2504 samples are never genotyped individually. African panel = AFR
super-population; European = EUR.

## Polarization

`CAnc` (chimp–human ancestor) from the archaic INFO where present; genome-wide fallback =
Ensembl EPO 6-primate ancestral, GRCh37 (`homo_sapiens_ancestor_GRCh37_e71.tar.bz2`, ~0.8 GB,
in `data/ancestral/`). Pairwise divergence (the primary CDDR statistic) does not require
polarization; derived/private-allele metrics do.

## Local caches produced by the pipeline (git-ignored)

- `data/cache/{genome}.chr{c}.variants.tsv.gz` — resolved variant-site genotypes + annotation
  (from `extract_variants.py`; ~50 s per genome per small chromosome).
- `data/cache/1000G.chr{c}.freq.tsv.gz` — modern super-population frequencies (`extract_modern.py`).

## Pilot scope acquired (this run)
chr21 + chr22 for all four archaics + masks + moderns + ancestral ≈ **7 GB**, downloaded in
~20 min at ~6–8 MB/s from EVA. Everything is resumable/idempotent (`src/download_pilot.sh`).
