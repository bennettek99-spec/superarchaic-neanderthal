# Superarchaic Introgression in Neanderthals

Hypothesis-testing project: **do high-coverage Neanderthal genomes contain regions consistent with
inheritance from the same deeply diverged ("superarchaic") ghost lineage that contributed ancestry
to Denisovans** — beyond what ordinary Neanderthal–Denisovan shared ancestry and incomplete lineage
sorting predict?

- **Model A (Neandersovan):** superarchaic → common ancestor of Nea + Den ⇒ both carry a signal.
- **Model B (Denisovan-only):** superarchaic → Denisovans after the Neanderthal split ⇒ Den-only.

**The objective is not to prove Model A.** A rigorous null / "cannot distinguish" result is a
success. See [`docs/LITERATURE_REVIEW.md`](docs/LITERATURE_REVIEW.md) for the full survey and the
predictions that separate A from B — the literature is genuinely split (Rogers 2020 ⇒ A; Hubisz
2020, TRACE 2026, DEEP 2026 ⇒ B-leaning).

## Design (whole-genome, high-coverage only — NOT sparse SNP-array data)

1. **Window scan** of the Denisovan genome (20/50/100 kb) → per-window divergence, private-allele
   density, coalescence-depth proxy, plus quality/mappability/callability/repeat/coverage metrics.
2. **Flag Candidate Deep Divergence Regions (CDDRs)** — conservatively; never auto-labelled "superarchaic."
3. **Neanderthal comparison** at every CDDR across **Altai + Vindija + Chagyrskaya**; classify
   Category 1–4 (shared-all / shared-some / Den-only / inconclusive). Replication across ≥2
   Neanderthals is mandatory.
4. **Local genealogies** and topology vs the species tree; deep-branch / ILS assessment.
5. **Alternative explanations** (ILS, ancient structure, mutation-rate, balancing selection, low
   recombination, reference bias, aDNA damage, mapping/seg-dup/low-complexity, error, contamination).
6. **Simulations** (Models 0–4 in msprime) through the identical pipeline → FPR, power, thresholds.
7. **Sensitivity** — block jackknife, bootstrap, leave-one-window/chromosome-out, multiple masks /
   window sizes / reference panels; robustness to injected missingness / error / damage / coverage loss.
8. **Final evidence ranking** per region: Strong / Moderate / Weak / Artifact-likely / Insufficient.

## Interpretation rules
Language is restricted to "consistent with", "compatible with", "cannot distinguish", "candidate",
"deeply divergent ancestry", "ghost lineage", "requires additional evidence." Observations,
statistics, demographic interpretation, and speculation are kept separate. No admixture % / date is
estimated unless the data genuinely support it.

## Data (see `docs/DATA.md` once populated; large files are git-ignored)
High-coverage archaics from the Max Planck EVA server (`cdna.eva.mpg.de/neandertal`, hg19/1000g
coordinates): Altai Neanderthal, Vindija 33.19, Chagyrskaya 8, Altai Denisovan — plus FilterBed
callability masks (mq25/mapab100). Moderns: 1000 Genomes phase3 (local) and/or SGDP. Polarization:
chimp allele / Ensembl EPO ancestral (GRCh37). **All-sites VCFs are large (~50–70 GB per genome
genome-wide); acquisition is staged and gated on explicit approval.**

## Environment
Windows 11, Python 3.14 (numpy/pandas/scipy/msprime/tskit/matplotlib/sklearn). **No
bcftools/cyvcf2/pysam/scikit-allel** → a pure-Python streaming VCF reader is used. Reuses the
validated statistics package from the main pipeline (`../archaic-introgression/archaic`,
f4/D/f4-ratio + block jackknife). Honest-framing ethos and conventions follow that pipeline.

## Status — chr21+chr22 pilot complete
- [x] **1. Literature review + comparison table** — `docs/LITERATURE_REVIEW.md`
- [x] **2. Dataset inventory** — `docs/DATA.md` (pilot ≈7 GB acquired; genome-wide ≈230–260 GB gated)
- [x] **3. CDDR catalog** — `results/cddr/candidate_catalog.win50000.tsv` (Denisovan-specific contrast)
- [x] **4. Neanderthal comparison** — `results/cddr/neanderthal_classified.win*.tsv` (mostly Category 3)
- [x] **5. Simulation report (Models 0–4)** — `results/sims/` + `sims/superarchaic_sim.py`
- [x] **6. Alternative explanations + sensitivity** — artifact covariates, common mask, 3 window sizes
- [x] **7–9. Findings + evidence ranking** — `docs/FINDINGS.md`, `results/figures/pilot_summary.png`

**Pilot answer to the success criterion: No Model-A (shared-superarchaic) signal on chr21+22.**
Deep-divergence tail is dominated by shared per-window rate variation + mapping artifacts; the
20 clean residual candidates are all Denisovan-specific (Model-B direction), not Neanderthal-
shared, and not yet significance-tested. See `docs/FINDINGS.md`. Genome-wide scan is the next step.

## Pipeline order
`src/download_pilot.sh` → `src/extract_variants.py` → `src/extract_modern.py` →
`src/scan_windows.py` → `src/neanderthal_compare.py` → `src/rank_candidates.py` →
`src/make_figures.py`; calibration: `sims/superarchaic_sim.py --mode discriminate`.
Tests: `python tests/test_pipeline.py`.

## Laptop-friendly staggered acquisition (`src/stagger.sh`)
Genome-wide raw VCFs total ~230–260 GB — too much to hold at once on a laptop. The scan only
needs the small per-chromosome **variant caches** (~1–2 MB/genome), not the VCFs, so acquisition
is staged: for each chromosome, one genome at a time, **download → extract cache → verify →
delete the VCF**. Peak disk stays at a single VCF (~6 GB worst case). Resumable (a genome is
skipped once its cache exists; a VCF is deleted only after its cache is confirmed).

```bash
bash src/stagger.sh 20 19        # process a couple of chromosomes now
bash src/stagger.sh $(seq 1 22)  # whole autosome set, deleting VCFs as it goes
```
After caches accumulate, run `scan_windows.py --chroms <list>` etc. on the cached chromosomes —
no VCFs required. Stop/resume anytime; completed chromosomes are skipped.
