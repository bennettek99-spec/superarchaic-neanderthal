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

## Why co-location is the wrong test (and what replaced it)

Under Model A, superarchaic ancestry entered the **common** Neanderthal+Denisovan ancestor and
then drifted and recombined independently for ~420 ky, so the two lineages end up carrying it at
**different loci**. Asking "are Neanderthals deep at *Denisovan* deep windows?" therefore has ~6%
power against a true Model A. Worse, the pilot's primary statistic

    den_excess = div_den_afr − mean(div_alt_afr, div_vin_afr, div_chag_afr)

**cancels Model A exactly** — both terms rise together — so Model A was structurally invisible to
it. See [`docs/FINDINGS.md §0`](docs/FINDINGS.md) for the correction, including a withdrawn
headline p-value that turned out to be a selection artefact.

The redesign tests Model A **without requiring the lineages to share loci**:

* an **external mutation-rate denominator** (`src/ratemap.py`) — human-vs-EPO-ancestral
  substitution density per window, which removes rate variation *without* contrasting the
  archaics against each other, so it leaves Model A intact;
* **polarized site patterns** (`src/scan_windows.py`) — in particular `pat_nea_all`, derived in
  all three Neanderthals and ancestral in the Denisovan, whose **spatial clustering** (not its
  mean) distinguishes introgressed haplotype blocks from incomplete lineage sorting;
* **locus-free statistics** (`src/modelA.py`) with weighted block-jackknife errors, run through
  the **same code path** on real and simulated data so power is calibrated rather than assumed.

## Design (whole-genome, high-coverage only — NOT sparse SNP-array data)

1. **Window scan** (20/50/100 kb) → per-window divergence, polarized derived-allele patterns,
   external mutation-rate denominator, plus quality/mappability/callability/repeat/coverage metrics.
2. **Model-A battery** — T1 clade clustering (the powered test), T2 tail symmetry (Model B),
   T3 cross-lineage covariance; block jackknife over 5 Mb blocks throughout.
3. **Flag Candidate Deep Divergence Regions (CDDRs)** — conservatively; never auto-labelled "superarchaic."
4. **Neanderthal comparison** at every CDDR across **Altai + Vindija + Chagyrskaya**; classify
   Category 1–4. Retained as a Model-B probe, with its ~6% Model-A power stated wherever it is reported.
5. **Alternative explanations** (ILS, ancient structure, mutation-rate, balancing selection, low
   recombination, reference bias, aDNA damage, mapping/seg-dup/low-complexity, error, contamination).
6. **Simulations** (Models 0–4 in msprime) through the identical statistics module, with per-window
   mutation rates resampled from the *measured* rate map and the African panel matched to 1000G
   → FPR, power, thresholds (`sims/power_calibration.py`).
7. **Sensitivity** — quality filters, choice of drift control, block jackknife, leave-one-chromosome-out,
   multiple masks / window sizes.
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

## Status — chr6 + chr13–22 (964 Mb, 33% of the autosomes, 15,616 usable 50-kb windows)
- [x] **1. Literature review + comparison table** — `docs/LITERATURE_REVIEW.md`
- [x] **2. Dataset inventory** — `docs/DATA.md`; chr13–22 cached, chr6 in acquisition
- [x] **3. External mutation-rate map** — `results/ratemap/`, median 0.0058–0.0068/bp
- [x] **4. Polarized site patterns** — in the window tables
- [x] **5. Model-A battery + sensitivity** — `results/modelA/`
- [x] **6. Simulation calibration with realistic overdispersion** — `results/sims/battery_summary.tsv`
- [x] **7. Selection-bias correction** — `docs/FINDINGS.md §0`
- [ ] **8. Replicate-level power/FPR** — `sims/power_calibration.py`, ready to run (~3 h)
- [ ] **9. Demography fitted to the observed pattern spectrum**

**Current answer to the success criterion.** *Co-location:* no — and it was the wrong question
(~6% power by construction). *Neanderthal-clade deep-ancestry clustering:* a **weak positive
signal**, α(`pat_nea_all`) − α(`pat_alt_only`) = **+0.114 ± 0.046 (z = 2.49)**, which simulation
puts at ≈0 under the null and under Model B and at +0.112 under Model A. It is stable across
data-quality filters — but **ancient structure (M4) reproduces it** (+0.152), it depends on which
drift control is used, and the demography is unfitted. **A candidate signal, not a finding**, and
explicitly not a claim about superarchaic ancestry. See [`docs/FINDINGS.md`](docs/FINDINGS.md).

**Updated with chr6: the signal does not hold.** The contrast falls to **+0.067 ± 0.061
(z = 1.09)** — adding 26% more data halved it, and leave-one-chromosome-out shows chr13–22 was
the outlier (dropping any other chromosome gives +0.03…+0.09). The value is consistent with the
null *and* with Model A, which are only 0.129 apart. **No significant signal, no claim in either
direction.**

**What changed is the power, and that is the point.** The old co-location test had ~6% power
under a true Model A *at any genome coverage* — more sequence could never have fixed it. The
replacement has **68% power now and 98% with all autosomes**, so finishing the genome is now the
highest-value next step, inverting the pilot's advice.

| coverage | windows | SE | z under Model A | power |
|---|---|---|---|---|
| now (chr6 + chr13–22, 33%) | 15,616 | 0.061 | 2.11 | 68% |
| all autosomes | 46,669 | 0.035 | 3.66 | 98% |

## Pipeline order
`src/stagger.sh <chroms>` (download → extract → verify → delete, now including the rate map) →
`src/scan_windows.py --jobs 4` → `src/modelA.py --sensitivity` → `src/rank_candidates.py` /
`src/neanderthal_compare.py` (Model-B probes) → `src/make_figures.py`.
Calibration: `sims/superarchaic_sim.py` battery, then `sims/power_calibration.py`.
Tests: `python tests/test_pipeline.py` (16 pass).

## Performance
The scan is vectorized end to end and parallel across chromosomes: **10 chromosomes × 3 window
sizes in ~90 s** on 4 workers. Per-window callability came from a Python loop over mask intervals
and is now an O(n) prefix-sum (~1200× faster on a real chromosome mask, verified against the old
path); base coding is a uint8 lookup instead of a per-element `.map()`; the four pandas outer
merges are a single sorted-union + `searchsorted` scatter. That last change also fixed a real bug —
the caches contain duplicate positions (multi-allelic records split across VCF lines) which the
outer merges **double-counted** in the annotation means.

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
no VCFs required. Stop/resume anytime; completed chromosomes are skipped. Each chromosome also
builds its mutation-rate map: the hg19 FASTA (~50–90 MB) is downloaded, reduced to per-window
substitution counts, and deleted, the same way the VCFs are.

## Overnight calibration (launch when ready)
```bash
python sims/power_calibration.py --reps 20 --jobs 6 --observed results/modelA/summary.chr13,14,15,16,17,18,19,20,21,22.win50000.tsv
```
~3 h. Produces the replicate-level false-positive rate, the power of each statistic under
Models 1–4, and where the observed value falls in each model's distribution — the calibrated
statement that a single-pass block jackknife cannot give.
