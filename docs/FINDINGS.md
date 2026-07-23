# Findings — Superarchaic Introgression in Neanderthals (chr21+chr22 pilot)

**Status: pilot (chromosomes 21 + 22 ≈ 3% of the genome).** The pipeline, statistics, and
simulation calibration are complete and validated; the genome-wide scan is the next step and is
gated on an explicit data-download decision (≈230–260 GB). Every statement below is a *pilot*
result. Language follows the project's interpretation rules ("candidate", "consistent with",
"cannot distinguish"); nothing is called "superarchaic".

---

## Bottom line (answer to the success-criterion question)

> *"Do Neanderthals show reproducible evidence of sharing the same candidate deeply divergent
> genomic regions observed in Denisovans, beyond what is expected under ordinary
> Neanderthal–Denisovan shared ancestry?"*

**On the pilot data: No.** Three independent lines converge:

1. **Focal-genome deep-divergence rates are equal across all four archaics** (Denisovan 3.5%,
   Altai 3.3%, Vindija 3.5%, Chagyrskaya 3.5% of 20-kb windows with z≥3). Equality is the
   signature of a *shared* confounder (per-window mutation-rate / coalescent-variance
   heterogeneity that lifts every archaic-vs-African divergence together), **not** lineage-
   specific deep ancestry.
2. **Candidate Deep Divergence Regions are Denisovan-specific, not Neanderthal-shared.** Using
   the confound-controlled contrast (below), 20 candidates survive artifact + multi-scale
   filtering on chr21+22 — all in the **Model-B direction (Denisovan-specific)**. The
   Neanderthal-comparison classifier put 8/10 (20 kb) and 2/2 (50 kb) CDDRs in **Category 3
   (Denisovan-only)**; essentially none show the shared-with-all-Neanderthals pattern.
3. **The strongest apparent Denisovan-specific windows are ≥65% artifacts** (repeat /
   low-mappability / segmental-duplication), i.e. paralog-mismapping inflating one genome's
   divergence — not ancestry.

So on chr21+22 there is **no Model-A (Neandersovan / shared-superarchaic) signal**, and the
residual Denisovan-specific candidates are (a) the *expected* Model-B direction (consistent with
Hubisz 2020 / Fu 2026 / DEEP 2026) and (b) not yet distinguishable from a heavy-tailed null on
3% of the genome. **A careful null is the honest pilot outcome.**

---

## 1. The pipeline is validated (correct population genetics)

On the common callable mask, median pairwise divergence per bp (chr21+22, 50 kb windows) is
textbook-ordered:

| Comparison | median divergence/bp | expectation |
|---|---|---|
| archaic vs African (den/alt/vin/chag) | 0.00145–0.00151, **all ≈ equal** | all archaics ~equidistant from Africans ✓ |
| Denisovan vs each Neanderthal | 0.00112–0.00115 | Nea–Den < archaic–modern ✓ |
| Neanderthal vs Neanderthal | 0.00018–0.00023 | within-clade smallest ✓ |

The archaic-vs-African values converge only **after** switching to a common intersection mask;
per-genome masks (2014 Altai/Denisovan vs 2016–18 Vindija/Chagyrskaya) gave a spurious
Vindija/Chagyrskaya deficit (0.0016 vs 0.0023) — a methodological trap this pipeline now avoids.

## 2. Raw divergence is confounded; the contrast is not

Because per-window rate variation inflates all archaic divergences together, a deep-divergence
window is **not** per se lineage-specific. The unconfounded statistic for Denisovan-specific
depth (the Model-B signature) is the contrast

    den_excess = div_den_afr − mean(div_alt_afr, div_vin_afr, div_chag_afr)

which cancels the shared component. **Simulation validates it**: across msprime Models 0–4
(20–30 Mb, seeded), the mean `den_excess` is largest under **M2 (Denisovan-only)** and smallest
under **M0 (null)** — the contrast tracks genuine Denisovan-specific ancestry.

## 3. Simulation calibration (Models 0–4) — and a key power limit

| Model | CDDR rate (z≥3 on den_afr) | %CDDRs classed "shared w/ Neanderthals" |
|---|---|---|
| M0 none (null) | 0.007 | 8% |
| M1 Neandersovan (**Model A**) | 0.020 | **6%** |
| M2 Denisovan-only (**Model B**) | 0.022 | 3% |
| M3 separate pulses | 0.027 | 4% |
| M4 ancient structure | 0.036 | 25% |

Two calibration lessons that shape the whole design:
- **The co-located-sharing test has low power for Model A.** Even when Model A is *true* (M1),
  only 6% of Denisovan CDDRs are classed "shared" — because superarchaic ancestry deposited in
  the Neandersovan ancestor is retained at *different loci* in Neanderthals vs Denisovans after
  ~400 ky of drift/recombination. So *co-location* is the wrong expectation for Model A; the
  right test is whether Neanderthals show their own genome-wide deep-divergence **excess** (they
  do not, beyond the shared confound — see Bottom line #1).
- **Ancient structure (M4) is the strongest confounder** (highest CDDR rate *and* highest
  apparent "sharing"), reproducing the long-standing Rogers-vs-structure ambiguity. Any positive
  claim must be defended against it.
- The msprime null (constant μ, constant recombination) **underestimates real overdispersion**
  (real z≥3 rate ~3.5% vs sim ~0.7%), so the sim gives *qualitative* calibration only; a
  genome-wide **empirical** null + block-jackknife is required for real significance.

## 4. Candidate Deep Divergence Region (CDDR) catalog — Denisovan-specific

`results/cddr/candidate_catalog.win50000.tsv` (57 windows, z_den_excess ≥ 3 on chr21+22):

| Evidence tier | n | meaning |
|---|---|---|
| Artifact-likely | 37 (65%) | repeat / low-mappability / segdup / low-callability dominated |
| Moderate | 20 | passes artifact filter **and** replicates across ≥1 other window size |
| Weak / Insufficient / Strong | 0 | — |

The 20 "Moderate" candidates are **Denisovan-specifically deep** and survive the alternative-
explanation filter — they are legitimate *Model-B-direction* candidates. **None are
Neanderthal-shared** (no Model-A candidates). They are **not** significance-tested: on 3% of the
genome with a heavy-tailed null and no block-jackknife, they remain candidates, not findings.

## 5. Alternative-explanation assessment (per the spec)

| Alternative | How addressed in the pilot | Verdict on pilot candidates |
|---|---|---|
| Mapping / paralog / segmental duplication | `Map20`, `RM`, `UR` covariates; ≥65% of top windows flagged | **dominant** explanation of the extreme tail |
| Low mappability / callability | mq25/mapab100 common mask; `callable_frac` gate | controlled at site level |
| Mutation-rate variation | the `den_excess` contrast cancels shared per-window rate | this is *why* the contrast is used |
| Ancient population structure | M4 simulation shows it mimics/inflates signal | **unresolved** — needs genome-wide + ARG |
| Incomplete lineage sorting | Nea–Nea and Den–Nea baselines; contrast design | baseline, not excess |
| aDNA damage | transversion-only pass available (`config.transversions_only_variant`) | pending in scaled run |
| Reference bias / low complexity / CpG | `CpG`, conservation covariates recorded | recorded, not yet decisive |
| Contamination | high-coverage archaics; not re-derived here | assumed handled upstream |

No candidate region "survives" as superarchaic: the extreme tail is artifact-dominated, and the
cleaner residual is (a) Denisovan-specific and (b) not distinguishable from ancient structure or
a heavy-tailed null on this data slice.

## 6. Final evidence ranking

- **Model A (superarchaic → Neanderthals too): no evidence** on the pilot. Not detected by any
  of the three lines above; simulation shows the co-location test is under-powered for it, so
  "not detected" is weaker than "excluded" — the genome-wide focal-rate test is the decisive one.
- **Model B (Denisovan-specific): weak, expected, unconfirmed candidates.** 20 Denisovan-specific
  regions survive artifact filtering, in the direction the literature predicts, but not
  significance-tested and not yet separated from ancient structure.
- **Dominant real signal on chr21+22: shared per-window rate variation + mapping artifacts** —
  i.e. the alternative explanations, not archaic introgression.

## 7. What the genome-wide analysis must add (to actually answer the question)

1. Scale to all autosomes (extraction ≈ 2 h; scan hours) → power + a real empirical null.
2. **Block jackknife / leave-one-chromosome-out** significance on `den_excess` and on the
   focal-genome CDDR-rate contrast (the decisive Model-A test).
3. Polarized **private-derived-allele density** and local-genealogy / TMRCA at surviving
   candidates (needs the ancestral FASTA, now downloaded) to separate deep-ancestry from rate.
4. Explicit **ancient-structure (M4) model fitting** — the one confounder the pilot cannot rule
   out — ideally with an ARG method (ARGweaver-D / SINGER) at candidate loci.
5. Genome-wide segmental-duplication / RepeatMasker tracks for exact per-window artifact masking.

**Reproduce:** `src/download_pilot.sh` → `extract_variants.py` → `extract_modern.py` →
`scan_windows.py` → `neanderthal_compare.py` → `rank_candidates.py`; `sims/superarchaic_sim.py
--mode discriminate` for calibration. Tests: `python tests/test_pipeline.py` (8 pass).
