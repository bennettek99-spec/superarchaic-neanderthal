# Findings — Superarchaic Introgression in Neanderthals

**Scope: chr6 + chr13–22 — 964 Mb, 33% of the autosomes, 15,616 usable 50-kb windows.**
Language follows the project's interpretation rules; nothing is called "superarchaic".

**Headline.** The pilot's reported Model-A null was a selection artefact and is
withdrawn (§0). Its replacement — a locus-free test that actually has power against
Model A — gives **+0.067 ± 0.061 (z = 1.09): no significant signal**, a value consistent
with the null *and* with Model A, because at 33% of the autosomes the two are only ~2.1
standard errors apart. The project's real output this round is therefore a **method with
quantified power (68% now, 98% at whole-genome) where the previous one had ~6% at any
coverage** — an unanswerable question turned into a sample-size one (§5).

---

## 0. Correction to the previous headline result

The earlier version of this document reported, as its headline:

> permutation p (Model-A co-location) = **1.000** (S_obs −3.35, null +0.10±0.14)

**That result was an artefact of the test's own selection rule and is withdrawn.**

`significance.py` selected candidate windows on

    z_den_excess ≥ 3,  where  den_excess = div_den_afr − mean(div_alt_afr, div_vin_afr, div_chag_afr)

and then used, as the test statistic, the mean Neanderthal residual depth at those same
windows. Selecting windows for having high Denisovan **and low Neanderthal** divergence
and then reporting that Neanderthal divergence is low there cannot come out any other
way. S_obs = −3.35 against a null s.d. of 0.14 is ~25 σ — a statement about the
selection rule, not about Neanderthals. The p = 1.000 was never evidence against Model A.

There was a second, deeper problem. The whole pilot was built on `den_excess`, which
**cancels Model A exactly**: under Model A superarchaic ancestry entered the common
Neanderthal+Denisovan ancestor, so both terms rise together and the contrast returns
zero. Model A was structurally invisible to the pilot's primary statistic. The reported
~6% power of the co-location test was a symptom of this, not the cause.

Both are fixed below. `src/significance.py` now selects on Denisovan depth alone and is
retained only as a Model-B probe; `--show-legacy` reproduces the old artefact on demand.

### The co-location test cannot be rescued at all

De-biasing it does not produce a usable test — it produces a smaller artefact. With
candidates selected on Denisovan depth alone and a residual taken in fine quantile
bins, S_obs = **+0.415, p = 0.004** — the *opposite sign* to the withdrawn result. But a
z ≥ 3 threshold cuts **through** a residual bin, so the candidates are the upper part of
the bins they fall in and inherit a positive residual by construction. Selecting the
same number of windows at **whole-bin boundaries** collapses S to **+0.03 … +0.10**:

| candidate set | S |
|---|---|
| z ≥ 3 threshold (126 windows) | **+0.415** |
| top 1 whole bin (50 windows) | +0.032 |
| top 2 whole bins (100 windows) | +0.037 |
| top 3 whole bins (150 windows) | +0.100 |

So the p = 0.004 is still mostly selection. And the problem is structural, not a matter
of finer binning: **condition fully on Denisovan depth and the test is vacuous by
construction** (the residual mean inside a stratum is zero by definition); **condition
less and it is dominated by the ~0.73 Denisovan–Neanderthal correlation that shared
ancestry produces under every model, including the null.** Its permutation null is
simply the wrong null.

This is why the redesign does not try to fix the co-location test. It replaces it.

---

## 1. What changed methodologically

### An external mutation-rate denominator (`src/ratemap.py`)
Per-window mutation-rate heterogeneity inflates every archaic-vs-modern divergence
together. The pilot removed it by contrasting the archaics against each other — which
also removed Model A. The replacement is a denominator that involves no archaic genome
at all: the density of substitutions on the human lineage since the human–chimp
ancestor, i.e. positions where hg19 differs from the Ensembl EPO ancestral sequence,
inside the common callable mask.

Measured median `sub_rate` is **0.0058–0.0068 per bp** across chr13–22 — the expected
value for ~6 Myr of human-lineage evolution, with ancestral coverage of 94–99% of the
callable mask. At 50 kb that is ~300 substitutions per window (~6% Poisson noise)
against ~72 for `div_den_afr`, so it is a near-noiseless local rate estimate.

It correlates **+0.40, +0.38, +0.37, +0.38** with Denisovan, Altai, Vindija and
Chagyrskaya divergence respectively — near-identical across all four, which is the
signature of a genuinely shared confound and confirms the mechanism.

**But it is not the whole story, and this corrects a claim in the pilot.** Normalizing
by `sub_rate` moves the cross-archaic correlation only from 0.72 to 0.70. Mutation rate
accounts for ~16% of the shared variance; the rest is shared *genealogy* — the archaics
share ancestry and are compared against the same African panel, so their divergences are
correlated under every model, including the null. "Per-window rate variation lifts all
archaics together" was true but minor. Consequently a rate denominator alone is not
enough: the null has to come from simulation, which is what §3 supplies.

A free by-product: hg19's soft-masking is RepeatMasker, so `repeat_frac` now gives the
genome-wide per-window repeat track that §7 of the previous version listed as missing.

### Polarized site patterns (`src/scan_windows.py`)
Each site is polarized against the ancestral allele (CAnc from the Altai/Denisovan INFO,
falling back to the 1000G EPO `AA`) and tallied into lineage-specific derived-allele
patterns, conditioned on all four archaics being callable and the derived allele being
absent from Africans:

| pattern | meaning |
|---|---|
| `pat_nea_all` | derived in **all three Neanderthals**, ancestral in Denisovan — the Model-A signature |
| `pat_den_only` | derived in Denisovan alone — the Model-B signature |
| `pat_all_arch` | derived in all four — ordinary Nea–Den shared ancestry |
| `pat_alt_only` / `pat_vin_only` / `pat_chag_only` | single-Neanderthal private — drift/error controls |

These count deep-branch mutations. Their **clustering**, not their mean, is what
separates introgressed haplotype blocks from incomplete lineage sorting — the mean is
confounded with the split times.

### Locus-free tests (`src/modelA.py`)
- **T1 clade clustering** — the headline. Model A predicts `pat_nea_all` clusters in
  blocks; the null predicts it is scattered. Reported as a *contrast* against control
  patterns so that rate, Ne and mappability, which cluster everything alike, cancel.
- **T2 tail symmetry** — deep-tail mass of rate-normalized depth per lineage. Rules
  Model B in or out. It **cannot** test Model A, which lifts all four lineages equally.
  The pilot's "all four focal rates are equal" observation was therefore never evidence
  against Model A, and is no longer presented as such.
- **T3 cross-lineage covariance** — the aggregate, unthresholded version of the
  co-location test.

Uncertainty throughout is a Busing–Meijer–van der Leeden weighted block jackknife over
5 Mb blocks (147 blocks on chr13–22), which respects linkage.

### Statistic design, twice corrected
Quasi-Poisson φ is **not** comparable across patterns of different abundance, since
Var = μ + αμ² gives φ = 1 + αμ. The contrast is therefore built on the mean-scale-free
α (unit-tested for invariance). A noise-corrected latent autocorrelation is reported
alongside it, because raw autocorrelation of counts is attenuated in proportion to how
rare a pattern is.

---

## 2. Observed results (chr6 + chr13–22, 50 kb, MHC excluded)

**T1 contrast: α(`pat_nea_all`) − α(`pat_alt_only`) = +0.067 ± 0.061 (z = +1.09).**

**T2:** Denisovan minus Neanderthal-mean tail mass = +0.0022 ± 0.0010 (z = +2.24) —
a weak Model-B-direction excess.

**T3:** corr(R_den, R_nea) = +0.725 ± 0.012.

### chr6 halved the T1 signal, and that matters more than the point estimate

On chr13–22 alone the contrast was +0.114 ± 0.046 (z = 2.49). Adding chr6 — **26% more
data** — dropped it to +0.067 ± 0.061 (z = 1.09). An effect that shrinks when data is
added is behaving like noise, and the leave-one-chromosome-out table makes the reason
plain: **every** chromosome dropped gives +0.03…+0.09, *except* dropping chr6, which
gives +0.114. The chr13–22 figure was the outlier, not the chr6 result.

| chromosome | n | contrast | z |  | drop | contrast | z |
|---|---|---|---|---|---|---|---|
| chr6 | 3,146 | +0.058 | +0.28 | | chr6 | **+0.114** | +2.49 |
| chr13 | 1,850 | +0.335 | +2.29 | | chr13 | +0.033 | +0.49 |
| chr14 | 1,692 | +0.196 | +2.14 | | chr14 | +0.059 | +0.92 |
| chr15 | 1,451 | +0.045 | +0.61 | | chr15 | +0.074 | +1.13 |
| chr16 | 1,381 | −0.011 | −0.10 | | chr16 | +0.079 | +1.15 |
| chr17 | 1,381 | −0.047 | −0.44 | | chr17 | +0.088 | +1.35 |
| chr18 | 1,452 | +0.300 | +1.47 | | chr18 | +0.039 | +0.60 |
| chr19 | 869 | +0.103 | +0.92 | | chr19 | +0.068 | +1.08 |
| chr20 | 1,156 | +0.068 | +0.55 | | chr20 | +0.066 | +1.01 |
| chr21 | 638 | +0.048 | +0.41 | | chr21 | +0.071 | +1.13 |
| chr22 | 600 | +0.316 | +1.40 | | chr22 | +0.057 | +0.91 |

Nine of eleven chromosomes lean positive (mean +0.128) but the scatter is large
(sd 0.136, range −0.047…+0.335). The between-chromosome scatter implies a standard error
on the mean of 0.041, against the pooled block-jackknife 0.061 — the same order, so the
jackknife is not materially understating the error.

### Polarized pattern spectrum (chr13–22 values, for reference)

| pattern | total | per window | α (relative dispersion) | lag-1 autocorr |
|---|---|---|---|---|
| `pat_nea_all` | 49,567 | 3.97 | 0.690 ± 0.037 | +0.425 |
| `pat_den_only` | 132,810 | 10.65 | 0.127 ± 0.006 | +0.289 |
| `pat_all_arch` | 14,912 | 1.20 | 2.982 ± 0.149 | +0.375 |
| `pat_alt_only` | 25,776 | 2.07 | 0.576 ± 0.031 | +0.312 |
| `pat_vin_only` | 4,091 | 0.33 | 3.330 ± 0.279 | +0.291 |
| `pat_chag_only` | 2,327 | 0.19 | 6.579 ± 2.094 | +0.247 |

### MHC positive control
Excluded from the primary scan and reported separately: 120 windows at chr6:28–34 Mb,
median `div_den_afr` 1.10× the genome-wide median. Elevated in the expected direction —
trans-species balancing selection produces genuinely ancient haplotypes — and modest at
50-kb resolution because the extreme signal sits in a narrower interval. It behaves as a
deep-but-shared region rather than an archaic-specific one, which is the correct
behaviour for a rate-normalized method.

---

## 3. Simulation calibration — what the T1 contrast means

All five models were run through the **identical** statistics module (10 seeds × 20 Mb
each), with per-window mutation rates resampled from the *measured* rate map rather than
held constant, and with the African panel matched to 1000G's 661 individuals (panel size
materially changes the "absent in Africans" filter).

| statistic | M0 null | **M1 = Model A** | M2 = Model B | M3 separate | M4 structure | **observed** |
|---|---|---|---|---|---|---|
| α(nea_all) − α(alt_only) | −0.017 | **+0.112** | −0.003 | −0.554 | +0.152 | **+0.114** |
| corr(R_den, R_nea) | 0.830 | 0.722 | 0.682 | 0.582 | 0.664 | 0.727 |
| legacy % co-located | 81% | 54% | 30% | 25% | 60% | 62% |

The T1 contrast behaves as designed: **≈0 under the null and under Model B, clearly
positive under Model A.** That is the important part — the statistic discriminates, which
`den_excess` provably could not.

The observed value does not. At **+0.067 ± 0.061** it is 1.38 SE from M0 (−0.017) and
0.74 SE from M1 (+0.112): **consistent with both.** The models are only 0.129 apart,
which is 2.1 SE at current coverage, so the data cannot yet choose between them. The
earlier chr13–22 value of +0.114 sat on M1, but as §2 shows that was a one-chromosome
accident.

### Why no claim is made

0. **The observed value is not significant** (z = 1.09) and is not stable across
   chromosomes (§2). Everything below would matter *if* it were.
1. **M4 (ancient structure) produces the same signal** (+0.152). Continuous deep
   structure without any discrete pulse mimics Model A on this statistic. This is the
   long-standing Rogers-versus-structure ambiguity, and nothing here resolves it — so
   even a significant result would not by itself support Model A over structure.
2. **The contrast depends on which drift control is used.** The three
   single-Neanderthal patterns are supposed to be exchangeable but give α = 0.58
   (Altai), 3.33 (Vindija), 6.58 (Chagyrskaya). Against Vindija or Chagyrskaya the
   contrast reverses sign (−2.64, −5.89). Altai is the only one with enough private
   sites (~26k vs ~4k and ~2k — it split from the Vindija/Chagyrskaya ancestor ~150 ka,
   giving a far longer private branch) for α to be well determined, and the simulated
   comparison uses the same control, so the contrast is at least internally consistent.
   It is still a choice, and it is reported rather than buried.
3. **The simulator does not yet reproduce the observed pattern spectrum.**
   `pat_all_arch` comes out ~8× too abundant, meaning the Neandersovan branch length
   and Ne are unfitted. The calibration is trustworthy in shape, not in absolute scale.
4. **One point estimate per model.** The z values come from a block jackknife inside a
   single simulated dataset, not from a distribution over independent replicates, so
   there is no false-positive rate or power figure yet.
   `sims/power_calibration.py` supplies that and is ready to run.

### What it is not
The contrast is **robust to data quality**, which removes the most obvious artefactual
explanation:

| subset | n | contrast | z |
|---|---|---|---|
| all windows | 12,470 | +0.114 | +2.49 |
| callable_frac ≥ 0.4 | 11,386 | +0.116 | +2.61 |
| callable_frac ≥ 0.5 | 8,752 | +0.101 | +2.17 |
| callable_frac ≥ 0.6 | 3,950 | +0.058 | +0.97 |
| balanced callability across genomes | 12,378 | +0.113 | +2.45 |
| non-artifact (Map20 ≥ 0.7, RM ≤ 0.6) | 5,199 | +0.123 | +2.05 |
| repeat_frac < 0.5 | 6,619 | +0.121 | +2.25 |

The attenuation in the strictest subset is consistent with the loss of two-thirds of the
windows rather than with the signal being quality-driven.

---

## 4. Bottom line

> *"Do Neanderthals show reproducible evidence of sharing the same candidate deeply
> divergent genomic regions observed in Denisovans, beyond what is expected under
> ordinary Neanderthal–Denisovan shared ancestry?"*

**Co-location: unanswerable, not merely negative — and it was always the wrong
question.** Under Model A the two lineages carry deep ancestry at *different* loci after
~420 ky, so the test has ~6% power by construction; and as §0 shows it has no valid null
at any conditioning level. The pilot's apparent significance was a selection artefact,
and the de-biased version is a smaller one.

**Neanderthal-clade deep-ancestry clustering: no significant signal.** The
Model-A-diagnostic statistic is +0.067 ± 0.061 (z = 1.09) — consistent with the null and
equally consistent with Model A, because at 33% of the autosomes those two hypotheses
are only 2.1 SE apart. The chr13–22 value of +0.114 (z = 2.49) did not survive the
addition of chr6.

**No claim is made in either direction.** The previous version reported a confident null
that was an artefact of its own selection rule; this version reports an honest
non-result with a known power. That is a smaller-sounding but far more useful position:
the question is now open and answerable, rather than closed on bad evidence.

---

## 5. What this bought: the question is now a sample-size problem

This is the substantive change. The old co-location test had **~6% power under a true
Model A at any genome coverage** — sequencing every chromosome would not have fixed it,
because co-location is the wrong expectation, not a sample-size shortfall. The redesigned
statistic separates M0 from M1 by 0.129 with a standard error that shrinks as √n:

| coverage | windows | SE | z under Model A | power |
|---|---|---|---|---|
| **now** (chr6 + chr13–22, 33%) | 15,616 | 0.061 | 2.11 | **68%** |
| + chr8–12 (53%) | 25,173 | 0.048 | 2.69 | 85% |
| all autosomes | 46,669 | 0.035 | 3.66 | **98%** |

So finishing the genome is now worth doing, and it was not before. Ranked by value:

1. **Acquire the remaining autosomes.** ~12 h of staged download per 5 chromosomes at
   the observed 12 min/genome, disk-bounded, fully resumable. Takes the test from 68% to
   98% power. This is now the highest-value action, which inverts the pilot's advice.
2. **`sims/power_calibration.py --reps 20 --jobs 6`** (~3 h) — replicate-level FPR and
   power. **This is now the most urgent item**, ahead of any demography work: the
   M1 − M0 separation of 0.129, on which the entire table above rests, is itself poorly
   determined. Independent 50 Mb runs put it anywhere from −0.033 to +0.189. Until the
   effect size is pinned down the power column is indicative only.
3. **Fix the AFRICAN model — not the Neandersovan branch.** `pat_all_arch` is ~8× too
   abundant in simulation, but *not* because the Neandersovan branch is mis-specified,
   and fitting that branch would be a mistake. The evidence:

   | check | result |
   |---|---|
   | `pat_all_arch`/`pat_den_only` vs `T_HN` (470/540/620/700 ka) | 0.224 / 0.408 / 0.650 / 0.828 — cleanly monotone, so the ratio *would* identify the branch |
   | observed ratio | **0.113** → extrapolates to `T_HN` ≈ 428 ka, an **~8 ky branch** against a 420 ka Nea–Den split. Impossible, and 120+ ky below the literature. |
   | aDNA damage / single-genome error | **ruled out** — transversion fraction is equal across all four patterns (0.291–0.313); a transversions-only pass moves the ratio only 0.113 → 0.122 |
   | the "absent in Africans" condition | **the cause** — it removes **94.3%** of `pat_all_arch` but only **38.3%** of `pat_den_only` |
   | ratio without that condition | **1.227**, *above* the simulated 0.650 |

   `pat_all_arch` variants arose on the Neandersovan branch (420–620 ka) and are old
   enough to also segregate in Africans; Denisovan-private variants are younger and
   mostly do not. The simulator's Africans are a single panmictic population at constant
   Ne = 15,000 and retain far less ancient variation than real African populations, so
   the filter bites far harder in the data than in the simulation. Either fit African
   structure/Ne (the survival fractions above are excellent targets), or — cheaper and
   more robust — use statistics that do not condition on an *exact-zero* African
   frequency.
4. **An explicit ancient-structure model.** M4 reproduces the Model-A signal, so even a
   98%-power result would not distinguish them. No amount of extra sequence fixes this;
   it needs a fitted structure model or an ARG method at candidate loci.
5. **Resolve the drift-control choice.** The three single-Neanderthal patterns are not
   exchangeable (α = 0.58 / 3.33 / 6.58); only Altai has enough private sites for α to
   be well determined. A control that does not depend on one genome's branch length
   would make T1 considerably more robust.

**Reproduce:** `src/stagger.sh <chroms>` → `src/scan_windows.py --jobs 4` →
`src/modelA.py --sensitivity` → `sims/superarchaic_sim.py` battery →
`sims/power_calibration.py`. Tests: `python tests/test_pipeline.py` (16 pass).
