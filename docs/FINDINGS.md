# Findings — Superarchaic Introgression in Neanderthals

**Scope: chromosomes 13–22 (~33% of the genome, 12,470 usable 50-kb windows), plus chr6
in acquisition.** Language follows the project's interpretation rules; nothing is called
"superarchaic".

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

## 2. Observed results (chr13–22, 50 kb, MHC excluded)

| pattern | total | per window | α (relative dispersion) | lag-1 autocorr |
|---|---|---|---|---|
| `pat_nea_all` | 49,567 | 3.97 | 0.690 ± 0.037 | +0.425 |
| `pat_den_only` | 132,810 | 10.65 | 0.127 ± 0.006 | +0.289 |
| `pat_all_arch` | 14,912 | 1.20 | 2.982 ± 0.149 | +0.375 |
| `pat_alt_only` | 25,776 | 2.07 | 0.576 ± 0.031 | +0.312 |
| `pat_vin_only` | 4,091 | 0.33 | 3.330 ± 0.279 | +0.291 |
| `pat_chag_only` | 2,327 | 0.19 | 6.579 ± 2.094 | +0.247 |

**T1 contrast: α(`pat_nea_all`) − α(`pat_alt_only`) = +0.114 ± 0.046 (z = +2.49).**

**T2:** Denisovan minus Neanderthal-mean tail mass = +0.0024 ± 0.0012 (z = +1.99) —
a weak Model-B-direction excess.

**T3:** corr(R_den, R_nea) = +0.727 ± 0.013.

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
positive under Model A.** The observed +0.114 sits essentially on M1 (+0.112) and away
from M0 (−0.017) and M2 (−0.003).

### Why this is *suggestive*, not a result

1. **M4 (ancient structure) produces the same signal** (+0.152). Continuous deep
   structure without any discrete pulse mimics Model A on this statistic. This is the
   long-standing Rogers-versus-structure ambiguity, and nothing here resolves it.
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

**Neanderthal-clade deep-ancestry clustering: a weak positive signal, consistent with
Model A and equally consistent with ancient structure.** The Model-A-diagnostic
statistic sits where a true Model A would put it and away from the null and Model B, is
stable across data-quality filters, but is z ≈ 2.5 on a third of the genome, depends on
the choice of drift control, rests on an unfitted demography, and does not separate
Model A from ancient population structure.

**This is a candidate signal that did not exist before, not a finding.** The previous
version reported a confident null that was an artefact; this version reports a weak
positive that is honestly bounded. Neither licenses a claim about superarchaic ancestry
in Neanderthals.

---

## 5. What would settle it

1. **`sims/power_calibration.py --reps 20 --jobs 6`** (~3 h) — replicate-level FPR and
   power, and where the observed value falls in each model's distribution. Until this
   runs, "z = 2.49" has no calibrated false-positive rate behind it.
2. **Fit the demography** so the simulated pattern spectrum matches the observed one
   (chiefly `pat_all_arch`). Power numbers are not quantitative until this is closed.
3. **An explicit ancient-structure model** fitted to the data. This is the one
   confounder that reproduces the signal, and no amount of extra sequence removes it.
4. **The remaining autosomes.** Note this is the *smallest* available gain: going from
   33% to 100% of the genome shrinks the standard error by only ~1.7×, taking z ≈ 2.5
   to ≈ 4.3 *if the effect is real*. The method changes above were worth far more than
   the data, which is why they came first.
5. **chr6/MHC as a positive control.** The MHC carries trans-species balancing
   selection — genuinely ancient haplotypes that are *not* archaic introgression. A
   correct rate-normalized method should flag it as ancient and shared, not
   archaic-specific. It is excluded from the primary scan and reported separately.

**Reproduce:** `src/stagger.sh <chroms>` → `src/scan_windows.py --jobs 4` →
`src/modelA.py --sensitivity` → `sims/superarchaic_sim.py` battery →
`sims/power_calibration.py`. Tests: `python tests/test_pipeline.py` (16 pass).
