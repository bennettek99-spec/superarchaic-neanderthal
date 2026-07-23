# Literature Review — Superarchaic ("ghost") Introgression in Neanderthals and Denisovans

**Project:** Investigating evidence for superarchaic introgression in Neanderthals using high-coverage archaic genomes.
**Deliverable 1 of 9.** Author: computational analysis for B. Kuhn. Date: 2026-07-23.

> **Scope of the question.** Multiple lines of evidence indicate that a *deeply diverged
> ("superarchaic") ghost lineage* — separated from the modern/Neanderthal/Denisovan
> ancestor on the order of ~1–2 million years ago, and often equated with *Homo erectus* or
> an *erectus*-like population — contributed ancestry to **Denisovans**. The open question this
> project targets is narrower and specific:
>
> **Did Neanderthals inherit material from the *same* superarchaic source, or is the superarchaic
> signal effectively Denisovan-specific?**
>
> - **Model A (Neandersovan):** superarchaic → the *common ancestor of Neanderthals and
>   Denisovans* ("Neandersovan") **before** they split ⇒ *both* later Neanderthals **and**
>   Denisovans should retain a detectable signal.
> - **Model B (Denisovan-only):** superarchaic → Denisovans (or the Denisovan lineage) **after**
>   the Neanderthal split ⇒ only Denisovans carry the deepest signal; apparent Neanderthal
>   similarity is ordinary shared Nea–Den ancestry / ILS.
>
> The literature is **genuinely divided** on exactly this point (Rogers et al. favour A; Hubisz
> et al. favour B), which is what makes the project worth doing. A null result (Model B, or
> "cannot distinguish") is a legitimate and valuable outcome.

---

## 1. Why this is answerable — and why it is hard

**Answerable, in principle.** Superarchaic ancestry deposits genomic segments whose lineages
coalesce with the rest of humanity *very* deeply (>1 Mya), far older than the Nea–Den split
(~400–700 kya) or ordinary incomplete lineage sorting (ILS) among Nea/Den/modern. Such segments
leave three converging signatures that high-coverage genomes can measure **directly** rather than
through modern-human proxies:

1. **Elevated pairwise sequence divergence** in windows (archaic-vs-modern, archaic-vs-archaic).
2. **Anomalously old local TMRCA / long internal branches** in local genealogies.
3. **Clusters of lineage-private derived alleles** at high local density.

A recent, directly relevant result (McAllister, Zöllner & Zhang 2026, *DEEP*) shows analytically
and by simulation that **Denisovan-mediated superarchaic introgression produces predictable shifts
in local coalescent depth that are well approximated by simple summary statistics — especially
pairwise sequence divergence** — *without* reconstructing full genealogies. This is direct support
for the window-scan design adopted here.

**Hard, in practice.** Every signature above is also produced, locally, by mundane processes:
ILS, ancient population structure, mutation-rate heterogeneity, balancing selection, low
recombination, reference bias, aDNA damage, mapping/alignment error in repetitive or duplicated
regions, and contamination. The entire methodological burden of this project is *separating a
genuine deep-ancestry signal from these alternatives* — hence the emphasis on masks, multiple
window sizes, block jackknife/bootstrap, simulation-calibrated false-positive rates, and a
conservative "Candidate Deep Divergence Region (CDDR)" label that is **never** equated with
"superarchaic" without surviving the alternative-explanation gauntlet.

---

## 2. Per-paper structured summaries

Grouped: (A) the archaic genomes themselves; (B) superarchaic/ghost inference; (C) ghost
introgression elsewhere & methodological analogs; (D) 2026 state-of-the-art directly on point;
(E) ARG/ILS method references.

### A. Foundational high-coverage archaic genomes

**A1 — Green et al. 2010, *Science* (draft Neanderthal genome).**
- *Question:* Did Neanderthals contribute to modern human genomes?
- *Data / archaic genomes:* low-coverage (~1.3×) composite Neanderthal (Vindija). Chimp outgroup.
- *Model / statistic:* introduced the **D-statistic (ABBA-BABA)**; four-population test.
- *Finding / stance:* 1–4% Neanderthal ancestry in non-Africans. **Foundational method, not about superarchaic.**
- *Limitations:* very low coverage; composite individual.

**A2 — Reich et al. 2010, *Nature* (Denisova discovery).**
- *Question:* identity/relationships of the Denisova phalanx individual.
- *Data:* ~1.9× Denisovan; established Denisovans as a sister group to Neanderthals.
- *Finding:* ~4–6% Denisovan ancestry in Melanesians. Not about superarchaic; sets up the Nea/Den/modern trichotomy.

**A3 — Meyer et al. 2012, *Science* (high-coverage Denisovan, ~30×).**
- *Question:* full-quality Denisovan genome and its relationships.
- *Data / archaic genomes:* ~30× Denisova 3; comparison to low-cov Neanderthal, modern panels.
- *Statistics:* D-statistics, F-statistics, heterozygosity, divergence.
- *Finding / stance:* first hints that Denisovans carry ancestry from a **more diverged, possibly
  archaic source**; noted the Denisovan lineage's unusual features. Seeds the "unknown archaic → Denisovan" idea.
- *Limitations:* single Denisovan; cannot yet localize deep segments confidently.

**A4 — Prüfer et al. 2014, *Nature* (Altai Neanderthal, ~52×).** — *pivotal for the ghost.*
- *Question:* high-quality Neanderthal genome; the network of archaic gene flow.
- *Data / archaic genomes:* **Altai Neanderthal (~52×)** + high-cov Denisovan + modern panels; chimp outgroup.
- *Model / statistics:* F/D-statistics, a coalescent model of four gene-flow edges; enhanced D.
- *Findings / stance:* (i) Neanderthal → modern non-Africans; (ii) Neanderthal → Denisovan;
  (iii) **gene flow into Denisovans from a population that diverged deeply — "possibly *Homo
  erectus*" — contributing on the order of a few percent of the Denisovan genome.**
  → **Supports a superarchaic → Denisovan event (Model-B-type).** Did *not* claim the same for Neanderthals.
- *Assumptions/limitations:* single genome per archaic; gene-flow proportions correlated/degenerate;
  "deeply diverged" source unsampled — inference is model-dependent.

**A5 — Prüfer et al. 2017, *Science* (Vindija 33.19 Neanderthal, ~30×).**
- *Data / archaic genomes:* **Vindija 33.19 (~30×)** + Altai + Denisovan.
- *Findings:* Vindija is closer to the Neanderthals that introgressed into modern humans than Altai
  is; revised Neanderthal contribution to non-Africans up (~1.8–2.6%). Refined Nea population history.
- *Stance on superarchaic:* not the focus; consistent with the Prüfer 2014 archaic-gene-flow network.

**A6 — Mafessoni et al. 2020, *PNAS* (Chagyrskaya 8 Neanderthal, ~27×).**
- *Data / archaic genomes:* **Chagyrskaya 8 (~27×)**, a second Altai-region Neanderthal.
- *Findings:* low heterozygosity / small Neanderthal Ne; gene flow *between* Neanderthal populations;
  Chagyrskaya closer to Vindija than to Altai. Neanderthal→modern estimates consistent with prior work.
- *Stance on superarchaic:* did not report a superarchaic-in-Neanderthal signal; **provides a
  third independent high-coverage Neanderthal — exactly the replication our Model A test needs.**

### B. Superarchaic / ghost inference in archaic hominins (the core debate)

**B1 — Rogers, Bohlender & Huff 2017, *PNAS*.**
- *Question:* early history (splits, sizes, admixture) of Neanderthals/Denisovans.
- *Data:* site-pattern frequency spectra from Altai Nea, Denisovan, Europeans, Africans.
- *Method:* **Legofit** — likelihood on the expected frequencies of nucleotide site patterns
  under an explicit demographic model.
- *Findings / stance:* introduced the **"Neandersovan"** ancestor; inferred deep splits and a
  bottleneck; groundwork for superarchaic admixture. Leans toward structure/admixture deep in the tree.

**B2 — Rogers, Harris & Achenbach 2020, *Science Advances*.** — **the Model-A anchor.**
- *Question:* did the Nea–Den ancestors interbreed with a superarchaic population?
- *Data / archaic genomes:* Altai Neanderthal + Denisovan + African + European (site-pattern spectrum).
- *Method / statistics:* **Legofit** on ~10 site patterns; model comparison (AIC-like), bootstrap.
- *Findings / stance:* **the ancestors of Neanderthals and Denisovans ("Neandersovans") interbred
  with a "superarchaic" population that separated ~2 Mya**; superarchaic Ne ≈ 20,000–50,000;
  confirms an *additional* superarchaic → Denisovan event; supports an early (~2 Mya) Neandersovan
  dispersal into Eurasia. → **Directly supports Model A: both Neanderthals and Denisovans should
  carry superarchaic ancestry.** DOI: 10.1126/sciadv.aay5483.
- *Assumptions/limitations:* strong reliance on a *pre-specified, discrete-pulse* demographic model
  with few site patterns; results depend on model choice; a single Neanderthal and single Denisovan;
  cannot by itself localize *which* Neanderthal regions are superarchaic. **Ancient structure (Model 4)
  can mimic pulse admixture in such spectra** — a key alternative this project must weigh.

**B3 — Hubisz, Williams & Siepel 2020, *PLoS Genetics*.** — **the Model-B anchor.**
- *Question:* map older/ghost gene-flow events directly in local genealogies.
- *Data / archaic genomes:* human + Altai Neanderthal + Denisovan genomes.
- *Method:* **ARGweaver-D** — demography-aware Bayesian sampling of ancestral recombination graphs;
  parses migrant lineages to give per-site introgression probabilities.
- *Findings / stance:* (i) **~3% of the *Neanderthal* genome introgressed from *ancient modern
  humans* ~200–300 kya** — note this is the *reverse* direction (H. sapiens → Neanderthal),
  **not** superarchaic → Neanderthal; (ii) **~1% of the *Denisovan* genome from an unsequenced,
  highly diverged "super-archaic" hominin**; (iii) ~15% of those superarchaic regions (≥~4 Mb)
  were passed on into modern humans. → **Supports Model B: superarchaic signal is Denisovan-specific;
  found no superarchaic-into-Neanderthal.** DOI: 10.1371/journal.pgen.1008895.
- *Assumptions/limitations:* ARG inference on few deep samples is hard; results depend on the
  demographic prior; "no evidence" for superarchaic-in-Neanderthal is not the same as evidence of
  absence — power at these depths is limited. **This is precisely the gap our direct, replicated,
  three-Neanderthal test aims to narrow.**

**B4 — Browning, Browning, Zhou, Tucci & Akey 2018, *Cell*.**
- *Question:* structure of Denisovan ancestry in living people.
- *Method:* **Sprime** (reference-free archaic-segment discovery from modern haplotypes).
- *Findings / stance:* **two distinct Denisovan components** in modern Asians/Oceanians, one more
  deeply diverged — evidence for *Denisovan internal structure/multiple sources*, which is easy to
  confuse with a superarchaic contribution. Relevant to disentangling "deep Denisovan" from "superarchaic."
- *Limitations:* modern-human-based; inferred segments, not archaic genomes directly.

**B5 — Villanea & Schraiber 2019, *Nature Ecology & Evolution*.**
- *Question:* one vs multiple Neanderthal→modern pulses.
- *Method:* **Approximate Bayesian Computation (ABC)** on the distribution of introgressed-tract counts.
- *Findings / stance:* favours **multiple episodes** of Neanderthal→modern gene flow. Not about
  superarchaic → Neanderthal, but a **methodological template** (ABC + tract statistics) for our
  simulation/inference layer, and a reminder that "one pulse" defaults are often wrong.

### C. Ghost introgression elsewhere & methodological analogs

**C1 — Hsieh et al. 2016, *Genome Research* (Central African Pygmies).**
- Whole-genome model-based inference (>60× data) → **archaic ("ghost") admixture into African
  AMH ancestors**. Establishes that ghost inference is feasible from divergence/LD summaries.
  DOI: 10.1101/gr.196634.115.

**C2 — Durvasula & Sankararaman 2020, *Science Advances* (ArchIE).**
- Reference-free ghost-segment detection (logistic/feature-based) → deeply diverged ghost ancestry
  in West Africans (~2–19% in some groups). Method analog for feature-based deep-ancestry detection.

**C3 — Teixeira et al. 2021, *Nature Ecology & Evolution* (ISEA).**
- Searched >400 modern genomes (>200 Island SE Asian) for superarchaic (H. luzonensis/floresiensis)
  admixture: **widespread Denisovan ancestry but NO substantial superarchaic signal** in living
  people beyond Denisovan. Bounds how much superarchaic reached *modern* humans; sharpens the
  distinction between "superarchaic in archaics" vs "superarchaic in moderns." DOI: 10.1038/s41559-021-01408-0.

**C4 — Pawar et al. 2023, *Nature Ecology & Evolution* (eastern gorillas).**
- **ABC + neural network** demographic inference → ~3% ghost archaic ancestry in eastern gorillas
  from a lineage diverged >3 Mya; X-chromosome depletion of introgression (parallels the archaic-human
  X-depletion). Direct methodological analog (ABC-NN) for our Models 0–4 inference. DOI: 10.1038/s41559-023-02145-2.

### D. 2026 state-of-the-art — directly on this project's question

**D1 — Zhang, Biddanda, Johnson, O'Dushlaine & Moorjani 2026, *bioRxiv* (TRACE).**
- *Method:* archaic-ancestry detection from **inferred ARGs of contemporary genomes alone** (no
  archaic reference, no unadmixed outgroup).
- *Findings / stance:* recovers known Neanderthal/Denisovan introgression; finds ghost admixture in
  Africans and non-Africans; ghost ancestry **persists in Neanderthal/Denisovan "ancestry deserts."**
  **In Oceanians, deep lineages are enriched in Denisovan — not Neanderthal — regions, "supporting a
  model of super-archaic gene flow."** → **Model-B-leaning** (superarchaic tracks Denisovan-associated
  regions). DOI: 10.64898/2026.03.03.709416. *(Preprint; not peer-reviewed.)*

**D2 — McAllister, Zöllner & Zhang 2026, *bioRxiv* (DEEP).** — **most methodologically relevant.**
- *Method:* ARG-free neural-network detection of **Denisovan-mediated superarchaic** ancestry, built
  on the analytic result that superarchaic introgression shifts **local coalescent depth**, well
  captured by **pairwise sequence divergence** in windows.
- *Findings / stance:* applied to Oceanians, Tibetans, Han → **~0.4–0.6% of genomic windows** with
  superarchaic evidence; **recurrent enrichment near the HLA locus** across populations. Frames the
  superarchaic contribution as *Denisovan-mediated* (Model-B framing for moderns). DOI: 10.64898/2026.06.25.734355.
  *(Preprint.)* → **Validates our summary-statistic choice (windowed divergence / coalescence depth)
  and flags HLA + immune loci as expected positive-control regions.**

**D3 — Fu et al. 2026, *Nature* (H. erectus enamel proteomes, China).**
- *Data:* ancient enamel proteins from six Middle-Pleistocene (~0.4 Mya) *H. erectus* specimens
  (Zhoukoudian, Hexian, Sunjiadong).
- *Finding / stance:* a shared enamel variant **AMBN(M273V) present in Denisovans and these H.
  erectus** specimens; authors conclude the **superarchaic-introgression regions in the Denisovan
  genome "are likely to have originated from H. erectus,"** and that late H. erectus may have
  coexisted/interacted with Denisovans in East Asia. → **Independent (proteomic/fossil) support that
  the Denisovan superarchaic source ≈ H. erectus; Denisovan-framed (Model-B-consistent).**
  DOI: 10.1038/s41586-026-10478-8.

### E. ARG / ILS method references (for the local-genealogy layer)

- **Rasmussen, Hubisz, Gronau & Siepel 2014** (ARGweaver) — MCMC ARG sampling; basis of ARGweaver-D.
- **Speidel et al. 2019** (Relate), **Kelleher et al. 2019 / Wohns et al. 2022** (tsinfer + tsdate),
  **Deng, Nielsen & Song 2024** (SINGER) — scalable ARG/tree-sequence inference and dating; candidate
  engines for local TMRCA at CDDRs.
- **ILS baseline:** under the multispecies coalescent, with Nea–Den split only a few hundred kya and
  large ancestral Ne, a substantial fraction of loci have genealogies discordant with the species
  tree **without any introgression** — the null every CDDR must be tested against.

---

## 3. Comparison table

| # | Paper (year) | Core question | Archaic genomes used | Method / statistic | Superarchaic → **Neanderthal** (Model A)? | Superarchaic → **Denisovan** (Model B)? | Ancient **structure** invoked? | Net stance |
|---|---|---|---|---|---|---|---|---|
| A4 | Prüfer 2014 (Nature) | Altai Nea genome; gene-flow network | Altai Nea 52×, Denisovan 30× | F/D-stats + coalescent edges | Not claimed | **Yes** (deep source, "possibly erectus") | Partly | Den-only superarchaic |
| A5 | Prüfer 2017 (Science) | Vindija genome; refine history | +Vindija 33.19 30× | F/D-stats | Not addressed | Consistent w/ 2014 | — | Neutral (Nea→modern focus) |
| A6 | Mafessoni 2020 (PNAS) | Chagyrskaya genome; Nea history | +Chagyrskaya 8 27× | het, F/D, Ne | **No signal reported** | — | — | No superarchaic-in-Nea reported |
| B2 | **Rogers 2020 (Sci Adv)** | Did Neandersovans admix w/ superarchaic? | Altai Nea, Denisovan | **Legofit** site-pattern spectrum | **YES (into Neandersovan ⇒ both)** | Yes (additional) | Alt. explanation | **Model A** |
| B3 | **Hubisz 2020 (PLoS Gen)** | Map deep gene flow in ARGs | Altai Nea, Denisovan | **ARGweaver-D** | **No** (found ancient *human*→Nea instead) | **Yes (~1%)** | — | **Model B** |
| B4 | Browning 2018 (Cell) | Structure of Denisovan ancestry | (modern haplotypes) | **Sprime** | — | 2 Denisovan sources (1 deep) | Denisovan structure | Deep Denisovan ≠ nec. superarchaic |
| B5 | Villanea 2019 (Nat E&E) | One vs many Nea pulses | (modern) | **ABC** tract counts | n/a (Nea→modern) | — | — | Methods template |
| C3 | Teixeira 2021 (Nat E&E) | Superarchaic in ISEA moderns? | (modern) | targeted ghost search | — | Denisovan yes; superarchaic **no** (in moderns) | — | Bounds superarchaic-in-moderns |
| C4 | Pawar 2023 (Nat E&E) | Ghost in gorillas | gorilla WGS | **ABC-NN** | analog | analog | analog | Methods analog (ghost feasible) |
| D1 | Zhang 2026 TRACE (bioRxiv) | Archaic ancestry from ARGs, no refs | (modern ARGs) | ARG features | deserts retain ghost | **Denisovan-region enriched** | — | **Model-B-leaning** |
| D2 | McAllister 2026 DEEP (bioRxiv) | Detect Denisovan-mediated superarchaic | (modern) | **windowed divergence / coalescent depth + NN** | — | **Yes**, 0.4–0.6% windows, HLA | — | **Model-B framing; validates our stat** |
| D3 | Fu 2026 (Nature) | Molecular ID of superarchaic source | H. erectus proteomes | paleoproteomics | — | Denisovan superarchaic ≈ **H. erectus** | — | Model-B-consistent (fossil) |

*(Preprints D1–D2 are not peer-reviewed and are cited as such.)*

---

## 4. Synthesis — where the field stands, and the gap this project fills

1. **A superarchaic → Denisovan contribution is broadly (not universally) accepted** (Prüfer 2014;
   Hubisz 2020; and now proteomic support that the source resembles *H. erectus*, Fu 2026).
2. **Whether Neanderthals share that same superarchaic source is unresolved and contested.**
   - **Rogers 2020** says *yes* (into the Neandersovan ancestor ⇒ Model A) — but from a *few site
     patterns under a pre-specified discrete-admixture model*, where **ancient structure (Model 4)
     is a well-known mimic**.
   - **Hubisz 2020** says *no* — it localizes superarchaic ancestry to Denisovans and finds the
     Neanderthal deep signal is instead *ancient-human → Neanderthal*. But ARG power at ~1–2 Mya
     depths with one Neanderthal is limited; "no evidence" ≠ "evidence of no."
3. **Modern-human-based 2026 methods (TRACE, DEEP)** consistently find the superarchaic signal
   **enriched in Denisovan-associated regions**, i.e. Model-B-leaning — but they infer through
   living genomes, not the archaic genomes themselves.

**The gap:** no study has done a *direct, replicated, region-level* test on the **high-coverage
archaic genomes** asking: *do candidate deep-divergence regions discovered in the Denisovan genome
recur in Altai, Vindija, and Chagyrskaya Neanderthals more than ILS/structure predict?* Three
independent high-coverage Neanderthals (Altai, Vindija, Chagyrskaya) now exist — enough to demand
**replication** of any Neanderthal signal and to separate a genuine shared-ancestry deposit from
one-genome artifacts. That is exactly this project's design.

## 5. Concrete, testable predictions that separate Model A from Model B

For a set of **Candidate Deep Divergence Regions (CDDRs)** discovered in the Denisovan genome:

| Observable | Model A (Neandersovan) | Model B (Denisovan-only) | Null (ILS / structure only) |
|---|---|---|---|
| CDDR overlap with elevated Nea–modern divergence | Present in **all 3** Neanderthals, correlated across them | Absent or ILS-level only | Scattered, not replicated |
| Sharing category (Cat 1–4) | Excess of **Cat 1** (Den + all Nea) beyond ILS expectation | Excess of **Cat 3** (Den-only) | Matches coalescent-simulation ILS rates |
| Local TMRCA at CDDRs (Nea lineages) | A subset **>1 Mya** deep, shared across Neanderthals | Nea TMRCA ≈ Nea–Den split depth | Old tail explained by ancestral Ne alone |
| Neanderthal-private derived-allele clusters in CDDRs | Present, replicated | Absent | Poisson-scatter background |
| Simulated FPR/power (Models 0–4 → same pipeline) | Distinguishes A from 0/4 above a power threshold | Distinguishes B | Sets the decision threshold |

**Decision rule (success criterion).** The project answers *yes* only if CDDRs show **reproducible,
multi-Neanderthal, simulation-calibrated** deep-divergence sharing **beyond** the ILS/structure null
and **beyond** what ordinary Nea–Den shared ancestry predicts — surviving all masks, window sizes,
and alternative-explanation checks. Otherwise the honest answer is *no* / *cannot distinguish*, which
would align with Hubisz 2020, TRACE, and DEEP.

## 6. Methodological takeaways adopted into the pipeline

- **Primary window statistic = pairwise sequence divergence + a coalescent-depth proxy** (per DEEP's
  analytic result), computed at **20/50/100 kb**, on **callable sites only** (FilterBed mq25/mapab100
  masks), polarized by a chimp/ancestral allele.
- **Positive controls we must recover:** HLA/immune loci (DEEP), and known adaptive-introgression /
  deep loci; **negative controls:** shared-archaic-fixed sites, and simulated Model-0 (no superarchaic).
- **Replication is mandatory:** any Neanderthal signal must appear in **≥2 of {Altai, Vindija,
  Chagyrskaya}** and correlate across them, or it is treated as artifact.
- **Simulation-first thresholds:** Models 0–4 (msprime) projected through the identical pipeline set
  the false-positive rate and the depth/divergence cutoffs *before* any real region is called.
- **ILS and ancient structure are the headline nulls,** not afterthoughts — Rogers-style pulse signals
  and structure signals can be nearly indistinguishable in low-dimensional summaries; only the
  *region-level, replicated* pattern can help separate them, and even then modestly.

---

### Sources & verification
Facts and figures above were drawn from the primary papers and cross-checked via **PubMed** and
**Consensus** (2026-07). Key DOIs: Rogers 2020 `10.1126/sciadv.aay5483`; Hubisz 2020
`10.1371/journal.pgen.1008895`; Teixeira 2021 `10.1038/s41559-021-01408-0`; Pawar 2023
`10.1038/s41559-023-02145-2`; Hsieh 2016 `10.1101/gr.196634.115`; TRACE 2026
`10.64898/2026.03.03.709416`; DEEP 2026 `10.64898/2026.06.25.734355`; Fu 2026 `10.1038/s41586-026-10478-8`.
Prüfer 2014 `10.1038/nature12886`; Prüfer 2017 `10.1126/science.aao1887`; Meyer 2012
`10.1126/science.1224344`; Mafessoni 2020 `10.1073/pnas.2004944117`. Preprints (TRACE, DEEP) are
explicitly flagged as non-peer-reviewed.
