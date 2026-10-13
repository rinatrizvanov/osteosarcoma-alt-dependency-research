# Literature Index — TCAB1 / 17p13.1 / ALT / BIR

**Purpose.** Context-efficient substitute for the PDFs in `Literature/`. Load this
instead of the papers. Built 2026-10-12 from full-text extraction of every readable
PDF in `Literature/` and `Literature/new/osteo-alt-tcab/`.

**Corpus.** 87 readable PDFs (80 in `new/`, 7 in the original folder) + 42 DOIs whose
PDFs could not be retrieved (listed in §"Unable to read"). Relative paths below are
from `Literature/`.

**How entries are tiered.** A flat 150-word entry for all 87 papers would be ~13k
words, which defeats the purpose. Tiers 1–2 get full entries; Tier 3 (off-axis corpus
— the scan's false positives, retrieved because a keyword matched, not because the
paper is about our question) gets one line each. Every paper appears exactly once.

**Relevance tags.** `(a)` explains the Stage 1 null · `(b)` suggests a variant
hypothesis · `(c)` constrains interpretation · `(d)` adjacent context.
`⚠ CONTRADICTS PRESPEC` = conflicts with a stated assumption in
`docs/analytical_prespec.md`; see §"Prespec conflicts" for the consolidated list.

**Stage 1 result this index is written against.** Primary 8v8 median split on WRAP53
log2(TPM+1) within 16 ALT+ DepMap osteosarcoma lines; BIR module (POLD3, RAD52,
SMARCAL1) minus HR module (BRCA1, BRCA2, RAD51) Chronos composite. Cliff's δ = **+0.281**
[−0.313, +0.844], p = 0.382 — *opposite* to the predicted direction. Verdict
`indeterminate`. Only 1 of the 16 ALT+ lines was ATRX-altered.

---

## Priority Reading

Six papers. Read in this order; they are sufficient to understand why Stage 1 came
back null and what the live alternatives are.

1. **Keegan et al. (2026) BMC Cancer** — the hypothesis' own source. TCAB1-inactivating
   SVs occur in ~40% of ALT+ osteosarcoma and are framed as *ALT-initiating*. If TCAB1
   loss causes ALT, conditioning on ALT+ conditions on the mediator. This is the single
   best explanation of the null.
2. **Taglialatela et al. (2026) Genes Dev** — the design Stage 1 was layered on top of,
   via the same DepMap. Crucially: "suppressing ALT abrogates their dependency on
   SMARCAL1" — the dependency tracks ALT *activity*, not upstream genotype.
3. **Bennett et al. (2026) bioRxiv (ALTitude)** — source of the silver ALT labels in our
   frozen manifest, and of the SMARCAL1 dependency claim Stage 1 tried to refine.
4. **Liu et al. (2023) Neuro-Oncology** — SMARCAL1 *loss* permits ALT, demonstrated in
   GBM **and osteosarcoma** cells. Inverts the polarity assumed for SMARCAL1 in our BIR
   module. ⚠ CONTRADICTS PRESPEC
5. **Henriksson et al. (2014) Genes Dev** — WRAP53β knockdown impairs HR **and** NHEJ via
   RNF8/ubiquitin signalling, cutting 53BP1, BRCA1 and RAD51 recruitment together. Kills
   the "selective HR crippling channels cells into BIR" mechanism. ⚠ CONTRADICTS PRESPEC
6. **Wu et al. (2025) Nat Commun (976 lines)** — telomere maintenance is not binary:
   TMM-double-negative and TMM-double-positive lines exist. Our ALT+/ALT− arm assignment
   is coarser than the biology, and this is a ready-made replication/extension substrate.

Runner-up worth knowing: **Bythell-Douglas et al. (2025) EMBO J** — FANCM is synthetic
lethal with BRCA1, SMARCAL1 *and* RAD52 and in ALT+ cancers. FANCM was one of only two
bootstrap-stable genes (with TERT) in the Aim 2 19-gene signature. Strongest single
candidate for a replacement endpoint.

---

## Prespec conflicts (consolidated)

Eight places where the corpus contradicts or materially constrains a stated assumption
in `docs/analytical_prespec.md`. Each is sourced to a paper read this session.

| # | Prespec assumption | What the literature says | Source |
|---|---|---|---|
| 1 | SMARCAL1 is a BIR effector whose essentiality rises with BIR reliance | SMARCAL1 *loss-of-function* permits/promotes ALT; re-expression suppresses ALT in GBM **and osteosarcoma** cells | Liu 2023; Brosnan-Cashman 2021; Sobieszek 2026 |
| 2 | WRAP53β loss selectively cripples classical HR, forcing cells onto BIR | WRAP53β knockdown impairs HR **and** NHEJ; it acts upstream at RNF8/ubiquitin signalling, reducing 53BP1, BRCA1 and RAD51 recruitment together — no selectivity | Henriksson 2014 |
| 3 | BRCA1/BRCA2/RAD51 are a neutral "canonical HR" comparator | BRCA2 loss *drives* ALT/BIR (G4 stabilisation → TERRA R-loops → LLPS → APBs → BIR telomere synthesis). BRCA2 KO can phenocopy ALT induction | Lee 2024 (NAR) |
| 4 | BIR and HR modules are orthogonal, so their difference is a clean contrast | RAD51, BRCA1 and BRCA2 actively *suppress* the POLQ/LIG3/PRIMPOL restart routes; BIR itself is MUS81–RAD52–POLD3. The arms are coupled with opposite sign | Batenburg 2023 |
| 5 | WRAP53 mRNA tertile is a usable proxy for TCAB1 functional status | TCAB1 protein level is set post-translationally by the CCT/TRiC chaperonin; HHS point mutants destabilise the protein with normal transcript | Jiang 2025; Bergstrand 2020 |
| 6 | TCAB1 status is the stratifier of interest within ALT+ lines | SMARCAL1 dependency tracks ALT activity and, in glioma, stratifies on **ATRX**. Our R5 ATRX check was vacuous (1/16 lines altered) | Taglialatela 2026; Brown 2026 |
| 7 | ALT status is binary (gold/silver positive vs negative) | TMM-DN and TMM-DP cell lines exist, including lines that are mixtures of ALT and TEL subclones | Wu 2025 (Nat Commun) |
| 8 | TCAB1-low vs TCAB1-intact is an exposure contrast within ALT+ | TCAB1 SV inactivation is framed as an **early, ALT-initiating** event in ~40% of ALT+ OS tumours — i.e. the mediator we conditioned on | Keegan 2026 (both versions) |

---

## Tier 1 — Core corpus

### 1. Keegan et al. (2026) BMC Cancer — TCAB1-inactivating SVs in osteosarcoma
https://doi.org/10.1186/s12885-026-16250-5
`new/osteo-alt-tcab/Keegan-2026-Complex_structural_variations_functionally_inactivate_the_te__s12885_026_16250_5/paper.pdf`
WGS of osteosarcoma samples found SVs in TP53 intron 1 that simultaneously inactivate
TP53 **and** TCAB1, exploiting the head-to-head 1,419-nt TP53/WRAP53 overlap at 17p13.1
(hg38 chr17:7,686,071–7,703,502). Validated in a 16-model GCCRI PDX cohort by TRAP
(telomerase-negative), TEL-FISH (ultrabright foci) and TP53/TCAB1 IHC (complete loss in
OS33, OS37, OS39, OS52, 1981; reduction in 0560-000).
**Quant.** ALT in 75% of pediatric OS; TCAB1 SVs in ~40% of ALT+ OS tumours; 16 PDX
models; SJSA1 as telomerase-positive / TCAB1-positive control.
**Relevance:** (a) frames TCAB1 loss as ALT-*initiating*, not an ALT modifier — so
conditioning on ALT+ conditions on the mediator. ⚠ CONTRADICTS PRESPEC (#8)

### 2. Keegan et al. (2026) bioRxiv — preprint of the above
https://doi.org/10.64898/2026.01.02.691864
`new/osteo-alt-tcab/Keegan-2026-Functional_inactivation_of_the_telomerase_chaperone_TCAB1_pr__2026_01_02_691864/paper.pdf`
· also at `4 TCAB.pdf`
Preprint version. Same content as the BMC Cancer paper; identical extracted text, zero
mentions of BIR, POLD3, RAD52 or WRAP53β-as-HR-factor. Cite the published version.
**Note.** `4 TCAB.pdf` and the `new/` copy are the same preprint downloaded twice —
identical text, differing PDF metadata (same byte length, different sha256).
**Relevance:** (a) the absence of any BIR framing in the source paper is itself evidence
that the BIR-module link was our inference, not Keegan's claim.

### 3. Taglialatela et al. (2026) Genes Dev — SMARCAL1 as an ALT target
https://doi.org/10.1101/gad.353377.125
`new/osteo-alt-tcab/Taglialatela-2026-SMARCAL1_is_a_candidate_therapeutic_target_for_ALT_positive__gad_353377_125/paper.pdf`
Used the Cancer Dependency Map to identify SMARCAL1 as a top selective dependency in
telomerase-negative tumours, then confirmed in ALT+/ALT− panels and osteosarcoma PDX
cells. Mechanism: SMARCAL1 loss → telomeric ssDNA accumulation (PRIMPOL repriming-
dependent); its annealing activity counteracts BLM unwinding. Depletion induces
senescence, creating a senolytic combination.
**Quant.** >50% of osteosarcomas are ALT+; OS 5-yr survival ~65%. Builds on Dharia 2021
pediatric CRISPR screens flagging SMARCAL1 in OS.
**Relevance:** (c) "suppressing ALT abrogates their dependency on SMARCAL1" — the
dependency is downstream of ALT activity, so a TCAB1 split within ALT+ lines should not
move it. Directly predicts our null. ⚠ CONTRADICTS PRESPEC (#6)

### 4. Taglialatela et al. (2026) bioRxiv — preprint of the above
https://doi.org/10.64898/2026.02.15.704070
`new/osteo-alt-tcab/Taglialatela-2026-SMARCAL1_is_a_candidate_therapeutic_target_for_ALT_positive__2026_02_15_704070/paper.pdf`
Preprint; near-identical term profile to the published version (different file, not a
byte duplicate). Cite Genes Dev.
**Relevance:** (d) provenance only.

### 5. Bennett et al. (2026) bioRxiv — ALTitude
https://doi.org/10.64898/2026.03.04.709438
`new/osteo-alt-tcab/Bennett-2026-Tabular_foundation_model_predicts_alternative_lengthening_of__2026_03_04_709438/paper.pdf`
· also at `5-important.pdf`
Tabular foundation model inferring ALT status from cell-line WGS without paired germline.
Deployed across DepMap, **doubling** the number of known ALT+ models. Integration with
CRISPR-Cas9 screens yields selective SMARCAL1 dependency in ALT+ lines; acute depletion
→ G2/M arrest, mitotic catastrophe, death.
**Quant.** SMARCAL1 mentioned 264×, ATRX 161×, osteosarcoma 64×; POLD3, RAD52, BIR and
WRAP53/TCAB1 appear **0 times**. Source of the `alt_silver` labels in our frozen manifest.
**Relevance:** (a) the published SMARCAL1 result is unstratified by TCAB1 *and makes no
BIR-module claim*. Our composite was testing a link this paper never asserts.
**Note.** `5-important.pdf` is the same preprint, different download instance.

### 6. Liu et al. (2023) Neuro-Oncology — SMARCAL1 LOF promotes ALT
https://doi.org/10.1093/neuonc/noad022
`new/osteo-alt-tcab/Liu-2023-Cancer_associated_SMARCAL1_loss_of_function_mutations_promot__neuonc_noad022/paper.pdf`
(PDF in folder is the 2022 bioRxiv version, 10.1101/2022.09.19.507319; same study.)
In a patient-derived ALT+ GBM line with native SMARCAL1 deficiency, inducible **rescue**
of SMARCAL1 suppressed ALT indicators and blocked de novo telomere synthesis in GBM
**and osteosarcoma cells**. SMARCAL1-deficient ALT+ cells propagate serially in vivo
without detectable telomerase.
**Quant.** ALT in ~10–15% of cancers; SMARCAL1 LOF found in ATRX-expressing ALT+ adult GBM.
**Relevance:** (c) SMARCAL1 polarity is the opposite of what the BIR module assumed — in
OS cells SMARCAL1 presence *suppresses* ALT. ⚠ CONTRADICTS PRESPEC (#1)

### 7. Brown et al. (2026) Neuro-Oncology — SMARCAL1 synthetic lethality, ATRX-stratified
https://doi.org/10.1093/neuonc/noaf300
`new/osteo-alt-tcab/Brown-2026-SMARCAL1_is_a_targetable_synthetic_lethal_therapeutic_vulner__neuonc_noaf300/paper.pdf`
DepMap CRISPR screens + patient-derived lines and xenografts identify SMARCAL1 as a
synthetic lethal vulnerability specifically in **ATRX-deficient** ALT+ glioma. SMARCAL1
localises to APBs; depletion hyperactivates ALT, causes G2-phase DSBs and mitotic
catastrophe; inducible depletion prolongs survival in intracranial IDH-mutant
astrocytoma xenografts.
**Quant.** ~10% of cancers use ALT; 28(4):895–910.
**Relevance:** (b) the already-published stratifier for SMARCAL1 dependency is **ATRX**,
not TCAB1 — and our R5 ATRX check was vacuous because only 1/16 ALT+ OS lines was
ATRX-altered. Explains why R5 carried no information. ⚠ CONTRADICTS PRESPEC (#6)

### 8. Henriksson et al. (2014) Genes Dev — WRAP53β orchestrates the DSB ubiquitin response
https://doi.org/10.1101/gad.246546.114
`new/osteo-alt-tcab/Henriksson-2014-The_scaffold_protein_WRAP53_orchestrates_the_ubiquitin_respo__gad_246546_114/paper.pdf`
WRAP53β localises to DSBs (ATM-, H2AX-, MDC1-dependent) and targets the E3 ligase RNF8
to lesions by bridging RNF8 and MDC1 through its WD40 domain. It thereby controls
ubiquitylation at damage sites and the downstream assembly of 53BP1, BRCA1 **and** RAD51.
Knockdown impairs DSB repair by **both HR and NHEJ**, causes spontaneous breaks and
delays recovery from radiation-induced arrest.
**Relevance:** (a) the mechanistic premise of Stage 1 — that WRAP53β loss selectively
disables classical HR and forces BIR — does not hold. The defect is upstream and global.
⚠ CONTRADICTS PRESPEC (#2)

### 9. Wu et al. (2025) Nat Commun — 976-line telomere maintenance resource
https://doi.org/10.1038/s41467-025-67190-w
`new/osteo-alt-tcab/Wu-2025-Large_scale_drug_sensitivity_gene_dependency_and_proteogenom__s41467_025_67190_w/paper.pdf`
Pan-cancer telomere-biology metrics for 976 cell lines across 28 tissue types and >60
cancer types, with DIA-MS proteomics plus integrated multi-omic, drug-sensitivity and
CRISPR-KO dependency data. Identifies ALT, TMM-double-negative (DN) and TMM-double-
positive (DP) lines — including DP lines that are mixtures of ALT and TEL subclones, and
EverShorter-Telomere lines with no TMM at all.
**Relevance:** (b)(c) binary ALT labelling is the wrong data model; this is the obvious
substrate for a continuous-ALT-activity reanalysis with proteomics, which also addresses
the WRAP53-mRNA-proxy problem. ⚠ CONTRADICTS PRESPEC (#7)

### 10. Lee et al. (2024) Nucleic Acids Res — BRCA2 loss drives ALT via G4/R-loop LLPS
https://doi.org/10.1093/nar/gkae251
`new/osteo-alt-tcab/Lee-2024-Disruption_of_G_quadruplex_dynamicity_by_BRCA2_abrogation_in__nar_gkae251/paper.pdf`
BRCA2 depletion stabilises telomeric G-quadruplexes, accumulating TERRA R-loops that
trigger liquid–liquid phase separation and APB assembly; disrupting R-loops abolishes
LLPS and impairs telomere synthesis, and engineered telomeric LLPS restores it. A subset
of BRCA2-null cells bypasses senescence into BIR-mediated telomere synthesis. TERRA
R-loops also recruit PRC2.
**Relevance:** (c) BRCA2 cannot serve as a neutral HR comparator inside an ALT+ panel —
its loss is itself pro-ALT/pro-BIR. Note BRCA2 was one of only two genes with a negative
δ in our result. ⚠ CONTRADICTS PRESPEC (#3)

### 11. Batenburg et al. (2023) Int J Mol Sci — CSB and BIR pathway choice
https://doi.org/10.3390/ijms241512419
`new/osteo-alt-tcab/Batenburg-2023-CSB_Regulates_Pathway_Choice_in_Response_to_DNA_Replication__ijms241512419/paper.pdf`
CSB inhibits PRIMPOL-dependent repriming at low camptothecin and, at high CPT, promotes
replication restart through **MUS81–RAD52–POLD3-dependent BIR**. Without CSB, restart
proceeds via POLQ–LIG3, LIG4 or PRIMPOL routes that are inhibited respectively by RAD51,
BRCA1 and BRCA2.
**Relevance:** (c) gives the canonical BIR module composition (validating POLD3/RAD52,
and arguing MUS81 belonged in the set), and shows the HR arm actively *suppresses*
alternative restart — so BIR−HR is a coupled, not orthogonal, difference.
⚠ CONTRADICTS PRESPEC (#4)

### 12. Jiang et al. (2025) Stem Cell Reports — CCT5/CCT7 control TCAB1 protein
https://doi.org/10.1016/j.stemcr.2025.102750
`new/osteo-alt-tcab/Jiang-2025-Chaperonin_proteins_CCT5_and_CCT7_epigenetically_restrict_th__j_stemcr_2025_102750/paper.pdf`
CCT5 or CCT7 depletion causes **TCAB1 protein loss** and impaired telomere maintenance;
in ESCs this disrupts telomere length homeostasis, triggers DDR, drives epigenetic
reprogramming and generates 2-cell-like cells.
**Relevance:** (c) TCAB1 abundance is set post-translationally by the CCT/TRiC
chaperonin. A WRAP53-transcript tertile is blind to this entire axis — a concrete
mechanism by which the Stage 1 stratifier could misclassify. ⚠ CONTRADICTS PRESPEC (#5)

### 13. Bergstrand et al. (2020) Cell Death Dis — biallelic WRAP53 mutations in HHS
https://doi.org/10.1038/s41419-020-2421-4
`new/osteo-alt-tcab/Bergstrand-2020-Biallelic_mutations_in_WRAP53_result_in_dysfunctional_telome__s41419_020_2421_4/paper.pdf`
A 3-year-old with Hoyeraal–Hreidarsson syndrome carried biallelic germline WRAP53
mutations (L283F, R398W) in the WD40 scaffold domain. Both destabilise and mislocalise
WRAP53β via misfolding by the TRiC chaperonin; mutants cannot elongate telomeres,
maintain Cajal bodies or repair DSBs.
**Quant.** n=1 index case; ~half of HHS cases have germline telomere-biology mutations.
**Relevance:** (c) demonstrates functional TCAB1 loss with intact transcript — the
clearest single counterexample to expression-as-proxy. ⚠ CONTRADICTS PRESPEC (#5)

### 14. Bythell-Douglas et al. (2025) EMBO J — FANCM structure and synthetic lethality
https://doi.org/10.1038/s44318-025-00468-3
`new/osteo-alt-tcab/BythellDouglas-2025-Structural_basis_of_Fanconi_anemia_pathway_activation_by_FAN__s44318_025_00468_3/paper.pdf`
Crystal structures of FANCM's N-terminal ATP-dependent translocase domain (2.2 Å) and
C-terminal region. Frames FANCM as a therapeutic target on the basis of synthetic lethal
interactions with **BRCA1, SMARCAL1 and RAD52**, and in ALT+ cancers.
**Relevance:** (b) FANCM and TERT were the only two bootstrap-stable genes in the Aim 2
19-gene signature. A FANCM-centred dependency endpoint is better grounded in both our
own prior result and this literature than the BIR−HR composite was.

### 15. Shah et al. (2024) Nat Commun — 53BP1 loss causes BIR hyperrecombination
https://doi.org/10.1038/s41467-024-52916-z
`new/osteo-alt-tcab/Shah-2024-53BP1_deficiency_leads_to_hyperrecombination_using_break_ind__s41467_024_52916_z/paper.pdf`
53BP1 suppresses BIR after end resection, distinct from its end-protection role. Loss of
53BP1 induces BIR-like hyperrecombination dependent on Polα-primase end fill-in on ssDNA
overhangs, leading to PCNA ubiquitination and PIF1 recruitment. At broken forks SMARCAD1
displaces 53BP1 to let ubiquitinated PCNA and PIF1 localise. Hyper-BIR produces template
switching and large deletions; 53BP1 and BIR are synthetic lethal.
**Relevance:** (b) Henriksson shows WRAP53β loss *reduces* 53BP1 assembly at DSBs; Shah
shows reduced 53BP1 *promotes* BIR and creates dependence on PIF1/PCNA-ubiquitination.
That is a mechanistically sharper variant endpoint than POLD3/RAD52/SMARCAL1 —
**PIF1, RAD18, POLA1, SMARCAD1** — and it is the strongest "right idea, wrong genes"
reading of the Stage 1 null.

### 16. Wu et al. (2023) Nat Commun — MiDAS requires Pol ζ then Pol δ
https://doi.org/10.1038/s41467-023-35992-5
`new/osteo-alt-tcab/Wu-2023-Mitotic_DNA_synthesis_in_response_to_replication_stress_requ__s41467_023_35992_5/paper.pdf`
MiDAS rescues under-replicated regions in mitosis. Prior work had it as POLD3-dependent
with BIR-like features; this paper shows human cells need both gap-filling (REV1, Pol ζ)
and replicative (Pol δ) polymerases, with a polymerase switch required for new synthesis
at under-replicated loci. Oncogene-activated cancer-cell survival depends on it.
**Quant.** POLD3/POLD1 mentioned 131×.
**Relevance:** (c) POLD3 essentiality reports a general oncogenic-replication-stress
requirement, not an ALT/BIR-specific one — a reason POLD3 dependency would not track
TCAB1 status within a uniformly stressed ALT+ panel.

### 17. Richards et al. (2023) Nat Commun — PARP1/2 and PolD3 ADP-ribosylation gate BIR
https://doi.org/10.1038/s41467-023-40071-w
`new/osteo-alt-tcab/Richards-2023-Regulation_of_Rad52_dependent_replication_fork_recovery_thro__s41467_023_40071_w/paper.pdf`
PARP1/PARP2 regulate Rad52-dependent fork repair to maintain viability **when HR is
dysfunctional**. Mre11 and ATM induce PARP activity under replication stress, promoting
BIR via RAD52 assembly at stalled forks. Site-specific ADP-ribosylation of PolD3 is
required for BIR activity and fork recovery.
**Relevance:** (b) the closest mechanistic match to the Stage 1 intuition (HR-deficient
cells lean on RAD52/POLD3 BIR) — but the switch is **post-translational**
(ADP-ribosylation), which CRISPR gene-effect of POLD3 would not resolve. Suggests a
PARP-sensitivity rather than essentiality readout.

### 18. Wu et al. (2025) Cell Rep — BIR repairs R-loop-associated deDSBs
https://doi.org/10.1016/j.celrep.2025.116386
`new/osteo-alt-tcab/Wu-2025-Break_induced_replication_is_activated_to_repair_R_loop_asso__j_celrep_2025_116386/paper.pdf`
Senataxin (SETX) loss induces BIR-mediated hyperrecombination at R-loop-associated
double-ended DSBs: non-canonical hyper-resection requiring **RAD52 and XPF**, stalled
Polα-primase fill-in from hybrid accumulation, then PCNA ubiquitination and PIF1 loading.
SETX dysfunction is synthetic lethal with loss of PIF1, RAD52 or XPF.
**Relevance:** (b) extends BIR beyond collapsed forks to R-loop lesions, and names a
concrete synthetic-lethal triplet (PIF1/RAD52/XPF). Pairs with Lee 2024's TERRA-R-loop
ALT mechanism.

### 19. Muoio et al. (2024) Nat Commun — PARP2 promotes BIR-mediated telomere fragility
https://doi.org/10.1038/s41467-024-47222-7
`new/osteo-alt-tcab/Muoio-2024-PARP2_promotes_Break_Induced_Replication_mediated_telomere_f__s41467_024_47222_7/paper.pdf`
PARP2 promotes replication-stress-induced telomere fragility and prevents telomere loss
under chronic oxidative lesions and BLM depletion. Fragility arises from BIR activity at
telomeres.
**Quant.** POLD3 97×, BIR 69×.
**Relevance:** (c)(b) telomeric BIR is PARP2- and BLM-context-dependent; supports a
drug-sensitivity (PARP) rather than gene-essentiality framing for the next step.

### 20. Broderick et al. (2023) Nat Commun — EXD2 dictates ALT BIR pathway choice
https://doi.org/10.1038/s41467-023-38029-z
`new/osteo-alt-tcab/Broderick-2023-Pathway_choice_in_the_alternative_telomere_lengthening_in_ne__s41467_023_38029_z/paper.pdf`
Establishes that ALT relies on **two distinct BIR processes** and that the EXD2 nuclease
is recruited to ALT telomeres to direct which one runs. EXD2 loss → telomere shortening,
elevated telomeric sister-chromatid exchange, C-circle formation and BIR-mediated
telomere synthesis.
**Quant.** BIR 64×, RAD52 78×, SMARCAL1 35×, C-circles 22×.
**Relevance:** (b)(a) "BIR" is not one module. A composite averaging POLD3/RAD52/SMARCAL1
mixes two pathways with potentially opposite dependencies — a direct candidate
explanation for a near-zero, wide-CI composite effect.

### 21. Cox, Maréchal & Flynn (2016) Cell Rep — SMARCAL1 resolves ALT telomere replication stress
https://doi.org/10.1016/j.celrep.2016.01.011
`new/osteo-alt-tcab/Cox-2016-SMARCAL1_Resolves_Replication_Stress_at_ALT_Telomeres__j_celrep_2016_01_011/paper.pdf`
SMARCAL1 associates with ALT telomeres to resolve replication stress and maintain
telomere stability; without it, stalled forks deteriorate into DSBs and chromosome
fusions. Establishes replication-stress resolution as a required step in ALT.
**Relevance:** (c) Flynn-lab origin of the SMARCAL1–ALT link, predating the DepMap
dependency framing. Read with Liu 2023 to see that "SMARCAL1 protects ALT telomeres" and
"SMARCAL1 loss permits ALT" are both supported — the polarity is context-dependent, which
is itself a reason a single-sign prediction was fragile.

### 22. Brosnan-Cashman et al. (2021) Mod Pathol — SMARCAL1 loss and ALT in giant-cell GBM
https://doi.org/10.1038/s41379-021-00841-7
`new/osteo-alt-tcab/BrosnanCashman-2021-SMARCAL1_loss_and_alternative_lengthening_of_telomeres_ALT_a__s41379_021_00841_7/paper.pdf`
Telomere-specific FISH on 44 cases (15 classic giant-cell GBM, 28 GBM with giant-cell
features, 1 anaplastic astrocytoma): **19/44 (43%) ALT+**. Of the ALT+ cases only
**9/19 (47.4%)** had ATRX loss by IHC, and **6/19 (31.5%)** showed loss or
mislocalisation of SMARCAL1.
**Relevance:** (c) quantifies SMARCAL1 *loss* in ALT+ tumours — i.e. ALT+ cells can lack
SMARCAL1 entirely, which is incompatible with treating SMARCAL1 essentiality as a
monotone BIR-reliance readout. ⚠ CONTRADICTS PRESPEC (#1)

### 23. Carson & Flynn (2023) Curr Opin Pharmacol — ALT vulnerabilities
https://doi.org/10.1016/j.coph.2023.102380
`3 Alt.pdf`
Review from the Flynn lab (BU) on therapeutic targeting of ALT. ALT is BIR-dependent and
relies on baseline replication stress to initiate elongation, giving elevated DNA damage
as an exploitable vulnerability; strategies tip the balance toward either inhibiting or
exacerbating ALT.
**Quant.** ALT 139×, ATRX 52×, RAD52 17×; "BRCA" appears **0 times**.
**Relevance:** (a) RAD51 and RAD52 are described as regulating *distinct branches* of ALT,
and the canonical BRCA axis is absent from the ALT-vulnerability frame entirely. The
BIR-vs-HR contrast Stage 1 built is not how this literature carves the problem.

### 24. Hoang & O'Sullivan (2020) Trends Cancer — ALT: building bridges
https://doi.org/10.1016/j.trecan.2019.12.009
`2 Alt.pdf`
Review of ALT telomere nucleoprotein architecture: chromatin mismanagement fosters cycles
of damage and replicative stress that activate homology-directed repair; chromatin
remodeller mutations determine ALT emergence and survival.
**Quant.** ALT 179×, ATRX 26×, BIR 39×, RAD52 25×, telomerase/TERT 47×.
**Relevance:** (d)(c) best single orientation to ALT mechanism; the emphasis on chromatin
remodellers (ATRX/DAXX) rather than repair-gene dosage is the field's default stratifier
and the one our panel could not test.

### 25. Folini et al. (2026) Res Sq preprint — C-circles as ALT biomarkers
https://doi.org/10.21203/rs.3.rs-11116403/v1
`new/osteo-alt-tcab/Folini-2026-Translational_relevance_of_Telomeric_C_circles_in_Cancer_fro__rs_3_rs_11116403_v1/paper.pdf`
Review positioning extrachromosomal telomeric C-circles as quantifiable ALT-activity
readouts with diagnostic, prognostic and therapeutic relevance, and proposing a unified
model in which replication stress drives C-circle biogenesis through two complementary
replication-coupled pathways shaped by telomeric structures and repair factors.
**Quant.** ALT 310×, C-circles 52×, ATRX 52×, SMARCAL1 40×; ALT in ~10–15% of cancers.
**Relevance:** (b) gives a **continuous** ALT-activity measure to replace the binary label
— the single highest-value design change available, and it converges with Broderick's
two-BIR-pathway model.

### 26. Armendáriz-Castillo et al. (2026) Biology — 28-gene ALT signature
https://doi.org/10.3390/biology15171531
`new/osteo-alt-tcab/ArmendarizCastillo-2026-Deconstructing_the_Alternative_Lengthening_of_Telomeres_Inte__biology15171531/paper.pdf`
Establishes and validates a 28-gene transcriptomic ALT signature, arguing that
ATRX/DAXX mutation screening fails to capture the macro-reprogramming sustaining ALT.
Multivariate Cox and time-dependent ROC show independent prediction of poor overall
survival in TCGA SARC and LGG, outperforming ATRX/DAXX mutation status.
**Relevance:** (b)(d) an independent, published comparator for our Aim 2 19-gene
signature; its explicit argument against genotype-based ALT stratification is the same
argument that undercuts a TCAB1-genotype-proxy split.

### 27. Diplas et al. (2018) Nat Commun — TERTp-wt/IDH-wt GBM genomic landscape
https://doi.org/10.1038/s41467-018-04448-6
`new/osteo-alt-tcab/Diplas-2018-The_genomic_landscape_of_TERT_promoter_wildtype_IDH_wildtype__s41467_018_04448_6/paper.pdf`
~20% of GBMs are TERTp-wt/IDH-wt with no established telomere-maintenance mechanism.
Genomic characterisation of this subgroup identified SMARCAL1 inactivating mutations as
an ATRX-independent route to ALT.
**Quant.** SMARCAL1 129×, ALT 94×, ATRX 40×, C-circles 34×.
**Relevance:** (c) origin of the SMARCAL1-LOF-as-ALT-driver line that Liu 2023 and
Brown 2026 develop; establishes that ALT+ tumours partition into ATRX-driven and
SMARCAL1-driven classes — a partition our 16-line panel could not resolve.

### 28. Mahmoudi et al. (2010) PLoS Biol — WRAP53 is essential for Cajal bodies
https://doi.org/10.1371/journal.pbio.1000521
`new/osteo-alt-tcab/Mahmoudi-2010-WRAP53_Is_Essential_for_Cajal_Body_Formation_and_for_Targeti__journal_pbio_1000521/paper.pdf`
WRAP53 protein maintains Cajal bodies and targets the SMN complex there; Cajal bodies
collapse on knockdown and cannot re-form. WRAP53 associates with coilin, SMN and
importin-β and is required for SMN–coilin and SMN–importin-β complex formation. Reduced
WRAP53–SMN binding in spinal muscular atrophy patients.
**Quant.** WRAP53 mentioned 376× — the densest WRAP53 source in the corpus.
**Relevance:** (d)(c) defines the non-telomeric WRAP53β function. Relevant because
WRAP53 depletion has pleiotropic consequences (splicing, SMN trafficking) that would
confound any dependency contrast built on low WRAP53 expression.

### 29. Mahmoudi et al. (2011) Cell Death Dis — WRAP53 promotes cancer cell survival
https://doi.org/10.1038/cddis.2010.90
`new/osteo-alt-tcab/Mahmoudi-2011-WRAP53_promotes_cancer_cell_survival_and_is_a_potential_targ__cddis_2010_90/paper.pdf`
WRAP53 protein is overexpressed across cancer cell lines of varied origin and its
overexpression promotes transformation. Knockdown triggers massive mitochondrial
apoptosis (Bax/Bak activation, loss of membrane potential, cytochrome c release),
blockable by Bcl-2. Tumour cells are more sensitive to depletion than normal cells.
**Relevance:** (a) WRAP53 is typically **over**expressed and is itself near-essential in
cancer lines. Lines surviving with genuinely low WRAP53 are a selected minority, which
both compresses the usable exposure range and biases the low arm toward compensated
genotypes — a plausible mechanical contributor to the null.

### 30. Niu et al. (2021) Cancer Cell Int — TCAB1 suppression induces senescence
https://doi.org/10.1186/s12935-020-01745-3
`new/osteo-alt-tcab/Niu-2021-Suppression_of_TCAB1_expression_induced_cellular_senescence__s12935_020_01745_3/paper.pdf`
TCAB1 knockdown in A549, H1299 and HeLa drove rapid senescence (SA-β-gal+), via increased
p21 **protein** without change in p21 mRNA — TCAB1 depletion shunts p21 away from
ubiquitin-mediated degradation. Supported by TCGA expression analysis.
**Quant.** TCAB1 157×; three cell lines.
**Relevance:** (a) TCAB1-low cells senesce. Among proliferating screened lines, the
TCAB1-low arm is survivorship-filtered, which attenuates any true dependency contrast.

### 31. Sobieszek & Idzik (2026) Gene Protein Dis — SMARCAL1 in Schimke dysplasia
https://doi.org/10.36922/GPD025480078
`new/osteo-alt-tcab/Sobieszek-2026-Genetic_alterations_in_Schimke_immunoosseous_dysplasia_Patho__gpd025480078/paper.pdf`
Review of SIOD, an autosomal recessive multisystem disorder from biallelic pathogenic
SMARCAL1 variants, covering the annealing-helicase role in fork remodelling and chromatin
dynamics plus therapeutic prospects.
**Quant.** SMARCAL1 132×.
**Relevance:** (c) germline biallelic SMARCAL1 loss is viable in humans — SMARCAL1 is not
constitutively essential, so a near-zero SMARCAL1 gene effect in some lines is expected
rather than informative. ⚠ CONTRADICTS PRESPEC (#1, supporting)

### 32. Sánchez-Morán et al. (2020) Sci Adv — nuclear WRAP53 in neuronal DSB repair
https://doi.org/10.1126/sciadv.abc5702
`new/osteo-alt-tcab/SanchezMoran-2020-Nuclear_WRAP53_promotes_neuronal_survival_and_functional_rec__sciadv_abc5702/paper.pdf`
In oxygen/glucose-deprived neurons and an in vivo stroke model, ROS induce DSBs and
WRAP53 nuclear translocation that activates the DSB repair response. Translocation is
faster in neurons carrying the WRAP53 SNP rs2287499 (c.202C>G).
**Quant.** WRAP53 214×.
**Relevance:** (c) WRAP53β's DDR function is regulated by **localisation**, not abundance
— a third route (after transcript-vs-protein and chaperonin folding) by which total
WRAP53 TPM fails as a functional proxy. ⚠ CONTRADICTS PRESPEC (#5, supporting)

### 33. Gadelha et al. (2022) Curr Issues Mol Biol — WRAP53 review
https://doi.org/10.3390/cimb44110372
`new/osteo-alt-tcab/Gadelha-2022-The_Role_of_WRAP53_in_Cell_Homeostasis_and_Carcinogenesis_On__cimb44110372/paper.pdf`
Review of the WRAP53 locus: antisense transcript WRAP53α stabilising TP53, and protein
WRAP53β (WDR79/TCAB1) scaffolding telomerase localisation, telomere assembly, Cajal body
integrity and DSB repair. Surveys WRAP53 overexpression as a cross-cancer biomarker.
**Relevance:** (d) the single best orientation document for the WRAP53 locus' three
products; useful for writing the Stage 1 limitation about which product WRAP53 TPM
actually measures.

### 34. Mueller et al. (2026) Haematologica — 17p11.2 breakpoints in TP53-altered AML
https://doi.org/10.3324/haematol.2026.301572
`new/osteo-alt-tcab/Mueller-2026-Beyond_17p13_1_recurrent_17p11_2_breakpoints_in_TP53_altered__haematol_2026_301572/paper.pdf`
Short Haematologica report (epub ahead of print) arguing that TP53-altered AML carries
recurrent breakpoints at **17p11.2**, beyond the canonical 17p13.1 locus.
**Quant.** ~2,700 words; 17p13 mentioned 7×.
**Relevance:** (c) the 17p architecture relevant to TP53 inactivation is not confined to
the TP53/WRAP53 overlap. If OS SVs behave similarly, a WRAP53-centred window
(chr17:7,560,000–7,710,000) may miss relevant rearrangements.

### 35. Showman et al. (2024) Cell Rep — centromeric array expansion by BIR in ALT+ cells
https://doi.org/10.1016/j.celrep.2024.113851
`new/osteo-alt-tcab/Showman-2024-Expansion_of_human_centromeric_arrays_in_cells_undergoing_br__j_celrep_2024_113851/paper.pdf`
ddPCR shows α-satellite centromeric arrays expand and contract within ~20 somatic
divisions of an **ALT+ cell line**; array variation among single-cell-derived subclones
ranges ~7% to ~100%, with expansion favoured over contraction.
**Relevance:** (c) BIR activity in ALT+ lines remodels the genome outside telomeres and
generates substantial subclonal heterogeneity — so a single bulk CRISPR screen per line
averages over genuinely heterogeneous BIR states.

### 36. Li et al. (2023) Nucleic Acids Res — MutSβ facilitates MiDAS at fragile sites
https://doi.org/10.1093/nar/gkad1112
`new/osteo-alt-tcab/Li-2023-MutS_protects_common_fragile_sites_by_facilitating_homology__nar_gkad1112/paper.pdf`
MutSβ (MSH2/MSH3) is required for MiDAS and common-fragile-site stability, for the
increased mitotic recombination induced by replication stress or FANCM loss at AT-rich
structure-prone sequences, and for HR when DSB ends carry secondary structures. MSH3 is
synthetic lethal with FANCM.
**Quant.** BIR 60×, RAD52 72×.
**Relevance:** (b) a second FANCM synthetic-lethal partner (with Bythell-Douglas); widens
the FANCM-centred alternative to include MSH3.

### 37. Sung (2025) Cells — ALT as a prognostic paradox
https://doi.org/10.3390/cells14201613
`new/osteo-alt-tcab/Sung-2025-Alternative_Lengthening_of_Telomeres_A_Prognostic_Paradox_in__cells14201613/paper.pdf`
Review: ALT is used by ~10–15% of tumours and is classically tied to genomic instability,
yet its prognostic direction is inconsistent across tumour types.
**Quant.** ALT 228×, ATRX 50×.
**Relevance:** (d) useful framing for a negative-result write-up — ALT-status effects are
not uniformly signed even for survival, so a null dependency contrast within ALT+ is not
anomalous.

### 38. Dueva & Iliakis (2020) NAR Cancer — RPA, a multifunctional ssDNA-binding protein
https://doi.org/10.1093/narcan/zcaa022
`new/osteo-alt-tcab/Dueva-2020-Replication_protein_A_a_multifunctional_protein_with_roles_i__narcan_zcaa022/paper.pdf`
Review of RPA as protector of transient ssDNA and recruitment platform for DDR factors,
with documented interactions including p53, RAD51, ATRIP/ETAA1, BLM, WRN and
**SMARCAL1/HARP**.
**Relevance:** (c) SMARCAL1 function is RPA-dependent, so SMARCAL1 gene effect is
partially a readout of ssDNA load — consistent with Taglialatela's mechanism, and another
reason SMARCAL1 essentiality reflects ALT activity rather than upstream genotype.

### 39. Dudragne et al. (2026) EMBO Rep — transient telomere uncapping drives rearrangements
https://doi.org/10.1038/s44319-026-00717-4
`new/osteo-alt-tcab/Dudragne-2026-Transient_telomere_uncapping_triggers_telomeric_and_subtelom__s44319_026_00717_4/paper.pdf`
In *S. cerevisiae* cdc13-1, transient telomere uncapping rapidly induces extensive
rearrangements despite an intact checkpoint. Two signatures in survivors: subtelomeric
recombination mostly involving Y′ elements, and telomeres elongated to ~10 kb
(~30-fold increase).
**Quant.** RAD52 40×.
**Relevance:** (d) yeast model of the recombinational telomere-maintenance route; useful
mechanistic background, not a mammalian constraint.

### 40. Zhang (2026) BioEssays — kinetic control of APB condensate formation
https://doi.org/10.1002/bies.70187
`new/osteo-alt-tcab/Zhang-2026-Kinetic_Control_of_Condensate_Function_How_the_Formation_Dyn__bies_70187/paper.pdf`
Essay: APBs form by two mechanistically and kinetically distinct pathways that modulate
the strength of ALT activity, which may explain extreme telomere-length heterogeneity in
ALT cells.
**Relevance:** (b) a second independent "ALT runs in two modes" argument (cf. Broderick,
Folini). Converges on the same conclusion: ALT activity is graded and bimodal, not a
label.

---

## Tier 2 — Supporting corpus

### 41. Gadelha et al. (2026) Cancer Rep — WRAP53 downregulated in AML
https://doi.org/10.1002/cnr2.70464
`new/osteo-alt-tcab/Gadelha-2026-WRAP53_is_Downregulated_in_Acute_Myeloid_Leukemia_Patients_a__cnr2_70464/paper.pdf`
qPCR of WRAP53 and hTERT in **110 AML patients** vs controls: WRAP53 significantly
decreased in AML, hTERT unchanged; moderate positive WRAP53–hTERT correlation.
**Relevance:** (c) a real-tumour precedent for low WRAP53 expression as a measurable
phenotype — and for its decoupling from telomerase transcript level.

### 42. Yuan et al. (2014) Oncol Lett — TCAB1 depletion blocks telomerase trafficking
https://doi.org/10.3892/ol.2014.2306
`new/osteo-alt-tcab/Yuan-2014-Telomerase_Cajal_body_protein_1_depletion_inhibits_telomeras__ol_2014_2306/paper.pdf`
TCAB1 knockdown in A549 was antiproliferative, inhibited telomerase trafficking to
telomeres and caused G1 arrest without apoptosis.
**Relevance:** (c) TCAB1 loss is anti-proliferative in telomerase-positive cells; the
expected phenotype differs by TMM, so pooling across a panel blurs it.

### 43. Sun et al. (2016) J Cell Mol Med — WDR79 overexpression in NSCLC
https://doi.org/10.1111/jcmm.12759
`new/osteo-alt-tcab/Sun-2016-Overexpression_of_WDR_79_in_nonsmall_cell_lung_cancer_is_lin__jcmm_12759/paper.pdf`
WDR79 (= TCAB1/WRAP53β) is frequently overexpressed in NSCLC lines and tissues;
knockdown inhibits proliferation.
**Relevance:** (d) part of the WDR79-in-NSCLC cluster (with #44–46) establishing TCAB1 as
an overexpressed pro-proliferative factor — i.e. the opposite exposure direction from ours.

### 44. Sun et al. (2017) Cell Death Dis — WDR79 via USP7–Mdm2–p53
https://doi.org/10.1038/cddis.2017.162
`new/osteo-alt-tcab/Sun-2017-WDR79_promotes_the_proliferation_of_non_small_cell_lung_canc__cddis_2017_162/paper.pdf`
WDR79 promotes NSCLC proliferation through USP7-mediated regulation of the Mdm2–p53 axis.
**Quant.** TCAB1/WDR79 154×.
**Relevance:** (c) couples TCAB1 to p53 stability *in trans*, independent of the 17p13.1
cis overlap — a confounder for any TP53-adjacent interpretation of TCAB1 effects.

### 45. Chen et al. (2018) J Cell Mol Med — WDR79 mediates NSCLC proliferation
https://doi.org/10.1111/jcmm.13580
`new/osteo-alt-tcab/Chen-2018-WDR79_mediates_the_proliferation_of_nonsmall_cell_lung_cance__jcmm_13580/paper.pdf`
Follow-up on the WDR79 overexpression phenotype in NSCLC.
**Relevance:** (d) cluster member; no independent constraint.

### 46. Wang et al. (2017) Sci Rep — EBV-induced TCAB1 upregulation in NPC
https://doi.org/10.1038/s41598-017-03156-3
`new/osteo-alt-tcab/Wang-2017-Epstein_Barr_virus_induced_up_regulation_of_TCAB1_is_involve__s41598_017_03156_3/paper.pdf`
EBV infection upregulates TCAB1 in nasopharyngeal carcinoma, implicating it in the DNA
damage response.
**Quant.** TCAB1 138×.
**Relevance:** (d) TCAB1 abundance is responsive to external (viral) signals — background
for treating TCAB1 expression as a stable line-level property.

### 47. Bonab et al. (2014) Med Oncol — WRAP53 rs2287499 and breast cancer
https://doi.org/10.1007/s12032-014-0168-4
`new/osteo-alt-tcab/Bonab-2014-Single_strand_conformational_polymorphism_analysis_of_a_comm__s12032_014_0168_4/paper.pdf`
SSCP association study of WRAP53 rs2287499 (C/G) with breast cancer risk and prognosis in
an Iranian cohort; WRAP53 framed as a candidate cancer susceptibility gene.
**Relevance:** (d) same SNP as Sánchez-Morán 2020. Germline-risk angle, not dependency.

### 48. Logan et al. (2018) Biol Open — scaRNA2/scaRNA9 dynamics under stress
https://doi.org/10.1242/bio.037101
`new/osteo-alt-tcab/Logan-2018-Altered_dynamics_of_scaRNA2_and_scaRNA9_in_response_to_stres__bio_037101/paper.pdf`
DNA-damaging agents alter the scaRNA9 full-length-to-fragment ratio and scaRNA2 levels;
processed fragments are proposed to form regulatory RNPs modulating snoRNP activity.
**Relevance:** (d) the scaRNA cargo WRAP53β traffics is itself stress-responsive —
peripheral mechanism, no bearing on the dependency contrast.

### 49. Mata-Garrido et al. (2016) Acta Neuropathol Commun — unrepaired DNA in neurons
https://doi.org/10.1186/s40478-016-0312-9
`new/osteo-alt-tcab/Mata_Garrido-2016-Neuronal_accumulation_of_unrepaired_DNA_in_a_novel_specific__s40478_016_0312_9/paper.pdf`
Irradiated rat sensory ganglion neurons accumulate persistent DNA lesions; characterises
the resulting chromatin, DDR-kinetic and transcriptional changes.
**Quant.** WRAP53 18×.
**Relevance:** (d) WRAP53-adjacent neuronal DDR work; off our axis.

### 50. Kotenko & Makovets (2024) EMBO Rep — RPA and PCNA recruit Pif1
https://doi.org/10.1038/s44319-024-00114-9
`new/osteo-alt-tcab/Kotenko-2024-The_functional_significance_of_the_RPA_and_PCNA_dependent_re__s44319_024_00114_9/paper.pdf`
Identifies RPA as the major recruiter of Pif1 in budding yeast, alongside the known
Pif1–PCNA interaction; the two modes act independently during telomerase inhibition.
**Relevance:** (b) mechanism for how PIF1 (the BIR helicase in Shah 2024 and Wu 2025)
gets to ssDNA — supporting detail for a PIF1-centred variant endpoint.

### 51. Lee et al. (2023) PLoS Genet — Pif1 nuclear localisation signal
https://doi.org/10.1371/journal.pgen.1010853
`new/osteo-alt-tcab/Lee-2023-Identification_of_the_nuclear_localization_signal_in_the_Sac__journal_pgen_1010853/paper.pdf`
Identifies a C-terminal NLS (781KKRK784) in *S. cerevisiae* Pif1; pif1-NLSΔ retains
mitochondrial function but is defective in nuclear roles.
**Quant.** BIR 117×.
**Relevance:** (d) yeast Pif1 cell biology; BIR-dense but not translatable to a human
dependency readout.

### 52. Mustafi et al. (2023) J Biol Chem — single-molecule Pif1 translocation
https://doi.org/10.1016/j.jbc.2023.104817
`new/osteo-alt-tcab/Mustafi-2023-Single_molecule_visualization_of_Pif1_helicase_translocation__j_jbc_2023_104817/paper.pdf`
Single-molecule imaging of Pif1 helicase translocation on ssDNA.
**Relevance:** (d) biophysics; background for the PIF1 axis only.

### 53. Touchman et al. (2000) Genome Res — 200-kb sequence around CTNS
https://doi.org/10.1101/gr.10.2.165
`new/osteo-alt-tcab/Touchman-2000-The_Genomic_Region_Encompassing_the_Nephropathic_Cystinosis__gr_10_2_165/paper.pdf`
Complete sequencing of a 200-kb segment at 17p13 encompassing CTNS, with discovery of a
novel gene inside the common cystinosis-causing deletion.
**Quant.** 17p13 mentioned 7×.
**Relevance:** (d) historical 17p13 architecture; predates and does not cover the
TP53/WRAP53 overlap.

### 54. Wang et al. (2022) Mol Oncol — pan-cancer telomerase-component signatures
https://doi.org/10.1002/1878-0261.13324
`new/osteo-alt-tcab/Wang-2022-Genomic_epigenomic_and_transcriptomic_signatures_for_telomer__1878_0261_13324/paper.pdf`
Pan-cancer genomic, epigenomic and transcriptomic analysis of telomerase complex
components.
**Quant.** TCAB1 61×, telomerase/TERT 299×.
**Relevance:** (d) reference for how TCAB1/WRAP53 behaves transcriptionally across TCGA —
a sanity check if the expression proxy is revisited.

### 55. Sung (2026) MedComm — stemness/telomere/metabolism nexus in GBM
https://doi.org/10.1002/mco2.70693
`new/osteo-alt-tcab/Sung-2026-The_Nexus_of_Stemness_Telomere_Stability_and_Metabolic_Repro__mco2_70693/paper.pdf`
Review proposing a "Triadic Nexus" framework linking glioma stem-like cells, telomere
stability and metabolic reprogramming.
**Relevance:** (d) conceptual; glioma-focused.

### 56. Tatar et al. (2026) Mol Biol Rep — redox regulation of telomerase and telomeres
https://doi.org/10.1007/s11033-026-12826-3
`new/osteo-alt-tcab/Tatar-2026-Redox_regulation_of_telomerase_and_telomeres_in_cancer_a_cri__s11033_026_12826_3/paper.pdf`
Critical review splitting the redox–telomere literature into five propositions —
including redox-associated ALT — and appraising each on exposure, model, causal test and
molecular directness, concluding much apparent disagreement is design heterogeneity.
**Relevance:** (d) methodologically interesting as a model for appraising a heterogeneous
literature; substantively off our axis.

### 57. Greiner-Mai et al. (2026) Clin Exp Pediatr — telomere biology disorders in childhood ILD
https://doi.org/10.3345/cep.2026.00290
`new/osteo-alt-tcab/Greiner_Mai-2026-Telomere_biology_disorders_associated_with_childhood_interst__cep_2026_00290/paper.pdf`
Clinical series of telomere biology disorders presenting as childhood interstitial lung
disease.
**Quant.** WRAP53 6×.
**Relevance:** (d) germline telomere-disorder context (cf. Bergstrand 2020).

### 58. Shay & Wright (2019) Nat Rev Genet — telomeres and telomerase: three decades
https://doi.org/10.1038/s41576-019-0099-1
`1 Telomeres.pdf`
Timeline review of the telomere/telomerase field: holoenzyme components, replicative
limits, human telomere-shortening disorders, telomerase inhibition as therapy, telomeric
regulation of gene expression.
**Quant.** telomerase/TERT 204×; ALT 15×.
**Relevance:** (d) general background; almost no ALT or BIR content.

### 59. Ahmed et al. (2024) Cancer Metastasis Rev — MCM10 in cancer
https://doi.org/10.1007/s10555-024-10209-3
`new/osteo-alt-tcab/Ahmed-2024-Multifaceted_role_of_the_DNA_replication_protein_MCM10_in_ma__s10555_024_10209_3/paper.pdf`
Review of MCM10 in replication initiation/elongation/termination and replisome assembly;
overexpression linked to aggressiveness in several cancers; disruption alters replication
timing with genome instability.
**Relevance:** (d) replication-stress background; MCM10 is not in our gene sets.

### 60. Liu et al. (2025) BMC Cancer — osteosarcoma single-cell + bulk atlas
https://doi.org/10.1186/s12885-025-13714-y
`6 single cell + good figures.pdf`
Human OS cell atlas from >110,000 single cells across 17 samples; ML tumour-purity
prediction models from transcriptome; Scissor + enrichment to link cell-intrinsic
features to prognosis; ATF6α loss-of-function tested genetically and pharmacologically;
CellChat for cell–cell communication.
**Relevance:** (d) OS-specific transcriptomic context and figure-design reference. No ALT,
TCAB1 or BIR content.

### 61. Yang et al. (2022) Front Immunol — OS immunogenic-cell-death 5-gene prognostic model
https://doi.org/10.3389/fimmu.2022.994034
`7.pdf`
scRNA-seq + bulk RNA-seq integration yielding a 5-gene immunogenic-cell-death panel as an
OS prognostic model.
**Relevance:** (d) OS prognostic-signature methodology; unrelated to telomere maintenance.

---

## Tier 3 — Off-axis corpus (keyword false positives)

Retrieved by the OpenAlex scan because a query term matched, but not about the
TCAB1/ALT/BIR question. Listed for completeness; none bears on Stage 1.

| Citation | DOI | Path | Why it's here |
|---|---|---|---|
| Garcia et al. (2025) bioRxiv — IDH-mutant astrocytoma cell lines | https://doi.org/10.64898/2025.12.30.696808 | `new/osteo-alt-tcab/Garcia-2025-Multi_omics_characterization_of_IDH_mutant_astrocytoma_deriv__2025_12_30_696808/paper.pdf` | Four new LGG lines; ALT mentioned 23× in passing |
| Qiu et al. (2026) Clin Cancer Res — PanNET cell states | https://doi.org/10.1158/1078-0432.ccr-25-2565 | `new/osteo-alt-tcab/Qiu-2026-Distinct_Malignant_Cell_States_and_Myeloid_Glutamate_Signali__1078_0432_ccr_25_2565/paper.pdf` | snRNA-seq of PanNET; ATRX/DAXX co-occurrence only |
| Avanthay et al. (2025) Endocr Pathol — PanNET molecular subtypes | https://doi.org/10.1007/s12022-025-09889-6 | `new/osteo-alt-tcab/Avanthay-2025-Molecular_Subtypes_of_Pancreatic_Neuroendocrine_Tumors_Mutat__s12022_025_09889_6/paper.pdf` | MEN1/DAXX/ATRX ADM subtypes; ATRX 48× |
| Labreche et al. (2018) Acta Neuropathol — glioma risk SNPs by subtype | https://doi.org/10.1007/s00401-018-1825-z | `new/osteo-alt-tcab/Labreche-2018-Diffuse_gliomas_classified_by_1p_19q_co_deletion_TERT_promot__s00401_018_1825_z/paper.pdf` | 1,659 tumours; TERTp/IDH/1p19q subgroups |
| Lulla et al. (2024) Neurooncol Adv — pediatric HGG phase II trial | https://doi.org/10.1093/noajnl/vdae035 | `new/osteo-alt-tcab/Lulla-2024-Vorinostat_temozolomide_or_bevacizumab_with_irradiation_and__noajnl_vdae035/paper.pdf` | Vorinostat/bevacizumab/TMZ trial; ALT incidental |
| Modrek et al. (2026) Nat Commun — persistent chromatin changes after DNA damage in GBM | https://doi.org/10.1038/s41467-026-77331-4 | `new/osteo-alt-tcab/Modrek-2026-Persistent_chromatin_alterations_and_gene_expression_reprogr__s41467_026_77331_4/paper.pdf` | No ALT/BIR/TCAB1 content |
| Alpizar-Sosa et al. (2025) BMC Genomics — *Leishmania* telomeres | https://doi.org/10.1186/s12864-025-12154-z | `new/osteo-alt-tcab/Alpizar_Sosa-2025-Leishmania_mexicana_telomeres_at_high_resolution_Ku80_TERT_a__s12864_025_12154_z/paper.pdf` | Parasite telomere maintenance; ALT 43× by homology |
| Wang et al. (2017) Nat Genet — DNA-repair genomics stratifies ovarian histotypes | https://doi.org/10.1038/ng.3849 | `new/osteo-alt-tcab/Wang-2017-Genomic_consequences_of_aberrant_DNA_repair_mechanisms_strat__ng_3849/paper.pdf` | "osteosarcoma" 51× in signature naming only |
| Spielmann et al. (2018) Nat Rev Genet — structural variation in the 3D genome | https://doi.org/10.1038/s41576-018-0007-0 | `new/osteo-alt-tcab/Spielmann-2018-Structural_variation_in_the_3D_genome__s41576_018_0007_0/paper.pdf` | SV/enhancer-hijacking review; no locus overlap |
| Weischenfeldt et al. (2016) Nat Genet — pan-cancer SCNA, IRS4/IGF2 | https://doi.org/10.1038/ng.3722 | `new/osteo-alt-tcab/Weischenfeldt-2016-Pan_cancer_analysis_of_somatic_copy_number_alterations_impli__ng_3722/paper.pdf` | Enhancer hijacking; not 17p |
| PCAWG Transcriptome Core Group (2020) Nature — genomic basis for RNA alterations | https://doi.org/10.1038/s41586-020-1970-0 | `new/osteo-alt-tcab/Group-2020-Genomic_basis_for_RNA_alterations_in_cancer__s41586_020_1970_0/paper.pdf` | Pan-cancer RNA alterations resource |
| Fehrmann et al. (2015) Nat Genet — global gene dosage sensitivity | https://doi.org/10.1038/ng.3173 | `new/osteo-alt-tcab/Fehrmann-2015-Gene_expression_analysis_identifies_global_gene_dosage_sensi__ng_3173/paper.pdf` | Dosage-sensitivity method |
| Gatza et al. (2014) Nat Genet — luminal breast cancer proliferation drivers | https://doi.org/10.1038/ng.3073 | `new/osteo-alt-tcab/Gatza-2014-An_integrated_genomics_approach_identifies_drivers_of_prolif__ng_3073/paper.pdf` | Breast integrative genomics |
| Seo et al. (2012) Genome Res — lung adenocarcinoma transcriptome | https://doi.org/10.1101/gr.145144.112 | `new/osteo-alt-tcab/Seo-2012-The_transcriptional_landscape_and_mutational_profile_of_lung__gr_145144_112/paper.pdf` | 200 lung adenocarcinomas |
| Deng et al. (2012) Gut — genomic alterations in gastric cancer | https://doi.org/10.1136/gutjnl-2011-301839 | `new/osteo-alt-tcab/Deng-2012-A_comprehensive_survey_of_genomic_alterations_in_gastric_can__gutjnl_2011_301839/paper.pdf` | Gastric target survey |
| Zhang et al. (2020) Cancers — global DNA hypomethylation in ovarian cancer | https://doi.org/10.3390/cancers12030764 | `new/osteo-alt-tcab/Zhang-2020-Global_DNA_Hypomethylation_in_Epithelial_Ovarian_Cancer_Pass__cancers12030764/paper.pdf` | GDHO and survival in EOC |
| Zheng et al. (2025) Mol Cancer — ALKBH5 in gastric cancer | https://doi.org/10.1186/s12943-024-02223-4 | `new/osteo-alt-tcab/Zheng-2025-ALKBH5_suppresses_gastric_cancer_tumorigenesis_and_metastasi__s12943_024_02223_4/paper.pdf` | m6A/ALKBH5; WRAP53 85× as an m6A target |
| He et al. (2020) Ann Transl Med — TP53 pathway in Wilms' tumour | https://doi.org/10.21037/atm-20-6047 | `new/osteo-alt-tcab/He-2020-Comprehensive_bioinformatics_analysis_of_the_TP53_signaling__atm_20_6047/paper.pdf` | Bioinformatics of TP53 signalling |
| Braicu et al. (2019) Cancers — non-coding RNAs in lung cancer | https://doi.org/10.3390/cancers11050605 | `new/osteo-alt-tcab/Braicu-2019-The_Function_of_Non_Coding_RNAs_in_Lung_Cancer_Tumorigenesis__cancers11050605/paper.pdf` | lncRNA review; WRAP53 13× (antisense transcript) |
| Goyal et al. (2016) Nucleic Acids Res — CRISPR/Cas9 for lncRNAs | https://doi.org/10.1093/nar/gkw883 | `new/osteo-alt-tcab/Goyal-2016-Challenges_of_CRISPR_Cas9_applications_for_long_non_coding_R__nar_gkw883/paper.pdf` | Methodological caveats for lncRNA editing |
| Canudas et al. (2019) Am J Clin Nutr — pistachio consumption and telomere genes | https://doi.org/10.1093/ajcn/nqz048 | `new/osteo-alt-tcab/Canudas-2019-Pistachio_consumption_modulates_DNA_oxidation_and_genes_rela__ajcn_nqz048/paper.pdf` | T2D nutrition RCT; telomere keyword match |
| Munquad & Das (2026) Front Bioinform — deep learning and omics review | https://doi.org/10.3389/fbinf.2026.1811813 | `new/osteo-alt-tcab/Munquad-2026-The_impact_of_deep_learning_and_omics_data_in_transforming_p__fbinf_2026_1811813/paper.pdf` | ML methods review |
| Cunningham (2026) Warwick thesis — particle Monte Carlo cluster analysis | https://doi.org/10.82444/warw.33820447 | `new/osteo-alt-tcab/Cunningham-2026-Particle_Monte_Carlo_methods_for_the_integrative_cluster_ana__warw_33820447/paper.pdf` | Statistical methodology thesis |
| Arjaans et al. (2026) Funct Integr Genomics — adipose/breast-tumour crosstalk | https://doi.org/10.1007/s10142-026-02012-w | `new/osteo-alt-tcab/Arjaans-2026-Exploring_the_transcriptional_crosstalk_between_adipose_tiss__s10142_026_02012_w/paper.pdf` | Breast LRR transcriptomics |

---

## ⚠ Unable to read

No PDF in `Literature/new/` is corrupt or scan-only — all 80 extracted clean full text
(minimum 2,710 words, Mueller 2026, a short Haematologica letter). The unreadable items
are papers the scan **failed to retrieve**.

**Five folders exist with supplementary material but no main text** (`paper.pdf` absent;
all Elsevier/paywalled). SI files are present where listed:

| Citation | DOI | Folder | Note |
|---|---|---|---|
| Lee et al. (2026) Mol Cells — divergent proteomic landscapes and functional requirements in **type I vs type II ALT** | https://doi.org/10.1016/j.mocell.2026.100405 | `new/osteo-alt-tcab/Lee-2026-Divergent_proteomic_landscapes_and_functional_requirements_f__j_mocell_2026_100405/` | ⚠ unable to read — 5 SI files only. **Highest-value gap**: directly tests whether ALT has mechanistically distinct subtypes |
| Audrey et al. (2024) Cell Rep — RAD52-dependent MiDAS in Cyclin E1-overexpressing cells | https://doi.org/10.1016/j.celrep.2024.114116 | `new/osteo-alt-tcab/Audrey-2024-RAD52_dependent_mitotic_DNA_synthesis_is_required_for_genome__j_celrep_2024_114116/` | ⚠ unable to read — 2 SI files only |
| Krishnamoorthy et al. (2021) Mol Cell — RADX confines fork reversal to stalled forks | https://doi.org/10.1016/j.molcel.2021.05.014 | `new/osteo-alt-tcab/Krishnamoorthy-2021-RADX_prevents_genome_instability_by_confining_replication_fo__j_molcel_2021_05_014/` | ⚠ unable to read — 4 SI files only |
| Carvalho et al. (2014) Am J Hum Genet — dosage changes of a **17p13.1** segment | https://doi.org/10.1016/j.ajhg.2014.10.006 | `new/osteo-alt-tcab/Carvalho-2014-Dosage_Changes_of_a_Segment_at_17p13_1_Lead_to_Intellectual__j_ajhg_2014_10_006/` | ⚠ unable to read — 2 SI files only |
| Li et al. (2013) Cell — eQTL analyses of breast cancer risk loci | https://doi.org/10.1016/j.cell.2012.12.034 | `new/osteo-alt-tcab/Li-2013-Integrative_eQTL_Based_Analyses_Reveal_the_Biology_of_Breast__j_cell_2012_12_034/` | ⚠ unable to read — 4 SI files only |

**Thirty-seven DOIs with empty `*__scihub` folders** (retrieval failed entirely; no text,
no SI). Marked ★ where the title indicates genuinely high relevance — these are worth
chasing by hand through institutional access:

- ★ https://doi.org/10.1016/j.cancergen.2026.06.001 — Beyond ATRX: the heterogeneous basis of ALT **in osteosarcoma** (2026)
- ★ https://doi.org/10.1080/15384101.2025.2609646 — TCAB1: a key promoter of tumour growth through telomere maintenance and senescence evasion (2026)
- ★ https://doi.org/10.1128/mcb.00379-12 — Telomerase recruitment requires both TCAB1 and Cajal bodies independently (2012)
- ★ https://doi.org/10.7554/elife.112326 — Break-induced replication drives large-scale genomic amplifications in cancer cells (2026)
- ★ https://doi.org/10.1073/pnas.1510750112 — SMARCAL1 maintains telomere integrity during DNA replication (2015)
- ★ https://doi.org/10.1016/j.dnarep.2017.06.015 — The role of SMARCAL1 in replication fork stability and telomere maintenance (2017)
- ★ https://doi.org/10.1080/19491034.2016.1179413 — SMARCAL1 and telomeres: replicating the troublesome ends (2016)
- ★ https://doi.org/10.1016/j.molcel.2016.11.027 — Wrap53, a natural p53 antisense transcript required for p53 induction upon DNA damage (2016)
- ★ https://doi.org/10.4161/cc.8.15.9223 — Wrap53, a novel regulator of p53 (2009)
- ★ https://doi.org/10.1093/neuonc/noae016 — Telomere maintenance mechanisms and therapeutic vulnerabilities in malignant gliomas (2024)
- ★ https://doi.org/10.1136/jcp-2023-209005 — ALT: mechanism and the pathogenesis of cancer (2023)
- ★ https://doi.org/10.1200/jco-26-00742 — Prognostic impact of telomere maintenance mechanisms in neuroblastoma (2026)
- ★ https://doi.org/10.1128/mcb.00151-10 — Specificity and stoichiometry of subunit interactions in the human telomerase holoenzyme (2010)
- ★ https://doi.org/10.17863/cam.129190 — Structural and biochemical characterisation of the human telomerase catalytic cycle (2025)
- https://doi.org/10.7554/elife.112326.1 · https://doi.org/10.7554/elife.112326.1.sa1 · https://doi.org/10.7554/elife.112326.1.sa2 · https://doi.org/10.7554/elife.112326.1.sa4 — versions and peer-review artifacts of the eLife BIR paper above, **not four separate works**
- https://doi.org/10.1080/15476286.2016.1211224 — processing elements and interactors of SMN/coilin/pseudogene scaRNAs (2016)
- https://doi.org/10.1016/j.critrevonc.2025.105017 — telomere maintenance meets cancer metabolism (2025)
- https://doi.org/10.1016/j.semcancer.2010.03.001 — non-LTR retrotransposons and genome stability (2010)
- https://doi.org/10.1002/gcc.20579 — molecular characterisation of 17p loss in B-CLL (2008)
- https://doi.org/10.1182/blood-2025-3506 — 17p cytogenetic changes in TP53-mutated AML (2025)
- https://doi.org/10.1182/blood-2025-3941 — modelling del(17p) in multiple myeloma (2025)
- https://doi.org/10.1038/modpathol.3800090 — cytogenetics of 43 aneurysmal bone cysts, 17p aberrations (2004)
- https://doi.org/10.1309/jqmr-323g-71y9-m7mb — cytogenetic findings in mantle cell lymphoma (2001)
- https://doi.org/10.1002/gcc.10158 — cytogenetics of cervical cancer cell lines (2003)
- https://doi.org/10.1016/j.ajpath.2017.04.022 — breast cancer molecular stratification (2017)
- https://doi.org/10.1158/1078-0432.ccr-19-3231 — single-cell characterisation of the gastric tumour microenvironment (2020)
- https://doi.org/10.1530/erc-24-0305 — management of small pancreatic neuroendocrine tumours (2025)
- https://doi.org/10.1086/341963 — pericentric inversion differing between chimpanzee and human (2002)
- https://doi.org/10.1073/pnas.050397097 — MLL-GAS7 translocation during topoisomerase-inhibitor chemotherapy (2000)
- https://doi.org/10.2174/0113816128421211260214133902 — multi-target therapy and prognostic modelling in liver cancer (2026)

---

## Gaps

Topics that matter to this hypothesis space and are **not covered** by any readable PDF
in `Literature/`.

1. **ALT subtypes (type I vs type II) — the largest gap.** Three readable papers argue ALT
   runs in two modes (Broderick 2023 EXD2; Zhang 2026 APB kinetics; Folini 2026 two
   C-circle pathways), but the one paper that directly profiles the subtypes
   (Lee 2026, Mol Cells, `10.1016/j.mocell.2026.100405`) could not be retrieved. If the
   BIR module splits across subtypes, that alone could produce our near-zero composite.
   **Get this paper first.**
2. **ALT heterogeneity in osteosarcoma specifically.** `10.1016/j.cancergen.2026.06.001`
   ("Beyond ATRX: the heterogeneous basis of ALT in osteosarcoma", 2026) is unretrieved.
   Every readable ALT-heterogeneity paper we have is glioma (Diplas, Liu, Brown,
   Brosnan-Cashman) or pan-cancer. We have **no** OS-specific ALT-subtype source.
3. **Axis G — expression vs structural-variant proxies — has no corpus at all.** This was
   already flagged in `openalex_scan/axis_g_expression_vs_sv_proxies.md`: the
   pre-registered query returned 6 off-topic works, and six broadenings drifted into
   maize, rice, yeast and Parkinson's. Nothing in `Literature/` validates (or invalidates)
   using gene expression as a proxy for a structural variant at a specific locus. The
   closest substitutes are indirect: Jiang 2025 and Bergstrand 2020 (protein-level
   decoupling) and Sánchez-Morán 2020 (localisation-level decoupling). **This remains the
   single weakest-evidenced assumption in the Stage 1 design and it is unresolved.**
4. **TCAB1 and BIR in the same paper: zero coverage.** Across all 87 readable PDFs, no
   paper mentions both TCAB1/WRAP53 and break-induced replication substantively. Keegan
   2026 has 0 BIR/POLD3/RAD52 mentions; Bennett 2026 has 0 TCAB1/WRAP53 mentions. The
   Stage 1 hypothesis joined two literatures that have not yet been joined by anyone.
   That is a reason the prediction was speculative, and also the reason a confirmed
   result would have been novel.
5. **DepMap strandedness / antisense-overlap quantification.** Nothing in the corpus
   addresses whether an unstranded protocol can resolve the 1,419-nt TP53/WRAP53
   head-to-head overlap. The notebook 02 branch decision ("direct, use WRAP53 TPM")
   still rests on our own empirical diagnostics, not on a citable source.
6. **Osteosarcoma CRISPR dependency methodology.** Taglialatela 2026 cites Dharia 2021
   (pediatric cancer dependencies) and Harris & Hawkins 2022; neither is in
   `Literature/`. If a Stage 2 dependency analysis is written up, those are the
   methodological precedents to hold.
7. **No negative/contradictory report on SMARCAL1 as an ALT dependency.** Every readable
   source agrees SMARCAL1 matters in ALT — they disagree only on *sign*. We have no
   published failure-to-replicate to calibrate our own null against.

---

## Other material in `Literature/`

`openalex_scan/` holds the seven-axis OpenAlex scan this corpus was drawn from: eight
axis files (a–g plus `axis_z_recombination_paradox_term.md`), `README.md` (synthesis) and
`next_hypotheses.md` (3–5 candidate Stage 2 hypotheses). Those files carry the search
provenance — queries, date, result counts, and the 132-DOI result list. This index
summarises the PDFs; the scan files record how they were found.

`new/osteo-alt-tcab/` also contains retrieval metadata: `_report.json`,
`_summary.csv`, `_scihub_manifest.csv` (42 second-pass DOIs, 37 retrieval failures).
