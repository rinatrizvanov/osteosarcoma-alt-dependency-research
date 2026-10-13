# Iteration 03 outcome — H4: do POLD3, RAD52 and SMARCAL1 behave as a module?

**Mechanical verdict: `null`** (T = −0.182, p = 0.229, n = 16, `rule = clean-null`).
**Verdict of record: `not-estimable` — design-insensitive by the frozen §6 gate**, which
fired and which forbids reading the primary as evidence against module coherence.

Prespec frozen `b24a8ed`; statistic and detection floor committed earlier at `cf73508` in a
module with no path to the outcome. Notebook idempotent.

---

## 1. The primary, and the gate that voids its interpretation

| | T | p | floor at this (n, k) | reaches floor |
|---|---|---|---|---|
| **Primary** — POLD3/RAD52/SMARCAL1 | **−0.1824** | **0.2285** | 0.2892 | no |
| Pol delta gate, raw | −0.1637 | 0.2825 | 0.2892 | no |
| Pol delta gate, standardised | −0.1961 | 0.1916 | 0.2892 | no |

Pairwise: POLD3–RAD52 **−0.268**, POLD3–SMARCAL1 **+0.024**, RAD52–SMARCAL1 **−0.303**.
Not one pair is positive. The observed T sits at the **12.5th percentile** of the
mean/SD-matched random-triplet null (q95 = +0.279), so the triplet is *less* coherent than
arbitrary genes with the same marginals. Four of five support conditions failed; the fifth
(no leave-one-out sign flip) passed, with T over 16 refits in [−0.248, −0.152].

**The §6 gate failed.** Pol delta reached −0.164 raw and −0.196 standardised against a
pre-registered threshold of +0.40 on both scales. The frozen text is categorical: *"If it
fails either, the primary is reported as 'design insensitive at this n', the verdict is
recorded as a feasibility failure, and the result is **not** read as evidence against the
BIR module."* That rule is applied as written.

## 2. But the other controls passed decisively, and that changes the reading

| Control module | k | T | p | reaches floor |
|---|---|---|---|---|
| **BTR** — BLM/TOP3A/RMI1/RMI2 | 4 | **+0.5255** | **0.0001** | **yes** |
| **HR** — BRCA1/BRCA2/RAD51 | 3 | **+0.4137** | **0.0039** | **yes** |
| Batenburg BIR — MUS81/RAD52/POLD3 | 3 | −0.2735 | 0.0665 | no |
| Extended fork restart (k = 5) | 5 | −0.0538 | 0.5196 | no |

BTR pairwise: BLM–TOP3A +0.682, RMI1–RMI2 +0.800, TOP3A–RMI2 +0.456, BLM–RMI2 +0.444,
TOP3A–RMI1 +0.388, BLM–RMI1 +0.382.

**This matters more than the primary.** The §6 gate was built on a worry that co-essentiality
cannot detect coherence among genes without dynamic range. BTR refutes that worry
empirically: it contains **TOP3A**, which is near-saturated (16/16 lines below −0.5), and
**BLM and RMI2**, which are essential in *no* line (0/16 below −0.5) — and it is the most
coherent module measured here, at p = 0.0001. The HR module passes too. So the design is
**demonstrably sensitive at n = 16**, across exactly the marginal regimes the gate was
designed to worry about.

So Pol delta's failure is not a property of the method. It is a property of **POLD3 and its
complex partners**: POLD1 (median −2.250), POLD2 (−1.058) and POLD3 (−1.666) are so deeply
essential that their cross-line variation carries no recoverable co-dependency signal, where
TOP3A at −1.078 still does.

**The conservative reading is the one of record**, because that is what was frozen. But the
informative reading is that the gate's failure is itself evidence **about POLD3** rather than
about the design: POLD3 cannot act as a module member in this cohort. Both readings are
stated; only the first governs the verdict.

## 3. What is established regardless of the gate

Two findings are deterministic given the 16 lines and do not depend on the primary test,
the gate, or any p-value.

**3.1 Two of the three Stage 1 BIR genes have no usable dynamic range** (pre-registered in
prespec §2.1, before the freeze):

| Gene | median | SD | lines < −0.5 | regime |
|---|---|---|---|---|
| POLD3 | −1.666 | 0.319 | **16/16** | ceiling: essential in every line |
| RAD52 | −0.010 | 0.173 | **0/16** | floor: essential in no line |
| SMARCAL1 | −0.462 | 0.395 | 8/16 | the only member with mid-range spread |

**3.2 The Stage 1 composite's variance came mostly from its HR arm.**

| Gene | arm | variance contribution | share |
|---|---|---|---|
| BRCA2 | HR | 0.03315 | **25.8%** |
| SMARCAL1 | BIR | 0.03007 | **23.4%** |
| BRCA1 | HR | 0.02809 | 21.9% |
| RAD51 | HR | 0.02624 | 20.5% |
| POLD3 | BIR | 0.01430 | 11.1% |
| RAD52 | BIR | −0.00351 | **−2.7%** |

The three HR genes supply **68.2%** of the composite's variance and the three BIR genes
**31.8%**, of which SMARCAL1 alone is 23.4% and **RAD52 contributes negatively**. The decomposition
is exact (contributions sum to the composite variance; asserted in code).

So Stage 1's "BIR-minus-HR" endpoint was, in variance terms, largely a measurement of its
*comparator* arm, with one informative BIR gene. This is the quantitative form of its own R8
observation, where dropping SMARCAL1 moved delta from +0.281 to +0.156.

**3.3 And the comparator arm was the coherent one.** HR passes the coherence test
(T = +0.414, p = 0.0039) while the BIR triplet does not. Stage 1's composite subtracted a
real module from a non-module.

## 4. Secondaries

- **S1b — Broderick's prediction is corroborated directionally.** Broderick 2023
  (https://doi.org/10.1038/s41467-023-38029-z) places SMARCAL1 in the RAD52-*independent*
  BIR branch, predicting rho(RAD52, SMARCAL1) <= 0. Observed **−0.303**, so the prediction
  holds in sign. It is **not** significant: at k = 2 the floor is 0.494 and p = 0.255, and
  the bootstrap CI is [−0.737, +0.214]. Suggestive, not established.
- **S4 — the literature's own BIR module is anti-coherent too.** Batenburg's
  MUS81–RAD52–POLD3 gives T = −0.274 (p = 0.067), driven by MUS81–RAD52 = **−0.512**. So
  substituting MUS81 for SMARCAL1, as the literature's actual module composition would have
  it, does not rescue coherence — it slightly worsens it.
- **S2b — standardisation changes nothing.** T = −0.188 (p = 0.213) on within-line z-scored
  effects, against −0.182 raw. The result is not an artefact of line-level screen scale.
- **S3 — the ALT− arm is the opposite sign.** T = **+0.378** (p = 0.058, n = 9, floor 0.384).
  Just short of its floor, and descriptive only at n = 9, but it runs *positive* where the
  ALT+ arm runs negative. Worth flagging and not over-reading.
- **S5 — subtype secondary: no defensible marker exists.** This is rule (c)'s
  explicit-statement branch, pre-registered in prespec §7 S5 and restated here as required.
  Type I/II ALT is defined in mouse ESCs via the ~7 kb mTALT template, and Lee et al. 2026's
  own supplementary material (https://doi.org/10.1016/j.mocell.2026.100405) reports that no
  compact human genomic segment resembling it was identified in GRCh38. The cohort carries no
  ALT-activity measure of any kind (`alt_assay_evidence` is `unknown` for all 35 rows;
  `telomerase_status` is constant within ALT+), and no candidate marker gene is bimodal over
  the 34 RNA-backed lines (largest guard-respecting interior gap below 1 SD for all 16
  tested). Any cut would be a researcher choice producing an unfalsifiable label. **No proxy
  was substituted.**
- **S6 — vacuous as implemented, and reported as such.** The batch-robustness refit was
  meant to drop commercially-deposited lines, but the depositor regex matched **0 of 16**
  lines, so the "refit" silently reproduced the primary. It is therefore **not** an
  independent check and contributes no information. Recorded as a failed check rather than
  counted as a passed one — the same treatment Stage 1 gave its vacuous R5.

## 5. Verdict, and why it is recorded conservatively

`src.verdict.evaluate_iteration` returns **`null`**: |aligned delta| = 0.182 < 0.200 with
p = 0.229. Four of five support conditions failed.

The frozen §6 gate overrides the *interpretation* of that null, so the **verdict of record is
`not-estimable` (design-insensitive)**, and the primary is not read as evidence against the
coherence of the BIR triplet.

This is deliberately the unfavourable bookkeeping. The data in §2 argue the gate was
mis-specified and that the null is in fact interpretable — but the rule was frozen before the
result was seen, it fired, and amending it now in light of BTR's passing would be exactly the
post-hoc rule change the loop forbids. Stage 1 set the precedent: its §10 clean-null ceiling
was left exactly as frozen even though it forced an awkward `indeterminate`, with the methods
lesson recorded rather than retrofitted.

**Methods lesson, recorded not retrofitted.** The gate was specified as a *single* designated
positive control, chosen because it contains POLD3. That was the wrong design: it conflates
"the method cannot detect coherence here" with "this particular gene carries no coherence
signal", which is precisely the distinction the gate existed to draw. A correct gate names
**several** control complexes spanning marginal regimes and passes if **any** of them clears
the threshold. Under that rule the design would have been declared sensitive — BTR at
p = 0.0001 — and the primary null interpretable. This is a defect in the iteration's own
pre-specification, found by the data, and it is written here rather than fixed silently.

**It is not fixed by re-running.** The primary T and p would be identical under a corrected
gate; only the interpretive label would change. Re-running the same test under a new
iteration number to obtain a more favourable label would be laundering, not science, so it is
not done. The corrected gate design belongs in any future prespec that uses this statistic.

## 6. Campaign accounting

**Tests run: 28** (the prespec enumerated 29; S6 is excluded because, as §4 records, it
reproduced the primary rather than running an independent refit). Enumerated: primary 1; Pol
delta raw and standardised 2; BTR 1; HR 1; S1 pairwise 3; S4 1; S4b 1; S2b 1; S3 1; LOO 16.
The matched-triplet calibration and the variance decomposition are not hypothesis tests.

**Campaign BH across all iteration primaries** (`src.verdict.campaign_fdr`): iteration 01
p = 0.280 and iteration 03 p = 0.2285 both give **q = 0.280**. Neither survives q < 0.1.
**No campaign positive.**

## 7. What this iteration does and does not license

**It does not** establish that POLD3, RAD52 and SMARCAL1 are mutually independent — the
frozen gate forbids that conclusion from this iteration, and co-essentiality across 16 lines
speaks to cross-line variance, not biochemistry. A pathway whose members sit at a ceiling or
a floor can be perfectly coupled mechanistically and show no correlation.

**It does** establish, deterministically and independently of the gate, that the Stage 1
endpoint was poorly constructed: two of its three BIR genes have no usable dynamic range in
this cohort, its variance was 68% comparator-arm, RAD52 contributed negatively, and the
comparator arm it subtracted was the only coherent module in the set. That is a real finding
about this project's prior work and it belongs in the preprint's discussion whatever else
happens.

**It also establishes** that the method works here — BTR and HR are both strongly coherent at
n = 16 — so "16 osteosarcoma lines is too few for co-dependency analysis" is not a defensible
excuse for any future null on this cohort.

## 8. Streak and stop conditions

Third consecutive non-positive iteration (01 `counter_signal`, 02 `not-estimable`, 03
`not-estimable`/`null`). Stop condition (b) fires at five.

**A new mechanistic angle did surface**, so the "no new mechanistic angle" clause does not
advance: the HR arm, not the BIR arm, is the coherent module in these lines, and BRCA2 is its
single largest variance contributor — which meets Lee 2024's finding that BRCA2 loss itself
drives ALT/BIR. A BRCA2/HR-centred endpoint is a literature-grounded candidate that did not
exist before this iteration.
