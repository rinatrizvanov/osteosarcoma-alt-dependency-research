# Post-mortem: when is a hypothesis wrong, and when is the test simply unable to see it?

**Filed 2026-10-13 at the user's instruction, for the record. Contains no new analysis.**
Every number below is quoted from a committed artifact in this repository; nothing was
computed to write this file. The campaign is closed (`FINAL.md`, stop condition (c),
tag `autoresearch-complete-v1.0`).

The question in the title is the one this project kept getting wrong, in both directions.
It is worth a document of its own because the campaign produced a clean worked example of
each error, and because the second one was mine, made after four iterations of being
careful about exactly this.

---

## 1. The distinction, stated once

A pre-registered test returns a number. That number can fail to support a hypothesis for
three different reasons, and they carry completely different consequences:

| | Meaning | What follows |
|---|---|---|
| **The hypothesis is wrong** | The effect is absent in this population | Evidence against the claim, bounded by the detection floor |
| **The test could not see it** | The design lacks power, range, or a valid measurement | No evidence either way. The hypothesis is *untested* |
| **The construct is broken** | The endpoint does not measure the thing named | No evidence either way, and the result is actively misleading if reported as a null |

The campaign's own vocabulary already encoded this: `src/verdict.py` separates `null`,
`refuted`, `counter_signal` and `inconclusive_underpowered`, and `NULL_SPACE.md` separates
`tested-and-null`, `not-estimable` and `skipped-missing-data`. Having the vocabulary did
not stop either error below. The taxonomy is necessary and not sufficient; what it needs
is a habit of asking which column a result belongs in **before** the label is chosen.

---

## 2. Error type A — calling a broken construct a null

**Where it happened: Stage 1.** The primary endpoint was
`mean(POLD3, RAD52, SMARCAL1) − mean(BRCA1, BRCA2, RAD51)`, reported as δ = +0.281,
p = 0.382, a non-confirmation. Iteration 03 established, deterministically, that this was
not a null of the hypothesis. It was a null of a construct that could not have carried it:

- POLD3 is essential in **16/16** lines (median −1.666), RAD52 in **0/16** (median −0.010).
  Two of the three BIR genes had no dependency variation to contribute.
- The HR arm supplied **68.2%** of the composite's variance; RAD52 contributed **−2.7%**.
  The decomposition is exact (`iteration_03/data/h4_variance_decomposition.csv`).
- The HR arm was a coherent module (T = +0.414, p = 0.0039); the BIR triplet was not
  (T = −0.182, no pair positive). **Stage 1 subtracted a module from a non-module.**

Nothing about the pre-registration was violated. The plan was frozen, the criteria were
evaluated by code, the verdict was not shaded. The endpoint was simply named after a
biological module that does not exist in this cohort, and no amount of procedural rigour
downstream of that choice could detect it.

**The transferable rule: check marginal regimes before composing any multi-gene score.**
That is now `NULL_SPACE.md` binding fact #2, and it is cheap — it is a table of medians and
counts, computable before any hypothesis is written.

---

## 3. Error type B — calling an untested hypothesis not-estimable

**Where it happened: my own screen that closed the campaign**, 2026-10-12, corrected the
same day by visible erratum (`a55e6c0`).

I closed Batenburg's cross-arm antagonism hypothesis on the following argument: LIG3 is
essential in 0/16 lines and PRIMPOL in 0/16, so they carry "no dependency signal", so a
cross-set statistic against them would have expectation ≈ 0 **whether or not the claim is
true** — uninformative by construction, therefore not-estimable.

The campaign's own positive control refutes this. BTR reached T = +0.526, p = 0.0001 while
containing two members in exactly that regime:

```
BLM   x TOP3A  rho +0.6824  p 0.0036    BLM essential in 0/16 lines
RMI1  x RMI2   rho +0.8000  p 0.0002    RMI2 essential in 0/16 lines
BLM   x RMI2   rho +0.4441  p 0.0848    BOTH essential in 0/16 lines
```

A gene that no line depends on can still have a graded profile tracking its obligate
complex partners at ρ +0.68. **Sub-threshold variation is recoverable biology here.**

The mechanism of the error is worth naming precisely, because it is subtle and it was
invited by the project's own documentation. Iteration 03's §7.3 class is labelled
"non-essential throughout" and glossed as spread that "is not dependency signal." That
gloss is correct for the use it was written for — whether a line *depends* on the gene, as
needed for a dependency *endpoint*. I carried it to a different use — whether the gene's
across-line *ranking* is informative, as needed for a correlational *exposure* — where it
does not hold. The label travelled; its scope condition did not.

**The transferable rule: a feasibility class is valid only for the use it was derived for.**
When reusing one, restate the use. "This gene is not a usable endpoint" and "this gene is
not a usable correlate" are different claims requiring different evidence.

---

## 4. The asymmetry that makes type B worse

Type A inflates the apparent evidence base: a broken construct reported as a null looks
like a result and enters the literature as one. Type B deflates it, and does so invisibly.

A `tested-and-null` entry publishes its numbers and invites contradiction. A
`not-estimable` entry publishes a *judgement* that the question cannot be asked here, and
nothing in the record prompts anyone to re-examine it. `NULL_SPACE.md` exists so that "no
future session repeats work recorded here" — which is exactly what makes a wrong
not-estimable entry self-sealing. Had the user not asked whether anything remained
explorable, that entry would have stood, and the preprint would have carried the sentence
"nothing further is testable on this cohort," which is false.

**The transferable rule: a not-estimable verdict needs a positive demonstration, not an
argument.** Iteration 02's kill is the correct pattern — it showed all four candidate ATRX
rules and the group sizes each produced (1, 2, 4, 1 against `MIN_CELL_SIZE = 4`), so a
reader can check it. My row-4 kill was an inference from a classification, with no
demonstration that the statistic in question actually fails. That is the difference between
the two, and it is checkable in advance by asking: *what would I show someone who doubted
this?*

---

## 5. What the campaign did get right about the distinction

Recorded so the lessons are not read as uniformly negative:

1. **Iteration 02 refused to run.** Killed at the feasibility gate before any outcome was
   read, with all four candidate exposure rules and their group sizes shown. Later
   corroborated by Bennett's published pan-cancer null (d = −0.248, P = 0.5, n = 20 v 14) —
   the right call on independent evidence, not just on our group sizes.
2. **Iteration 03 applied a frozen rule that fired against the evidence.** The §6 gate's
   single designated control failed while BTR and HR passed decisively, showing the design
   *was* sensitive. The rule was applied anyway, the defect recorded as found-by-data, and
   the identical test deliberately **not** re-run under a corrected label — T and p would
   be unchanged and only the label would differ.
3. **Iteration 04 refused to bank vacuous failures.** Three of four support conditions
   failed arithmetically because they were written relative to a primary that came back at
   exactly zero; the negative control's own p was 0.540, behaving exactly as intended. The
   memo says so rather than reporting "3 of 4 checks failed."
4. **Iteration 04 pre-committed to not over-reading its own null.** The frozen prespec
   stated that a PRIMPOL null would not refute Taglialatela's model — floor |ρ| ≥ 0.500,
   transcript rather than repriming activity, and the mechanism's stronger arm (BLM)
   untestable here. Written before the outcome, against interest.
5. **The detection floor was made genuinely pre-registerable.** The coherence statistic
   depends on values only through within-gene ranks, so its null depends only on (n, k) —
   proven, tested, and therefore immune to being moved by any gene-effect value.
6. **Co-essentiality was calibrated, not assumed.** BTR at p = 0.0001 on these very 16
   lines means "n = 16 is too small" was removed as an available excuse for any null in the
   campaign. This is what separates these nulls from uninformative ones — and, as §3 shows,
   it is also the fact that refuted my own kill. A calibration control earns its keep twice.

---

## 6. The checklist this leaves behind

For any future dependency-screen hypothesis, in order, before a prespec is frozen:

1. **Marginals first.** Median, SD, and n-below-threshold for every gene in the design.
   Costs one table. Catches error type A.
2. **Name the use of every feasibility class you cite.** Endpoint or correlate. They have
   different requirements. Catches error type B.
3. **Distinguish "no dependency" from "no information."** A gene essential in no line may
   still rank lines informatively. Do not infer the second from the first.
4. **A not-estimable verdict owes a demonstration.** Show the configurations that fail and
   the numbers they produce, as iteration 02 did. An argument from classification is not
   enough.
5. **Write support conditions as absolute tests**, never relative to the primary's
   magnitude — a primary of zero makes relative conditions undefined, which is the most
   likely outcome in an underpowered campaign.
6. **Segment-validate any gene-level copy number** before using it as an exposure. Third
   defect of this kind in this project.
7. **Ask what you would show a doubter.** If the answer is "my reasoning," the verdict is
   not ready.

---

## 7. Standing status

Closed. `FINAL.md` is the campaign's conclusion of record and its erratum is part of it:
the screen is **twelve of thirteen**, with Batenburg's row 4 standing as `untested` on data
already held, and the Pol ζ axis (Wu 2023) untouched and estimable. Neither is a resumption
of this loop; both would need a new mandate and their own prespec. No analytical work
continues on this project.
