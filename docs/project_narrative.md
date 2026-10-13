# The whole arc, in plain English

Written for a scientifically literate reader who has not been following this project. No
prior knowledge of the specific genes is assumed. Roughly three pages.

---

## 1. The problem

Cells cannot divide forever. Every division shortens the telomeres that cap each
chromosome, and when they get too short the cell stops. Cancers have to solve this, and
most do it by switching telomerase back on. But a substantial minority — including a large
share of osteosarcomas, a bone cancer that mostly affects teenagers — use a completely
different route called **ALT**, for *alternative lengthening of telomeres*. ALT cells
rebuild their telomeres by copying them from each other using DNA-repair machinery.

This matters therapeutically. ALT depends on repair enzymes that normal cells do not lean
on nearly as hard, so in principle you could kill ALT cancer cells by blocking those
enzymes while leaving healthy tissue alone. To look for such targets you need two things:
a reliable way to tell which tumours or cell lines are using ALT, and a way to find out
what those cells depend on.

Both exist in public form. The **DepMap** project has systematically knocked out nearly
every gene in over a thousand human cancer cell lines and measured how much each knockout
hurt. That gives a dependency score per gene per cell line. What DepMap does *not* ship is
a reliable ALT label. Hence Part 1.

## 2. Part 1 — teaching a model to recognise ALT

The first project built a classifier. Starting from manually curated ALT status for a set
of osteosarcoma cell lines, it trained a regularised regression model on the expression
levels of 19 genes already implicated in telomere biology and DNA repair, and asked how
well those 19 genes predict ALT status.

**It worked, by the usual measure.** On the primary, most stringently labelled set of 20
cell lines the model reached an AUROC of 0.948 — meaning if you pick one ALT-positive and
one ALT-negative line at random, the model ranks them correctly about 95% of the time. A
permutation test (reshuffling the labels a thousand times) put that well beyond chance. On
a looser label definition with 31 lines it reached 0.829.

**And then a check made the result more complicated.** Rather than stopping at the headline
number, the project asked whether the *specific genes* the model relied on were stable:
refit the model on 200 bootstrap resamples and see which genes keep getting picked. Only
three of the nineteen survived that test at N = 20. Fourteen were selected less than half
the time. Across both cohorts only two genes — `FANCM` and `TERT` — were consistently
retained.

That is a real and somewhat deflating finding, and it is reported in full. The
interpretation is that at twenty samples the model can tell the classes apart, but the
particular combination of genes it uses to do so is largely a product of which twenty
samples you happened to have. The discrimination is probably real; the gene list is not a
signature anyone should carry to a new cohort without re-validating it. Encouragingly, the
one gene that was always retained, `TERT`, carries a negative weight — lower telomerase,
more likely ALT — which is exactly the biology you would hope to see.

## 3. The hypothesis that motivated Part 2

A 2026 paper described something striking in osteosarcoma genomes. At a spot on chromosome
17 called 17p13.1, two genes sit head to head and overlap: `TP53`, the famous tumour
suppressor, and `WRAP53`, which encodes a protein called TCAB1. Structural rearrangements
in that region can knock out both at once. TCAB1 is needed to deliver telomerase to
telomeres, and the same protein has a second job helping organise DNA double-strand-break
repair.

That suggests a clean prediction. A cell that has lost TCAB1 has lost its telomerase route,
so if it is alive and dividing it must be using ALT. And it has partially lost the
conventional repair scaffold. So it should be *unusually* dependent on the specific
non-conventional repair machinery that ALT runs on — a process called **break-induced
replication**, or BIR.

Concretely: among ALT-positive osteosarcoma lines, the ones with low TCAB1 should be more
dependent on BIR genes (`POLD3`, `RAD52`, `SMARCAL1`) relative to conventional repair genes
(`BRCA1`, `BRCA2`, `RAD51`). That is a testable, directional claim, and all the data to test
it already existed.

## 4. Part 2 — the test, pre-registered

Because the prediction was appealing, the sample size small, and the number of ways to slice
the analysis large, the whole plan was written down and committed to version control
*before* any dependency value was read: which lines, which split, which statistic, what
would count as support, what would count as failure, and what would count as uninformative.
Code enforced it — one module was the only thing allowed to load the outcome data, and the
analysis notebooks assert that module is not loaded.

**The result did not confirm the prediction.** The effect size was +0.281 where the
hypothesis predicted a negative number, with p = 0.382. In other words, not only was there
no detectable effect, the point estimate leaned the opposite way. That reversal held up
across every robustness check: all sixteen leave-one-out refits, both alternative split
definitions, and the continuous version of the analysis all agreed on the wrong sign. Low
statistical power explains a failure to detect something; it does not explain a consistent
reversal.

## 5. What came next, and why it is the most useful part

Rather than stop, the project ran four further pre-registered iterations testing
alternatives drawn from the recent literature — stratifying on `ATRX` instead, using `BLM`
or `PRIMPOL` expression as the exposure, and asking whether the gene sets behave as
coordinated modules at all. Every iteration was frozen before its outcome was read, and
critically, the significance threshold was corrected across *all* of them together. A
result only counted if it survived at q < 0.1 across the whole campaign, not p < 0.05 within
its own iteration. Nothing came close: the final corrected values were 0.420, 0.420 and
1.000.

The fourth iteration is where the project learned something. It asked whether the three
"BIR genes" actually behave as a module in these cells — and they do not. Worse, two of
them cannot carry the analysis at all. `POLD3` is essential in every single line (knock it
out and any cell dies, ALT or not), and `RAD52` is essential in none. Only `SMARCAL1` varies
in an informative way. When the composite score was decomposed, **68% of its variance came
from the comparison arm** rather than the BIR arm, and `RAD52` contributed slightly
negatively.

So the original endpoint was mostly measuring its own control group, and subtracting a
genuinely coordinated module (the HR genes do cohere) from a set of genes that is not a
module at all. The experiment could not have detected the effect it was designed to detect.
No amount of statistical care downstream of that choice would have revealed it — which is
precisely why the check was worth running.

One more control matters: a known protein complex (BTR) was tested on the same sixteen cell
lines and came out strongly coherent, p = 0.0001. That rules out the easy excuse. Sixteen
lines is enough to detect this kind of structure when it is there, so the nulls carry
information rather than reflecting a cohort too small to say anything.

## 6. The honest conclusion

**The hypothesis has not been refuted.** It has not been properly tested. The exposure
variable — TCAB1 expression — turned out to be a poor stand-in for the structural damage the
original paper described, because only two of sixteen lines had anything like the relevant
genomic architecture and the rest sat in one undifferentiated clump. And a 2025 paper showed
TCAB1 protein levels are controlled after the gene is transcribed, which means the
transcript was arguably never the right measurement in the first place.

What genuinely remains open is a real observation the project could not explain. Dependency
on `SMARCAL1` varies a lot among ALT-positive lines — eight of sixteen depend on it
substantially, eight do not. Four separate explanations were tested and none accounted for
it. An independent group reported the same unexplained heterogeneity, so this is a real open
problem rather than a local artefact.

The campaign stopped because every remaining literature-backed idea needed something the
data does not contain: a direct measurement of ALT activity, structural-variant calls,
protein levels, or a subtype label that would have to be invented. Inventing one was exactly
the mistake that started this, so the project stopped instead of manufacturing a fifth
hypothesis to stay alive.

One correction belongs here, because it is part of the record. The first version of that
closing judgement declared one remaining hypothesis untestable on a flawed argument. A later
check against the project's own positive control showed the argument was wrong, and the
hypothesis is genuinely untested rather than closed. The erratum is published in the
repository rather than quietly folded in.
