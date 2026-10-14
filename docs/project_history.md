# How this repository was consolidated

This repository was assembled on **2026-10-13** from two separate working repositories. Git
history was **not** carried across: this repo begins with a fresh initial commit. That was a
deliberate choice, and this file is the record of what it means.

## The two parent repositories

| Part here | Parent repository (local) | Commits | Final state |
|---|---|---|---|
| `part1_alt_classification/` | `tcab1-alt-osteo pipeline` | 13 | HEAD `5af4d0f` "notebooks 13", 2026-06-12 |
| `part2_bir_dependency_test/` | `tcab1-bir-dependency` | 30 | HEAD `f372c40`, tags `stage1-v1.0-null`, `autoresearch-complete-v1.0` |

Both parents remain intact on local disk and are unmodified by this consolidation. **They
are the authoritative record of development order**, which matters here more than usual: the
pre-registration claims in Part 2 rest on commit ordering — prespec committed before
analysis code, analysis code committed before execution — and that evidence lives in the
parent repository's `git log`, not in this one.

If the ordering ever needs to be demonstrated to a reviewer, it must be shown from
`tcab1-bir-dependency`, or that repository must be published alongside this one. The freeze
commit hashes are cited inside the individual prespec and outcome-memo files, so the claims
are at least checkable against the parent.

## Why history was flattened

Three reasons, in order of weight:

1. **Clean paths.** Preserving history across a directory move without rewriting commits
   (`git subtree`) leaves historical commits pointing at the old top-level paths, so
   `git log -- part2_bir_dependency_test/...` does not reach back past the merge. Rewriting
   them properly needs `git-filter-repo`, which was not installed, and which would have
   changed every commit hash — including the two tagged states.
2. **Tag integrity.** `stage1-v1.0-null` and `autoresearch-complete-v1.0` keep pointing at
   their original commits in the parent repository. A rewrite would have orphaned them.
3. **Clean slate for release.** A flat start made it straightforward to guarantee that no
   excluded file exists anywhere in this repository's history (below), which matters for a
   repository that may later be made public.

## What was deliberately excluded

Content, not just working-tree state — none of this exists in this repository's history:

- **Raw DepMap data** (~1.2 GB) and the parent projects' local caches (~3 GB total). Not
  redistributed. The five required files, their release and their exact SHA256 checksums are
  specified in the top-level `README.md` under "How to reproduce".
- **Publisher PDFs and supplementary material** (~1.2 GB of literature). Only
  `Literature/INDEX.md`, the OpenAlex scan notes and two small retrieval manifests are kept
  as the provenance record.
- **`_scihub_manifest.csv`** — a 37-row log of failed retrieval attempts. It contained no
  copyrighted text, but the filename documents paywall circumvention and this repository may
  become public. The same 37 DOIs are listed in `Literature/INDEX.md` under "Unable to read",
  so nothing is lost.
- **`docs/.Rhistory`** from Part 1 — stray shell history, no content value.
- **Run-varying provenance logs** (`provenance.jsonl`), which change on every execution.
- **Preprint scaffolding.** A `preprint/` directory (outline, figure plan, abstract draft,
  bibliography skeleton) was drafted in the parent repository before the decision was taken
  not to publish. It is removed here rather than kept with a disclaimer, which would have
  added noise rather than clarity. **No preprint is planned**; the result does not warrant
  formal publication, and that was a deliberate decision rather than a default. The
  scaffolding remains in the parent repository's history for anyone who wants it. Note that
  `part2_bir_dependency_test/docs/preprint_outline.md` is a *different* and older document —
  the Stage 1 closeout outline — and is retained as part of the original record.

## What was added here and did not exist in either parent

- `README.md`, `LICENSE` (MIT), `CITATION.cff`
- `docs/project_narrative.md`, `docs/methodology_summary.md`, `docs/PREREGISTRATION_INDEX.md`,
  this file
- `part1_alt_classification/README.md` — Part 1 had no README of its own
- A parent-context banner prepended to Part 2's existing README

## Intentional duplication, do not "fix"

Eight of Part 1's files appear twice: once in place, and once vendored into
`part2_bir_dependency_test/references/` as frozen read-only copies pinned by
`references/MANIFEST.sha256`. Part 2's entire design depends on consuming
`master_manifest.csv` as an immutable input that it never re-derives, and the checksums are
how that is enforced. Deduplicating the copies would break the pins and void the guarantee.
The vendored copies also record their own provenance — the source project path, the source
git HEAD (`5af4d0f`), the copy timestamp and the original filenames.

## Content carried across, for the record

Part 1 contributed 78 tracked files plus 7 documentation files that were untracked in the
parent at consolidation time (slides, a one-pager, a presentation, and the results summary
that bridges to Part 2); these were included here and so enter version control for the first
time. Part 2 contributed 157 files after the exclusions above. With the 9 new top-level and
`docs/` files, the repository tracks **251 files, 5.5 MB**.

## Public release — 2026-10-14

The repository was audited before being made public. Four classes of change were made; no
analysis code, result, statistic or pre-specified decision rule was touched.

1. **Personal contact address removed.** A personal email appeared in eight OpenAlex query
   URLs recorded in `part2_bir_dependency_test/Literature/openalex_scan/` (as the `mailto`
   politeness parameter the API asks for) and on the closing slide of Part 1's deck. It is
   now `<CONTACT_EMAIL>`, matching the already-redacted `<API_KEY>` in the same URLs. The
   OpenAlex API keys were never committed in cleartext.

2. **Third-party name redacted.** An academic advisor was shown Part 1 in June 2026 and is
   named in the parent repositories as the recipient of the presentation, as the holder of a
   PDX cohort the project considered requesting, and in internal notes about how and when to
   approach them. They have not reviewed or approved this public release and did not
   consent to being named in it, so those references now read as role descriptions
   ("an academic advisor", "a collaborator-held PDX cohort", "the collaborating lab").
   **Citations to their published work are untouched** — Cox, Maréchal & Flynn (2016) and
   Carson & Flynn (2023) are cited by name throughout `Literature/INDEX.md` and the
   iteration scans, exactly as any other source. This redaction edits two frozen documents
   (`part2_bir_dependency_test/docs/analytical_prespec.md` §on falsification, and
   `autoresearch/PROGRAM.md` stop condition (d) and scope list). Neither edit changes a
   threshold, an endpoint, a stop rule or a decision; both change only who the
   data-access counterparty is called. The unredacted originals remain in the two private
   parent repositories, which are the authoritative record of pre-registration ordering.

3. **Slide binaries regenerated.** `part1_alt_classification/docs/project_slides.pptx` and
   `project_slides.pdf` embedded both redacted strings, so they were rebuilt from
   `docs/make_slides.py` (python-pptx) and re-exported to PDF rather than edited in place.

4. **Two README claims tightened.** The "does not redistribute DepMap" bullet now states
   precisely what derived DepMap-sourced values the small output tables do contain and under
   which licence, and the file count in this section was corrected from 245 to the actual
   251.

Deliberately **not** changed: absolute local paths such as
`/Users/<user>/Desktop/Projects/...` that appear inside executed notebook output cells, two
generated reports and the header of `references/MANIFEST.sha256`. Those are execution
records and provenance headers; rewriting them would misrepresent what was run. They expose
no credential and no information beyond a directory name.
