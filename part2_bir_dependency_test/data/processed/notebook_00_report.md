# Notebook 00 — Setup and reference verification

Verifies the environment and the frozen `references/` snapshot. No analysis.

## Environment

- Python: 3.12.13

| Library | Version |
|---|---|
| pandas | 2.3.3 |
| numpy | 2.5.4 |
| scipy | 1.18.1 |
| sklearn | 1.9.1 |
| matplotlib | 3.11.2 |

## Reference snapshot

All 8 files present, read-only (mode 444), sha256 matching `references/MANIFEST.sha256`.

Note: git does not track the write bit, so `chmod -w` does not survive `git clone`. This notebook restores mode 444 when it finds a reference file writable; the sha256 check is the hard integrity guarantee. Whether anything needed restoring on a given run is recorded in `provenance.jsonl`, not here, so this report stays byte-identical across runs.

| File | sha256 match | Mode | Bytes |
|---|---|---|---|
| aim2_signature.py | yes | 444 | 28320 |
| aim2_utils.py | yes | 444 | 5380 |
| feasibility_decision_memo.md | yes | 444 | 8293 |
| label_rules.py | yes | 444 | 16708 |
| label_spec_v1.0.md | yes | 444 | 5718 |
| manifest_builder.py | yes | 444 | 35068 |
| master_manifest.csv | yes | 444 | 64029 |
| notebook_13_report.md | yes | 444 | 4432 |

## Source provenance

```
Source project : /Users/rinatrizvanov/Desktop/Projects/Project 1 - Comp Bio for Osteosarcoma/tcab1-alt-osteo pipeline
Source git HEAD: 5af4d0f31ed0d777c2daca32c4cfd1e2d8d02ea8
Source git log : 5af4d0f 2026-06-12 notebooks 13
Copied         : 2026-10-11T22:54:32Z
```

`master_manifest.csv` is the frozen label input. Stage 1 never re-derives ALT labels.
