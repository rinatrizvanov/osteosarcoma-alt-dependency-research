# Data Provenance

DepMap release: **26Q1**. Files copied read-only from the Aim 2 project
and sha256-verified byte-identical to the originals.

| File | Release | sha256 | Bytes |
|---|---|---|---|
| `data/raw/depmap/CRISPRGeneEffect.csv` | 26Q1 | `e610a4cefb13a82b...` | 440646050 |
| `data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 26Q1 | `0377be80c525fde9...` | 305007605 |
| `data/raw/depmap/Model.csv` | 26Q1 | `ea4e0b2a3bc806f8...` | 697455 |
| `data/raw/depmap/OmicsCNGeneWGS.csv` | 26Q1 | `d2fe65b516802fb2...` | 409005582 |
| `data/raw/depmap/OmicsCNSegmentsWGS.csv` | 26Q1 | `690e861b3a35db89...` | 63718976 |

## Re-fetch policy

If a file must be replaced, fetch it from the DepMap **26Q1 archive
specifically**, never from "latest". The frozen ALT labels in
`references/master_manifest.csv` were derived against this release; a different
release changes which lines are screened and would silently break the label join.

## Frozen label input

| File | sha256 | Role |
|---|---|---|
| `references/master_manifest.csv` | `d51c5f180391ffbc...` | ALT labels (`alt_call_direction`, `alt_status`); never re-derived |
