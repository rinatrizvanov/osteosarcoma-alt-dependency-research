# Notebook 05 — Pass A Report
Run: 2026-05-30T21:05:14Z

Recommended definition: strict

## Counts
_universe_ = `include_aim1_primary_{set}` true; _primary_contrast_ = `disrupted + intact` (go/no-go uses `min(disrupted, intact)`); _uncertain_ = in universe but neither arm (e.g. Broad RNA-proxy with `tcab1_genomic_status` unknown).

- strict: universe=40, primary_contrast=40, disrupted=8, intact=32, uncertain=0, aim2_pos=60, aim2_neg=19
- pragmatic: universe=40, primary_contrast=40, disrupted=8, intact=32, uncertain=0, aim2_pos=66, aim2_neg=24
- broad: universe=57, primary_contrast=41, disrupted=9, intact=32, uncertain=16, aim2_pos=66, aim2_neg=24

## Status: PIPELINE PAUSED — awaiting human review and label_spec_v1.0 freeze