# Notebook 10 — ALTitude × DepMap ALT Merge
_Generated: 2026-05-20T14:34:35Z_

Integrates ALTitude curated ALT calls into the manifest for RNA-backed DepMap OS cell lines. First mutating notebook since the audit freeze.

## Priority hierarchy
| Tier | Source | alt_status written | source_confidence |
|---|---|---|---|
| 1 (kept) | Existing manifest gold | unchanged | unchanged |
| 2 | ALTitude experimental | alt_gold_positive / alt_negative_high_conf | peer_reviewed |
| 2★ | HOS special case (PMID:33344901) | alt_negative_high_conf | peer_reviewed |
| 3 | ALTitude prediction-only | alt_silver | computational |
| — | No call | unchanged (alt_unknown) | unchanged |

## Invariant checks
- I1 (no downgrade):             **PASS**
- I2 (no peer_reviewed/unknown): **PASS**
- I3 (no gold→silver):           **PASS**

## Ground-truth checks
- ACH-000364 (U2OS): `alt_gold_positive` [PASS]
- ACH-000410 (SAOS2): `alt_gold_positive` [PASS]
- ACH-000748 (SJSA1): `alt_negative_high_conf` [PASS]
- ACH-000359 (MG63): `alt_negative_high_conf` [PASS]
- ACH-000613 (HOS, special case): `alt_negative_high_conf` [NaN pred_class; C-circle negative; PMID:33344901]

## Merge actions
- Rows written:  **0**
- skip_non_os: 1
- skip_same: 31

## alt_status distribution — overall
| alt_status | pre | post |
|---|---|---|
| alt_gold_positive | 60 | 60 |
| alt_negative_high_conf | 19 | 19 |
| alt_silver | 11 | 11 |
| alt_unknown | 25 | 25 |

## alt_status distribution — DepMap_26Q1 cohort
| alt_status | pre | post |
|---|---|---|
| alt_gold_positive | 12 | 12 |
| alt_negative_high_conf | 8 | 8 |
| alt_silver | 11 | 11 |
| alt_unknown | 4 | 4 |

## include_aim2 flag counts — before and after
| flag | pre | post |
|---|---|---|
| include_aim2_strict | 79 | 79 |
| include_aim2_pragmatic | 90 | 90 |
| include_aim2_broad | 90 | 90 |

## alt_silver direction split (DepMap_26Q1)
- positive (ALT+): **6**
- negative (ALT-): **5**

| sample_id | alt_call_direction |
|---|---|
| ACH-001275 | negative |
| ACH-001813 | negative |
| ACH-001817 | negative |
| ACH-002433 | positive |
| ACH-002779 | negative |
| ACH-003171 | positive |
| ACH-003176 | negative |
| ACH-003178 | positive |
| ACH-003185 | positive |
| ACH-003573 | positive |
| ACH-003587 | positive |

## alt_status × alt_call_direction (all cohorts)
| alt_status | positive | negative | unknown |
|---|---|---|---|
| alt_gold_positive | 60 | 0 | 0 |
| alt_negative_high_conf | 0 | 19 | 0 |
| alt_silver | 6 | 5 | 0 |
| alt_unknown | 0 | 0 | 25 |

## Idempotence
Second-pass writes (must be 0): see cell below.

## Next steps
1. Re-run notebooks 05 (Pass A) and 04 (Pass B) to update audit tables.
2. No DE or modeling code has been written.
