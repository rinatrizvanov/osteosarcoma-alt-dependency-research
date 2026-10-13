Prespec freeze: 5f2f039

PRIMARY  rho(PRIMPOL, SMARCAL1) = +0.0000  p_perm = 1.0000  n = 16  floor = 0.4943
         aligned = -0.0000  (predicted_sign = -1)
VERDICT  null  (rule: clean-null)
         Clean null: |aligned delta| = 0.000 < 0.2 (p = 1). No effect of interpretable size in either direction.

Exposure gate (evaluated before the outcome was read):
   gene  n     sd  n_zero  guarded_gap  shapiro_p  sd_ok  no_zeros_ok  gap_ok  shapiro_ok  gate_passed
PRIMPOL 16 0.7610       0       0.4536     0.1174   True         True    True        True         True
 ZRANB3 16 0.4269       0       0.2059     0.9419  False         True    True        True        False
    BLM 16 0.8947       0       0.2196     0.0059   True         True    True       False        False

Analyses:
                analysis exposure                  outcome  n     rho  p_perm  floor_p05  reaches_floor  aligned
                 primary  PRIMPOL                 SMARCAL1 16  0.0000  1.0000     0.4943          False  -0.0000
     S1 negative control   ZRANB3                 SMARCAL1 16  0.1647  0.5403     0.4943          False  -0.1647
     S2 screen amplitude  PRIMPOL  genome-wide mean effect 16 -0.2794  0.2967     0.4943          False   0.2794
        S2b standardised  PRIMPOL SMARCAL1 (z within line) 16  0.0118  0.9703     0.4943          False  -0.0118
     S3 ALT- specificity  PRIMPOL                 SMARCAL1  9 -0.1333  0.7501     0.6833          False   0.1333
S4 BLM (UNINTERPRETABLE)      BLM                 SMARCAL1 16  0.1882  0.4813     0.4943          False  -0.1882

LOO rho range: [-0.0929, +0.1857], no sign flip: False

Support conditions:
  PASS exposure_variance_gate: sd=0.761, shapiro p=0.117, guarded gap=0.454
  FAIL negative_control_null: ZRANB3 rho=+0.1647, p=0.5403
  FAIL not_screen_confounded: screen-amplitude rho=-0.2794 vs primary +0.0000
  FAIL loo_no_sign_flip: range [-0.0929, +0.1857]

S4 BLM is UNINTERPRETABLE by prespec §2/§5: it fails the exposure gate
(Shapiro p = 0.0059); its distribution in these lines is 2 low
outliers plus 14 lines inside a 1.496-unit band. A large value is not evidence.
