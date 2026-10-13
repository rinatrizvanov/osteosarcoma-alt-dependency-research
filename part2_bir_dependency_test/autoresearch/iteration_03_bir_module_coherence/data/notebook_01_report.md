Prespec freeze: b24a8ed | statistic commit: cf73508

PRIMARY  T = -0.1824  p_perm = 0.2285  (n = 16, floor 0.2892)
VERDICT  null  (rule: clean-null)
         Clean null: |aligned delta| = 0.182 < 0.2 (p = 0.2285). No effect of interpretable size in either direction.

Pairwise rho:
  POLD3-RAD52: -0.2676
  POLD3-SMARCAL1: +0.0235
  RAD52-SMARCAL1: -0.3029

Gate (Pol delta): raw -0.1637, standardised -0.1961, threshold 0.4 -> FAILED
Matched triplets: observed at the 12.5th percentile, q95 +0.2794
S2b standardised T: -0.1882
S3 ALT- T: +0.3778 (n = 9, floor 0.384)
S1b Broderick rho(RAD52,SMARCAL1) = -0.3029 (supported)
LOO T range: [-0.2476, -0.1524], no sign flip: True

Support conditions:
  FAIL all_pairs_positive
  FAIL min_pair_magnitude
  FAIL exceeds_matched_triplets
  FAIL assay_sensitive
  PASS loo_no_sign_flip

Module table:
                           module                           genes  k  n_lines       T  p_perm  floor_p05  reaches_floor  min_pair_rho  max_pair_rho  has_ties
            BIR triplet (Stage 1)            POLD3+RAD52+SMARCAL1  3       16 -0.1824  0.2285     0.2892          False       -0.3029        0.0235     False
            Pol delta (gate, raw)               POLD1+POLD2+POLD3  3       16 -0.1637  0.2825     0.2892          False       -0.2235       -0.0735     False
   Pol delta (gate, standardised)               POLD1+POLD2+POLD3  3       16 -0.1961  0.1916     0.2892          False       -0.2471       -0.1324     False
        BTR (BLM/TOP3A/RMI1/RMI2)             BLM+TOP3A+RMI1+RMI2  4       16  0.5255  0.0001     0.2029           True        0.3824        0.8000     False
   HR module (Stage 1 subtrahend)               BRCA1+BRCA2+RAD51  3       16  0.4137  0.0039     0.2892           True        0.2618        0.6882     False
Batenburg BIR (MUS81/RAD52/POLD3)               MUS81+RAD52+POLD3  3       16 -0.2735  0.0665     0.2892          False       -0.5118       -0.0412     False
      Extended fork restart (k=5) MUS81+RAD52+POLD3+SMARCAL1+PIF1  5       16 -0.0538  0.5196     0.1556          False       -0.5118        0.4559     False

Variance decomposition of the Stage 1 composite:
    gene arm  weight  own_sd  variance_contribution  share_of_total
   POLD3 BIR  0.3333  0.3191                0.01430          0.1114
   RAD52 BIR  0.3333  0.1727               -0.00351         -0.0274
SMARCAL1 BIR  0.3333  0.3950                0.03007          0.2343
   BRCA1  HR -0.3333  0.4395                0.02809          0.2189
   BRCA2  HR -0.3333  0.3749                0.03315          0.2583
   RAD51  HR -0.3333  0.3168                0.02624          0.2045
