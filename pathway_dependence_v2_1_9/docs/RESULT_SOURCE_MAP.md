# Result sources and mathematical objects

| Claim or figure | Object | Builder | Outputs |
| --- | --- | --- | --- |
| Main CFD budget, tolerance and edge-cost profiles; Figures 7.1/7.2 | Functional inside one common partition, then optimize across partitions | run_real_cfd_common_partition_profiles.py | 37_REAL_CFD_COMMON_PARTITION_PROFILES.csv/json |
| Main CFD compact summary and shared threshold interval | Summary of the preceding common-partition profiles | build_authoritative_cfd_summary.py | 42_REAL_CFD_COMMON_PARTITION_SUMMARY.csv/json |
| Weight reoptimization and precommitment | Fixed-target qualifying cut family | analyze_weight_stability.py | weight_stability_results.json |
| Auxiliary union budget and threshold analyses | Functional of the union of all affordable assignments | Four builders in auxiliary/robust_union | 28/29/33/35_ROBUST_UNION files |

Only the first five CSV/JSON files are regenerated and compared by default.
`--include-robust-union` adds eight separately grouped auxiliary files.
The distinction concerns quantifier order, regardless of whether particular
computed numbers happen to be close. All raw vectors and figure assets are
unchanged. Historical precomputation records remain under provenance.
