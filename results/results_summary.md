# Paper 2 results (v2, ch4_c2h6; main transport ZISCALE 0.80, range 0.8-1)
Aircraft-implied ZISCALE (05c, lower bound): 0.78

## T1 main result (ZISCALE 0.80, CH4+C2H6)
| prior | region | prior_t_hr | posterior_t_hr | ratio_post_prior | prior_fossil_share | posterior_fossil_share | fossil_t_hr |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EPA GHGI 2020 | Paper 1 box (urban) | 2.8 | 3.6 [2.6-5.1] | 1.29 | 0.55 | 0.43 [0.20-0.64] | 1.5 |
| GRA2PES v1.1 | Paper 1 box (urban) | 1.7 | 5.2 [3.8-7.0] | 3.06 | 0.70 | 0.90 [0.83-0.95] | 4.6 |
| GRA2PES v2.0b | Paper 1 box (urban) | 16.5 | 6.6 [4.8-9.2] | 0.40 | 0.23 | 0.39 [0.27-0.52] | 2.5 |
| EPA GHGI 2020 | Observed domain | 25.9 | 35.3 [27.9-43.5] | 1.36 | 0.39 | 0.30 [0.22-0.40] | 10.5 |
| GRA2PES v1.1 | Observed domain | 34.6 | 34.4 [27.5-43.3] | 1.00 | 0.38 | 0.41 [0.29-0.55] | 14.0 |
| GRA2PES v2.0b | Observed domain | 59.0 | 32.5 [26.6-40.0] | 0.55 | 0.42 | 0.37 [0.30-0.45] | 12.0 |
| EPA GHGI 2020 | DJ Basin | 18.8 | 26.9 [20.9-33.5] | 1.43 | 0.40 | 0.29 [0.21-0.40] | 7.9 |
| GRA2PES v1.1 | DJ Basin | 29.1 | 24.7 [19.3-31.9] | 0.85 | 0.36 | 0.29 [0.18-0.44] | 7.1 |
| GRA2PES v2.0b | DJ Basin | 34.6 | 20.8 [16.8-25.8] | 0.60 | 0.53 | 0.39 [0.31-0.50] | 8.2 |

## T2 transport sweep (values EPA / v1.1 / v2)
| ziscale | box_EPA_v11_v2 | obs_box | djb | fossil_box | fossil_obs_box | fossil_djb | rho_fit | urban_signal_share |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.37 | 1.9 / 1.9 / 3.5 | 21.8 / 20.8 / 22.1 | 17.1 / 16.4 / 15.3 | 0.44 / 0.82 / 0.33 | 0.38 / 0.35 / 0.39 | 0.38 / 0.29 / 0.43 | 0.52 / 0.49 / 0.52 | 67 / 53 / 70 |
| 0.5 | 3.6 / 3.5 / 5.6 | 33.4 / 30.7 / 33.4 | 25.3 / 23.0 / 21.7 | 0.42 / 0.84 / 0.40 | 0.32 / 0.29 / 0.33 | 0.32 / 0.18 / 0.34 | 0.53 / 0.51 / 0.53 | 78 / 60 / 81 |
| 0.8 | 3.6 / 5.2 / 6.6 | 35.3 / 34.4 / 32.5 | 26.9 / 24.7 / 20.8 | 0.43 / 0.90 / 0.39 | 0.30 / 0.41 / 0.37 | 0.29 / 0.29 / 0.39 | 0.60 / 0.56 / 0.59 | 72 / 61 / 72 |
| 1 | 5.1 / 5.5 / 7.1 | 31.2 / 34.4 / 31.8 | 21.7 / 24.2 / 19.6 | 0.20 / 0.90 / 0.39 | 0.26 / 0.41 / 0.39 | 0.28 / 0.28 / 0.42 | 0.57 / 0.51 / 0.54 | 69 / 56 / 65 |
| 1.2 | 2.9 / 2.7 / 6.0 | 26.0 / 27.9 / 26.1 | 19.6 / 21.8 / 16.7 | 0.27 / 0.78 / 0.36 | 0.29 / 0.24 / 0.41 | 0.30 / 0.16 / 0.45 | 0.45 / 0.39 / 0.41 | 51 / 38 / 54 |

## T3 tracer sets at the main transport
| prior | tracers | box_t_hr | fossil_share | DOFS | info_gain_bits | dIG_vs_ch4 | sd_fossil_share | regimes_well_partial_not |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EPA GHGI 2020 | ch4 | 4.1 [2.8-5.9] | 0.61 [0.24-0.81] | 3.2 | 5.5 | 0.0 | 0.177 | 1/2/6 |
| EPA GHGI 2020 | ch4_c2h6 | 3.6 [2.6-5.1] | 0.43 [0.20-0.64] | 4.3 | 8.2 | 2.7 | 0.133 | 3/2/4 |
| GRA2PES v1.1 | ch4 | 6.7 [4.9-8.8] | 0.95 [0.90-0.97] | 3.8 | 6.6 | 0.0 | 0.021 | 2/3/0 |
| GRA2PES v1.1 | ch4_c2h6 | 5.2 [3.8-7.0] | 0.90 [0.83-0.95] | 4.2 | 7.7 | 1.1 | 0.036 | 3/2/0 |
| GRA2PES v1.1 | ch4_c2h6_co | 5.6 [4.3-7.2] | 0.88 [0.82-0.93] | 4.2 | 8.2 | 1.6 | 0.034 | 3/2/0 |
| GRA2PES v2.0b | ch4 | 6.8 [4.7-9.6] | 0.40 [0.18-0.60] | 5.3 | 8.7 | 0.0 | 0.129 | 3/5/1 |
| GRA2PES v2.0b | ch4_c2h6 | 6.6 [4.8-9.2] | 0.39 [0.27-0.52] | 7.0 | 12.6 | 3.9 | 0.077 | 5/4/0 |
| GRA2PES v2.0b | ch4_c2h6_co | 6.1 [4.6-8.3] | 0.45 [0.32-0.59] | 7.0 | 12.9 | 4.3 | 0.081 | 5/4/0 |

## T4 uncertainty budget: totals (log sd; 1 sigma)
| quantity | posterior_1sigma_log | across_priors_log_sd | across_transport_range_log_sd | across_transport_full_sweep_log_sd | across_background_log_sd | drift_over_base | dominant |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Paper 1 box total | 0.197 | 0.301 | 0.061 | 0.365 | 0.132 | 0.72 | prior |
| Observed-domain total | 0.135 | 0.042 | 0.016 | 0.198 | 0.079 | 0.84 | posterior |
| DJ Basin total | 0.143 | 0.131 | 0.043 | 0.168 | 0.071 | 0.88 | posterior |

## T4 uncertainty budget: fossil shares
| quantity | posterior_1sigma | across_priors_range | across_transport_range_sd | across_background_sd |
| --- | --- | --- | --- | --- |
| Paper 1 box fossil share | 0.076 | 0.39-0.90 | 0.006 | 0.016 |
| Observed-domain fossil share | 0.054 | 0.30-0.41 | 0.018 | NA |
| DJ Basin fossil share | 0.059 | 0.29-0.39 | 0.007 | 0.017 |

## T5 headline ensemble (priors x transport range x background base/drift)
| region | n_members | total_median | total_range | total_90pct_members | fossil_share_median | fossil_share_range | fossil_share_excl_v11 | biogenic_share_median |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Paper 1 box (urban) | 12 | 4.9 | 2.6-7.1 | 2.6-9.9 | 0.41 | 0.20-0.90 | 0.39 (0.20-0.43) | 0.59 |
| Observed domain | 12 | 30.6 | 26.3-35.3 | 24.5-43.5 | 0.38 | 0.26-0.41 | 0.33 (0.26-0.39) | 0.62 |
| DJ Basin | 12 | 21.9 | 17.3-26.9 | 15.8-33.5 | 0.29 | 0.28-0.42 | 0.34 (0.28-0.42) | 0.71 |

## T6 prior dependence by region (main transport, CH4+C2H6)
| region | prior_spread_log | post_spread_log | total_prior_dependence | post_spread_over_post_sigma | fossil_prior_range | fossil_post_range | fossil_prior_dependence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Paper 1 box (urban) | 1.197 | 0.301 | 0.25 | 1.56 | 0.23-0.70 | 0.39-0.90 | 1.09 |
| Observed domain | 0.417 | 0.042 | 0.1 | 0.32 | 0.38-0.42 | 0.30-0.41 | NA |
| DJ Basin | 0.315 | 0.131 | 0.42 | 0.92 | 0.36-0.53 | 0.29-0.39 | 0.6 |
