# Data dictionary

| Field | Meaning |
|---|---|
| `patient_id` | TCGA or GEO patient/sample identifier |
| `OS_time` | Overall-survival follow-up time in days unless explicitly stated |
| `OS_event` | 1 = death, 0 = censored |
| `split` / `dataset` | Frozen Training or Testing membership |
| `risk_score` | Frozen gene-model linear predictor |
| `risk_group` | Classification using the TCGA-training cutoff |
| `clinical_lp` | Linear predictor from the clinical-only Cox model |
| `combined_lp` | Linear predictor from the clinical-plus-gene Cox model |
| `proliferation_score` | Prespecified MKI67/PCNA/TOP2A/CCNB1/CDC20 benchmark |
| `predicted_survival` | Model-estimated survival probability at a stated horizon |
| `brier_score` | IPCW prediction error; lower is better |
| `net_benefit` | Decision-curve net benefit at a stated threshold |
| `selection_frequency` | Fraction of nested runs selecting a gene by LASSO |
| `outer_screen_frequency` | Fraction of outer folds passing a gene through Cox screening |
| `connectivity_score` | LINCS/CMap signature connectivity, expected range -100 to 100 |

Raw and processed tables are tab-delimited. Missing values are represented as `NA`.
