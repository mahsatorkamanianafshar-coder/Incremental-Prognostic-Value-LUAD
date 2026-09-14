# Exact methods and reproducibility details for manuscript preparation

This document records the implemented settings in the LUAD project as of
2026-07-18. It is intended for Methods writing and reviewer responses. The code
and cited artifacts remain the source of truth.

## 1. Cohort flow and outcome construction

### TCGA expression cohort

The local TCGA-LUAD `SummarizedExperiment` contains 60,660 gene features and 599
samples: 540 Primary Tumor and 59 Solid Tissue Normal samples. TCGA barcodes were
reduced to the first 12 characters to define patient IDs.

### Overall-survival outcome

The clinical source table contains 585 unique participants. Overall survival was
constructed in `scripts/02_download_TCGA_LUAD_clinical.R` as follows:

- `OS_event=1` for case-insensitive `vital_status == "Dead"`.
- `OS_event=0` for case-insensitive `vital_status == "Alive"`.
- For deaths, `OS_time = days_to_death`.
- For living/censored participants, `OS_time = days_to_last_follow_up`.
- Other vital-status values, nonnumeric time, and missing time were set to missing.
- Duplicate participant rows are ranked deterministically: valid recognized
  outcome first, death over alive if statuses conflict, longest matching endpoint
  next, and smallest original source-row number last. The source currently has no
  duplicated participant IDs, so this rule does not change its 585 rows; the
  audit is saved in `TCGA_LUAD_clinical_row_selection_audit.tsv`.

The intermediate survival table contains 513 participants: 184 deaths and 329
censored observations. Four have zero follow-up and are retained in the
intermediate table for transparency but excluded before splitting/modeling by the
prespecified `OS_time > 0` criterion. There are 509 positive-time participants,
including 183 deaths. Five positive-time clinical participants have no eligible
Primary Tumor expression sample; the final expression-matched modeling cohort
therefore contains 504 patients and 182 deaths.

### Frozen split procedure

`scripts/00_create_leakage_safe_split.R` was run before feature discovery.
Eligibility required finite positive OS time, binary event status, a unique
patient ID, and at least one Primary Tumor sample. With random seed 123, patients
were grouped by `OS_event`; a uniform random ordering was generated within each
event stratum, and the first `floor(0.70*n)` in each stratum was assigned to
Training. The remainder was assigned to Testing.

| Partition | Censored | Deaths | Total |
|---|---:|---:|---:|
| Training | 225 | 127 | 352 |
| Testing | 97 | 55 | 152 |
| Total | 322 | 182 | 504 |

Patient IDs are unique and the partitions are disjoint. The Testing partition is
not used for DEG discovery, survival screening, Top-100 selection, lambda tuning,
coefficient estimation, imputation estimation, or cutoff estimation.

## 2. Normal-sample handling

> Final revision status (2026-07-18): full-training DEG, all 50 outer folds,
> final refit and every model-dependent downstream analysis were regenerated
> under the paired design.

### Full-training discovery

Training-only DEG analysis uses raw unstranded counts from Primary Tumor and Solid
Tissue Normal samples whose 12-character patient IDs belong to Training. For each
patient and sample type, the lexicographically smallest complete TCGA barcode is
retained. Only patients having both selected tumor and normal samples are used.
This yields 41 matched pairs (82 samples) and eliminates multiple-aliquot
pseudoreplication. DESeq2 uses `design = ~ patient_block + sample_type`, with a
syntactically safe factor copy of patient ID as the
blocking and Tumor-versus-Normal contrast. Count >=10 in at least three samples is
required. The regenerated analysis tested 27,722 genes, identified 7,407 genes at
adjusted P<0.05 and |log2FC|>=1, and retained 1,888 strict protein-coding DEGs at
adjusted P<0.01, |log2FC|>=1.5, and baseMean>=50.

### Normal samples inside each outer fold

For every outer fold, `scripts/lib_outer_fold_feature_discovery.R` receives only
the outer-analysis patient IDs. It then:

1. Selects only patients with both Tumor and Normal samples in that
   outer-analysis set.
2. Within patient and sample type, retains the lexicographically smallest full
   TCGA barcode.
3. Excludes all samples, including normal samples, belonging to outer-assessment
   patients.
4. Excludes normal-only patients who are not members of the survival-eligible
   frozen training cohort.
5. Repeats raw-count filtering and paired DESeq2 with
   `design = ~ patient_block + sample_type`.
6. Uses all outer-analysis tumor patients for survival screening, retaining the
   lexicographically smallest Primary Tumor barcode per patient.

The completed paired outer-fold run used 281-282 analysis patients and 27-38
matched pairs per fold. Between 26,361 and 27,477 genes were analyzed and
1,789-1,980 passed strict filtering. The shared fail-closed
`gene_type`/`gene_biotype` rule was used in full training and outer folds, with adjusted
DEG P<0.01, absolute log2 fold change >=1.5, and
baseMean >=50. The 500 leading strict DEGs were tested with univariate Cox models
after outer-analysis-only VST. Genes with nominal Cox P<0.05 were ranked by
`abs(log2FC) * -log10(Cox P)` and up to 100 were retained. Forty-eight folds had
100 eligible genes; two had 86 and 99. No nonsignificant genes were added merely
to fill the list.

All folds tested 500 survival candidates. Nominal survival-positive counts ranged
from 66 to 287. Forty-six folds selected 100 genes; the remaining four selected
96, 90, 74 and 66. Evidence: `tables/model/outer_fold_discovery_summary.tsv` and
`results/model/nested_cv/outer_fold_top100_discovery.tsv`.

## 3. Complete settings for extended validation (script 08b)

### Input populations and models

`scripts/08b_extended_model_validation.R` reads
`results/model/clinical_model_predictions.tsv` and separates the frozen Training
(n=352, 127 events) and Testing (n=152, 55 events) partitions. Four prespecified
scores are evaluated:

- Gene-only: `risk_score`.
- Clinical-only: `clinical_lp`.
- Clinical plus gene: `combined_lp`.
- Proliferation-only: `proliferation_score`.

For each score, a univariable Cox model is fitted in Training and is used to
generate test predictions at 365, 1,095, and 1,825 days. Baseline cumulative
hazard is obtained using `basehaz(centered=FALSE)`, and survival is calculated as
`exp[-H0(t)*exp(beta*score)]`.

### Calibration

- Test predictions are divided into five groups with `dplyr::ntile`.
- Mean predicted survival is compared with Kaplan-Meier survival in each group at
  each horizon; `extend=TRUE` is used when extracting KM survival.
- The reported cloglog calibration-in-the-large statistic is
  `cloglog(KM observed survival) - cloglog(mean predicted survival)`, with
  probabilities truncated to [1e-6, 1-1e-6].
- The overall calibration slope is obtained by fitting test survival against the
  linear predictor produced by the corresponding Training Cox score model. It is
  horizon-independent and repeated in the horizon table for joining.

Current slopes are 0.6599 (gene), 0.5200 (clinical), 0.4976 (combined), and 0.6217
(proliferation). These indicate miscalibration/overfitting for several models and
must not be presented as ideal calibration.

### IPCW Brier score and integrated Brier score

The censoring distribution is estimated in the Testing cohort with
`survfit(Surv(time, 1-event) ~ 1)`. IPCW weights are:

- `1/G(T_i-)` for an observed event by the horizon;
- `1/G(t)` for a participant known to survive beyond the horizon;
- zero for other observations.

Estimated censoring survival is floored at 1e-6. The outcome is survival beyond
the horizon. The integrated Brier score is a trapezoidal average over 102 unique
grid points from day 30 to day 1,825; 1, 3, and 5 years are forced into the
otherwise evenly spaced 100-point grid. It is not described as a time-zero IBS.

### Decision-curve analysis

- Threshold probabilities: 0.05 to 0.50 in increments of 0.01.
- Horizons: 365, 1,095, and 1,825 days.
- The full Testing denominator (n=152) is retained at every horizon.
- Cases by the horizon receive `1/G(T_i-)`, event-free controls receive `1/G(t)`,
  and patients censored before the horizon contribute zero weight.
- There are 9, 68, and 82 early-censored patients at 1, 3, and 5 years.
- Outcome is death by the horizon.
- Net benefit is `TP/n - FP/n * threshold/(1-threshold)`.
- Treat-all and treat-none comparators are saved.

This is an IPCW horizon-specific DCA, not a complete-case DCA.

### Bootstrap

- Population: untouched Testing set only.
- Replicates: 1,000 by default via `BOOTSTRAP_REPS`.
- Seed: 20260714.
- Sampling: n=152 patient rows with replacement per replicate.
- Statistic: Harrell concordance from `survival::concordance`, with
  `reverse=TRUE` for higher-risk scores.
- Metrics: all four model C-indices and combined-minus-clinical C-index.
- Intervals: empirical percentile 2.5th and 97.5th quantiles.
- Failed concordance calculations are returned as missing and summaries use
  `na.rm=TRUE`.

The bootstrap interval for combined-minus-clinical C-index includes zero.

### Proportional hazards and subgroup analysis

`cox.zph` is applied to each Training score model. Prespecified binary subgroups
are:

- Age: Older versus Younger using the Training median age (66 years).
- Stage: Advanced (III/IV) versus Early (I/II).
- Sex: Male versus Female.
- Smoking: Ever if quantitative smoking exposure >0, otherwise Never.

Training interaction models are `Surv ~ combined_lp * subgroup` with
`coxph.control(iter.max=50)`. A subgroup with only one observed level is marked
not estimable; this applies to smoking. Testing subgroup C-indices are reported
only when a subgroup has at least 20 participants and at least five events. No
multiplicity correction is applied, so the nominal sex interaction is exploratory.

### Output plots

Calibration is saved at 9x6 inches and DCA at 10x6 inches, both 300 dpi. The
calibration plot includes the 45-degree line and horizon facets. DCA uses separate
horizon facets with free y scales.

## 4. Clinical preprocessing

The actual clinical predictors are age, pathological stage, sex, and quantitative
smoking exposure.

- Age: `age_at_index` in years is preferred. Only if that entire source is absent
  is `age_at_diagnosis/365.25` used.
- Stage: strings are uppercased, an initial `STAGE` prefix is removed, and leading
  Roman numerals map to I-IV. Stage is modeled categorically with frozen levels
  I, II, III, and IV.
- Sex: the first available of `sex_at_birth`, `gender`, and `sex` is mapped to
  male=1 and female=0.
- Smoking: the first available quantitative field among `pack_years_smoked`,
  `cigarettes_per_day`, and `years_smoked` is used. In this dataset, pack-years
  is used as a continuous predictor. A binary tobacco-status fallback is used only
  if no quantitative smoking column has finite values.

Predictors must have at least two finite Training values/levels. Missing values
are imputed using Training only and the same values are applied to Testing:

| Variable | Imputation rule | Frozen value |
|---|---|---:|
| Age | Median | 66 |
| Stage | Mode | I (numeric code 1 before factor conversion) |
| Male indicator | Mode | 0 (female) |
| Smoking exposure | Median | 40 pack-years |

No missingness indicator is added. The clinical Cox model includes the four
predictors; the combined model adds the frozen gene risk score. Clinical and
combined linear predictors are learned in Training and applied unchanged to
Testing. Internal effect comparisons standardize each evaluated score to one SD
within the evaluation dataset; external analyses instead report HR per 0.1 frozen
risk-score unit.

The proliferation benchmark is the mean of fixed within-sample reference-panel
percentiles for available MKI67, PCNA, TOP2A, CCNB1, and CDC20; at least three are
required.

## 5. STRING settings

`scripts/07_PPI_network_STRING_hub_genes.R` uses:

- STRINGdb database version 11.5.
- Species 9606 (Homo sapiens).
- Minimum combined score 700 (high confidence).
- Input: 100 unique gene symbols from the full-Training Top-100 list.
- Mapping: `removeUnmappedRows=TRUE`; 100/100 mapped in the current run.
- Interactions are restricted to pairs in the input set.
- Network: undirected igraph; loops removed; duplicate edges simplified using the
  maximum combined score.
- Current network: 100 nodes and 1,154 edges; STRING combined scores retain the
  prespecified threshold of at least 700.
- Strength weights: `combined_score/1000`.
- Centralities: degree, weighted strength, normalized betweenness, and normalized
  closeness.
- Hub score: sum of z-scaled degree, strength, betweenness, and closeness; a
  constant metric is assigned a zero z-score.
- Explainable network priority: z-scaled expression/survival priority plus
  z-scaled hub score.
- Plot seed: 123; Fruchterman-Reingold layout; only the 20 highest-degree nodes are
  labeled; node size is `4 + 2*log1p(degree)` and edge width scales with STRING
  combined score.

STRING is used for biological prioritization/explainability, not as an
outer-fold outcome-selection step in the fully nested performance estimate.

## 6. DGIdb settings and provenance limitation

### Direction-aware signature

`scripts/11_prepare_drug_repurposing_signature.R` defines:

- Upregulated risk gene -> decrease/inhibit.
- Downregulated protective gene -> increase/activate.
- Discordant expression/survival directions -> caution.
- Final-model genes are labeled separately.
- High network priority is degree at or above the 75th percentile.
- The 20 leading risk hubs by hub score form the broad inhibition query.

### Versioned DGIdb query

`scripts/12a_query_DGIdb_graphql.R` submits the exact regenerated 20-gene list to
the official `https://dgidb.org/api/graphql` endpoint. The frozen query used
DGIdb 5.0.11 on 2026-07-18 UTC. Requested and returned gene lists, endpoint,
version, timestamp, export method and coverage are written to
`DGIdb_query_metadata.tsv` and `DGIdb_query_gene_coverage.tsv`.

All 20 intended genes were sent. DGIdb recognized 14 (70%) and reported at least
one interaction for seven. The export contains 516 distinct interaction rows.
The six unrecognized genes and seven recognized genes without interactions are
explicit database-coverage results, rather than omissions from the request.
After direction-aware joining, all 516 interaction rows remain relevant and
yield 419 candidate drug names.

### Drug ranking

Only `drug_desired_action == "Decrease_or_inhibit"` is retained. For each drug:

- Target genes are collapsed to a semicolon-separated list.
- The number of distinct signature targets is counted.
- Mean gene HR, mean network priority, and maximum hub score are calculated.
- Infinite values are converted to missing.
- Priority score = number of targets + mean network priority (or zero if missing)
  + maximum hub score (or zero if missing).

The top 20 are saved. `scripts/13_filter_candidate_drugs.R` then uses a manually
specified annotation table for mechanism and manuscript category. Unlisted drugs
are marked for manual review. This manual categorization is heuristic, not an
independent validation layer. Regulatory-approval, indication, and DGIdb
interaction-score columns are retained in clean interaction data but are not part
of the implemented drug-priority formula.

## 7. External microarray mapping and frozen scoring

### Shared procedure

Both scripts preferentially combine the archived local series-matrix and platform
annotation files; `getGEO(GSEMatrix=TRUE, AnnotGPL=TRUE)` is a network fallback.
If multiple matrices are returned, the largest sample count is selected.
GEO `characteristics` fields formatted as `key: value` are parsed into cleaned
columns. The supplied normalized expression matrices are used; no new MAS5, RMA,
IRON, or batch correction is fitted in this project.

Candidate annotation columns match `gene.?symbol|symbol|hgnc|gene.?name`. The
column maximizing overlap with the 9 frozen model genes is selected; expression
row names are a fallback. Multi-symbol annotations are split on `///`, `//`,
semicolon, comma, or pipe. Frozen model/reference symbols are preferred over the
first arbitrary token. The frozen alias `PCLAF -> KIAA0101` is applied, and this
also permits selection of CLEC3B when it is the second symbol in a multi-symbol
GPL570 annotation. Multiple probes mapping to one selected symbol are collapsed
by arithmetic mean.

Sensitivity outputs compare arithmetic mean with median and the single
highest-IQR probe, and repeat scoring after removing each reference gene in turn.
In GSE31210, correlations with primary scores were 0.847 for median and 0.858 for
highest-IQR mapping. In GSE72094 they were 0.993 and 0.908, respectively. All
sensitivity Cox effects retained the same harmful-score direction.

All 9 model genes and all 15 reference genes must be present; otherwise scoring
stops. For model gene g and sample i, the frozen feature is:

`[1 + number of the 15 reference genes with expression <= expression(g,i)] / 16`.

The risk score is the dot product of these 9 within-sample percentile features
and the frozen coefficients. No external cohort mean, SD, coefficient, reference
panel, or cutoff is estimated. The frozen TCGA-training median cutoff 1.546816 is
applied only for exploratory binary grouping. Primary external inference uses the
continuous score divided by 0.1, so HRs are per 0.1 score unit.

### GSE31210

- Platform: GPL570.
- Source matrix: 54,675 probes/features and 246 samples.
- GEO-reported processing: MAS5 normalization.
- Selected outcome fields: `days_before_death_censor_ch1` and `death_ch1`.
- Numeric strings are cleaned; string event labels map dead/deceased/death/died to
  1 and alive/living/censor to 0. Numeric 0/1 is preserved; numeric 1/2 maps to
  alive/dead.
- Eligibility: finite positive time and binary nonmissing event; 226 samples and
  35 deaths remain.
- Cox model: `Surv ~ risk_score/0.1` (coxph default Efron ties).
- Time-dependent ROC: 365, 1,095, and 1,825 days with `timeROC`, cause=1,
  marginal weighting, and `iid=TRUE`.
- Frozen cutoff grouping is 226 low-risk and zero high-risk; KM inference is not
  estimable and continuous-score analysis is primary.

### GSE72094

- Platform: GPL15048.
- Source matrix: 60,607 probes/features and 442 samples.
- GEO-reported processing: IRON normalization followed by PLS batch correction
  for RNA quality.
- Selected outcome fields: `survival_time_in_days_ch1` and `vital_status_ch1`.
- At least 50 positive times, 50 valid events, five deaths, and five censored
  observations are required before proceeding.
- Eligibility: finite positive time and binary nonmissing event; 398 samples and
  113 deaths remain after matching expression and phenotype.
- Cox model: `Surv ~ risk_score/0.1`, ties=`efron`, with design matrix retained.
- A Cox PH diagnostic is saved for the external continuous score.
- Frozen cutoff grouping is 388 low-risk versus 10 high-risk; binary/KM results
  are exploratory because of the imbalance.

## 8. Manuscript wording constraints

- Call the nested estimate “fully outer-fold feature-discovery nested CV.”
- Describe Top-100 as “up to 100 nominally survival-associated genes.”
- State that DEG is matched-pair, patient-blocked, and uses one deterministic
  barcode per patient/sample type.
- Call 08b DCA IPCW horizon-specific DCA with the full test denominator.
- State that IBS integrates 102 points from day 30 to day 1,825.
- Describe DGIdb as an exact versioned GraphQL query, while reporting 14/20 API
  recognition and 7/20 genes with interactions as database-coverage limitations.
- Do not claim external platform renormalization; scoring uses supplied normalized
  expression plus the frozen within-sample percentile transformation.
- Keep continuous external Cox analysis primary because frozen-cutoff groups are
  severely imbalanced.
