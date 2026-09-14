# Incremental prognostic-model redesign: locked analysis plan

## Primary-analysis migration amendment (2026-08-12)

Following reviewer-methodology review, the corrected M=20
smoking-constrained PMM pipeline replaces the historical M=5
legacy-unconstrained pipeline as the canonical primary analysis. This change
is method-driven: known never-smokers with missing pack-years are assigned the
prospectively defined structural zero, and known ever-smokers draw only from
known-ever donor pools. M=50 remains a computational sensitivity because it
produced essentially the same estimates and identical final genes as M=20.

The primary discrimination estimand is amended to fold-stratified Harrell
concordance. Comparable patient pairs are formed only within each outer
assessment fold; fold-specific concordances are weighted by their numbers of
comparable pairs within each repeat, and the ten repeat estimates are averaged.
This avoids comparisons between linear predictors from independently fitted
outer-fold models and requires no between-fold LP rescaling. Patient-cluster
bootstrap intervals recalculate this same estimand but remain conditional on
the fitted OOF predictions rather than rerunning the full modeling pipeline.

The canonical frozen model is version
`incremental_redesign_v3_M20_smoking_constrained_primary`. It contains ten
genes, with CCNA2 replacing CBX7 relative to the archived historical M=5
model. The historical model and overwritten canonical outputs are retained in
`archive/M5_frozen_primary_before_M20_migration_2026-08-12`.

Locked before executing the redesigned model on 2026-07-21.

Methodological audit amendment implemented on 2026-07-23:
the candidate universe, all inner-fold preprocessing, and absolute-risk
calibration were made symmetric between resampling and final refit. This
amendment responds to audit findings and is not outcome-driven.

## Cohort roles

- Development cohort: the 352 patients labelled `Training` in the frozen split.
- Secondary internal cohort: the 152 patients labelled `Testing`; this cohort has
  been observed previously and is not described as untouched, independent, or
  confirmatory.
- External cohorts: GSE72094 (primary incremental assessment) and GSE31210
  (limited exploratory incremental assessment because it has 35 events).

## Primary estimand and model-selection criterion

The primary criterion is the paired repeat-level difference in Harrell C-index
between the joint combined model and the clinical-only model in 10 repeats of
5-fold outer cross-validation within the 352-patient development cohort.
Secondary criteria are paired differences in IBS, 1/3/5-year Brier score and
time-dependent AUC, calibration, model size, and gene-selection stability.

No minimum improvement is imposed. A negative or uncertain result is retained.

## Locked modeling grid

- Outer CV: 10 repeats x 5 folds, event-stratified.
- Inner CV: 5 folds, event-stratified.
- Candidate input: the first 20 genes after the locked screening/ranking
  sequence; candidate-set size is no longer tuned.
- Elastic-net alpha: 0, 0.25, 0.50, 0.75, 1.
- Lambda candidates: `lambda.min` and `lambda.1se`.
- Complexity constraint: at most 15 nonzero gene coefficients.
- Near-tie rule: within 0.005 C-index of the best valid inner result, select
  by lower IBS, fewer genes, `lambda.1se`, then larger alpha.
- Current fifteen-gene reference panel and within-sample percentile
  transformation are retained for the first redesign.

## Clinical preprocessing

- Stage is categorical (I-IV).
- Sex is binary using the recorded birth-sex field.
- Age and stage missingness are imputed using analysis-fold medians/modes.
- The primary comparator uses age, stage, sex, and log1p(pack-years).
- Pack-years is multiply imputed by predictive mean matching using donors from
  the corresponding analysis fold only; outer-assessment values and outcomes
  are not used by the imputer.
- Sensitivities omit smoking or use raw ever/never smoking status.
- Linear and restricted-cubic-spline representations of age and pack-years are
  compared inside development resampling only.

## Feature screening and models

- Every outer-analysis and inner-training partition independently reruns
  deterministic paired tumor-normal selection, paired DESeq2 with
  `design = ~ patient_block + sample_type`, strict protein-coding filtering,
  and strict Top-500 selection.
- Within each inner-training partition, the Top-500 is followed by univariate
  survival screening (`P < 0.05`), clinical-adjusted ranking, and a maximum
  Top-20 Elastic Net input. No archived Top-100 universe is used.
- Imputation, survival screening, adjusted ranking, glmnet standardization,
  and alpha/lambda fitting are learned only from the relevant inner-training
  partition. The inner-assessment partition is used only for prediction and
  performance calculation.
- Univariate and clinical-adjusted survival screens are stored separately.
- Clinical-adjusted screening uses
  `Surv(time,event) ~ gene_z + age + stage + sex + smoking` and Rubin pooling
  across imputed analysis-fold datasets.
- Gene-only models use univariate-ranked candidates.
- Combined models use clinical-adjusted-ranked candidates.
- Combined Elastic Net uses penalty factor zero for all clinical design columns
  and one for all gene columns.

## Absolute-risk calibration

- Clinical-only and combined models use the same development-only calibration
  form: `Surv(time, event) ~ linear_predictor`.
- Each frozen model stores the calibration slope and baseline cumulative
  hazard from that calibration fit. External cohorts never refit either
  quantity.

## Reporting restriction

The 152-patient cohort is reported as “previously observed secondary internal
validation cohort.” It is not used to select preprocessing, covariates,
features, alpha, lambda, candidate limit, or cutoff.

## Reviewer M3 shared-calibration diagnostic (post hoc amendment; 2026-08-12)

- This diagnostic uses only the canonical M20 smoking-constrained model.
- Clinical-only and combined development linear predictors are separately
  standardized with their full-development means and standard deviations.
- The two representations are stacked with weight 0.5 per row, and one Cox
  calibration model estimates a single slope and a single baseline cumulative
  hazard from the development cohort only.
- The frozen development transformations, shared slope, and shared baseline
  hazard are applied unchanged to the secondary and external cohorts. No
  cohort-specific recalibration is permitted.
- The diagnostic is supplementary and does not replace the locked primary
  model-specific absolute-risk calibration. Its purpose is to test whether
  conclusions about relative IBS depend on using separate calibration maps.
