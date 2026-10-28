# Supplementary File S25 — corrected M=20 inner-crossfit workflow

This folder hosts the source-traceable **Supplementary File S25** associated with:

**Exploratory Evaluation of the Incremental Prognostic Value of a Transcriptomic Model in Lung Adenocarcinoma**

## Download

**Archive:** [Supplementary_File_S25_M20_InnerCrossfit_2026-10-18.zip](Supplementary_File_S25_M20_InnerCrossfit_2026-10-18.zip)

**Size:** 18,494,034 bytes (about 18.5 MB)

**SHA-256:**
```text
de0bead86eff2c5c983904f4432350d50f7b9592c8ffa8163492f97640b78939
```

To check the downloaded archive on Windows PowerShell:

```powershell
Get-FileHash ".\Supplementary_File_S25_M20_InnerCrossfit_2026-10-18.zip" -Algorithm SHA256
```

Compare its `Hash` field with the SHA-256 above.

## What the archive contains

This ZIP is the **October 18, 2026 reissue** of the October 16 corrective inner-crossfit handoff. It was repackaged with provenance/integrity records; the analysis was **not rerun** for this reissue.

- 10 repeats × 5 outer folds, using 352 development patients, with 3,520 out-of-fold prediction rows.
- Corrected internal fold-stratified Harrell C-index results: clinical **0.655165**, combined **0.675445**, ΔC **+0.020280**.
- Conditional patient-cluster bootstrap 95% CI for ΔC: **−0.005836 to +0.044854** (1,000 replicates). The interval includes zero.
- Frozen joint model, coefficients, clinical comparator and selected prediction outputs.
- External outputs for GSE72094, GSE31210 eligible, and GSE31210 protocol-aligned analyses.
- Selected R scripts, PowerShell launch scripts, available logs, manifests and validation/audit reports.

## Where to look after extracting

All paths below are relative to the extracted ZIP root.

| File or folder | Purpose |
| --- | --- |
| `README_S25.md` | Reissue provenance, authoritative scope and limitations |
| `MANIFEST_SHA256.tsv` | Reissued archive integrity manifest |
| `ORIGINAL_SOURCE_SHA256.txt` | Hash of the original source ZIP |
| `M20_inner_crossfit_2026-10-16/README_FINAL.md` | Final original handoff summary |
| `M20_inner_crossfit_2026-10-16/output_audit/primary_fold_stratified_summary.tsv` | Primary C-index point estimates |
| `M20_inner_crossfit_2026-10-16/results/model/nested_cv/paired_oof_predictions.tsv` | Held-out predictions |
| `M20_inner_crossfit_2026-10-16/results/model/redesign/frozen_incremental_joint_model.rds` | Frozen fitted model |
| `M20_inner_crossfit_2026-10-16/final_audit/VALIDATION.tsv` | Run-level validation summary |
| `M20_inner_crossfit_2026-10-16/scripts/` | Selected R code accompanying the handoff |

## Scope and cautions

- **Not a complete end-to-end rerun:** the original raw expression matrices, large caches and installed R library are not included. Source checksums and indexes help trace required external inputs.
- **Not an independent replication:** this reissue preserved existing model outputs; integrity checks do not establish independent reproduction of the full modeling pipeline.
- The bootstrap confidence interval is **conditional on the existing out-of-fold predictions**; the complete selection pipeline was not refitted in every bootstrap replicate.
- Some earlier stage-specific notes inside the ZIP were superseded by `README_FINAL.md` and `final_audit/FINAL_AUDIT_REPORT.md`.
- Historical 30–730-day GSE72094 IBS results referenced in Supplementary Table S37 are not reproduced by this package.
- The recorded uncertainty must be retained; this archive does **not** establish clinically confirmed superiority.

For the study's broader documentation, see the [repository README](../../README.md). Source datasets are accessible from the [data access guide](../README.md).
