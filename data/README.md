# Data access

Raw TCGA-LUAD and GEO expression matrices, and large intermediate expression/discovery caches, are **not included** in this repository.

## Public source datasets

- [TCGA-LUAD — NCI Genomic Data Commons](https://portal.gdc.cancer.gov/projects/TCGA-LUAD)
- [GSE72094 — NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72094)
- [GSE31210 — NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210)

Use the study's associated manuscript and methodological documents for cohort selection, outcome definitions, and preprocessing details.

## Supplementary audit archive

[Supplementary File S25 and its documentation](supplementary/README.md) provide selected scripts, predictions, model objects, audit results, and integrity manifests from the corrected M=20 workflow.

**Important:** The repository has no top-level `scripts/` folder. Selected scripts are packaged *inside* the S25 ZIP. The ZIP is not a complete raw-data-to-results workflow; original source data, caches and a compatible R environment are needed for a full rerun.
