# Exploratory Evaluation of the Incremental Prognostic Value of a Transcriptomic Model in Lung Adenocarcinoma

This repository provides **selected research outputs and reproducibility documentation** associated with the study:

**“Exploratory Evaluation of the Incremental Prognostic Value of a Transcriptomic Model in Lung Adenocarcinoma.”**

## Study overview

The study examines whether transcriptomic information adds prognostic value beyond a clinical comparator in lung adenocarcinoma (LUAD). The modeling workflow includes repeated nested cross-validation, Elastic Net Cox regression, leakage-controlled feature selection, and assessment in external microarray cohorts. The reported improvement is modest, with uncertainty that includes no improvement. The findings are exploratory and should not be interpreted as a clinically validated prediction tool.

## Source cohorts

- [TCGA-LUAD (NCI Genomic Data Commons)](https://portal.gdc.cancer.gov/projects/TCGA-LUAD)
- [GSE72094 (NCBI GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72094)
- [GSE31210 (NCBI GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210)

The source expression datasets are **not redistributed** here. Obtain them from their original repositories and follow the documentation for study-specific preprocessing and endpoint definitions.

## Repository contents

| Location | Contents |
| --- | --- |
| [data/](data/) | Source-dataset guidance and supplementary archive |
| [data/supplementary/](data/supplementary/) | **Supplementary File S25**, with audit outputs, selected R scripts, model objects and predictions |
| [figures/](figures/) | Selected model and sensitivity-analysis figures |
| [DATA_DICTIONARY.md](DATA_DICTIONARY.md) | Variable and output definitions |
| [INCREMENTAL_REDESIGN_DESIGN_LOCK.md](INCREMENTAL_REDESIGN_DESIGN_LOCK.md) | Versioned methodological design notes, including later amendments |
| [METHODS_REPRODUCIBILITY_DETAILS.md](METHODS_REPRODUCIBILITY_DETAILS.md) | Historical methodological documentation (dated July 18, 2026; see note below) |
| [EXPERIMENTAL_VALIDATION_PROTOCOL.md](EXPERIMENTAL_VALIDATION_PROTOCOL.md) | Proposed experimental follow-up; **not** completed wet-lab validation |
| [MANIFEST.csv](MANIFEST.csv) | Repository project manifest |
| [renv.lock](renv.lock) | Recorded R dependency lockfile |
| [CITATION.cff](CITATION.cff) | Citation metadata |

### Supplementary File S25

**File:** [Supplementary_File_S25_M20_InnerCrossfit_2026-10-18.zip](data/supplementary/Supplementary_File_S25_M20_InnerCrossfit_2026-10-18.zip)

**SHA-256:**
`de0bead86eff2c5c983904f4432350d50f7b9592c8ffa8163492f97640b78939`

S25 documents the corrected M=20 inner-crossfit analysis and associated outputs. For navigation, contents, integrity checks and limitations, see the [supplementary archive README](data/supplementary/README.md) and `README_S25.md` inside the archive.

## Reproducibility and scope

This repository is **not a complete, stand-alone raw-data-to-results reproduction environment**:

- There is **no top-level `scripts/` directory** in this repository. Selected R scripts and PowerShell launch scripts are included **inside the S25 ZIP** under `M20_inner_crossfit_2026-10-16/`.
- The archive does not contain the full raw count matrices, large expression/discovery caches, or the original package library. An exact rerun requires the corresponding inputs and environment.
- The S25 ZIP was repackaged and checked for integrity, not independently rerun from raw data.
- Historical files in this repository describe different analysis stages. In particular, `METHODS_REPRODUCIBILITY_DETAILS.md` is a July 18, 2026 snapshot and should **not** be treated as a fully updated description of the later corrected M=20 primary analysis. The associated final manuscript, supplementary manuscript and the later S25 handoff govern final result reporting.

## Citation

Repository citation metadata are available in [CITATION.cff](CITATION.cff). The journal article citation and DOI can be added after publication.

## Author

Mahsa Torkamanian-Afshar  
Department of Software Engineering  
Faculty of Engineering and Natural Sciences  
İstanbul Topkapı University  
Istanbul, Türkiye  

ORCID: [0000-0002-8658-4013](https://orcid.org/0000-0002-8658-4013)

## License

See [LICENSE](LICENSE) for the repository's stated license. Third-party source datasets remain subject to their own access and reuse terms.
