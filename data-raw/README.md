# Raw data

This directory contains the clinical files used for the METABRIC survival analysis.

## Data source

The data were obtained from the METABRIC study on cBioPortal (https://www.cbioportal.org/study/summary?id=brca_metabric), study identifier `brca_metabric`.

Access date: August 22, 2026.

## Required files 

The analysis requires the following files:

* `data_clinical_patient.txt`
* `data_clinical_sample.txt`
* `meta_clinical_patient.txt`
* `meta_clinical_sample.txt`

The files must retain these exact names because the analysis scripts refer to them using project-relative paths.

## Version-control policy

Raw clinical files are deliberately excluded from Git version control through `.gitignore`. They must be downloaded separately and placed in this directory before running the analysis.

The repository contains only analysis code, documentation, and non-identifying derived outputs.

## License

The source data are available under the Open Data Commons Open Database License. The original licensing statement is retained in `DATA_LICENSE.txt` in the repository root.
