# r-analysis-portfolio
This is the repository to my new project conducted in R using publicly available health data, more details will come soon.

### *PROJECT UPDATE* 23.08: The primary overall survival analysis is complete. The initial full follow-up Cox model showed substantial violations of the proportional hazards assumption, so the final analysis used interval-specific Cox models for 0-2, 2-5, 5-10, and beyond 10 years. Each model adjusted for age at diagnosis and Nottingham Prognostic Index using natural cubic splines and was stratified by source cohort. Holm adjustment was applied separately to the four global subtype tests and to the 20 comparisons between individual subtypes and Luminal A across the four follow-up intervals. Age- and NPI-stratified sensitivity analyses produced similar estimates, with all changes in subtype hazard ratios below 10%. The final tables, proportional hazards diagnostics, sensitivity analyses, and interval-specific hazard ratio figures are generated reproducibly by scripts `06`–`12`.

### P.S: The originally planned C-index comparison was removed because the final analysis no longer relies on a single full follow-up model with constant hazard ratios, making one overall discrimination estimate difficult to interpret.

# Investigation into the prognostic contribution of molecular subtype in breast cancer

## A survival analysis of the METABRIC cohort

This project investigates whether the PAM50 + Claudin-low molecular subtype classification provides prognostic information for survival beyond age at diagnosis and the Nottingham Prognostic Index (NPI) in women with primary breast cancer.

The project is a reproducible observational survival analysis conducted in R. It examines prognostic associations and model discrimination. It does not estimate causal treatment effects or produce a model intended for clinical decision-making.

## 1. Clinical question and objectives

### Primary research question

Among women with primary breast carcinoma in METABRIC, does the PAM50 + Claudin-low subtype classification provide additional prognostic information for overall survival beyond age at diagnosis and the Nottingham Prognostic Index (NPI)?

### Primary objective

To compare:

1. A clinical Cox model containing age at diagnosis and Nottingham Prognostic Index (NPI)
2. An extended Cox model containing age, Nottingham Prognostic Index (NPI), and molecular subtype

The analysis will assess whether adding molecular subtype:

- Improves the statistical fit of the survival model
- Improves the model’s ability to rank patients according to mortality risk
- Remains associated with overall survival after adjustment for established clinical prognosis

### Secondary objective

To examine associations between molecular subtype and relapse-free survival endpoint using the same clinical adjustment strategy.

### Interpretation

The study evaluates ***prognostic associations***. It does not attempt to show that molecular subtype or any recorded treatment causes a survival outcome.

## 2. Study design and data source

This is a retrospective cohort study using patient-level clinical data from the Molecular Taxonomy of Breast Cancer International Consortium.

The cBioPortal release contains 2,509 patients with corresponding primary breast tumor samples. The original METABRIC study included large discovery and validation collections of primary breast tumors with long-term clinical follow-up.

The analysis uses:

* `data_clinical_patient.txt`
* `data_clinical_sample.txt`
* `meta_clinical_patient.txt`
* `meta_clinical_sample.txt`

The patient and sample files are joined using `PATIENT_ID`.

Data inspection confirmed the following:

* 2,509 unique patients
* 2,509 unique samples
* One sample per patient
* No duplicated patient or sample identifiers
* All samples identified as primary tumors
* All participants recorded as female

Data references:

* [METABRIC on cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric)
* [Curtis et al., Nature 2012](https://pubmed.ncbi.nlm.nih.gov/22522925/) - original METABRIC publication describing the main discovery and validation cohorts, molecular profiling, and subtype framework.
* [Pereira et al., Nature Communications 2016](https://www.nature.com/articles/ncomms11479) - describes the expanded METABRIC profiling reflected in the current cBioPortal release.

## 3. Study population and eligibility criteria

### Source population

The source population contains 2,509 women with primary breast cancer. Age at diagnosis ranges from approximately 22 to 96 years. No additional age restriction is required.

### Overall-survival cohort

Patients are included if they have:

* A primary breast carcinoma
* An interpretable PAM50 + Claudin-low subtype
* A recorded overall-survival time
* A known overall-survival status
* Usable survival information, defined as observable follow-up after diagnosis or a death recorded at the start of follow-up

Patients are excluded for the following reasons:

* 529 have no molecular-subtype classification
* 6 are coded `NC` and cannot be assigned to an interpretable subtype
* 2 have breast angiosarcoma, a non-epithelial malignancy outside the target population
* 1 is alive and censored at month zero, contributing neither follow-up time nor an event

Among the 529 patients without subtype information, 528 also lack overall-survival information. These patients are counted once under missing subtype in the participant flow.

The final overall-survival cohort contains:

* 1,971 patients
* 1,138 deaths
* 833 censored observations

### Relapse-free survival cohort

The relapse-free survival (RFS) cohort will include patients who meet all of the following criteria:

* Eligible epithelial primary breast tumor.
* Available PAM50/Claudin-low subtype classification.
* Available RFS status and RFS follow-up time.
* No recorded evidence of stage IV or distant metastatic disease at baseline.
* RFS information that is consistent with the defined study baseline and subsequent follow-up.

Patients will be excluded from the RFS cohort if they have:

* A non-epithelial breast malignancy.
* Missing PAM50/Claudin-low subtype classification.
* Missing RFS status or RFS follow-up time.
* Known stage IV or distant metastatic disease at baseline.
* Documented relapse before or at the study baseline.
* An invalid or irreconcilable sequence between baseline, follow-up, and relapse.

Patients will ***not*** be excluded solely because relapse or death occurred shortly after baseline. Excluding valid early events would preferentially retain patients who remained event-free for longer and could bias the survival estimates. Events recorded at zero months will be reviewed individually. They will be retained if the event occurred after baseline and the zero value reflects rounding or limited time precision. They will be excluded if the record indicates that relapse or metastatic disease was already present at baseline or if the event chronology cannot be reconciled.

Missing tumor stage alone will ***not*** be an exclusion criterion. Patients with missing stage will be retained unless another variable indicates metastatic disease at baseline. This avoids excluding a substantial proportion of the cohort based only on incomplete stage documentation. However, because metastatic status cannot be confirmed for these patients, this uncertainty will be reported as a study limitation.

Patients will also ***not*** be excluded solely because they were censored, had short but observable follow-up, received or did not receive a particular treatment, had missing nonessential covariates, or belonged to a small but valid molecular subtype group.

## 4. Outcomes

### Primary outcome: overall survival

Overall survival is measured in months from initial breast cancer diagnosis until death from any cause.

* Time variable: `OS_MONTHS`
* Event: `OS_STATUS = 1:DECEASED`
* Censored: `OS_STATUS = 0:LIVING`

Patients alive at their last recorded follow-up are censored at that time.

A patient censored at month zero contributes neither observable follow-up nor an event and is therefore excluded. A death recorded at month zero would be retained because it represents an observed event.

### Secondary outcome: relapse-free survival

The METABRIC data dictionary defines a relapse-free-survival event as:

* Locoregional relapse
* Distant relapse
* Disease-specific death

Patients without one of these events are censored at their recorded follow-up time.

* Time variable: `RFS_MONTHS`
* Event: `RFS_STATUS = 1:Recurred`
* Censored: `RFS_STATUS = 0:Not Recurred`

Three patients have an RFS event recorded at month zero. Review of the available clinical fields found no recorded stage IV disease or other evidence establishing that relapse preceded baseline; all three also have recorded breast surgery and subsequent overall-survival follow-up. They will therefore be retained as observed early events, while acknowledging that the data do not provide exact event dates within the month. Patients censored at month zero will be excluded because they contribute neither observable follow-up nor an event.

## 5. Primary exposure

The primary exposure is `CLAUDIN_SUBTYPE`, described in the data dictionary as the PAM50 + Claudin-low subtype classification.

The six analysis categories are:

* Luminal A
* Luminal B
* HER2-enriched
* Basal-like
* Normal-like
* Claudin-low

Luminal A is the prespecified reference group because it is the largest category and provides a clinically interpretable comparison.

The variable will not be described as ordinary PAM50 because Claudin-low is included as an additional category.

The global contribution of molecular subtype will be evaluated before individual subtype comparisons are interpreted.

## 6. Clinical predictors

### Primary clinical predictors

The primary clinical model includes:

* Age at diagnosis
* Nottingham Prognostic Index (NPI)
* Cohort as a stratification variable

Age and NPI are complete for all 1,971 patients in the overall-survival cohort.

NPI combines information from:

* Tumor size
* Histologic grade
* Lymph node involvement

NPI will therefore not be entered into the same model as all three of its components.

### Continuous variable modeling

Age and NPI will remain continuous. They will not be divided into arbitrary categories.

Restricted cubic splines with four knots at standard quantile locations will be used to allow nonlinear associations between:

* Age and the log hazard
* NPI and the log hazard

The fitted relationships will be presented graphically to support interpretation.

### Cohort stratification

The eligible patients belong to component cohorts numbered 1-5. These cohorts differ in size and follow-up patterns.

The Cox models will be stratified by `COHORT`. This allows each cohort to have its own baseline hazard while estimating common associations for age, NPI, and molecular subtype.

The cohort variable is included to account for differences between component cohorts and will not be interpreted as a biological predictor.

### Variables not included in the primary model

Tumor stage is not included because it is missing for approximately 26% of the eligible cohort.

Tumor size, grade, and lymph node involvement are not entered separately because their prognostic information is represented by NPI.

ER, PR, and HER2 status are not added to the primary model because they overlap biologically with the molecular subtype classification.

Chemotherapy, hormone therapy, and radiation therapy are described as cohort characteristics but are not used to estimate treatment effects.

## 7. Descriptive analysis

The analysis will report:

* Participant characteristics overall
* Participant characteristics by molecular subtype
* Age and NPI distributions
* Tumor size, grade, and lymph node involvement
* Receptor status
* Recorded treatment indicators
* Numbers of deaths, relapse-free survival events, and censored observations
* Missingness for every reported variable
* Follow-up duration using the reverse Kaplan-Meier method
* Numbers of patients and events within each molecular subtype

Categorical variables will be summarized using counts and percentages.

Continuous variables will be summarized using means and standard deviations or medians and interquartile ranges, depending on their distributions.

Routine hypothesis testing will not be used in the descriptive table because molecular subtype groups are not randomized groups.

## 8. Unadjusted survival analysis

### Kaplan-Meier estimation

Kaplan-Meier curves will describe overall survival across the six molecular subtype categories.

The figure will include:

* Survival probabilities
* 95 % confidence intervals
* Censoring marks
* Numbers at risk
* Clearly labeled time units
* A color-blind-accessible palette

Median survival will be reported only for subtype groups in which the estimated survival curve reaches 50%.

Survival probabilities at five and ten years will be reported with confidence intervals.

### Log-rank test

A global log-rank test will compare the unadjusted survival distributions across all six molecular subtypes.

The global comparison will be interpreted before examining individual subtype differences.

The log-rank test will be treated as an unadjusted comparison and will not replace the multivariable Cox analysis.

## 9. Multivariable Cox analysis

### Clinical model

The clinical Cox model will contain age, NPI, and cohort stratification:

```text
Overall survival ~ spline(age) + spline(NPI) + strata(cohort)
```

### Extended molecular model

The extended Cox model will add molecular subtype:

```text
Overall survival ~ spline(age) + spline(NPI) +
                   molecular subtype + strata(cohort)
```

### Model estimates

The analysis will report:

* Adjusted hazard ratios
* 95 % confidence intervals
* A global test for molecular subtype
* Individual subtype contrasts using Luminal A as the reference
* Corresponding p-values
* Graphical presentation of the adjusted subtype estimates

The number of model parameters is small relative to the 1,138 observed deaths, providing sufficient information for the prespecified model.

No automated stepwise variable selection will be used.

### Relapse-free survival model

The same clinical and molecular model structure will be applied to the relapse-free-survival outcome:

```text
Relapse-free survival ~ spline(age) + spline(NPI) +
                        molecular subtype + strata(cohort)
```

The relapse-free survival analysis is secondary and will be interpreted in light of the dataset’s composite event definition.

## 10. Added prognostic information and discrimination

### Likelihood-ratio test

The clinical and extended overall-survival models will be compared using a likelihood-ratio test.

This test evaluates whether adding molecular subtype improves model fit beyond age, NPI, and cohort stratification.

The global likelihood-ratio test is the primary assessment of the added statistical contribution of molecular subtype.

### Harrell’s concordance index

Harrell’s concordance index will be calculated for:

* The clinical model
* The extended molecular model

The C-index measures how well each model ranks patients according to mortality risk while accounting for censored observations.

The analysis will report:

* The C-index for each model
* Ninety-five percent confidence intervals
* The absolute difference in C-index after adding molecular subtype

The C-index comparison will be interpreted as a measure of model discrimination, not as proof of clinical usefulness.

Because the C-index is calculated using the same cohort in which the model is fitted, it will be described as apparent model discrimination.

## 11. Model diagnostics

### Proportional-hazards assumption

The proportional-hazards assumption will be assessed using:

* Scaled Schoenfeld residuals
* Variable-specific tests
* A global proportional-hazards test
* Diagnostic plots against follow-up time

The statistical tests and plots will be interpreted together. A small p-value alone will not determine whether a violation is practically important.

A constant hazard ratio will be interpreted only when the corresponding proportional-hazards assumption is adequately supported.

Any material violation will be reported explicitly, and the affected coefficient will not be described as having a constant effect throughout follow-up.

### Continuous-variable functional form

The fitted spline relationships for age and NPI will be inspected for:

* Nonlinearity
* Implausible behavior at the boundaries
* Sparse data in extreme ranges
* Excessively complex patterns unsupported by the data

### Subtype sample size

The numbers of patients and events within each molecular subtype will be examined before interpreting subtype-specific estimates.

Hazard ratios with wide confidence intervals will be described as imprecise rather than treated as evidence of no association.

## 12. Prespecified sensitivity analyses

### Stage restriction

The overall-survival analysis will be repeated after excluding patients recorded as:

* Stage 0
* Stage IV

This analysis assesses whether the primary findings are affected by patients whose recorded stage is inconsistent with the main early-to-locally-advanced invasive breast cancer population.

Patients with missing stage will not be excluded solely because stage information is unavailable.

### Normal-like subtype

The analysis will be repeated after excluding the Normal-like category.

This assesses whether conclusions are influenced by a subtype that may partly reflect low tumor cellularity or normal-tissue contamination.

### Month-zero relapse events

The relapse-free-survival analysis will be repeated after excluding the three patients with events recorded at month zero.

This assesses whether immediate recorded events materially affect the secondary results.

## 13. Participant flow

The overall-survival cohort will be presented as:

```text
2,509 patients in the source dataset
  - 529 missing molecular subtype
  - 6 coded as unclassified subtype
  - 2 breast angiosarcomas
  - 1 alive and censored at month zero
= 1,971 patients in the overall-survival cohort
```

The relapse-free-survival cohort will be presented as:

```text
The relapse-free-survival cohort will be presented as:

```text
1,972 patients meeting tumor and subtype criteria
  - 1 missing relapse-free survival event status
  - 1 censored at month zero
  - 10 with known stage IV disease at baseline
= 1,960 patients in the relapse-free survival cohort
```

The final relapse-free-survival cohort contains 790 events and 1,170 censored observations. The three patients with events recorded at month zero remain included.

Every exclusion will be implemented in code and recorded in a participant-flow table.

## 14. Planned results presentation

No numerical survival results will be reported before the analysis is completed.

### Main Figure 1

Kaplan-Meier overall-survival curves by molecular subtype with confidence intervals and numbers at risk.

### Main Figure 2

Forest plot of adjusted molecular-subtype hazard ratios from the extended Cox model.

This figure will be used only for coefficients whose proportional-hazards assumptions support interpretation as constant hazard ratios.

### Main Figure 3

Comparison of the clinical and extended model C-indices, including confidence intervals and the absolute change after adding molecular subtype.

### Secondary figure

Kaplan-Meier relapse-free-survival curves by molecular subtype.

### Main tables

* Table 1: Participant characteristics overall and by molecular subtype
* Table 2: Clinical and extended overall-survival model estimates
* Table 3: Model comparison, C-index results, and sensitivity analyses

## 15. Limitations

* METABRIC is a retrospective research cohort rather than a population-representative cancer registry.
* The source tumors were collected between 1977 and 2005, and survival patterns may not represent patients receiving contemporary treatment.
* Much of the HER2-positive cohort predates routine trastuzumab treatment, which is important when interpreting HER2-enriched survival.
* Molecular subtype is unavailable for 529 patients, and 528 of these patients also lack overall-survival information. The analyzed molecular cohort may therefore differ systematically from the complete source cohort.
* The supplied subtype variable combines PAM50 and Claudin-low classification and should not be described as standard PAM50 alone.
* The Normal-like category may partly reflect low tumor cellularity or normal-tissue contamination.
* Historical pathology terminology and clinical practices differed across contributing centers and time periods.
* Overall survival includes deaths unrelated to breast cancer.
* The relapse-free-survival outcome is a composite endpoint and does not separately identify every type of event.
* Tumor stage is missing for approximately one-quarter of the eligible cohort.
* Treatment indicators do not contain sufficient information about treatment timing, regimen, adherence, or subsequent therapy to support causal treatment comparisons.
* The C-index estimates describe performance within the analyzed METABRIC cohort and may overstate performance in other patient populations.
* This project is intended to demonstrate reproducible statistical practice and must not be presented as a clinical decision tool.

## 16. Reproducibility

The analysis will be conducted entirely in R.

Reproducibility measures will include:

* Preserving the original source files without manual editing
* Recording the data source and download date
* Keeping raw data outside version control
* Joining files using documented patient identifiers
* Implementing all eligibility criteria in code
* Recording each exclusion and its reason
* Checking uniqueness of patient and sample identifiers
* Verifying survival-time and event-status coding
* Recording all recoding decisions
* Using project-relative file paths
* Recording package versions with `renv`
* Rendering the final analysis using Quarto
* Generating every table and figure directly from code
* Retaining the original ODbL attribution and license information

The repository will contain analysis code, documentation, and non-identifying derived outputs.

## 17. Planned repository structure

```text
.
├── README.md
├── analysis.qmd
├── R/
│   ├── 01_download_data.R
│   ├── 02_prepare_data.R
│   ├── 03_define_cohorts.R
│   ├── 04_descriptive_analysis.R
│   ├── 05_survival_analysis.R
│   ├── 06_model_diagnostics.R
│   └── 07_sensitivity_analysis.R
├── data-raw/
├── data-derived/
├── output/
│   ├── figures/
│   └── tables/
├── renv.lock
└── .gitignore
```

Raw data files in `data-raw/` will not be committed to GitHub. The repository will include instructions or code for obtaining the original cBioPortal data.
