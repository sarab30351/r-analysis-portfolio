
# 11_npi_stratified_sensitivity.R
# Purpose: NPI showed evidence of a changing association with mortality during years 2 to 5. This sensitivity analysis stratifies the baseline hazard by established NPI prognostic group: excellent ≤2.4, good >2.4–3.4, moderate >3.4–5.4, and poor >5.4.
#
# Age remains adjusted using the same natural cubic spline as in the primary model. METABRIC source cohort and NPI group are used as stratification variables.


# Check packages and project location 

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' must be installed.")
}

if (!requireNamespace("splines", quietly = TRUE)) {
  stop("Package 'splines' must be installed.")
}

suppressPackageStartupMessages(
  library(survival)
)

if (!file.exists("r-analysis-portfolio.Rproj")) {
  stop(
    paste(
      "Open r-analysis-portfolio.Rproj before running this script.",
      "The working directory must be the repository root."
    )
  )
}


# Load the OS cohort and primary model results 

os_cohort_path <- file.path(
  "data-derived",
  "metabric_os_cohort.rds"
)

primary_model_path <- file.path(
  "data-derived",
  "os_refined_interval_cox_models.rds"
)

if (!file.exists(os_cohort_path)) {
  stop("The overall-survival cohort could not be found.")
}

if (!file.exists(primary_model_path)) {
  stop("Run Script 08 before running this sensitivity analysis.")
}

os_model_data <- readRDS(os_cohort_path)
primary_results <- readRDS(primary_model_path)


# Create years of follow-up if needed 

if (!"os_years" %in% names(os_model_data)) {
  
  if (!"os_months" %in% names(os_model_data)) {
    stop("Neither os_years nor os_months was found.")
  }
  
  os_model_data$os_years <- os_model_data$os_months / 12
}


# Validate the required variables 
required_variables <- c(
  "os_years",
  "os_event",
  "age_at_diagnosis",
  "npi",
  "cohort",
  "molecular_subtype"
)

missing_variables <- setdiff(
  required_variables,
  names(os_model_data)
)

if (length(missing_variables) > 0) {
  stop(
    paste(
      "Missing required variable(s):",
      paste(missing_variables, collapse = ", ")
    )
  )
}

if (!all(stats::complete.cases(
  os_model_data[, required_variables]
))) {
  stop("At least one model variable contains missing values.")
}


# Standardize factor levels 

subtype_levels <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

os_model_data$molecular_subtype <- factor(
  os_model_data$molecular_subtype,
  levels = subtype_levels
)

os_model_data$cohort <- factor(
  os_model_data$cohort
)


# Define established NPI prognostic groups

os_model_data$npi_group <- cut(
  os_model_data$npi,
  breaks = c(
    -Inf,
    2.4,
    3.4,
    5.4,
    Inf
  ),
  labels = c(
    "Excellent: NPI <= 2.4",
    "Good: NPI > 2.4 to 3.4",
    "Moderate: NPI > 3.4 to 5.4",
    "Poor: NPI > 5.4"
  ),
  right = TRUE
)

if (any(is.na(os_model_data$npi_group))) {
  stop("At least one participant could not be assigned an NPI group.")
}


# Recreate the age spline knots 

age_knots <- as.numeric(
  stats::quantile(
    os_model_data$age_at_diagnosis,
    probs = c(0.05, 0.35, 0.65, 0.95),
    names = FALSE
  )
)


# Construct the 2-to-5-year risk set 

two_to_five_data <- os_model_data[
  os_model_data$os_years > 2,
  ,
  drop = FALSE
]

two_to_five_data$entry_time <- 2

two_to_five_data$exit_time <- pmin(
  two_to_five_data$os_years,
  5
)

two_to_five_data$interval_event <- as.integer(
  two_to_five_data$os_event == 1L &
    two_to_five_data$os_years > 2 &
    two_to_five_data$os_years <= 5
)

if (any(
  two_to_five_data$exit_time <=
  two_to_five_data$entry_time
)) {
  stop("Invalid follow-up time found in the 2-to-5-year risk set.")
}


# Summarize the NPI groups in this risk set 

npi_levels <- levels(
  two_to_five_data$npi_group
)

npi_group_summary <- data.frame(
  npi_group = npi_levels,
  
  participants = vapply(
    npi_levels,
    function(group) {
      sum(two_to_five_data$npi_group == group)
    },
    FUN.VALUE = numeric(1)
  ),
  
  deaths = vapply(
    npi_levels,
    function(group) {
      sum(
        two_to_five_data$interval_event[
          two_to_five_data$npi_group == group
        ]
      )
    },
    FUN.VALUE = numeric(1)
  )
)


# Fit the NPI-stratified models 

clinical_model <- survival::coxph(
  survival::Surv(
    entry_time,
    exit_time,
    interval_event
  ) ~
    splines::ns(
      age_at_diagnosis,
      knots = age_knots[c(2, 3)],
      Boundary.knots = age_knots[c(1, 4)]
    ) +
    strata(cohort, npi_group),
  data = two_to_five_data,
  ties = "efron",
  model = TRUE,
  x = TRUE,
  y = TRUE
)

extended_model <- survival::coxph(
  survival::Surv(
    entry_time,
    exit_time,
    interval_event
  ) ~
    splines::ns(
      age_at_diagnosis,
      knots = age_knots[c(2, 3)],
      Boundary.knots = age_knots[c(1, 4)]
    ) +
    molecular_subtype +
    strata(cohort, npi_group),
  data = two_to_five_data,
  ties = "efron",
  model = TRUE,
  x = TRUE,
  y = TRUE
)


# Test the added contribution of molecular subtype 

model_comparison <- stats::anova(
  clinical_model,
  extended_model,
  test = "LRT"
)

model_comparison_summary <- data.frame(
  interval = "2 to 5 years",
  participants = nrow(two_to_five_data),
  deaths = sum(two_to_five_data$interval_event),
  likelihood_ratio_chisq =
    model_comparison[2, "Chisq"],
  degrees_of_freedom =
    model_comparison[2, "Df"],
  p_value =
    model_comparison[2, "Pr(>|Chi|)"]
)


# Extract subtype estimates 

coefficient_table <- summary(
  extended_model
)$coefficients

confidence_intervals <- stats::confint(
  extended_model
)

subtype_rows <- grep(
  "^molecular_subtype",
  rownames(coefficient_table)
)

subtype_names <- rownames(
  coefficient_table
)[subtype_rows]

sensitivity_estimates <- data.frame(
  interval = "2 to 5 years",
  
  comparison = paste(
    sub(
      "^molecular_subtype",
      "",
      subtype_names
    ),
    "vs Luminal A"
  ),
  
  sensitivity_hazard_ratio = exp(
    coefficient_table[subtype_rows, "coef"]
  ),
  
  sensitivity_lower_95_ci = exp(
    confidence_intervals[subtype_names, 1]
  ),
  
  sensitivity_upper_95_ci = exp(
    confidence_intervals[subtype_names, 2]
  ),
  
  sensitivity_p_value =
    coefficient_table[subtype_rows, "Pr(>|z|)"],
  
  row.names = NULL
)


# Compare against the primary 2-to-5-year model 

primary_estimates <- primary_results$subtype_estimates

primary_estimates <- primary_estimates[
  primary_estimates$interval == "2 to 5 years",
  c(
    "interval",
    "comparison",
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci"
  )
]

names(primary_estimates)[3:5] <- c(
  "primary_hazard_ratio",
  "primary_lower_95_ci",
  "primary_upper_95_ci"
)

hazard_ratio_comparison <- merge(
  primary_estimates,
  sensitivity_estimates,
  by = c(
    "interval",
    "comparison"
  ),
  all = TRUE,
  sort = FALSE
)

hazard_ratio_comparison$percent_change <- 100 * (
  hazard_ratio_comparison$sensitivity_hazard_ratio -
    hazard_ratio_comparison$primary_hazard_ratio
) /
  hazard_ratio_comparison$primary_hazard_ratio


# Recheck proportional hazards 

ph_terms <- survival::cox.zph(
  extended_model,
  transform = "km",
  terms = TRUE,
  global = TRUE
)

ph_term_table <- data.frame(
  term = rownames(ph_terms$table),
  ph_terms$table,
  row.names = NULL,
  check.names = FALSE
)

ph_coefficients <- survival::cox.zph(
  extended_model,
  transform = "km",
  terms = FALSE,
  global = TRUE
)

subtype_ph_rows <- grepl(
  "^molecular_subtype",
  rownames(ph_coefficients$table)
)

subtype_ph_table <- data.frame(
  coefficient = rownames(
    ph_coefficients$table
  )[subtype_ph_rows],
  
  ph_coefficients$table[
    subtype_ph_rows,
    ,
    drop = FALSE
  ],
  
  row.names = NULL,
  check.names = FALSE
)


# Print results 

cat("\nNPI-group composition during years 2 to 5:\n")
print(
  npi_group_summary,
  row.names = FALSE
)

cat("\nLikelihood-ratio comparison:\n")
print(model_comparison)

cat("\nModel-comparison summary:\n")
print(
  model_comparison_summary,
  row.names = FALSE
)

cat("\nComparison with the primary subtype estimates:\n")
print(
  hazard_ratio_comparison,
  row.names = FALSE
)

cat("\nTerm-level proportional-hazards tests:\n")
print(
  ph_term_table,
  row.names = FALSE
)

cat("\nSubtype coefficient-level PH tests:\n")
print(
  subtype_ph_table,
  row.names = FALSE
)


# Save sensitivity results locally 

sensitivity_path <- file.path(
  "data-derived",
  "os_npi_stratified_sensitivity.rds"
)

saveRDS(
  list(
    npi_group_summary = npi_group_summary,
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_comparison = model_comparison_summary,
    sensitivity_estimates = sensitivity_estimates,
    hazard_ratio_comparison = hazard_ratio_comparison,
    ph_terms = ph_term_table,
    subtype_ph_results = subtype_ph_table
  ),
  sensitivity_path
)

cat(
  "\nNPI-stratified sensitivity results saved locally to: ",
  normalizePath(sensitivity_path),
  "\n",
  sep = ""
)

### This sensitivity analysis resolved the remaining concerns. I now feel more confident in these primary findings:
### Subtype provides substantial added prognostic information during 0-2 and 2-5 years.
### There is no clear overall added contribution during 5-10 years.
### Subtype differences reappear beyond 10 years, primarily through the lower late hazard among the remaining Basal-like risk set.
### Sensitivity analyses using age group strata and NPI prognostic group strata produced similar subtype hazard-ratio estimates, with absolute changes below 10%. (see documented 08_refined_interval_cox_models.R) for primary model results.