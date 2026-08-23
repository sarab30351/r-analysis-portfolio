
# 10_age_stratified_sensitivity.R
# Purpose:
# The refined interval models showed that the association between age and mortality changes during follow-up. In this sensitivity analysis,
# participants are divided into five equally sized age groups.

# Age group and source cohort are used as stratification variables, allowing separate baseline hazards without estimating a constant hazard ratio for age.

# This is a sensitivity analysis, primary model is in script 08_refined_interval_cox_models.R.


# Check packages and project location 

suppressPackageStartupMessages(
  library(survival)
)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' must be installed.")
}

if (!requireNamespace("splines", quietly = TRUE)) {
  stop("Package 'splines' must be installed.")
}

if (!file.exists("r-analysis-portfolio.Rproj")) {
  stop(
    paste(
      "Open r-analysis-portfolio.Rproj before running this script.",
      "The working directory must be the repository root."
    )
  )
}


# Load data and primary results 
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


# Validate variables 

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


# Define five balanced age strata 
age_stratum_breaks <- as.numeric(
  stats::quantile(
    os_model_data$age_at_diagnosis,
    probs = seq(0, 1, by = 0.20),
    names = FALSE
  )
)

if (anyDuplicated(age_stratum_breaks)) {
  stop(
    paste(
      "Age quintile boundaries are not unique.",
      "Age strata cannot be constructed reliably."
    )
  )
}

os_model_data$age_stratum <- cut(
  os_model_data$age_at_diagnosis,
  breaks = age_stratum_breaks,
  include.lowest = TRUE,
  labels = paste(
    "Age quintile",
    1:5
  )
)

if (any(is.na(os_model_data$age_stratum))) {
  stop("At least one participant could not be assigned an age stratum.")
}


# Recreate the NPI spline knots 

npi_knots <- as.numeric(
  stats::quantile(
    os_model_data$npi,
    probs = c(0.05, 0.50, 0.95),
    names = FALSE
  )
)


# Create an interval-specific risk set 

create_interval_data <- function(
    data,
    start_time,
    end_time,
    interval_label
) {
  
  interval_data <- data[
    data$os_years > start_time,
    ,
    drop = FALSE
  ]
  
  interval_data$entry_time <- start_time
  
  if (is.infinite(end_time)) {
    interval_data$exit_time <- interval_data$os_years
  } else {
    interval_data$exit_time <- pmin(
      interval_data$os_years,
      end_time
    )
  }
  
  interval_data$interval_event <- as.integer(
    interval_data$os_event == 1L &
      interval_data$os_years > start_time &
      (
        is.infinite(end_time) |
          interval_data$os_years <= end_time
      )
  )
  
  interval_data$interval_label <- interval_label
  
  interval_data
}


# Fit one age-stratified interval model 

fit_age_stratified_models <- function(
    interval_data,
    interval_label
) {
  
  clinical_model <- survival::coxph(
    survival::Surv(
      entry_time,
      exit_time,
      interval_event
    ) ~
      splines::ns(
        npi,
        knots = npi_knots[2],
        Boundary.knots = npi_knots[c(1, 3)]
      ) +
      strata(cohort, age_stratum),
    data = interval_data,
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
        npi,
        knots = npi_knots[2],
        Boundary.knots = npi_knots[c(1, 3)]
      ) +
      molecular_subtype +
      strata(cohort, age_stratum),
    data = interval_data,
    ties = "efron",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  model_comparison <- stats::anova(
    clinical_model,
    extended_model,
    test = "LRT"
  )
  
  comparison_summary <- data.frame(
    interval = interval_label,
    participants = nrow(interval_data),
    deaths = sum(interval_data$interval_event),
    likelihood_ratio_chisq =
      model_comparison[2, "Chisq"],
    degrees_of_freedom =
      model_comparison[2, "Df"],
    p_value =
      model_comparison[2, "Pr(>|Chi|)"]
  )
  
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
  
  subtype_estimates <- data.frame(
    interval = interval_label,
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
  
  ph_test <- survival::cox.zph(
    extended_model,
    transform = "km",
    terms = TRUE,
    global = TRUE
  )
  
  ph_table <- data.frame(
    interval = interval_label,
    term = rownames(ph_test$table),
    ph_test$table,
    row.names = NULL,
    check.names = FALSE
  )
  
  list(
    interval = interval_label,
    clinical_model = clinical_model,
    extended_model = extended_model,
    comparison_summary = comparison_summary,
    subtype_estimates = subtype_estimates,
    ph_test = ph_test,
    ph_table = ph_table
  )
}


# Construct the same four intervals used in script 08_refined_interval_cox_models.R 

sensitivity_interval_data <- list(
  zero_to_two = create_interval_data(
    os_model_data,
    start_time = 0,
    end_time = 2,
    interval_label = "0 to 2 years"
  ),
  
  two_to_five = create_interval_data(
    os_model_data,
    start_time = 2,
    end_time = 5,
    interval_label = "2 to 5 years"
  ),
  
  five_to_ten = create_interval_data(
    os_model_data,
    start_time = 5,
    end_time = 10,
    interval_label = "5 to 10 years"
  ),
  
  beyond_ten = create_interval_data(
    os_model_data,
    start_time = 10,
    end_time = Inf,
    interval_label = "Beyond 10 years"
  )
)


# Fit all sensitivity models 

age_stratified_results <- Map(
  f = function(interval_data, interval_name) {
    fit_age_stratified_models(
      interval_data = interval_data,
      interval_label = interval_name
    )
  },
  interval_data = sensitivity_interval_data,
  interval_name = vapply(
    sensitivity_interval_data,
    function(data) {
      unique(data$interval_label)
    },
    FUN.VALUE = character(1)
  )
)


# Combine results 

sensitivity_comparisons <- do.call(
  rbind,
  lapply(
    age_stratified_results,
    function(result) {
      result$comparison_summary
    }
  )
)

sensitivity_subtype_estimates <- do.call(
  rbind,
  lapply(
    age_stratified_results,
    function(result) {
      result$subtype_estimates
    }
  )
)

sensitivity_ph_results <- do.call(
  rbind,
  lapply(
    age_stratified_results,
    function(result) {
      result$ph_table
    }
  )
)


# Compare sensitivity and primary hazard ratios 

primary_subtype_estimates <-
  primary_results$subtype_estimates

primary_subtype_estimates <- primary_subtype_estimates[
  ,
  c(
    "interval",
    "comparison",
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci"
  )
]

names(primary_subtype_estimates)[3:5] <- c(
  "primary_hazard_ratio",
  "primary_lower_95_ci",
  "primary_upper_95_ci"
)

hazard_ratio_comparison <- merge(
  primary_subtype_estimates,
  sensitivity_subtype_estimates,
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


# Print age-stratum boundaries and results 

cat("\nAge-stratum boundaries:\n")

print(
  data.frame(
    boundary = c(
      "Minimum",
      "20th percentile",
      "40th percentile",
      "60th percentile",
      "80th percentile",
      "Maximum"
    ),
    age_years = age_stratum_breaks
  ),
  row.names = FALSE
)

cat("\nAge-stratified model comparisons:\n")

print(
  sensitivity_comparisons,
  row.names = FALSE
)

cat("\nComparison with the primary subtype estimates:\n")

print(
  hazard_ratio_comparison,
  row.names = FALSE
)

cat("\nProportional-hazards tests after age stratification:\n")

print(
  sensitivity_ph_results,
  row.names = FALSE
)


# Save sensitivity results locally 

sensitivity_model_path <- file.path(
  "data-derived",
  "os_age_stratified_sensitivity.rds"
)

saveRDS(
  list(
    age_stratum_breaks = age_stratum_breaks,
    interval_results = age_stratified_results,
    model_comparisons = sensitivity_comparisons,
    subtype_estimates = sensitivity_subtype_estimates,
    hazard_ratio_comparison = hazard_ratio_comparison,
    proportional_hazards_results = sensitivity_ph_results
  ),
  sensitivity_model_path
)

cat(
  "\nSensitivity-analysis results saved locally to: ",
  normalizePath(sensitivity_model_path),
  "\n",
  sep = ""
)

### It is now evident the changing effect of age was not driving subtype findings.
### It's worth investigating into NPI as well.