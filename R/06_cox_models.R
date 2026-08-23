
# 06_cox_models.R
#
# Purpose: Fit and compare clinical and extended Cox models for overall survival.


# Load required packages

suppressPackageStartupMessages({
  library(survival)
  library(splines)
})


# Import the overall-survival cohort

os_file <- here::here(
  "data-derived",
  "metabric_os_cohort.rds"
)

if (!file.exists(os_file)) {
  stop(
    "The overall-survival cohort is missing. ",
    "Run R/03_define_cohorts.R first."
  )
}

subtype_levels <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

os_model_data <- readRDS(os_file) |>
  dplyr::mutate(
    os_years = os_months / 12,
    molecular_subtype = factor(
      molecular_subtype,
      levels = subtype_levels
    ),
    cohort = factor(cohort)
  )


# Validate the model data

required_model_variables <- c(
  "os_years",
  "os_event",
  "age_at_diagnosis",
  "npi",
  "cohort",
  "molecular_subtype"
)

missing_variables <- setdiff(
  required_model_variables,
  names(os_model_data)
)

if (length(missing_variables) > 0) {
  stop(
    "Required model variables are missing: ",
    paste(missing_variables, collapse = ", ")
  )
}

if (anyNA(os_model_data[, required_model_variables])) {
  stop(
    "The Cox-model dataset contains missing values in required variables."
  )
}

if (!all(os_model_data$os_event %in% c(0L, 1L))) {
  stop("Overall-survival event values must be coded as 0 or 1.")
}


# Define spline knots

age_knot_probabilities <- c(
  0.05,
  0.35,
  0.65,
  0.95
)

npi_knot_probabilities <- c(
  0.05,
  0.50,
  0.95
)

age_knots <- unname(
  stats::quantile(
    os_model_data$age_at_diagnosis,
    probs = age_knot_probabilities
  )
)

npi_knots <- unname(
  stats::quantile(
    os_model_data$npi,
    probs = npi_knot_probabilities
  )
)

if (any(diff(age_knots) <= 0)) {
  stop("The age spline knots are not distinct.")
}

if (any(diff(npi_knots) <= 0)) {
  stop("The NPI spline knots are not distinct.")
}


# Fit the clinical Cox model

clinical_model <- survival::coxph(
  survival::Surv(
    os_years,
    os_event
  ) ~
    splines::ns(
      age_at_diagnosis,
      knots = age_knots[c(2, 3)],
      Boundary.knots = age_knots[c(1, 4)]
    ) +
    splines::ns(
      npi,
      knots = npi_knots[2],
      Boundary.knots = npi_knots[c(1, 3)]
    ) +
    strata(cohort),
  data = os_model_data,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)


# Fit the extended Cox model

extended_model <- survival::coxph(
  survival::Surv(
    os_years,
    os_event
  ) ~
    splines::ns(
      age_at_diagnosis,
      knots = age_knots[c(2, 3)],
      Boundary.knots = age_knots[c(1, 4)]
    ) +
    splines::ns(
      npi,
      knots = npi_knots[2],
      Boundary.knots = npi_knots[c(1, 3)]
    ) +
    molecular_subtype +
    strata(cohort),
  data = os_model_data,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)


# Confirm that both models used the same observations

if (
  clinical_model$n != nrow(os_model_data) ||
  extended_model$n != nrow(os_model_data)
) {
  stop("One or both Cox models did not use the complete model cohort.")
}


# Compare the clinical and extended models

model_comparison <- stats::anova(
  clinical_model,
  extended_model,
  test = "LRT"
)


# Extract adjusted molecular-subtype estimates

extended_model_summary <- summary(extended_model)

subtype_rows <- grepl(
  "^molecular_subtype",
  rownames(extended_model_summary$coefficients)
)

subtype_coefficients <- extended_model_summary$coefficients[
  subtype_rows,
  ,
  drop = FALSE
]

subtype_confidence_intervals <- extended_model_summary$conf.int[
  subtype_rows,
  ,
  drop = FALSE
]

subtype_results <- data.frame(
  comparison = sub(
    "^molecular_subtype",
    "",
    rownames(subtype_coefficients)
  ),
  hazard_ratio = subtype_confidence_intervals[, "exp(coef)"],
  lower_95_ci = subtype_confidence_intervals[, "lower .95"],
  upper_95_ci = subtype_confidence_intervals[, "upper .95"],
  p_value = subtype_coefficients[, "Pr(>|z|)"],
  row.names = NULL
)


# Save model objects for later diagnostics

model_bundle <- list(
  clinical_model = clinical_model,
  extended_model = extended_model,
  age_knots = age_knots,
  npi_knots = npi_knots,
  model_participants = nrow(os_model_data),
  model_events = sum(os_model_data$os_event)
)

saveRDS(
  model_bundle,
  here::here(
    "data-derived",
    "os_cox_models.rds"
  )
)


# Display initial results

print(model_comparison)
print(subtype_results)

message(
  "Cox models fitted using ",
  nrow(os_model_data),
  " participants and ",
  sum(os_model_data$os_event),
  " deaths."
)


### One adjustment is the supposed 4-knot spline for NPI would place its two inner knots at 4.016 and 4.080, which is almost the same location. I instead used these:
### Age: 4 knots at the 5th, 35th, 65th, and 95th percentiles.
### NPI: 3 knots at the 5th, 50th, and 95th percentiles.
### Cohort stratification 

### Warning: Basal-like, HER2-enriched and Luminal B subtypes show strong violations. 
### Basal- like: (p ≈ 2 × 10⁻¹⁶) Its relative hazard is high early and declines sharply.
### HER2-enriched: (p = 0.00033). Its disadvantage is concentrated earlier.
### Luminal B: (p = 0.021) weaker, but still a violation.

### Thinking of applying interval-specific models, 0–5 years, 5–10 years, conditional on surviving beyond 5, beyond 10 years, conditional on surviving beyond 10 (see documented 07_interval_cox_models.R)