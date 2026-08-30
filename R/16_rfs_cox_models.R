# 16_rfs_cox_models.R

# Purpose: Fit full follow-up relapse-free survival Cox models.

# Models:
# 1. Clinical model: age spline + NPI spline + cohort stratification
# 2. Extended model: clinical model + molecular subtype

# The subtype hazard ratios remain provisional until proportional hazards diagnostics are completed.


# 1. Package and project checks

required_packages <- c(
  "survival",
  "splines"
)

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0L) {
  stop(
    paste(
      "Install the following packages before continuing:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages(
  library(survival)
)

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

project_file <- file.path(
  project_root,
  "r-analysis-portfolio.Rproj"
)

if (!file.exists(project_file)) {
  stop(
    paste(
      "The working directory is not the project root:",
      project_root
    )
  )
}


# 2. Load the RFS analysis cohort

rfs_data_path <- file.path(
  project_root,
  "data-derived",
  "metabric_rfs_cohort.rds"
)

if (!file.exists(rfs_data_path)) {
  stop(
    paste(
      "RFS cohort file not found:",
      rfs_data_path
    )
  )
}

rfs_model_data <- as.data.frame(
  readRDS(rfs_data_path)
)

required_columns <- c(
  "patient_id",
  "rfs_months",
  "rfs_event",
  "age_at_diagnosis",
  "npi",
  "molecular_subtype",
  "cohort"
)

missing_columns <- setdiff(
  required_columns,
  names(rfs_model_data)
)

if (length(missing_columns) > 0L) {
  stop(
    paste(
      "The following required columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# 3. Prepare and validate the model variables

subtype_levels <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

rfs_model_data$molecular_subtype <- factor(
  as.character(rfs_model_data$molecular_subtype),
  levels = subtype_levels
)

contrasts(
  rfs_model_data$molecular_subtype
) <- stats::contr.treatment(
  subtype_levels,
  base = 1
)

rfs_model_data$cohort <- factor(
  rfs_model_data$cohort
)

rfs_model_data$rfs_years <- (
  rfs_model_data$rfs_months / 12
)

model_columns <- c(
  "rfs_years",
  "rfs_event",
  "age_at_diagnosis",
  "npi",
  "molecular_subtype",
  "cohort"
)

if (anyNA(rfs_model_data[model_columns])) {
  stop("The RFS model variables contain missing values.")
}

if (
  !all(
    rfs_model_data$rfs_event %in% c(0, 1)
  )
) {
  stop("RFS event values must be coded as zero or one.")
}

if (any(rfs_model_data$rfs_years < 0)) {
  stop("RFS follow-up time cannot be negative.")
}

if (nrow(rfs_model_data) != 1960L) {
  stop(
    paste(
      "Expected 1,960 participants but found",
      nrow(rfs_model_data)
    )
  )
}

if (sum(rfs_model_data$rfs_event) != 790L) {
  stop(
    paste(
      "Expected 790 RFS events but found",
      sum(rfs_model_data$rfs_event)
    )
  )
}

month_zero_event_count <- sum(
  rfs_model_data$rfs_years == 0 &
    rfs_model_data$rfs_event == 1
)

if (month_zero_event_count != 3L) {
  stop(
    paste(
      "Expected three month-zero RFS events but found",
      month_zero_event_count
    )
  )
}

subtype_counts <- table(
  rfs_model_data$molecular_subtype
)

if (any(subtype_counts == 0L)) {
  stop("At least one molecular subtype has no participants.")
}

if (
  levels(rfs_model_data$molecular_subtype)[1] !=
  "Luminal A"
) {
  stop("Luminal A is not the molecular subtype reference level.")
}


# 4. Define the spline knots

# Age uses four knots at the 5th, 35th, 65th, and 95th percentiles.
# The 35th and 65th percentiles are the internal knots.

# NPI uses three knots at the 5th, 50th, and 95th percentiles.
# The 50th percentile is the internal knot.

# This matches the spline structure used in the OS models while calculating the numerical knot values from the RFS analytical cohort.

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
    rfs_model_data$age_at_diagnosis,
    probs = age_knot_probabilities,
    type = 7
  )
)

npi_knots <- unname(
  stats::quantile(
    rfs_model_data$npi,
    probs = npi_knot_probabilities,
    type = 7
  )
)

if (
  anyDuplicated(age_knots) > 0L ||
  anyDuplicated(npi_knots) > 0L
) {
  stop("The calculated spline knots are not unique.")
}

rfs_spline_knots <- data.frame(
  variable = c(
    rep("Age at diagnosis", 4),
    rep("NPI", 3)
  ),
  knot_role = c(
    "Lower boundary",
    "Internal knot 1",
    "Internal knot 2",
    "Upper boundary",
    "Lower boundary",
    "Internal knot",
    "Upper boundary"
  ),
  percentile = c(
    age_knot_probabilities,
    npi_knot_probabilities
  ),
  value = c(
    age_knots,
    npi_knots
  ),
  stringsAsFactors = FALSE
)


# 5. Define the Cox model formulas

rfs_clinical_formula <- survival::Surv(
  rfs_years,
  rfs_event
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
  strata(cohort)

rfs_extended_formula <- stats::update(
  rfs_clinical_formula,
  . ~ . + molecular_subtype
)


# 6. Fit the clinical and extended models

rfs_clinical_model <- survival::coxph(
  formula = rfs_clinical_formula,
  data = rfs_model_data,
  ties = "efron",
  na.action = stats::na.fail,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

rfs_extended_model <- survival::coxph(
  formula = rfs_extended_formula,
  data = rfs_model_data,
  ties = "efron",
  na.action = stats::na.fail,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

if (
  rfs_clinical_model$n != 1960L ||
  rfs_extended_model$n != 1960L
) {
  stop("The Cox models did not use all 1,960 participants.")
}

if (
  rfs_clinical_model$nevent != 790L ||
  rfs_extended_model$nevent != 790L
) {
  stop("The Cox models did not use all 790 RFS events.")
}


# 7. Compare the models using a likelihood-ratio test

clinical_loglik_object <- stats::logLik(
  rfs_clinical_model
)

extended_loglik_object <- stats::logLik(
  rfs_extended_model
)

clinical_loglik <- as.numeric(
  clinical_loglik_object
)

extended_loglik <- as.numeric(
  extended_loglik_object
)

clinical_parameter_count <- attr(
  clinical_loglik_object,
  "df"
)

extended_parameter_count <- attr(
  extended_loglik_object,
  "df"
)

rfs_lrt_chisq <- 2 * (
  extended_loglik -
    clinical_loglik
)

rfs_lrt_df <- (
  extended_parameter_count -
    clinical_parameter_count
)

rfs_lrt_p_value <- stats::pchisq(
  rfs_lrt_chisq,
  df = rfs_lrt_df,
  lower.tail = FALSE
)

if (rfs_lrt_df != 5L) {
  stop(
    paste(
      "Expected five degrees of freedom for subtype but found",
      rfs_lrt_df
    )
  )
}

rfs_global_subtype_test <- data.frame(
  endpoint = "Relapse-free survival",
  comparison = paste(
    "Extended molecular subtype model",
    "versus clinical model"
  ),
  participants = rfs_extended_model$n,
  rfs_events = rfs_extended_model$nevent,
  clinical_log_likelihood = clinical_loglik,
  extended_log_likelihood = extended_loglik,
  likelihood_ratio_chisq = rfs_lrt_chisq,
  degrees_of_freedom = rfs_lrt_df,
  p_value = rfs_lrt_p_value,
  stringsAsFactors = FALSE
)


# 8. Extract the subtype hazard ratios

rfs_extended_summary <- summary(
  rfs_extended_model
)

coefficient_table <- (
  rfs_extended_summary$coefficients
)

confidence_table <- (
  rfs_extended_summary$conf.int
)

subtype_terms <- grep(
  "^molecular_subtype",
  rownames(coefficient_table),
  value = TRUE
)

if (length(subtype_terms) != 5L) {
  stop(
    paste(
      "Expected five subtype coefficients but found",
      length(subtype_terms)
    )
  )
}

comparison_subtypes <- sub(
  "^molecular_subtype",
  "",
  subtype_terms
)

rfs_subtype_hazard_ratios <- data.frame(
  term = subtype_terms,
  comparison = paste(
    comparison_subtypes,
    "vs Luminal A"
  ),
  log_hazard_ratio = coefficient_table[
    subtype_terms,
    "coef"
  ],
  standard_error = coefficient_table[
    subtype_terms,
    "se(coef)"
  ],
  hazard_ratio = confidence_table[
    subtype_terms,
    "exp(coef)"
  ],
  lower_95_ci = confidence_table[
    subtype_terms,
    "lower .95"
  ],
  upper_95_ci = confidence_table[
    subtype_terms,
    "upper .95"
  ],
  z_statistic = coefficient_table[
    subtype_terms,
    "z"
  ],
  p_value = coefficient_table[
    subtype_terms,
    "Pr(>|z|)"
  ],
  stringsAsFactors = FALSE,
  row.names = NULL
)


# 9. Display the provisional results

rfs_model_overview <- data.frame(
  model = c(
    "Clinical model",
    "Extended molecular subtype model"
  ),
  participants = c(
    rfs_clinical_model$n,
    rfs_extended_model$n
  ),
  rfs_events = c(
    rfs_clinical_model$nevent,
    rfs_extended_model$nevent
  ),
  estimated_parameters = c(
    clinical_parameter_count,
    extended_parameter_count
  ),
  fitted_log_likelihood = c(
    clinical_loglik,
    extended_loglik
  ),
  stringsAsFactors = FALSE
)

cat("\nRFS spline knot specification:\n")

rfs_spline_knots_print <- rfs_spline_knots

rfs_spline_knots_print$value <- round(
  rfs_spline_knots_print$value,
  3
)

print(
  rfs_spline_knots_print,
  row.names = FALSE
)

cat("\nRFS Cox model overview:\n")

rfs_model_overview_print <- rfs_model_overview

rfs_model_overview_print$fitted_log_likelihood <- round(
  rfs_model_overview_print$fitted_log_likelihood,
  3
)

print(
  rfs_model_overview_print,
  row.names = FALSE
)

cat("\nGlobal likelihood ratio test for molecular subtype:\n")

rfs_global_test_print <- rfs_global_subtype_test

rfs_global_test_print[
  c(
    "clinical_log_likelihood",
    "extended_log_likelihood",
    "likelihood_ratio_chisq"
  )
] <- round(
  rfs_global_test_print[
    c(
      "clinical_log_likelihood",
      "extended_log_likelihood",
      "likelihood_ratio_chisq"
    )
  ],
  3
)

rfs_global_test_print$p_value <- format.pval(
  rfs_global_subtype_test$p_value,
  digits = 4,
  eps = 0.001
)

print(
  rfs_global_test_print,
  row.names = FALSE
)

cat(
  paste0(
    "\nProvisional full follow-up subtype hazard ratios ",
    "before proportional hazards testing:\n"
  )
)

rfs_hr_print <- rfs_subtype_hazard_ratios

rfs_hr_print[
  c(
    "log_hazard_ratio",
    "standard_error",
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci",
    "z_statistic"
  )
] <- round(
  rfs_hr_print[
    c(
      "log_hazard_ratio",
      "standard_error",
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci",
      "z_statistic"
    )
  ],
  3
)

rfs_hr_print$p_value <- format.pval(
  rfs_subtype_hazard_ratios$p_value,
  digits = 4,
  eps = 0.001
)

print(
  rfs_hr_print[
    ,
    c(
      "comparison",
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci",
      "p_value"
    )
  ],
  row.names = FALSE
)

cat(
  paste0(
    "\nThese full follow-up hazard ratios are provisional. ",
    "Do not interpret them as constant effects until the ",
    "proportional hazards diagnostics are evaluated.\n"
  )
)


# 10. Export the full follow-up model results

rfs_table_directory <- file.path(
  project_root,
  "output",
  "tables"
)

rfs_data_directory <- file.path(
  project_root,
  "data-derived"
)

dir.create(
  rfs_table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  rfs_data_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

rfs_full_followup_tables <- list(
  rfs_full_followup_spline_knots = rfs_spline_knots,
  rfs_full_followup_model_overview = rfs_model_overview,
  rfs_full_followup_global_subtype_test = rfs_global_subtype_test,
  rfs_full_followup_subtype_hazard_ratios =
    rfs_subtype_hazard_ratios
)

rfs_full_followup_table_paths <- file.path(
  rfs_table_directory,
  paste0(
    names(rfs_full_followup_tables),
    ".csv"
  )
)

invisible(
  Map(
    function(data, path) {
      utils::write.csv(
        data,
        path,
        row.names = FALSE,
        na = ""
      )
    },
    rfs_full_followup_tables,
    rfs_full_followup_table_paths
  )
)

rfs_full_followup_model_results <- list(
  clinical_model = rfs_clinical_model,
  extended_model = rfs_extended_model,
  clinical_formula = rfs_clinical_formula,
  extended_formula = rfs_extended_formula,
  spline_knots = rfs_spline_knots,
  model_overview = rfs_model_overview,
  global_subtype_test = rfs_global_subtype_test,
  subtype_hazard_ratios = rfs_subtype_hazard_ratios
)

rfs_full_followup_results_path <- file.path(
  rfs_data_directory,
  "rfs_full_followup_cox_models.rds"
)

saveRDS(
  rfs_full_followup_model_results,
  rfs_full_followup_results_path
)


# 11. Validate the exported model results

expected_table_rows <- c(
  rfs_full_followup_spline_knots = 7L,
  rfs_full_followup_model_overview = 2L,
  rfs_full_followup_global_subtype_test = 1L,
  rfs_full_followup_subtype_hazard_ratios = 5L
)

observed_table_rows <- vapply(
  rfs_full_followup_table_paths,
  function(path) {
    nrow(
      utils::read.csv(path)
    )
  },
  FUN.VALUE = integer(1)
)

if (
  !identical(
    unname(observed_table_rows),
    unname(expected_table_rows)
  )
) {
  stop("One or more exported tables have unexpected row counts.")
}

cat("\nScript 16 completed successfully.\n")

cat(
  "Model results:",
  rfs_full_followup_results_path,
  "\n"
)

cat("Tables:\n")

cat(
  paste0(
    "  ",
    rfs_full_followup_table_paths,
    collapse = "\n"
  ),
  "\n"
)