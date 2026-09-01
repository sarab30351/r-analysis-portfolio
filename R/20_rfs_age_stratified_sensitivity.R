# 20_rfs_age_stratified_sensitivity.R

# Purpose: Test whether replacing the age spline with age-stratified baseline hazards materially changes the 5-to-10 year RFS subtype estimates.

# This targeted sensitivity analysis responds to the residual age signal found in Script 19. It does not replace the primary four-interval analysis and does not apply multiplicity adjustment.


# 1. Check packages and project location

required_packages <- c("survival", "splines")

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0L) {
  stop(
    "Install before continuing: ",
    paste(missing_packages, collapse = ", ")
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

if (!file.exists(file.path(project_root, "r-analysis-portfolio.Rproj"))) {
  stop("Run this script from the project root: ", project_root)
}


# 2. Load the RFS cohort and saved results

input_paths <- c(
  cohort = file.path(
    project_root,
    "data-derived",
    "metabric_rfs_cohort.rds"
  ),
  primary = file.path(
    project_root,
    "data-derived",
    "rfs_interval_cox_models.rds"
  ),
  diagnostics = file.path(
    project_root,
    "data-derived",
    "rfs_interval_diagnostics.rds"
  )
)

missing_inputs <- input_paths[!file.exists(input_paths)]

if (length(missing_inputs) > 0L) {
  stop(
    "Run Scripts 18 and 19 first. Missing: ",
    paste(missing_inputs, collapse = ", ")
  )
}

rfs_data <- as.data.frame(readRDS(input_paths["cohort"]))
primary_results <- readRDS(input_paths["primary"])
diagnostic_results <- readRDS(input_paths["diagnostics"])

required_primary_components <- c(
  "spline_knots",
  "model_overview",
  "subtype_hazard_ratios"
)

missing_primary_components <- setdiff(
  required_primary_components,
  names(primary_results)
)

if (length(missing_primary_components) > 0L) {
  stop(
    "Script 18 results are missing: ",
    paste(missing_primary_components, collapse = ", ")
  )
}

if (!"term_tests" %in% names(diagnostic_results)) {
  stop("Script 19 term-level diagnostics were not found.")
}


# 3. Prepare and validate the model variables

required_columns <- c(
  "patient_id",
  "rfs_months",
  "rfs_event",
  "age_at_diagnosis",
  "npi",
  "molecular_subtype",
  "cohort"
)

missing_columns <- setdiff(required_columns, names(rfs_data))

if (length(missing_columns) > 0L) {
  stop(
    "The RFS cohort is missing: ",
    paste(missing_columns, collapse = ", ")
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

rfs_data$molecular_subtype <- factor(
  as.character(rfs_data$molecular_subtype),
  levels = subtype_levels
)

contrasts(rfs_data$molecular_subtype) <- stats::contr.treatment(
  subtype_levels,
  base = 1
)

rfs_data$cohort <- factor(rfs_data$cohort)
rfs_data$rfs_years <- rfs_data$rfs_months / 12

model_columns <- c(
  "patient_id",
  "rfs_years",
  "rfs_event",
  "age_at_diagnosis",
  "npi",
  "molecular_subtype",
  "cohort"
)

if (anyNA(rfs_data[model_columns])) {
  stop("The RFS model variables contain missing values.")
}

if (
  nrow(rfs_data) != 1960L ||
  sum(rfs_data$rfs_event) != 790L ||
  any(rfs_data$rfs_years < 0) ||
  !all(rfs_data$rfs_event %in% c(0, 1)) ||
  levels(rfs_data$molecular_subtype)[1] != "Luminal A"
) {
  stop("The data do not match the validated RFS cohort.")
}


# 4. Recover the primary NPI spline knots

spline_knots <- primary_results$spline_knots

if (!all(c("variable", "percentile", "value") %in% names(spline_knots))) {
  stop("The saved spline-knot table has unexpected columns.")
}

npi_rows <- spline_knots$variable == "NPI"

npi_knots <- spline_knots$value[npi_rows][
  order(spline_knots$percentile[npi_rows])
]

if (length(npi_knots) != 3L || anyDuplicated(npi_knots) > 0L) {
  stop("The primary NPI spline-knot specification is invalid.")
}


# 5. Confirm the age diagnostic signal

age_signal <- diagnostic_results$term_tests[
  diagnostic_results$term_tests$interval == "5 to 10 years" &
    diagnostic_results$term_tests$term == "Age spline",
  ,
  drop = FALSE
]

if (nrow(age_signal) != 1L) {
  stop("The expected 5-to-10-year age diagnostic was not found.")
}

if (age_signal$p_value >= 0.05) {
  warning(
    paste(
      "The age diagnostic no longer matches the signal",
      "that motivated this sensitivity analysis."
    )
  )
}


# 6. Define age strata and the 5-to-10-year risk set

age_breaks <- as.numeric(
  stats::quantile(
    rfs_data$age_at_diagnosis,
    probs = seq(0, 1, by = 0.20),
    names = FALSE,
    type = 7
  )
)

if (anyDuplicated(age_breaks) > 0L) {
  stop("The age-quintile boundaries are not unique.")
}

age_levels <- paste("Age quintile", 1:5)

rfs_data$age_stratum <- cut(
  rfs_data$age_at_diagnosis,
  breaks = age_breaks,
  include.lowest = TRUE,
  labels = age_levels
)

if (anyNA(rfs_data$age_stratum)) {
  stop("At least one participant was not assigned an age stratum.")
}

interval_data <- rfs_data[
  rfs_data$rfs_years > 5,
  ,
  drop = FALSE
]

interval_data$interval_time <- pmin(interval_data$rfs_years, 10) - 5

interval_data$interval_event <- as.integer(
  interval_data$rfs_event == 1L &
    interval_data$rfs_years <= 10
)

if (
  nrow(interval_data) != 1332L ||
  sum(interval_data$interval_event) != 180L ||
  any(interval_data$interval_time <= 0)
) {
  stop("The 5-to-10-year risk set failed validation.")
}

age_stratum_summary <- data.frame(
  age_stratum = age_levels,
  lower_age_years = age_breaks[-length(age_breaks)],
  upper_age_years = age_breaks[-1],
  participants_at_risk = vapply(
    age_levels,
    function(level) {
      sum(interval_data$age_stratum == level)
    },
    integer(1)
  ),
  rfs_events = vapply(
    age_levels,
    function(level) {
      sum(
        interval_data$interval_event[
          interval_data$age_stratum == level
        ]
      )
    },
    integer(1)
  ),
  stringsAsFactors = FALSE,
  row.names = NULL
)

if (
  sum(age_stratum_summary$participants_at_risk) != nrow(interval_data) ||
  sum(age_stratum_summary$rfs_events) != sum(interval_data$interval_event)
) {
  stop("The age-stratum summary failed validation.")
}


# 7. Fit the age-stratified models

clinical_formula <- survival::Surv(
  interval_time,
  interval_event
) ~
  splines::ns(
    npi,
    knots = npi_knots[2],
    Boundary.knots = npi_knots[c(1, 3)]
  ) +
  strata(cohort, age_stratum)

extended_formula <- stats::update(
  clinical_formula,
  . ~ . + molecular_subtype
)

fit_model <- function(formula) {
  survival::coxph(
    formula = formula,
    data = interval_data,
    ties = "efron",
    na.action = stats::na.fail,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
}

clinical_model <- fit_model(clinical_formula)
extended_model <- fit_model(extended_formula)

if (
  clinical_model$n != nrow(interval_data) ||
  extended_model$n != nrow(interval_data) ||
  clinical_model$nevent != 180L ||
  extended_model$nevent != 180L
) {
  stop("A sensitivity model did not use the complete risk set.")
}

clinical_loglik <- stats::logLik(clinical_model)
extended_loglik <- stats::logLik(extended_model)

clinical_parameters <- attr(clinical_loglik, "df")
extended_parameters <- attr(extended_loglik, "df")
likelihood_ratio_df <- extended_parameters - clinical_parameters

if (
  clinical_parameters != 2L ||
  extended_parameters != 7L ||
  likelihood_ratio_df != 5L
) {
  stop("The sensitivity models have unexpected parameter counts.")
}

likelihood_ratio_chisq <- 2 * (
  as.numeric(extended_loglik) - as.numeric(clinical_loglik)
)

model_overview <- data.frame(
  analysis = "Age-stratified sensitivity",
  interval = "5 to 10 years",
  participants_at_risk = nrow(interval_data),
  rfs_events = sum(interval_data$interval_event),
  clinical_parameters = clinical_parameters,
  extended_parameters = extended_parameters,
  likelihood_ratio_chisq = likelihood_ratio_chisq,
  degrees_of_freedom = likelihood_ratio_df,
  p_value = stats::pchisq(
    likelihood_ratio_chisq,
    df = likelihood_ratio_df,
    lower.tail = FALSE
  ),
  stringsAsFactors = FALSE
)


# 8. Compare subtype estimates with the primary model

model_summary <- summary(extended_model)
coefficient_table <- model_summary$coefficients
confidence_table <- model_summary$conf.int

subtype_terms <- grep(
  "^molecular_subtype",
  rownames(coefficient_table),
  value = TRUE
)

if (length(subtype_terms) != 5L) {
  stop("The sensitivity model does not contain five subtype coefficients.")
}

comparisons <- paste(
  sub("^molecular_subtype", "", subtype_terms),
  "vs Luminal A"
)

sensitivity_estimates <- data.frame(
  comparison = comparisons,
  sensitivity_hazard_ratio = confidence_table[subtype_terms, "exp(coef)"],
  sensitivity_lower_95_ci = confidence_table[subtype_terms, "lower .95"],
  sensitivity_upper_95_ci = confidence_table[subtype_terms, "upper .95"],
  sensitivity_p_value = coefficient_table[subtype_terms, "Pr(>|z|)"],
  stringsAsFactors = FALSE,
  row.names = NULL
)

primary_estimates <- primary_results$subtype_hazard_ratios[
  primary_results$subtype_hazard_ratios$interval == "5 to 10 years",
  c("comparison", "hazard_ratio", "lower_95_ci", "upper_95_ci")
]

primary_rows <- match(
  sensitivity_estimates$comparison,
  primary_estimates$comparison
)

if (nrow(primary_estimates) != 5L || anyNA(primary_rows)) {
  stop("The primary and sensitivity estimates could not be matched.")
}

primary_estimates <- primary_estimates[primary_rows, , drop = FALSE]

subtype_comparison <- data.frame(
  analysis = "Age-stratified sensitivity",
  interval = "5 to 10 years",
  comparison = sensitivity_estimates$comparison,
  primary_hazard_ratio = primary_estimates$hazard_ratio,
  primary_lower_95_ci = primary_estimates$lower_95_ci,
  primary_upper_95_ci = primary_estimates$upper_95_ci,
  sensitivity_hazard_ratio = sensitivity_estimates$sensitivity_hazard_ratio,
  sensitivity_lower_95_ci = sensitivity_estimates$sensitivity_lower_95_ci,
  sensitivity_upper_95_ci = sensitivity_estimates$sensitivity_upper_95_ci,
  sensitivity_p_value = sensitivity_estimates$sensitivity_p_value,
  stringsAsFactors = FALSE,
  row.names = NULL
)

subtype_comparison$percent_hazard_ratio_change <- 100 * (
  subtype_comparison$sensitivity_hazard_ratio -
    subtype_comparison$primary_hazard_ratio
) / subtype_comparison$primary_hazard_ratio

subtype_comparison$absolute_percent_change <- abs(
  subtype_comparison$percent_hazard_ratio_change
)

if (anyNA(subtype_comparison)) {
  stop("The subtype comparison contains missing values.")
}


# 9. Diagnose the age-stratified model

ph_terms <- survival::cox.zph(
  extended_model,
  transform = "km",
  terms = TRUE,
  singledf = FALSE,
  global = TRUE
)

ph_coefficients <- survival::cox.zph(
  extended_model,
  transform = "km",
  terms = FALSE,
  global = TRUE
)

find_one <- function(pattern, values, label) {
  matches <- grep(pattern, values)
  
  if (length(matches) != 1L) {
    stop("Expected one diagnostic row for ", label, ".")
  }
  
  matches
}

term_matrix <- ph_terms$table
term_names <- rownames(term_matrix)

term_rows <- c(
  find_one("npi", term_names, "NPI"),
  find_one("molecular_subtype", term_names, "molecular subtype"),
  find_one("^GLOBAL$", term_names, "the global model")
)

term_tests <- data.frame(
  diagnostic_level = "Model term",
  term = c("NPI spline", "Molecular subtype", "Global model"),
  comparison = NA_character_,
  chi_square = term_matrix[term_rows, "chisq"],
  degrees_of_freedom = term_matrix[term_rows, "df"],
  p_value = term_matrix[term_rows, "p"],
  stringsAsFactors = FALSE,
  row.names = NULL
)

if (!identical(as.integer(term_tests$degrees_of_freedom), c(2L, 5L, 7L))) {
  stop("The term-level diagnostics have unexpected degrees of freedom.")
}

coefficient_matrix <- ph_coefficients$table
subtype_rows <- match(subtype_terms, rownames(coefficient_matrix))

if (anyNA(subtype_rows)) {
  stop("The subtype coefficients could not be matched to their diagnostics.")
}

subtype_tests <- data.frame(
  diagnostic_level = "Subtype coefficient",
  term = subtype_terms,
  comparison = comparisons,
  chi_square = coefficient_matrix[subtype_rows, "chisq"],
  degrees_of_freedom = coefficient_matrix[subtype_rows, "df"],
  p_value = coefficient_matrix[subtype_rows, "p"],
  stringsAsFactors = FALSE,
  row.names = NULL
)

if (!all(subtype_tests$degrees_of_freedom == 1)) {
  stop("The subtype diagnostics have unexpected degrees of freedom.")
}

ph_tests <- rbind(term_tests, subtype_tests)
ph_tests$analysis <- "Age-stratified sensitivity"
ph_tests$interval <- "5 to 10 years"

ph_tests$diagnostic_flag <- ifelse(
  ph_tests$p_value < 0.05,
  "Evidence against proportional hazards",
  "No evidence against proportional hazards"
)

ph_tests <- ph_tests[
  ,
  c(
    "analysis",
    "interval",
    "diagnostic_level",
    "term",
    "comparison",
    "chi_square",
    "degrees_of_freedom",
    "p_value",
    "diagnostic_flag"
  )
]


# 10. Display the sensitivity results

cat("\nAge-stratum composition during years 5 to 10:\n")
print(age_stratum_summary, row.names = FALSE)

cat("\nAge-stratified model overview:\n")

model_overview_print <- model_overview

model_overview_print$likelihood_ratio_chisq <- round(
  model_overview_print$likelihood_ratio_chisq,
  3
)

model_overview_print$p_value <- format.pval(
  model_overview$p_value,
  digits = 4,
  eps = 0.001
)

print(model_overview_print, row.names = FALSE)

cat("\nComparison with the primary subtype estimates:\n")

comparison_print <- subtype_comparison

comparison_print[
  c(
    "primary_hazard_ratio",
    "sensitivity_hazard_ratio",
    "percent_hazard_ratio_change",
    "absolute_percent_change"
  )
] <- round(
  comparison_print[
    c(
      "primary_hazard_ratio",
      "sensitivity_hazard_ratio",
      "percent_hazard_ratio_change",
      "absolute_percent_change"
    )
  ],
  3
)

print(
  comparison_print[
    ,
    c(
      "comparison",
      "primary_hazard_ratio",
      "sensitivity_hazard_ratio",
      "percent_hazard_ratio_change",
      "absolute_percent_change"
    )
  ],
  row.names = FALSE
)

cat("\nProportional hazards diagnostics after age stratification:\n")

ph_print <- ph_tests
ph_print$chi_square <- round(ph_print$chi_square, 3)

ph_print$p_value <- format.pval(
  ph_tests$p_value,
  digits = 4,
  eps = 0.001
)

print(ph_print, row.names = FALSE)

cat(
  paste0(
    "\nThis sensitivity model tests whether different handling of age ",
    "materially changes the 5-to-10-year subtype estimates. It does ",
    "not replace the primary model automatically.\n"
  )
)


# 11. Export and validate the results

table_directory <- file.path(project_root, "output", "tables")
data_directory <- file.path(project_root, "data-derived")

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)

tables <- list(
  rfs_age_stratum_summary = age_stratum_summary,
  rfs_age_stratified_sensitivity_model_overview = model_overview,
  rfs_age_stratified_sensitivity = subtype_comparison,
  rfs_age_stratified_sensitivity_ph_tests = ph_tests
)

table_paths <- file.path(
  table_directory,
  paste0(names(tables), ".csv")
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
    tables,
    table_paths
  )
)

results_path <- file.path(
  data_directory,
  "rfs_age_stratified_sensitivity.rds"
)

saveRDS(
  list(
    age_breaks = age_breaks,
    age_stratum_summary = age_stratum_summary,
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_overview = model_overview,
    subtype_comparison = subtype_comparison,
    ph_terms = ph_terms,
    ph_coefficients = ph_coefficients,
    ph_tests = ph_tests
  ),
  results_path
)

expected_rows <- c(
  rfs_age_stratum_summary = 5L,
  rfs_age_stratified_sensitivity_model_overview = 1L,
  rfs_age_stratified_sensitivity = 5L,
  rfs_age_stratified_sensitivity_ph_tests = 8L
)

observed_rows <- vapply(
  table_paths,
  function(path) {
    nrow(utils::read.csv(path))
  },
  integer(1)
)

if (!identical(unname(observed_rows), unname(expected_rows))) {
  stop("At least one exported table has an unexpected row count.")
}

if (!file.exists(results_path)) {
  stop("The age-stratified sensitivity RDS file was not created.")
}

cat("\nScript 20 completed successfully.\n")
cat("Sensitivity results:", normalizePath(results_path), "\n")
cat("Tables:\n")

cat(
  paste0(
    "  ",
    normalizePath(table_paths),
    collapse = "\n"
  ),
  "\n"
)