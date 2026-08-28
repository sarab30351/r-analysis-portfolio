# 13_os_additional_sensitivity.R
#
# Purpose:
# Assess whether the primary overall survival results are sensitive to excluding Stage 0 or Stage IV disease and to excluding the Normal-like subtype.


# Check packages and project location

required_packages <- c(
  "survival",
  "splines"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

if (!file.exists("r-analysis-portfolio.Rproj")) {
  stop(
    "Open r-analysis-portfolio.Rproj before running this script."
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
  stop("Run Script 03 before this script.")
}

if (!file.exists(primary_model_path)) {
  stop("Run Script 08 before this script.")
}

os_model_data <- readRDS(os_cohort_path)
primary_results <- readRDS(primary_model_path)

if (!"os_years" %in% names(os_model_data)) {
  os_model_data$os_years <- os_model_data$os_months / 12
}

required_variables <- c(
  "os_years",
  "os_event",
  "age_at_diagnosis",
  "npi",
  "cohort",
  "molecular_subtype",
  "tumor_stage"
)

missing_variables <- setdiff(
  required_variables,
  names(os_model_data)
)

if (length(missing_variables) > 0) {
  stop(
    "Missing required variable(s): ",
    paste(missing_variables, collapse = ", ")
  )
}

if (!all(stats::complete.cases(
  os_model_data[, setdiff(required_variables, "tumor_stage")]
))) {
  stop("At least one primary model variable is missing.")
}


# Preserve primary factor levels and spline knots

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

os_model_data$cohort <- factor(os_model_data$cohort)

age_knots <- primary_results$age_knots
npi_knots <- primary_results$npi_knots


# Define the two sensitivity cohorts

stage_sensitivity_data <- os_model_data[
  is.na(os_model_data$tumor_stage) |
    !(os_model_data$tumor_stage %in% c(0, 4)),
  ,
  drop = FALSE
]

normal_like_sensitivity_data <- droplevels(
  os_model_data[
    os_model_data$molecular_subtype != "Normal-like",
    ,
    drop = FALSE
  ]
)


# Verify cohort changes

sensitivity_cohort_summary <- data.frame(
  analysis = c(
    "Primary OS cohort",
    "Exclude recorded Stage 0 or IV",
    "Exclude Normal-like"
  ),
  participants = c(
    nrow(os_model_data),
    nrow(stage_sensitivity_data),
    nrow(normal_like_sensitivity_data)
  ),
  deaths = c(
    sum(os_model_data$os_event),
    sum(stage_sensitivity_data$os_event),
    sum(normal_like_sensitivity_data$os_event)
  )
)

expected_participants <- c(1971L, 1950L, 1823L)
expected_deaths <- c(1138L, 1127L, 1059L)

if (
  !all(sensitivity_cohort_summary$participants == expected_participants) ||
  !all(sensitivity_cohort_summary$deaths == expected_deaths) ||
  sum(is.na(stage_sensitivity_data$tumor_stage)) != 510L
) {
  stop("Sensitivity-cohort counts differ from the validated counts.")
}

print(sensitivity_cohort_summary, row.names = FALSE)
print(table(os_model_data$tumor_stage, useNA = "always"))

# Attach survival so coxph recognizes strata() in the model formulas

suppressPackageStartupMessages(
  library(survival)
)

if (length(age_knots) != 4L || length(npi_knots) != 3L) {
  stop("The saved primary spline knots are incomplete.")
}


# Define the primary follow-up intervals

interval_definitions <- data.frame(
  interval = c(
    "0 to 2 years",
    "2 to 5 years",
    "5 to 10 years",
    "Beyond 10 years"
  ),
  start_time = c(0, 2, 5, 10),
  end_time = c(2, 5, 10, Inf)
)


# Create one interval-specific risk set

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
  
  if (any(interval_data$exit_time <= interval_data$entry_time)) {
    stop("Invalid follow-up time in interval: ", interval_label)
  }
  
  interval_data
}


# Fit one interval-specific sensitivity model

fit_sensitivity_interval <- function(
    data,
    analysis_label,
    interval_label,
    start_time,
    end_time
) {
  interval_data <- create_interval_data(
    data = data,
    start_time = start_time,
    end_time = end_time,
    interval_label = interval_label
  )
  
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
      splines::ns(
        npi,
        knots = npi_knots[2],
        Boundary.knots = npi_knots[c(1, 3)]
      ) +
      strata(cohort),
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
  
  global_result <- data.frame(
    analysis = analysis_label,
    interval = interval_label,
    participants = nrow(interval_data),
    deaths = sum(interval_data$interval_event),
    likelihood_ratio_chisq = model_comparison[2, "Chisq"],
    degrees_of_freedom = model_comparison[2, "Df"],
    p_value = model_comparison[2, "Pr(>|Chi|)"]
  )
  
  coefficient_table <- summary(extended_model)$coefficients
  confidence_intervals <- stats::confint(extended_model)
  
  subtype_rows <- grep(
    "^molecular_subtype",
    rownames(coefficient_table)
  )
  
  subtype_names <- rownames(coefficient_table)[subtype_rows]
  
  subtype_estimates <- data.frame(
    analysis = analysis_label,
    interval = interval_label,
    comparison = paste(
      sub("^molecular_subtype", "", subtype_names),
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
    sensitivity_p_value = coefficient_table[
      subtype_rows,
      "Pr(>|z|)"
    ],
    row.names = NULL
  )
  
  list(
    clinical_model = clinical_model,
    extended_model = extended_model,
    global_result = global_result,
    subtype_estimates = subtype_estimates
  )
}


# Fit all four intervals for one sensitivity cohort

fit_sensitivity_analysis <- function(data, analysis_label) {
  interval_results <- lapply(
    seq_len(nrow(interval_definitions)),
    function(i) {
      fit_sensitivity_interval(
        data = data,
        analysis_label = analysis_label,
        interval_label = interval_definitions$interval[i],
        start_time = interval_definitions$start_time[i],
        end_time = interval_definitions$end_time[i]
      )
    }
  )
  
  names(interval_results) <- interval_definitions$interval
  
  list(
    interval_results = interval_results,
    global_results = do.call(
      rbind,
      lapply(interval_results, function(x) x$global_result)
    ),
    subtype_estimates = do.call(
      rbind,
      lapply(interval_results, function(x) x$subtype_estimates)
    )
  )
}


# Fit both sensitivity analyses

stage_sensitivity_results <- fit_sensitivity_analysis(
  data = stage_sensitivity_data,
  analysis_label = "Exclude recorded Stage 0 or IV"
)

normal_like_sensitivity_results <- fit_sensitivity_analysis(
  data = normal_like_sensitivity_data,
  analysis_label = "Exclude Normal-like"
)

additional_sensitivity_results <- list(
  stage = stage_sensitivity_results,
  normal_like = normal_like_sensitivity_results
)

additional_global_results <- do.call(
  rbind,
  lapply(
    additional_sensitivity_results,
    function(x) x$global_results
  )
)

additional_subtype_estimates <- do.call(
  rbind,
  lapply(
    additional_sensitivity_results,
    function(x) x$subtype_estimates
  )
)


# Compare sensitivity and primary hazard ratios

primary_subtype_estimates <- primary_results$subtype_estimates[
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
  additional_subtype_estimates,
  primary_subtype_estimates,
  by = c("interval", "comparison"),
  all.x = TRUE,
  sort = FALSE
)

hazard_ratio_comparison$percent_change <- 100 * (
  hazard_ratio_comparison$sensitivity_hazard_ratio /
    hazard_ratio_comparison$primary_hazard_ratio - 1
)


# Print results for review

cat("\nGlobal subtype tests:\n")
print(additional_global_results, row.names = FALSE)

cat("\nComparison with primary subtype hazard ratios:\n")
print(hazard_ratio_comparison, row.names = FALSE)

cat("\nLargest absolute percentage change by analysis:\n")
print(
  aggregate(
    abs(percent_change) ~ analysis,
    data = hazard_ratio_comparison,
    FUN = max
  ),
  row.names = FALSE
)


# Prepare ordered export tables

analysis_levels <- c(
  "Exclude recorded Stage 0 or IV",
  "Exclude Normal-like"
)

comparison_levels <- paste(
  subtype_levels[-1],
  "vs Luminal A"
)

additional_global_results$analysis <- factor(
  additional_global_results$analysis,
  levels = analysis_levels
)

additional_global_results$interval <- factor(
  additional_global_results$interval,
  levels = interval_definitions$interval
)

global_results_export <- additional_global_results[
  order(
    additional_global_results$analysis,
    additional_global_results$interval
  ),
]

global_results_export$analysis <- as.character(
  global_results_export$analysis
)

global_results_export$interval <- as.character(
  global_results_export$interval
)

global_results_export$likelihood_ratio_chisq <- round(
  global_results_export$likelihood_ratio_chisq,
  3
)

global_results_export$p_value_display <- ifelse(
  global_results_export$p_value < 0.001,
  "<0.001",
  sprintf("%.3f", global_results_export$p_value)
)

hazard_ratio_comparison$analysis <- factor(
  hazard_ratio_comparison$analysis,
  levels = analysis_levels
)

hazard_ratio_comparison$interval <- factor(
  hazard_ratio_comparison$interval,
  levels = interval_definitions$interval
)

hazard_ratio_comparison$comparison <- factor(
  hazard_ratio_comparison$comparison,
  levels = comparison_levels
)

hazard_ratio_export <- hazard_ratio_comparison[
  order(
    hazard_ratio_comparison$analysis,
    hazard_ratio_comparison$interval,
    hazard_ratio_comparison$comparison
  ),
]

hazard_ratio_export$analysis <- as.character(
  hazard_ratio_export$analysis
)

hazard_ratio_export$interval <- as.character(
  hazard_ratio_export$interval
)

hazard_ratio_export$comparison <- as.character(
  hazard_ratio_export$comparison
)

hazard_ratio_columns <- c(
  "primary_hazard_ratio",
  "primary_lower_95_ci",
  "primary_upper_95_ci",
  "sensitivity_hazard_ratio",
  "sensitivity_lower_95_ci",
  "sensitivity_upper_95_ci"
)

hazard_ratio_export[hazard_ratio_columns] <- lapply(
  hazard_ratio_export[hazard_ratio_columns],
  round,
  digits = 3
)

hazard_ratio_export$percent_change <- round(
  hazard_ratio_export$percent_change,
  1
)

hazard_ratio_export$sensitivity_p_value_display <- ifelse(
  hazard_ratio_export$sensitivity_p_value < 0.001,
  "<0.001",
  sprintf("%.3f", hazard_ratio_export$sensitivity_p_value)
)


# Save model objects and exported tables

model_output_path <- file.path(
  "data-derived",
  "os_additional_sensitivity.rds"
)

table_directory <- file.path("output", "tables")
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

cohort_output_path <- file.path(
  table_directory,
  "os_additional_sensitivity_cohorts.csv"
)

global_output_path <- file.path(
  table_directory,
  "os_additional_sensitivity_global_tests.csv"
)

hazard_ratio_output_path <- file.path(
  table_directory,
  "os_additional_sensitivity_hazard_ratios.csv"
)

saveRDS(
  list(
    cohort_summary = sensitivity_cohort_summary,
    sensitivity_results = additional_sensitivity_results,
    global_results = additional_global_results,
    hazard_ratio_comparison = hazard_ratio_comparison
  ),
  model_output_path
)

utils::write.csv(
  sensitivity_cohort_summary,
  cohort_output_path,
  row.names = FALSE
)

utils::write.csv(
  global_results_export,
  global_output_path,
  row.names = FALSE
)

utils::write.csv(
  hazard_ratio_export,
  hazard_ratio_output_path,
  row.names = FALSE
)

message("")
message("Script 13 completed successfully.")
message("Model results: ", normalizePath(model_output_path))
message("Tables:")
message("  ", normalizePath(cohort_output_path))
message("  ", normalizePath(global_output_path))
message("  ", normalizePath(hazard_ratio_output_path))