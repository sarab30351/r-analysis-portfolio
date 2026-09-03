# 25_rfs_additional_sensitivity.R

# Purpose: Assess whether the final five-interval RFS results are sensitive to excluding the three events recorded at month zero or the Normal-like subtype.

# Both analyses use the final interval structure and model specification from Script 22. They are sensitivity analyses and do not replace the final models. Their p-values are not included in the final Holm-adjustment families.


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


# 2. Load the cohort and final RFS results

input_paths <- c(
  cohort = file.path(
    project_root,
    "data-derived",
    "metabric_rfs_cohort.rds"
  ),
  primary_models = file.path(
    project_root,
    "data-derived",
    "rfs_interval_cox_models.rds"
  ),
  final_results = file.path(
    project_root,
    "data-derived",
    "rfs_final_interval_analysis.rds"
  )
)

missing_inputs <- input_paths[!file.exists(input_paths)]

if (length(missing_inputs) > 0L) {
  stop(
    "Run Scripts 18 and 22 first. Missing: ",
    paste(missing_inputs, collapse = ", ")
  )
}

rfs_data <- as.data.frame(readRDS(input_paths["cohort"]))
primary_results <- readRDS(input_paths["primary_models"])
final_results <- readRDS(input_paths["final_results"])

required_primary_components <- c(
  "spline_knots"
)

required_final_components <- c(
  "interval_definitions",
  "model_overview",
  "subtype_hazard_ratios"
)

if (
  any(!required_primary_components %in% names(primary_results)) ||
  any(!required_final_components %in% names(final_results))
) {
  stop("A saved RFS analysis has an unexpected structure.")
}


# 3. Prepare and validate the primary cohort

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

month_zero_event <- (
  rfs_data$rfs_years == 0 &
    rfs_data$rfs_event == 1L
)

if (
  nrow(rfs_data) != 1960L ||
  sum(rfs_data$rfs_event) != 790L ||
  sum(month_zero_event) != 3L ||
  any(rfs_data$rfs_years < 0) ||
  !all(rfs_data$rfs_event %in% c(0, 1)) ||
  levels(rfs_data$molecular_subtype)[1] != "Luminal A"
) {
  stop("The data do not match the validated RFS cohort.")
}


# 4. Recover the final interval and spline specifications

interval_definitions <- final_results$interval_definitions

expected_interval_keys <- c(
  "zero_to_two",
  "two_to_three_point_five",
  "three_point_five_to_five",
  "five_to_ten",
  "beyond_ten"
)

expected_interval_labels <- c(
  "0 to 2 years",
  "2 to 3.5 years",
  "3.5 to 5 years",
  "5 to 10 years",
  "Beyond 10 years"
)

if (
  !all(
    c(
      "interval_key",
      "interval",
      "start_year",
      "end_year"
    ) %in% names(interval_definitions)
  ) ||
  !identical(
    as.character(interval_definitions$interval_key),
    expected_interval_keys
  ) ||
  !identical(
    as.character(interval_definitions$interval),
    expected_interval_labels
  )
) {
  stop("The final interval specification is unexpected.")
}

spline_knots <- primary_results$spline_knots

if (!all(c("variable", "percentile", "value") %in% names(spline_knots))) {
  stop("The saved spline-knot table has unexpected columns.")
}

extract_knots <- function(variable_label) {
  rows <- spline_knots$variable == variable_label
  
  spline_knots$value[rows][
    order(spline_knots$percentile[rows])
  ]
}

age_knots <- extract_knots("Age at diagnosis")
npi_knots <- extract_knots("NPI")

if (
  length(age_knots) != 4L ||
  length(npi_knots) != 3L ||
  anyDuplicated(age_knots) > 0L ||
  anyDuplicated(npi_knots) > 0L
) {
  stop("The primary spline-knot specification is invalid.")
}


# 5. Create the two sensitivity cohorts

month_zero_sensitivity_data <- rfs_data[
  !month_zero_event,
  ,
  drop = FALSE
]

normal_like_sensitivity_data <- droplevels(
  rfs_data[
    rfs_data$molecular_subtype != "Normal-like",
    ,
    drop = FALSE
  ]
)

contrasts(
  normal_like_sensitivity_data$molecular_subtype
) <- stats::contr.treatment(
  levels(normal_like_sensitivity_data$molecular_subtype),
  base = 1
)

sensitivity_cohort_summary <- data.frame(
  analysis = c(
    "Primary five-interval RFS analysis",
    "Exclude three month-zero RFS events",
    "Exclude Normal-like subtype"
  ),
  participants = c(
    nrow(rfs_data),
    nrow(month_zero_sensitivity_data),
    nrow(normal_like_sensitivity_data)
  ),
  rfs_events = c(
    sum(rfs_data$rfs_event),
    sum(month_zero_sensitivity_data$rfs_event),
    sum(normal_like_sensitivity_data$rfs_event)
  ),
  participants_removed = c(0L, 3L, 146L),
  rfs_events_removed = c(0L, 3L, 63L),
  stringsAsFactors = FALSE
)

if (
  !identical(
    sensitivity_cohort_summary$participants,
    c(1960L, 1957L, 1814L)
  ) ||
  !identical(
    sensitivity_cohort_summary$rfs_events,
    c(790L, 787L, 727L)
  ) ||
  any(month_zero_sensitivity_data$rfs_years == 0) ||
  "Normal-like" %in%
  levels(normal_like_sensitivity_data$molecular_subtype)
) {
  stop("A sensitivity cohort failed validation.")
}

sensitivity_data <- list(
  month_zero = month_zero_sensitivity_data,
  normal_like = normal_like_sensitivity_data
)

analysis_labels <- c(
  month_zero = "Exclude three month-zero RFS events",
  normal_like = "Exclude Normal-like subtype"
)


# 6. Define the interval risk-set function

create_interval_data <- function(data, start_year, end_year, interval_label) {
  at_risk <- if (start_year == 0) {
    data$rfs_years >= start_year
  } else {
    data$rfs_years > start_year
  }
  
  interval_data <- data[
    at_risk,
    ,
    drop = FALSE
  ]
  
  interval_exit <- if (is.infinite(end_year)) {
    interval_data$rfs_years
  } else {
    pmin(interval_data$rfs_years, end_year)
  }
  
  interval_data$time_since_interval_start <- (
    interval_exit - start_year
  )
  
  event_within_interval <- if (is.infinite(end_year)) {
    rep(TRUE, nrow(interval_data))
  } else {
    interval_data$rfs_years <= end_year
  }
  
  interval_data$interval_event <- as.integer(
    interval_data$rfs_event == 1L &
      event_within_interval
  )
  
  interval_data$interval <- interval_label
  
  if (
    nrow(interval_data) == 0L ||
    sum(interval_data$interval_event) == 0L ||
    any(interval_data$time_since_interval_start < 0) ||
    any(
      interval_data$time_since_interval_start == 0 &
      interval_data$interval_event == 0
    ) ||
    any(table(interval_data$molecular_subtype) == 0L)
  ) {
    stop("The risk set failed validation for: ", interval_label)
  }
  
  interval_data
}


# 7. Fit one interval-specific sensitivity model

fit_interval <- function(
    interval_data,
    analysis_label,
    interval_label
) {
  clinical_formula <- survival::Surv(
    time_since_interval_start,
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
    strata(cohort)
  
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
  expected_events <- sum(interval_data$interval_event)
  
  if (
    clinical_model$n != nrow(interval_data) ||
    extended_model$n != nrow(interval_data) ||
    clinical_model$nevent != expected_events ||
    extended_model$nevent != expected_events
  ) {
    stop("A model did not use the complete risk set for: ", interval_label)
  }
  
  clinical_loglik <- stats::logLik(clinical_model)
  extended_loglik <- stats::logLik(extended_model)
  clinical_parameters <- attr(clinical_loglik, "df")
  extended_parameters <- attr(extended_loglik, "df")
  likelihood_ratio_df <- extended_parameters - clinical_parameters
  expected_subtype_parameters <- nlevels(interval_data$molecular_subtype) - 1L
  
  if (
    clinical_parameters != 5L ||
    extended_parameters != 5L + expected_subtype_parameters ||
    likelihood_ratio_df != expected_subtype_parameters
  ) {
    stop("Unexpected parameter counts for: ", interval_label)
  }
  
  likelihood_ratio_chisq <- 2 * (
    as.numeric(extended_loglik) - as.numeric(clinical_loglik)
  )
  
  model_overview <- data.frame(
    analysis = analysis_label,
    interval = interval_label,
    participants_at_risk = nrow(interval_data),
    rfs_events = expected_events,
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
  
  model_summary <- summary(extended_model)
  coefficient_table <- model_summary$coefficients
  confidence_table <- model_summary$conf.int
  
  subtype_terms <- grep(
    "^molecular_subtype",
    rownames(coefficient_table),
    value = TRUE
  )
  
  if (length(subtype_terms) != expected_subtype_parameters) {
    stop("The subtype coefficients are incomplete for: ", interval_label)
  }
  
  comparison_subtypes <- sub(
    "^molecular_subtype",
    "",
    subtype_terms
  )
  
  subtype_estimates <- data.frame(
    analysis = analysis_label,
    interval = interval_label,
    term = subtype_terms,
    comparison = paste(
      comparison_subtypes,
      "vs Luminal A"
    ),
    sensitivity_hazard_ratio = confidence_table[
      subtype_terms,
      "exp(coef)"
    ],
    sensitivity_lower_95_ci = confidence_table[
      subtype_terms,
      "lower .95"
    ],
    sensitivity_upper_95_ci = confidence_table[
      subtype_terms,
      "upper .95"
    ],
    sensitivity_p_value = coefficient_table[
      subtype_terms,
      "Pr(>|z|)"
    ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  list(
    interval = interval_label,
    interval_data = interval_data,
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_overview = model_overview,
    subtype_estimates = subtype_estimates
  )
}


# 8. Fit all five intervals for one sensitivity cohort

fit_sensitivity_analysis <- function(data, analysis_label) {
  interval_data <- setNames(
    lapply(
      seq_len(nrow(interval_definitions)),
      function(index) {
        create_interval_data(
          data = data,
          start_year = interval_definitions$start_year[index],
          end_year = interval_definitions$end_year[index],
          interval_label = interval_definitions$interval[index]
        )
      }
    ),
    interval_definitions$interval_key
  )
  
  participant_counts <- vapply(interval_data, nrow, integer(1))
  event_counts <- vapply(
    interval_data,
    function(interval) sum(interval$interval_event),
    integer(1)
  )
  
  if (
    participant_counts[1] != nrow(data) ||
    any(diff(participant_counts) >= 0) ||
    sum(event_counts) != sum(data$rfs_event)
  ) {
    stop("The five risk sets failed validation for: ", analysis_label)
  }
  
  interval_results <- setNames(
    lapply(
      seq_len(nrow(interval_definitions)),
      function(index) {
        interval_key <- interval_definitions$interval_key[index]
        
        fit_interval(
          interval_data = interval_data[[interval_key]],
          analysis_label = analysis_label,
          interval_label = interval_definitions$interval[index]
        )
      }
    ),
    interval_definitions$interval_key
  )
  
  combine_component <- function(component) {
    result <- do.call(
      rbind,
      lapply(
        interval_results,
        function(interval_result) interval_result[[component]]
      )
    )
    
    row.names(result) <- NULL
    result
  }
  
  list(
    interval_results = interval_results,
    model_overview = combine_component("model_overview"),
    subtype_estimates = combine_component("subtype_estimates")
  )
}


# 9. Fit both sensitivity analyses and compare with the final models

sensitivity_results <- setNames(
  lapply(
    names(sensitivity_data),
    function(analysis_key) {
      fit_sensitivity_analysis(
        data = sensitivity_data[[analysis_key]],
        analysis_label = analysis_labels[[analysis_key]]
      )
    }
  ),
  names(sensitivity_data)
)

model_overview <- do.call(
  rbind,
  lapply(
    sensitivity_results,
    function(result) result$model_overview
  )
)

sensitivity_estimates <- do.call(
  rbind,
  lapply(
    sensitivity_results,
    function(result) result$subtype_estimates
  )
)

row.names(model_overview) <- NULL
row.names(sensitivity_estimates) <- NULL

primary_estimate_columns <- c(
  "interval",
  "comparison",
  "hazard_ratio",
  "lower_95_ci",
  "upper_95_ci"
)

if (
  any(
    !primary_estimate_columns %in%
    names(final_results$subtype_hazard_ratios)
  )
) {
  stop("The final subtype estimate table has unexpected columns.")
}

primary_estimates <- final_results$subtype_hazard_ratios[
  ,
  primary_estimate_columns,
  drop = FALSE
]

names(primary_estimates)[3:5] <- c(
  "primary_hazard_ratio",
  "primary_lower_95_ci",
  "primary_upper_95_ci"
)

primary_rows <- match(
  paste(
    sensitivity_estimates$interval,
    sensitivity_estimates$comparison
  ),
  paste(
    primary_estimates$interval,
    primary_estimates$comparison
  )
)

if (anyNA(primary_rows)) {
  stop("Sensitivity estimates could not be matched to final estimates.")
}

hazard_ratio_comparison <- cbind(
  sensitivity_estimates,
  primary_estimates[
    primary_rows,
    3:5,
    drop = FALSE
  ]
)

hazard_ratio_comparison$percent_change <- 100 * (
  hazard_ratio_comparison$sensitivity_hazard_ratio /
    hazard_ratio_comparison$primary_hazard_ratio - 1
)

hazard_ratio_comparison$absolute_percent_change <- abs(
  hazard_ratio_comparison$percent_change
)

change_summary <- do.call(
  rbind,
  lapply(
    analysis_labels,
    function(analysis_label) {
      rows <- hazard_ratio_comparison$analysis == analysis_label
      
      data.frame(
        analysis = analysis_label,
        comparisons = sum(rows),
        largest_absolute_percent_change = max(
          hazard_ratio_comparison$absolute_percent_change[rows]
        ),
        median_absolute_percent_change = stats::median(
          hazard_ratio_comparison$absolute_percent_change[rows]
        ),
        comparisons_changing_null_direction = sum(
          sign(
            hazard_ratio_comparison$sensitivity_hazard_ratio[rows] - 1
          ) !=
            sign(
              hazard_ratio_comparison$primary_hazard_ratio[rows] - 1
            )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

row.names(change_summary) <- NULL

if (
  nrow(model_overview) != 10L ||
  nrow(hazard_ratio_comparison) != 45L ||
  nrow(change_summary) != 2L
) {
  stop("The combined sensitivity results have unexpected dimensions.")
}


# 10. Display the sensitivity results

cat("\nSensitivity-cohort summary:\n")
print(sensitivity_cohort_summary, row.names = FALSE)

cat("\nSensitivity-model overview:\n")

overview_print <- model_overview
overview_print$likelihood_ratio_chisq <- round(
  overview_print$likelihood_ratio_chisq,
  3
)
overview_print$p_value <- format.pval(
  model_overview$p_value,
  digits = 4,
  eps = 0.001
)

print(
  overview_print[
    ,
    c(
      "analysis",
      "interval",
      "participants_at_risk",
      "rfs_events",
      "likelihood_ratio_chisq",
      "degrees_of_freedom",
      "p_value"
    )
  ],
  row.names = FALSE
)

cat("\nComparison with the final subtype estimates:\n")

comparison_print <- hazard_ratio_comparison
comparison_print[
  c(
    "sensitivity_hazard_ratio",
    "sensitivity_lower_95_ci",
    "sensitivity_upper_95_ci",
    "primary_hazard_ratio",
    "percent_change",
    "absolute_percent_change"
  )
] <- round(
  comparison_print[
    c(
      "sensitivity_hazard_ratio",
      "sensitivity_lower_95_ci",
      "sensitivity_upper_95_ci",
      "primary_hazard_ratio",
      "percent_change",
      "absolute_percent_change"
    )
  ],
  3
)
comparison_print$sensitivity_p_value <- format.pval(
  hazard_ratio_comparison$sensitivity_p_value,
  digits = 4,
  eps = 0.001
)

print(
  comparison_print[
    ,
    c(
      "analysis",
      "interval",
      "comparison",
      "primary_hazard_ratio",
      "sensitivity_hazard_ratio",
      "sensitivity_lower_95_ci",
      "sensitivity_upper_95_ci",
      "sensitivity_p_value",
      "percent_change",
      "absolute_percent_change"
    )
  ],
  row.names = FALSE
)

cat("\nChange summary:\n")

change_print <- change_summary
change_print$largest_absolute_percent_change <- round(
  change_print$largest_absolute_percent_change,
  3
)
change_print$median_absolute_percent_change <- round(
  change_print$median_absolute_percent_change,
  3
)
print(change_print, row.names = FALSE)

cat(
  paste0(
    "\nThese analyses test robustness to two prespecified exclusions. ",
    "Their p-values are descriptive sensitivity results and are not added ",
    "to the final Holm-adjustment families.\n"
  )
)


# 11. Export and validate the sensitivity results

table_directory <- file.path(project_root, "output", "tables")
data_directory <- file.path(project_root, "data-derived")

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)

tables <- list(
  rfs_additional_sensitivity_cohorts = sensitivity_cohort_summary,
  rfs_additional_sensitivity_model_overview = model_overview,
  rfs_additional_sensitivity_hazard_ratios = hazard_ratio_comparison,
  rfs_additional_sensitivity_change_summary = change_summary
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
  "rfs_additional_sensitivity.rds"
)

saveRDS(
  list(
    sensitivity_cohort_summary = sensitivity_cohort_summary,
    interval_definitions = interval_definitions,
    sensitivity_results = sensitivity_results,
    model_overview = model_overview,
    hazard_ratio_comparison = hazard_ratio_comparison,
    change_summary = change_summary,
    multiplicity_note = paste(
      "Sensitivity p-values are descriptive and are not included",
      "in the final Holm-adjustment families."
    )
  ),
  results_path
)

expected_rows <- c(
  rfs_additional_sensitivity_cohorts = 3L,
  rfs_additional_sensitivity_model_overview = 10L,
  rfs_additional_sensitivity_hazard_ratios = 45L,
  rfs_additional_sensitivity_change_summary = 2L
)

observed_rows <- vapply(
  table_paths,
  function(path) nrow(utils::read.csv(path)),
  integer(1)
)

if (!identical(unname(observed_rows), unname(expected_rows))) {
  stop("At least one sensitivity table has an unexpected row count.")
}

if (!file.exists(results_path)) {
  stop("The additional sensitivity RDS file was not created.")
}

cat("\nScript 25 completed successfully.\n")
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