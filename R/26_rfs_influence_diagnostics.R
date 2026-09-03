# 26_rfs_influence_diagnostics.R

# Purpose: Assess whether individual participants substantially influence the final five-interval RFS subtype hazard ratios assembled in Script 22.

# DFBETAS provide a fast approximation to the change caused by deleting one participant. For each of the 25 subtype coefficients, this script screens the three participants with the largest absolute DFBETAS and then performs exact leave-one-participant-out refits. These are diagnostics, not new outcome tests, and their results are not added to the Holm-adjustment families.


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


# 2. Load the cohort and saved final RFS analysis

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

required_primary_components <- "spline_knots"

required_final_components <- c(
  "interval_definitions",
  "models",
  "model_overview",
  "subtype_hazard_ratios",
  "subtype_event_counts"
)

if (
  any(!required_primary_components %in% names(primary_results)) ||
  any(!required_final_components %in% names(final_results))
) {
  stop("A saved RFS analysis has an unexpected structure.")
}


# 3. Prepare and validate the RFS cohort

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

if (
  anyNA(rfs_data[model_columns]) ||
  anyDuplicated(rfs_data$patient_id) > 0L ||
  nrow(rfs_data) != 1960L ||
  sum(rfs_data$rfs_event) != 790L ||
  sum(rfs_data$rfs_years == 0 & rfs_data$rfs_event == 1L) != 3L ||
  any(rfs_data$rfs_years < 0) ||
  !all(rfs_data$rfs_event %in% c(0, 1)) ||
  levels(rfs_data$molecular_subtype)[1] != "Luminal A"
) {
  stop("The data do not match the validated RFS cohort.")
}


# 4. Recover and validate the final specifications

interval_definitions <- final_results$interval_definitions
final_models <- final_results$models
final_overview <- final_results$model_overview

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

required_interval_columns <- c(
  "interval_key",
  "interval",
  "start_year",
  "end_year"
)

required_overview_columns <- c(
  "interval",
  "participants_at_risk",
  "rfs_events"
)

if (
  any(!required_interval_columns %in% names(interval_definitions)) ||
  any(!required_overview_columns %in% names(final_overview)) ||
  !identical(
    as.character(interval_definitions$interval_key),
    expected_interval_keys
  ) ||
  !identical(
    as.character(interval_definitions$interval),
    expected_interval_labels
  ) ||
  !identical(names(final_models), expected_interval_keys) ||
  nrow(final_overview) != 5L
) {
  stop("The final five-interval specification is unexpected.")
}

valid_models <- vapply(
  final_models,
  function(model) {
    inherits(model, "coxph") &&
      !is.null(model$x) &&
      !is.null(model$y)
  },
  logical(1)
)

if (!all(valid_models)) {
  stop("At least one final model is incomplete or is not a Cox model.")
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


# 5. Recreate the five final risk sets

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

interval_data_sets <- setNames(
  lapply(
    seq_len(nrow(interval_definitions)),
    function(index) {
      create_interval_data(
        data = rfs_data,
        start_year = interval_definitions$start_year[index],
        end_year = interval_definitions$end_year[index],
        interval_label = interval_definitions$interval[index]
      )
    }
  ),
  interval_definitions$interval_key
)

recreated_participant_counts <- vapply(
  interval_data_sets,
  nrow,
  integer(1)
)

recreated_event_counts <- vapply(
  interval_data_sets,
  function(data) sum(data$interval_event),
  integer(1)
)

overview_rows <- match(
  expected_interval_labels,
  final_overview$interval
)

if (
  anyNA(overview_rows) ||
  !all(
    recreated_participant_counts ==
    final_overview$participants_at_risk[overview_rows]
  ) ||
  !all(
    recreated_event_counts ==
    final_overview$rfs_events[overview_rows]
  )
) {
  stop("The recreated risk sets do not match the final model overview.")
}


# 6. Confirm that recreated data align with every saved final model

extended_formula <- survival::Surv(
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
  molecular_subtype +
  strata(cohort)

for (interval_key in expected_interval_keys) {
  saved_model <- final_models[[interval_key]]
  interval_data <- interval_data_sets[[interval_key]]
  interval_label <- interval_data_sets[[interval_key]]$interval[1]
  
  coefficient_names <- names(stats::coef(saved_model))
  subtype_terms <- grep(
    "^molecular_subtype",
    coefficient_names,
    value = TRUE
  )
  
  responses_match <- isTRUE(
    all.equal(
      unname(as.matrix(saved_model$y)),
      unname(
        cbind(
          time = interval_data$time_since_interval_start,
          status = interval_data$interval_event
        )
      ),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )
  
  row_order_matches <- identical(
    row.names(saved_model$model),
    row.names(interval_data)
  )
  
  if (
    is.null(saved_model$model) ||
    saved_model$n != nrow(interval_data) ||
    saved_model$nevent != sum(interval_data$interval_event) ||
    length(coefficient_names) != 10L ||
    length(subtype_terms) != 5L ||
    !row_order_matches ||
    !responses_match
  ) {
    stop(
      "The recreated data do not align with the saved model for: ",
      interval_label
    )
  }
}


# 7. Validate subtype-specific event information

count_subtype_events <- function(interval_data, interval_label) {
  do.call(
    rbind,
    lapply(
      subtype_levels,
      function(subtype) {
        subtype_data <- interval_data[
          interval_data$molecular_subtype == subtype,
          ,
          drop = FALSE
        ]
        
        data.frame(
          analysis = "Final five-interval RFS influence diagnostics",
          interval = interval_label,
          molecular_subtype = subtype,
          participants_at_risk = nrow(subtype_data),
          rfs_events = sum(subtype_data$interval_event),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

subtype_event_counts <- do.call(
  rbind,
  lapply(
    seq_along(expected_interval_keys),
    function(index) {
      count_subtype_events(
        interval_data = interval_data_sets[[
          expected_interval_keys[index]
        ]],
        interval_label = expected_interval_labels[index]
      )
    }
  )
)

row.names(subtype_event_counts) <- NULL

saved_event_counts <- final_results$subtype_event_counts

required_event_columns <- c(
  "interval",
  "molecular_subtype",
  "participants_at_risk",
  "rfs_events"
)

if (any(!required_event_columns %in% names(saved_event_counts))) {
  stop("The saved subtype-event table has unexpected columns.")
}

saved_event_rows <- match(
  paste(
    subtype_event_counts$interval,
    subtype_event_counts$molecular_subtype
  ),
  paste(
    saved_event_counts$interval,
    saved_event_counts$molecular_subtype
  )
)

if (
  anyNA(saved_event_rows) ||
  !all(
    subtype_event_counts$participants_at_risk ==
    saved_event_counts$participants_at_risk[saved_event_rows]
  ) ||
  !all(
    subtype_event_counts$rfs_events ==
    saved_event_counts$rfs_events[saved_event_rows]
  )
) {
  stop("The subtype-specific event counts do not match Script 22.")
}


# 8. Screen the three largest absolute DFBETAS per subtype coefficient

extract_dfbeta_screen <- function(
    model,
    interval_data,
    interval_key,
    interval_label,
    top_n = 3L
) {
  dfbetas <- stats::residuals(
    model,
    type = "dfbetas"
  )
  
  if (is.null(dim(dfbetas))) {
    dfbetas <- matrix(dfbetas, ncol = 1L)
  }
  
  coefficient_names <- names(stats::coef(model))
  
  if (
    ncol(dfbetas) != length(coefficient_names) ||
    nrow(dfbetas) != nrow(interval_data)
  ) {
    stop("DFBETAS dimensions differ for: ", interval_label)
  }
  
  colnames(dfbetas) <- coefficient_names
  
  subtype_terms <- grep(
    "^molecular_subtype",
    coefficient_names,
    value = TRUE
  )
  
  if (length(subtype_terms) != 5L) {
    stop("Expected five subtype coefficients for: ", interval_label)
  }
  
  if (any(!is.finite(dfbetas[, subtype_terms, drop = FALSE]))) {
    stop("Non-finite subtype DFBETAS found for: ", interval_label)
  }
  
  do.call(
    rbind,
    lapply(
      subtype_terms,
      function(coefficient_name) {
        values <- dfbetas[, coefficient_name]
        
        selected_rows <- head(
          order(abs(values), decreasing = TRUE),
          top_n
        )
        
        selected_data <- interval_data[
          selected_rows,
          ,
          drop = FALSE
        ]
        
        data.frame(
          interval_key = interval_key,
          interval = interval_label,
          term = coefficient_name,
          comparison = paste(
            sub("^molecular_subtype", "", coefficient_name),
            "vs Luminal A"
          ),
          influence_rank = seq_along(selected_rows),
          patient_id = as.character(selected_data$patient_id),
          participant_subtype = as.character(
            selected_data$molecular_subtype
          ),
          interval_event = selected_data$interval_event,
          dfbetas = values[selected_rows],
          absolute_dfbetas = abs(values[selected_rows]),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    )
  )
}

dfbeta_screen <- do.call(
  rbind,
  lapply(
    seq_along(expected_interval_keys),
    function(index) {
      interval_key <- expected_interval_keys[index]
      
      extract_dfbeta_screen(
        model = final_models[[interval_key]],
        interval_data = interval_data_sets[[interval_key]],
        interval_key = interval_key,
        interval_label = expected_interval_labels[index]
      )
    }
  )
)

row.names(dfbeta_screen) <- NULL

if (
  nrow(dfbeta_screen) != 75L ||
  anyDuplicated(
    dfbeta_screen[
      ,
      c("interval", "comparison", "influence_rank")
    ]
  ) > 0L
) {
  stop("The DFBETAS screen has unexpected dimensions or duplicates.")
}


# 9. Refit after deleting each screened participant

fit_leave_one_out <- function(screen_row) {
  interval_key <- as.character(screen_row$interval_key)
  interval_label <- as.character(screen_row$interval)
  comparison_label <- as.character(screen_row$comparison)
  coefficient_name <- as.character(screen_row$term)
  omitted_patient <- as.character(screen_row$patient_id)
  
  interval_data <- interval_data_sets[[interval_key]]
  primary_model <- final_models[[interval_key]]
  
  if (is.null(interval_data) || is.null(primary_model)) {
    stop("Missing interval data or model for: ", interval_label)
  }
  
  if (sum(interval_data$patient_id == omitted_patient) != 1L) {
    stop(
      "Expected one matching participant for deletion: ",
      omitted_patient
    )
  }
  
  leave_one_out_data <- interval_data[
    interval_data$patient_id != omitted_patient,
    ,
    drop = FALSE
  ]
  
  leave_one_out_data$molecular_subtype <- factor(
    leave_one_out_data$molecular_subtype,
    levels = subtype_levels
  )
  
  contrasts(
    leave_one_out_data$molecular_subtype
  ) <- stats::contr.treatment(
    subtype_levels,
    base = 1
  )
  
  leave_one_out_data$cohort <- factor(
    leave_one_out_data$cohort,
    levels = levels(rfs_data$cohort)
  )
  
  leave_one_out_model <- survival::coxph(
    formula = extended_formula,
    data = leave_one_out_data,
    ties = "efron",
    na.action = stats::na.fail
  )
  
  primary_coefficients <- stats::coef(primary_model)
  leave_one_out_coefficients <- stats::coef(leave_one_out_model)
  
  if (
    !coefficient_name %in% names(primary_coefficients) ||
    !coefficient_name %in% names(leave_one_out_coefficients) ||
    anyNA(leave_one_out_coefficients)
  ) {
    stop(
      "A subtype coefficient is missing after deleting: ",
      omitted_patient,
      " in ",
      interval_label
    )
  }
  
  primary_beta <- unname(
    primary_coefficients[coefficient_name]
  )
  
  leave_one_out_beta <- unname(
    leave_one_out_coefficients[coefficient_name]
  )
  
  primary_standard_error <- sqrt(
    stats::vcov(primary_model)[
      coefficient_name,
      coefficient_name
    ]
  )
  
  leave_one_out_standard_error <- sqrt(
    stats::vcov(leave_one_out_model)[
      coefficient_name,
      coefficient_name
    ]
  )
  
  primary_hazard_ratio <- exp(primary_beta)
  leave_one_out_hazard_ratio <- exp(leave_one_out_beta)
  
  percent_change <- 100 * (
    leave_one_out_hazard_ratio /
      primary_hazard_ratio - 1
  )
  
  data.frame(
    interval_key = interval_key,
    interval = interval_label,
    term = coefficient_name,
    comparison = comparison_label,
    influence_rank = as.integer(screen_row$influence_rank),
    omitted_patient_id = omitted_patient,
    participant_subtype = as.character(
      screen_row$participant_subtype
    ),
    interval_event = as.integer(screen_row$interval_event),
    dfbetas = as.numeric(screen_row$dfbetas),
    primary_hazard_ratio = primary_hazard_ratio,
    leave_one_out_hazard_ratio = leave_one_out_hazard_ratio,
    leave_one_out_lower_95_ci = exp(
      leave_one_out_beta -
        stats::qnorm(0.975) *
        leave_one_out_standard_error
    ),
    leave_one_out_upper_95_ci = exp(
      leave_one_out_beta +
        stats::qnorm(0.975) *
        leave_one_out_standard_error
    ),
    percent_hazard_ratio_change = percent_change,
    absolute_percent_change = abs(percent_change),
    standardized_log_hr_change = (
      leave_one_out_beta - primary_beta
    ) / primary_standard_error,
    null_direction_changed = (
      sign(leave_one_out_beta) !=
        sign(primary_beta)
    ),
    material_change_10_percent = (
      abs(percent_change) >= 10
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

leave_one_out_results <- do.call(
  rbind,
  lapply(
    seq_len(nrow(dfbeta_screen)),
    function(index) {
      fit_leave_one_out(
        dfbeta_screen[
          index,
          ,
          drop = FALSE
        ]
      )
    }
  )
)

leave_one_out_results$absolute_standardized_change <- abs(
  leave_one_out_results$standardized_log_hr_change
)

if (nrow(leave_one_out_results) != 75L) {
  stop("The exact leave-one-out results should contain 75 rows.")
}


# 10. Summarize the largest exact change for each subtype comparison

comparison_key <- interaction(
  leave_one_out_results$interval_key,
  leave_one_out_results$comparison,
  drop = TRUE
)

comparison_groups <- split(
  leave_one_out_results,
  comparison_key
)

maximum_change_by_comparison <- do.call(
  rbind,
  lapply(
    comparison_groups,
    function(results) {
      maximum_result <- results[
        which.max(results$absolute_percent_change),
        ,
        drop = FALSE
      ]
      
      maximum_result$any_null_direction_changed <- any(
        results$null_direction_changed
      )
      
      maximum_result$any_material_change_10_percent <- any(
        results$material_change_10_percent
      )
      
      maximum_result
    }
  )
)

row.names(maximum_change_by_comparison) <- NULL

interval_order <- match(
  maximum_change_by_comparison$interval_key,
  expected_interval_keys
)

comparison_order <- match(
  maximum_change_by_comparison$comparison,
  paste(subtype_levels[-1], "vs Luminal A")
)

maximum_change_by_comparison <- maximum_change_by_comparison[
  order(interval_order, comparison_order),
  ,
  drop = FALSE
]

if (nrow(maximum_change_by_comparison) != 25L) {
  stop("The maximum-change summary should contain 25 rows.")
}

required_hazard_ratio_columns <- c(
  "interval",
  "comparison",
  "holm_significant"
)

if (
  any(
    !required_hazard_ratio_columns %in%
    names(final_results$subtype_hazard_ratios)
  )
) {
  stop("The final subtype table lacks Holm significance indicators.")
}

holm_significant_rows <- final_results$subtype_hazard_ratios[
  final_results$subtype_hazard_ratios$holm_significant,
  c("interval", "comparison"),
  drop = FALSE
]

headline_keys <- paste(
  holm_significant_rows$interval,
  holm_significant_rows$comparison
)

headline_results <- leave_one_out_results[
  paste(
    leave_one_out_results$interval,
    leave_one_out_results$comparison
  ) %in% headline_keys,
  ,
  drop = FALSE
]

if (
  nrow(holm_significant_rows) != 4L ||
  nrow(headline_results) != 12L
) {
  stop("The Holm-significant headline comparison set is unexpected.")
}

overall_summary <- data.frame(
  screened_deletions = nrow(leave_one_out_results),
  subtype_comparisons = nrow(maximum_change_by_comparison),
  comparisons_with_10_percent_change = sum(
    maximum_change_by_comparison$
      any_material_change_10_percent
  ),
  comparisons_changing_null_direction = sum(
    maximum_change_by_comparison$
      any_null_direction_changed
  ),
  largest_absolute_percent_change = max(
    maximum_change_by_comparison$
      absolute_percent_change
  ),
  largest_absolute_standardized_change = max(
    maximum_change_by_comparison$
      absolute_standardized_change
  ),
  stringsAsFactors = FALSE
)


# 11. Display the influence diagnostics

cat("\nSubtype-specific event information:\n")

print(
  subtype_event_counts,
  row.names = FALSE
)

screen_print <- dfbeta_screen

screen_print$dfbetas <- round(
  screen_print$dfbetas,
  4
)

screen_print$absolute_dfbetas <- round(
  screen_print$absolute_dfbetas,
  4
)

cat(
  "\nThree largest absolute DFBETAS per subtype coefficient:\n"
)

print(
  screen_print[
    ,
    c(
      "interval",
      "comparison",
      "influence_rank",
      "patient_id",
      "participant_subtype",
      "interval_event",
      "dfbetas"
    )
  ],
  row.names = FALSE
)

summary_print <- maximum_change_by_comparison

summary_print[
  c(
    "primary_hazard_ratio",
    "leave_one_out_hazard_ratio",
    "percent_hazard_ratio_change",
    "absolute_percent_change",
    "absolute_standardized_change"
  )
] <- round(
  summary_print[
    c(
      "primary_hazard_ratio",
      "leave_one_out_hazard_ratio",
      "percent_hazard_ratio_change",
      "absolute_percent_change",
      "absolute_standardized_change"
    )
  ],
  3
)

display_columns <- c(
  "interval",
  "comparison",
  "omitted_patient_id",
  "participant_subtype",
  "interval_event",
  "primary_hazard_ratio",
  "leave_one_out_hazard_ratio",
  "percent_hazard_ratio_change",
  "absolute_standardized_change",
  "any_null_direction_changed",
  "any_material_change_10_percent"
)

cat(
  "\nLargest exact deletion change for each comparison:\n"
)

print(
  summary_print[
    ,
    display_columns
  ],
  row.names = FALSE
)

headline_summary <- maximum_change_by_comparison[
  paste(
    maximum_change_by_comparison$interval,
    maximum_change_by_comparison$comparison
  ) %in% headline_keys,
  ,
  drop = FALSE
]

headline_summary[
  c(
    "primary_hazard_ratio",
    "leave_one_out_hazard_ratio",
    "percent_hazard_ratio_change",
    "absolute_standardized_change"
  )
] <- round(
  headline_summary[
    c(
      "primary_hazard_ratio",
      "leave_one_out_hazard_ratio",
      "percent_hazard_ratio_change",
      "absolute_standardized_change"
    )
  ],
  3
)

cat(
  "\nLargest changes for the four Holm-significant comparisons:\n"
)

print(
  headline_summary[
    ,
    display_columns
  ],
  row.names = FALSE
)

overall_print <- overall_summary

overall_print$largest_absolute_percent_change <- round(
  overall_print$largest_absolute_percent_change,
  3
)

overall_print$largest_absolute_standardized_change <- round(
  overall_print$largest_absolute_standardized_change,
  3
)

cat("\nOverall influence summary:\n")

print(
  overall_print,
  row.names = FALSE
)

cat(
  paste0(
    "\nThe 10% flag is a pragmatic descriptive threshold, not a formal ",
    "statistical test. DFBETAS identify participants for exact checking; ",
    "the leave-one-out refits determine the observed change. These ",
    "diagnostics are not added to the Holm-adjustment families.\n"
  )
)


# 12. Export and validate the diagnostic results

table_directory <- file.path(
  project_root,
  "output",
  "tables"
)

data_directory <- file.path(
  project_root,
  "data-derived"
)

dir.create(
  table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  data_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dfbeta_export <- dfbeta_screen

dfbeta_export$dfbetas <- round(
  dfbeta_export$dfbetas,
  4
)

dfbeta_export$absolute_dfbetas <- round(
  dfbeta_export$absolute_dfbetas,
  4
)

leave_one_out_export <- leave_one_out_results
maximum_change_export <- maximum_change_by_comparison

hazard_ratio_columns <- c(
  "primary_hazard_ratio",
  "leave_one_out_hazard_ratio",
  "leave_one_out_lower_95_ci",
  "leave_one_out_upper_95_ci"
)

change_columns <- c(
  "percent_hazard_ratio_change",
  "absolute_percent_change"
)

standardized_columns <- c(
  "standardized_log_hr_change",
  "absolute_standardized_change"
)

leave_one_out_export[
  hazard_ratio_columns
] <- lapply(
  leave_one_out_export[hazard_ratio_columns],
  round,
  digits = 3
)

leave_one_out_export[
  change_columns
] <- lapply(
  leave_one_out_export[change_columns],
  round,
  digits = 1
)

leave_one_out_export[
  standardized_columns
] <- lapply(
  leave_one_out_export[standardized_columns],
  round,
  digits = 3
)

maximum_change_export[
  hazard_ratio_columns
] <- lapply(
  maximum_change_export[hazard_ratio_columns],
  round,
  digits = 3
)

maximum_change_export[
  change_columns
] <- lapply(
  maximum_change_export[change_columns],
  round,
  digits = 1
)

maximum_change_export[
  standardized_columns
] <- lapply(
  maximum_change_export[standardized_columns],
  round,
  digits = 3
)

tables <- list(
  rfs_influence_subtype_event_counts = subtype_event_counts,
  rfs_influence_dfbeta_screen = dfbeta_export,
  rfs_influence_leave_one_out = leave_one_out_export,
  rfs_influence_leave_one_out_summary = maximum_change_export
)

table_paths <- file.path(
  table_directory,
  paste0(
    names(tables),
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
    tables,
    table_paths
  )
)

results_path <- file.path(
  data_directory,
  "rfs_influence_diagnostics.rds"
)

saveRDS(
  list(
    interval_definitions = interval_definitions,
    subtype_event_counts = subtype_event_counts,
    dfbeta_screen = dfbeta_screen,
    leave_one_out_results = leave_one_out_results,
    maximum_change_by_comparison =
      maximum_change_by_comparison,
    holm_significant_comparisons =
      holm_significant_rows,
    overall_summary = overall_summary,
    screening_rule = paste(
      "Three largest absolute DFBETAS for each final interval-specific",
      "subtype coefficient, followed by exact leave-one-participant-out",
      "refitting"
    ),
    interpretation_note = paste(
      "The 10% change flag is descriptive. Influence diagnostics are not",
      "included in the Holm-adjustment families."
    )
  ),
  results_path
)

expected_rows <- c(
  30L,
  75L,
  75L,
  25L
)

written_rows <- vapply(
  table_paths,
  function(path) {
    nrow(
      utils::read.csv(path)
    )
  },
  integer(1)
)

if (
  !file.exists(results_path) ||
  !all(file.exists(table_paths)) ||
  !identical(
    unname(written_rows),
    expected_rows
  )
) {
  stop(
    "One or more influence-diagnostic outputs failed validation."
  )
}

message("")
message("Script 26 completed successfully.")
message(
  "Diagnostic results: ",
  normalizePath(results_path)
)
message("Tables:")

for (path in table_paths) {
  message(
    "  ",
    normalizePath(path)
  )
}