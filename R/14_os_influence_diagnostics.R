# 14_os_influence_diagnostics.R
#
# Purpose: Assess subtype-specific event information and whether individual participants substantively influence the primary interval-specific subtype hazard ratios.


# Check requirements

required_packages <- c("survival", "splines")

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
  stop("Open r-analysis-portfolio.Rproj before running this script.")
}

suppressPackageStartupMessages(
  library(survival)
)


# Load the OS cohort and primary interval models

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
  "patient_id",
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
    "Missing required variable(s): ",
    paste(missing_variables, collapse = ", ")
  )
}

if (!all(stats::complete.cases(os_model_data[, required_variables]))) {
  stop("At least one required variable is missing.")
}

if (anyDuplicated(os_model_data$patient_id)) {
  stop("Patient identifiers are not unique in the OS cohort.")
}


# Preserve the primary model settings

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

if (length(age_knots) != 4L || length(npi_knots) != 3L) {
  stop("The saved primary spline knots are incomplete.")
}


# Recreate the four primary risk sets

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
  
  interval_data$exit_time <- if (is.infinite(end_time)) {
    interval_data$os_years
  } else {
    pmin(interval_data$os_years, end_time)
  }
  
  interval_data$interval_event <- as.integer(
    interval_data$os_event == 1L &
      interval_data$os_years > start_time &
      (
        is.infinite(end_time) |
          interval_data$os_years <= end_time
      )
  )
  
  interval_data$interval <- interval_label
  interval_data
}

interval_data_sets <- Map(
  f = function(start_time, end_time, interval_label) {
    create_interval_data(
      data = os_model_data,
      start_time = start_time,
      end_time = end_time,
      interval_label = interval_label
    )
  },
  start_time = interval_definitions$start_time,
  end_time = interval_definitions$end_time,
  interval_label = interval_definitions$interval
)

names(interval_data_sets) <- interval_definitions$interval


# Count participants and deaths by interval and subtype

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
          interval = interval_label,
          molecular_subtype = subtype,
          participants_at_risk = nrow(subtype_data),
          deaths = sum(subtype_data$interval_event)
        )
      }
    )
  )
}

subtype_event_counts <- do.call(
  rbind,
  Map(
    count_subtype_events,
    interval_data_sets,
    names(interval_data_sets)
  )
)

row.names(subtype_event_counts) <- NULL


# Verify interval totals against the primary results

interval_totals <- data.frame(
  interval = names(interval_data_sets),
  participants = vapply(
    interval_data_sets,
    nrow,
    FUN.VALUE = integer(1)
  ),
  deaths = vapply(
    interval_data_sets,
    function(data) as.integer(sum(data$interval_event)),
    FUN.VALUE = integer(1)
  )
)

primary_index <- match(
  interval_totals$interval,
  primary_results$model_comparisons$interval
)

if (anyNA(primary_index)) {
  stop("At least one primary interval result is missing.")
}

primary_totals <- primary_results$model_comparisons[
  primary_index,
  c("participants", "deaths")
]

if (
  !all(interval_totals$participants == primary_totals$participants) ||
  !all(interval_totals$deaths == primary_totals$deaths)
) {
  stop("Recreated interval totals do not match the primary models.")
}

cat("\nSubtype-specific event information:\n")
print(subtype_event_counts, row.names = FALSE)



# Match the saved primary models to the recreated intervals

primary_interval_labels <- vapply(
  primary_results$interval_results,
  function(result) result$interval,
  FUN.VALUE = character(1)
)

primary_model_index <- match(
  names(interval_data_sets),
  primary_interval_labels
)

if (anyNA(primary_model_index)) {
  stop("At least one saved primary interval model is missing.")
}

primary_models <- lapply(
  primary_results$interval_results[primary_model_index],
  function(result) result$extended_model
)

names(primary_models) <- names(interval_data_sets)


# Identify the three largest absolute DFBETAS per subtype coefficient

extract_dfbeta_screen <- function(
    model,
    interval_data,
    interval_label,
    top_n = 3L
) {
  dfbetas <- stats::residuals(
    model,
    type = "dfbetas"
  )
  
  coefficient_names <- names(
    stats::coef(model)
  )
  
  if (ncol(dfbetas) != length(coefficient_names)) {
    stop(
      "DFBETAS columns do not match model coefficients in interval: ",
      interval_label
    )
  }
  
  colnames(dfbetas) <- coefficient_names
  
  subtype_columns <- grep(
    "^molecular_subtype",
    coefficient_names
  )
  
  if (length(subtype_columns) != 5L) {
    stop(
      "Expected five subtype coefficients in interval: ",
      interval_label
    )
  }
  
  # Cox models preserve the input order but renumber subsetted rows.
  
  if (nrow(dfbetas) != nrow(interval_data)) {
    stop(
      "DFBETAS row count differs in interval: ",
      interval_label
    )
  }
  
  data_index <- seq_len(nrow(interval_data))
  
  do.call(
    rbind,
    lapply(
      subtype_columns,
      function(column) {
        coefficient_name <- colnames(dfbetas)[column]
        values <- dfbetas[, column]
        
        selected_rows <- head(
          order(abs(values), decreasing = TRUE),
          top_n
        )
        
        selected_data <- interval_data[
          data_index[selected_rows],
          ,
          drop = FALSE
        ]
        
        data.frame(
          interval = interval_label,
          comparison = paste(
            sub(
              "^molecular_subtype",
              "",
              coefficient_name
            ),
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
          row.names = NULL
        )
      }
    )
  )
}

dfbeta_screen <- do.call(
  rbind,
  Map(
    extract_dfbeta_screen,
    primary_models,
    interval_data_sets,
    names(interval_data_sets)
  )
)

row.names(dfbeta_screen) <- NULL

cat("\nLargest absolute DFBETAS for subtype coefficients:\n")
print(dfbeta_screen, row.names = FALSE)


# Refit after deleting each screened participant

fit_leave_one_out <- function(screen_row) {
  interval_label <- as.character(screen_row$interval)
  comparison_label <- as.character(screen_row$comparison)
  omitted_patient <- as.character(screen_row$patient_id)
  
  interval_data <- interval_data_sets[[interval_label]]
  primary_model <- primary_models[[interval_label]]
  
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
  
  leave_one_out_data$cohort <- factor(
    leave_one_out_data$cohort,
    levels = levels(os_model_data$cohort)
  )
  
  leave_one_out_model <- survival::coxph(
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
    data = leave_one_out_data,
    ties = "efron"
  )
  
  subtype_name <- sub(
    " vs Luminal A$",
    "",
    comparison_label
  )
  
  coefficient_name <- paste0(
    "molecular_subtype",
    subtype_name
  )
  
  primary_coefficients <- stats::coef(primary_model)
  leave_one_out_coefficients <- stats::coef(
    leave_one_out_model
  )
  
  if (
    !coefficient_name %in% names(primary_coefficients) ||
    !coefficient_name %in%
    names(leave_one_out_coefficients)
  ) {
    stop(
      "Missing subtype coefficient after deleting: ",
      omitted_patient
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
  leave_one_out_hazard_ratio <- exp(
    leave_one_out_beta
  )
  
  data.frame(
    interval = interval_label,
    comparison = comparison_label,
    influence_rank = screen_row$influence_rank,
    omitted_patient_id = omitted_patient,
    participant_subtype = as.character(
      screen_row$participant_subtype
    ),
    interval_event = screen_row$interval_event,
    dfbetas = screen_row$dfbetas,
    primary_hazard_ratio = primary_hazard_ratio,
    leave_one_out_hazard_ratio =
      leave_one_out_hazard_ratio,
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
    percent_hazard_ratio_change = 100 * (
      leave_one_out_hazard_ratio /
        primary_hazard_ratio - 1
    ),
    standardized_log_hr_change = (
      leave_one_out_beta - primary_beta
    ) / primary_standard_error,
    row.names = NULL
  )
}

leave_one_out_results <- do.call(
  rbind,
  lapply(
    seq_len(nrow(dfbeta_screen)),
    function(i) {
      fit_leave_one_out(
        dfbeta_screen[i, , drop = FALSE]
      )
    }
  )
)

leave_one_out_results$absolute_percent_change <- abs(
  leave_one_out_results$percent_hazard_ratio_change
)

leave_one_out_results$absolute_standardized_change <- abs(
  leave_one_out_results$standardized_log_hr_change
)


# Find the largest exact change for each comparison

comparison_key <- interaction(
  leave_one_out_results$interval,
  leave_one_out_results$comparison,
  drop = TRUE
)

maximum_change_by_comparison <- do.call(
  rbind,
  lapply(
    split(leave_one_out_results, comparison_key),
    function(results) {
      results[
        which.max(results$absolute_percent_change),
        ,
        drop = FALSE
      ]
    }
  )
)

row.names(maximum_change_by_comparison) <- NULL

interval_order <- match(
  maximum_change_by_comparison$interval,
  interval_definitions$interval
)

comparison_order <- match(
  maximum_change_by_comparison$comparison,
  paste(subtype_levels[-1], "vs Luminal A")
)

maximum_change_by_comparison <-
  maximum_change_by_comparison[
    order(interval_order, comparison_order),
    ,
    drop = FALSE
  ]

headline_basal_results <- leave_one_out_results[
  leave_one_out_results$comparison ==
    "Basal-like vs Luminal A" &
    leave_one_out_results$interval %in% c(
      "0 to 2 years",
      "Beyond 10 years"
    ),
  ,
  drop = FALSE
]

display_columns <- c(
  "interval",
  "comparison",
  "omitted_patient_id",
  "primary_hazard_ratio",
  "leave_one_out_hazard_ratio",
  "percent_hazard_ratio_change",
  "absolute_standardized_change"
)

cat(
  "\nLargest exact deletion change for each comparison:\n"
)

print(
  maximum_change_by_comparison[, display_columns],
  row.names = FALSE
)

cat(
  "\nExact deletion results for headline Basal-like findings:\n"
)

print(
  headline_basal_results[, display_columns],
  row.names = FALSE
)

# Prepare export tables

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

leave_one_out_export[hazard_ratio_columns] <- lapply(
  leave_one_out_export[hazard_ratio_columns],
  round,
  digits = 3
)

leave_one_out_export[change_columns] <- lapply(
  leave_one_out_export[change_columns],
  round,
  digits = 1
)

leave_one_out_export[standardized_columns] <- lapply(
  leave_one_out_export[standardized_columns],
  round,
  digits = 3
)

maximum_change_export[hazard_ratio_columns] <- lapply(
  maximum_change_export[hazard_ratio_columns],
  round,
  digits = 3
)

maximum_change_export[change_columns] <- lapply(
  maximum_change_export[change_columns],
  round,
  digits = 1
)

maximum_change_export[standardized_columns] <- lapply(
  maximum_change_export[standardized_columns],
  round,
  digits = 3
)


# Save diagnostic results and tables

model_output_path <- file.path(
  "data-derived",
  "os_influence_diagnostics.rds"
)

table_directory <- file.path("output", "tables")

dir.create(
  table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

event_count_output_path <- file.path(
  table_directory,
  "os_influence_subtype_event_counts.csv"
)

dfbeta_output_path <- file.path(
  table_directory,
  "os_influence_dfbeta_screen.csv"
)

leave_one_out_output_path <- file.path(
  table_directory,
  "os_influence_leave_one_out.csv"
)

summary_output_path <- file.path(
  table_directory,
  "os_influence_leave_one_out_summary.csv"
)

saveRDS(
  list(
    subtype_event_counts = subtype_event_counts,
    dfbeta_screen = dfbeta_screen,
    leave_one_out_results = leave_one_out_results,
    maximum_change_by_comparison =
      maximum_change_by_comparison
  ),
  model_output_path
)

utils::write.csv(
  subtype_event_counts,
  event_count_output_path,
  row.names = FALSE
)

utils::write.csv(
  dfbeta_export,
  dfbeta_output_path,
  row.names = FALSE
)

utils::write.csv(
  leave_one_out_export,
  leave_one_out_output_path,
  row.names = FALSE
)

utils::write.csv(
  maximum_change_export,
  summary_output_path,
  row.names = FALSE
)

message("")
message("Script 14 completed successfully.")
message(
  "Diagnostic results: ",
  normalizePath(model_output_path)
)
message("Tables:")
message("  ", normalizePath(event_count_output_path))
message("  ", normalizePath(dfbeta_output_path))
message("  ", normalizePath(leave_one_out_output_path))
message("  ", normalizePath(summary_output_path))