# 22_rfs_final_interval_results.R

# Purpose: Assemble the final five-interval relapse-free survival analysis and apply Holm adjustment to the final results.

# The initial 2-to-5-year interval retained evidence of time-varying molecular-subtype associations. Script 21 divided that interval at its neutral midpoint of 3.5 years. The midpoint was selected independently of any estimated change point, and both refined intervals had adequate proportional-hazards diagnostics.

# This script does not refit models. It combines the three retained models from Script 18 with the two refined models from Script 21. Holm adjustment is applied separately to the five global subtype tests and the 25 subtype comparisons. Diagnostic and sensitivity p-values from Scripts 19 through 21 remain separate and are not multiplicity adjusted here.


# 1. Check the project location

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

if (!file.exists(file.path(project_root, "r-analysis-portfolio.Rproj"))) {
  stop("Run this script from the project root: ", project_root)
}


# 2. Load the saved RFS analyses

input_paths <- c(
  primary_models = file.path(
    project_root,
    "data-derived",
    "rfs_interval_cox_models.rds"
  ),
  midpoint_refinement = file.path(
    project_root,
    "data-derived",
    "rfs_two_to_five_split_sensitivity.rds"
  )
)

missing_inputs <- input_paths[!file.exists(input_paths)]

if (length(missing_inputs) > 0L) {
  stop(
    "Run Scripts 18 and 21 first. Missing: ",
    paste(missing_inputs, collapse = ", ")
  )
}

primary_results <- readRDS(input_paths["primary_models"])
midpoint_results <- readRDS(input_paths["midpoint_refinement"])


# 3. Validate the saved object structures

require_components <- function(object, components, label) {
  missing_components <- setdiff(components, names(object))
  
  if (length(missing_components) > 0L) {
    stop(
      label,
      " is missing: ",
      paste(missing_components, collapse = ", ")
    )
  }
}

require_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  
  if (length(missing_columns) > 0L) {
    stop(
      label,
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

require_components(
  primary_results,
  c(
    "interval_results",
    "model_overview",
    "global_subtype_tests",
    "subtype_hazard_ratios"
  ),
  "Script 18 results"
)

require_components(
  midpoint_results,
  c(
    "split_definitions",
    "split_results",
    "model_overview",
    "subtype_comparison",
    "ph_tests"
  ),
  "Script 21 results"
)

overview_columns <- c(
  "interval",
  "participants_at_risk",
  "rfs_events",
  "clinical_parameters",
  "extended_parameters"
)

global_test_columns <- c(
  "likelihood_ratio_chisq",
  "degrees_of_freedom",
  "p_value"
)

combined_overview_columns <- c(
  overview_columns,
  global_test_columns
)

hazard_ratio_columns <- c(
  "interval",
  "comparison",
  "hazard_ratio",
  "lower_95_ci",
  "upper_95_ci",
  "p_value"
)

require_columns(
  primary_results$model_overview,
  overview_columns,
  "Script 18 model overview"
)

require_columns(
  primary_results$global_subtype_tests,
  c(
    "interval",
    "participants_at_risk",
    "rfs_events",
    global_test_columns
  ),
  "Script 18 global subtype tests"
)

require_columns(
  midpoint_results$model_overview,
  combined_overview_columns,
  "Script 21 model overview"
)

require_columns(
  primary_results$subtype_hazard_ratios,
  hazard_ratio_columns,
  "Script 18 subtype estimates"
)

require_columns(
  midpoint_results$subtype_comparison,
  hazard_ratio_columns,
  "Script 21 subtype estimates"
)


# 4. Define the final five intervals

final_interval_definitions <- data.frame(
  interval_key = c(
    "zero_to_two",
    "two_to_three_point_five",
    "three_point_five_to_five",
    "five_to_ten",
    "beyond_ten"
  ),
  interval = c(
    "0 to 2 years",
    "2 to 3.5 years",
    "3.5 to 5 years",
    "5 to 10 years",
    "Beyond 10 years"
  ),
  start_year = c(0, 2, 3.5, 5, 10),
  end_year = c(2, 3.5, 5, 10, Inf),
  model_source = c(
    "Script 18",
    "Script 21",
    "Script 21",
    "Script 18",
    "Script 18"
  ),
  boundary_basis = c(
    "Retained from the initial interval specification",
    "Neutral midpoint refinement after the 2-to-5-year subtype diagnostic",
    "Neutral midpoint refinement after the 2-to-5-year subtype diagnostic",
    "Retained from the initial interval specification",
    "Retained from the initial interval specification"
  ),
  stringsAsFactors = FALSE
)

final_interval_keys <- final_interval_definitions$interval_key
final_interval_labels <- final_interval_definitions$interval

retained_interval_labels <- c(
  "0 to 2 years",
  "5 to 10 years",
  "Beyond 10 years"
)

refined_interval_labels <- c(
  "2 to 3.5 years",
  "3.5 to 5 years"
)

if (
  !identical(
    as.character(midpoint_results$split_definitions$interval),
    refined_interval_labels
  )
) {
  stop("The saved midpoint refinement has unexpected intervals.")
}


# 5. Assemble and validate the final models

extract_primary_model <- function(interval_key) {
  result <- primary_results$interval_results[[interval_key]]
  
  if (is.null(result) || is.null(result$extended_model)) {
    stop("The primary extended model is missing for: ", interval_key)
  }
  
  result$extended_model
}

extract_refined_model <- function(interval_key) {
  result <- midpoint_results$split_results[[interval_key]]
  
  if (is.null(result) || is.null(result$extended_model)) {
    stop("The refined extended model is missing for: ", interval_key)
  }
  
  result$extended_model
}

final_models <- list(
  zero_to_two = extract_primary_model("zero_to_two"),
  two_to_three_point_five = extract_refined_model(
    "two_to_three_point_five"
  ),
  three_point_five_to_five = extract_refined_model(
    "three_point_five_to_five"
  ),
  five_to_ten = extract_primary_model("five_to_ten"),
  beyond_ten = extract_primary_model("beyond_ten")
)

if (!identical(names(final_models), final_interval_keys)) {
  stop("The final model order is incorrect.")
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


# 6. Assemble the final global subtype tests

primary_model_summary <- primary_results$model_overview[
  primary_results$model_overview$interval %in% retained_interval_labels,
  overview_columns,
  drop = FALSE
]

primary_global_tests <- primary_results$global_subtype_tests[
  primary_results$global_subtype_tests$interval %in% retained_interval_labels,
  ,
  drop = FALSE
]

primary_global_rows <- match(
  primary_model_summary$interval,
  primary_global_tests$interval
)

if (
  nrow(primary_model_summary) != 3L ||
  nrow(primary_global_tests) != 3L ||
  anyNA(primary_global_rows) ||
  any(
    primary_model_summary$participants_at_risk !=
    primary_global_tests$participants_at_risk[primary_global_rows]
  ) ||
  any(
    primary_model_summary$rfs_events !=
    primary_global_tests$rfs_events[primary_global_rows]
  )
) {
  stop(
    paste(
      "Script 18 model summaries and global subtype tests",
      "could not be matched."
    )
  )
}

primary_overview <- cbind(
  primary_model_summary,
  primary_global_tests[
    primary_global_rows,
    global_test_columns,
    drop = FALSE
  ]
)

refined_overview <- midpoint_results$model_overview[
  midpoint_results$model_overview$interval %in% refined_interval_labels,
  combined_overview_columns,
  drop = FALSE
]

final_model_overview <- rbind(
  primary_overview,
  refined_overview
)

overview_order <- match(
  final_interval_labels,
  final_model_overview$interval
)

if (
  nrow(final_model_overview) != 5L ||
  anyNA(overview_order)
) {
  stop("The five final model summaries could not be assembled.")
}

final_model_overview <- final_model_overview[
  overview_order,
  ,
  drop = FALSE
]

row.names(final_model_overview) <- NULL
final_model_overview$analysis <- "Final five-interval RFS model"
final_model_overview$multiplicity_family <- (
  "Five global molecular-subtype likelihood-ratio tests"
)
final_model_overview$holm_adjusted_p_value <- stats::p.adjust(
  final_model_overview$p_value,
  method = "holm"
)
final_model_overview$holm_significant <- (
  final_model_overview$holm_adjusted_p_value < 0.05
)

final_model_overview <- final_model_overview[
  ,
  c(
    "analysis",
    "interval",
    "participants_at_risk",
    "rfs_events",
    "clinical_parameters",
    "extended_parameters",
    "likelihood_ratio_chisq",
    "degrees_of_freedom",
    "p_value",
    "multiplicity_family",
    "holm_adjusted_p_value",
    "holm_significant"
  )
]

if (
  any(final_model_overview$clinical_parameters != 5L) ||
  any(final_model_overview$extended_parameters != 10L) ||
  any(final_model_overview$degrees_of_freedom != 5L) ||
  sum(final_model_overview$rfs_events) != 790L ||
  any(diff(final_model_overview$participants_at_risk) >= 0)
) {
  stop("The final model overview failed validation.")
}

for (index in seq_along(final_models)) {
  model <- final_models[[index]]
  
  if (
    model$n != final_model_overview$participants_at_risk[index] ||
    model$nevent != final_model_overview$rfs_events[index]
  ) {
    stop(
      "A final model does not match its overview: ",
      final_interval_labels[index]
    )
  }
}


# 7. Assemble the final subtype estimates

extract_estimates <- function(data, interval_labels) {
  data[
    data$interval %in% interval_labels,
    hazard_ratio_columns,
    drop = FALSE
  ]
}

primary_estimates <- extract_estimates(
  primary_results$subtype_hazard_ratios,
  retained_interval_labels
)

refined_estimates <- extract_estimates(
  midpoint_results$subtype_comparison,
  refined_interval_labels
)

final_subtype_estimates <- rbind(
  primary_estimates,
  refined_estimates
)

comparison_labels <- c(
  "Luminal B vs Luminal A",
  "HER2-enriched vs Luminal A",
  "Basal-like vs Luminal A",
  "Normal-like vs Luminal A",
  "Claudin-low vs Luminal A"
)

estimate_order <- order(
  match(final_subtype_estimates$interval, final_interval_labels),
  match(final_subtype_estimates$comparison, comparison_labels)
)

final_subtype_estimates <- final_subtype_estimates[
  estimate_order,
  ,
  drop = FALSE
]

row.names(final_subtype_estimates) <- NULL

interval_row_counts <- table(
  factor(
    final_subtype_estimates$interval,
    levels = final_interval_labels
  )
)

if (
  nrow(final_subtype_estimates) != 25L ||
  any(interval_row_counts != 5L) ||
  anyNA(match(final_subtype_estimates$comparison, comparison_labels)) ||
  any(final_subtype_estimates$hazard_ratio <= 0) ||
  any(final_subtype_estimates$lower_95_ci <= 0) ||
  any(final_subtype_estimates$upper_95_ci <= 0) ||
  any(
    final_subtype_estimates$lower_95_ci >
    final_subtype_estimates$hazard_ratio
  ) ||
  any(
    final_subtype_estimates$upper_95_ci <
    final_subtype_estimates$hazard_ratio
  )
) {
  stop("The final subtype estimates failed validation.")
}

final_subtype_estimates$analysis <- "Final five-interval RFS model"
final_subtype_estimates$molecular_subtype <- sub(
  " vs Luminal A$",
  "",
  final_subtype_estimates$comparison
)
final_subtype_estimates$reference <- "Luminal A"
final_subtype_estimates$multiplicity_family <- (
  "Twenty-five interval-specific subtype comparisons"
)
final_subtype_estimates$holm_adjusted_p_value <- stats::p.adjust(
  final_subtype_estimates$p_value,
  method = "holm"
)
final_subtype_estimates$holm_significant <- (
  final_subtype_estimates$holm_adjusted_p_value < 0.05
)

final_subtype_estimates <- final_subtype_estimates[
  ,
  c(
    "analysis",
    "interval",
    "molecular_subtype",
    "reference",
    "comparison",
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci",
    "p_value",
    "multiplicity_family",
    "holm_adjusted_p_value",
    "holm_significant"
  )
]


# 8. Derive subtype-specific participant and event counts

extract_event_counts <- function(model, interval_label) {
  subtype_columns <- grep(
    "^molecular_subtype",
    colnames(model$x),
    value = TRUE
  )
  
  if (length(subtype_columns) != 5L) {
    stop("The model matrix lacks subtype coefficients for: ", interval_label)
  }
  
  subtype_matrix <- model$x[
    ,
    subtype_columns,
    drop = FALSE
  ]
  
  if (
    !all(subtype_matrix %in% c(0, 1)) ||
    any(rowSums(subtype_matrix) > 1)
  ) {
    stop("The subtype indicators are invalid for: ", interval_label)
  }
  
  subtype_labels <- sub(
    "^molecular_subtype",
    "",
    subtype_columns
  )
  
  subtype_character <- rep("Luminal A", model$n)
  
  for (index in seq_along(subtype_columns)) {
    subtype_character[subtype_matrix[, index] == 1] <- (
      subtype_labels[index]
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
  
  subtype <- factor(
    subtype_character,
    levels = subtype_levels
  )
  
  event <- as.integer(model$y[, ncol(model$y)])
  
  if (
    length(subtype) != model$n ||
    length(event) != model$n ||
    anyNA(subtype) ||
    !all(event %in% c(0L, 1L))
  ) {
    stop("Model data could not be recovered for: ", interval_label)
  }
  
  counts <- data.frame(
    analysis = "Final five-interval RFS model",
    interval = interval_label,
    molecular_subtype = subtype_levels,
    participants_at_risk = vapply(
      subtype_levels,
      function(level) {
        sum(subtype == level)
      },
      integer(1)
    ),
    rfs_events = vapply(
      subtype_levels,
      function(level) {
        sum(event[subtype == level])
      },
      integer(1)
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  if (
    sum(counts$participants_at_risk) != model$n ||
    sum(counts$rfs_events) != model$nevent
  ) {
    stop("Recovered event counts failed validation for: ", interval_label)
  }
  
  counts
}

final_event_counts <- do.call(
  rbind,
  Map(
    extract_event_counts,
    final_models,
    final_interval_labels
  )
)

row.names(final_event_counts) <- NULL

if (nrow(final_event_counts) != 30L) {
  stop("The final subtype event-count table should contain 30 rows.")
}


# 9. Display the final results

cat("\nFinal five-interval model overview:\n")

overview_print <- final_model_overview
overview_print$likelihood_ratio_chisq <- round(
  overview_print$likelihood_ratio_chisq,
  3
)
overview_print$p_value <- format.pval(
  final_model_overview$p_value,
  digits = 4,
  eps = 0.001
)
overview_print$holm_adjusted_p_value <- format.pval(
  final_model_overview$holm_adjusted_p_value,
  digits = 4,
  eps = 0.001
)

print(
  overview_print[
    ,
    c(
      "interval",
      "participants_at_risk",
      "rfs_events",
      "likelihood_ratio_chisq",
      "degrees_of_freedom",
      "p_value",
      "holm_adjusted_p_value",
      "holm_significant"
    )
  ],
  row.names = FALSE
)

cat("\nFinal interval-specific subtype estimates:\n")

estimate_print <- final_subtype_estimates

estimate_print[
  c("hazard_ratio", "lower_95_ci", "upper_95_ci")
] <- round(
  estimate_print[
    c("hazard_ratio", "lower_95_ci", "upper_95_ci")
  ],
  3
)

estimate_print$p_value <- format.pval(
  final_subtype_estimates$p_value,
  digits = 4,
  eps = 0.001
)

estimate_print$holm_adjusted_p_value <- format.pval(
  final_subtype_estimates$holm_adjusted_p_value,
  digits = 4,
  eps = 0.001
)

print(
  estimate_print[
    ,
    c(
      "interval",
      "comparison",
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci",
      "p_value",
      "holm_adjusted_p_value",
      "holm_significant"
    )
  ],
  row.names = FALSE
)

cat("\nSubtype-specific event counts:\n")
print(final_event_counts, row.names = FALSE)

cat(
  "\nHolm-significant global subtype tests:",
  sum(final_model_overview$holm_significant),
  "of 5.\n"
)

cat(
  "Holm-significant subtype comparisons:",
  sum(final_subtype_estimates$holm_significant),
  "of 25.\n"
)

cat(
  paste0(
    "The 3.5-year boundary is a transparent post-diagnostic midpoint ",
    "refinement. It must not be described as prespecified or as an ",
    "externally established clinical cutoff.\n"
  )
)


# 10. Export and validate the final results

table_directory <- file.path(project_root, "output", "tables")
data_directory <- file.path(project_root, "data-derived")

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)

tables <- list(
  rfs_final_interval_definitions = final_interval_definitions,
  rfs_final_model_overview = final_model_overview,
  rfs_final_subtype_hazard_ratios = final_subtype_estimates,
  rfs_final_subtype_event_counts = final_event_counts
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
  "rfs_final_interval_analysis.rds"
)

saveRDS(
  list(
    interval_definitions = final_interval_definitions,
    models = final_models,
    model_overview = final_model_overview,
    subtype_hazard_ratios = final_subtype_estimates,
    subtype_event_counts = final_event_counts,
    multiplicity_method = "Holm",
    global_test_family_size = 5L,
    subtype_comparison_family_size = 25L,
    boundary_note = paste(
      "The 3.5-year boundary was selected as the neutral midpoint",
      "of the diagnostically inadequate 2-to-5-year interval."
    )
  ),
  results_path
)

expected_rows <- c(
  rfs_final_interval_definitions = 5L,
  rfs_final_model_overview = 5L,
  rfs_final_subtype_hazard_ratios = 25L,
  rfs_final_subtype_event_counts = 30L
)

observed_rows <- vapply(
  table_paths,
  function(path) {
    nrow(utils::read.csv(path))
  },
  integer(1)
)

if (!identical(unname(observed_rows), unname(expected_rows))) {
  stop("At least one final table has an unexpected row count.")
}

if (!file.exists(results_path)) {
  stop("The final RFS interval-analysis RDS file was not created.")
}

cat("\nScript 22 completed successfully.\n")
cat("Final results:", normalizePath(results_path), "\n")
cat("Tables:\n")
cat(
  paste0(
    "  ",
    normalizePath(table_paths),
    collapse = "\n"
  ),
  "\n"
)