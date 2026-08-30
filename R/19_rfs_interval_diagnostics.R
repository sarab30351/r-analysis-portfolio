# 19_rfs_interval_diagnostics.R

# Purpose: Assess the proportional hazards assumption within each interval-specific relapse-free survival model.

# This script diagnoses the saved extended models.


# 1. Check packages and project location

required_packages <- "survival"

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


# 2. Load and validate the interval-specific models

model_results_path <- file.path(
  project_root,
  "data-derived",
  "rfs_interval_cox_models.rds"
)

if (!file.exists(model_results_path)) {
  stop(
    paste(
      "Interval-specific RFS models not found:",
      model_results_path,
      "Run Script 18 first."
    )
  )
}

rfs_interval_model_results <- readRDS(
  model_results_path
)

required_components <- c(
  "interval_definitions",
  "interval_results",
  "model_overview"
)

missing_components <- setdiff(
  required_components,
  names(rfs_interval_model_results)
)

if (length(missing_components) > 0L) {
  stop(
    paste(
      "The saved model results are missing:",
      paste(missing_components, collapse = ", ")
    )
  )
}

interval_definitions <- (
  rfs_interval_model_results$interval_definitions
)

interval_results <- (
  rfs_interval_model_results$interval_results
)

model_overview <- (
  rfs_interval_model_results$model_overview
)

expected_interval_keys <- c(
  "zero_to_two",
  "two_to_five",
  "five_to_ten",
  "beyond_ten"
)

expected_interval_labels <- c(
  "0 to 2 years",
  "2 to 5 years",
  "5 to 10 years",
  "Beyond 10 years"
)

if (
  !identical(
    names(interval_results),
    expected_interval_keys
  ) ||
  !identical(
    as.character(interval_definitions$interval),
    expected_interval_labels
  ) ||
  nrow(model_overview) != 4L
) {
  stop("The saved interval structure is unexpected.")
}

expected_participant_counts <- setNames(
  model_overview$participants_at_risk,
  model_overview$interval
)

expected_event_counts <- setNames(
  model_overview$rfs_events,
  model_overview$interval
)

validate_interval_model <- function(
    result,
    interval_key,
    interval_label
) {
  required_result_components <- c(
    "interval",
    "extended_model"
  )
  
  missing_result_components <- setdiff(
    required_result_components,
    names(result)
  )
  
  if (length(missing_result_components) > 0L) {
    stop(
      paste(
        "The saved result for",
        interval_label,
        "is missing:",
        paste(missing_result_components, collapse = ", ")
      )
    )
  }
  
  model <- result$extended_model
  
  if (!inherits(model, "coxph")) {
    stop(
      paste(
        "The extended model is not a Cox model in interval:",
        interval_label
      )
    )
  }
  
  if (
    result$interval != interval_label ||
    model$n != expected_participant_counts[[interval_label]] ||
    model$nevent != expected_event_counts[[interval_label]]
  ) {
    stop(
      paste(
        "The saved model does not match the overview in interval:",
        interval_label
      )
    )
  }
  
  if (is.null(model$x) || is.null(model$y)) {
    stop(
      paste(
        "The model lacks data required for diagnostics in interval:",
        interval_label
      )
    )
  }
  
  coefficient_names <- names(
    stats::coef(model)
  )
  
  subtype_terms <- grep(
    "^molecular_subtype",
    coefficient_names,
    value = TRUE
  )
  
  if (
    length(coefficient_names) != 10L ||
    length(subtype_terms) != 5L
  ) {
    stop(
      paste(
        "The model has unexpected coefficients in interval:",
        interval_label
      )
    )
  }
  
  list(
    interval_key = interval_key,
    interval = interval_label,
    model = model,
    subtype_terms = subtype_terms
  )
}

validated_models <- setNames(
  lapply(
    seq_along(expected_interval_keys),
    function(index) {
      interval_key <- expected_interval_keys[index]
      
      validate_interval_model(
        result = interval_results[[interval_key]],
        interval_key = interval_key,
        interval_label = expected_interval_labels[index]
      )
    }
  ),
  expected_interval_keys
)


# 3. Define diagnostic helper functions

find_single_match <- function(
    pattern,
    values,
    diagnostic_label
) {
  matching_indices <- grep(
    pattern,
    values
  )
  
  if (length(matching_indices) != 1L) {
    stop(
      paste(
        "Expected one match for",
        diagnostic_label,
        "but found",
        length(matching_indices)
      )
    )
  }
  
  matching_indices
}

extract_interval_diagnostics <- function(
    validated_model
) {
  model <- validated_model$model
  interval_label <- validated_model$interval
  subtype_terms <- validated_model$subtype_terms
  
  ph_terms <- survival::cox.zph(
    model,
    transform = "km",
    terms = TRUE,
    singledf = FALSE,
    global = TRUE
  )
  
  ph_coefficients <- survival::cox.zph(
    model,
    transform = "km",
    terms = FALSE,
    global = TRUE
  )
  
  term_test_matrix <- ph_terms$table
  term_test_names <- rownames(term_test_matrix)
  
  term_row_indices <- c(
    find_single_match(
      "age_at_diagnosis",
      term_test_names,
      paste(interval_label, "age spline")
    ),
    find_single_match(
      "npi",
      term_test_names,
      paste(interval_label, "NPI spline")
    ),
    find_single_match(
      "molecular_subtype",
      term_test_names,
      paste(interval_label, "molecular subtype")
    ),
    find_single_match(
      "^GLOBAL$",
      term_test_names,
      paste(interval_label, "global model")
    )
  )
  
  term_tests <- data.frame(
    interval = interval_label,
    term = c(
      "Age spline",
      "NPI spline",
      "Molecular subtype",
      "Global model"
    ),
    chi_square = term_test_matrix[
      term_row_indices,
      "chisq"
    ],
    degrees_of_freedom = term_test_matrix[
      term_row_indices,
      "df"
    ],
    p_value = term_test_matrix[
      term_row_indices,
      "p"
    ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  term_tests$diagnostic_flag <- ifelse(
    term_tests$p_value < 0.05,
    "Evidence against proportional hazards",
    "No evidence against proportional hazards"
  )
  
  expected_term_degrees_of_freedom <- c(
    3L,
    2L,
    5L,
    10L
  )
  
  if (
    !identical(
      as.integer(term_tests$degrees_of_freedom),
      expected_term_degrees_of_freedom
    )
  ) {
    stop(
      paste(
        "The term diagnostics have unexpected degrees of freedom in interval:",
        interval_label
      )
    )
  }
  
  coefficient_test_matrix <- ph_coefficients$table
  coefficient_test_names <- rownames(
    coefficient_test_matrix
  )
  
  subtype_test_rows <- match(
    subtype_terms,
    coefficient_test_names
  )
  
  if (anyNA(subtype_test_rows)) {
    stop(
      paste(
        "Subtype coefficients could not be matched in interval:",
        interval_label
      )
    )
  }
  
  comparison_subtypes <- sub(
    "^molecular_subtype",
    "",
    subtype_terms
  )
  
  subtype_tests <- data.frame(
    interval = interval_label,
    term = subtype_terms,
    comparison = paste(
      comparison_subtypes,
      "vs Luminal A"
    ),
    chi_square = coefficient_test_matrix[
      subtype_test_rows,
      "chisq"
    ],
    degrees_of_freedom = coefficient_test_matrix[
      subtype_test_rows,
      "df"
    ],
    p_value = coefficient_test_matrix[
      subtype_test_rows,
      "p"
    ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  subtype_tests$diagnostic_flag <- ifelse(
    subtype_tests$p_value < 0.05,
    "Evidence against proportional hazards",
    "No evidence against proportional hazards"
  )
  
  if (!all(subtype_tests$degrees_of_freedom == 1)) {
    stop(
      paste(
        "Subtype diagnostics have unexpected degrees of freedom in interval:",
        interval_label
      )
    )
  }
  
  list(
    interval_key = validated_model$interval_key,
    interval = interval_label,
    ph_terms = ph_terms,
    ph_coefficients = ph_coefficients,
    term_tests = term_tests,
    subtype_tests = subtype_tests,
    subtype_terms = subtype_terms,
    comparison_subtypes = comparison_subtypes
  )
}


# 4. Run the diagnostics for all four intervals

rfs_interval_diagnostics <- setNames(
  lapply(
    validated_models,
    extract_interval_diagnostics
  ),
  expected_interval_keys
)

rfs_interval_ph_terms <- do.call(
  rbind,
  lapply(
    rfs_interval_diagnostics,
    function(result) {
      result$term_tests
    }
  )
)

rfs_interval_ph_subtypes <- do.call(
  rbind,
  lapply(
    rfs_interval_diagnostics,
    function(result) {
      result$subtype_tests
    }
  )
)

row.names(rfs_interval_ph_terms) <- NULL
row.names(rfs_interval_ph_subtypes) <- NULL


# 5. Display the diagnostic results

format_diagnostic_table <- function(data) {
  data_print <- data
  
  data_print$chi_square <- round(
    data_print$chi_square,
    3
  )
  
  data_print$p_value <- format.pval(
    data$p_value,
    digits = 4,
    eps = 0.001
  )
  
  data_print
}

for (interval_label in expected_interval_labels) {
  cat(
    "\n",
    interval_label,
    "\n",
    sep = ""
  )
  
  interval_term_tests <- rfs_interval_ph_terms[
    rfs_interval_ph_terms$interval == interval_label,
    ,
    drop = FALSE
  ]
  
  interval_subtype_tests <- rfs_interval_ph_subtypes[
    rfs_interval_ph_subtypes$interval == interval_label,
    ,
    drop = FALSE
  ]
  
  cat("\nTerm-level proportional hazards diagnostics:\n")
  
  print(
    format_diagnostic_table(
      interval_term_tests
    ),
    row.names = FALSE
  )
  
  cat("\nIndividual subtype proportional hazards diagnostics:\n")
  
  print(
    format_diagnostic_table(
      interval_subtype_tests
    )[
      ,
      c(
        "interval",
        "comparison",
        "chi_square",
        "degrees_of_freedom",
        "p_value",
        "diagnostic_flag"
      )
    ],
    row.names = FALSE
  )
}

term_signal_count <- sum(
  rfs_interval_ph_terms$p_value < 0.05
)

subtype_signal_count <- sum(
  rfs_interval_ph_subtypes$p_value < 0.05
)

cat(
  paste0(
    "\nDiagnostic signals were found in ",
    term_signal_count,
    " of 16 term-level tests and ",
    subtype_signal_count,
    " of 20 individual subtype tests.\n"
  )
)

cat(
  paste0(
    "These p-values are diagnostic signals rather than final ",
    "outcome tests. The curves, event counts and size of any ",
    "time trend must be assessed before deciding whether the ",
    "four intervals are adequate.\n"
  )
)


# 6. Create one diagnostic figure for each interval

figure_directory <- file.path(
  project_root,
  "output",
  "figures"
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_paths <- setNames(
  file.path(
    figure_directory,
    paste0(
      "rfs_interval_ph_",
      expected_interval_keys,
      ".png"
    )
  ),
  expected_interval_keys
)

plot_ph_panel <- function(
    ph_object,
    plot_column,
    panel_title
) {
  graphics::plot(
    ph_object,
    var = plot_column,
    resid = TRUE,
    se = TRUE,
    xlab = "Years since interval start",
    ylab = "Time-varying coefficient",
    main = panel_title,
    lwd = 2
  )
  
  graphics::abline(
    h = 0,
    col = "gray55",
    lty = 2,
    lwd = 1.5
  )
}

create_interval_figure <- function(
    diagnostic_result,
    path
) {
  term_plot_names <- colnames(
    diagnostic_result$ph_terms$y
  )
  
  term_plot_columns <- c(
    find_single_match(
      "age_at_diagnosis",
      term_plot_names,
      paste(diagnostic_result$interval, "age plot")
    ),
    find_single_match(
      "npi",
      term_plot_names,
      paste(diagnostic_result$interval, "NPI plot")
    ),
    find_single_match(
      "molecular_subtype",
      term_plot_names,
      paste(diagnostic_result$interval, "subtype term plot")
    )
  )
  
  subtype_plot_columns <- match(
    diagnostic_result$subtype_terms,
    colnames(diagnostic_result$ph_coefficients$y)
  )
  
  if (anyNA(subtype_plot_columns)) {
    stop(
      paste(
        "Subtype coefficients could not be matched to plots in interval:",
        diagnostic_result$interval
      )
    )
  }
  
  grDevices::png(
    filename = path,
    width = 3200,
    height = 3200,
    res = 200
  )
  
  old_graphics_parameters <- graphics::par(
    no.readonly = TRUE
  )
  
  on.exit(
    {
      graphics::par(old_graphics_parameters)
      grDevices::dev.off()
    },
    add = TRUE
  )
  
  graphics::par(
    mfrow = c(3, 3),
    mar = c(4.2, 4.2, 3.3, 1.3),
    oma = c(0, 0, 4, 0)
  )
  
  plot_ph_panel(
    diagnostic_result$ph_terms,
    term_plot_columns[1],
    "Age spline"
  )
  
  plot_ph_panel(
    diagnostic_result$ph_terms,
    term_plot_columns[2],
    "NPI spline"
  )
  
  plot_ph_panel(
    diagnostic_result$ph_terms,
    term_plot_columns[3],
    "Molecular subtype"
  )
  
  for (index in seq_along(subtype_plot_columns)) {
    plot_ph_panel(
      diagnostic_result$ph_coefficients,
      subtype_plot_columns[index],
      paste(
        diagnostic_result$comparison_subtypes[index],
        "vs Luminal A"
      )
    )
  }
  
  graphics::plot.new()
  
  graphics::mtext(
    paste(
      "RFS proportional hazards diagnostics:",
      diagnostic_result$interval
    ),
    outer = TRUE,
    side = 3,
    line = 1.5,
    font = 2,
    cex = 1.4
  )
}

invisible(
  Map(
    create_interval_figure,
    rfs_interval_diagnostics,
    figure_paths
  )
)


# 7. Export the diagnostic results

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

rfs_interval_diagnostic_tables <- list(
  rfs_interval_ph_terms = rfs_interval_ph_terms,
  rfs_interval_ph_subtypes = rfs_interval_ph_subtypes
)

rfs_interval_diagnostic_table_paths <- file.path(
  table_directory,
  paste0(
    names(rfs_interval_diagnostic_tables),
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
    rfs_interval_diagnostic_tables,
    rfs_interval_diagnostic_table_paths
  )
)

rfs_interval_diagnostic_results <- list(
  interval_definitions = interval_definitions,
  diagnostics = rfs_interval_diagnostics,
  term_tests = rfs_interval_ph_terms,
  subtype_tests = rfs_interval_ph_subtypes,
  transform = "km",
  figure_paths = figure_paths
)

rfs_interval_diagnostic_results_path <- file.path(
  data_directory,
  "rfs_interval_diagnostics.rds"
)

saveRDS(
  rfs_interval_diagnostic_results,
  rfs_interval_diagnostic_results_path
)


# 8. Validate the exported results

expected_table_rows <- c(
  rfs_interval_ph_terms = 16L,
  rfs_interval_ph_subtypes = 20L
)

observed_table_rows <- vapply(
  rfs_interval_diagnostic_table_paths,
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
  stop(
    "One or more diagnostic tables have unexpected row counts."
  )
}

figure_information <- file.info(
  figure_paths
)

if (
  any(!file.exists(figure_paths)) ||
  any(is.na(figure_information$size)) ||
  any(figure_information$size <= 0)
) {
  stop(
    "One or more diagnostic figures were not created correctly."
  )
}

if (!file.exists(rfs_interval_diagnostic_results_path)) {
  stop(
    "The diagnostic results file was not created."
  )
}

cat("\nScript 19 completed successfully.\n")

cat(
  "Diagnostic results:",
  rfs_interval_diagnostic_results_path,
  "\n"
)

cat("Figures:\n")

cat(
  paste0(
    "  ",
    figure_paths,
    collapse = "\n"
  ),
  "\n"
)

cat("Tables:\n")

cat(
  paste0(
    "  ",
    rfs_interval_diagnostic_table_paths,
    collapse = "\n"
  ),
  "\n"
)