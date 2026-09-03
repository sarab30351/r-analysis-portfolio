# 24_rfs_boundary_diagnostics.R

# Purpose: Assess the proportional-hazards assumption in the four alternative boundary models fitted in Script 23.

# This script diagnoses the saved 2-to-3, 3-to-5, 2-to-4 and 4-to-5-year models. It does not refit them, select a boundary or modify the final five-interval analysis from Script 22.


# 1. Check packages and project location

required_packages <- "survival"

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


# 2. Load the saved alternative boundary models

model_results_path <- file.path(
  project_root,
  "data-derived",
  "rfs_boundary_sensitivity.rds"
)

if (!file.exists(model_results_path)) {
  stop(
    "Script 23 results were not found: ",
    model_results_path
  )
}

boundary_results_object <- readRDS(
  model_results_path
)

required_components <- c(
  "reference_boundary_year",
  "alternative_boundary_years",
  "boundary_definitions",
  "boundary_results",
  "model_overview"
)

missing_components <- setdiff(
  required_components,
  names(boundary_results_object)
)

if (length(missing_components) > 0L) {
  stop(
    "Script 23 results are missing: ",
    paste(missing_components, collapse = ", ")
  )
}

boundary_definitions <- (
  boundary_results_object$boundary_definitions
)

boundary_results <- (
  boundary_results_object$boundary_results
)

model_overview <- (
  boundary_results_object$model_overview
)

expected_interval_keys <- c(
  "two_to_three",
  "three_to_five",
  "two_to_four",
  "four_to_five"
)

expected_interval_labels <- c(
  "2 to 3 years",
  "3 to 5 years",
  "2 to 4 years",
  "4 to 5 years"
)

if (
  boundary_results_object$reference_boundary_year != 3.5 ||
  !identical(
    as.numeric(
      boundary_results_object$alternative_boundary_years
    ),
    c(3, 4)
  ) ||
  !identical(
    names(boundary_results),
    expected_interval_keys
  ) ||
  !identical(
    as.character(boundary_definitions$interval),
    expected_interval_labels
  ) ||
  nrow(model_overview) != 4L
) {
  stop(
    "The saved alternative boundary structure is unexpected."
  )
}


# 3. Validate the four saved models

validate_model <- function(
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
      "The saved result for ",
      interval_label,
      " is missing: ",
      paste(
        missing_result_components,
        collapse = ", "
      )
    )
  }
  
  model <- result$extended_model
  
  if (!inherits(model, "coxph")) {
    stop(
      "The saved object is not a Cox model for: ",
      interval_label
    )
  }
  
  overview_row <- model_overview[
    model_overview$interval == interval_label,
    ,
    drop = FALSE
  ]
  
  if (
    nrow(overview_row) != 1L ||
    result$interval != interval_label ||
    model$n != overview_row$participants_at_risk ||
    model$nevent != overview_row$rfs_events
  ) {
    stop(
      "The saved model does not match its overview for: ",
      interval_label
    )
  }
  
  if (is.null(model$x) || is.null(model$y)) {
    stop(
      "The model lacks data required for diagnostics: ",
      interval_label
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
      "The model has unexpected coefficients for: ",
      interval_label
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
      
      validate_model(
        result = boundary_results[[interval_key]],
        interval_key = interval_key,
        interval_label = expected_interval_labels[index]
      )
    }
  ),
  expected_interval_keys
)


# 4. Define diagnostic helpers

find_one <- function(
    pattern,
    values,
    label
) {
  matches <- grep(
    pattern,
    values
  )
  
  if (length(matches) != 1L) {
    stop(
      "Expected one match for ",
      label,
      " but found ",
      length(matches),
      "."
    )
  }
  
  matches
}

diagnose_model <- function(
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
  
  term_matrix <- ph_terms$table
  term_names <- rownames(
    term_matrix
  )
  
  term_rows <- c(
    find_one(
      "age_at_diagnosis",
      term_names,
      paste(
        interval_label,
        "age at diagnosis"
      )
    ),
    find_one(
      "npi",
      term_names,
      paste(
        interval_label,
        "NPI"
      )
    ),
    find_one(
      "molecular_subtype",
      term_names,
      paste(
        interval_label,
        "molecular subtype"
      )
    ),
    find_one(
      "^GLOBAL$",
      term_names,
      paste(
        interval_label,
        "global model"
      )
    )
  )
  
  term_tests <- data.frame(
    analysis = (
      "Alternative RFS boundary diagnostics"
    ),
    interval = interval_label,
    term = c(
      "Age-at-diagnosis spline",
      "NPI spline",
      "Molecular subtype",
      "Global model"
    ),
    chi_square = term_matrix[
      term_rows,
      "chisq"
    ],
    degrees_of_freedom = term_matrix[
      term_rows,
      "df"
    ],
    p_value = term_matrix[
      term_rows,
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
  
  if (
    !identical(
      as.integer(
        term_tests$degrees_of_freedom
      ),
      c(3L, 2L, 5L, 10L)
    )
  ) {
    stop(
      paste(
        "Unexpected term diagnostic degrees",
        "of freedom for:",
        interval_label
      )
    )
  }
  
  coefficient_matrix <- (
    ph_coefficients$table
  )
  
  subtype_rows <- match(
    subtype_terms,
    rownames(coefficient_matrix)
  )
  
  if (anyNA(subtype_rows)) {
    stop(
      "Subtype diagnostics could not be matched for: ",
      interval_label
    )
  }
  
  comparison_subtypes <- sub(
    "^molecular_subtype",
    "",
    subtype_terms
  )
  
  subtype_tests <- data.frame(
    analysis = (
      "Alternative RFS boundary diagnostics"
    ),
    interval = interval_label,
    term = subtype_terms,
    comparison = paste(
      comparison_subtypes,
      "vs Luminal A"
    ),
    chi_square = coefficient_matrix[
      subtype_rows,
      "chisq"
    ],
    degrees_of_freedom = coefficient_matrix[
      subtype_rows,
      "df"
    ],
    p_value = coefficient_matrix[
      subtype_rows,
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
  
  if (
    !all(
      subtype_tests$degrees_of_freedom == 1
    )
  ) {
    stop(
      paste(
        "Unexpected subtype diagnostic degrees",
        "of freedom for:",
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


# 5. Diagnose all four alternative models

boundary_diagnostics <- setNames(
  lapply(
    validated_models,
    diagnose_model
  ),
  expected_interval_keys
)

term_tests <- do.call(
  rbind,
  lapply(
    boundary_diagnostics,
    function(result) {
      result$term_tests
    }
  )
)

subtype_tests <- do.call(
  rbind,
  lapply(
    boundary_diagnostics,
    function(result) {
      result$subtype_tests
    }
  )
)

row.names(term_tests) <- NULL
row.names(subtype_tests) <- NULL

if (
  nrow(term_tests) != 16L ||
  nrow(subtype_tests) != 20L
) {
  stop(
    "The combined diagnostic tables have unexpected dimensions."
  )
}


# 6. Display the diagnostic results

format_diagnostic_table <- function(data) {
  result <- data
  
  result$chi_square <- round(
    result$chi_square,
    3
  )
  
  result$p_value <- format.pval(
    data$p_value,
    digits = 4,
    eps = 0.001
  )
  
  result
}

for (interval_label in expected_interval_labels) {
  cat(
    "\n",
    interval_label,
    "\n",
    sep = ""
  )
  
  interval_term_tests <- term_tests[
    term_tests$interval == interval_label,
    ,
    drop = FALSE
  ]
  
  interval_subtype_tests <- subtype_tests[
    subtype_tests$interval == interval_label,
    ,
    drop = FALSE
  ]
  
  cat(
    "\nTerm-level proportional hazards diagnostics:\n"
  )
  
  print(
    format_diagnostic_table(
      interval_term_tests
    )[
      ,
      c(
        "interval",
        "term",
        "chi_square",
        "degrees_of_freedom",
        "p_value",
        "diagnostic_flag"
      )
    ],
    row.names = FALSE
  )
  
  cat(
    paste0(
      "\nIndividual subtype ",
      "proportional hazards diagnostics:\n"
    )
  )
  
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
  term_tests$p_value < 0.05
)

subtype_signal_count <- sum(
  subtype_tests$p_value < 0.05
)

cat(
  "\nDiagnostic signals were found in ",
  term_signal_count,
  " of 16 term-level tests and ",
  subtype_signal_count,
  " of 20 individual subtype tests.\n",
  sep = ""
)

cat(
  paste0(
    "These are sensitivity model diagnostics. ",
    "They do not determine the final boundary ",
    "and are not added to the Holm-adjustment ",
    "families.\n"
  )
)


# 7. Create one diagnostic figure for each alternative model

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
      "rfs_boundary_ph_",
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

create_diagnostic_figure <- function(
    diagnostic_result,
    path
) {
  term_plot_names <- colnames(
    diagnostic_result$ph_terms$y
  )
  
  term_plot_columns <- c(
    find_one(
      "age_at_diagnosis",
      term_plot_names,
      paste(
        diagnostic_result$interval,
        "age plot"
      )
    ),
    find_one(
      "npi",
      term_plot_names,
      paste(
        diagnostic_result$interval,
        "NPI plot"
      )
    ),
    find_one(
      "molecular_subtype",
      term_plot_names,
      paste(
        diagnostic_result$interval,
        "subtype plot"
      )
    )
  )
  
  subtype_plot_columns <- match(
    diagnostic_result$subtype_terms,
    colnames(
      diagnostic_result$ph_coefficients$y
    )
  )
  
  if (anyNA(subtype_plot_columns)) {
    stop(
      "Subtype plots could not be matched for: ",
      diagnostic_result$interval
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
      graphics::par(
        old_graphics_parameters
      )
      
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
    "Age-at-diagnosis spline"
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
  
  for (
    index in seq_along(
      subtype_plot_columns
    )
  ) {
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
      "RFS alternative boundary diagnostics:",
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
    create_diagnostic_figure,
    boundary_diagnostics,
    figure_paths
  )
)


# 8. Export and validate the diagnostic results

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

tables <- list(
  rfs_boundary_ph_terms = term_tests,
  rfs_boundary_ph_subtypes = subtype_tests
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
  "rfs_boundary_diagnostics.rds"
)

saveRDS(
  list(
    reference_boundary_year = 3.5,
    alternative_boundary_years = c(3, 4),
    boundary_definitions = boundary_definitions,
    diagnostics = boundary_diagnostics,
    term_tests = term_tests,
    subtype_tests = subtype_tests,
    transform = "km",
    figure_paths = figure_paths,
    multiplicity_note = paste(
      "Diagnostic p-values are not included",
      "in the final Holm-adjustment families."
    )
  ),
  results_path
)

expected_rows <- c(
  rfs_boundary_ph_terms = 16L,
  rfs_boundary_ph_subtypes = 20L
)

observed_rows <- vapply(
  table_paths,
  function(path) {
    nrow(
      utils::read.csv(path)
    )
  },
  integer(1)
)

if (
  !identical(
    unname(observed_rows),
    unname(expected_rows)
  )
) {
  stop(
    "At least one diagnostic table has an unexpected row count."
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
    "At least one diagnostic figure was not created correctly."
  )
}

if (!file.exists(results_path)) {
  stop(
    "The boundary diagnostic RDS file was not created."
  )
}

cat(
  "\nScript 24 completed successfully.\n"
)

cat(
  "Diagnostic results:",
  normalizePath(results_path),
  "\n"
)

cat("Figures:\n")

cat(
  paste0(
    "  ",
    normalizePath(figure_paths),
    collapse = "\n"
  ),
  "\n"
)

cat("Tables:\n")

cat(
  paste0(
    "  ",
    normalizePath(table_paths),
    collapse = "\n"
  ),
  "\n"
)