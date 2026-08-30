# 17_rfs_full_followup_diagnostics.R

# Purpose: Assess the proportional hazards assumption for the full follow-up relapse-free survival model fitted in Script 16.

# This script diagnoses the existing extended model. It does not refit the model or treat the full follow-up hazard ratios as final.


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


# 2. Load and validate the full follow-up models

model_results_path <- file.path(
  project_root,
  "data-derived",
  "rfs_full_followup_cox_models.rds"
)

if (!file.exists(model_results_path)) {
  stop(
    paste(
      "Full follow-up RFS models not found:",
      model_results_path,
      "Run Script 16 first."
    )
  )
}

rfs_model_results <- readRDS(
  model_results_path
)

required_components <- c(
  "clinical_model",
  "extended_model"
)

missing_components <- setdiff(
  required_components,
  names(rfs_model_results)
)

if (length(missing_components) > 0L) {
  stop(
    paste(
      "The saved model results are missing:",
      paste(missing_components, collapse = ", ")
    )
  )
}

rfs_clinical_model <- (
  rfs_model_results$clinical_model
)

rfs_extended_model <- (
  rfs_model_results$extended_model
)

if (
  !inherits(rfs_clinical_model, "coxph") ||
  !inherits(rfs_extended_model, "coxph")
) {
  stop("The saved objects are not Cox proportional hazards models.")
}

if (
  rfs_extended_model$n != 1960L ||
  rfs_extended_model$nevent != 790L
) {
  stop(
    "The extended model does not match the validated RFS cohort."
  )
}

if (
  is.null(rfs_extended_model$x) ||
  is.null(rfs_extended_model$y)
) {
  stop(
    paste(
      "The extended model does not contain the data",
      "required for proportional hazards diagnostics."
    )
  )
}

coefficient_names <- names(
  stats::coef(rfs_extended_model)
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
    "The extended model does not contain the expected coefficients."
  )
}


# 3. Run the proportional hazards diagnostics

rfs_ph_terms <- survival::cox.zph(
  rfs_extended_model,
  transform = "km",
  terms = TRUE,
  singledf = FALSE,
  global = TRUE
)

rfs_ph_coefficients <- survival::cox.zph(
  rfs_extended_model,
  transform = "km",
  terms = FALSE,
  global = TRUE
)


# 4. Create the diagnostic tables

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

term_test_matrix <- (
  rfs_ph_terms$table
)

term_test_names <- rownames(
  term_test_matrix
)

term_row_indices <- c(
  find_single_match(
    "age_at_diagnosis",
    term_test_names,
    "the age spline"
  ),
  find_single_match(
    "npi",
    term_test_names,
    "the NPI spline"
  ),
  find_single_match(
    "molecular_subtype",
    term_test_names,
    "molecular subtype"
  ),
  find_single_match(
    "^GLOBAL$",
    term_test_names,
    "the global model"
  )
)

rfs_ph_term_tests <- data.frame(
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

rfs_ph_term_tests$diagnostic_flag <- ifelse(
  rfs_ph_term_tests$p_value < 0.05,
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
    as.integer(
      rfs_ph_term_tests$degrees_of_freedom
    ),
    expected_term_degrees_of_freedom
  )
) {
  stop(
    "The term-level diagnostics have unexpected degrees of freedom."
  )
}

coefficient_test_matrix <- (
  rfs_ph_coefficients$table
)

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
      "The subtype coefficients could not be matched",
      "to their diagnostic tests."
    )
  )
}

comparison_subtypes <- sub(
  "^molecular_subtype",
  "",
  subtype_terms
)

rfs_ph_subtype_tests <- data.frame(
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

rfs_ph_subtype_tests$diagnostic_flag <- ifelse(
  rfs_ph_subtype_tests$p_value < 0.05,
  "Evidence against proportional hazards",
  "No evidence against proportional hazards"
)

if (
  !all(
    rfs_ph_subtype_tests$degrees_of_freedom == 1
  )
) {
  stop(
    paste(
      "At least one subtype diagnostic has",
      "unexpected degrees of freedom."
    )
  )
}


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

cat(
  "\nTerm-level proportional hazards diagnostics:\n"
)

print(
  format_diagnostic_table(
    rfs_ph_term_tests
  ),
  row.names = FALSE
)

cat(
  "\nIndividual subtype proportional hazards diagnostics:\n"
)

print(
  format_diagnostic_table(
    rfs_ph_subtype_tests
  )[
    ,
    c(
      "comparison",
      "chi_square",
      "degrees_of_freedom",
      "p_value",
      "diagnostic_flag"
    )
  ],
  row.names = FALSE
)

subtype_ph_signal <- (
  rfs_ph_term_tests$p_value[
    rfs_ph_term_tests$term == "Molecular subtype"
  ] < 0.05
)

global_ph_signal <- (
  rfs_ph_term_tests$p_value[
    rfs_ph_term_tests$term == "Global model"
  ] < 0.05
)

if (subtype_ph_signal || global_ph_signal) {
  cat(
    paste0(
      "\nThe full follow-up model does not support constant ",
      "subtype hazard ratios across follow-up. The full ",
      "follow-up estimates remain provisional, and ",
      "interval-specific RFS models are required.\n"
    )
  )
}


# 6. Create the diagnostic figures

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

term_figure_path <- file.path(
  figure_directory,
  "rfs_full_followup_ph_terms.png"
)

subtype_figure_path <- file.path(
  figure_directory,
  "rfs_full_followup_ph_subtypes.png"
)

term_plot_names <- colnames(
  rfs_ph_terms$y
)

term_plot_columns <- c(
  find_single_match(
    "age_at_diagnosis",
    term_plot_names,
    "the age spline plot"
  ),
  find_single_match(
    "npi",
    term_plot_names,
    "the NPI spline plot"
  ),
  find_single_match(
    "molecular_subtype",
    term_plot_names,
    "the molecular subtype plot"
  )
)

subtype_plot_columns <- match(
  subtype_terms,
  colnames(rfs_ph_coefficients$y)
)

if (anyNA(subtype_plot_columns)) {
  stop(
    paste(
      "The subtype coefficients could not be matched",
      "to their diagnostic plots."
    )
  )
}

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
    xlab = "Years since diagnosis",
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

create_term_figure <- function(path) {
  grDevices::png(
    filename = path,
    width = 3000,
    height = 2200,
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
    mfrow = c(2, 2),
    mar = c(4.5, 4.5, 3.5, 1.5),
    oma = c(0, 0, 3.5, 0)
  )
  
  plot_ph_panel(
    rfs_ph_terms,
    term_plot_columns[1],
    "Age spline"
  )
  
  plot_ph_panel(
    rfs_ph_terms,
    term_plot_columns[2],
    "NPI spline"
  )
  
  plot_ph_panel(
    rfs_ph_terms,
    term_plot_columns[3],
    "Molecular subtype"
  )
  
  graphics::plot.new()
  
  graphics::mtext(
    paste(
      "RFS proportional hazards diagnostics:",
      "complete model terms"
    ),
    outer = TRUE,
    side = 3,
    line = 1,
    font = 2,
    cex = 1.4
  )
}

create_subtype_figure <- function(path) {
  grDevices::png(
    filename = path,
    width = 3000,
    height = 2400,
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
    mfrow = c(3, 2),
    mar = c(4.5, 4.5, 3.5, 1.5),
    oma = c(0, 0, 3.5, 0)
  )
  
  for (
    index in seq_along(
      subtype_plot_columns
    )
  ) {
    plot_ph_panel(
      rfs_ph_coefficients,
      subtype_plot_columns[index],
      paste(
        comparison_subtypes[index],
        "vs Luminal A"
      )
    )
  }
  
  graphics::plot.new()
  
  graphics::mtext(
    paste(
      "RFS proportional hazards diagnostics:",
      "subtype comparisons"
    ),
    outer = TRUE,
    side = 3,
    line = 1,
    font = 2,
    cex = 1.4
  )
}

create_term_figure(
  term_figure_path
)

create_subtype_figure(
  subtype_figure_path
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

rfs_diagnostic_tables <- list(
  rfs_full_followup_ph_terms =
    rfs_ph_term_tests,
  rfs_full_followup_ph_subtypes =
    rfs_ph_subtype_tests
)

rfs_diagnostic_table_paths <- file.path(
  table_directory,
  paste0(
    names(rfs_diagnostic_tables),
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
    rfs_diagnostic_tables,
    rfs_diagnostic_table_paths
  )
)

rfs_diagnostic_results <- list(
  ph_terms = rfs_ph_terms,
  ph_coefficients = rfs_ph_coefficients,
  term_tests = rfs_ph_term_tests,
  subtype_tests = rfs_ph_subtype_tests,
  transform = "km",
  term_figure_path = term_figure_path,
  subtype_figure_path = subtype_figure_path
)

rfs_diagnostic_results_path <- file.path(
  data_directory,
  "rfs_full_followup_diagnostics.rds"
)

saveRDS(
  rfs_diagnostic_results,
  rfs_diagnostic_results_path
)


# 8. Validate the exported results

expected_table_rows <- c(
  rfs_full_followup_ph_terms = 4L,
  rfs_full_followup_ph_subtypes = 5L
)

observed_table_rows <- vapply(
  rfs_diagnostic_table_paths,
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

figure_paths <- c(
  term_figure_path,
  subtype_figure_path
)

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

if (!file.exists(rfs_diagnostic_results_path)) {
  stop(
    "The diagnostic results file was not created."
  )
}

cat(
  "\nScript 17 completed successfully.\n"
)

cat(
  "Diagnostic results:",
  rfs_diagnostic_results_path,
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
    rfs_diagnostic_table_paths,
    collapse = "\n"
  ),
  "\n"
)