# 21_rfs_two_to_five_split_sensitivity.R

# Purpose: Explore the residual molecular-subtype proportional-hazards signal within the primary 2-to-5-year RFS interval.

# The interval is divided at 3.5 years, its neutral midpoint. This boundary was not selected from an estimated change point or an individual subtype curve.

# The two models are exploratory sensitivity analyses. They do not replace the primary 2-to-5-year model automatically and their p-values are not included in the primary multiplicity-adjustment families.


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

required_diagnostic_components <- c(
  "term_tests",
  "subtype_tests"
)

missing_diagnostic_components <- setdiff(
  required_diagnostic_components,
  names(diagnostic_results)
)

if (length(missing_diagnostic_components) > 0L) {
  stop(
    "Script 19 results are missing: ",
    paste(missing_diagnostic_components, collapse = ", ")
  )
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


# 4. Recover the primary spline knots

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


# 5. Confirm the residual 2-to-5-year subtype signal

subtype_signal <- diagnostic_results$term_tests[
  diagnostic_results$term_tests$interval == "2 to 5 years" &
    diagnostic_results$term_tests$term == "Molecular subtype",
  ,
  drop = FALSE
]

if (nrow(subtype_signal) != 1L) {
  stop("The expected 2-to-5-year subtype diagnostic was not found.")
}

if (subtype_signal$p_value >= 0.05) {
  warning(
    paste(
      "The subtype diagnostic no longer matches the signal",
      "that motivated this sensitivity analysis."
    )
  )
}

flagged_comparisons <- diagnostic_results$subtype_tests$comparison[
  diagnostic_results$subtype_tests$interval == "2 to 5 years" &
    diagnostic_results$subtype_tests$p_value < 0.05
]

expected_flagged_comparisons <- c(
  "HER2-enriched vs Luminal A",
  "Normal-like vs Luminal A"
)

if (!setequal(flagged_comparisons, expected_flagged_comparisons)) {
  warning(
    paste(
      "The flagged individual comparisons differ from the",
      "HER2-enriched and Normal-like signals expected from Script 19."
    )
  )
}


# 6. Create the two exploratory risk sets

split_definitions <- data.frame(
  interval_key = c(
    "two_to_three_point_five",
    "three_point_five_to_five"
  ),
  interval = c(
    "2 to 3.5 years",
    "3.5 to 5 years"
  ),
  start_year = c(2, 3.5),
  end_year = c(3.5, 5),
  stringsAsFactors = FALSE
)

create_interval_data <- function(data, start_year, end_year, interval_label) {
  interval_data <- data[
    data$rfs_years > start_year,
    ,
    drop = FALSE
  ]
  
  interval_data$interval_time <- (
    pmin(interval_data$rfs_years, end_year) - start_year
  )
  
  interval_data$interval_event <- as.integer(
    interval_data$rfs_event == 1L &
      interval_data$rfs_years <= end_year
  )
  
  interval_data$interval <- interval_label
  
  if (
    nrow(interval_data) == 0L ||
    sum(interval_data$interval_event) == 0L ||
    any(interval_data$interval_time <= 0) ||
    any(table(interval_data$molecular_subtype) == 0L)
  ) {
    stop("The risk set failed validation for: ", interval_label)
  }
  
  interval_data
}

split_data <- setNames(
  lapply(
    seq_len(nrow(split_definitions)),
    function(index) {
      create_interval_data(
        data = rfs_data,
        start_year = split_definitions$start_year[index],
        end_year = split_definitions$end_year[index],
        interval_label = split_definitions$interval[index]
      )
    }
  ),
  split_definitions$interval_key
)

# An event recorded exactly at 3.5 years belongs to the first interval. Only
# participants observed event-free beyond 3.5 years enter the second interval.

primary_overview <- primary_results$model_overview[
  primary_results$model_overview$interval == "2 to 5 years",
  ,
  drop = FALSE
]

if (nrow(primary_overview) != 1L) {
  stop("The primary 2-to-5-year model overview was not found.")
}

split_participants <- vapply(split_data, nrow, integer(1))

split_events <- vapply(
  split_data,
  function(data) {
    sum(data$interval_event)
  },
  integer(1)
)

if (
  split_participants[1] != primary_overview$participants_at_risk ||
  split_participants[2] >= split_participants[1] ||
  sum(split_events) != primary_overview$rfs_events
) {
  stop("The split risk sets do not reproduce the primary interval.")
}


# 7. Define the shared model formulas and diagnostic helper

clinical_formula <- survival::Surv(
  interval_time,
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

find_one <- function(pattern, values, label) {
  matches <- grep(pattern, values)
  
  if (length(matches) != 1L) {
    stop("Expected one diagnostic row for ", label, ".")
  }
  
  matches
}


# 8. Fit and diagnose one exploratory interval

fit_interval <- function(interval_data, interval_label) {
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
  
  if (
    clinical_parameters != 5L ||
    extended_parameters != 10L ||
    likelihood_ratio_df != 5L
  ) {
    stop("Unexpected parameter counts for: ", interval_label)
  }
  
  likelihood_ratio_chisq <- 2 * (
    as.numeric(extended_loglik) - as.numeric(clinical_loglik)
  )
  
  model_overview <- data.frame(
    analysis = "2-to-5-year midpoint split",
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
  
  if (
    length(subtype_terms) != 5L ||
    anyNA(coefficient_table[subtype_terms, "coef"])
  ) {
    stop("The subtype coefficients are incomplete for: ", interval_label)
  }
  
  comparison_subtypes <- sub(
    "^molecular_subtype",
    "",
    subtype_terms
  )
  
  subtype_estimates <- data.frame(
    analysis = "2-to-5-year midpoint split",
    interval = interval_label,
    term = subtype_terms,
    comparison = paste(comparison_subtypes, "vs Luminal A"),
    hazard_ratio = confidence_table[subtype_terms, "exp(coef)"],
    lower_95_ci = confidence_table[subtype_terms, "lower .95"],
    upper_95_ci = confidence_table[subtype_terms, "upper .95"],
    p_value = coefficient_table[subtype_terms, "Pr(>|z|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  subtype_event_counts <- data.frame(
    analysis = "2-to-5-year midpoint split",
    interval = interval_label,
    molecular_subtype = subtype_levels,
    participants_at_risk = vapply(
      subtype_levels,
      function(subtype) {
        sum(interval_data$molecular_subtype == subtype)
      },
      integer(1)
    ),
    rfs_events = vapply(
      subtype_levels,
      function(subtype) {
        sum(
          interval_data$interval_event[
            interval_data$molecular_subtype == subtype
          ]
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  if (
    sum(subtype_event_counts$participants_at_risk) != nrow(interval_data) ||
    sum(subtype_event_counts$rfs_events) != expected_events
  ) {
    stop("Subtype event counts failed validation for: ", interval_label)
  }
  
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
  
  term_matrix <- ph_terms$table
  term_names <- rownames(term_matrix)
  
  term_rows <- c(
    find_one("age_at_diagnosis", term_names, "age at diagnosis"),
    find_one("npi", term_names, "NPI"),
    find_one("molecular_subtype", term_names, "molecular subtype"),
    find_one("^GLOBAL$", term_names, "the global model")
  )
  
  term_tests <- data.frame(
    diagnostic_level = "Model term",
    term = c(
      "Age-at-diagnosis spline",
      "NPI spline",
      "Molecular subtype",
      "Global model"
    ),
    comparison = NA_character_,
    chi_square = term_matrix[term_rows, "chisq"],
    degrees_of_freedom = term_matrix[term_rows, "df"],
    p_value = term_matrix[term_rows, "p"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  if (
    !identical(
      as.integer(term_tests$degrees_of_freedom),
      c(3L, 2L, 5L, 10L)
    )
  ) {
    stop(
      "Unexpected term-level diagnostic degrees of freedom for: ",
      interval_label
    )
  }
  
  coefficient_matrix <- ph_coefficients$table
  subtype_rows <- match(subtype_terms, rownames(coefficient_matrix))
  
  if (anyNA(subtype_rows)) {
    stop(
      "The subtype diagnostics could not be matched for: ",
      interval_label
    )
  }
  
  subtype_tests <- data.frame(
    diagnostic_level = "Subtype coefficient",
    term = subtype_terms,
    comparison = paste(comparison_subtypes, "vs Luminal A"),
    chi_square = coefficient_matrix[subtype_rows, "chisq"],
    degrees_of_freedom = coefficient_matrix[subtype_rows, "df"],
    p_value = coefficient_matrix[subtype_rows, "p"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  if (!all(subtype_tests$degrees_of_freedom == 1)) {
    stop(
      "Unexpected subtype diagnostic degrees of freedom for: ",
      interval_label
    )
  }
  
  ph_tests <- rbind(term_tests, subtype_tests)
  ph_tests$analysis <- "2-to-5-year midpoint split"
  ph_tests$interval <- interval_label
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
  
  list(
    interval = interval_label,
    interval_data = interval_data,
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_overview = model_overview,
    subtype_estimates = subtype_estimates,
    subtype_event_counts = subtype_event_counts,
    ph_terms = ph_terms,
    ph_coefficients = ph_coefficients,
    ph_tests = ph_tests
  )
}


# 9. Fit both exploratory intervals and combine the results

split_results <- setNames(
  lapply(
    seq_len(nrow(split_definitions)),
    function(index) {
      interval_key <- split_definitions$interval_key[index]
      
      fit_interval(
        interval_data = split_data[[interval_key]],
        interval_label = split_definitions$interval[index]
      )
    }
  ),
  split_definitions$interval_key
)

combine_component <- function(component) {
  do.call(
    rbind,
    lapply(
      split_results,
      function(result) {
        result[[component]]
      }
    )
  )
}

model_overview <- combine_component("model_overview")
subtype_estimates <- combine_component("subtype_estimates")
subtype_event_counts <- combine_component("subtype_event_counts")
ph_tests <- combine_component("ph_tests")

row.names(model_overview) <- NULL
row.names(subtype_estimates) <- NULL
row.names(subtype_event_counts) <- NULL
row.names(ph_tests) <- NULL

primary_estimates <- primary_results$subtype_hazard_ratios[
  primary_results$subtype_hazard_ratios$interval == "2 to 5 years",
  c("comparison", "hazard_ratio", "lower_95_ci", "upper_95_ci")
]

primary_rows <- match(
  subtype_estimates$comparison,
  primary_estimates$comparison
)

if (nrow(primary_estimates) != 5L || anyNA(primary_rows)) {
  stop(
    "The primary and exploratory subtype estimates could not be matched."
  )
}

subtype_comparison <- subtype_estimates
subtype_comparison$primary_2_to_5_hazard_ratio <- (
  primary_estimates$hazard_ratio[primary_rows]
)
subtype_comparison$primary_2_to_5_lower_95_ci <- (
  primary_estimates$lower_95_ci[primary_rows]
)
subtype_comparison$primary_2_to_5_upper_95_ci <- (
  primary_estimates$upper_95_ci[primary_rows]
)


# 10. Display the exploratory results

cat("\nMidpoint-split model overview:\n")

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
print(overview_print, row.names = FALSE)

cat("\nSubtype-specific events in each exploratory interval:\n")
print(subtype_event_counts, row.names = FALSE)

cat("\nExploratory subtype hazard ratios:\n")

subtype_print <- subtype_comparison
subtype_print[
  c(
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci",
    "primary_2_to_5_hazard_ratio"
  )
] <- round(
  subtype_print[
    c(
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci",
      "primary_2_to_5_hazard_ratio"
    )
  ],
  3
)
subtype_print$p_value <- format.pval(
  subtype_comparison$p_value,
  digits = 4,
  eps = 0.001
)
print(
  subtype_print[
    ,
    c(
      "interval",
      "comparison",
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci",
      "p_value",
      "primary_2_to_5_hazard_ratio"
    )
  ],
  row.names = FALSE
)

cat("\nProportional hazards diagnostics after the midpoint split:\n")

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
    "\nThe two hazard-ratio sets describe different follow-up periods. ",
    "They should not be judged against the primary 2-to-5-year average ",
    "using a percentage-change threshold. These analyses assess the ",
    "residual time trend and remain exploratory.\n"
  )
)


# 11. Export and validate the results

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
  rfs_two_to_five_split_model_overview = model_overview,
  rfs_two_to_five_split_sensitivity = subtype_comparison,
  rfs_two_to_five_split_event_counts = subtype_event_counts,
  rfs_two_to_five_split_ph_tests = ph_tests
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
  "rfs_two_to_five_split_sensitivity.rds"
)

saveRDS(
  list(
    split_definitions = split_definitions,
    split_results = split_results,
    model_overview = model_overview,
    subtype_comparison = subtype_comparison,
    subtype_event_counts = subtype_event_counts,
    ph_tests = ph_tests
  ),
  results_path
)

expected_rows <- c(
  rfs_two_to_five_split_model_overview = 2L,
  rfs_two_to_five_split_sensitivity = 10L,
  rfs_two_to_five_split_event_counts = 12L,
  rfs_two_to_five_split_ph_tests = 18L
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
  stop("The midpoint-split sensitivity RDS file was not created.")
}

cat("\nScript 21 completed successfully.\n")
cat(
  "Sensitivity results:",
  normalizePath(results_path),
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