# 18_rfs_interval_cox_models.R

# Purpose: Fit interval-specific relapse-free survival Cox models 

# Intervals:
# 1. 0 to 2 years
# 2. 2 to 5 years
# 3. 5 to 10 years
# 4. Beyond 10 years

# The estimates remain provisional until the interval-specific diagnostics are completed in Script 19.


# 1. Check packages and project location

required_packages <- c(
  "survival",
  "splines"
)

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


# 2. Load the RFS cohort and Script 16 model specification

rfs_cohort_path <- file.path(
  project_root,
  "data-derived",
  "metabric_rfs_cohort.rds"
)

full_followup_model_path <- file.path(
  project_root,
  "data-derived",
  "rfs_full_followup_cox_models.rds"
)

if (!file.exists(rfs_cohort_path)) {
  stop(
    paste(
      "RFS cohort not found:",
      rfs_cohort_path
    )
  )
}

if (!file.exists(full_followup_model_path)) {
  stop(
    paste(
      "Script 16 model results not found:",
      full_followup_model_path
    )
  )
}

rfs_model_data <- as.data.frame(
  readRDS(rfs_cohort_path)
)

full_followup_results <- readRDS(
  full_followup_model_path
)

if (!"spline_knots" %in% names(full_followup_results)) {
  stop("The Script 16 results do not contain the spline knots.")
}

rfs_spline_knots <- (
  full_followup_results$spline_knots
)

required_knot_columns <- c(
  "variable",
  "percentile",
  "value"
)

if (
  !all(
    required_knot_columns %in%
    names(rfs_spline_knots)
  )
) {
  stop("The saved spline knot table has unexpected columns.")
}

age_knot_rows <- (
  rfs_spline_knots$variable == "Age at diagnosis"
)

npi_knot_rows <- (
  rfs_spline_knots$variable == "NPI"
)

age_knots <- rfs_spline_knots$value[
  age_knot_rows
][
  order(
    rfs_spline_knots$percentile[age_knot_rows]
  )
]

npi_knots <- rfs_spline_knots$value[
  npi_knot_rows
][
  order(
    rfs_spline_knots$percentile[npi_knot_rows]
  )
]

if (
  length(age_knots) != 4L ||
  length(npi_knots) != 3L ||
  anyDuplicated(age_knots) > 0L ||
  anyDuplicated(npi_knots) > 0L
) {
  stop("The saved spline knot specification is invalid.")
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

missing_columns <- setdiff(
  required_columns,
  names(rfs_model_data)
)

if (length(missing_columns) > 0L) {
  stop(
    paste(
      "The RFS cohort is missing:",
      paste(missing_columns, collapse = ", ")
    )
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

rfs_model_data$molecular_subtype <- factor(
  as.character(
    rfs_model_data$molecular_subtype
  ),
  levels = subtype_levels
)

contrasts(
  rfs_model_data$molecular_subtype
) <- stats::contr.treatment(
  subtype_levels,
  base = 1
)

rfs_model_data$cohort <- factor(
  rfs_model_data$cohort
)

rfs_model_data$rfs_years <- (
  rfs_model_data$rfs_months / 12
)

model_columns <- c(
  "patient_id",
  "rfs_years",
  "rfs_event",
  "age_at_diagnosis",
  "npi",
  "molecular_subtype",
  "cohort"
)

if (anyNA(rfs_model_data[model_columns])) {
  stop("The RFS model variables contain missing values.")
}

if (
  nrow(rfs_model_data) != 1960L ||
  sum(rfs_model_data$rfs_event) != 790L
) {
  stop("The data do not match the validated RFS cohort.")
}

if (
  any(rfs_model_data$rfs_years < 0) ||
  !all(rfs_model_data$rfs_event %in% c(0, 1))
) {
  stop("The RFS time or event coding is invalid.")
}

month_zero_events <- (
  rfs_model_data$rfs_years == 0 &
    rfs_model_data$rfs_event == 1
)

month_zero_censored <- (
  rfs_model_data$rfs_years == 0 &
    rfs_model_data$rfs_event == 0
)

if (
  sum(month_zero_events) != 3L ||
  any(month_zero_censored)
) {
  stop("The month zero RFS records differ from the validated cohort.")
}

if (
  levels(rfs_model_data$molecular_subtype)[1] !=
  "Luminal A"
) {
  stop("Luminal A is not the subtype reference level.")
}


# 4. Define the follow-up intervals

interval_definitions <- data.frame(
  interval_key = c(
    "zero_to_two",
    "two_to_five",
    "five_to_ten",
    "beyond_ten"
  ),
  interval = c(
    "0 to 2 years",
    "2 to 5 years",
    "5 to 10 years",
    "Beyond 10 years"
  ),
  start_year = c(
    0,
    2,
    5,
    10
  ),
  end_year = c(
    2,
    5,
    10,
    Inf
  ),
  stringsAsFactors = FALSE
)


# 5. Create the interval-specific risk sets

# Time is measured from the beginning of each interval. All participants in a later interval survived event-free beyond that interval's lower boundary.

# Resetting the time origin preserves the interval-specific risk sets and allows the three valid month zero events to remain at time zero without assigning artificial event times.

create_interval_data <- function(
    data,
    start_year,
    end_year,
    interval_label
) {
  if (start_year == 0) {
    at_risk <- data$rfs_years >= start_year
  } else {
    at_risk <- data$rfs_years > start_year
  }
  
  interval_data <- data[
    at_risk,
    ,
    drop = FALSE
  ]
  
  if (is.infinite(end_year)) {
    interval_exit <- interval_data$rfs_years
  } else {
    interval_exit <- pmin(
      interval_data$rfs_years,
      end_year
    )
  }
  
  interval_data$time_since_interval_start <- (
    interval_exit - start_year
  )
  
  interval_data$interval_event <- as.integer(
    interval_data$rfs_event == 1L &
      (
        is.infinite(end_year) |
          interval_data$rfs_years <= end_year
      )
  )
  
  interval_data$interval <- interval_label
  interval_data$interval_start_year <- start_year
  interval_data$interval_end_year <- end_year
  
  if (
    start_year == 0 &&
    any(
      interval_data$time_since_interval_start < 0
    )
  ) {
    stop(
      paste(
        "Negative follow-up found in interval:",
        interval_label
      )
    )
  }
  
  if (
    start_year > 0 &&
    any(
      interval_data$time_since_interval_start <= 0
    )
  ) {
    stop(
      paste(
        "Non-positive follow-up found in interval:",
        interval_label
      )
    )
  }
  
  if (
    any(
      interval_data$time_since_interval_start == 0 &
      interval_data$interval_event == 0
    )
  ) {
    stop(
      paste(
        "A zero-time censored record was found in interval:",
        interval_label
      )
    )
  }
  
  if (
    any(
      table(
        interval_data$molecular_subtype
      ) == 0L
    )
  ) {
    stop(
      paste(
        "At least one subtype is absent from interval:",
        interval_label
      )
    )
  }
  
  interval_data
}

interval_data_sets <- setNames(
  lapply(
    seq_len(
      nrow(interval_definitions)
    ),
    function(index) {
      create_interval_data(
        data = rfs_model_data,
        start_year = interval_definitions$start_year[index],
        end_year = interval_definitions$end_year[index],
        interval_label = interval_definitions$interval[index]
      )
    }
  ),
  interval_definitions$interval_key
)

interval_participant_counts <- vapply(
  interval_data_sets,
  nrow,
  FUN.VALUE = integer(1)
)

interval_event_counts <- vapply(
  interval_data_sets,
  function(data) {
    sum(data$interval_event)
  },
  FUN.VALUE = integer(1)
)

if (
  interval_participant_counts[1] != 1960L ||
  sum(interval_event_counts) != 790L ||
  any(diff(interval_participant_counts) >= 0)
) {
  stop("The interval-specific risk sets failed validation.")
}

first_interval_data <- (
  interval_data_sets$zero_to_two
)

if (
  sum(
    first_interval_data$time_since_interval_start == 0 &
    first_interval_data$interval_event == 1
  ) != 3L
) {
  stop("The first interval did not retain all month zero events.")
}


# 6. Define the interval-specific model formulas

rfs_interval_clinical_formula <- survival::Surv(
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

rfs_interval_extended_formula <- stats::update(
  rfs_interval_clinical_formula,
  . ~ . + molecular_subtype
)


# 7. Fit the models for one interval

fit_interval_models <- function(
    interval_data,
    interval_label,
    start_year,
    end_year
) {
  clinical_model <- survival::coxph(
    formula = rfs_interval_clinical_formula,
    data = interval_data,
    ties = "efron",
    na.action = stats::na.fail,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  extended_model <- survival::coxph(
    formula = rfs_interval_extended_formula,
    data = interval_data,
    ties = "efron",
    na.action = stats::na.fail,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  expected_event_count <- sum(
    interval_data$interval_event
  )
  
  if (
    clinical_model$n != nrow(interval_data) ||
    extended_model$n != nrow(interval_data) ||
    clinical_model$nevent != expected_event_count ||
    extended_model$nevent != expected_event_count
  ) {
    stop(
      paste(
        "A model did not use the full risk set in interval:",
        interval_label
      )
    )
  }
  
  clinical_loglik_object <- stats::logLik(
    clinical_model
  )
  
  extended_loglik_object <- stats::logLik(
    extended_model
  )
  
  clinical_loglik <- as.numeric(
    clinical_loglik_object
  )
  
  extended_loglik <- as.numeric(
    extended_loglik_object
  )
  
  clinical_parameter_count <- attr(
    clinical_loglik_object,
    "df"
  )
  
  extended_parameter_count <- attr(
    extended_loglik_object,
    "df"
  )
  
  likelihood_ratio_chisq <- 2 * (
    extended_loglik - clinical_loglik
  )
  
  likelihood_ratio_df <- (
    extended_parameter_count -
      clinical_parameter_count
  )
  
  likelihood_ratio_p_value <- stats::pchisq(
    likelihood_ratio_chisq,
    df = likelihood_ratio_df,
    lower.tail = FALSE
  )
  
  if (
    clinical_parameter_count != 5L ||
    extended_parameter_count != 10L ||
    likelihood_ratio_df != 5L
  ) {
    stop(
      paste(
        "Unexpected parameter count in interval:",
        interval_label
      )
    )
  }
  
  extended_summary <- summary(
    extended_model
  )
  
  coefficient_table <- (
    extended_summary$coefficients
  )
  
  confidence_table <- (
    extended_summary$conf.int
  )
  
  subtype_terms <- grep(
    "^molecular_subtype",
    rownames(coefficient_table),
    value = TRUE
  )
  
  if (length(subtype_terms) != 5L) {
    stop(
      paste(
        "Expected five subtype coefficients in interval:",
        interval_label
      )
    )
  }
  
  comparison_subtypes <- sub(
    "^molecular_subtype",
    "",
    subtype_terms
  )
  
  model_overview <- data.frame(
    interval = interval_label,
    start_year = start_year,
    end_year = end_year,
    participants_at_risk = nrow(interval_data),
    rfs_events = expected_event_count,
    clinical_parameters = clinical_parameter_count,
    extended_parameters = extended_parameter_count,
    stringsAsFactors = FALSE
  )
  
  global_subtype_test <- data.frame(
    interval = interval_label,
    participants_at_risk = nrow(interval_data),
    rfs_events = expected_event_count,
    clinical_log_likelihood = clinical_loglik,
    extended_log_likelihood = extended_loglik,
    likelihood_ratio_chisq = likelihood_ratio_chisq,
    degrees_of_freedom = likelihood_ratio_df,
    p_value = likelihood_ratio_p_value,
    stringsAsFactors = FALSE
  )
  
  subtype_hazard_ratios <- data.frame(
    interval = interval_label,
    term = subtype_terms,
    comparison = paste(
      comparison_subtypes,
      "vs Luminal A"
    ),
    log_hazard_ratio = coefficient_table[
      subtype_terms,
      "coef"
    ],
    standard_error = coefficient_table[
      subtype_terms,
      "se(coef)"
    ],
    hazard_ratio = confidence_table[
      subtype_terms,
      "exp(coef)"
    ],
    lower_95_ci = confidence_table[
      subtype_terms,
      "lower .95"
    ],
    upper_95_ci = confidence_table[
      subtype_terms,
      "upper .95"
    ],
    z_statistic = coefficient_table[
      subtype_terms,
      "z"
    ],
    p_value = coefficient_table[
      subtype_terms,
      "Pr(>|z|)"
    ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  subtype_event_counts <- do.call(
    rbind,
    lapply(
      subtype_levels,
      function(subtype) {
        subtype_rows <- (
          interval_data$molecular_subtype == subtype
        )
        
        data.frame(
          interval = interval_label,
          molecular_subtype = subtype,
          participants_at_risk = sum(subtype_rows),
          rfs_events = sum(
            interval_data$interval_event[subtype_rows]
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  if (
    sum(subtype_event_counts$participants_at_risk) !=
    nrow(interval_data) ||
    sum(subtype_event_counts$rfs_events) !=
    expected_event_count
  ) {
    stop(
      paste(
        "Subtype event counts failed validation in interval:",
        interval_label
      )
    )
  }
  
  list(
    interval = interval_label,
    start_year = start_year,
    end_year = end_year,
    interval_data = interval_data,
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_overview = model_overview,
    global_subtype_test = global_subtype_test,
    subtype_hazard_ratios = subtype_hazard_ratios,
    subtype_event_counts = subtype_event_counts
  )
}


# 8. Fit all four interval-specific models

rfs_interval_results <- setNames(
  lapply(
    seq_len(
      nrow(interval_definitions)
    ),
    function(index) {
      interval_key <- (
        interval_definitions$interval_key[index]
      )
      
      fit_interval_models(
        interval_data = interval_data_sets[[interval_key]],
        interval_label = interval_definitions$interval[index],
        start_year = interval_definitions$start_year[index],
        end_year = interval_definitions$end_year[index]
      )
    }
  ),
  interval_definitions$interval_key
)

rfs_interval_model_overview <- do.call(
  rbind,
  lapply(
    rfs_interval_results,
    function(result) {
      result$model_overview
    }
  )
)

rfs_interval_global_tests <- do.call(
  rbind,
  lapply(
    rfs_interval_results,
    function(result) {
      result$global_subtype_test
    }
  )
)

rfs_interval_subtype_hazard_ratios <- do.call(
  rbind,
  lapply(
    rfs_interval_results,
    function(result) {
      result$subtype_hazard_ratios
    }
  )
)

rfs_interval_subtype_event_counts <- do.call(
  rbind,
  lapply(
    rfs_interval_results,
    function(result) {
      result$subtype_event_counts
    }
  )
)

row.names(rfs_interval_model_overview) <- NULL
row.names(rfs_interval_global_tests) <- NULL
row.names(rfs_interval_subtype_hazard_ratios) <- NULL
row.names(rfs_interval_subtype_event_counts) <- NULL


# 9. Display the provisional interval-specific results

cat("\nRFS interval-specific model overview:\n")

print(
  rfs_interval_model_overview,
  row.names = FALSE
)

cat("\nSubtype-specific RFS event information:\n")

print(
  rfs_interval_subtype_event_counts,
  row.names = FALSE
)

cat("\nGlobal likelihood-ratio tests for molecular subtype:\n")

rfs_global_tests_print <- (
  rfs_interval_global_tests
)

rfs_global_tests_print[
  c(
    "clinical_log_likelihood",
    "extended_log_likelihood",
    "likelihood_ratio_chisq"
  )
] <- round(
  rfs_global_tests_print[
    c(
      "clinical_log_likelihood",
      "extended_log_likelihood",
      "likelihood_ratio_chisq"
    )
  ],
  3
)

rfs_global_tests_print$p_value <- format.pval(
  rfs_interval_global_tests$p_value,
  digits = 4,
  eps = 0.001
)

print(
  rfs_global_tests_print,
  row.names = FALSE
)

cat("\nProvisional interval-specific subtype hazard ratios:\n")

rfs_interval_hr_print <- (
  rfs_interval_subtype_hazard_ratios
)

rfs_interval_hr_print[
  c(
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci"
  )
] <- round(
  rfs_interval_hr_print[
    c(
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci"
    )
  ],
  3
)

rfs_interval_hr_print$p_value <- format.pval(
  rfs_interval_subtype_hazard_ratios$p_value,
  digits = 4,
  eps = 0.001
)

print(
  rfs_interval_hr_print[
    ,
    c(
      "interval",
      "comparison",
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci",
      "p_value"
    )
  ],
  row.names = FALSE
)

cat(
  paste0(
    "\nThese estimates are provisional. Interval-specific ",
    "proportional hazards diagnostics and multiplicity adjustment ",
    "must be completed before the results are interpreted as final.\n"
  )
)


# 10. Export the interval-specific model results

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

rfs_interval_tables <- list(
  rfs_interval_model_overview =
    rfs_interval_model_overview,
  rfs_interval_global_subtype_tests =
    rfs_interval_global_tests,
  rfs_interval_subtype_hazard_ratios =
    rfs_interval_subtype_hazard_ratios,
  rfs_interval_subtype_event_counts =
    rfs_interval_subtype_event_counts
)

rfs_interval_table_paths <- file.path(
  table_directory,
  paste0(
    names(rfs_interval_tables),
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
    rfs_interval_tables,
    rfs_interval_table_paths
  )
)

rfs_interval_model_results <- list(
  interval_definitions = interval_definitions,
  spline_knots = rfs_spline_knots,
  clinical_formula = rfs_interval_clinical_formula,
  extended_formula = rfs_interval_extended_formula,
  interval_results = rfs_interval_results,
  model_overview = rfs_interval_model_overview,
  global_subtype_tests = rfs_interval_global_tests,
  subtype_hazard_ratios =
    rfs_interval_subtype_hazard_ratios,
  subtype_event_counts =
    rfs_interval_subtype_event_counts
)

rfs_interval_results_path <- file.path(
  data_directory,
  "rfs_interval_cox_models.rds"
)

saveRDS(
  rfs_interval_model_results,
  rfs_interval_results_path
)


# 11. Validate the exported results

expected_table_rows <- c(
  rfs_interval_model_overview = 4L,
  rfs_interval_global_subtype_tests = 4L,
  rfs_interval_subtype_hazard_ratios = 20L,
  rfs_interval_subtype_event_counts = 24L
)

observed_table_rows <- vapply(
  rfs_interval_table_paths,
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
  stop("One or more exported tables have unexpected row counts.")
}

if (!file.exists(rfs_interval_results_path)) {
  stop("The interval-specific model results file was not created.")
}

cat("\nScript 18 completed successfully.\n")

cat(
  "Model results:",
  rfs_interval_results_path,
  "\n"
)

cat("Tables:\n")

cat(
  paste0(
    "  ",
    rfs_interval_table_paths,
    collapse = "\n"
  ),
  "\n"
)