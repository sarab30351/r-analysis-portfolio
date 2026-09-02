# 23_rfs_boundary_sensitivity.R

# Purpose: Assess whether the final RFS conclusions depend strongly on placing the post-diagnostic boundary at exactly 3.5 years.

# The 2-to-5-year period is refitted using alternative boundaries at 3.0 and 4.0 years. These models are sensitivity analyses; they do not replace the final five-interval analysis from Script 22. Their p-values are not included in the Holm-adjustment families and no boundary is selected by its p-value.


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


# 2. Load the cohort and saved analyses

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
  "spline_knots",
  "model_overview"
)

required_final_components <- c(
  "interval_definitions",
  "subtype_hazard_ratios"
)

if (
  any(!required_primary_components %in% names(primary_results)) ||
  any(!required_final_components %in% names(final_results))
) {
  stop("A saved RFS analysis has an unexpected structure.")
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

comparison_labels <- paste(
  subtype_levels[-1],
  "vs Luminal A"
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


# 4. Recover the primary spline specification

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


# 5. Define the alternative interval boundaries

boundary_definitions <- data.frame(
  interval_key = c(
    "two_to_three",
    "three_to_five",
    "two_to_four",
    "four_to_five"
  ),
  boundary_year = c(3, 3, 4, 4),
  phase = c(
    "Earlier portion",
    "Later portion",
    "Earlier portion",
    "Later portion"
  ),
  interval = c(
    "2 to 3 years",
    "3 to 5 years",
    "2 to 4 years",
    "4 to 5 years"
  ),
  start_year = c(2, 3, 2, 4),
  end_year = c(3, 5, 4, 5),
  stringsAsFactors = FALSE
)

create_interval_data <- function(
    data,
    start_year,
    end_year,
    interval_label
) {
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

boundary_data <- setNames(
  lapply(
    seq_len(nrow(boundary_definitions)),
    function(index) {
      create_interval_data(
        data = rfs_data,
        start_year = boundary_definitions$start_year[index],
        end_year = boundary_definitions$end_year[index],
        interval_label = boundary_definitions$interval[index]
      )
    }
  ),
  boundary_definitions$interval_key
)

# An event recorded exactly at an alternative boundary belongs to the earlier
# interval. Only participants observed event-free beyond it enter the later one.

primary_overview <- primary_results$model_overview[
  primary_results$model_overview$interval == "2 to 5 years",
  ,
  drop = FALSE
]

if (nrow(primary_overview) != 1L) {
  stop("The primary 2-to-5-year model overview was not found.")
}

for (boundary in c(3, 4)) {
  rows <- boundary_definitions$boundary_year == boundary
  keys <- boundary_definitions$interval_key[rows]
  
  participant_counts <- vapply(
    boundary_data[keys],
    nrow,
    integer(1)
  )
  
  event_counts <- vapply(
    boundary_data[keys],
    function(data) {
      sum(data$interval_event)
    },
    integer(1)
  )
  
  if (
    participant_counts[1] !=
    primary_overview$participants_at_risk ||
    participant_counts[2] >= participant_counts[1] ||
    sum(event_counts) != primary_overview$rfs_events
  ) {
    stop("The risk sets failed validation for boundary: ", boundary)
  }
}


# 6. Define the shared models

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


# 7. Fit one alternative interval

fit_interval <- function(interval_data, definition) {
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
    stop(
      "A model did not use the complete risk set for: ",
      definition$interval
    )
  }
  
  clinical_loglik <- stats::logLik(clinical_model)
  extended_loglik <- stats::logLik(extended_model)
  
  clinical_parameters <- attr(
    clinical_loglik,
    "df"
  )
  
  extended_parameters <- attr(
    extended_loglik,
    "df"
  )
  
  likelihood_ratio_df <- (
    extended_parameters - clinical_parameters
  )
  
  if (
    clinical_parameters != 5L ||
    extended_parameters != 10L ||
    likelihood_ratio_df != 5L
  ) {
    stop(
      "Unexpected parameter counts for: ",
      definition$interval
    )
  }
  
  likelihood_ratio_chisq <- 2 * (
    as.numeric(extended_loglik) -
      as.numeric(clinical_loglik)
  )
  
  identifying_columns <- data.frame(
    analysis = "Alternative RFS boundary sensitivity",
    boundary_year = definition$boundary_year,
    phase = definition$phase,
    interval = definition$interval,
    stringsAsFactors = FALSE
  )
  
  model_overview <- cbind(
    identifying_columns,
    data.frame(
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
    stop(
      "The subtype coefficients are incomplete for: ",
      definition$interval
    )
  }
  
  comparison_subtypes <- sub(
    "^molecular_subtype",
    "",
    subtype_terms
  )
  
  subtype_estimates <- cbind(
    identifying_columns[
      rep(1, 5),
      ,
      drop = FALSE
    ],
    data.frame(
      term = subtype_terms,
      comparison = paste(
        comparison_subtypes,
        "vs Luminal A"
      ),
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
      p_value = coefficient_table[
        subtype_terms,
        "Pr(>|z|)"
      ],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  )
  
  subtype_event_counts <- cbind(
    identifying_columns[
      rep(1, 6),
      ,
      drop = FALSE
    ],
    data.frame(
      molecular_subtype = subtype_levels,
      participants_at_risk = vapply(
        subtype_levels,
        function(subtype) {
          sum(
            interval_data$molecular_subtype ==
              subtype
          )
        },
        integer(1)
      ),
      rfs_events = vapply(
        subtype_levels,
        function(subtype) {
          sum(
            interval_data$interval_event[
              interval_data$molecular_subtype ==
                subtype
            ]
          )
        },
        integer(1)
      ),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  )
  
  subtype_event_counts$precision_flag <- ifelse(
    subtype_event_counts$rfs_events < 10L,
    "Fewer than 10 events",
    "10 or more events"
  )
  
  if (
    sum(
      subtype_event_counts$participants_at_risk
    ) != nrow(interval_data) ||
    sum(
      subtype_event_counts$rfs_events
    ) != expected_events
  ) {
    stop(
      "Subtype event counts failed validation for: ",
      definition$interval
    )
  }
  
  list(
    interval = definition$interval,
    interval_data = interval_data,
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_overview = model_overview,
    subtype_estimates = subtype_estimates,
    subtype_event_counts = subtype_event_counts
  )
}


# 8. Fit all four alternative intervals and combine the results

boundary_results <- setNames(
  lapply(
    seq_len(nrow(boundary_definitions)),
    function(index) {
      interval_key <- (
        boundary_definitions$interval_key[index]
      )
      
      fit_interval(
        interval_data = boundary_data[[interval_key]],
        definition = boundary_definitions[
          index,
          ,
          drop = FALSE
        ]
      )
    }
  ),
  boundary_definitions$interval_key
)

combine_component <- function(component) {
  result <- do.call(
    rbind,
    lapply(
      boundary_results,
      function(interval_result) {
        interval_result[[component]]
      }
    )
  )
  
  row.names(result) <- NULL
  
  result
}

model_overview <- combine_component(
  "model_overview"
)

subtype_estimates <- combine_component(
  "subtype_estimates"
)

subtype_event_counts <- combine_component(
  "subtype_event_counts"
)

if (
  nrow(model_overview) != 4L ||
  nrow(subtype_estimates) != 20L ||
  nrow(subtype_event_counts) != 24L
) {
  stop(
    paste(
      "The combined boundary-sensitivity results",
      "have unexpected dimensions."
    )
  )
}


# 9. Compare the alternative patterns with the final 3.5-year split

required_estimate_columns <- c(
  "interval",
  "comparison",
  "hazard_ratio"
)

if (
  any(
    !required_estimate_columns %in%
    names(
      final_results$subtype_hazard_ratios
    )
  )
) {
  stop(
    "The final subtype-estimate table has unexpected columns."
  )
}

reference_estimates <- (
  final_results$subtype_hazard_ratios[
    final_results$subtype_hazard_ratios$interval %in%
      c(
        "2 to 3.5 years",
        "3.5 to 5 years"
      ),
    required_estimate_columns,
    drop = FALSE
  ]
)

if (
  nrow(reference_estimates) != 10L ||
  any(
    table(reference_estimates$interval) != 5L
  )
) {
  stop(
    "The final 3.5-year reference estimates could not be recovered."
  )
}

classify_null_direction <- function(
    hazard_ratio
) {
  ifelse(
    hazard_ratio > 1,
    "Above 1",
    ifelse(
      hazard_ratio < 1,
      "Below 1",
      "Equal to 1"
    )
  )
}

classify_time_pattern <- function(
    early_hazard_ratio,
    late_hazard_ratio
) {
  ifelse(
    late_hazard_ratio > early_hazard_ratio,
    "Higher HR later",
    ifelse(
      late_hazard_ratio < early_hazard_ratio,
      "Lower HR later",
      "No change"
    )
  )
}

build_pattern_comparison <- function(boundary) {
  final_early <- reference_estimates[
    reference_estimates$interval ==
      "2 to 3.5 years",
    ,
    drop = FALSE
  ]
  
  final_late <- reference_estimates[
    reference_estimates$interval ==
      "3.5 to 5 years",
    ,
    drop = FALSE
  ]
  
  alternative <- subtype_estimates[
    subtype_estimates$boundary_year ==
      boundary,
    ,
    drop = FALSE
  ]
  
  alternative_early <- alternative[
    alternative$phase == "Earlier portion",
    ,
    drop = FALSE
  ]
  
  alternative_late <- alternative[
    alternative$phase == "Later portion",
    ,
    drop = FALSE
  ]
  
  final_early_rows <- match(
    comparison_labels,
    final_early$comparison
  )
  
  final_late_rows <- match(
    comparison_labels,
    final_late$comparison
  )
  
  alternative_early_rows <- match(
    comparison_labels,
    alternative_early$comparison
  )
  
  alternative_late_rows <- match(
    comparison_labels,
    alternative_late$comparison
  )
  
  if (
    anyNA(final_early_rows) ||
    anyNA(final_late_rows) ||
    anyNA(alternative_early_rows) ||
    anyNA(alternative_late_rows)
  ) {
    stop(
      "Subtype patterns could not be matched for boundary: ",
      boundary
    )
  }
  
  final_early_hr <- (
    final_early$hazard_ratio[
      final_early_rows
    ]
  )
  
  final_late_hr <- (
    final_late$hazard_ratio[
      final_late_rows
    ]
  )
  
  alternative_early_hr <- (
    alternative_early$hazard_ratio[
      alternative_early_rows
    ]
  )
  
  alternative_late_hr <- (
    alternative_late$hazard_ratio[
      alternative_late_rows
    ]
  )
  
  final_pattern <- classify_time_pattern(
    final_early_hr,
    final_late_hr
  )
  
  alternative_pattern <- classify_time_pattern(
    alternative_early_hr,
    alternative_late_hr
  )
  
  data.frame(
    alternative_boundary_year = boundary,
    comparison = comparison_labels,
    final_3_5_early_hazard_ratio = (
      final_early_hr
    ),
    final_3_5_late_hazard_ratio = (
      final_late_hr
    ),
    alternative_early_hazard_ratio = (
      alternative_early_hr
    ),
    alternative_late_hazard_ratio = (
      alternative_late_hr
    ),
    final_late_to_early_ratio = (
      final_late_hr / final_early_hr
    ),
    alternative_late_to_early_ratio = (
      alternative_late_hr /
        alternative_early_hr
    ),
    final_time_pattern = final_pattern,
    alternative_time_pattern = (
      alternative_pattern
    ),
    time_pattern_agreement = (
      final_pattern == alternative_pattern
    ),
    early_null_direction_agreement = (
      classify_null_direction(
        final_early_hr
      ) ==
        classify_null_direction(
          alternative_early_hr
        )
    ),
    late_null_direction_agreement = (
      classify_null_direction(
        final_late_hr
      ) ==
        classify_null_direction(
          alternative_late_hr
        )
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

pattern_comparison <- do.call(
  rbind,
  lapply(
    c(3, 4),
    build_pattern_comparison
  )
)

row.names(pattern_comparison) <- NULL

if (nrow(pattern_comparison) != 10L) {
  stop(
    "The boundary-pattern comparison should contain 10 rows."
  )
}


# 10. Display the sensitivity results

cat(
  "\nAlternative-boundary model overview:\n"
)

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
      "boundary_year",
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

cat(
  "\nAlternative-boundary subtype estimates:\n"
)

estimate_print <- subtype_estimates

estimate_print[
  c(
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci"
  )
] <- round(
  estimate_print[
    c(
      "hazard_ratio",
      "lower_95_ci",
      "upper_95_ci"
    )
  ],
  3
)

estimate_print$p_value <- format.pval(
  subtype_estimates$p_value,
  digits = 4,
  eps = 0.001
)

print(
  estimate_print[
    ,
    c(
      "boundary_year",
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
  "\nComparison with the final 3.5-year pattern:\n"
)

pattern_print <- pattern_comparison

ratio_columns <- c(
  "final_3_5_early_hazard_ratio",
  "final_3_5_late_hazard_ratio",
  "alternative_early_hazard_ratio",
  "alternative_late_hazard_ratio",
  "final_late_to_early_ratio",
  "alternative_late_to_early_ratio"
)

pattern_print[ratio_columns] <- round(
  pattern_print[ratio_columns],
  3
)

print(
  pattern_print,
  row.names = FALSE
)

sparse_counts <- subtype_event_counts[
  subtype_event_counts$rfs_events < 10L,
  ,
  drop = FALSE
]

cat(
  "\nSubtype cells with fewer than 10 RFS events:\n"
)

if (nrow(sparse_counts) == 0L) {
  cat("None.\n")
} else {
  print(
    sparse_counts[
      ,
      c(
        "boundary_year",
        "interval",
        "molecular_subtype",
        "participants_at_risk",
        "rfs_events"
      )
    ],
    row.names = FALSE
  )
}

cat(
  paste0(
    "\nThese alternative boundaries test whether the ",
    "substantive pattern depends on choosing exactly ",
    "3.5 years. Their p-values are descriptive ",
    "sensitivity results and are not added to the ",
    "final Holm families.\n"
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
  rfs_boundary_sensitivity_model_overview = (
    model_overview
  ),
  rfs_boundary_sensitivity_subtype_estimates = (
    subtype_estimates
  ),
  rfs_boundary_sensitivity_event_counts = (
    subtype_event_counts
  ),
  rfs_boundary_sensitivity_pattern_comparison = (
    pattern_comparison
  )
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
  "rfs_boundary_sensitivity.rds"
)

saveRDS(
  list(
    reference_boundary_year = 3.5,
    alternative_boundary_years = c(3, 4),
    boundary_definitions = boundary_definitions,
    boundary_results = boundary_results,
    model_overview = model_overview,
    subtype_estimates = subtype_estimates,
    subtype_event_counts = (
      subtype_event_counts
    ),
    pattern_comparison = pattern_comparison,
    multiplicity_note = paste(
      "Sensitivity p-values are descriptive",
      "and are not included in the final",
      "Holm-adjustment families."
    )
  ),
  results_path
)

expected_rows <- c(
  rfs_boundary_sensitivity_model_overview = 4L,
  rfs_boundary_sensitivity_subtype_estimates = 20L,
  rfs_boundary_sensitivity_event_counts = 24L,
  rfs_boundary_sensitivity_pattern_comparison = 10L
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
    paste(
      "At least one boundary-sensitivity table",
      "has an unexpected row count."
    )
  )
}

if (!file.exists(results_path)) {
  stop(
    "The boundary-sensitivity RDS file was not created."
  )
}

cat(
  "\nScript 23 completed successfully.\n"
)

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