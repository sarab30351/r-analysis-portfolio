
# 08_refined_interval_cox_models.R
# Purpose:
# The original Cox model and the first refined version showed problems as documented in each of them (06_cox_models.R) and (07_interval_cox_models.R)
# I therefore decided to refine the intervals:
#   0 to 2 years
#   2 to 5 years
#   5 to 10 years
#   Beyond 10 years
#
# Models are adjusted for age at diagnosis and the Nottingham Prognostic Index using natural cubic splines. Source cohort is handled through stratification. Luminal A is the reference subtype.


# Check packages 

suppressPackageStartupMessages(
  library(survival)
)

required_packages <- c(
  "survival",
  "splines"
)

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
    paste(
      "Install the following package(s) before running this script:",
      paste(missing_packages, collapse = ", ")
    )
  )
}


# Check project location 

if (!file.exists("r-analysis-portfolio.Rproj")) {
  stop(
    paste(
      "Open r-analysis-portfolio.Rproj before running this script.",
      "The working directory must be the repository root."
    )
  )
}

os_cohort_path <- file.path(
  "data-derived",
  "metabric_os_cohort.rds"
)

if (!file.exists(os_cohort_path)) {
  stop(
    paste(
      "Overall-survival cohort not found at:",
      os_cohort_path
    )
  )
}


# Load and validate the overall-survival cohort 

os_cohort <- readRDS(os_cohort_path)

# Convert overall-survival follow-up from months to years 

if (!"os_years" %in% names(os_cohort)) {
  
  if (!"os_months" %in% names(os_cohort)) {
    stop(
      paste(
        "Neither os_years nor os_months was found",
        "in the overall-survival cohort."
      )
    )
  }
  
  os_cohort$os_years <- os_cohort$os_months / 12
}

required_model_variables <- c(
  "os_years",
  "os_event",
  "age_at_diagnosis",
  "npi",
  "cohort",
  "molecular_subtype"
)

missing_variables <- setdiff(
  required_model_variables,
  names(os_cohort)
)

if (length(missing_variables) > 0) {
  stop(
    paste(
      "Required variable(s) missing from the OS cohort:",
      paste(missing_variables, collapse = ", ")
    )
  )
}

complete_model_records <- stats::complete.cases(
  os_cohort[, required_model_variables]
)

if (!all(complete_model_records)) {
  stop(
    paste(
      sum(!complete_model_records),
      "record(s) have missing model variables."
    )
  )
}

os_model_data <- os_cohort

subtype_levels <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

unexpected_subtypes <- setdiff(
  unique(as.character(os_model_data$molecular_subtype)),
  subtype_levels
)

if (length(unexpected_subtypes) > 0) {
  stop(
    paste(
      "Unexpected molecular subtype(s):",
      paste(unexpected_subtypes, collapse = ", ")
    )
  )
}

os_model_data$molecular_subtype <- factor(
  os_model_data$molecular_subtype,
  levels = subtype_levels
)

os_model_data$cohort <- factor(
  os_model_data$cohort
)


# Recreate the spline knots used in Script 06 

age_knots <- as.numeric(
  stats::quantile(
    os_model_data$age_at_diagnosis,
    probs = c(0.05, 0.35, 0.65, 0.95),
    names = FALSE
  )
)

npi_knots <- as.numeric(
  stats::quantile(
    os_model_data$npi,
    probs = c(0.05, 0.50, 0.95),
    names = FALSE
  )
)


# Function to create an interval-specific risk set 
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
  
  if (is.infinite(end_time)) {
    interval_data$exit_time <- interval_data$os_years
  } else {
    interval_data$exit_time <- pmin(
      interval_data$os_years,
      end_time
    )
  }
  
  interval_data$interval_event <- as.integer(
    interval_data$os_event == 1L &
      interval_data$os_years > start_time &
      (
        is.infinite(end_time) |
          interval_data$os_years <= end_time
      )
  )
  
  interval_data$interval_label <- interval_label
  
  if (any(interval_data$exit_time <= interval_data$entry_time)) {
    stop(
      paste(
        "Invalid follow-up time found in interval:",
        interval_label
      )
    )
  }
  
  interval_data
}


# Function to fit and evaluate one interval 

fit_interval_models <- function(
    interval_data,
    interval_label
) {
  
  clinical_model <- survival::coxph(
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
      strata(cohort),
    data = interval_data,
    ties = "efron",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  extended_model <- survival::coxph(
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
    data = interval_data,
    ties = "efron",
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  model_comparison <- stats::anova(
    clinical_model,
    extended_model,
    test = "LRT"
  )
  
  model_comparison_summary <- data.frame(
    interval = interval_label,
    participants = nrow(interval_data),
    deaths = sum(interval_data$interval_event),
    clinical_log_likelihood = as.numeric(
      stats::logLik(clinical_model)
    ),
    extended_log_likelihood = as.numeric(
      stats::logLik(extended_model)
    ),
    likelihood_ratio_chisq = unname(
      model_comparison[2, "Chisq"]
    ),
    degrees_of_freedom = unname(
      model_comparison[2, "Df"]
    ),
    p_value = unname(
      model_comparison[2, "Pr(>|Chi|)"]
    )
  )
  
  coefficient_table <- summary(
    extended_model
  )$coefficients
  
  confidence_intervals <- stats::confint(
    extended_model
  )
  
  subtype_rows <- grep(
    "^molecular_subtype",
    rownames(coefficient_table)
  )
  
  if (length(subtype_rows) != 5) {
    stop(
      paste(
        "Expected five subtype coefficients in interval:",
        interval_label
      )
    )
  }
  
  subtype_coefficient_names <- rownames(
    coefficient_table
  )[subtype_rows]
  
  subtype_estimates <- data.frame(
    interval = interval_label,
    comparison = paste(
      sub(
        "^molecular_subtype",
        "",
        subtype_coefficient_names
      ),
      "vs Luminal A"
    ),
    hazard_ratio = exp(
      coefficient_table[subtype_rows, "coef"]
    ),
    lower_95_ci = exp(
      confidence_intervals[
        subtype_coefficient_names,
        1
      ]
    ),
    upper_95_ci = exp(
      confidence_intervals[
        subtype_coefficient_names,
        2
      ]
    ),
    p_value = coefficient_table[
      subtype_rows,
      "Pr(>|z|)"
    ],
    row.names = NULL
  )
  
  ph_terms <- survival::cox.zph(
    extended_model,
    transform = "km",
    terms = TRUE,
    global = TRUE
  )
  
  ph_term_table <- data.frame(
    term = rownames(ph_terms$table),
    ph_terms$table,
    row.names = NULL,
    check.names = FALSE
  )
  
  ph_coefficients <- survival::cox.zph(
    extended_model,
    transform = "km",
    terms = FALSE,
    global = TRUE
  )
  
  subtype_ph_rows <- grepl(
    "^molecular_subtype",
    rownames(ph_coefficients$table)
  )
  
  subtype_ph_table <- data.frame(
    coefficient = rownames(
      ph_coefficients$table
    )[subtype_ph_rows],
    ph_coefficients$table[
      subtype_ph_rows,
      ,
      drop = FALSE
    ],
    row.names = NULL,
    check.names = FALSE
  )
  
  list(
    interval = interval_label,
    participants = nrow(interval_data),
    deaths = sum(interval_data$interval_event),
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_comparison = model_comparison,
    model_comparison_summary = model_comparison_summary,
    subtype_estimates = subtype_estimates,
    ph_terms = ph_terms,
    ph_term_table = ph_term_table,
    ph_coefficients = ph_coefficients,
    subtype_ph_table = subtype_ph_table
  )
}


# Create the four refined risk sets 

interval_data_sets <- list(
  zero_to_two = create_interval_data(
    data = os_model_data,
    start_time = 0,
    end_time = 2,
    interval_label = "0 to 2 years"
  ),
  
  two_to_five = create_interval_data(
    data = os_model_data,
    start_time = 2,
    end_time = 5,
    interval_label = "2 to 5 years"
  ),
  
  five_to_ten = create_interval_data(
    data = os_model_data,
    start_time = 5,
    end_time = 10,
    interval_label = "5 to 10 years"
  ),
  
  beyond_ten = create_interval_data(
    data = os_model_data,
    start_time = 10,
    end_time = Inf,
    interval_label = "Beyond 10 years"
  )
)


# Fit all four interval-specific models 

refined_interval_results <- Map(
  f = function(interval_data, interval_name) {
    fit_interval_models(
      interval_data = interval_data,
      interval_label = interval_name
    )
  },
  interval_data = interval_data_sets,
  interval_name = vapply(
    interval_data_sets,
    function(data) {
      unique(data$interval_label)
    },
    FUN.VALUE = character(1)
  )
)


# Combine the main results 

refined_model_comparisons <- do.call(
  rbind,
  lapply(
    refined_interval_results,
    function(result) {
      result$model_comparison_summary
    }
  )
)

refined_subtype_estimates <- do.call(
  rbind,
  lapply(
    refined_interval_results,
    function(result) {
      result$subtype_estimates
    }
  )
)

refined_ph_results <- do.call(
  rbind,
  lapply(
    refined_interval_results,
    function(result) {
      
      table <- result$ph_term_table
      table$interval <- result$interval
      
      table[
        ,
        c(
          "interval",
          setdiff(names(table), "interval")
        )
      ]
    }
  )
)


# Print results 
for (result in refined_interval_results) {
  
  cat(
    "\n\n",
    result$interval,
    "\n",
    sep = ""
  )
  
  cat(
    "Participants: ",
    result$participants,
    "; deaths: ",
    result$deaths,
    "\n",
    sep = ""
  )
  
  cat("\nLikelihood-ratio comparison:\n")
  print(result$model_comparison)
  
  cat("\nAdjusted subtype estimates:\n")
  print(
    result$subtype_estimates,
    row.names = FALSE
  )
  
  cat("\nTerm-level proportional-hazards tests:\n")
  print(
    result$ph_term_table,
    row.names = FALSE
  )
  
  cat("\nSubtype coefficient-level PH tests:\n")
  print(
    result$subtype_ph_table,
    row.names = FALSE
  )
}


# Save model-development results locally 

refined_model_path <- file.path(
  "data-derived",
  "os_refined_interval_cox_models.rds"
)

saveRDS(
  list(
    age_knots = age_knots,
    npi_knots = npi_knots,
    interval_results = refined_interval_results,
    model_comparisons = refined_model_comparisons,
    subtype_estimates = refined_subtype_estimates,
    proportional_hazards_results = refined_ph_results
  ),
  refined_model_path
)

cat(
  "\n\nRefined interval models saved locally to: ",
  normalizePath(refined_model_path),
  "\n",
  sep = ""
)

### Refinement worked in the most critical areas i.e the severe violations, however there are still smaller problems to address. Update: I stopped refining the model for now (after sensitivty analysis in script 11).  Refining it further would introduce uneccessary complexity.
### In the first two years, Basal-like disease has an estimated HR of 7.18 and HER2-enriched disease an HR of 3.93 versus Luminal A. Between two and five years, those HRs decline to 2.05 and 2.29. Between five and ten years, there is no clear overall subtype contribution beyond age and NPI. Beyond ten years, the apparent lower Basal-like hazard applies only to ten-year survivors. This late hazard reversal does not represent an overall survival advantage from diagnosis and may partly reflect selection of long-term survivors and the small number of late Basal-like deaths.
