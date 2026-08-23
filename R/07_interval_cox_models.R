
# 07_interval_cox_models.R
# Purpose: Fit interval-specific Cox models after the full follow-up model demonstrated violations of the proportional-hazards assumption (documented in 06_cox_models.R)


# Load packages

suppressPackageStartupMessages({
  library(survival)
  library(splines)
})


# Import data and previously defined spline knots

os_file <- here::here(
  "data-derived",
  "metabric_os_cohort.rds"
)

model_file <- here::here(
  "data-derived",
  "os_cox_models.rds"
)

if (!file.exists(os_file)) {
  stop("Run R/03_define_cohorts.R first.")
}

if (!file.exists(model_file)) {
  stop("Run R/06_cox_models.R first.")
}

model_bundle <- readRDS(model_file)

age_knots <- model_bundle$age_knots
npi_knots <- model_bundle$npi_knots

subtype_levels <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

os_model_data <- readRDS(os_file) |>
  dplyr::mutate(
    os_years = os_months / 12,
    molecular_subtype = factor(
      molecular_subtype,
      levels = subtype_levels
    ),
    cohort = factor(cohort)
  )


# Construct interval-specific risk sets

interval_0_to_5 <- os_model_data |>
  dplyr::mutate(
    entry_time = 0,
    exit_time = pmin(os_years, 5),
    interval_event = as.integer(
      os_event == 1L & os_years <= 5
    )
  )

interval_5_to_10 <- os_model_data |>
  dplyr::filter(
    os_years > 5
  ) |>
  dplyr::mutate(
    entry_time = 5,
    exit_time = pmin(os_years, 10),
    interval_event = as.integer(
      os_event == 1L & os_years <= 10
    )
  )

interval_after_10 <- os_model_data |>
  dplyr::filter(
    os_years > 10
  ) |>
  dplyr::mutate(
    entry_time = 10,
    exit_time = os_years,
    interval_event = as.integer(
      os_event == 1L
    )
  )


# Function for fitting one interval

fit_interval_models <- function(interval_data, interval_label) {
  
  if (any(interval_data$exit_time <= interval_data$entry_time)) {
    stop(
      "Invalid follow-up times in interval: ",
      interval_label
    )
  }
  
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
    x = TRUE,
    y = TRUE,
    model = TRUE
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
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  
  model_comparison <- stats::anova(
    clinical_model,
    extended_model,
    test = "LRT"
  )
  
  extended_summary <- summary(extended_model)
  
  subtype_rows <- grepl(
    "^molecular_subtype",
    rownames(extended_summary$coefficients)
  )
  
  subtype_coefficients <- extended_summary$coefficients[
    subtype_rows,
    ,
    drop = FALSE
  ]
  
  subtype_confidence_intervals <- extended_summary$conf.int[
    subtype_rows,
    ,
    drop = FALSE
  ]
  
  subtype_results <- data.frame(
    interval = interval_label,
    comparison = sub(
      "^molecular_subtype",
      "",
      rownames(subtype_coefficients)
    ),
    hazard_ratio = subtype_confidence_intervals[
      ,
      "exp(coef)"
    ],
    lower_95_ci = subtype_confidence_intervals[
      ,
      "lower .95"
    ],
    upper_95_ci = subtype_confidence_intervals[
      ,
      "upper .95"
    ],
    p_value = subtype_coefficients[
      ,
      "Pr(>|z|)"
    ],
    row.names = NULL
  )
  
  ph_test <- survival::cox.zph(
    extended_model,
    transform = "km",
    terms = TRUE,
    global = TRUE
  )
  
  list(
    interval = interval_label,
    participants = nrow(interval_data),
    events = sum(interval_data$interval_event),
    clinical_model = clinical_model,
    extended_model = extended_model,
    model_comparison = model_comparison,
    subtype_results = subtype_results,
    ph_test = ph_test
  )
}


# Fit the three intervals

interval_results <- list(
  "0 to 5 years" = fit_interval_models(
    interval_0_to_5,
    "0 to 5 years"
  ),
  "5 to 10 years" = fit_interval_models(
    interval_5_to_10,
    "5 to 10 years"
  ),
  "Beyond 10 years" = fit_interval_models(
    interval_after_10,
    "Beyond 10 years"
  )
)


# Save models for subsequent reporting

saveRDS(
  interval_results,
  here::here(
    "data-derived",
    "os_interval_cox_models.rds"
  )
)


# Display model comparisons, subtype estimates, and diagnostics

for (result_name in names(interval_results)) {
  
  current_result <- interval_results[[result_name]]
  
  cat(
    "\n\n",
    result_name,
    "\nParticipants: ",
    current_result$participants,
    "; deaths: ",
    current_result$events,
    "\n",
    sep = ""
  )
  
  cat("\nLikelihood-ratio comparison:\n")
  print(current_result$model_comparison)
  
  cat("\nAdjusted subtype estimates:\n")
  print(current_result$subtype_results)
  
  cat("\nProportional-hazards tests:\n")
  print(current_result$ph_test)
}

### Warning: During 0–5 years, Basal-like and HER2-enriched tumors have particularly high early hazards, but their effects are not constant throughout those five years. Therefore, HR 3.02 and HR 2.61 are still averages across a period where the effects change.
### During 5–10 years, the overall subtype test is not significant. Although Luminal B has p = 0.014 individually, we should not emphasize that result when the  global test is p = 0.105.
### Beyond 10 years, Basal-like has HR 0.42. This does not mean Basal-like disease is generally protective though. It applies only to the selected group who survived and remained observed beyond ten years. There were also only 14 Basal-like deaths in this period.

### For refinement, I will split the early period into 0–2 and 2–5 years.
### I will keep 5–10 and beyond 10 years.
### I will address age separately because its effect changes over time.
### I will continue testing NPI and molecular subtype within every interval.
### see documented (08_refined_interval_cox_models.R)