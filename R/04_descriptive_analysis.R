
# 04_descriptive_analysis.R
#
# Purpose: Describe the primary overall-survival cohort and summarize outcome information for both the OS and RFS cohorts.



# Import analytic cohorts

os_file <- here::here(
  "data-derived",
  "metabric_os_cohort.rds"
)

rfs_file <- here::here(
  "data-derived",
  "metabric_rfs_cohort.rds"
)

if (!file.exists(os_file) || !file.exists(rfs_file)) {
  stop(
    "One or both analytic cohort files are missing. ",
    "Run R/03_define_cohorts.R first."
  )
}

os_cohort <- readRDS(os_file) |>
  dplyr::mutate(
    os_years = os_months / 12
  )

rfs_cohort <- readRDS(rfs_file) |>
  dplyr::mutate(
    rfs_years = rfs_months / 12
  )


# Variable labels

variable_labels <- c(
  age_at_diagnosis = "Age at diagnosis (years)",
  npi = "Nottingham Prognostic Index",
  tumor_size = "Tumor size",
  lymph_nodes_examined_positive = "Positive lymph nodes",
  tumor_stage = "Tumor stage",
  grade = "Histologic grade",
  er_status = "ER status",
  pr_status = "PR status",
  her2_status = "HER2 status",
  chemotherapy = "Chemotherapy",
  hormone_therapy = "Hormone therapy",
  radio_therapy = "Radiation therapy",
  breast_surgery = "Breast surgery",
  histological_subtype = "Histologic subtype",
  cohort = "Component cohort"
)


# Continuous-variable summary function 

summarize_continuous_variable <- function(
    data,
    variable,
    variable_label
) {
  
  analysis_groups <- c(
    "Overall",
    levels(data$molecular_subtype)
  )
  
  group_summaries <- lapply(
    analysis_groups,
    function(group_name) {
      
      group_data <- if (group_name == "Overall") {
        data
      } else {
        data[
          data$molecular_subtype == group_name,
          ,
          drop = FALSE
        ]
      }
      
      values <- group_data[[variable]]
      observed_values <- values[!is.na(values)]
      
      data.frame(
        variable = variable,
        variable_label = variable_label,
        group = group_name,
        n_total = nrow(group_data),
        n_observed = length(observed_values),
        n_missing = sum(is.na(values)),
        mean = if (length(observed_values) > 0) {
          mean(observed_values)
        } else {
          NA_real_
        },
        standard_deviation = if (length(observed_values) > 1) {
          stats::sd(observed_values)
        } else {
          NA_real_
        },
        median = if (length(observed_values) > 0) {
          stats::median(observed_values)
        } else {
          NA_real_
        },
        first_quartile = if (length(observed_values) > 0) {
          unname(stats::quantile(observed_values, 0.25))
        } else {
          NA_real_
        },
        third_quartile = if (length(observed_values) > 0) {
          unname(stats::quantile(observed_values, 0.75))
        } else {
          NA_real_
        },
        minimum = if (length(observed_values) > 0) {
          min(observed_values)
        } else {
          NA_real_
        },
        maximum = if (length(observed_values) > 0) {
          max(observed_values)
        } else {
          NA_real_
        }
      )
    }
  )
  
  dplyr::bind_rows(group_summaries)
}


# Categorical variable summary function 
summarize_categorical_variable <- function(
    data,
    variable,
    variable_label
) {
  
  analysis_groups <- c(
    "Overall",
    levels(data$molecular_subtype)
  )
  
  group_summaries <- lapply(
    analysis_groups,
    function(group_name) {
      
      group_data <- if (group_name == "Overall") {
        data
      } else {
        data[
          data$molecular_subtype == group_name,
          ,
          drop = FALSE
        ]
      }
      
      values <- as.character(group_data[[variable]])
      
      values[
        is.na(values) |
          trimws(values) == ""
      ] <- "Missing"
      
      counts <- as.data.frame(
        table(values),
        stringsAsFactors = FALSE
      )
      
      data.frame(
        variable = variable,
        variable_label = variable_label,
        group = group_name,
        level = counts$values,
        n = counts$Freq,
        denominator = nrow(group_data),
        percent = 100 * counts$Freq / nrow(group_data)
      )
    }
  )
  
  dplyr::bind_rows(group_summaries)
}


# Summarize continuous characteristics 

continuous_variables <- c(
  "age_at_diagnosis",
  "npi",
  "tumor_size",
  "lymph_nodes_examined_positive"
)

continuous_summaries <- lapply(
  continuous_variables,
  function(variable) {
    summarize_continuous_variable(
      data = os_cohort,
      variable = variable,
      variable_label = unname(variable_labels[[variable]])
    )
  }
) |>
  dplyr::bind_rows()


# Summarize categorical characteristics

categorical_variables <- c(
  "tumor_stage",
  "grade",
  "er_status",
  "pr_status",
  "her2_status",
  "chemotherapy",
  "hormone_therapy",
  "radio_therapy",
  "breast_surgery",
  "histological_subtype",
  "cohort"
)

categorical_summaries <- lapply(
  categorical_variables,
  function(variable) {
    summarize_categorical_variable(
      data = os_cohort,
      variable = variable,
      variable_label = unname(variable_labels[[variable]])
    )
  }
) |>
  dplyr::bind_rows()


# Summarize missingness in the primary OS cohort 

missingness_variables <- c(
  continuous_variables,
  categorical_variables
)

count_missing <- function(variable) {
  
  values <- os_cohort[[variable]]
  
  if (is.character(values)) {
    sum(is.na(values) | trimws(values) == "")
  } else {
    sum(is.na(values))
  }
}

missingness_summary <- data.frame(
  variable = missingness_variables,
  variable_label = unname(
    variable_labels[missingness_variables]
  ),
  n_total = nrow(os_cohort),
  n_missing = vapply(
    missingness_variables,
    count_missing,
    numeric(1)
  )
) |>
  dplyr::mutate(
    percent_missing = 100 * n_missing / n_total
  )


# Summarize participants and events by subtype 

os_subtype_summary <- os_cohort |>
  dplyr::group_by(
    molecular_subtype,
    .drop = FALSE
  ) |>
  dplyr::summarise(
    participants = dplyr::n(),
    events = sum(os_event),
    censored = sum(os_event == 0L),
    event_percent = 100 * events / participants,
    .groups = "drop"
  ) |>
  dplyr::transmute(
    endpoint = "Overall survival",
    molecular_subtype = as.character(molecular_subtype),
    participants,
    events,
    censored,
    event_percent
  )

rfs_subtype_summary <- rfs_cohort |>
  dplyr::group_by(
    molecular_subtype,
    .drop = FALSE
  ) |>
  dplyr::summarise(
    participants = dplyr::n(),
    events = sum(rfs_event),
    censored = sum(rfs_event == 0L),
    event_percent = 100 * events / participants,
    .groups = "drop"
  ) |>
  dplyr::transmute(
    endpoint = "Relapse-free survival",
    molecular_subtype = as.character(molecular_subtype),
    participants,
    events,
    censored,
    event_percent
  )

subtype_outcomes <- dplyr::bind_rows(
  os_subtype_summary,
  rfs_subtype_summary
)


# Estimate follow-up using reverse Kaplan-Meier 

summarize_reverse_km <- function(
    data,
    time_variable,
    event_variable,
    endpoint_name
) {
  
  analysis_groups <- c(
    "Overall",
    levels(data$molecular_subtype)
  )
  
  followup_summaries <- lapply(
    analysis_groups,
    function(group_name) {
      
      group_data <- if (group_name == "Overall") {
        data
      } else {
        data[
          data$molecular_subtype == group_name,
          ,
          drop = FALSE
        ]
      }
      
      # Reverse Kaplan-Meier:
      # censored observations become events for the follow-up distribution while outcome events are treated as censored.
      followup_object <- survival::Surv(
        time = group_data[[time_variable]],
        event = 1L - group_data[[event_variable]]
      )
      
      followup_fit <- survival::survfit(
        followup_object ~ 1
      )
      
      fit_table <- summary(followup_fit)$table
      
      data.frame(
        endpoint = endpoint_name,
        group = group_name,
        participants = nrow(group_data),
        total_observed_person_years =
          sum(group_data[[time_variable]]),
        median_followup_years =
          unname(fit_table["median"]),
        lower_95_ci =
          unname(fit_table["0.95LCL"]),
        upper_95_ci =
          unname(fit_table["0.95UCL"])
      )
    }
  )
  
  dplyr::bind_rows(followup_summaries)
}

os_followup <- summarize_reverse_km(
  data = os_cohort,
  time_variable = "os_years",
  event_variable = "os_event",
  endpoint_name = "Overall survival"
)

rfs_followup <- summarize_reverse_km(
  data = rfs_cohort,
  time_variable = "rfs_years",
  event_variable = "rfs_event",
  endpoint_name = "Relapse-free survival"
)

followup_summary <- dplyr::bind_rows(
  os_followup,
  rfs_followup
)


# Save aggregate descriptive tables 

output_directory <- here::here(
  "output",
  "tables"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  continuous_summaries,
  here::here(
    "output",
    "tables",
    "descriptive_continuous_os.csv"
  )
)

readr::write_csv(
  categorical_summaries,
  here::here(
    "output",
    "tables",
    "descriptive_categorical_os.csv"
  )
)

readr::write_csv(
  missingness_summary,
  here::here(
    "output",
    "tables",
    "missingness_os.csv"
  )
)

readr::write_csv(
  subtype_outcomes,
  here::here(
    "output",
    "tables",
    "subtype_outcomes.csv"
  )
)

readr::write_csv(
  followup_summary,
  here::here(
    "output",
    "tables",
    "followup_reverse_km.csv"
  )
)


# Report completion 

message("Descriptive-analysis checks completed.")
message("OS cohort participants: ", nrow(os_cohort))
message("RFS cohort participants: ", nrow(rfs_cohort))
message("Five aggregate descriptive tables were saved.")