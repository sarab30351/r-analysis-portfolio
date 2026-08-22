
# 03_define_cohorts.R
#
# Purpose: Apply the eligibility criteria for the overall-survival and relapse-free survival analyses.
#
# Each excluded patient is assigned one primary exclusion reason so that the participant-flow counts are mutually exclusive.


# Import prepared data 

input_file <- here::here(
  "data-derived",
  "metabric_clinical_prepared.rds"
)

if (!file.exists(input_file)) {
  stop(
    "The prepared clinical dataset is missing. ",
    "Run R/02_prepare_data.R first."
  )
}

prepared_data <- readRDS(input_file)


# Define exclusion-reason order 

os_reason_levels <- c(
  "Missing molecular subtype",
  "Unclassified molecular subtype",
  "Breast angiosarcoma",
  "Missing or invalid overall-survival time",
  "Missing overall-survival status",
  "Censored at month zero",
  "Included"
)

rfs_reason_levels <- c(
  "Missing molecular subtype",
  "Unclassified molecular subtype",
  "Breast angiosarcoma",
  "Missing or invalid relapse-free-survival time",
  "Missing relapse-free-survival status",
  "Censored at month zero",
  "Known stage IV disease at baseline",
  "Included"
)


# Assign one primary exclusion reason per patient 

cohort_data <- prepared_data |>
  dplyr::mutate(
    
    os_exclusion_reason = dplyr::case_when(
      is.na(claudin_subtype) ~
        "Missing molecular subtype",
      
      claudin_subtype == "NC" ~
        "Unclassified molecular subtype",
      
      non_epithelial_angiosarcoma ~
        "Breast angiosarcoma",
      
      !os_time_valid ~
        "Missing or invalid overall-survival time",
      
      is.na(os_event) ~
        "Missing overall-survival status",
      
      os_zero_censored ~
        "Censored at month zero",
      
      TRUE ~
        "Included"
    ),
    
    rfs_exclusion_reason = dplyr::case_when(
      is.na(claudin_subtype) ~
        "Missing molecular subtype",
      
      claudin_subtype == "NC" ~
        "Unclassified molecular subtype",
      
      non_epithelial_angiosarcoma ~
        "Breast angiosarcoma",
      
      !rfs_time_valid ~
        "Missing or invalid relapse-free-survival time",
      
      is.na(rfs_event) ~
        "Missing relapse-free-survival status",
      
      rfs_zero_censored ~
        "Censored at month zero",
      
      stage_iv_recorded ~
        "Known stage IV disease at baseline",
      
      TRUE ~
        "Included"
    ),
    
    os_exclusion_reason = factor(
      os_exclusion_reason,
      levels = os_reason_levels
    ),
    
    rfs_exclusion_reason = factor(
      rfs_exclusion_reason,
      levels = rfs_reason_levels
    )
  )


# Create the final analytic cohorts 

os_cohort <- cohort_data |>
  dplyr::filter(os_exclusion_reason == "Included") |>
  droplevels()

rfs_cohort <- cohort_data |>
  dplyr::filter(rfs_exclusion_reason == "Included") |>
  droplevels()


# Validate final cohort counts 

check_count <- function(actual, expected, description) {
  
  if (actual != expected) {
    stop(
      description,
      ": observed ",
      actual,
      ", expected ",
      expected,
      "."
    )
  }
}

if (anyNA(os_cohort$os_event)) {
  stop("The OS cohort contains missing event indicators.")
}

if (anyNA(rfs_cohort$rfs_event)) {
  stop("The RFS cohort contains missing event indicators.")
}

if (
  anyNA(os_cohort$age_at_diagnosis) ||
  anyNA(os_cohort$npi)
) {
  stop("Age or NPI is missing in the OS cohort.")
}

if (
  anyNA(rfs_cohort$age_at_diagnosis) ||
  anyNA(rfs_cohort$npi)
) {
  stop("Age or NPI is missing in the RFS cohort.")
}

if (any(rfs_cohort$stage_iv_recorded)) {
  stop("A known stage IV patient remains in the RFS cohort.")
}

check_count(
  nrow(os_cohort),
  1971L,
  "Overall-survival cohort size"
)

check_count(
  sum(os_cohort$os_event),
  1138L,
  "Overall-survival event count"
)

check_count(
  sum(os_cohort$os_event == 0L),
  833L,
  "Overall-survival censored count"
)

check_count(
  nrow(rfs_cohort),
  1960L,
  "Relapse-free-survival cohort size"
)

check_count(
  sum(rfs_cohort$rfs_event),
  790L,
  "Relapse-free-survival event count"
)

check_count(
  sum(rfs_cohort$rfs_event == 0L),
  1170L,
  "Relapse-free-survival censored count"
)

check_count(
  sum(rfs_cohort$rfs_zero_event),
  3L,
  "Relapse-free-survival month-zero event count"
)


# Create participant-flow counts 

os_flow <- cohort_data |>
  dplyr::count(
    os_exclusion_reason,
    .drop = FALSE,
    name = "n"
  ) |>
  dplyr::filter(n > 0) |>
  dplyr::transmute(
    endpoint = "Overall survival",
    reason = as.character(os_exclusion_reason),
    n = n
  )

rfs_flow <- cohort_data |>
  dplyr::count(
    rfs_exclusion_reason,
    .drop = FALSE,
    name = "n"
  ) |>
  dplyr::filter(n > 0) |>
  dplyr::transmute(
    endpoint = "Relapse-free survival",
    reason = as.character(rfs_exclusion_reason),
    n = n
  )

participant_flow <- dplyr::bind_rows(
  os_flow,
  rfs_flow
)

cohort_summary <- data.frame(
  endpoint = c(
    "Overall survival",
    "Relapse-free survival"
  ),
  participants = c(
    nrow(os_cohort),
    nrow(rfs_cohort)
  ),
  events = c(
    sum(os_cohort$os_event),
    sum(rfs_cohort$rfs_event)
  ),
  censored = c(
    sum(os_cohort$os_event == 0L),
    sum(rfs_cohort$rfs_event == 0L)
  )
)


# Save local analytic datasets 

saveRDS(
  os_cohort,
  here::here(
    "data-derived",
    "metabric_os_cohort.rds"
  )
)

saveRDS(
  rfs_cohort,
  here::here(
    "data-derived",
    "metabric_rfs_cohort.rds"
  )
)


# Save non-identifying aggregate outputs 

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
  participant_flow,
  here::here(
    "output",
    "tables",
    "participant_flow.csv"
  )
)

readr::write_csv(
  cohort_summary,
  here::here(
    "output",
    "tables",
    "cohort_summary.csv"
  )
)

# Report cohort results

message("Cohort-definition checks passed.")

message(
  "Overall-survival cohort: ",
  nrow(os_cohort),
  " patients; ",
  sum(os_cohort$os_event),
  " deaths; ",
  sum(os_cohort$os_event == 0L),
  " censored."
)

message(
  "Relapse-free-survival cohort: ",
  nrow(rfs_cohort),
  " patients; ",
  sum(rfs_cohort$rfs_event),
  " events; ",
  sum(rfs_cohort$rfs_event == 0L),
  " censored."
)

message(
  "RFS patients retained with missing stage: ",
  sum(rfs_cohort$stage_missing)
)

message(
  "RFS month-zero events retained: ",
  sum(rfs_cohort$rfs_zero_event)
)

