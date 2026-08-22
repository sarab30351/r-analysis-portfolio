
# 02_prepare_data.R
#
# Purpose:
#   Prepare analysis variables and eligibility indicators for the METABRIC
#   overall-survival and relapse-free survival cohorts.
#
# Important:
#   This script does not exclude patients. Cohorts will be defined separately.



# Locate and import the joined dataset 
input_file <- here::here(
  "data-derived",
  "metabric_clinical_joined.rds"
)

if (!file.exists(input_file)) {
  stop(
    "The joined clinical dataset is missing. ",
    "Run R/01_import_validate_data.R first."
  )
}

clinical_data <- readRDS(input_file) |>
  dplyr::rename_with(tolower)


# Validate the observed categorical values 

check_expected_values <- function(variable, expected_values, variable_name) {
  
  observed_values <- unique(variable[!is.na(variable)])
  unexpected_values <- setdiff(observed_values, expected_values)
  
  if (length(unexpected_values) > 0) {
    stop(
      "Unexpected value(s) in ",
      variable_name,
      ": ",
      paste(unexpected_values, collapse = ", ")
    )
  }
}

check_expected_values(
  clinical_data$os_status,
  c("0:LIVING", "1:DECEASED"),
  "os_status"
)

check_expected_values(
  clinical_data$rfs_status,
  c("0:Not Recurred", "1:Recurred"),
  "rfs_status"
)

check_expected_values(
  clinical_data$claudin_subtype,
  c(
    "LumA",
    "LumB",
    "Her2",
    "Basal",
    "Normal",
    "claudin-low",
    "NC"
  ),
  "claudin_subtype"
)


# Prepare variables and eligibility indicators 

prepared_data <- clinical_data |>
  dplyr::mutate(
    
    # Cohort is categorical and will later be used for stratification.
    cohort = factor(cohort),
    
    # Ensure that stage is stored numerically.
    tumor_stage = as.numeric(as.character(tumor_stage)),
    
    # Recode overall-survival status.
    os_event = dplyr::case_when(
      os_status == "1:DECEASED" ~ 1L,
      os_status == "0:LIVING"   ~ 0L,
      TRUE                      ~ NA_integer_
    ),
    
    # Recode relapse-free-survival status.
    rfs_event = dplyr::case_when(
      rfs_status == "1:Recurred"     ~ 1L,
      rfs_status == "0:Not Recurred" ~ 0L,
      TRUE                           ~ NA_integer_
    ),
    
    # Recode the PAM50 + Claudin-low classification.
    molecular_subtype = dplyr::case_when(
      claudin_subtype == "LumA"         ~ "Luminal A",
      claudin_subtype == "LumB"         ~ "Luminal B",
      claudin_subtype == "Her2"         ~ "HER2-enriched",
      claudin_subtype == "Basal"        ~ "Basal-like",
      claudin_subtype == "Normal"       ~ "Normal-like",
      claudin_subtype == "claudin-low"  ~ "Claudin-low",
      claudin_subtype == "NC"           ~ NA_character_,
      is.na(claudin_subtype)            ~ NA_character_
    ),
    
    molecular_subtype = factor(
      molecular_subtype,
      levels = c(
        "Luminal A",
        "Luminal B",
        "HER2-enriched",
        "Basal-like",
        "Normal-like",
        "Claudin-low"
      )
    ),
    
    # Identify records with an interpretable molecular classification.
    subtype_interpretable =
      !is.na(claudin_subtype) &
      claudin_subtype != "NC",
    
    # Identify the non-epithelial tumors excluded from the target population.
    non_epithelial_angiosarcoma =
      dplyr::coalesce(
        cancer_type_detailed == "Breast Angiosarcoma",
        FALSE
      ),
    
    # Record stage availability and known stage IV disease.
    stage_missing = is.na(tumor_stage),
    
    stage_iv_recorded =
      dplyr::coalesce(tumor_stage == 4, FALSE),
    
    # Check that survival times are nonmissing and nonnegative.
    os_time_valid =
      !is.na(os_months) &
      os_months >= 0,
    
    rfs_time_valid =
      !is.na(rfs_months) &
      rfs_months >= 0,
    
    # Identify month-zero observations.
    os_zero_censored =
      dplyr::coalesce(
        os_months == 0 & os_event == 0L,
        FALSE
      ),
    
    os_zero_event =
      dplyr::coalesce(
        os_months == 0 & os_event == 1L,
        FALSE
      ),
    
    rfs_zero_censored =
      dplyr::coalesce(
        rfs_months == 0 & rfs_event == 0L,
        FALSE
      ),
    
    rfs_zero_event =
      dplyr::coalesce(
        rfs_months == 0 & rfs_event == 1L,
        FALSE
      ),
    
    # An event at month zero is informative. A censored record at month zero is
    # not informative because it contributes neither follow-up nor an event.
    os_information_available =
      os_time_valid &
      !is.na(os_event) &
      (os_months > 0 | os_event == 1L),
    
    rfs_information_available =
      rfs_time_valid &
      !is.na(rfs_event) &
      (rfs_months > 0 | rfs_event == 1L)
  )


# Validate the recoded variables 

# Validate the recoded variables ------------------------------------------------

os_events_missing_time <- sum(
  prepared_data$os_event == 1L &
    is.na(prepared_data$os_months),
  na.rm = TRUE
)

rfs_events_missing_time <- sum(
  prepared_data$rfs_event == 1L &
    is.na(prepared_data$rfs_months),
  na.rm = TRUE
)

if (any(prepared_data$os_months < 0, na.rm = TRUE)) {
  stop("Negative overall-survival time was detected.")
}

if (any(prepared_data$rfs_months < 0, na.rm = TRUE)) {
  stop("Negative relapse-free-survival time was detected.")
}

# Save the prepared full dataset

output_file <- here::here(
  "data-derived",
  "metabric_clinical_prepared.rds"
)

saveRDS(
  object = prepared_data,
  file = output_file
)


# Report preparation results 

message("Data preparation checks passed.")
message("Prepared records: ", nrow(prepared_data))
message(
  "Interpretable molecular subtypes: ",
  sum(prepared_data$subtype_interpretable)
)
message(
  "Breast angiosarcomas: ",
  sum(prepared_data$non_epithelial_angiosarcoma)
)
message(
  "Known stage IV records: ",
  sum(prepared_data$stage_iv_recorded)
)
message(
  "Overall-survival month-zero censored records: ",
  sum(prepared_data$os_zero_censored)
)
message(
  "Relapse-free-survival month-zero events: ",
  sum(prepared_data$rfs_zero_event)
)
message(
  "Relapse-free-survival month-zero censored records: ",
  sum(prepared_data$rfs_zero_censored)
)
message("Prepared dataset saved to: ", output_file)

message(
  "Overall-survival events with missing time: ",
  os_events_missing_time
)

message(
  "Relapse-free-survival events with missing time: ",
  rfs_events_missing_time
)
