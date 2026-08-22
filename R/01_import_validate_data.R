
# 01_import_validate_data.R
#
# Purpose:
#   Import the METABRIC patient and sample clinical files, verify their basic
#   structure, join them one-to-one by PATIENT_ID, and save the joined dataset.


# Define source-file paths 

patient_file <- here::here(
  "data-raw",
  "data_clinical_patient.txt"
)

sample_file <- here::here(
  "data-raw",
  "data_clinical_sample.txt"
)


# Confirm that the required files exist 
required_files <- c(patient_file, sample_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}


# Import the clinical files 
# Lines beginning with # contain cBioPortal metadata rather than observations.
patient_data <- readr::read_tsv(
  file = patient_file,
  comment = "#",
  na = c("", "NA"),
  trim_ws = TRUE,
  show_col_types = FALSE
)

sample_data <- readr::read_tsv(
  file = sample_file,
  comment = "#",
  na = c("", "NA"),
  trim_ws = TRUE,
  show_col_types = FALSE
)


# Check required variables 
required_patient_variables <- c(
  "PATIENT_ID",
  "AGE_AT_DIAGNOSIS",
  "NPI",
  "OS_MONTHS",
  "OS_STATUS",
  "CLAUDIN_SUBTYPE",
  "RFS_MONTHS",
  "RFS_STATUS"
)

required_sample_variables <- c(
  "PATIENT_ID",
  "SAMPLE_ID",
  "CANCER_TYPE_DETAILED",
  "SAMPLE_TYPE",
  "TUMOR_STAGE"
)

missing_patient_variables <- setdiff(
  required_patient_variables,
  names(patient_data)
)

missing_sample_variables <- setdiff(
  required_sample_variables,
  names(sample_data)
)

if (length(missing_patient_variables) > 0) {
  stop(
    "Missing required patient variables: ",
    paste(missing_patient_variables, collapse = ", ")
  )
}

if (length(missing_sample_variables) > 0) {
  stop(
    "Missing required sample variables: ",
    paste(missing_sample_variables, collapse = ", ")
  )
}


# Validate identifiers and dataset structure 

expected_patient_count <- 2509L

if (nrow(patient_data) != expected_patient_count) {
  stop(
    "Unexpected number of patient records: ",
    nrow(patient_data),
    ". Expected: ",
    expected_patient_count
  )
}

if (nrow(sample_data) != expected_patient_count) {
  stop(
    "Unexpected number of sample records: ",
    nrow(sample_data),
    ". Expected: ",
    expected_patient_count
  )
}

if (anyDuplicated(patient_data$PATIENT_ID) > 0) {
  stop("Duplicated PATIENT_ID values were found in the patient file.")
}

if (anyDuplicated(sample_data$PATIENT_ID) > 0) {
  stop("More than one sample was found for at least one patient.")
}

if (anyDuplicated(sample_data$SAMPLE_ID) > 0) {
  stop("Duplicated SAMPLE_ID values were found in the sample file.")
}

if (!setequal(patient_data$PATIENT_ID, sample_data$PATIENT_ID)) {
  stop("The patient identifiers do not match between the two files.")
}

if (
  anyNA(sample_data$SAMPLE_TYPE) ||
  !all(sample_data$SAMPLE_TYPE == "Primary")
) {
  stop("Not all samples are recorded as primary tumors.")
}

if (
  anyNA(patient_data$SEX) ||
  !all(patient_data$SEX == "Female")
) {
  stop("Not all patients are recorded as female.")
}


# Join patient-level and sample-level information 
clinical_data <- patient_data |>
  dplyr::inner_join(
    sample_data,
    by = "PATIENT_ID",
    relationship = "one-to-one"
  ) |>
  dplyr::arrange(PATIENT_ID)

if (nrow(clinical_data) != expected_patient_count) {
  stop("The joined dataset does not contain the expected number of patients.")
}


# Save the joined dataset locally 

derived_file <- here::here(
  "data-derived",
  "metabric_clinical_joined.rds"
)

saveRDS(
  object = clinical_data,
  file = derived_file
)


# Report completed checks 

message("Import and validation checks passed.")
message("Patient records: ", nrow(patient_data))
message("Sample records: ", nrow(sample_data))
message("Joined records: ", nrow(clinical_data))
message("Joined variables: ", ncol(clinical_data))
message("Derived dataset saved to: ", derived_file)