
# 12_os_model_results.R
# Purpose: Export the finalized overall survival Cox model results and create an interval-specific adjusted hazard ratio forest plot.
#
# Primary analysis: Refined interval Cox models from Script 08 (08_refined_interval_cox_models.R).
#
# Sensitivity analyses: Age-stratified models from Script 10 (10_age_stratified_sensitivity.R). NPI-stratified 2-to-5-year model from Script 11 (11_npi_stratified_sensitivity.R).


# Check required packages 

required_packages <- c(
  "ggplot2",
  "scales"
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
      "Install the following package(s) before running Script 12:",
      paste(missing_packages, collapse = ", ")
    )
  )
}


# Define paths 

primary_results_path <- file.path(
  "data-derived",
  "os_refined_interval_cox_models.rds"
)

age_sensitivity_path <- file.path(
  "data-derived",
  "os_age_stratified_sensitivity.rds"
)

npi_sensitivity_path <- file.path(
  "data-derived",
  "os_npi_stratified_sensitivity.rds"
)

required_files <- c(
  primary_results_path,
  age_sensitivity_path,
  npi_sensitivity_path
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Required result file(s) missing:",
      paste(missing_files, collapse = ", ")
    )
  )
}

table_directory <- file.path(
  "output",
  "tables"
)

figure_directory <- file.path(
  "output",
  "figures"
)

dir.create(
  table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# Read results 

primary_results <- readRDS(
  primary_results_path
)

age_sensitivity_results <- readRDS(
  age_sensitivity_path
)

npi_sensitivity_results <- readRDS(
  npi_sensitivity_path
)


# Check required result objects 

required_primary_objects <- c(
  "model_comparisons",
  "subtype_estimates",
  "proportional_hazards_results"
)

missing_primary_objects <- setdiff(
  required_primary_objects,
  names(primary_results)
)

if (length(missing_primary_objects) > 0) {
  stop(
    paste(
      "The primary results file is missing:",
      paste(missing_primary_objects, collapse = ", ")
    )
  )
}

if (!"hazard_ratio_comparison" %in%
    names(age_sensitivity_results)) {
  stop(
    paste(
      "The age-stratified sensitivity file does not contain",
      "`hazard_ratio_comparison`."
    )
  )
}

if (!"hazard_ratio_comparison" %in%
    names(npi_sensitivity_results)) {
  stop(
    paste(
      "The NPI-stratified sensitivity file does not contain",
      "`hazard_ratio_comparison`."
    )
  )
}


# Helper functions 

format_p_value <- function(p_value) {
  ifelse(
    is.na(p_value),
    NA_character_,
    ifelse(
      p_value < 0.001,
      "<0.001",
      formatC(
        p_value,
        format = "f",
        digits = 3
      )
    )
  )
}

check_columns <- function(
    data,
    required_columns,
    object_name
) {
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      paste(
        object_name,
        "is missing the following column(s):",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
}


# Extract primary results 

primary_model_comparisons <- primary_results$model_comparisons

primary_subtype_estimates <- primary_results$subtype_estimates

primary_ph_results <- primary_results$proportional_hazards_results


check_columns(
  primary_model_comparisons,
  c(
    "interval",
    "participants",
    "deaths",
    "likelihood_ratio_chisq",
    "degrees_of_freedom",
    "p_value"
  ),
  "Primary model-comparison table"
)

check_columns(
  primary_subtype_estimates,
  c(
    "interval",
    "comparison",
    "hazard_ratio",
    "lower_95_ci",
    "upper_95_ci",
    "p_value"
  ),
  "Primary subtype-estimate table"
)

check_columns(
  primary_ph_results,
  c(
    "interval",
    "term",
    "chisq",
    "df",
    "p"
  ),
  "Primary proportional hazards table"
)

# Adjust primary p-values for multiple testing 

# Four interval-level likelihood-ratio tests
primary_model_comparisons$p_value_holm_4_intervals <-
  stats::p.adjust(
    primary_model_comparisons$p_value,
    method = "holm"
  )

# Twenty subtype-vs-Luminal-A contrasts:
# five comparisons across four intervals
primary_subtype_estimates$p_value_holm_20_contrasts <-
  stats::p.adjust(
    primary_subtype_estimates$p_value,
    method = "holm"
  )


# Define interval order 

interval_order <- c(
  "0 to 2 years",
  "2 to 5 years",
  "5 to 10 years",
  "Beyond 10 years"
)


# Prepare primary model comparison table 

model_comparison_export <- primary_model_comparisons

model_comparison_export$interval <- factor(
  model_comparison_export$interval,
  levels = interval_order
)

model_comparison_export <- model_comparison_export[
  order(model_comparison_export$interval),
  ,
  drop = FALSE
]

model_comparison_export$interval <- as.character(
  model_comparison_export$interval
)

model_comparison_export$p_value_display <- format_p_value(
  model_comparison_export$p_value
)


model_comparison_export$p_value_holm_4_intervals_display <-
  format_p_value(
    model_comparison_export$p_value_holm_4_intervals
  )


model_comparison_export$likelihood_ratio_chisq <- round(
  model_comparison_export$likelihood_ratio_chisq,
  digits = 3
)

model_comparison_export$p_value <- signif(
  model_comparison_export$p_value,
  digits = 6
)


model_comparison_export$p_value_holm_4_intervals <- signif(
  model_comparison_export$p_value_holm_4_intervals,
  digits = 6
)


# Prepare primary subtype estimate table 
subtype_estimate_export <- primary_subtype_estimates

subtype_estimate_export$p_value_display <- format_p_value(
  subtype_estimate_export$p_value
)

subtype_estimate_export$p_value_holm_20_contrasts_display <-
  format_p_value(
    subtype_estimate_export$p_value_holm_20_contrasts
  )


subtype_estimate_export$hazard_ratio <- round(
  subtype_estimate_export$hazard_ratio,
  digits = 3
)

subtype_estimate_export$lower_95_ci <- round(
  subtype_estimate_export$lower_95_ci,
  digits = 3
)

subtype_estimate_export$upper_95_ci <- round(
  subtype_estimate_export$upper_95_ci,
  digits = 3
)

subtype_estimate_export$p_value <- signif(
  subtype_estimate_export$p_value,
  digits = 6
)

subtype_estimate_export$p_value_holm_20_contrasts <- signif(
  subtype_estimate_export$p_value_holm_20_contrasts,
  digits = 6
)

subtype_estimate_export$interval <- factor(
  subtype_estimate_export$interval,
  levels = interval_order
)

subtype_estimate_export <- subtype_estimate_export[
  order(
    subtype_estimate_export$interval,
    subtype_estimate_export$comparison
  ),
  ,
  drop = FALSE
]

subtype_estimate_export$interval <- as.character(
  subtype_estimate_export$interval
)


# Prepare proportional hazards table 
ph_results_export <- primary_ph_results

ph_results_export$p_value_display <- format_p_value(
  ph_results_export$p
)

ph_results_export$chisq <- round(
  ph_results_export$chisq,
  digits = 3
)

ph_results_export$p <- signif(
  ph_results_export$p,
  digits = 6
)

ph_results_export$interval <- factor(
  ph_results_export$interval,
  levels = interval_order
)

ph_results_export <- ph_results_export[
  order(ph_results_export$interval),
  ,
  drop = FALSE
]

ph_results_export$interval <- as.character(
  ph_results_export$interval
)


# Prepare age stratified sensitivity table 

age_sensitivity_export <-
  age_sensitivity_results$hazard_ratio_comparison

check_columns(
  age_sensitivity_export,
  c(
    "interval",
    "comparison",
    "primary_hazard_ratio",
    "sensitivity_hazard_ratio",
    "sensitivity_lower_95_ci",
    "sensitivity_upper_95_ci",
    "sensitivity_p_value",
    "percent_change"
  ),
  "Age-stratified sensitivity table"
)

age_sensitivity_export$sensitivity_p_value_display <-
  format_p_value(
    age_sensitivity_export$sensitivity_p_value
  )

columns_to_round_age <- c(
  "primary_hazard_ratio",
  "primary_lower_95_ci",
  "primary_upper_95_ci",
  "sensitivity_hazard_ratio",
  "sensitivity_lower_95_ci",
  "sensitivity_upper_95_ci"
)

age_sensitivity_export[
  columns_to_round_age
] <- lapply(
  age_sensitivity_export[columns_to_round_age],
  round,
  digits = 3
)

age_sensitivity_export$percent_change <- round(
  age_sensitivity_export$percent_change,
  digits = 1
)

age_sensitivity_export$sensitivity_p_value <- signif(
  age_sensitivity_export$sensitivity_p_value,
  digits = 6
)


# Prepare NPI stratified sensitivity table 

npi_sensitivity_export <-
  npi_sensitivity_results$hazard_ratio_comparison

check_columns(
  npi_sensitivity_export,
  c(
    "interval",
    "comparison",
    "primary_hazard_ratio",
    "sensitivity_hazard_ratio",
    "sensitivity_lower_95_ci",
    "sensitivity_upper_95_ci",
    "sensitivity_p_value",
    "percent_change"
  ),
  "NPI-stratified sensitivity table"
)

npi_sensitivity_export$sensitivity_p_value_display <-
  format_p_value(
    npi_sensitivity_export$sensitivity_p_value
  )

columns_to_round_npi <- c(
  "primary_hazard_ratio",
  "primary_lower_95_ci",
  "primary_upper_95_ci",
  "sensitivity_hazard_ratio",
  "sensitivity_lower_95_ci",
  "sensitivity_upper_95_ci"
)

npi_sensitivity_export[
  columns_to_round_npi
] <- lapply(
  npi_sensitivity_export[columns_to_round_npi],
  round,
  digits = 3
)

npi_sensitivity_export$percent_change <- round(
  npi_sensitivity_export$percent_change,
  digits = 1
)

npi_sensitivity_export$sensitivity_p_value <- signif(
  npi_sensitivity_export$sensitivity_p_value,
  digits = 6
)


# Export tables 

model_comparison_output_path <- file.path(
  table_directory,
  "os_interval_model_comparisons.csv"
)

subtype_estimate_output_path <- file.path(
  table_directory,
  "os_interval_subtype_hazard_ratios.csv"
)

ph_results_output_path <- file.path(
  table_directory,
  "os_interval_ph_tests.csv"
)

age_sensitivity_output_path <- file.path(
  table_directory,
  "os_age_stratified_sensitivity.csv"
)

npi_sensitivity_output_path <- file.path(
  table_directory,
  "os_npi_stratified_sensitivity.csv"
)

utils::write.csv(
  model_comparison_export,
  model_comparison_output_path,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  subtype_estimate_export,
  subtype_estimate_output_path,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  ph_results_export,
  ph_results_output_path,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  age_sensitivity_export,
  age_sensitivity_output_path,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  npi_sensitivity_export,
  npi_sensitivity_output_path,
  row.names = FALSE,
  na = ""
)


# Prepare forest-plot data 

plot_data <- primary_subtype_estimates

plot_data$subtype <- sub(
  " vs Luminal A$",
  "",
  plot_data$comparison
)

subtype_order <- c(
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

plot_data$subtype <- factor(
  plot_data$subtype,
  levels = rev(subtype_order)
)

interval_p_values <- setNames(
  primary_model_comparisons$p_value_holm_4_intervals,
  primary_model_comparisons$interval
)

interval_facet_labels <- paste0(
  interval_order,
  "\nHolm-adjusted likelihood-ratio p ",
  format_p_value(
    interval_p_values[interval_order]
  )
)

names(interval_facet_labels) <- interval_order

plot_data$facet_label <- unname(
  interval_facet_labels[plot_data$interval]
)

plot_data$facet_label <- factor(
  plot_data$facet_label,
  levels = interval_facet_labels
)

subtype_colors <- c(
  "Luminal B" = "#E69F00",
  "HER2-enriched" = "#D55E00",
  "Basal-like" = "#CC79A7",
  "Normal-like" = "#009E73",
  "Claudin-low" = "#56B4E9"
)


# Create adjusted hazard-ratio forest plot 

adjusted_hr_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = hazard_ratio,
    y = subtype,
    color = subtype
  )
) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "gray45"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = lower_95_ci,
      xend = upper_95_ci,
      yend = subtype
    ),
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::facet_wrap(
    ~ facet_label,
    ncol = 2
  ) +
  ggplot2::scale_color_manual(
    values = subtype_colors,
    guide = "none"
  ) +
  ggplot2::scale_x_log10(
    limits = c(0.2, 16),
    breaks = c(
      0.25,
      0.5,
      1,
      2,
      4,
      8,
      16
    ),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  ggplot2::labs(
    title = paste(
      "Adjusted mortality hazard ratios",
      "by molecular subtype"
    ),
    subtitle = paste(
      "Reference: Luminal A;",
      "models adjusted for age and NPI",
      "and stratified by METABRIC source cohort"
    ),
    x = "Adjusted hazard ratio, logarithmic scale",
    y = NULL,
    caption = paste(
      paste(
        "Points show hazard ratios and horizontal lines show",
        "pointwise 95% confidence intervals."
      ),
      paste(
        "Facet p-values use Holm adjustment across the four",
        "interval-level likelihood-ratio tests."
      ),
      "Estimates apply within each follow-up interval.",
      sep = "\n"
    )
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    axis.text.y = ggplot2::element_text(
      color = "black"
    ),
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 9,
      lineheight = 1.1,
      margin = ggplot2::margin(t = 8)
    )
  )
# Save forest plot 

forest_plot_output_path <- file.path(
  figure_directory,
  "os_adjusted_hazard_ratios_by_interval.png"
)

ggplot2::ggsave(
  filename = forest_plot_output_path,
  plot = adjusted_hr_plot,
  width = 12,
  height = 9,
  units = "in",
  dpi = 320,
  bg = "white"
)


# Print verification output 

print(adjusted_hr_plot)

message("")
message("Script 12 completed successfully.")
message("")
message("Tables created:")

for (output_path in c(
  model_comparison_output_path,
  subtype_estimate_output_path,
  ph_results_output_path,
  age_sensitivity_output_path,
  npi_sensitivity_output_path
)) {
  message(
    "  ",
    normalizePath(output_path)
  )
}

message("")
message("Figure created:")
message(
  "  ",
  normalizePath(forest_plot_output_path)
)

### Main interpretations will be in script 08 and 11 ### 