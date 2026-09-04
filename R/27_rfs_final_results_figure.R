# 27_rfs_final_results_figure.R

# Purpose: Create the final interval-specific adjusted hazard-ratio figure for relapse-free survival using the five-interval results assembled in Script 22.

# This script does not refit models, change interval boundaries, repeat diagnostics or perform additional hypothesis tests. It presents the final estimates and multiplicity-adjusted results already saved by Script 22.


# 1. Check packages and project location

required_packages <- c(
  "ggplot2",
  "scales"
)

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0L) {
  stop(
    "Install before continuing: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

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
    "Run this script from the project root: ",
    project_root
  )
}


# 2. Load the final RFS results

results_path <- file.path(
  project_root,
  "data-derived",
  "rfs_final_interval_analysis.rds"
)

if (!file.exists(results_path)) {
  stop(
    "Final RFS results not found: ",
    results_path,
    ". Run Script 22 first."
  )
}

final_results <- readRDS(
  results_path
)

required_components <- c(
  "interval_definitions",
  "model_overview",
  "subtype_hazard_ratios",
  "multiplicity_method",
  "global_test_family_size",
  "subtype_comparison_family_size",
  "boundary_note"
)

missing_components <- setdiff(
  required_components,
  names(final_results)
)

if (length(missing_components) > 0L) {
  stop(
    "The final RFS results are missing: ",
    paste(
      missing_components,
      collapse = ", "
    )
  )
}

interval_definitions <- (
  final_results$interval_definitions
)

model_overview <- (
  final_results$model_overview
)

subtype_estimates <- (
  final_results$subtype_hazard_ratios
)


# 3. Validate the final results

require_columns <- function(
    data,
    columns,
    label
) {
  missing_columns <- setdiff(
    columns,
    names(data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      label,
      " is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
}

interval_columns <- c(
  "interval_key",
  "interval",
  "start_year",
  "end_year"
)

overview_columns <- c(
  "interval",
  "participants_at_risk",
  "rfs_events",
  "likelihood_ratio_chisq",
  "degrees_of_freedom",
  "p_value",
  "holm_adjusted_p_value",
  "holm_significant"
)

estimate_columns <- c(
  "interval",
  "comparison",
  "hazard_ratio",
  "lower_95_ci",
  "upper_95_ci",
  "p_value",
  "holm_adjusted_p_value",
  "holm_significant"
)

require_columns(
  interval_definitions,
  interval_columns,
  "Final interval definitions"
)

require_columns(
  model_overview,
  overview_columns,
  "Final model overview"
)

require_columns(
  subtype_estimates,
  estimate_columns,
  "Final subtype estimates"
)

interval_order <- c(
  "0 to 2 years",
  "2 to 3.5 years",
  "3.5 to 5 years",
  "5 to 10 years",
  "Beyond 10 years"
)

subtype_order <- c(
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

comparison_order <- paste(
  subtype_order,
  "vs Luminal A"
)

if (
  !identical(
    as.character(
      interval_definitions$interval
    ),
    interval_order
  ) ||
  !setequal(
    as.character(
      model_overview$interval
    ),
    interval_order
  ) ||
  nrow(model_overview) != 5L ||
  nrow(subtype_estimates) != 25L ||
  !identical(
    final_results$multiplicity_method,
    "Holm"
  ) ||
  final_results$global_test_family_size != 5L ||
  final_results$subtype_comparison_family_size != 25L
) {
  stop(
    "The final RFS result structure is unexpected."
  )
}

expected_result_keys <- as.vector(
  outer(
    interval_order,
    comparison_order,
    paste
  )
)

observed_result_keys <- paste(
  subtype_estimates$interval,
  subtype_estimates$comparison
)

if (
  !setequal(
    observed_result_keys,
    expected_result_keys
  ) ||
  anyDuplicated(
    observed_result_keys
  ) > 0L
) {
  stop(
    paste(
      "The final subtype comparisons",
      "are incomplete or duplicated."
    )
  )
}

numeric_estimate_columns <- c(
  "hazard_ratio",
  "lower_95_ci",
  "upper_95_ci",
  "p_value",
  "holm_adjusted_p_value"
)

if (
  any(
    !is.finite(
      as.matrix(
        subtype_estimates[
          numeric_estimate_columns
        ]
      )
    )
  ) ||
  any(
    subtype_estimates$lower_95_ci <= 0
  ) ||
  any(
    subtype_estimates$hazard_ratio <= 0
  ) ||
  any(
    subtype_estimates$upper_95_ci <= 0
  ) ||
  any(
    subtype_estimates$lower_95_ci >
    subtype_estimates$hazard_ratio
  ) ||
  any(
    subtype_estimates$upper_95_ci <
    subtype_estimates$hazard_ratio
  ) ||
  any(
    subtype_estimates$p_value < 0 |
    subtype_estimates$p_value > 1
  ) ||
  any(
    subtype_estimates$holm_adjusted_p_value < 0 |
    subtype_estimates$holm_adjusted_p_value > 1
  ) ||
  any(
    subtype_estimates$holm_adjusted_p_value +
    sqrt(.Machine$double.eps) <
    subtype_estimates$p_value
  ) ||
  !identical(
    as.logical(
      subtype_estimates$holm_significant
    ),
    subtype_estimates$holm_adjusted_p_value < 0.05
  )
) {
  stop(
    "At least one final subtype estimate failed validation."
  )
}

if (
  any(
    !is.finite(
      model_overview$holm_adjusted_p_value
    )
  ) ||
  any(
    model_overview$holm_adjusted_p_value < 0 |
    model_overview$holm_adjusted_p_value > 1
  ) ||
  !identical(
    as.logical(
      model_overview$holm_significant
    ),
    model_overview$holm_adjusted_p_value < 0.05
  ) ||
  sum(
    model_overview$holm_significant
  ) != 4L ||
  sum(
    subtype_estimates$holm_significant
  ) != 4L
) {
  stop(
    "The final Holm-adjusted results failed validation."
  )
}


# 4. Prepare the forest-plot data

format_p_value <- function(
    p_value
) {
  ifelse(
    p_value < 0.001,
    "< 0.001",
    paste0(
      "= ",
      formatC(
        p_value,
        format = "f",
        digits = 3
      )
    )
  )
}

overview_rows <- match(
  interval_order,
  model_overview$interval
)

if (anyNA(overview_rows)) {
  stop(
    "The final interval overview could not be ordered."
  )
}

ordered_overview <- model_overview[
  overview_rows,
  ,
  drop = FALSE
]

interval_facet_labels <- paste0(
  ordered_overview$interval,
  "\nn = ",
  ordered_overview$participants_at_risk,
  "; RFS events = ",
  ordered_overview$rfs_events,
  "\nHolm-adjusted global p ",
  format_p_value(
    ordered_overview$holm_adjusted_p_value
  )
)

names(
  interval_facet_labels
) <- interval_order

plot_data <- subtype_estimates

plot_data$subtype <- sub(
  " vs Luminal A$",
  "",
  plot_data$comparison
)

plot_data$subtype <- factor(
  plot_data$subtype,
  levels = rev(
    subtype_order
  )
)

plot_data$facet_label <- unname(
  interval_facet_labels[
    plot_data$interval
  ]
)

plot_data$facet_label <- factor(
  plot_data$facet_label,
  levels = interval_facet_labels
)

plot_data$comparison_significance <- ifelse(
  plot_data$holm_significant,
  "Holm-adjusted p < 0.05",
  "Not Holm-significant"
)

plot_data$comparison_significance <- factor(
  plot_data$comparison_significance,
  levels = c(
    "Holm-adjusted p < 0.05",
    "Not Holm-significant"
  )
)

if (
  anyNA(
    plot_data$subtype
  ) ||
  anyNA(
    plot_data$facet_label
  ) ||
  anyNA(
    plot_data$comparison_significance
  )
) {
  stop(
    "The forest-plot labels could not be constructed."
  )
}

plot_data <- plot_data[
  order(
    match(
      plot_data$interval,
      interval_order
    ),
    match(
      plot_data$comparison,
      comparison_order
    )
  ),
  ,
  drop = FALSE
]

subtype_colors <- c(
  "Luminal B" = "#E69F00",
  "HER2-enriched" = "#D55E00",
  "Basal-like" = "#CC79A7",
  "Normal-like" = "#009E73",
  "Claudin-low" = "#56B4E9"
)

x_limits <- c(
  0.05,
  8
)

if (
  min(
    plot_data$lower_95_ci
  ) < x_limits[1] ||
  max(
    plot_data$upper_95_ci
  ) > x_limits[2]
) {
  stop(
    paste(
      "A confidence interval falls outside",
      "the planned figure range."
    )
  )
}


# 5. Create the final RFS forest plot

rfs_forest_plot <- ggplot2::ggplot(
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
    ggplot2::aes(
      shape = comparison_significance
    ),
    size = 3.2,
    stroke = 1
  ) +
  ggplot2::facet_wrap(
    ~ facet_label,
    ncol = 3,
    axes = "all_x",
    axis.labels = "all_x",
    drop = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = subtype_colors,
    guide = "none"
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Holm-adjusted p < 0.05" = 16,
      "Not Holm-significant" = 1
    ),
    name = NULL
  ) +
  ggplot2::scale_x_log10(
    limits = x_limits,
    breaks = c(
      0.05,
      0.1,
      0.25,
      0.5,
      1,
      2,
      4,
      8
    ),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  ggplot2::labs(
    title = paste(
      "Adjusted hazard ratios for relapse-free survival events",
      "by molecular subtype"
    ),
    subtitle = paste(
      paste(
        "Reference: Luminal A;",
        "models adjusted for age at diagnosis and NPI"
      ),
      "and stratified by METABRIC source cohort"
    ),
    x = paste(
      "Adjusted hazard ratio,",
      "logarithmic scale"
    ),
    y = NULL,
    caption = paste(
      paste(
        "Points show hazard ratios and horizontal lines show",
        "pointwise 95% confidence intervals. Filled points",
        "identify comparisons that remain significant after",
        "Holm adjustment across all 25 subtype comparisons."
      ),
      paste(
        "Facet p-values are likelihood-ratio tests adjusted",
        "by Holm's method across the five interval-level",
        "global subtype tests."
      ),
      paste(
        "Estimates are conditional on remaining event-free",
        "to the start of each interval. The 3.5-year boundary",
        "is a post-diagnostic midpoint refinement."
      ),
      sep = "\n"
    )
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.spacing = grid::unit(
      1.2,
      "lines"
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 10.5,
      lineheight = 1.05
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
      size = 8.5,
      lineheight = 1.1,
      margin = ggplot2::margin(
        t = 10
      )
    ),
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = ggplot2::margin(
      t = 2,
      b = 2
    )
  )


# 6. Save and validate the final figure

figure_directory <- file.path(
  project_root,
  "output",
  "figures"
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_path <- file.path(
  figure_directory,
  "rfs_adjusted_hazard_ratios_by_interval.png"
)

ggplot2::ggsave(
  filename = figure_path,
  plot = rfs_forest_plot,
  width = 15,
  height = 8.5,
  units = "in",
  dpi = 320,
  bg = "white"
)

figure_information <- file.info(
  figure_path
)

if (
  !file.exists(
    figure_path
  ) ||
  is.na(
    figure_information$size
  ) ||
  figure_information$size <= 0
) {
  stop(
    "The final RFS figure was not created correctly."
  )
}

figure_summary <- data.frame(
  final_intervals = length(
    interval_order
  ),
  subtype_comparisons = nrow(
    subtype_estimates
  ),
  holm_significant_global_tests = sum(
    model_overview$holm_significant
  ),
  holm_significant_subtype_comparisons = sum(
    subtype_estimates$holm_significant
  ),
  stringsAsFactors = FALSE
)

print(
  rfs_forest_plot
)

cat(
  "\nFinal RFS figure summary:\n"
)

print(
  figure_summary,
  row.names = FALSE
)

cat(
  paste0(
    "\nThe figure presents the final results from Script 22. ",
    "It does not refit the models or add any tests ",
    "to the Holm families.\n"
  )
)

message("")
message(
  "Script 27 completed successfully."
)
message(
  "Figure: ",
  normalizePath(
    figure_path
  )
)