# 15_rfs_descriptive_survival.R
#
# Purpose: Validate the relapse-free survival cohort and produce unadjusted descriptive survival results by molecular subtype.
#
# The METABRIC RFS endpoint combines locoregional relapse (breast, chest wall, lymph nodes), distant relapse (any other organ or tissue in the body) and disease-specific death (only deaths caused by breast cancer, but not necessarily with recorded recurrence).


# Check packages and project location

required_packages <- c(
  "survival",
  "ggplot2",
  "patchwork",
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
    "Install the following package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

if (!file.exists("r-analysis-portfolio.Rproj")) {
  stop(
    "Open r-analysis-portfolio.Rproj before running this script."
  )
}


# Load the final RFS cohort

rfs_cohort_path <- file.path(
  "data-derived",
  "metabric_rfs_cohort.rds"
)

if (!file.exists(rfs_cohort_path)) {
  stop("Run Script 03 before this script.")
}

rfs_data <- readRDS(rfs_cohort_path)

if (!"rfs_years" %in% names(rfs_data)) {
  rfs_data$rfs_years <- rfs_data$rfs_months / 12
}


# Check required variables

required_variables <- c(
  "patient_id",
  "rfs_months",
  "rfs_years",
  "rfs_event",
  "rfs_zero_event",
  "rfs_zero_censored",
  "age_at_diagnosis",
  "npi",
  "cohort",
  "molecular_subtype",
  "tumor_stage",
  "stage_missing",
  "stage_iv_recorded",
  "breast_surgery",
  "os_months",
  "os_event"
)

missing_variables <- setdiff(
  required_variables,
  names(rfs_data)
)

if (length(missing_variables) > 0) {
  stop(
    "Missing required variable(s): ",
    paste(missing_variables, collapse = ", ")
  )
}


# Validate identifiers and model variables

core_variables <- c(
  "patient_id",
  "rfs_months",
  "rfs_years",
  "rfs_event",
  "age_at_diagnosis",
  "npi",
  "cohort",
  "molecular_subtype"
)

if (!all(stats::complete.cases(
  rfs_data[, core_variables]
))) {
  stop("At least one required RFS variable is missing.")
}

if (anyDuplicated(rfs_data$patient_id)) {
  stop("Patient identifiers are not unique.")
}

if (any(rfs_data$rfs_months < 0)) {
  stop("A negative RFS follow-up time was detected.")
}

if (any(!rfs_data$rfs_event %in% c(0L, 1L))) {
  stop("RFS event coding must contain only zero and one.")
}


# Preserve the molecular-subtype order

subtype_levels <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

unexpected_subtypes <- setdiff(
  unique(as.character(rfs_data$molecular_subtype)),
  subtype_levels
)

if (length(unexpected_subtypes) > 0) {
  stop(
    "Unexpected molecular subtype(s): ",
    paste(unexpected_subtypes, collapse = ", ")
  )
}

rfs_data$molecular_subtype <- factor(
  rfs_data$molecular_subtype,
  levels = subtype_levels
)

rfs_data$cohort <- factor(rfs_data$cohort)


# Validate month-zero and stage indicators

zero_event_check <- (
  rfs_data$rfs_months == 0 &
    rfs_data$rfs_event == 1L
)

zero_censored_check <- (
  rfs_data$rfs_months == 0 &
    rfs_data$rfs_event == 0L
)

if (!all(rfs_data$rfs_zero_event == zero_event_check)) {
  stop("The month-zero RFS event indicator is inconsistent.")
}

if (!all(
  rfs_data$rfs_zero_censored == zero_censored_check
)) {
  stop("The month-zero censored indicator is inconsistent.")
}

if (!all(
  rfs_data$stage_missing == is.na(rfs_data$tumor_stage)
)) {
  stop("The missing-stage indicator is inconsistent.")
}


# Verify the final cohort counts

if (
  nrow(rfs_data) != 1960L ||
  sum(rfs_data$rfs_event == 1L) != 790L ||
  sum(rfs_data$rfs_event == 0L) != 1170L ||
  sum(zero_event_check) != 3L ||
  sum(zero_censored_check) != 0L ||
  any(rfs_data$stage_iv_recorded)
) {
  stop("The RFS cohort differs from the validated cohort.")
}


# Create the cohort-audit summary

rfs_cohort_audit <- data.frame(
  check = c(
    "Participants",
    "RFS events",
    "Censored observations",
    "Month-zero RFS events retained",
    "Month-zero censored observations",
    "Recorded Stage IV cases",
    "Participants with missing stage",
    "Missing age values",
    "Missing NPI values",
    "Duplicated patient identifiers"
  ),
  value = c(
    nrow(rfs_data),
    sum(rfs_data$rfs_event == 1L),
    sum(rfs_data$rfs_event == 0L),
    sum(zero_event_check),
    sum(zero_censored_check),
    sum(rfs_data$stage_iv_recorded),
    sum(rfs_data$stage_missing),
    sum(is.na(rfs_data$age_at_diagnosis)),
    sum(is.na(rfs_data$npi)),
    anyDuplicated(rfs_data$patient_id)
  )
)


# Review the three retained month-zero events

month_zero_events <- rfs_data[
  zero_event_check,
  c(
    "patient_id",
    "molecular_subtype",
    "tumor_stage",
    "stage_missing",
    "breast_surgery",
    "rfs_months",
    "os_months",
    "os_event"
  ),
  drop = FALSE
]

month_zero_events$molecular_subtype <- as.character(
  month_zero_events$molecular_subtype
)


# Print the audit results

cat("\nRFS cohort audit:\n")
print(rfs_cohort_audit, row.names = FALSE)

cat("\nRetained month-zero RFS events:\n")
print(month_zero_events, row.names = FALSE)

cat("\nRecorded tumor-stage distribution:\n")
print(
  table(
    rfs_data$tumor_stage,
    useNA = "always"
  )
)


# Summarize RFS outcomes by molecular subtype

rfs_subtype_summary <- do.call(
  rbind,
  lapply(
    subtype_levels,
    function(subtype) {
      subtype_data <- rfs_data[
        rfs_data$molecular_subtype == subtype,
        ,
        drop = FALSE
      ]
      
      participants <- nrow(subtype_data)
      events <- sum(subtype_data$rfs_event == 1L)
      
      data.frame(
        molecular_subtype = subtype,
        participants = participants,
        rfs_events = events,
        censored = sum(subtype_data$rfs_event == 0L),
        event_percent = 100 * events / participants,
        row.names = NULL
      )
    }
  )
)

if (
  sum(rfs_subtype_summary$participants) != nrow(rfs_data) ||
  sum(rfs_subtype_summary$rfs_events) !=
  sum(rfs_data$rfs_event)
) {
  stop("Subtype outcome totals do not match the RFS cohort.")
}


# Estimate follow-up using reverse Kaplan-Meier

summarize_reverse_km <- function(
    data,
    group_label
) {
  reverse_survival_object <- survival::Surv(
    time = data$rfs_years,
    event = 1L - data$rfs_event
  )
  
  reverse_km_fit <- survival::survfit(
    reverse_survival_object ~ 1
  )
  
  fit_table <- summary(reverse_km_fit)$table
  
  data.frame(
    group = group_label,
    participants = nrow(data),
    total_observed_person_years = sum(data$rfs_years),
    median_followup_years = unname(
      fit_table["median"]
    ),
    lower_95_ci = unname(
      fit_table["0.95LCL"]
    ),
    upper_95_ci = unname(
      fit_table["0.95UCL"]
    ),
    row.names = NULL
  )
}

followup_groups <- c(
  "Overall",
  subtype_levels
)

rfs_followup_summary <- do.call(
  rbind,
  lapply(
    followup_groups,
    function(group_label) {
      group_data <- if (group_label == "Overall") {
        rfs_data
      } else {
        rfs_data[
          rfs_data$molecular_subtype == group_label,
          ,
          drop = FALSE
        ]
      }
      
      summarize_reverse_km(
        data = group_data,
        group_label = group_label
      )
    }
  )
)


# Print descriptive results

rfs_subtype_display <- rfs_subtype_summary
rfs_subtype_display$event_percent <- round(
  rfs_subtype_display$event_percent,
  1
)

rfs_followup_display <- rfs_followup_summary

followup_numeric_columns <- c(
  "total_observed_person_years",
  "median_followup_years",
  "lower_95_ci",
  "upper_95_ci"
)

rfs_followup_display[followup_numeric_columns] <- lapply(
  rfs_followup_display[followup_numeric_columns],
  round,
  digits = 2
)

cat("\nRFS outcomes by molecular subtype:\n")
print(rfs_subtype_display, row.names = FALSE)

cat("\nReverse Kaplan-Meier follow-up summary:\n")
print(rfs_followup_display, row.names = FALSE)

# Fit unadjusted Kaplan-Meier RFS curves

rfs_survival_object <- survival::Surv(
  time = rfs_data$rfs_years,
  event = rfs_data$rfs_event
)

rfs_km_fit <- survival::survfit(
  rfs_survival_object ~ molecular_subtype,
  data = rfs_data,
  conf.type = "log-log"
)


# Perform the global unadjusted log-rank test

rfs_logrank_test <- survival::survdiff(
  rfs_survival_object ~ molecular_subtype,
  data = rfs_data
)

logrank_degrees_freedom <- length(
  rfs_logrank_test$n
) - 1L

logrank_p_value <- stats::pchisq(
  rfs_logrank_test$chisq,
  df = logrank_degrees_freedom,
  lower.tail = FALSE
)

rfs_logrank_results <- data.frame(
  endpoint = "Relapse-free survival",
  test = "Global log-rank test",
  chi_square = unname(rfs_logrank_test$chisq),
  degrees_freedom = logrank_degrees_freedom,
  p_value = logrank_p_value,
  row.names = NULL
)


# Estimate RFS probabilities at 5 and 10 years

rfs_time_summary <- summary(
  rfs_km_fit,
  times = c(5, 10),
  extend = FALSE
)

rfs_probabilities <- data.frame(
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    as.character(rfs_time_summary$strata)
  ),
  time_years = rfs_time_summary$time,
  n_risk = rfs_time_summary$n.risk,
  rfs_probability = rfs_time_summary$surv,
  lower_95_ci = rfs_time_summary$lower,
  upper_95_ci = rfs_time_summary$upper,
  row.names = NULL
)


# Estimate median RFS where estimable

median_table <- summary(rfs_km_fit)$table

median_rfs <- data.frame(
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    rownames(median_table)
  ),
  participants = unname(
    median_table[, "records"]
  ),
  rfs_events = unname(
    median_table[, "events"]
  ),
  median_rfs_years = unname(
    median_table[, "median"]
  ),
  lower_95_ci = unname(
    median_table[, "0.95LCL"]
  ),
  upper_95_ci = unname(
    median_table[, "0.95UCL"]
  ),
  row.names = NULL
)


# Calculate numbers at risk at five-year intervals

maximum_observed_years <- max(
  rfs_data$rfs_years
)

last_risk_table_year <- floor(
  maximum_observed_years / 5
) * 5

risk_table_times <- seq(
  from = 0,
  to = last_risk_table_year,
  by = 5
)

risk_summary <- summary(
  rfs_km_fit,
  times = risk_table_times,
  extend = TRUE
)

rfs_risk_table <- data.frame(
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    as.character(risk_summary$strata)
  ),
  time_years = risk_summary$time,
  n_risk = risk_summary$n.risk,
  row.names = NULL
)


# Prepare readable console tables

rfs_logrank_display <- rfs_logrank_results

rfs_logrank_display$chi_square <- round(
  rfs_logrank_display$chi_square,
  3
)

rfs_logrank_display$p_value <- format.pval(
  rfs_logrank_display$p_value,
  digits = 4,
  eps = 0.001
)

rfs_probability_display <- rfs_probabilities

probability_columns <- c(
  "rfs_probability",
  "lower_95_ci",
  "upper_95_ci"
)

rfs_probability_display[probability_columns] <- lapply(
  rfs_probability_display[probability_columns],
  round,
  digits = 3
)

median_rfs_display <- median_rfs

median_columns <- c(
  "median_rfs_years",
  "lower_95_ci",
  "upper_95_ci"
)

median_rfs_display[median_columns] <- lapply(
  median_rfs_display[median_columns],
  round,
  digits = 2
)

rfs_risk_matrix <- stats::xtabs(
  n_risk ~ molecular_subtype + time_years,
  data = rfs_risk_table
)


# Print unadjusted RFS results

cat("\nGlobal unadjusted RFS log-rank test:\n")
print(rfs_logrank_display, row.names = FALSE)

cat("\nFive- and ten-year RFS probabilities:\n")
print(rfs_probability_display, row.names = FALSE)

cat("\nMedian RFS by molecular subtype:\n")
print(median_rfs_display, row.names = FALSE)

cat("\nNumbers at risk by year:\n")
print(rfs_risk_matrix)


# Define accessible subtype colors

subtype_colors <- c(
  "Luminal A" = "#0072B2",
  "Luminal B" = "#E69F00",
  "HER2-enriched" = "#D55E00",
  "Basal-like" = "#CC79A7",
  "Normal-like" = "#009E73",
  "Claudin-low" = "#56B4E9"
)


# Extract Kaplan-Meier curve coordinates

rfs_km_summary <- summary(rfs_km_fit)

rfs_km_coordinates <- data.frame(
  time = rfs_km_summary$time,
  survival = rfs_km_summary$surv,
  lower = rfs_km_summary$lower,
  upper = rfs_km_summary$upper,
  n_risk = rfs_km_summary$n.risk,
  n_event = rfs_km_summary$n.event,
  n_censor = rfs_km_summary$n.censor,
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    as.character(rfs_km_summary$strata)
  ),
  plot_order = 1L,
  row.names = NULL
)

rfs_km_ci_data <- rfs_km_coordinates[
  !is.na(rfs_km_coordinates$lower) &
    !is.na(rfs_km_coordinates$upper),
  ,
  drop = FALSE
]


# Add the initial survival value for each subtype

initial_km_rows <- data.frame(
  time = 0,
  survival = 1,
  lower = 1,
  upper = 1,
  n_risk = rfs_subtype_summary$participants,
  n_event = 0,
  n_censor = 0,
  molecular_subtype =
    rfs_subtype_summary$molecular_subtype,
  plot_order = 0L,
  row.names = NULL
)

rfs_km_plot_data <- rbind(
  initial_km_rows,
  rfs_km_coordinates
)

rfs_km_plot_data$molecular_subtype <- factor(
  rfs_km_plot_data$molecular_subtype,
  levels = subtype_levels
)

rfs_km_plot_data <- rfs_km_plot_data[
  order(
    rfs_km_plot_data$molecular_subtype,
    rfs_km_plot_data$time,
    rfs_km_plot_data$plot_order
  ),
  ,
  drop = FALSE
]

rfs_km_ci_data$molecular_subtype <- factor(
  rfs_km_ci_data$molecular_subtype,
  levels = subtype_levels
)


# Restrict the displayed figure to 20 years

plot_risk_times <- c(0, 5, 10, 15, 20)

plot_risk_table_data <- rfs_risk_table[
  rfs_risk_table$time_years %in% plot_risk_times,
  ,
  drop = FALSE
]

plot_risk_table_data$molecular_subtype <- factor(
  plot_risk_table_data$molecular_subtype,
  levels = subtype_levels
)


# Create the Kaplan-Meier panel

rfs_km_plot <- ggplot2::ggplot(
  rfs_km_plot_data,
  ggplot2::aes(
    x = time,
    color = molecular_subtype
  )
) +
  ggplot2::geom_step(
    data = rfs_km_ci_data,
    ggplot2::aes(y = lower),
    linewidth = 0.3,
    alpha = 0.25,
    show.legend = FALSE
  ) +
  ggplot2::geom_step(
    data = rfs_km_ci_data,
    ggplot2::aes(y = upper),
    linewidth = 0.3,
    alpha = 0.25,
    show.legend = FALSE
  ) +
  ggplot2::geom_step(
    ggplot2::aes(y = survival),
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    data = rfs_km_plot_data[
      rfs_km_plot_data$n_censor > 0,
      ,
      drop = FALSE
    ],
    ggplot2::aes(y = survival),
    shape = 3,
    size = 0.8,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = subtype_colors,
    breaks = subtype_levels
  ) +
  ggplot2::scale_x_continuous(
    breaks = plot_risk_times,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-0.3, 20.3),
    expand = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Relapse-free survival by molecular subtype",
    subtitle = paste(
      "METABRIC RFS cohort;",
      "composite endpoint"
    ),
    x = NULL,
    y = "Relapse-free survival probability",
    color = "Molecular subtype",
    caption = paste(
      "Thin lines show 95% confidence limits;",
      "plus signs indicate censoring.",
      "Display restricted to 20 years because",
      "risk sets are sparse thereafter."
    )
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title.position = "plot"
  )


# Create the numbers-at-risk panel

rfs_risk_table_plot <- ggplot2::ggplot(
  plot_risk_table_data,
  ggplot2::aes(
    x = time_years,
    y = molecular_subtype,
    label = n_risk,
    color = molecular_subtype
  )
) +
  ggplot2::geom_text(
    size = 3.2,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = subtype_colors
  ) +
  ggplot2::scale_x_continuous(
    breaks = plot_risk_times,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-0.3, 20.3)
  ) +
  ggplot2::scale_y_discrete(
    limits = rev(subtype_levels)
  ) +
  ggplot2::labs(
    title = "Number at risk",
    x = "Years since diagnosis",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(
      color = "black"
    ),
    plot.title = ggplot2::element_text(
      size = 10,
      face = "bold"
    )
  )


# Combine and save the figure

combined_rfs_km_plot <- patchwork::wrap_plots(
  rfs_km_plot,
  rfs_risk_table_plot,
  ncol = 1,
  heights = c(3.6, 1.4)
)

figure_directory <- file.path(
  "output",
  "figures"
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

rfs_figure_path <- file.path(
  figure_directory,
  "km_relapse_free_survival_by_subtype.png"
)

ggplot2::ggsave(
  filename = rfs_figure_path,
  plot = combined_rfs_km_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

message(
  "RFS Kaplan-Meier figure saved to: ",
  normalizePath(rfs_figure_path)
)

# Identify and verify the project root
project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

if (
  !file.exists(
    file.path(
      project_root,
      "r-analysis-portfolio.Rproj"
    )
  )
) {
  stop(
    paste(
      "The working directory is not the project root:",
      project_root
    )
  )
}



# Export descriptive RFS results

rfs_output_table_dir <- file.path(
  project_root,
  "output",
  "tables"
)

rfs_output_data_dir <- file.path(
  project_root,
  "data-derived"
)

dir.create(
  rfs_output_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  rfs_output_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

rfs_export_tables <- list(
  rfs_cohort_audit = as.data.frame(rfs_cohort_audit),
  rfs_month_zero_events = as.data.frame(month_zero_events),
  rfs_subtype_outcomes = as.data.frame(rfs_subtype_summary),
  rfs_followup_summary = as.data.frame(rfs_followup_summary),
  rfs_logrank_test = as.data.frame(rfs_logrank_results),
  rfs_survival_probabilities = as.data.frame(rfs_probabilities),
  rfs_median_survival = as.data.frame(median_rfs),
  rfs_numbers_at_risk = as.data.frame(rfs_risk_table)
)

rfs_table_paths <- file.path(
  rfs_output_table_dir,
  paste0(names(rfs_export_tables), ".csv")
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
    rfs_export_tables,
    rfs_table_paths
  )
)

rfs_descriptive_results <- list(
  cohort_audit = rfs_cohort_audit,
  month_zero_events = month_zero_events,
  subtype_outcomes = rfs_subtype_summary,
  followup_summary = rfs_followup_summary,
  logrank_test = rfs_logrank_results,
  survival_probabilities = rfs_probabilities,
  median_survival = median_rfs,
  numbers_at_risk = rfs_risk_table,
  kaplan_meier_fit = rfs_km_fit,
  figure_file = file.path(
    "output",
    "figures",
    "km_relapse_free_survival_by_subtype.png"
  )
)

rfs_results_path <- file.path(
  rfs_output_data_dir,
  "rfs_descriptive_survival.rds"
)

saveRDS(
  rfs_descriptive_results,
  rfs_results_path
)

expected_table_rows <- c(
  rfs_cohort_audit = 10L,
  rfs_month_zero_events = 3L,
  rfs_subtype_outcomes = 6L,
  rfs_followup_summary = 7L,
  rfs_logrank_test = 1L,
  rfs_survival_probabilities = 12L,
  rfs_median_survival = 6L,
  rfs_numbers_at_risk = 36L
)

observed_table_rows <- vapply(
  rfs_table_paths,
  function(path) {
    nrow(utils::read.csv(path))
  },
  FUN.VALUE = integer(1)
)

stopifnot(
  identical(
    unname(observed_table_rows),
    unname(expected_table_rows)
  )
)

message("\nScript 15 completed successfully.")
message("Descriptive results: ", rfs_results_path)
message("Figure: ", rfs_figure_path)
message("Tables:")
message(
  paste0(
    "  ",
    rfs_table_paths,
    collapse = "\n"
  )
)