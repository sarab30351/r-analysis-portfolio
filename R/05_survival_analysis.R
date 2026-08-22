
# 05_survival_analysis.R
#
# Purpose:
# Perform unadjusted overall-survival analysis by molecular subtype.


# Import the overall-survival cohort 

os_file <- here::here(
  "data-derived",
  "metabric_os_cohort.rds"
)

if (!file.exists(os_file)) {
  stop(
    "The overall-survival cohort is missing. ",
    "Run R/03_define_cohorts.R first."
  )
}

os_cohort <- readRDS(os_file) |>
  dplyr::mutate(
    os_years = os_months / 12
  )


# Define subtype order and accessible colors 

subtype_levels <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal-like",
  "Normal-like",
  "Claudin-low"
)

subtype_colors <- c(
  "Luminal A" = "#0072B2",
  "Luminal B" = "#E69F00",
  "HER2-enriched" = "#D55E00",
  "Basal-like" = "#CC79A7",
  "Normal-like" = "#009E73",
  "Claudin-low" = "#56B4E9"
)

os_cohort <- os_cohort |>
  dplyr::mutate(
    molecular_subtype = factor(
      molecular_subtype,
      levels = subtype_levels
    )
  )


# Fit unadjusted Kaplan-Meier curves 

os_survival_object <- survival::Surv(
  time = os_cohort$os_years,
  event = os_cohort$os_event
)

os_km_fit <- survival::survfit(
  os_survival_object ~ molecular_subtype,
  data = os_cohort,
  conf.type = "log-log"
)


# Extract Kaplan-Meier curve coordinates 

os_km_summary <- summary(os_km_fit)

os_km_data <- data.frame(
  time = os_km_summary$time,
  survival = os_km_summary$surv,
  lower = os_km_summary$lower,
  upper = os_km_summary$upper,
  n_risk = os_km_summary$n.risk,
  n_event = os_km_summary$n.event,
  n_censor = os_km_summary$n.censor,
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    as.character(os_km_summary$strata)
  )
)

os_km_ci_data <- os_km_data |>
  dplyr::filter(
    !is.na(lower),
    !is.na(upper)
  )

initial_rows <- os_cohort |>
  dplyr::count(
    molecular_subtype,
    .drop = FALSE,
    name = "n_risk"
  ) |>
  dplyr::transmute(
    time = 0,
    survival = 1,
    lower = 1,
    upper = 1,
    n_risk = n_risk,
    n_event = 0,
    n_censor = 0,
    molecular_subtype = as.character(molecular_subtype)
  )

os_km_data <- dplyr::bind_rows(
  initial_rows,
  os_km_data
) |>
  dplyr::mutate(
    molecular_subtype = factor(
      molecular_subtype,
      levels = subtype_levels
    )
  ) |>
  dplyr::arrange(
    molecular_subtype,
    time
  )


# Select x-axis and risk-table times 

maximum_observed_years <- max(
  os_cohort$os_years,
  na.rm = TRUE
)

last_risk_table_year <- floor(
  maximum_observed_years / 5
) * 5

risk_table_times <- seq(
  from = 0,
  to = last_risk_table_year,
  by = 5
)

# Create risk-table data 
risk_summary <- summary(
  os_km_fit,
  times = risk_table_times,
  extend = TRUE
)

risk_table_data <- data.frame(
  time = risk_summary$time,
  n_risk = risk_summary$n.risk,
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    as.character(risk_summary$strata)
  )
) |>
  dplyr::mutate(
    molecular_subtype = factor(
      molecular_subtype,
      levels = subtype_levels
    )
  )

# Create Kaplan-Meier plot

os_km_plot <- ggplot2::ggplot(
  os_km_data,
  ggplot2::aes(
    x = time,
    color = molecular_subtype
  )
) +
  ggplot2::geom_step(
    data = os_km_ci_data,
    ggplot2::aes(y = lower),
    linewidth = 0.3,
    alpha = 0.25,
    show.legend = FALSE
  ) +
  ggplot2::geom_step(
    data = os_km_ci_data,
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
    data = dplyr::filter(
      os_km_data,
      n_censor > 0
    ),
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
    breaks = risk_table_times,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-0.4, max(risk_table_times) + 0.4),
    expand = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Overall survival by molecular subtype",
    subtitle = "METABRIC overall survival cohort",
    x = NULL,
    y = "Overall survival probability",
    color = "Molecular subtype",
    caption = paste(
      "Thin lines show 95% confidence limits; plus signs indicate censoring.",
      paste0(
        "Display restricted to ",
        max(risk_table_times),
        " years because risk sets are sparse thereafter."
      )
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


# Create the numbers at risk table

os_risk_table_plot <- ggplot2::ggplot(
  risk_table_data,
  ggplot2::aes(
    x = time,
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
    breaks = risk_table_times,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-0.4, max(risk_table_times) + 0.4)
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
  
# Combine and save the Kaplan-Meier figure 

combined_os_km_plot <- patchwork::wrap_plots(
  os_km_plot,
  os_risk_table_plot,
  ncol = 1,
  heights = c(3.6, 1.4)
)

figure_directory <- here::here(
  "output",
  "figures"
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

ggplot2::ggsave(
  filename = here::here(
    "output",
    "figures",
    "km_overall_survival_by_subtype.png"
  ),
  plot = combined_os_km_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)


# Perform the global log-rank test 

os_logrank_test <- survival::survdiff(
  os_survival_object ~ molecular_subtype,
  data = os_cohort
)

logrank_degrees_freedom <- length(
  os_logrank_test$n
) - 1L

logrank_p_value <- stats::pchisq(
  os_logrank_test$chisq,
  df = logrank_degrees_freedom,
  lower.tail = FALSE
)

logrank_results <- data.frame(
  endpoint = "Overall survival",
  test = "Global log-rank test",
  chi_square = unname(os_logrank_test$chisq),
  degrees_freedom = logrank_degrees_freedom,
  p_value = logrank_p_value
)


# Estimate survival probabilities at 5 and 10 years 

survival_time_summary <- summary(
  os_km_fit,
  times = c(5, 10),
  extend = FALSE
)

survival_probabilities <- data.frame(
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    as.character(survival_time_summary$strata)
  ),
  time_years = survival_time_summary$time,
  n_risk = survival_time_summary$n.risk,
  survival_probability = survival_time_summary$surv,
  lower_95_ci = survival_time_summary$lower,
  upper_95_ci = survival_time_summary$upper
)


# Estimate median survival  

median_table <- summary(os_km_fit)$table

median_survival <- data.frame(
  molecular_subtype = sub(
    "^molecular_subtype=",
    "",
    rownames(median_table)
  ),
  participants = unname(
    median_table[, "records"]
  ),
  events = unname(
    median_table[, "events"]
  ),
  median_survival_years = unname(
    median_table[, "median"]
  ),
  lower_95_ci = unname(
    median_table[, "0.95LCL"]
  ),
  upper_95_ci = unname(
    median_table[, "0.95UCL"]
  )
)


# Save results 

table_directory <- here::here(
  "output",
  "tables"
)

dir.create(
  table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  logrank_results,
  here::here(
    "output",
    "tables",
    "os_logrank_test.csv"
  )
)

readr::write_csv(
  survival_probabilities,
  here::here(
    "output",
    "tables",
    "os_survival_probabilities.csv"
  )
)

readr::write_csv(
  median_survival,
  here::here(
    "output",
    "tables",
    "os_median_survival.csv"
  )
)


# Report completion 

message("Unadjusted overall-survival analysis completed.")
message(
  "Global log-rank p-value: ",
  format.pval(
    logrank_p_value,
    digits = 3,
    eps = 0.001
  )
)
message(
  "Kaplan-Meier figure saved to ",
  here::here(
    "output",
    "figures",
    "km_overall_survival_by_subtype.png"
  )
)

# Three late follow-up points were warnings because they were at the extreme end of follow-up i.e very few participants at risk (confidence limits NA) so I cut off the confidence lines right before those points.

# os_km_ci_data <- os_km_data |>
  # dplyr::filter(
   # !is.na(lower),
   # !is.na(upper)
   # )
# data = os_km_ci_data