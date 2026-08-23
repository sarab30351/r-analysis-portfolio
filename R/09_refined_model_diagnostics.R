#
# 09_refined_model_diagnostics.R
#
# Purpose: Examine whether the remaining statistically significant PH tests represent meaningful changes over time or minor deviations detected because of the sample size.


suppressPackageStartupMessages(
  library(survival)
)

# Check project location 
if (!file.exists("r-analysis-portfolio.Rproj")) {
  stop(
    paste(
      "Open r-analysis-portfolio.Rproj before running this script.",
      "The working directory must be the repository root."
    )
  )
}


# Load refined model results 

model_path <- file.path(
  "data-derived",
  "os_refined_interval_cox_models.rds"
)

if (!file.exists(model_path)) {
  stop(
    paste(
      "Refined interval models not found at:",
      model_path
    )
  )
}

refined_models <- readRDS(model_path)
interval_results <- refined_models$interval_results

required_intervals <- c(
  "two_to_five",
  "five_to_ten",
  "beyond_ten"
)

missing_intervals <- setdiff(
  required_intervals,
  names(interval_results)
)

if (length(missing_intervals) > 0) {
  stop(
    paste(
      "Expected interval result(s) not found:",
      paste(missing_intervals, collapse = ", ")
    )
  )
}


# Function for selecting and plotting one PH curve 

plot_selected_ph_curve <- function(
    ph_object,
    coefficient_pattern,
    panel_title
) {
  
  matching_columns <- grep(
    coefficient_pattern,
    colnames(ph_object$y)
  )
  
  if (length(matching_columns) != 1) {
    stop(
      paste(
        "Expected one matching PH coefficient for:",
        panel_title,
        "but found",
        length(matching_columns)
      )
    )
  }
  
  graphics::plot(
    ph_object,
    var = matching_columns,
    resid = FALSE,
    se = TRUE,
    xlab = "Years since diagnosis",
    ylab = "Estimated coefficient",
    main = panel_title,
    lwd = 2
  )
  
  graphics::abline(
    h = 0,
    col = "gray60",
    lty = 3
  )
}


# Create diagnostic figure 

figure_directory <- file.path(
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
  "os_refined_ph_diagnostics.png"
)

grDevices::png(
  filename = figure_path,
  width = 2600,
  height = 1800,
  res = 200
)

old_graphics_parameters <- graphics::par(
  no.readonly = TRUE
)

graphics::par(
  mfrow = c(2, 3),
  mar = c(4.5, 4.5, 3.5, 1.5),
  oma = c(0, 0, 3, 0)
)


# Basal-like coefficient during 2 to 5 years

plot_selected_ph_curve(
  ph_object =
    interval_results$two_to_five$ph_coefficients,
  coefficient_pattern =
    "^molecular_subtypeBasal-like$",
  panel_title =
    "Basal-like vs Luminal A: 2 to 5 years"
)


# NPI term during 2 to 5 years

plot_selected_ph_curve(
  ph_object =
    interval_results$two_to_five$ph_terms,
  coefficient_pattern =
    "^splines::ns\\(npi,",
  panel_title =
    "NPI: 2 to 5 years"
)


# Age term during 5 to 10 years

plot_selected_ph_curve(
  ph_object =
    interval_results$five_to_ten$ph_terms,
  coefficient_pattern =
    "^splines::ns\\(age_at_diagnosis,",
  panel_title =
    "Age at diagnosis: 5 to 10 years"
)


# Age term beyond 10 years

plot_selected_ph_curve(
  ph_object =
    interval_results$beyond_ten$ph_terms,
  coefficient_pattern =
    "^splines::ns\\(age_at_diagnosis,",
  panel_title =
    "Age at diagnosis: beyond 10 years"
)


# Claudin-low coefficient beyond 10 years

plot_selected_ph_curve(
  ph_object =
    interval_results$beyond_ten$ph_coefficients,
  coefficient_pattern =
    "^molecular_subtypeClaudin-low$",
  panel_title =
    "Claudin-low vs Luminal A: beyond 10 years"
)


# Leave the final panel empty

graphics::plot.new()

graphics::mtext(
  "Selected proportional hazards diagnostics after interval refinement",
  outer = TRUE,
  side = 3,
  line = 1,
  font = 2,
  cex = 1.3
)

graphics::par(old_graphics_parameters)
grDevices::dev.off()


# Report output 

cat(
  "Refined diagnostic figure saved to: ",
  normalizePath(figure_path),
  "\n",
  sep = ""
)