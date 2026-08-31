# ==============================================================================
# Script 166 FINAL
# Supplementary Figure S8
# GSE185263 external sensitivity and secondary analyses
#
# Purpose:
#   Publication packaging of FROZEN external-validation results from:
#     - Script 142b: GSE185263 external validation
#     - Script 143: location-specific confidence intervals / descriptive pooling
#
# Figure S8 panels:
#   A. Five-gene score by in-hospital mortality
#   B. Five-gene score by collection site: ICU vs Emergency Room
#   C. Location-specific score-SOFA correlations with 95% CI
#   D. Contextual sepsis-versus-healthy contrast
#
# IMPORTANT:
#   - NO new inferential statistical tests are performed.
#   - All P values, q values, AUCs and confidence intervals are frozen.
#   - Sepsis-vs-healthy discrimination is contextual only.
#   - Collection locations are sensitivity strata within one dataset,
#     NOT independent external validation cohorts.
#   - The pooled location estimate is descriptive only.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. START
# ------------------------------------------------------------------------------

cat("\n")
cat("====================================================================\n")
cat("Running Script 166 FINAL\n")
cat("Supplementary Figure S8\n")
cat("GSE185263 external sensitivity and secondary analyses\n")
cat("====================================================================\n\n")


# ------------------------------------------------------------------------------
# 1. PROJECT DIRECTORY
# ------------------------------------------------------------------------------

project_candidates <- c(
  Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2")),
  path.expand("~/Sepsis_DESeq2")
)

project_dir <- project_candidates[file.exists(project_candidates)][1]

if (is.na(project_dir) || length(project_dir) == 0) {
  stop(
    "Could not locate Sepsis_DESeq2 project directory.\n",
    "Checked:\n",
    paste(project_candidates, collapse = "\n")
  )
}

setwd(project_dir)

cat("Project directory:\n")
print(normalizePath(project_dir, winslash = "\\", mustWork = TRUE))
cat("\n")


# ------------------------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "ggplot2",
  "openxlsx",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
})

cat("Required packages loaded successfully.\n\n")


# ------------------------------------------------------------------------------
# 3. INPUT DIRECTORIES
# ------------------------------------------------------------------------------

dir_142b <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "142b_external_validation_GSE185263",
  "tables"
)

dir_143 <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "143_multicohort_integration",
  "tables"
)

dir_150 <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "150_Figure5_GSE185263_severity_validation",
  "tables"
)


# ------------------------------------------------------------------------------
# 4. INPUT FILES
# ------------------------------------------------------------------------------

file_sample_scores <- file.path(
  dir_142b,
  "142b_GSE185263_sample_scores.csv"
)

file_secondary <- file.path(
  dir_142b,
  "142b_secondary_score_associations.csv"
)

file_location_142b <- file.path(
  dir_142b,
  "142b_location_specific_SOFA_correlations.csv"
)

file_location_143 <- file.path(
  dir_143,
  "143_GSE185263_location_correlations_with_CI.csv"
)

file_location_meta_143 <- file.path(
  dir_143,
  "143_GSE185263_location_random_effects_meta.csv"
)

file_figure5_audit <- file.path(
  dir_150,
  "150_Figure5_numerical_audit.csv"
)


required_files <- c(
  file_sample_scores,
  file_secondary,
  file_location_142b,
  file_location_143,
  file_location_meta_143
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("Frozen input files found successfully.\n\n")


# ------------------------------------------------------------------------------
# 5. OUTPUT DIRECTORIES
# ------------------------------------------------------------------------------

output_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "166_FigureS8_GSE185263_external_sensitivity"
)

figure_dir <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")
text_dir   <- file.path(output_dir, "text")
audit_dir  <- file.path(output_dir, "audit")
log_dir    <- file.path(output_dir, "logs")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,    recursive = TRUE, showWarnings = FALSE)

cat("Output folder:\n")
cat(normalizePath(output_dir, winslash = "\\", mustWork = TRUE), "\n\n")


# ------------------------------------------------------------------------------
# 6. READ FROZEN SOURCES
# ------------------------------------------------------------------------------

sample_scores <- readr::read_csv(
  file_sample_scores,
  show_col_types = FALSE
)

secondary <- readr::read_csv(
  file_secondary,
  show_col_types = FALSE
)

location_142b <- readr::read_csv(
  file_location_142b,
  show_col_types = FALSE
)

location_143 <- readr::read_csv(
  file_location_143,
  show_col_types = FALSE
)

location_meta_143 <- readr::read_csv(
  file_location_meta_143,
  show_col_types = FALSE
)

figure5_audit <- NULL

if (file.exists(file_figure5_audit)) {
  figure5_audit <- readr::read_csv(
    file_figure5_audit,
    show_col_types = FALSE
  )
}


# ------------------------------------------------------------------------------
# 7. EXACT SCHEMA CHECKS
# ------------------------------------------------------------------------------

required_sample_columns <- c(
  "geo_accession",
  "sample_title",
  "disease_state",
  "age_numeric",
  "sex_standardized",
  "collection_location",
  "collection_site",
  "mortality",
  "sofa",
  "five_gene_score",
  "five_gene_score_all_reference",
  "five_gene_score_healthy_reference",
  "sofa_group"
)

required_secondary_columns <- c(
  "analysis",
  "score",
  "case",
  "control",
  "n_case",
  "n_control",
  "median_case",
  "median_control",
  "median_difference",
  "W",
  "p_value",
  "AUC",
  "CI_low",
  "CI_high",
  "BH_secondary"
)

required_location143_columns <- c(
  "collection_location",
  "direction_concordant",
  "n",
  "rho",
  "p_value",
  "BH_location",
  "CI_low",
  "CI_high"
)

required_meta_columns <- c(
  "k",
  "total_n",
  "Q",
  "Q_df",
  "Q_p",
  "I2_percent",
  "tau2_fisher_z",
  "fixed_rho",
  "fixed_CI_low",
  "fixed_CI_high",
  "fixed_p",
  "random_rho",
  "random_CI_low",
  "random_CI_high",
  "random_p"
)


check_columns <- function(data, required, object_name) {
  
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    stop(
      object_name,
      " is missing required columns:\n",
      paste(missing, collapse = ", ")
    )
  }
  
  TRUE
}


schema_audit <- data.frame(
  object = c(
    "sample_scores",
    "secondary_associations",
    "location_CI",
    "location_meta"
  ),
  schema_match = c(
    check_columns(
      sample_scores,
      required_sample_columns,
      "sample_scores"
    ),
    check_columns(
      secondary,
      required_secondary_columns,
      "secondary"
    ),
    check_columns(
      location_143,
      required_location143_columns,
      "location_143"
    ),
    check_columns(
      location_meta_143,
      required_meta_columns,
      "location_meta_143"
    )
  )
)

cat("SCHEMA AUDIT\n")
cat("------------\n")
print(schema_audit)
cat("\n")


# ------------------------------------------------------------------------------
# 8. FROZEN SECONDARY ENDPOINT ROWS
# ------------------------------------------------------------------------------

mortality_result <- secondary %>%
  dplyr::filter(analysis == "Died_vs_Survived")

icu_result <- secondary %>%
  dplyr::filter(analysis == "ICU_vs_Emergency_Room")

context_result <- secondary %>%
  dplyr::filter(analysis == "Sepsis_vs_healthy_contextual")


if (nrow(mortality_result) != 1) {
  stop("Expected exactly one Died_vs_Survived row.")
}

if (nrow(icu_result) != 1) {
  stop("Expected exactly one ICU_vs_Emergency_Room row.")
}

if (nrow(context_result) != 1) {
  stop("Expected exactly one Sepsis_vs_healthy_contextual row.")
}


# ------------------------------------------------------------------------------
# 9. SAMPLE-LEVEL DATASETS FOR PLOTTING
# ------------------------------------------------------------------------------

sepsis_samples <- sample_scores %>%
  dplyr::filter(tolower(disease_state) == "sepsis")

healthy_samples <- sample_scores %>%
  dplyr::filter(tolower(disease_state) == "healthy")


mortality_data <- sepsis_samples %>%
  dplyr::filter(mortality %in% c("Survived", "Died")) %>%
  dplyr::mutate(
    mortality_plot = factor(
      mortality,
      levels = c("Survived", "Died")
    )
  )


icu_data <- sepsis_samples %>%
  dplyr::filter(collection_site %in% c("Emergency Room", "ICU")) %>%
  dplyr::mutate(
    site_plot = factor(
      collection_site,
      levels = c("Emergency Room", "ICU")
    )
  )


context_data <- sample_scores %>%
  dplyr::filter(
    tolower(disease_state) %in% c("healthy", "sepsis")
  ) %>%
  dplyr::mutate(
    context_group = factor(
      ifelse(
        tolower(disease_state) == "healthy",
        "Healthy controls",
        "Sepsis"
      ),
      levels = c(
        "Healthy controls",
        "Sepsis"
      )
    )
  )


# ------------------------------------------------------------------------------
# 10. SAMPLE COUNT AUDIT
# ------------------------------------------------------------------------------

count_audit <- data.frame(
  metric = c(
    "Total samples",
    "Sepsis samples",
    "Healthy samples",
    "Mortality Died",
    "Mortality Survived",
    "ICU",
    "Emergency Room"
  ),
  expected = c(
    392,
    348,
    44,
    mortality_result$n_case,
    mortality_result$n_control,
    icu_result$n_case,
    icu_result$n_control
  ),
  observed = c(
    nrow(sample_scores),
    nrow(sepsis_samples),
    nrow(healthy_samples),
    sum(mortality_data$mortality == "Died"),
    sum(mortality_data$mortality == "Survived"),
    sum(icu_data$collection_site == "ICU"),
    sum(icu_data$collection_site == "Emergency Room")
  )
)

count_audit$pass <-
  abs(count_audit$expected - count_audit$observed) < 1e-12


cat("SAMPLE COUNT AUDIT\n")
cat("------------------\n")
print(count_audit)
cat("\n")

if (!all(count_audit$pass)) {
  stop("Sample count audit failed.")
}


# ------------------------------------------------------------------------------
# 11. MEDIAN CONSISTENCY AUDIT
# ------------------------------------------------------------------------------

observed_mort_died <- median(
  mortality_data$five_gene_score[
    mortality_data$mortality == "Died"
  ],
  na.rm = TRUE
)

observed_mort_survived <- median(
  mortality_data$five_gene_score[
    mortality_data$mortality == "Survived"
  ],
  na.rm = TRUE
)

observed_icu <- median(
  icu_data$five_gene_score[
    icu_data$collection_site == "ICU"
  ],
  na.rm = TRUE
)

observed_er <- median(
  icu_data$five_gene_score[
    icu_data$collection_site == "Emergency Room"
  ],
  na.rm = TRUE
)

observed_context_sepsis <- median(
  context_data$five_gene_score_all_reference[
    context_data$context_group == "Sepsis"
  ],
  na.rm = TRUE
)

observed_context_healthy <- median(
  context_data$five_gene_score_all_reference[
    context_data$context_group == "Healthy controls"
  ],
  na.rm = TRUE
)


median_audit <- data.frame(
  metric = c(
    "Mortality Died median",
    "Mortality Survived median",
    "ICU median",
    "Emergency Room median",
    "Contextual Sepsis median",
    "Contextual Healthy median"
  ),
  expected = c(
    mortality_result$median_case,
    mortality_result$median_control,
    icu_result$median_case,
    icu_result$median_control,
    context_result$median_case,
    context_result$median_control
  ),
  observed = c(
    observed_mort_died,
    observed_mort_survived,
    observed_icu,
    observed_er,
    observed_context_sepsis,
    observed_context_healthy
  )
)

median_audit$absolute_difference <-
  abs(median_audit$expected - median_audit$observed)

median_audit$pass <-
  median_audit$absolute_difference < 1e-8


cat("MEDIAN CONSISTENCY AUDIT\n")
cat("------------------------\n")
print(median_audit)
cat("\n")

if (!all(median_audit$pass)) {
  stop("Median consistency audit failed.")
}


# ------------------------------------------------------------------------------
# 12. NUMERICAL ANCHOR AUDIT
# ------------------------------------------------------------------------------

location_pooled <- location_meta_143[1, ]


anchor_audit <- data.frame(
  metric = c(
    "Mortality AUC",
    "Mortality P",
    "ICU vs ER AUC",
    "ICU vs ER P",
    "Contextual sepsis vs healthy AUC",
    "Location pooled rho",
    "Location total n"
  ),
  observed = c(
    mortality_result$AUC,
    mortality_result$p_value,
    icu_result$AUC,
    icu_result$p_value,
    context_result$AUC,
    location_pooled$fixed_rho,
    sum(location_143$n)
  ),
  expected = c(
    0.6272644,
    0.003447767,
    0.6438658,
    8.173562e-05,
    0.9494514,
    0.2407839,
    345
  ),
  tolerance = c(
    1e-6,
    1e-8,
    1e-6,
    1e-9,
    1e-6,
    1e-6,
    0
  )
)

anchor_audit$absolute_difference <-
  abs(anchor_audit$observed - anchor_audit$expected)

anchor_audit$pass <-
  anchor_audit$absolute_difference <= anchor_audit$tolerance


cat("FROZEN NUMERICAL ANCHOR AUDIT\n")
cat("-----------------------------\n")
print(anchor_audit)
cat("\n")

if (!all(anchor_audit$pass)) {
  stop("Frozen numerical anchor audit failed.")
}


# ------------------------------------------------------------------------------
# 13. LOCATION-SPECIFIC SENSITIVITY AUDIT
# ------------------------------------------------------------------------------

if (nrow(location_143) != 5) {
  stop("Expected exactly five collection-location strata.")
}

if (!all(location_143$direction_concordant)) {
  stop(
    "At least one location-specific association is not directionally concordant."
  )
}

if (!all(location_143$rho > 0)) {
  stop(
    "Expected all five frozen location-specific rho estimates to be positive."
  )
}


location_plot_data <- location_143 %>%
  dplyr::mutate(
    location_label = tools::toTitleCase(collection_location)
  ) %>%
  dplyr::arrange(rho) %>%
  dplyr::mutate(
    location_label = factor(
      location_label,
      levels = location_label
    )
  )


cat("LOCATION-SPECIFIC SENSITIVITY\n")
cat("-----------------------------\n")
print(location_143)
cat("\n")

cat(
  "All five location-specific rho estimates positive = ",
  all(location_143$rho > 0),
  "\n",
  sep = ""
)

cat(
  "Descriptive fixed-effect pooled rho = ",
  signif(location_pooled$fixed_rho, 7),
  "\n",
  sep = ""
)

cat(
  "95% CI = ",
  signif(location_pooled$fixed_CI_low, 7),
  " to ",
  signif(location_pooled$fixed_CI_high, 7),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 14. FORMATTING HELPERS
# ------------------------------------------------------------------------------

format_p <- function(x) {
  
  if (is.na(x)) {
    return("NA")
  }
  
  if (x == 0) {
    return("<1e-16")
  }
  
  if (x < 0.001) {
    return(formatC(
      x,
      format = "e",
      digits = 2
    ))
  }
  
  sprintf("%.4f", x)
}


format_auc <- function(x) {
  sprintf("%.3f", x)
}


# ------------------------------------------------------------------------------
# 15. PUBLICATION PALETTE
# ------------------------------------------------------------------------------

col_control <- "#4C78A8"
col_case    <- "#E45756"
col_context <- "#59A14F"
col_forest  <- "#3A6EA5"
col_pooled  <- "#7A5195"
col_neutral <- "#6B6B6B"


# ------------------------------------------------------------------------------
# 16. COMMON THEME
# ------------------------------------------------------------------------------

theme_publication <- function() {
  
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11.5,
        margin = ggplot2::margin(b = 7)
      ),
      axis.title = ggplot2::element_text(
        size = 10.5
      ),
      axis.text = ggplot2::element_text(
        size = 9.5
      ),
      axis.text.x = ggplot2::element_text(
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        color = "black"
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(
        t = 8,
        r = 10,
        b = 8,
        l = 8
      )
    )
}


# ------------------------------------------------------------------------------
# 17. PANEL A — IN-HOSPITAL MORTALITY
# ------------------------------------------------------------------------------

mort_label <- paste0(
  "AUC = ",
  format_auc(mortality_result$AUC),
  "\nq = ",
  format_p(mortality_result$BH_secondary)
)


panel_A <- ggplot2::ggplot(
  mortality_data,
  ggplot2::aes(
    x = mortality_plot,
    y = five_gene_score,
    fill = mortality_plot
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    alpha = 0.55,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Survived" = col_control,
      "Died" = col_case
    )
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = max(mortality_data$five_gene_score, na.rm = TRUE),
    label = mort_label,
    vjust = -0.55,
    size = 3.3
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.20)
    )
  ) +
  ggplot2::labs(
    title = "In-hospital mortality",
    x = NULL,
    y = "Five-gene host-response score"
  ) +
  theme_publication()


# ------------------------------------------------------------------------------
# 18. PANEL B — ICU vs EMERGENCY ROOM
# ------------------------------------------------------------------------------

icu_label <- paste0(
  "AUC = ",
  format_auc(icu_result$AUC),
  "\nq = ",
  format_p(icu_result$BH_secondary)
)


panel_B <- ggplot2::ggplot(
  icu_data,
  ggplot2::aes(
    x = site_plot,
    y = five_gene_score,
    fill = site_plot
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    alpha = 0.55,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Emergency Room" = col_control,
      "ICU" = col_case
    )
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = max(icu_data$five_gene_score, na.rm = TRUE),
    label = icu_label,
    vjust = -0.55,
    size = 3.3
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.20)
    )
  ) +
  ggplot2::labs(
    title = "Collection site",
    x = NULL,
    y = "Five-gene host-response score"
  ) +
  theme_publication()


# ------------------------------------------------------------------------------
# 19. PANEL C — LOCATION-SPECIFIC SCORE–SOFA ASSOCIATIONS
# ------------------------------------------------------------------------------

pooled_label <- paste0(
  "Descriptive pooled rho = ",
  sprintf("%.3f", location_pooled$fixed_rho)
)


panel_C <- ggplot2::ggplot(
  location_plot_data,
  ggplot2::aes(
    x = rho,
    y = location_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.55,
    color = col_neutral
  ) +
  ggplot2::geom_vline(
    xintercept = location_pooled$fixed_rho,
    linetype = "dashed",
    linewidth = 0.7,
    color = col_pooled
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.18,
    linewidth = 0.65,
    color = col_forest
  ) +
  ggplot2::geom_point(
    size = 2.7,
    color = col_forest
  ) +
  ggplot2::annotate(
    "text",
    x = location_pooled$fixed_rho,
    y = length(levels(location_plot_data$location_label)) + 0.58,
    label = pooled_label,
    color = col_pooled,
    hjust = 0.5,
    size = 3.1
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      min(
        -0.45,
        min(location_plot_data$CI_low, na.rm = TRUE) - 0.03
      ),
      max(
        0.80,
        max(location_plot_data$CI_high, na.rm = TRUE) + 0.03
      )
    ),
    breaks = seq(
      -0.4,
      0.8,
      by = 0.2
    )
  ) +
  ggplot2::coord_cartesian(
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Location-specific score–SOFA associations",
    x = "Spearman rho with 24-h SOFA",
    y = NULL
  ) +
  theme_publication() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 23,
      r = 10,
      b = 8,
      l = 8
    )
  )


# ------------------------------------------------------------------------------
# 20. PANEL D — CONTEXTUAL SEPSIS vs HEALTHY CONTRAST
# ------------------------------------------------------------------------------

context_label <- paste0(
  "AUC = ",
  format_auc(context_result$AUC),
  "\ncontextual only"
)


panel_D <- ggplot2::ggplot(
  context_data,
  ggplot2::aes(
    x = context_group,
    y = five_gene_score_all_reference,
    fill = context_group
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    alpha = 0.55,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Healthy controls" = col_control,
      "Sepsis" = col_context
    )
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = max(
      context_data$five_gene_score_all_reference,
      na.rm = TRUE
    ),
    label = context_label,
    vjust = -0.55,
    size = 3.3
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.20)
    )
  ) +
  ggplot2::labs(
    title = "Contextual sepsis–healthy contrast",
    x = NULL,
    y = "Five-gene score\n(all-sample reference)"
  ) +
  theme_publication()


# ------------------------------------------------------------------------------
# 21. COMBINE FIGURE
# ------------------------------------------------------------------------------

figure_S8 <- (
  panel_A | panel_B
) / (
  panel_C | panel_D
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Supplementary Figure S8. External sensitivity analyses in GSE185263",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      )
    )
  )


# ------------------------------------------------------------------------------
# 22. SAVE FIGURE
# ------------------------------------------------------------------------------

png_file <- file.path(
  figure_dir,
  "166_FigureS8_GSE185263_external_sensitivity.png"
)

pdf_file <- file.path(
  figure_dir,
  "166_FigureS8_GSE185263_external_sensitivity.pdf"
)

tiff_file <- file.path(
  figure_dir,
  "166_FigureS8_GSE185263_external_sensitivity.tiff"
)


ggplot2::ggsave(
  filename = png_file,
  plot = figure_S8,
  width = 11.5,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = pdf_file,
  plot = figure_S8,
  width = 11.5,
  height = 8.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggplot2::ggsave(
  filename = tiff_file,
  plot = figure_S8,
  width = 11.5,
  height = 8.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------------------------
# 23. FIGURE SOURCE-DATA WORKBOOK
# ------------------------------------------------------------------------------

provenance_table <- data.frame(
  figure_component = c(
    "Panel A: mortality",
    "Panel B: ICU vs ER",
    "Panel C: location-specific SOFA",
    "Panel C: descriptive pooled location estimate",
    "Panel D: contextual sepsis vs healthy"
  ),
  frozen_source = c(
    "142b_secondary_score_associations.csv + 142b_GSE185263_sample_scores.csv",
    "142b_secondary_score_associations.csv + 142b_GSE185263_sample_scores.csv",
    "143_GSE185263_location_correlations_with_CI.csv",
    "143_GSE185263_location_random_effects_meta.csv",
    "142b_secondary_score_associations.csv + 142b_GSE185263_sample_scores.csv"
  ),
  new_inferential_testing = rep("NO", 5),
  interpretation = c(
    "Secondary outcome context",
    "Secondary clinical-setting context",
    "Sensitivity analysis within one external dataset",
    "Descriptive only; not independent-cohort meta-analysis",
    "Contextual discrimination only; not primary clinical validation endpoint"
  ),
  stringsAsFactors = FALSE
)


readme_table <- data.frame(
  item = c(
    "Figure",
    "Dataset",
    "Purpose",
    "Primary external endpoint",
    "Panels",
    "New hypothesis testing",
    "Location interpretation",
    "Contextual discrimination interpretation"
  ),
  value = c(
    "Supplementary Figure S8",
    "GSE185263",
    "Secondary and sensitivity visualization supporting frozen external severity evaluation",
    "Continuous 24-h SOFA association is reported in Main Figure 5 and Supplementary Table S10",
    "A mortality; B ICU vs Emergency Room; C location-specific SOFA associations; D contextual sepsis vs healthy",
    "NO",
    "Five collection locations are strata within one dataset, not independent validation cohorts",
    "Sepsis-vs-healthy AUC is contextual only and does not establish a calibrated diagnostic assay"
  ),
  stringsAsFactors = FALSE
)


wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  wb,
  "S8_ReadMe"
)

openxlsx::writeData(
  wb,
  "S8_ReadMe",
  readme_table
)


openxlsx::addWorksheet(
  wb,
  "Secondary_frozen"
)

openxlsx::writeData(
  wb,
  "Secondary_frozen",
  secondary
)


openxlsx::addWorksheet(
  wb,
  "Mortality_samples"
)

openxlsx::writeData(
  wb,
  "Mortality_samples",
  mortality_data
)


openxlsx::addWorksheet(
  wb,
  "ICU_ER_samples"
)

openxlsx::writeData(
  wb,
  "ICU_ER_samples",
  icu_data
)


openxlsx::addWorksheet(
  wb,
  "Location_SOFA"
)

openxlsx::writeData(
  wb,
  "Location_SOFA",
  location_143
)


openxlsx::addWorksheet(
  wb,
  "Location_pooled"
)

openxlsx::writeData(
  wb,
  "Location_pooled",
  location_meta_143
)


openxlsx::addWorksheet(
  wb,
  "Contextual_samples"
)

openxlsx::writeData(
  wb,
  "Contextual_samples",
  context_data
)


openxlsx::addWorksheet(
  wb,
  "Provenance"
)

openxlsx::writeData(
  wb,
  "Provenance",
  provenance_table
)


openxlsx::addWorksheet(
  wb,
  "Anchor_audit"
)

openxlsx::writeData(
  wb,
  "Anchor_audit",
  anchor_audit
)


openxlsx::addWorksheet(
  wb,
  "Count_audit"
)

openxlsx::writeData(
  wb,
  "Count_audit",
  count_audit
)


openxlsx::addWorksheet(
  wb,
  "Median_audit"
)

openxlsx::writeData(
  wb,
  "Median_audit",
  median_audit
)


for (sheet_name in names(wb)) {
  
  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  openxlsx::setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:50,
    widths = "auto"
  )
}


source_data_file <- file.path(
  table_dir,
  "166_FigureS8_source_data.xlsx"
)

openxlsx::saveWorkbook(
  wb,
  source_data_file,
  overwrite = TRUE
)


# ------------------------------------------------------------------------------
# 24. FIGURE CAPTION — ENGLISH
# ------------------------------------------------------------------------------

caption_EN <- paste0(
  "Supplementary Figure S8. Secondary and sensitivity analyses of the ",
  "five-gene host-response score in GSE185263. ",
  "(A) Distribution of the frozen five-gene score according to in-hospital ",
  "mortality among patients with sepsis. Patients who died had higher scores ",
  "than survivors (AUC = ",
  sprintf("%.3f", mortality_result$AUC),
  "; Benjamini-Hochberg-adjusted q = ",
  format_p(mortality_result$BH_secondary),
  "). ",
  "(B) Distribution of the score according to collection site. ICU samples ",
  "showed higher scores than Emergency Room samples (AUC = ",
  sprintf("%.3f", icu_result$AUC),
  "; q = ",
  format_p(icu_result$BH_secondary),
  "). ",
  "(C) Location-specific Spearman correlations between the five-gene score ",
  "and continuous 24-h SOFA. Points indicate frozen correlation estimates ",
  "and horizontal bars indicate 95% confidence intervals. All five estimates ",
  "were positive, although magnitude and precision varied among locations. ",
  "The purple dashed line indicates the descriptive fixed-effect pooled ",
  "estimate (rho = ",
  sprintf("%.3f", location_pooled$fixed_rho),
  "); the collection locations are sensitivity strata within a single dataset ",
  "and are not independent external validation cohorts. ",
  "(D) Contextual comparison of sepsis and healthy samples using the ",
  "all-sample-reference version of the frozen five-gene score ",
  "(AUC = ",
  sprintf("%.3f", context_result$AUC),
  "). This contrast is shown for biological context only and was not the ",
  "primary clinically relevant external validation endpoint. ",
  "No new inferential statistical analyses were performed for this figure; ",
  "all statistics were taken from the frozen Script 142b/143 analyses."
)


caption_EN_file <- file.path(
  text_dir,
  "166_FigureS8_caption_EN.txt"
)

writeLines(
  caption_EN,
  con = caption_EN_file,
  useBytes = TRUE
)


# ------------------------------------------------------------------------------
# 25. FIGURE CAPTION — RUSSIAN
# ------------------------------------------------------------------------------

caption_RU <- paste0(
  "Дополнительный рисунок S8. Вторичные анализы и анализы чувствительности ",
  "пятигенного показателя ответа организма во внешнем наборе GSE185263. ",
  "(A) Распределение замороженного пятигенного показателя в зависимости от ",
  "внутрибольничной летальности у пациентов с сепсисом. У умерших пациентов ",
  "значения показателя были выше, чем у выживших (AUC = ",
  sprintf("%.3f", mortality_result$AUC),
  "; скорректированное по Бенджамини–Хохбергу q = ",
  format_p(mortality_result$BH_secondary),
  "). ",
  "(B) Распределение показателя в зависимости от места забора образца. ",
  "В образцах из отделения интенсивной терапии значения показателя были выше, ",
  "чем в образцах из отделения неотложной помощи (AUC = ",
  sprintf("%.3f", icu_result$AUC),
  "; q = ",
  format_p(icu_result$BH_secondary),
  "). ",
  "(C) Корреляции Спирмена между пятигенным показателем и непрерывным ",
  "24-часовым SOFA в отдельных географических подгруппах. Точки показывают ",
  "замороженные оценки корреляции, горизонтальные линии — 95% доверительные ",
  "интервалы. Все пять оценок имели положительное направление, хотя величина ",
  "эффекта и точность различались между подгруппами. Фиолетовая пунктирная ",
  "линия показывает описательную фиксированную объединенную оценку (rho = ",
  sprintf("%.3f", location_pooled$fixed_rho),
  "). Эти географические подгруппы принадлежат одному набору данных и не ",
  "рассматриваются как независимые внешние когорты. ",
  "(D) Контекстное сравнение пациентов с сепсисом и здоровых контролей с ",
  "использованием версии пятигенного показателя, стандартизированной по всем ",
  "образцам (AUC = ",
  sprintf("%.3f", context_result$AUC),
  "). Этот анализ представлен только как биологический контекст и не является ",
  "основной клинически значимой конечной точкой внешней проверки. ",
  "Для рисунка не выполнялись новые статистические проверки; все показатели ",
  "взяты из замороженных анализов Scripts 142b/143."
)


caption_RU_file <- file.path(
  text_dir,
  "166_FigureS8_caption_RU.txt"
)

writeLines(
  caption_RU,
  con = caption_RU_file,
  useBytes = TRUE
)


# ------------------------------------------------------------------------------
# 26. RESULTS PLACEMENT TEXT
# ------------------------------------------------------------------------------

results_placement_EN <- paste0(
  "Secondary and sensitivity analyses were directionally consistent with the ",
  "primary external severity result (Supplementary Fig. S8). The five-gene ",
  "score was higher among patients who died in hospital than among survivors ",
  "(AUC = ",
  sprintf("%.3f", mortality_result$AUC),
  ", P = ",
  format_p(mortality_result$p_value),
  ") and among ICU compared with Emergency Room samples (AUC = ",
  sprintf("%.3f", icu_result$AUC),
  ", P = ",
  format_p(icu_result$p_value),
  "). Location-stratified analyses yielded positive score-SOFA correlations ",
  "in all five collection locations, although precision varied among strata. ",
  "The five locations represent sensitivity strata within GSE185263 rather ",
  "than independent validation cohorts. Sepsis-versus-healthy discrimination ",
  "was high (AUC = ",
  sprintf("%.3f", context_result$AUC),
  ") but was considered contextual rather than a primary external validation ",
  "endpoint."
)


results_file <- file.path(
  text_dir,
  "166_Results_3.9_supplementary_sentence_EN.txt"
)

writeLines(
  results_placement_EN,
  con = results_file,
  useBytes = TRUE
)


# ------------------------------------------------------------------------------
# 27. REPORTING GUARDRAILS
# ------------------------------------------------------------------------------

guardrails <- c(
  "Supplementary Figure S8 reporting guardrails",
  "",
  "- No new inferential statistical analysis was performed.",
  "- Mortality and ICU-versus-Emergency-Room analyses are secondary.",
  "- Sepsis-versus-healthy discrimination is contextual only.",
  "- The primary external endpoint remains continuous 24-h SOFA in sepsis.",
  "- Main Figure 5 contains the primary SOFA analyses.",
  "- Figure S8 complements rather than duplicates Main Figure 5.",
  "- Collection locations are sensitivity strata within one GSE185263 dataset.",
  "- Collection locations must not be described as independent validation cohorts.",
  "- The fixed-effect pooled location estimate is descriptive only.",
  "- No random-effects location meta-analysis is required for the main claim.",
  "- GSE185263 supports replication of an organ-dysfunction-severity association.",
  "- GSE185263 does not establish calibration of a clinical diagnostic or prognostic assay."
)


guardrail_file <- file.path(
  text_dir,
  "166_FigureS8_reporting_guardrails.txt"
)

writeLines(
  guardrails,
  con = guardrail_file,
  useBytes = TRUE
)


# ------------------------------------------------------------------------------
# 28. SAVE INTERNAL AUDIT WORKBOOK
# ------------------------------------------------------------------------------

audit_wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  audit_wb,
  "Schema_audit"
)
openxlsx::writeData(
  audit_wb,
  "Schema_audit",
  schema_audit
)

openxlsx::addWorksheet(
  audit_wb,
  "Count_audit"
)
openxlsx::writeData(
  audit_wb,
  "Count_audit",
  count_audit
)

openxlsx::addWorksheet(
  audit_wb,
  "Median_audit"
)
openxlsx::writeData(
  audit_wb,
  "Median_audit",
  median_audit
)

openxlsx::addWorksheet(
  audit_wb,
  "Anchor_audit"
)
openxlsx::writeData(
  audit_wb,
  "Anchor_audit",
  anchor_audit
)

openxlsx::addWorksheet(
  audit_wb,
  "Secondary_frozen"
)
openxlsx::writeData(
  audit_wb,
  "Secondary_frozen",
  secondary
)

openxlsx::addWorksheet(
  audit_wb,
  "Location_142b"
)
openxlsx::writeData(
  audit_wb,
  "Location_142b",
  location_142b
)

openxlsx::addWorksheet(
  audit_wb,
  "Location_143_CI"
)
openxlsx::writeData(
  audit_wb,
  "Location_143_CI",
  location_143
)

openxlsx::addWorksheet(
  audit_wb,
  "Location_meta_143"
)
openxlsx::writeData(
  audit_wb,
  "Location_meta_143",
  location_meta_143
)

if (!is.null(figure5_audit)) {
  openxlsx::addWorksheet(
    audit_wb,
    "Figure5_audit"
  )
  openxlsx::writeData(
    audit_wb,
    "Figure5_audit",
    figure5_audit
  )
}

audit_file <- file.path(
  audit_dir,
  "166_INTERNAL_AUDIT_FigureS8_GSE185263.xlsx"
)

openxlsx::saveWorkbook(
  audit_wb,
  audit_file,
  overwrite = TRUE
)


# ------------------------------------------------------------------------------
# 29. SESSION INFO
# ------------------------------------------------------------------------------

session_file <- file.path(
  log_dir,
  "166_sessionInfo.txt"
)

capture.output(
  sessionInfo(),
  file = session_file
)


# ------------------------------------------------------------------------------
# 30. FINAL CONSOLE SUMMARY
# ------------------------------------------------------------------------------

cat("\n")
cat("====================================================================\n")
cat("Script 166 FINAL completed successfully.\n")
cat("====================================================================\n\n")


cat("FIGURE S8 PANELS\n")
cat("----------------\n")
cat("A. In-hospital mortality\n")
cat("B. ICU vs Emergency Room\n")
cat("C. Location-specific score-SOFA correlations\n")
cat("D. Contextual sepsis vs healthy contrast\n\n")


cat("SECONDARY FROZEN RESULTS\n")
cat("------------------------\n")

cat(
  "Mortality: n = ",
  mortality_result$n_case,
  " vs ",
  mortality_result$n_control,
  "; AUC = ",
  signif(mortality_result$AUC, 7),
  "; P = ",
  signif(mortality_result$p_value, 7),
  "; BH q = ",
  signif(mortality_result$BH_secondary, 7),
  "\n",
  sep = ""
)

cat(
  "ICU vs Emergency Room: n = ",
  icu_result$n_case,
  " vs ",
  icu_result$n_control,
  "; AUC = ",
  signif(icu_result$AUC, 7),
  "; P = ",
  signif(icu_result$p_value, 7),
  "; BH q = ",
  signif(icu_result$BH_secondary, 7),
  "\n",
  sep = ""
)

cat(
  "Contextual sepsis vs healthy: n = ",
  context_result$n_case,
  " vs ",
  context_result$n_control,
  "; AUC = ",
  signif(context_result$AUC, 7),
  "\n\n",
  sep = ""
)


cat("LOCATION SENSITIVITY\n")
cat("--------------------\n")

print(
  location_143 %>%
    dplyr::select(
      collection_location,
      n,
      rho,
      CI_low,
      CI_high,
      p_value,
      BH_location
    )
)

cat("\n")

cat(
  "All five location-specific rho estimates positive = ",
  all(location_143$rho > 0),
  "\n",
  sep = ""
)

cat(
  "Descriptive pooled rho = ",
  signif(location_pooled$fixed_rho, 7),
  "; 95% CI = ",
  signif(location_pooled$fixed_CI_low, 7),
  " to ",
  signif(location_pooled$fixed_CI_high, 7),
  "; P = ",
  signif(location_pooled$fixed_p, 7),
  "\n\n",
  sep = ""
)


cat("OUTPUT FILES\n")
cat("------------\n")

cat(
  "PNG:\n",
  normalizePath(png_file, winslash = "\\", mustWork = TRUE),
  "\n\n",
  sep = ""
)

cat(
  "PDF:\n",
  normalizePath(pdf_file, winslash = "\\", mustWork = TRUE),
  "\n\n",
  sep = ""
)

cat(
  "TIFF:\n",
  normalizePath(tiff_file, winslash = "\\", mustWork = TRUE),
  "\n\n",
  sep = ""
)

cat(
  "Figure source data:\n",
  normalizePath(source_data_file, winslash = "\\", mustWork = TRUE),
  "\n\n",
  sep = ""
)

cat(
  "Internal audit:\n",
  normalizePath(audit_file, winslash = "\\", mustWork = TRUE),
  "\n\n",
  sep = ""
)

cat(
  "English caption:\n",
  normalizePath(caption_EN_file, winslash = "\\", mustWork = TRUE),
  "\n\n",
  sep = ""
)

cat(
  "Russian caption:\n",
  normalizePath(caption_RU_file, winslash = "\\", mustWork = TRUE),
  "\n\n",
  sep = ""
)


cat("REPORTING GUARDRAILS\n")
cat("--------------------\n")
cat("- No new inferential statistical tests.\n")
cat("- Mortality and ICU/ER analyses are secondary.\n")
cat("- Sepsis/healthy discrimination is contextual only.\n")
cat("- Locations are sensitivity strata, not independent cohorts.\n")
cat("- Pooled location rho is descriptive only.\n")
cat("- Primary external endpoint remains continuous 24-h SOFA.\n")
cat("- Figure S8 complements frozen Main Figure 5.\n\n")

cat("Done.\n")