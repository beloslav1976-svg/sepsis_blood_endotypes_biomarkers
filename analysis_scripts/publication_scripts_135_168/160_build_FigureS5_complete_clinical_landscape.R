################################################################################
# Script 160
# FINAL v3
#
# Supplementary Figure S5
#
# Complete exploratory clinical-association landscape
#
# Project:
#   Sepsis_DESeq2
#
#
# PURPOSE
# -------
#
# Visualize the complete frozen 60-test exploratory clinical-association
# family from Supplementary Table S7.
#
#
# FIGURE STRUCTURE
# ----------------
#
# A. Continuous molecular-clinical correlations
#
#       Primary five-gene score
#       SRSq
#
#    Effect metric:
#       Spearman rho
#
#
# B. Molecular scores across categorical clinical groups
#
#       Primary five-gene score
#       SRSq
#
#    Effect metric:
#       frozen median score difference
#
#
# C. Continuous clinical variables across transcriptomic endotypes
#
#       SRS
#       CTS
#
#    Because the frozen SRS and CTS analyses used different effect metrics
#    and the SRS raw effects are expressed in heterogeneous clinical units,
#    Panel C displays:
#
#       -log10(global BH q)
#
#
# D. Categorical clinical variables across transcriptomic endotypes
#
#       SRS
#       CTS
#
#    Frozen Fisher exact tests did not contain effect-size estimates.
#    Panel D therefore also displays:
#
#       -log10(global BH q)
#
#
# Panels C and D use the SAME evidence scale.
#
#
# TEST ACCOUNTING
# ---------------
#
#   Panel A = 20 tests
#   Panel B = 10 tests
#   Panel C = 20 tests
#   Panel D = 10 tests
#
#   TOTAL   = 60 tests
#
#
# THIS SCRIPT DOES NOT
# --------------------
#
#   - calculate any new association
#   - calculate any new P value
#   - recalculate BH/FDR
#   - calculate any post hoc clinical effect size
#   - refit any model
#   - alter the frozen 60-test multiplicity family
#
#
# FROZEN GLOBAL-BH SIGNIFICANT RESULTS
# ------------------------------------
#
# Primary five-gene score vs CRP:
#
#   rho = 0.5743504
#   P   = 0.0003085409
#   q   = 0.01851246
#
#
# SRSq vs CRP:
#
#   rho = 0.5260209
#   P   = 0.001172438
#   q   = 0.03517314
#
#
# INTERPRETATION
# --------------
#
# These results support association of the blood transcriptomic
# host-response axis with systemic inflammatory activity.
#
# They do NOT establish:
#
#   - prognosis
#   - causality
#   - clinical diagnostic validity
#   - treatment response
#
################################################################################


cat("====================================================================\n")
cat("Running Script 160 FINAL v3\n")
cat("Supplementary Figure S5\n")
cat("Complete exploratory clinical-association landscape\n")
cat("====================================================================\n\n")


# =============================================================================
# 1. PROJECT DIRECTORY
# =============================================================================

project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))

if (!dir.exists(project_dir)) {
  project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))
}

if (!dir.exists(project_dir)) {
  stop("Sepsis_DESeq2 project directory not found.")
}

setwd(project_dir)


cat("Project directory:\n")

print(
  normalizePath(
    getwd(),
    winslash = "\\",
    mustWork = TRUE
  )
)


# =============================================================================
# 2. PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr",
  "stringr",
  "ggplot2",
  "readxl",
  "openxlsx",
  "gridExtra"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "Missing required package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(gridExtra)
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "160_FigureS5_complete_clinical_landscape"
)


figures_dir <- file.path(
  output_dir,
  "figures"
)


tables_dir <- file.path(
  output_dir,
  "tables"
)


text_dir <- file.path(
  output_dir,
  "text"
)


audit_dir <- file.path(
  output_dir,
  "audit"
)


for (
  one_dir in c(
    output_dir,
    figures_dir,
    tables_dir,
    text_dir,
    audit_dir
  )
) {
  
  dir.create(
    one_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


cat("\nOutput folder:\n")

cat(
  normalizePath(
    output_dir,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


# =============================================================================
# 4. INPUT — FROZEN TABLE S7
# =============================================================================

tableS7_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "159_TableS7_complete_clinical_associations",
  "tables",
  "159_TableS7_complete_exploratory_clinical_associations.xlsx"
)


if (!file.exists(tableS7_file)) {
  
  stop(
    paste0(
      "Frozen Supplementary Table S7 not found:\n",
      tableS7_file
    )
  )
}


cat("\nFrozen Table S7:\n")

cat(
  normalizePath(
    tableS7_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


s7_sheets <- readxl::excel_sheets(
  tableS7_file
)


cat("\nTable S7 sheets:\n")

print(
  s7_sheets
)


if (!("Complete_60_tests" %in% s7_sheets)) {
  
  stop(
    "Sheet Complete_60_tests not found in frozen Table S7."
  )
}


clinical <- readxl::read_excel(
  tableS7_file,
  sheet = "Complete_60_tests"
) %>%
  
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat("\nFrozen clinical table dimensions:\n")

cat(
  nrow(clinical),
  " x ",
  ncol(clinical),
  "\n",
  sep = ""
)


# =============================================================================
# 5. REQUIRED COLUMN AUDIT
# =============================================================================

required_columns <- c(
  "framework",
  "clinical_variable",
  "clinical_label",
  "test_family",
  "test",
  "n",
  "statistic",
  "effect",
  "effect_name",
  "p_value",
  "group_summary",
  "note",
  "test_id",
  "BH_global",
  "BH_within_framework",
  "BH_within_test_family",
  "significance_global",
  "Global_BH_significant"
)


missing_columns <- setdiff(
  required_columns,
  names(clinical)
)


if (length(missing_columns) > 0) {
  
  cat("\nObserved columns:\n")
  
  print(
    names(clinical)
  )
  
  
  stop(
    paste0(
      "Missing required Table S7 column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


if (nrow(clinical) != 60) {
  
  stop(
    paste0(
      "Expected exactly 60 frozen tests; observed ",
      nrow(clinical),
      "."
    )
  )
}


# =============================================================================
# 6. STANDARDIZE NUMERIC FIELDS
# =============================================================================

clinical <- clinical %>%
  
  dplyr::mutate(
    
    n =
      suppressWarnings(
        as.numeric(n)
      ),
    
    statistic =
      suppressWarnings(
        as.numeric(statistic)
      ),
    
    effect =
      suppressWarnings(
        as.numeric(effect)
      ),
    
    p_value =
      suppressWarnings(
        as.numeric(p_value)
      ),
    
    BH_global =
      suppressWarnings(
        as.numeric(BH_global)
      ),
    
    BH_within_framework =
      suppressWarnings(
        as.numeric(BH_within_framework)
      ),
    
    BH_within_test_family =
      suppressWarnings(
        as.numeric(BH_within_test_family)
      ),
    
    Global_BH_significant =
      as.logical(
        Global_BH_significant
      )
  )


# =============================================================================
# 7. EFFECT-SIZE AVAILABILITY AUDIT
# =============================================================================

missing_effect <- is.na(
  clinical$effect
)


missing_effect_table <- clinical[
  missing_effect,
  c(
    "framework",
    "clinical_variable",
    "clinical_label",
    "test_family",
    "test",
    "n",
    "effect",
    "effect_name",
    "p_value",
    "BH_global"
  ),
  drop = FALSE
]


cat("\nFROZEN EFFECT-SIZE AVAILABILITY\n")
cat("-------------------------------\n")


cat(
  "Tests with defined frozen effect = ",
  sum(
    !missing_effect
  ),
  "\n",
  sep = ""
)


cat(
  "Tests without frozen effect = ",
  sum(
    missing_effect
  ),
  "\n",
  sep = ""
)


if (sum(missing_effect) != 10) {
  
  print(
    missing_effect_table,
    row.names = FALSE
  )
  
  stop(
    paste0(
      "Expected exactly 10 frozen tests without effect size; observed ",
      sum(missing_effect),
      "."
    )
  )
}


if (
  !all(
    clinical$test_family[
      missing_effect
    ] ==
    "categorical_endotype"
  )
) {
  
  stop(
    "Unexpected missing effect outside categorical_endotype tests."
  )
}


if (
  !all(
    clinical$test[
      missing_effect
    ] ==
    "Fisher_exact"
  )
) {
  
  stop(
    "Frozen missing-effect rows are not exclusively Fisher exact tests."
  )
}


cat(
  "Missing-effect audit passed: 10/10 are categorical-endotype Fisher exact tests.\n"
)


# =============================================================================
# 8. P VALUE / BH AUDIT
# =============================================================================

if (
  any(
    !is.finite(
      clinical$p_value
    )
  )
) {
  
  stop(
    "At least one frozen test has a non-finite nominal P value."
  )
}


if (
  any(
    !is.finite(
      clinical$BH_global
    )
  )
) {
  
  stop(
    "At least one frozen test has a non-finite global BH value."
  )
}


if (
  any(
    clinical$p_value < 0 |
    clinical$p_value > 1
  )
) {
  
  stop(
    "At least one nominal P value is outside [0,1]."
  )
}


if (
  any(
    clinical$BH_global < 0 |
    clinical$BH_global > 1
  )
) {
  
  stop(
    "At least one global BH value is outside [0,1]."
  )
}


# =============================================================================
# 9. GLOBAL MULTIPLICITY AUDIT
# =============================================================================

if (
  sum(
    clinical$Global_BH_significant
  ) !=
  2
) {
  
  stop(
    paste0(
      "Expected exactly two globally significant tests; observed ",
      sum(
        clinical$Global_BH_significant
      ),
      "."
    )
  )
}


if (
  !all(
    clinical$Global_BH_significant ==
    (
      clinical$BH_global <
      0.05
    )
  )
) {
  
  stop(
    "Stored global-significance flag is inconsistent with BH_global."
  )
}


cat("\nGlobal frozen audit passed:\n")

cat(
  "60 tests; 2 global-BH significant.\n"
)


# =============================================================================
# 10. TEST-FAMILY STRUCTURE
# =============================================================================

family_counts <- clinical %>%
  
  dplyr::count(
    framework,
    test_family,
    name = "n_tests"
  ) %>%
  
  dplyr::arrange(
    framework,
    test_family
  )


cat("\nTEST-FAMILY STRUCTURE\n")
cat("---------------------\n")


print(
  family_counts,
  row.names = FALSE
)


expected_family_counts <- data.frame(
  
  framework = c(
    "CTS_class",
    "CTS_class",
    "Primary_5gene_score",
    "Primary_5gene_score",
    "SRS_class",
    "SRS_class",
    "SRSq",
    "SRSq"
  ),
  
  test_family = c(
    "categorical_endotype",
    "continuous_by_endotype",
    "continuous_correlation",
    "molecular_score_by_clinical_group",
    "categorical_endotype",
    "continuous_by_endotype",
    "continuous_correlation",
    "molecular_score_by_clinical_group"
  ),
  
  n_tests = c(
    5,
    10,
    10,
    5,
    5,
    10,
    10,
    5
  ),
  
  stringsAsFactors = FALSE
)


family_audit <- dplyr::full_join(
  expected_family_counts,
  family_counts,
  by = c(
    "framework",
    "test_family"
  ),
  suffix = c(
    "_expected",
    "_observed"
  )
) %>%
  
  dplyr::mutate(
    match =
      n_tests_expected ==
      n_tests_observed
  )


if (
  any(
    is.na(
      family_audit$match
    )
  ) ||
  !all(
    family_audit$match
  )
) {
  
  print(
    family_audit,
    row.names = FALSE
  )
  
  stop(
    "Frozen 60-test family structure does not match expected architecture."
  )
}


cat(
  "Frozen 20 + 10 + 20 + 10 test-family structure audit passed.\n"
)


# =============================================================================
# 11. EFFECT METRIC AUDIT
# =============================================================================

effect_name_summary <- clinical %>%
  
  dplyr::mutate(
    
    effect_name_display =
      ifelse(
        is.na(
          effect_name
        ) |
          trimws(
            as.character(
              effect_name
            )
          ) ==
          "",
        "<not defined>",
        as.character(
          effect_name
        )
      )
  ) %>%
  
  dplyr::distinct(
    framework,
    test_family,
    effect_name_display
  ) %>%
  
  dplyr::arrange(
    test_family,
    framework,
    effect_name_display
  )


cat("\nEFFECT METRICS BY TEST FAMILY\n")
cat("-----------------------------\n")


print(
  effect_name_summary,
  row.names = FALSE
)


# =============================================================================
# 12. DISPLAY VARIABLES
# =============================================================================

clinical <- clinical %>%
  
  dplyr::mutate(
    
    framework_display =
      dplyr::case_when(
        
        framework ==
          "Primary_5gene_score" ~
          "Five-gene score",
        
        framework ==
          "SRSq" ~
          "SRSq",
        
        framework ==
          "SRS_class" ~
          "SRS",
        
        framework ==
          "CTS_class" ~
          "CTS",
        
        TRUE ~
          framework
      ),
    
    clinical_display =
      as.character(
        clinical_label
      ),
    
    significance_label =
      factor(
        ifelse(
          Global_BH_significant,
          "q < 0.05",
          "q \u2265 0.05"
        ),
        levels = c(
          "q \u2265 0.05",
          "q < 0.05"
        )
      ),
    
    minus_log10_global_BH =
      -log10(
        BH_global
      )
  )


# =============================================================================
# 13. PANEL A DATA
# =============================================================================

panel_A_data <- clinical %>%
  
  dplyr::filter(
    test_family ==
      "continuous_correlation"
  )


if (nrow(panel_A_data) != 20) {
  
  stop(
    paste0(
      "Panel A expected 20 tests; observed ",
      nrow(panel_A_data),
      "."
    )
  )
}


if (
  any(
    !is.finite(
      panel_A_data$effect
    )
  )
) {
  
  stop(
    "Panel A contains a non-finite frozen Spearman effect."
  )
}


if (
  !all(
    grepl(
      "rho|spearman",
      panel_A_data$effect_name,
      ignore.case = TRUE
    )
  )
) {
  
  stop(
    "Panel A contains an unexpected effect metric."
  )
}


panel_A_levels <- unique(
  panel_A_data$clinical_display
)


panel_A_data$clinical_display <- factor(
  panel_A_data$clinical_display,
  levels = rev(
    panel_A_levels
  )
)


panel_A_data$framework_display <- factor(
  panel_A_data$framework_display,
  levels = c(
    "Five-gene score",
    "SRSq"
  )
)


# =============================================================================
# 14. PANEL B DATA
# =============================================================================

panel_B_data <- clinical %>%
  
  dplyr::filter(
    test_family ==
      "molecular_score_by_clinical_group"
  )


if (nrow(panel_B_data) != 10) {
  
  stop(
    paste0(
      "Panel B expected 10 tests; observed ",
      nrow(panel_B_data),
      "."
    )
  )
}


if (
  any(
    !is.finite(
      panel_B_data$effect
    )
  )
) {
  
  stop(
    "Panel B contains an unexpected non-finite frozen effect."
  )
}


panel_B_levels <- unique(
  panel_B_data$clinical_display
)


panel_B_data$clinical_display <- factor(
  panel_B_data$clinical_display,
  levels = rev(
    panel_B_levels
  )
)


panel_B_data$framework_display <- factor(
  panel_B_data$framework_display,
  levels = c(
    "Five-gene score",
    "SRSq"
  )
)


# =============================================================================
# 15. PANEL C DATA
# =============================================================================

panel_C_data <- clinical %>%
  
  dplyr::filter(
    test_family ==
      "continuous_by_endotype"
  )


if (nrow(panel_C_data) != 20) {
  
  stop(
    paste0(
      "Panel C expected 20 tests; observed ",
      nrow(panel_C_data),
      "."
    )
  )
}


if (
  any(
    !is.finite(
      panel_C_data$minus_log10_global_BH
    )
  )
) {
  
  stop(
    "Panel C contains a non-finite transformed global BH value."
  )
}


panel_C_levels <- unique(
  panel_C_data$clinical_display
)


panel_C_data$clinical_display <- factor(
  panel_C_data$clinical_display,
  levels = rev(
    panel_C_levels
  )
)


panel_C_data$framework_display <- factor(
  panel_C_data$framework_display,
  levels = c(
    "SRS",
    "CTS"
  )
)


# =============================================================================
# 16. PANEL D DATA
# =============================================================================

panel_D_data <- clinical %>%
  
  dplyr::filter(
    test_family ==
      "categorical_endotype"
  )


if (nrow(panel_D_data) != 10) {
  
  stop(
    paste0(
      "Panel D expected 10 tests; observed ",
      nrow(panel_D_data),
      "."
    )
  )
}


if (
  !all(
    panel_D_data$test ==
    "Fisher_exact"
  )
) {
  
  stop(
    "Panel D contains a test other than Fisher exact."
  )
}


if (
  any(
    !is.na(
      panel_D_data$effect
    )
  )
) {
  
  stop(
    "Panel D unexpectedly contains a frozen effect-size estimate."
  )
}


if (
  any(
    !is.finite(
      panel_D_data$minus_log10_global_BH
    )
  )
) {
  
  stop(
    "Panel D contains a non-finite transformed global BH value."
  )
}


panel_D_levels <- unique(
  panel_D_data$clinical_display
)


panel_D_data$clinical_display <- factor(
  panel_D_data$clinical_display,
  levels = rev(
    panel_D_levels
  )
)


panel_D_data$framework_display <- factor(
  panel_D_data$framework_display,
  levels = c(
    "SRS",
    "CTS"
  )
)


# =============================================================================
# 17. SHARED C/D EVIDENCE SCALE
# =============================================================================

global_bh_threshold <- -log10(
  0.05
)


max_CD <- max(
  c(
    panel_C_data$minus_log10_global_BH,
    panel_D_data$minus_log10_global_BH
  ),
  na.rm = TRUE
)


evidence_x_max <- max(
  1.45,
  ceiling(
    max_CD *
      10
  ) /
    10 +
    0.05
)


cat("\nPANELS C/D EVIDENCE SCALE\n")
cat("-------------------------\n")


cat(
  "Global BH threshold (-log10 0.05) = ",
  global_bh_threshold,
  "\n",
  sep = ""
)


cat(
  "Shared x-axis maximum = ",
  evidence_x_max,
  "\n",
  sep = ""
)


# =============================================================================
# 18. COLORS
# =============================================================================

framework_colors <- c(
  
  "Five-gene score" =
    "#7A0177",
  
  "SRSq" =
    "#2C7FB8",
  
  "SRS" =
    "#41AB5D",
  
  "CTS" =
    "#D95F0E"
)


# =============================================================================
# 19. COMMON THEME
# =============================================================================

clinical_theme <- ggplot2::theme_bw(
  base_size = 10
) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 9.3
      ),
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 9.5
      ),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    panel.grid.major.y =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom"
  )


# =============================================================================
# 20. PANEL A
# =============================================================================

panel_A <- ggplot2::ggplot(
  
  panel_A_data,
  
  ggplot2::aes(
    x = effect,
    y = clinical_display,
    color = framework_display,
    shape = significance_label
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#777777"
  ) +
  
  ggplot2::geom_point(
    size = 3.3,
    alpha = 0.95
  ) +
  
  ggplot2::facet_wrap(
    ~ framework_display,
    nrow = 1
  ) +
  
  ggplot2::scale_color_manual(
    values = framework_colors,
    drop = FALSE
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "q \u2265 0.05" = 16,
      "q < 0.05" = 8
    ),
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      -1,
      1
    ),
    breaks = seq(
      -1,
      1,
      by = 0.25
    )
  ) +
  
  ggplot2::labs(
    title = "A  Continuous clinical correlations",
    subtitle = "Spearman correlations; star symbols denote global BH q < 0.05",
    x = "Spearman rho",
    y = NULL,
    color = NULL,
    shape = "Global BH"
  ) +
  
  ggplot2::guides(
    
    color =
      ggplot2::guide_legend(
        order = 1
      ),
    
    shape =
      ggplot2::guide_legend(
        order = 2
      )
  ) +
  
  clinical_theme


# =============================================================================
# 21. PANEL B
# =============================================================================

panel_B <- ggplot2::ggplot(
  
  panel_B_data,
  
  ggplot2::aes(
    x = effect,
    y = clinical_display,
    color = framework_display
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#777777"
  ) +
  
  ggplot2::geom_point(
    size = 3.2,
    alpha = 0.95
  ) +
  
  ggplot2::facet_wrap(
    ~ framework_display,
    nrow = 1,
    scales = "free_x"
  ) +
  
  ggplot2::scale_color_manual(
    values = framework_colors,
    drop = FALSE
  ) +
  
  ggplot2::labs(
    title = "B  Molecular scores across categorical clinical groups",
    subtitle = "Frozen effect = median score difference between prespecified groups",
    x = "Median score difference",
    y = NULL,
    color = NULL
  ) +
  
  clinical_theme


# =============================================================================
# 22. PANEL C
# =============================================================================

panel_C <- ggplot2::ggplot(
  
  panel_C_data,
  
  ggplot2::aes(
    x = minus_log10_global_BH,
    y = clinical_display,
    color = framework_display
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = global_bh_threshold,
    linetype = "dashed",
    linewidth = 0.55,
    color = "#444444"
  ) +
  
  ggplot2::geom_point(
    size = 3.2,
    alpha = 0.95
  ) +
  
  ggplot2::facet_wrap(
    ~ framework_display,
    nrow = 1
  ) +
  
  ggplot2::scale_color_manual(
    values = framework_colors,
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      evidence_x_max
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  ggplot2::labs(
    title = "C  Continuous clinical variables across transcriptomic endotypes",
    subtitle = "Multiplicity-adjusted evidence; dashed line denotes global BH q = 0.05",
    x = expression(-log[10]("global BH q")),
    y = NULL,
    color = NULL
  ) +
  
  clinical_theme


# =============================================================================
# 23. PANEL D
# =============================================================================

panel_D <- ggplot2::ggplot(
  
  panel_D_data,
  
  ggplot2::aes(
    x = minus_log10_global_BH,
    y = clinical_display,
    color = framework_display
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = global_bh_threshold,
    linetype = "dashed",
    linewidth = 0.55,
    color = "#444444"
  ) +
  
  ggplot2::geom_point(
    size = 3.2,
    alpha = 0.95
  ) +
  
  ggplot2::facet_wrap(
    ~ framework_display,
    nrow = 1
  ) +
  
  ggplot2::scale_color_manual(
    values = framework_colors,
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      evidence_x_max
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  ggplot2::labs(
    title = "D  Categorical clinical variables across transcriptomic endotypes",
    subtitle = "Fisher exact tests; dashed line denotes global BH q = 0.05",
    x = expression(-log[10]("global BH q")),
    y = NULL,
    color = NULL
  ) +
  
  clinical_theme


# =============================================================================
# 24. COMBINE FIGURE
# =============================================================================

figure_S5_grob <- gridExtra::arrangeGrob(
  
  panel_A,
  panel_B,
  panel_C,
  panel_D,
  
  ncol = 1,
  
  heights = c(
    1.25,
    0.9,
    1.18,
    0.88
  )
)


# =============================================================================
# 25. SAVE COMBINED FIGURE
# =============================================================================

figure_pdf <- file.path(
  figures_dir,
  "160_FigureS5_complete_clinical_association_landscape.pdf"
)


figure_png <- file.path(
  figures_dir,
  "160_FigureS5_complete_clinical_association_landscape.png"
)


figure_tiff <- file.path(
  figures_dir,
  "160_FigureS5_complete_clinical_association_landscape.tiff"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure_S5_grob,
  width = 13,
  height = 20,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure_S5_grob,
  width = 13,
  height = 20,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure_S5_grob,
  width = 13,
  height = 20,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 26. SAVE INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  
  A_continuous_correlations =
    panel_A,
  
  B_scores_by_clinical_group =
    panel_B,
  
  C_continuous_by_endotype =
    panel_C,
  
  D_categorical_by_endotype =
    panel_D
)


individual_heights <- c(
  7.2,
  5.5,
  6.8,
  5.5
)


for (
  i in seq_along(
    individual_panels
  )
) {
  
  panel_name <- names(
    individual_panels
  )[i]
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      figures_dir,
      paste0(
        "160_FigureS5_",
        panel_name,
        ".pdf"
      )
    ),
    
    plot =
      individual_panels[[panel_name]],
    
    width =
      11,
    
    height =
      individual_heights[i],
    
    units =
      "in",
    
    device =
      grDevices::cairo_pdf,
    
    bg =
      "white"
  )
}


# =============================================================================
# 27. SOURCE DATA
# =============================================================================

source_data_file <- file.path(
  tables_dir,
  "160_FigureS5_source_data.xlsx"
)


wb <- openxlsx::createWorkbook()


source_objects <- list(
  
  FigureS5A =
    panel_A_data,
  
  FigureS5B =
    panel_B_data,
  
  FigureS5C =
    panel_C_data,
  
  FigureS5D =
    panel_D_data,
  
  Missing_effect_audit =
    missing_effect_table,
  
  Family_structure =
    family_counts,
  
  Effect_metrics =
    effect_name_summary
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE
)


for (
  sheet_name in names(
    source_objects
  )
) {
  
  data_object <- source_objects[[sheet_name]]
  
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  
  openxlsx::writeData(
    wb,
    sheet_name,
    data_object,
    withFilter = TRUE
  )
  
  
  openxlsx::addStyle(
    wb,
    sheet_name,
    header_style,
    rows = 1,
    cols = seq_len(
      ncol(
        data_object
      )
    ),
    gridExpand = TRUE
  )
  
  
  openxlsx::freezePane(
    wb,
    sheet_name,
    firstActiveRow = 2
  )
  
  
  openxlsx::setColWidths(
    wb,
    sheet_name,
    cols = seq_len(
      ncol(
        data_object
      )
    ),
    widths = "auto"
  )
}


openxlsx::saveWorkbook(
  wb,
  source_data_file,
  overwrite = TRUE
)


# =============================================================================
# 28. GLOBAL-BH SIGNIFICANT ROWS
# =============================================================================

global_sig <- clinical %>%
  
  dplyr::filter(
    Global_BH_significant
  )


if (nrow(global_sig) != 2) {
  
  stop(
    "Expected exactly two globally significant rows."
  )
}


primary_crp <- global_sig %>%
  
  dplyr::filter(
    framework ==
      "Primary_5gene_score",
    clinical_variable ==
      "crp_numeric"
  )


srsq_crp <- global_sig %>%
  
  dplyr::filter(
    framework ==
      "SRSq",
    clinical_variable ==
      "crp_numeric"
  )


if (
  nrow(primary_crp) != 1 ||
  nrow(srsq_crp) != 1
) {
  
  stop(
    "Could not uniquely recover both frozen CRP anchors."
  )
}


# =============================================================================
# 29. FIGURE LEGEND
# =============================================================================

figure_legend <- paste0(
  
  "Supplementary Figure S5. Complete exploratory clinical-association ",
  "landscape for blood transcriptomic host-response measures. ",
  
  "(A) Spearman correlations of the primary five-gene host-response score ",
  "and continuous SRSq output with continuous clinical variables. Star-shaped ",
  "symbols denote associations remaining significant after Benjamini-Hochberg ",
  "correction across the complete 60-test exploratory family. ",
  
  "(B) Frozen median differences in the primary five-gene score and SRSq ",
  "between prespecified categorical clinical groups. ",
  
  "(C) Multiplicity-adjusted evidence for associations of continuous clinical ",
  "variables with SRS and CTS transcriptomic classes. Because the corresponding ",
  "frozen analyses used heterogeneous effect metrics and clinical units, the ",
  "panel displays -log10 of the global Benjamini-Hochberg-adjusted P value ",
  "rather than raw effect magnitude. ",
  
  "(D) Multiplicity-adjusted evidence for Fisher exact tests relating ",
  "categorical clinical variables to SRS and CTS classes. No effect-size ",
  "statistic was calculated for these tests in the frozen analysis; therefore, ",
  "the panel likewise displays -log10 of the global adjusted P value. ",
  
  "Panels C and D use the same horizontal evidence scale, with the dashed line ",
  "indicating global BH q=0.05. Together, panels A-D represent all 60 evaluable ",
  "tests in the exploratory multiplicity family. Only two associations remained ",
  "significant after global correction: the primary five-gene score versus ",
  "C-reactive protein (Spearman rho=",
  sprintf(
    "%.3f",
    primary_crp$effect
  ),
  ", q=",
  sprintf(
    "%.3f",
    primary_crp$BH_global
  ),
  ") and SRSq versus C-reactive protein (rho=",
  sprintf(
    "%.3f",
    srsq_crp$effect
  ),
  ", q=",
  sprintf(
    "%.3f",
    srsq_crp$BH_global
  ),
  "). No statistical association, effect estimate, or multiple-testing ",
  "adjustment was recalculated for figure generation. These analyses provide ",
  "exploratory molecular-clinical context rather than independent diagnostic ",
  "or prognostic validation."
)


legend_file <- file.path(
  text_dir,
  "160_FigureS5_legend_EN.txt"
)


writeLines(
  figure_legend,
  legend_file
)


# =============================================================================
# 30. PROPOSED RESULTS 3.6
# =============================================================================

results_3_6 <- paste0(
  
  "We next examined relationships between transcriptomic host-response ",
  "measures and available clinical characteristics across a complete ",
  "exploratory family of 60 tests (Supplementary Fig. S5 and Supplementary ",
  "Table S7). After Benjamini-Hochberg correction across the entire test ",
  "family, only two associations remained significant. The primary five-gene ",
  "score correlated positively with C-reactive protein (Spearman rho=",
  sprintf(
    "%.3f",
    primary_crp$effect
  ),
  ", P=",
  format(
    primary_crp$p_value,
    scientific = TRUE,
    digits = 3
  ),
  ", q=",
  sprintf(
    "%.3f",
    primary_crp$BH_global
  ),
  "), and SRSq showed a concordant positive association with C-reactive ",
  "protein (rho=",
  sprintf(
    "%.3f",
    srsq_crp$effect
  ),
  ", P=",
  format(
    srsq_crp$p_value,
    scientific = TRUE,
    digits = 3
  ),
  ", q=",
  sprintf(
    "%.3f",
    srsq_crp$BH_global
  ),
  "). No other molecular-clinical association remained significant after ",
  "global correction. The complete association landscape showed no additional ",
  "endotype-clinical relationship that survived correction across the full ",
  "multiplicity family. These findings therefore support a relationship ",
  "between the blood transcriptomic host-response axis and systemic ",
  "inflammatory activity, while providing no evidence of independent ",
  "prognostic validation."
)


results_file <- file.path(
  text_dir,
  "160_proposed_Results_3.6_clinical_associations_EN.txt"
)


writeLines(
  results_3_6,
  results_file
)


# =============================================================================
# 31. AUDIT WORKBOOK
# =============================================================================

audit_file <- file.path(
  audit_dir,
  "160_FigureS5_audit.xlsx"
)


audit_summary <- data.frame(
  
  Metric = c(
    "Total frozen tests",
    "Tests with frozen effect",
    "Tests without frozen effect",
    "Panel A tests",
    "Panel B tests",
    "Panel C tests",
    "Panel D tests",
    "Total represented",
    "Global-BH significant tests",
    "Primary score vs CRP rho",
    "Primary score vs CRP P",
    "Primary score vs CRP global BH",
    "SRSq vs CRP rho",
    "SRSq vs CRP P",
    "SRSq vs CRP global BH",
    "Global-BH evidence threshold",
    "Panels C/D shared x maximum"
  ),
  
  Value = c(
    nrow(clinical),
    sum(!is.na(clinical$effect)),
    sum(is.na(clinical$effect)),
    nrow(panel_A_data),
    nrow(panel_B_data),
    nrow(panel_C_data),
    nrow(panel_D_data),
    nrow(panel_A_data) +
      nrow(panel_B_data) +
      nrow(panel_C_data) +
      nrow(panel_D_data),
    sum(clinical$Global_BH_significant),
    primary_crp$effect,
    primary_crp$p_value,
    primary_crp$BH_global,
    srsq_crp$effect,
    srsq_crp$p_value,
    srsq_crp$BH_global,
    global_bh_threshold,
    evidence_x_max
  ),
  
  stringsAsFactors = FALSE
)


wb_audit <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wb_audit,
  "Audit_summary"
)


openxlsx::writeData(
  wb_audit,
  "Audit_summary",
  audit_summary
)


openxlsx::addStyle(
  wb_audit,
  "Audit_summary",
  header_style,
  rows = 1,
  cols = seq_len(
    ncol(
      audit_summary
    )
  ),
  gridExpand = TRUE
)


openxlsx::setColWidths(
  wb_audit,
  "Audit_summary",
  cols = seq_len(
    ncol(
      audit_summary
    )
  ),
  widths = "auto"
)


openxlsx::addWorksheet(
  wb_audit,
  "Missing_effect_tests"
)


openxlsx::writeData(
  wb_audit,
  "Missing_effect_tests",
  missing_effect_table
)


openxlsx::addStyle(
  wb_audit,
  "Missing_effect_tests",
  header_style,
  rows = 1,
  cols = seq_len(
    ncol(
      missing_effect_table
    )
  ),
  gridExpand = TRUE
)


openxlsx::setColWidths(
  wb_audit,
  "Missing_effect_tests",
  cols = seq_len(
    ncol(
      missing_effect_table
    )
  ),
  widths = "auto"
)


openxlsx::saveWorkbook(
  wb_audit,
  audit_file,
  overwrite = TRUE
)


# =============================================================================
# 32. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "160_sessionInfo.txt"
  )
)


# =============================================================================
# 33. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 160 FINAL v3 completed successfully.\n")
cat("====================================================================\n\n")


cat("COMPLETE CLINICAL LANDSCAPE\n")
cat("---------------------------\n")


cat(
  "Total frozen tests = ",
  nrow(clinical),
  "\n",
  sep = ""
)


cat(
  "Frozen tests with defined effect = ",
  sum(
    !is.na(
      clinical$effect
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Frozen tests without effect = ",
  sum(
    is.na(
      clinical$effect
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Panel A = ",
  nrow(panel_A_data),
  "\n",
  sep = ""
)


cat(
  "Panel B = ",
  nrow(panel_B_data),
  "\n",
  sep = ""
)


cat(
  "Panel C = ",
  nrow(panel_C_data),
  "\n",
  sep = ""
)


cat(
  "Panel D = ",
  nrow(panel_D_data),
  "\n",
  sep = ""
)


cat(
  "Total represented = ",
  nrow(panel_A_data) +
    nrow(panel_B_data) +
    nrow(panel_C_data) +
    nrow(panel_D_data),
  "\n",
  sep = ""
)


cat("\nGLOBAL-BH SIGNIFICANT RESULTS\n")
cat("-----------------------------\n")


print(
  global_sig[
    ,
    c(
      "framework",
      "clinical_label",
      "effect",
      "effect_name",
      "p_value",
      "BH_global"
    )
  ],
  row.names = FALSE
)


cat("\nPANELS C/D GLOBAL-BH EVIDENCE RANGE\n")
cat("-----------------------------------\n")


cat(
  "q=0.05 threshold = ",
  global_bh_threshold,
  "\n",
  sep = ""
)


cat(
  "Panel C max -log10(q) = ",
  max(
    panel_C_data$minus_log10_global_BH
  ),
  "\n",
  sep = ""
)


cat(
  "Panel D max -log10(q) = ",
  max(
    panel_D_data$minus_log10_global_BH
  ),
  "\n",
  sep = ""
)


cat(
  "Shared C/D x-axis max = ",
  evidence_x_max,
  "\n",
  sep = ""
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Figure S5 PDF:\n  ",
  normalizePath(
    figure_pdf,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S5 PNG:\n  ",
  normalizePath(
    figure_png,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S5 TIFF:\n  ",
  normalizePath(
    figure_tiff,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure source data:\n  ",
  normalizePath(
    source_data_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure legend:\n  ",
  normalizePath(
    legend_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Proposed Results 3.6 text:\n  ",
  normalizePath(
    results_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")


cat(
  "- Figure S5 represents all 60 frozen clinical-association tests.\n"
)


cat(
  "- No association, P value, effect estimate, or BH/FDR is recalculated.\n"
)


cat(
  "- Panel A retains directly comparable Spearman rho effect sizes.\n"
)


cat(
  "- Panel B retains within-score median differences.\n"
)


cat(
  "- Panel C does not compare heterogeneous raw clinical effect units.\n"
)


cat(
  "- Panels C and D display the same -log10(global BH q) evidence scale.\n"
)


cat(
  "- Dashed lines in C/D indicate global BH q=0.05.\n"
)


cat(
  "- Only Primary score vs CRP and SRSq vs CRP are globally significant.\n"
)


cat(
  "- Clinical associations remain exploratory.\n"
)


cat(
  "- Do not interpret Figure S5 as prognostic validation.\n"
)


cat("\nDone.\n")