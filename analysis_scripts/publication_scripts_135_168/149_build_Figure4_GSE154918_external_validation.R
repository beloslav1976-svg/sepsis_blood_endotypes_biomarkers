################################################################################
# Script 149
# FINAL Main Figure 4
# External evaluation of the frozen five-gene host-response signature
# Dataset: GSE154918
#
# Project: Sepsis_DESeq2
#
#
# FIGURE STRUCTURE
# ----------------
# A. Frozen score across four baseline clinical states
# B. Prespecified primary external comparison:
#       Sepsis + septic shock versus uncomplicated infection
# C. ROC curve for the prespecified primary external comparison
# D. Secondary score-level external contrasts
#
#
# PURPOSE
# -------
# Publication packaging of FINALIZED Script 141 results.
#
#
# IMPORTANT
# ---------
# This script does NOT:
#   - download or reprocess GEO expression data;
#   - recalculate the five-gene score;
#   - re-standardize expression;
#   - perform feature selection;
#   - alter the frozen five-gene composition;
#   - refit coefficients;
#   - optimize a cutoff;
#   - reverse score direction;
#   - rerun Wilcoxon tests;
#   - rerun BH correction;
#   - rerun ordered-status tests;
#   - rerun component-gene inference;
#   - create new hypothesis tests.
#
# For Panel C only, ROC coordinates are generated from the already frozen
# sample-level score table for visualization. The resulting AUC is audited
# against the frozen AUC stored in Script 141 output.
#
#
# FROZEN PANEL
# ------------
# UP:
#   CD177
#   HK3
#   IRAK3
#
# DOWN:
#   CARD11
#   IKZF2
#
#
# FROZEN SCORE
# ------------
# mean[z(CD177), z(HK3), z(IRAK3)] -
# mean[z(CARD11), z(IKZF2)]
#
#
# PRIMARY EXTERNAL COMPARISON
# ---------------------------
# Baseline Seps_P + Shock_P versus Inf1_P
#
# n:
#   case = 39
#   control = 12
#
# median:
#   case = 1.766097488
#   control = 0.998897965
#
# Wilcoxon:
#   P = 0.1074155362
#
# BH across five score-level comparisons:
#   q = 0.1404956452
#
# fixed-direction AUC:
#   0.655982906
#
# 95% CI:
#   0.505583980 - 0.806381832
#
#
# BASELINE GROUPS
# ---------------
# Hlty    = 40
# Inf1_P  = 12
# Seps_P  = 20
# Shock_P = 19
#
# Total baseline n = 91
#
#
# ORDERED BASELINE ANALYSIS
# -------------------------
# Hlty -> Inf1_P -> Seps_P -> Shock_P
#
# Spearman rho = 0.8418590062
# P = 1.4462160066e-25
#
# Kruskal-Wallis:
# chi-square = 68.45317977
# P = 9.150159001e-15
#
#
# KEY SECONDARY RESULT
# --------------------
# Shock_P versus Inf1_P
#
# median:
#   shock = 1.952711539
#   uncomplicated infection = 0.998897965
#
# P = 0.00890133605
# BH q = 0.02225334012
#
# AUC = 0.7850877193
# 95% CI = 0.6177149010 - 0.9524605376
#
#
# INTERPRETATION
# --------------
# External transcriptomic evaluation/replication of a frozen molecular score.
#
# NOT:
#   - validation of a calibrated clinical diagnostic assay;
#   - validation of a clinical cutoff;
#   - proof of diagnostic superiority.
#
#
# OUTPUT
# ------
# results/blood_endotypes_biomarkers/
#   149_Figure4_GSE154918_external_validation/
#
################################################################################


cat("====================================================================\n")
cat("Running Script 149\n")
cat("FINAL Main Figure 4: GSE154918 external evaluation\n")
cat("====================================================================\n\n")


# =============================================================================
# 1. PROJECT DIRECTORY
# =============================================================================

project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))

if (!dir.exists(project_dir)) {
  project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))
}

if (!dir.exists(project_dir)) {
  stop(
    "Sepsis_DESeq2 project directory not found."
  )
}

setwd(
  project_dir
)


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
  "tidyr",
  "ggplot2",
  "patchwork",
  "openxlsx",
  "pROC",
  "stringr",
  "forcats",
  "scales"
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
  library(tidyr)
  
  library(ggplot2)
  library(patchwork)
  
  library(openxlsx)
  library(pROC)
  
  library(stringr)
  library(forcats)
  library(scales)
  
})


# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================


find_project_file <- function(
    candidates,
    recursive_pattern = NULL,
    description = "file"
) {
  
  candidates <- unique(
    candidates
  )
  
  
  existing <- candidates[
    file.exists(
      candidates
    )
  ]
  
  
  if (length(existing) > 0) {
    
    return(
      existing[1]
    )
  }
  
  
  if (!is.null(recursive_pattern)) {
    
    hits <- list.files(
      path = "results",
      pattern = recursive_pattern,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    
    hits <- sort(
      unique(
        hits
      )
    )
    
    
    if (length(hits) > 0) {
      
      if (length(hits) > 1) {
        
        cat(
          "\nMultiple candidate files found for ",
          description,
          ":\n",
          sep = ""
        )
        
        
        for (h in hits) {
          
          cat(
            "  ",
            h,
            "\n",
            sep = ""
          )
        }
        
        
        cat(
          "Using first candidate:\n  ",
          hits[1],
          "\n",
          sep = ""
        )
      }
      
      
      return(
        hits[1]
      )
    }
  }
  
  
  return(
    NA_character_
  )
}


theme_publication <- function(
    base_size = 10
) {
  
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    
    ggplot2::theme(
      
      plot.title =
        ggplot2::element_text(
          face = "bold",
          size = base_size + 1
        ),
      
      plot.subtitle =
        ggplot2::element_text(
          size = base_size - 0.5,
          color = "grey25"
        ),
      
      plot.tag =
        ggplot2::element_text(
          face = "bold",
          size = base_size + 2
        ),
      
      axis.title =
        ggplot2::element_text(
          face = "bold"
        ),
      
      legend.title =
        ggplot2::element_text(
          face = "bold"
        ),
      
      plot.margin =
        ggplot2::margin(
          7,
          7,
          7,
          7
        )
    )
}


format_p <- function(x) {
  
  x <- as.numeric(
    x
  )
  
  
  if (!is.finite(x)) {
    
    return(
      "NA"
    )
  }
  
  
  if (x < 0.001) {
    
    return(
      format(
        x,
        scientific = TRUE,
        digits = 3
      )
    )
  }
  
  
  sprintf(
    "%.4f",
    x
  )
}


format_q <- function(x) {
  
  x <- as.numeric(
    x
  )
  
  
  if (!is.finite(x)) {
    
    return(
      "NA"
    )
  }
  
  
  if (x < 0.001) {
    
    return(
      format(
        x,
        scientific = TRUE,
        digits = 3
      )
    )
  }
  
  
  sprintf(
    "%.4f",
    x
  )
}


# =============================================================================
# 4. LOCATE FROZEN SCRIPT 141 WORKBOOK
# =============================================================================

source_workbook <- find_project_file(
  
  candidates = c(
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "141_external_validation_GSE154918",
      "tables",
      "141_GSE154918_external_validation.xlsx"
    ),
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "141_external_validation_GSE154918",
      "141_GSE154918_external_validation.xlsx"
    )
  ),
  
  recursive_pattern =
    "^141_GSE154918_external_validation\\.xlsx$",
  
  description =
    "Script 141 GSE154918 external-validation workbook"
)


if (
  is.na(
    source_workbook
  ) ||
  !file.exists(
    source_workbook
  )
) {
  
  stop(
    paste0(
      "Could not locate frozen Script 141 workbook:\n",
      "141_GSE154918_external_validation.xlsx"
    )
  )
}


cat("\n====================================================================\n")
cat("SOURCE WORKBOOK\n")
cat("====================================================================\n")


cat(
  normalizePath(
    source_workbook,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


# =============================================================================
# 5. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "149_Figure4_GSE154918_external_validation"
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


logs_dir <- file.path(
  output_dir,
  "logs"
)


for (d in c(
  output_dir,
  figures_dir,
  tables_dir,
  text_dir,
  logs_dir
)) {
  
  dir.create(
    d,
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
# 6. WORKBOOK SHEET AUDIT
# =============================================================================

required_sheets <- c(
  "00_run_info",
  "04_sample_scores",
  "05_score_comparisons",
  "07_gene_direction_audit",
  "08_ordered_tests",
  "09_pairwise_tests",
  "10_score_group_summary",
  "12_scaling_sensitivity"
)


available_sheets <- openxlsx::getSheetNames(
  source_workbook
)


cat("\nAvailable Script 141 sheets:\n")

print(
  available_sheets
)


missing_sheets <- setdiff(
  required_sheets,
  available_sheets
)


if (length(missing_sheets) > 0) {
  
  stop(
    paste0(
      "Required Script 141 sheet(s) missing:\n",
      paste(
        missing_sheets,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 7. LOAD FROZEN SOURCE TABLES
# =============================================================================

run_info <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "00_run_info"
)


scores_df <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "04_sample_scores"
)


comparisons <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "05_score_comparisons"
)


gene_direction <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "07_gene_direction_audit"
)


ordered_tests <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "08_ordered_tests"
)


pairwise_tests <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "09_pairwise_tests"
)


group_summary <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "10_score_group_summary"
)


scaling_sensitivity <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "12_scaling_sensitivity"
)


# =============================================================================
# 8. REQUIRED COLUMN AUDIT
# =============================================================================

required_score_columns <- c(
  "geo_accession",
  "sample_title",
  "status",
  "five_gene_score",
  "five_gene_score_healthy_reference",
  "severity_order"
)


missing_score_columns <- setdiff(
  required_score_columns,
  names(
    scores_df
  )
)


if (length(missing_score_columns) > 0) {
  
  stop(
    paste0(
      "Missing columns in 04_sample_scores:\n",
      paste(
        missing_score_columns,
        collapse = ", "
      )
    )
  )
}


required_comparison_columns <- c(
  "comparison",
  "score",
  "case_status",
  "control_status",
  "n_case",
  "n_control",
  "case_median",
  "control_median",
  "median_difference_case_minus_control",
  "p_value",
  "auc_fixed_direction",
  "auc_ci_low",
  "auc_ci_high",
  "p_BH_across_score_comparisons"
)


missing_comparison_columns <- setdiff(
  required_comparison_columns,
  names(
    comparisons
  )
)


if (length(missing_comparison_columns) > 0) {
  
  stop(
    paste0(
      "Missing columns in 05_score_comparisons:\n",
      paste(
        missing_comparison_columns,
        collapse = ", "
      )
    )
  )
}


required_ordered_columns <- c(
  "analysis",
  "statistic",
  "statistic_name",
  "p_value"
)


missing_ordered_columns <- setdiff(
  required_ordered_columns,
  names(
    ordered_tests
  )
)


if (length(missing_ordered_columns) > 0) {
  
  stop(
    paste0(
      "Missing columns in 08_ordered_tests:\n",
      paste(
        missing_ordered_columns,
        collapse = ", "
      )
    )
  )
}


required_summary_columns <- c(
  "status",
  "n",
  "median",
  "q1",
  "q3",
  "mean",
  "sd"
)


missing_summary_columns <- setdiff(
  required_summary_columns,
  names(
    group_summary
  )
)


if (length(missing_summary_columns) > 0) {
  
  stop(
    paste0(
      "Missing columns in 10_score_group_summary:\n",
      paste(
        missing_summary_columns,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 9. NORMALIZE TYPES
# =============================================================================

scores_df <- scores_df %>%
  
  dplyr::mutate(
    
    status =
      as.character(
        status
      ),
    
    five_gene_score =
      as.numeric(
        five_gene_score
      ),
    
    five_gene_score_healthy_reference =
      as.numeric(
        five_gene_score_healthy_reference
      ),
    
    severity_order =
      as.numeric(
        severity_order
      )
  )


comparisons <- comparisons %>%
  
  dplyr::mutate(
    
    n_case =
      as.numeric(
        n_case
      ),
    
    n_control =
      as.numeric(
        n_control
      ),
    
    case_median =
      as.numeric(
        case_median
      ),
    
    control_median =
      as.numeric(
        control_median
      ),
    
    median_difference_case_minus_control =
      as.numeric(
        median_difference_case_minus_control
      ),
    
    p_value =
      as.numeric(
        p_value
      ),
    
    auc_fixed_direction =
      as.numeric(
        auc_fixed_direction
      ),
    
    auc_ci_low =
      as.numeric(
        auc_ci_low
      ),
    
    auc_ci_high =
      as.numeric(
        auc_ci_high
      ),
    
    p_BH_across_score_comparisons =
      as.numeric(
        p_BH_across_score_comparisons
      )
  )


ordered_tests <- ordered_tests %>%
  
  dplyr::mutate(
    
    statistic =
      as.numeric(
        statistic
      ),
    
    p_value =
      as.numeric(
        p_value
      )
  )


group_summary <- group_summary %>%
  
  dplyr::mutate(
    
    n =
      as.numeric(
        n
      ),
    
    median =
      as.numeric(
        median
      ),
    
    q1 =
      as.numeric(
        q1
      ),
    
    q3 =
      as.numeric(
        q3
      ),
    
    mean =
      as.numeric(
        mean
      ),
    
    sd =
      as.numeric(
        sd
      )
  )


# =============================================================================
# 10. STATUS DEFINITIONS
# =============================================================================

status_levels <- c(
  "Hlty",
  "Inf1_P",
  "Seps_P",
  "Shock_P"
)


status_labels <- c(
  
  "Hlty" =
    "Healthy",
  
  "Inf1_P" =
    "Uncomplicated\ninfection",
  
  "Seps_P" =
    "Sepsis",
  
  "Shock_P" =
    "Septic shock"
)


scores_df$status <- factor(
  scores_df$status,
  levels = status_levels
)


group_summary$status <- factor(
  group_summary$status,
  levels = status_levels
)


# =============================================================================
# 11. EXPECTED FROZEN VALUES
# =============================================================================

expected_status_n <- c(
  Hlty = 40,
  Inf1_P = 12,
  Seps_P = 20,
  Shock_P = 19
)


expected_group_median <- c(
  Hlty = -1.85800069834073,
  Inf1_P = 0.998897965331387,
  Seps_P = 1.49268704995069,
  Shock_P = 1.95271153926247
)


expected_primary <- list(
  
  n_case =
    39,
  
  n_control =
    12,
  
  case_median =
    1.76609748815924,
  
  control_median =
    0.998897965331387,
  
  p =
    0.107415536186119,
  
  q =
    0.140495645247327,
  
  auc =
    0.655982905982906,
  
  ci_low =
    0.505583979979148,
  
  ci_high =
    0.806381831986664
)


expected_shock_inf <- list(
  
  n_case =
    19,
  
  n_control =
    12,
  
  case_median =
    1.95271153926247,
  
  control_median =
    0.998897965331387,
  
  p =
    0.00890133604882126,
  
  q =
    0.0222533401220532,
  
  auc =
    0.785087719298246,
  
  ci_low =
    0.617714901041834,
  
  ci_high =
    0.952460537554658
)


expected_ordered_rho <- 0.841859006173835
expected_ordered_p <- 1.44621600658424e-25

expected_kw_chisq <- 68.4531797721729
expected_kw_p <- 9.15015900066155e-15

expected_scaling_rho <- 0.99386845039019


# =============================================================================
# 12. EXTRACT FROZEN PRIMARY / SECONDARY RESULTS
# =============================================================================

primary_row <- comparisons %>%
  
  dplyr::filter(
    comparison ==
      "Sepsis_or_shock_vs_uncomplicated"
  )


shock_inf_row <- comparisons %>%
  
  dplyr::filter(
    comparison ==
      "Shock_vs_uncomplicated"
  )


if (nrow(primary_row) != 1) {
  
  stop(
    "Primary external comparison not uniquely identified."
  )
}


if (nrow(shock_inf_row) != 1) {
  
  stop(
    "Shock vs uncomplicated infection comparison not uniquely identified."
  )
}


ordered_spearman_row <- ordered_tests %>%
  
  dplyr::filter(
    analysis ==
      "Spearman_score_vs_ordered_status"
  )


ordered_kw_row <- ordered_tests %>%
  
  dplyr::filter(
    analysis ==
      "Kruskal_Wallis_across_four_baseline_groups"
  )


if (nrow(ordered_spearman_row) != 1) {
  
  stop(
    "Ordered Spearman result not uniquely identified."
  )
}


if (nrow(ordered_kw_row) != 1) {
  
  stop(
    "Ordered Kruskal-Wallis result not uniquely identified."
  )
}


# =============================================================================
# 13. NUMERICAL AUDIT
# =============================================================================

observed_status_n <- table(
  scores_df$status
)


status_audit <- data.frame(
  
  status =
    status_levels,
  
  observed =
    as.numeric(
      observed_status_n[
        status_levels
      ]
    ),
  
  expected =
    as.numeric(
      expected_status_n[
        status_levels
      ]
    ),
  
  stringsAsFactors = FALSE
)


status_audit$match <-
  status_audit$observed ==
  status_audit$expected


if (!all(
  status_audit$match
)) {
  
  stop(
    "Baseline status-count audit failed."
  )
}


# -----------------------------------------------------------------------------
# Group median audit
# -----------------------------------------------------------------------------

group_median_audit <- group_summary %>%
  
  dplyr::mutate(
    
    expected_median =
      unname(
        expected_group_median[
          as.character(
            status
          )
        ]
      ),
    
    difference =
      median -
      expected_median
  )


if (
  any(
    abs(
      group_median_audit$difference
    ) >
    1e-10
  )
) {
  
  stop(
    "Baseline group-median audit failed."
  )
}


# -----------------------------------------------------------------------------
# Primary result audit
# -----------------------------------------------------------------------------

primary_checks <- c(
  
  primary_row$n_case[1] ==
    expected_primary$n_case,
  
  primary_row$n_control[1] ==
    expected_primary$n_control,
  
  abs(
    primary_row$case_median[1] -
      expected_primary$case_median
  ) <
    1e-10,
  
  abs(
    primary_row$control_median[1] -
      expected_primary$control_median
  ) <
    1e-10,
  
  abs(
    primary_row$p_value[1] -
      expected_primary$p
  ) <
    1e-12,
  
  abs(
    primary_row$p_BH_across_score_comparisons[1] -
      expected_primary$q
  ) <
    1e-12,
  
  abs(
    primary_row$auc_fixed_direction[1] -
      expected_primary$auc
  ) <
    1e-12,
  
  abs(
    primary_row$auc_ci_low[1] -
      expected_primary$ci_low
  ) <
    1e-12,
  
  abs(
    primary_row$auc_ci_high[1] -
      expected_primary$ci_high
  ) <
    1e-12
)


if (!all(
  primary_checks
)) {
  
  stop(
    "Frozen primary external-comparison audit failed."
  )
}


# -----------------------------------------------------------------------------
# Key secondary result audit
# -----------------------------------------------------------------------------

shock_inf_checks <- c(
  
  shock_inf_row$n_case[1] ==
    expected_shock_inf$n_case,
  
  shock_inf_row$n_control[1] ==
    expected_shock_inf$n_control,
  
  abs(
    shock_inf_row$p_value[1] -
      expected_shock_inf$p
  ) <
    1e-12,
  
  abs(
    shock_inf_row$p_BH_across_score_comparisons[1] -
      expected_shock_inf$q
  ) <
    1e-12,
  
  abs(
    shock_inf_row$auc_fixed_direction[1] -
      expected_shock_inf$auc
  ) <
    1e-12,
  
  abs(
    shock_inf_row$auc_ci_low[1] -
      expected_shock_inf$ci_low
  ) <
    1e-12,
  
  abs(
    shock_inf_row$auc_ci_high[1] -
      expected_shock_inf$ci_high
  ) <
    1e-12
)


if (!all(
  shock_inf_checks
)) {
  
  stop(
    "Frozen Shock vs uncomplicated-infection audit failed."
  )
}


# -----------------------------------------------------------------------------
# Ordered result audit
# -----------------------------------------------------------------------------

if (
  abs(
    ordered_spearman_row$statistic[1] -
    expected_ordered_rho
  ) >
  1e-12
) {
  
  stop(
    "Ordered baseline Spearman rho mismatch."
  )
}


if (
  abs(
    ordered_spearman_row$p_value[1] -
    expected_ordered_p
  ) >
  1e-30
) {
  
  stop(
    "Ordered baseline Spearman P mismatch."
  )
}


if (
  abs(
    ordered_kw_row$statistic[1] -
    expected_kw_chisq
  ) >
  1e-10
) {
  
  stop(
    "Ordered baseline Kruskal-Wallis statistic mismatch."
  )
}


if (
  abs(
    ordered_kw_row$p_value[1] -
    expected_kw_p
  ) >
  1e-20
) {
  
  stop(
    "Ordered baseline Kruskal-Wallis P mismatch."
  )
}


# -----------------------------------------------------------------------------
# Scaling sensitivity audit
# -----------------------------------------------------------------------------

if (
  nrow(
    scaling_sensitivity
  ) != 1
) {
  
  stop(
    "Scaling sensitivity result not uniquely identified."
  )
}


if (
  abs(
    as.numeric(
      scaling_sensitivity$spearman_rho[1]
    ) -
    expected_scaling_rho
  ) >
  1e-10
) {
  
  stop(
    "Scaling-sensitivity audit failed."
  )
}


# -----------------------------------------------------------------------------
# Component-gene direction audit
# -----------------------------------------------------------------------------

n_direction_concordant <- sum(
  as.logical(
    gene_direction$direction_concordant
  ),
  na.rm = TRUE
)


n_gene_BH_significant <- sum(
  as.numeric(
    gene_direction$p_BH_five_genes
  ) <
    0.05,
  na.rm = TRUE
)


if (
  n_direction_concordant != 5
) {
  
  stop(
    "Expected all five component genes to retain directional concordance."
  )
}


if (
  n_gene_BH_significant != 0
) {
  
  stop(
    "Unexpected BH-significant individual component gene detected."
  )
}


# -----------------------------------------------------------------------------
# Consolidated audit table
# -----------------------------------------------------------------------------

audit_table <- data.frame(
  
  metric = c(
    "Baseline samples",
    "Healthy n",
    "Uncomplicated infection n",
    "Sepsis n",
    "Septic shock n",
    "Primary case n",
    "Primary control n",
    "Primary case median",
    "Primary control median",
    "Primary Wilcoxon P",
    "Primary BH q",
    "Primary AUC",
    "Primary AUC CI low",
    "Primary AUC CI high",
    "Ordered Spearman rho",
    "Ordered Spearman P",
    "Ordered Kruskal-Wallis chi-square",
    "Ordered Kruskal-Wallis P",
    "Shock vs uncomplicated P",
    "Shock vs uncomplicated BH q",
    "Shock vs uncomplicated AUC",
    "Directionally concordant genes",
    "BH-significant component genes",
    "Scaling sensitivity rho"
  ),
  
  observed = c(
    nrow(
      scores_df
    ),
    observed_status_n["Hlty"],
    observed_status_n["Inf1_P"],
    observed_status_n["Seps_P"],
    observed_status_n["Shock_P"],
    primary_row$n_case[1],
    primary_row$n_control[1],
    primary_row$case_median[1],
    primary_row$control_median[1],
    primary_row$p_value[1],
    primary_row$p_BH_across_score_comparisons[1],
    primary_row$auc_fixed_direction[1],
    primary_row$auc_ci_low[1],
    primary_row$auc_ci_high[1],
    ordered_spearman_row$statistic[1],
    ordered_spearman_row$p_value[1],
    ordered_kw_row$statistic[1],
    ordered_kw_row$p_value[1],
    shock_inf_row$p_value[1],
    shock_inf_row$p_BH_across_score_comparisons[1],
    shock_inf_row$auc_fixed_direction[1],
    n_direction_concordant,
    n_gene_BH_significant,
    scaling_sensitivity$spearman_rho[1]
  ),
  
  expected = c(
    91,
    40,
    12,
    20,
    19,
    expected_primary$n_case,
    expected_primary$n_control,
    expected_primary$case_median,
    expected_primary$control_median,
    expected_primary$p,
    expected_primary$q,
    expected_primary$auc,
    expected_primary$ci_low,
    expected_primary$ci_high,
    expected_ordered_rho,
    expected_ordered_p,
    expected_kw_chisq,
    expected_kw_p,
    expected_shock_inf$p,
    expected_shock_inf$q,
    expected_shock_inf$auc,
    5,
    0,
    expected_scaling_rho
  ),
  
  stringsAsFactors = FALSE
)


audit_table$difference <-
  audit_table$observed -
  audit_table$expected


cat("\n====================================================================\n")
cat("NUMERICAL AUDIT\n")
cat("====================================================================\n")

print(
  audit_table
)


cat("\nNumerical audit passed successfully.\n")


# =============================================================================
# 14. PUBLICATION COLORS
# =============================================================================

col_healthy <- "#56B4E9"
col_infection <- "#7F8C8D"
col_sepsis <- "#E69F00"
col_shock <- "#D55E00"

col_primary_case <- "#D55E00"
col_primary_control <- "#7F8C8D"

col_roc <- "#0072B2"

col_forest <- "#0072B2"
col_reference <- "#7F8C8D"


status_colors <- c(
  
  "Hlty" =
    col_healthy,
  
  "Inf1_P" =
    col_infection,
  
  "Seps_P" =
    col_sepsis,
  
  "Shock_P" =
    col_shock
)


# =============================================================================
# 15. PANEL A — SCORE ACROSS BASELINE STATES
# =============================================================================

panel_A_annotation <- paste0(
  
  "Ordered Spearman \u03c1 = ",
  sprintf(
    "%.3f",
    ordered_spearman_row$statistic[1]
  ),
  
  "\nP = ",
  format_p(
    ordered_spearman_row$p_value[1]
  )
)


p_A <- ggplot2::ggplot(
  
  scores_df,
  
  ggplot2::aes(
    x = status,
    y = five_gene_score,
    fill = status,
    color = status
  )
  
) +
  
  ggplot2::geom_boxplot(
    width = 0.60,
    alpha = 0.40,
    outlier.shape = NA,
    linewidth = 0.55
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.0,
    alpha = 0.78
  ) +
  
  ggplot2::scale_fill_manual(
    values = status_colors
  ) +
  
  ggplot2::scale_color_manual(
    values = status_colors
  ) +
  
  ggplot2::scale_x_discrete(
    labels = status_labels
  ) +
  
  ggplot2::annotate(
    "label",
    x = 1.05,
    y = Inf,
    label = panel_A_annotation,
    hjust = 0,
    vjust = 1.10,
    size = 3.1,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    axis.text.x =
      ggplot2::element_text(
        size = 8.2
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "A",
    
    title =
      "Frozen score across baseline clinical states",
    
    subtitle =
      "GSE154918 baseline samples; n=91",
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 16. PANEL B — PRESPECIFIED PRIMARY COMPARISON
# =============================================================================

primary_plot_data <- scores_df %>%
  
  dplyr::filter(
    status %in%
      c(
        "Inf1_P",
        "Seps_P",
        "Shock_P"
      )
  ) %>%
  
  dplyr::mutate(
    
    primary_group =
      dplyr::if_else(
        status ==
          "Inf1_P",
        "Uncomplicated infection",
        "Sepsis / septic shock"
      ),
    
    primary_group =
      factor(
        primary_group,
        levels = c(
          "Uncomplicated infection",
          "Sepsis / septic shock"
        )
      )
  )


if (
  nrow(
    primary_plot_data
  ) != 51
) {
  
  stop(
    "Primary plotting dataset should contain 51 baseline samples."
  )
}


primary_colors <- c(
  
  "Uncomplicated infection" =
    col_primary_control,
  
  "Sepsis / septic shock" =
    col_primary_case
)


panel_B_annotation <- paste0(
  
  "Median: ",
  sprintf(
    "%.3f",
    primary_row$control_median[1]
  ),
  " vs ",
  sprintf(
    "%.3f",
    primary_row$case_median[1]
  ),
  
  "\nWilcoxon P = ",
  format_p(
    primary_row$p_value[1]
  ),
  
  "\nBH q = ",
  format_q(
    primary_row$p_BH_across_score_comparisons[1]
  )
)


p_B <- ggplot2::ggplot(
  
  primary_plot_data,
  
  ggplot2::aes(
    x = primary_group,
    y = five_gene_score,
    fill = primary_group,
    color = primary_group
  )
  
) +
  
  ggplot2::geom_boxplot(
    width = 0.58,
    alpha = 0.42,
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  
  ggplot2::geom_jitter(
    width = 0.11,
    height = 0,
    size = 2.3,
    alpha = 0.82
  ) +
  
  ggplot2::scale_fill_manual(
    values = primary_colors
  ) +
  
  ggplot2::scale_color_manual(
    values = primary_colors
  ) +
  
  ggplot2::annotate(
    "label",
    x = 1.02,
    y = Inf,
    label = panel_B_annotation,
    hjust = 0,
    vjust = 1.10,
    size = 3.1,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    axis.text.x =
      ggplot2::element_text(
        size = 8.3
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "B",
    
    title =
      "Prespecified primary external comparison",
    
    subtitle =
      paste0(
        "Uncomplicated infection (n=",
        primary_row$n_control[1],
        ") vs sepsis/septic shock (n=",
        primary_row$n_case[1],
        ")"
      ),
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 17. PANEL C — PRIMARY ROC
# =============================================================================
#
# ROC coordinates are recreated from FROZEN sample-level scores solely
# for visualization.
#
# The AUC is then checked against the frozen Script 141 AUC.
#
# =============================================================================

roc_object <- pROC::roc(
  
  response =
    primary_plot_data$primary_group,
  
  predictor =
    primary_plot_data$five_gene_score,
  
  levels = c(
    "Uncomplicated infection",
    "Sepsis / septic shock"
  ),
  
  direction = "<",
  
  quiet = TRUE
)


roc_auc_visual <- as.numeric(
  pROC::auc(
    roc_object
  )
)


if (
  abs(
    roc_auc_visual -
    primary_row$auc_fixed_direction[1]
  ) >
  1e-12
) {
  
  stop(
    paste0(
      "ROC visualization AUC does not reproduce frozen AUC. ",
      "Observed ",
      roc_auc_visual,
      "; frozen ",
      primary_row$auc_fixed_direction[1]
    )
  )
}


roc_coords <- as.data.frame(
  
  pROC::coords(
    
    roc_object,
    
    x = "all",
    
    ret = c(
      "specificity",
      "sensitivity"
    ),
    
    transpose = FALSE
  )
) %>%
  
  dplyr::mutate(
    
    false_positive_rate =
      1 -
      specificity
  )


panel_C_annotation <- paste0(
  
  "AUC = ",
  sprintf(
    "%.3f",
    primary_row$auc_fixed_direction[1]
  ),
  
  "\n95% CI ",
  sprintf(
    "%.3f",
    primary_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    primary_row$auc_ci_high[1]
  )
)


p_C <- ggplot2::ggplot(
  roc_coords,
  ggplot2::aes(
    x = false_positive_rate,
    y = sensitivity
  )
) +
  
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey55"
  ) +
  
  ggplot2::geom_step(
    linewidth = 1.05,
    color = col_roc
  ) +
  
  ggplot2::annotate(
    "label",
    x = 0.05,
    y = 0.95,
    label = panel_C_annotation,
    hjust = 0,
    vjust = 1,
    size = 3.4,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  ggplot2::coord_equal(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    ),
    expand = FALSE
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    
    tag =
      "C",
    
    title =
      "Primary fixed-direction ROC",
    
    subtitle =
      "Sepsis/septic shock versus uncomplicated infection",
    
    x =
      "1 \u2212 specificity",
    
    y =
      "Sensitivity"
  )


# =============================================================================
# 18. PANEL D — SECONDARY SCORE-LEVEL CONTRASTS
# =============================================================================

secondary_comparisons <- comparisons %>%
  
  dplyr::filter(
    comparison !=
      "Sepsis_or_shock_vs_uncomplicated"
  ) %>%
  
  dplyr::mutate(
    
    figure_label =
      dplyr::case_when(
        
        comparison ==
          "Sepsis_or_shock_vs_healthy" ~
          "Sepsis / shock vs healthy",
        
        comparison ==
          "Sepsis_vs_uncomplicated" ~
          "Sepsis vs uncomplicated infection",
        
        comparison ==
          "Shock_vs_uncomplicated" ~
          "Septic shock vs uncomplicated infection",
        
        comparison ==
          "Shock_vs_sepsis" ~
          "Septic shock vs sepsis",
        
        TRUE ~
          comparison
      )
  )


secondary_order <- c(
  "Sepsis / shock vs healthy",
  "Sepsis vs uncomplicated infection",
  "Septic shock vs uncomplicated infection",
  "Septic shock vs sepsis"
)


secondary_comparisons$figure_label <- factor(
  secondary_comparisons$figure_label,
  levels =
    rev(
      secondary_order
    )
)


secondary_comparisons <- secondary_comparisons %>%
  
  dplyr::mutate(
    
    auc_text =
      sprintf(
        "%.3f",
        auc_fixed_direction
      ),
    
    q_text =
      paste0(
        "q=",
        vapply(
          p_BH_across_score_comparisons,
          format_q,
          character(1)
        )
      )
  )


p_D <- ggplot2::ggplot(
  
  secondary_comparisons,
  
  ggplot2::aes(
    x = auc_fixed_direction,
    y = figure_label
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey55"
  ) +
  
  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = auc_ci_low,
      xmax = auc_ci_high
    ),
    width = 0.13,
    linewidth = 0.75,
    color = col_forest
  ) +
  
  ggplot2::geom_point(
    size = 3.8,
    color = col_forest
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = auc_text
    ),
    hjust = -0.35,
    size = 3.0
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = q_text
    ),
    x = 0.31,
    hjust = 0,
    size = 2.8,
    color = "grey30"
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      0.30,
      1.05
    ),
    breaks = seq(
      0.4,
      1.0,
      by = 0.2
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    axis.title.y =
      ggplot2::element_blank(),
    
    axis.text.y =
      ggplot2::element_text(
        size = 8.4
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "D",
    
    title =
      "Secondary external contrasts",
    
    subtitle =
      "Fixed-direction AUC with DeLong 95% CI; BH q shown at left",
    
    x =
      "Area under the ROC curve"
  )


# =============================================================================
# 19. ASSEMBLE FINAL FIGURE 4
# =============================================================================
#
# No global title inside the image.
#
# Scientific Figure 4 title is retained in the figure legend.
#
# =============================================================================

top_row <- (
  p_A |
    p_B
) +
  
  patchwork::plot_layout(
    widths = c(
      1,
      1
    )
  )


bottom_row <- (
  p_C |
    p_D
) +
  
  patchwork::plot_layout(
    widths = c(
      0.85,
      1.15
    )
  )


figure4 <- (
  top_row /
    bottom_row
) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1.02
    )
  )


# =============================================================================
# 20. EXPORT FINAL FIGURE 4
# =============================================================================

figure_png <- file.path(
  figures_dir,
  "149_Figure4_GSE154918_external_validation.png"
)


figure_pdf <- file.path(
  figures_dir,
  "149_Figure4_GSE154918_external_validation.pdf"
)


figure_tiff <- file.path(
  figures_dir,
  "149_Figure4_GSE154918_external_validation.tiff"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure4,
  width = 13.5,
  height = 9.2,
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure4,
  width = 13.5,
  height = 9.2,
  device =
    if (
      capabilities(
        "cairo"
      )
    ) {
      
      grDevices::cairo_pdf
      
    } else {
      
      grDevices::pdf
    },
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure4,
  width = 13.5,
  height = 9.2,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 21. EXPORT INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  
  A_baseline_states =
    p_A,
  
  B_primary_comparison =
    p_B,
  
  C_primary_ROC =
    p_C,
  
  D_secondary_contrasts =
    p_D
)


for (nm in names(
  individual_panels
)) {
  
  ggplot2::ggsave(
    
    filename =
      file.path(
        figures_dir,
        paste0(
          "149_panel_",
          nm,
          ".png"
        )
      ),
    
    plot =
      individual_panels[[nm]],
    
    width = 6.7,
    
    height = 5.2,
    
    dpi = 600,
    
    bg = "white"
  )
}


# =============================================================================
# 22. FIGURE SOURCE WORKBOOK
# =============================================================================

run_info_149 <- data.frame(
  
  item = c(
    "script",
    "figure",
    "analysis_mode",
    "source_script",
    "source_dataset",
    "score_composition",
    "score_definition",
    "primary_external_comparison",
    "feature_selection",
    "coefficient_refitting",
    "cutoff_optimization",
    "score_direction_flipping",
    "followup_samples_primary",
    "ROC_coordinates",
    "interpretation"
  ),
  
  value = c(
    "149_build_Figure4_GSE154918_external_validation.R",
    "Main Figure 4",
    "publication packaging of frozen Script 141 results",
    "141_external_validation_GSE154918.R",
    "GSE154918",
    "CD177; HK3; IRAK3; CARD11; IKZF2",
    "mean z(CD177,HK3,IRAK3) - mean z(CARD11,IKZF2)",
    "Seps_P + Shock_P versus Inf1_P",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "reconstructed from frozen sample-level scores for visualization only; AUC audited against frozen result",
    "external transcriptomic replication/evaluation, not calibrated clinical-assay validation"
  ),
  
  stringsAsFactors = FALSE
)


source_output_workbook <- file.path(
  tables_dir,
  "149_Figure4_source_data.xlsx"
)


openxlsx::write.xlsx(
  
  list(
    
    Run_info =
      run_info_149,
    
    Numerical_audit =
      audit_table,
    
    Status_audit =
      status_audit,
    
    Group_summary =
      group_summary,
    
    Sample_scores =
      scores_df,
    
    Primary_comparison =
      primary_row,
    
    Secondary_comparisons =
      secondary_comparisons,
    
    Ordered_tests =
      ordered_tests,
    
    Pairwise_tests =
      pairwise_tests,
    
    Component_gene_audit =
      gene_direction,
    
    Scaling_sensitivity =
      scaling_sensitivity,
    
    ROC_coordinates =
      roc_coords
    
  ),
  
  source_output_workbook,
  
  overwrite = TRUE
)


# =============================================================================
# 23. PROVENANCE TABLE
# =============================================================================

provenance <- data.frame(
  
  Panel = c(
    "A",
    "B",
    "C",
    "D"
  ),
  
  Content = c(
    "Frozen score across four baseline clinical states",
    "Prespecified primary external comparison",
    "ROC curve for primary external comparison",
    "Secondary score-level external contrasts"
  ),
  
  Upstream_source = c(
    "141 workbook sheets 04 + 08 + 10",
    "141 workbook sheets 04 + 05",
    "141 workbook sheets 04 + 05",
    "141 workbook sheet 05"
  ),
  
  New_hypothesis_testing = c(
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  
  Visualization_only_calculation = c(
    "NO",
    "NO",
    "ROC coordinates only",
    "NO"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  provenance,
  file.path(
    tables_dir,
    "149_Figure4_provenance.csv"
  ),
  row.names = FALSE
)


write.csv(
  audit_table,
  file.path(
    tables_dir,
    "149_Figure4_numerical_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 24. FINAL FIGURE LEGEND — ENGLISH
# =============================================================================

caption_en <- paste0(
  
  "Figure 4. External evaluation in GSE154918 shows directional replication ",
  "of the five-gene score but limited discrimination between sepsis and ",
  "uncomplicated infection. ",
  
  "(A) Distribution of the frozen five-gene host-response score across ",
  "baseline healthy controls (n=40), uncomplicated infection (n=12), sepsis ",
  "(n=20), and septic shock (n=19). The five-gene score was calculated using ",
  "the prespecified composition CD177, HK3, IRAK3, CARD11, and IKZF2, with ",
  "gene-wise standardization across all 91 baseline samples and without ",
  "feature selection, coefficient refitting, cutoff optimization, or post hoc ",
  "direction reversal. Across the ordered baseline states healthy control to ",
  "uncomplicated infection to sepsis to septic shock, the score showed a ",
  "positive rank association (Spearman rho=",
  sprintf(
    "%.3f",
    ordered_spearman_row$statistic[1]
  ),
  ", P=",
  format(
    ordered_spearman_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "). ",
  
  "(B) Prespecified primary external comparison of baseline sepsis/septic ",
  "shock versus uncomplicated infection. Median scores were ",
  sprintf(
    "%.3f",
    primary_row$case_median[1]
  ),
  " and ",
  sprintf(
    "%.3f",
    primary_row$control_median[1]
  ),
  ", respectively (Wilcoxon P=",
  sprintf(
    "%.3f",
    primary_row$p_value[1]
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.3f",
    primary_row$p_BH_across_score_comparisons[1]
  ),
  "). The prespecified primary comparison therefore did not reach statistical ",
  "significance. ",
  
  "(C) Fixed-direction receiver-operating-characteristic curve for the same ",
  "primary comparison. The AUC was ",
  sprintf(
    "%.3f",
    primary_row$auc_fixed_direction[1]
  ),
  " (95% CI ",
  sprintf(
    "%.3f",
    primary_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    primary_row$auc_ci_high[1]
  ),
  "). ",
  
  "(D) Secondary score-level external comparisons summarized by fixed-direction ",
  "AUC and DeLong 95% confidence intervals. Among the within-infection ",
  "comparisons, septic shock showed greater separation from uncomplicated ",
  "infection (AUC=",
  sprintf(
    "%.3f",
    shock_inf_row$auc_fixed_direction[1]
  ),
  ", 95% CI ",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_high[1]
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.4f",
    shock_inf_row$p_BH_across_score_comparisons[1]
  ),
  "). All five component genes retained the prespecified expression direction ",
  "in the primary external contrast, although none remained individually ",
  "significant after correction across the five-gene family. Follow-up samples ",
  "were excluded from the primary analysis. These results support directional ",
  "external transcriptomic replication of the frozen molecular score but do ",
  "not constitute validation of a pre-calibrated clinical diagnostic assay or ",
  "decision threshold."
)


writeLines(
  caption_en,
  file.path(
    text_dir,
    "149_Figure4_caption_EN.txt"
  )
)


# =============================================================================
# 25. FINAL FIGURE LEGEND — RUSSIAN
# =============================================================================

caption_ru <- paste0(
  
  "Рисунок 4. Внешняя оценка в GSE154918 демонстрирует направленную ",
  "репликацию пятигенного score, но ограниченную способность различать сепсис ",
  "и неосложненную инфекцию. ",
  
  "(A) Распределение замороженного пятигенного host-response score среди ",
  "исходных здоровых контролей (n=40), пациентов с неосложненной инфекцией ",
  "(n=12), сепсисом (n=20) и септическим шоком (n=19). Пятигенный score ",
  "рассчитывали на основе заранее зафиксированных CD177, HK3, IRAK3, CARD11 ",
  "и IKZF2 с gene-wise стандартизацией по всем 91 baseline samples без нового ",
  "feature selection, refitting коэффициентов, оптимизации cutoff или post hoc ",
  "смены направления score. В последовательности healthy control -> ",
  "uncomplicated infection -> sepsis -> septic shock наблюдалась положительная ",
  "ранговая ассоциация (Spearman rho=",
  sprintf(
    "%.3f",
    ordered_spearman_row$statistic[1]
  ),
  ", P=",
  format(
    ordered_spearman_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "). ",
  
  "(B) Заранее определенное основное внешнее сравнение sepsis/septic shock ",
  "против uncomplicated infection. Медианы score составляли ",
  sprintf(
    "%.3f",
    primary_row$case_median[1]
  ),
  " и ",
  sprintf(
    "%.3f",
    primary_row$control_median[1]
  ),
  " соответственно (Wilcoxon P=",
  sprintf(
    "%.3f",
    primary_row$p_value[1]
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.3f",
    primary_row$p_BH_across_score_comparisons[1]
  ),
  "). Таким образом, основной prespecified endpoint статистической значимости ",
  "не достиг. ",
  
  "(C) Fixed-direction ROC для того же основного сравнения. AUC составила ",
  sprintf(
    "%.3f",
    primary_row$auc_fixed_direction[1]
  ),
  " (95% ДИ ",
  sprintf(
    "%.3f",
    primary_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    primary_row$auc_ci_high[1]
  ),
  "). ",
  
  "(D) Вторичные score-level сравнения, представленные как fixed-direction ",
  "AUC с 95% ДИ DeLong. Среди сравнений внутри инфекционных состояний наиболее ",
  "выраженное разделение наблюдалось между septic shock и uncomplicated ",
  "infection (AUC=",
  sprintf(
    "%.3f",
    shock_inf_row$auc_fixed_direction[1]
  ),
  ", 95% ДИ ",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_high[1]
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.4f",
    shock_inf_row$p_BH_across_score_comparisons[1]
  ),
  "). Все пять компонентов сохранили заранее ожидаемое направление экспрессии ",
  "в основном внешнем сравнении, однако ни один отдельный ген не сохранял ",
  "значимость после коррекции внутри пятигенной family. Follow-up samples не ",
  "включались в primary analysis. Результаты подтверждают направленную внешнюю ",
  "транскриптомную репликацию замороженного molecular score, но не являются ",
  "валидацией заранее калиброванного клинического диагностического теста или ",
  "decision threshold."
)


writeLines(
  caption_ru,
  file.path(
    text_dir,
    "149_Figure4_caption_RU.txt"
  )
)


writeLines(
  c(
    caption_en,
    "",
    caption_ru
  ),
  file.path(
    text_dir,
    "149_Figure4_caption_EN_RU.txt"
  )
)


# =============================================================================
# 26. RESULTS SECTION 3.8 — FINAL ENGLISH DRAFT
# =============================================================================

results_38_en <- paste0(
  
  "The frozen five-gene host-response signature was next evaluated in the ",
  "independent whole-blood RNA-seq dataset GSE154918 without gene substitution, ",
  "feature selection, coefficient refitting, cutoff optimization, or post hoc ",
  "score-direction reversal. Baseline samples comprised 40 healthy controls, ",
  "12 patients with uncomplicated infection, 20 with sepsis, and 19 with ",
  "septic shock; follow-up samples were excluded from the primary analysis. ",
  "Across these ordered baseline states, the score increased progressively ",
  "(Spearman rho=",
  sprintf(
    "%.3f",
    ordered_spearman_row$statistic[1]
  ),
  ", P=",
  format(
    ordered_spearman_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; Fig. 4A). ",
  
  "However, the prespecified primary external comparison of sepsis/septic ",
  "shock versus uncomplicated infection was formally negative. Median scores ",
  "were ",
  sprintf(
    "%.3f",
    primary_row$case_median[1]
  ),
  " and ",
  sprintf(
    "%.3f",
    primary_row$control_median[1]
  ),
  ", respectively (Wilcoxon P=",
  sprintf(
    "%.3f",
    primary_row$p_value[1]
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.3f",
    primary_row$p_BH_across_score_comparisons[1]
  ),
  "), with a fixed-direction AUC of ",
  sprintf(
    "%.3f",
    primary_row$auc_fixed_direction[1]
  ),
  " (95% CI ",
  sprintf(
    "%.3f",
    primary_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    primary_row$auc_ci_high[1]
  ),
  ") (Fig. 4B,C). ",
  
  "Secondary analyses showed stronger separation between septic shock and ",
  "uncomplicated infection (AUC=",
  sprintf(
    "%.3f",
    shock_inf_row$auc_fixed_direction[1]
  ),
  ", 95% CI ",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_high[1]
  ),
  "; P=",
  format(
    shock_inf_row$p_value[1],
    scientific = FALSE,
    digits = 4
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.4f",
    shock_inf_row$p_BH_across_score_comparisons[1]
  ),
  ") (Fig. 4D). All five component genes retained their prespecified direction ",
  "in the primary comparison, although individual-gene evidence was weak: only ",
  "CARD11 was nominally significant and no component remained significant after ",
  "correction across the five-gene family. Cohort-standardized and ",
  "healthy-reference implementations were highly concordant (Spearman rho=",
  sprintf(
    "%.3f",
    scaling_sensitivity$spearman_rho[1]
  ),
  "). Thus, GSE154918 provided directional external transcriptomic replication ",
  "and evidence of increasing molecular signal with greater clinical severity, ",
  "but did not support significant discrimination of the prespecified combined ",
  "sepsis/septic-shock group from uncomplicated infection."
)


writeLines(
  results_38_en,
  file.path(
    text_dir,
    "149_Results_3.8_EN.txt"
  )
)


# =============================================================================
# 27. RESULTS SECTION 3.8 — RUSSIAN WORKING DRAFT
# =============================================================================

results_38_ru <- paste0(
  
  "Замороженная пятигенная host-response сигнатура была далее оценена во ",
  "внешнем whole-blood RNA-seq наборе GSE154918 без замены генов, нового ",
  "feature selection, refitting коэффициентов, оптимизации cutoff или post hoc ",
  "смены направления score. Baseline cohort включала 40 здоровых контролей, ",
  "12 пациентов с неосложненной инфекцией, 20 с сепсисом и 19 с септическим ",
  "шоком; follow-up samples исключались из primary analysis. В последовательности ",
  "этих baseline состояний score возрастал (Spearman rho=",
  sprintf(
    "%.3f",
    ordered_spearman_row$statistic[1]
  ),
  ", P=",
  format(
    ordered_spearman_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; рис. 4A). ",
  
  "Однако заранее определенное основное внешнее сравнение sepsis/septic shock ",
  "против uncomplicated infection было формально отрицательным. Медианы score ",
  "составляли ",
  sprintf(
    "%.3f",
    primary_row$case_median[1]
  ),
  " и ",
  sprintf(
    "%.3f",
    primary_row$control_median[1]
  ),
  " соответственно (Wilcoxon P=",
  sprintf(
    "%.3f",
    primary_row$p_value[1]
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.3f",
    primary_row$p_BH_across_score_comparisons[1]
  ),
  "), fixed-direction AUC = ",
  sprintf(
    "%.3f",
    primary_row$auc_fixed_direction[1]
  ),
  " (95% ДИ ",
  sprintf(
    "%.3f",
    primary_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    primary_row$auc_ci_high[1]
  ),
  ") (рис. 4B,C). ",
  
  "Во вторичных анализах более выраженное разделение наблюдалось между septic ",
  "shock и uncomplicated infection (AUC=",
  sprintf(
    "%.3f",
    shock_inf_row$auc_fixed_direction[1]
  ),
  ", 95% ДИ ",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    shock_inf_row$auc_ci_high[1]
  ),
  "; P=",
  format(
    shock_inf_row$p_value[1],
    scientific = FALSE,
    digits = 4
  ),
  "; BH-adjusted P=",
  sprintf(
    "%.4f",
    shock_inf_row$p_BH_across_score_comparisons[1]
  ),
  ") (рис. 4D). Все пять генов сохранили заранее ожидаемое направление эффекта ",
  "в основном сравнении, однако доказательства на уровне отдельных genes были ",
  "слабыми: только CARD11 был nominally significant, и ни один компонент не ",
  "сохранял значимость после коррекции внутри five-gene family. ",
  "Cohort-standardized и healthy-reference implementations были практически ",
  "полностью согласованы (Spearman rho=",
  sprintf(
    "%.3f",
    scaling_sensitivity$spearman_rho[1]
  ),
  "). Таким образом, GSE154918 подтверждает направленную внешнюю ",
  "транскриптомную репликацию и усиление molecular signal с увеличением ",
  "клинической тяжести, но не подтверждает статистически значимое различение ",
  "заранее определенной объединенной группы sepsis/septic shock от ",
  "uncomplicated infection."
)


writeLines(
  results_38_ru,
  file.path(
    text_dir,
    "149_Results_3.8_RU.txt"
  )
)


# =============================================================================
# 28. FIGURE-TO-RESULTS PLACEMENT
# =============================================================================

placement <- data.frame(
  
  Results_section = c(
    "3.8",
    "3.8",
    "3.8",
    "3.8"
  ),
  
  Panel = c(
    "4A",
    "4B",
    "4C",
    "4D"
  ),
  
  Content = c(
    "Frozen score across baseline clinical states",
    "Prespecified primary external comparison",
    "Primary fixed-direction ROC",
    "Secondary external contrasts"
  ),
  
  Recommended_first_citation = c(
    "Sentence reporting ordered baseline-state association",
    "Sentence reporting formally negative primary comparison",
    "Sentence reporting primary AUC and 95% CI",
    "Sentence reporting secondary shock-versus-uncomplicated-infection result"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  placement,
  file.path(
    tables_dir,
    "149_Figure4_Results_placement.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 29. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "149_sessionInfo.txt"
  )
)


# =============================================================================
# 30. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 149 completed successfully.\n")
cat("====================================================================\n\n")


cat("BASELINE GSE154918 COHORT\n")
cat("------------------------\n")


cat(
  "Total baseline n = ",
  nrow(
    scores_df
  ),
  "\n",
  sep = ""
)


print(
  table(
    scores_df$status
  )
)


cat("\nORDERED BASELINE ANALYSIS\n")
cat("-------------------------\n")


cat(
  "Spearman rho = ",
  sprintf(
    "%.6f",
    ordered_spearman_row$statistic[1]
  ),
  "\n",
  sep = ""
)


cat(
  "P = ",
  format(
    ordered_spearman_row$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat(
  "Kruskal-Wallis chi-square = ",
  sprintf(
    "%.6f",
    ordered_kw_row$statistic[1]
  ),
  "\n",
  sep = ""
)


cat(
  "P = ",
  format(
    ordered_kw_row$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat("\nPRIMARY EXTERNAL COMPARISON\n")
cat("---------------------------\n")


cat(
  "Sepsis/septic shock vs uncomplicated infection\n"
)


cat(
  "n = ",
  primary_row$n_case[1],
  " vs ",
  primary_row$n_control[1],
  "\n",
  sep = ""
)


cat(
  "Median = ",
  sprintf(
    "%.6f",
    primary_row$case_median[1]
  ),
  " vs ",
  sprintf(
    "%.6f",
    primary_row$control_median[1]
  ),
  "\n",
  sep = ""
)


cat(
  "Wilcoxon P = ",
  format(
    primary_row$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat(
  "BH q = ",
  sprintf(
    "%.8f",
    primary_row$p_BH_across_score_comparisons[1]
  ),
  "\n",
  sep = ""
)


cat(
  "AUC = ",
  sprintf(
    "%.6f",
    primary_row$auc_fixed_direction[1]
  ),
  "\n",
  sep = ""
)


cat(
  "95% CI = [",
  sprintf(
    "%.6f",
    primary_row$auc_ci_low[1]
  ),
  ", ",
  sprintf(
    "%.6f",
    primary_row$auc_ci_high[1]
  ),
  "]\n",
  sep = ""
)


cat("\nKEY SECONDARY COMPARISON\n")
cat("------------------------\n")


cat(
  "Septic shock vs uncomplicated infection\n"
)


cat(
  "Wilcoxon P = ",
  format(
    shock_inf_row$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat(
  "BH q = ",
  sprintf(
    "%.8f",
    shock_inf_row$p_BH_across_score_comparisons[1]
  ),
  "\n",
  sep = ""
)


cat(
  "AUC = ",
  sprintf(
    "%.6f",
    shock_inf_row$auc_fixed_direction[1]
  ),
  "\n",
  sep = ""
)


cat(
  "95% CI = [",
  sprintf(
    "%.6f",
    shock_inf_row$auc_ci_low[1]
  ),
  ", ",
  sprintf(
    "%.6f",
    shock_inf_row$auc_ci_high[1]
  ),
  "]\n",
  sep = ""
)


cat("\nCOMPONENT-GENE AUDIT\n")
cat("--------------------\n")


cat(
  "Directionally concordant genes = ",
  n_direction_concordant,
  "/5\n",
  sep = ""
)


cat(
  "BH-significant individual genes = ",
  n_gene_BH_significant,
  "/5\n",
  sep = ""
)


cat("\nSCALING SENSITIVITY\n")
cat("-------------------\n")


cat(
  "Cohort-standardized vs healthy-reference rho = ",
  sprintf(
    "%.6f",
    scaling_sensitivity$spearman_rho[1]
  ),
  "\n",
  sep = ""
)


cat("\nMAIN FIGURE\n")
cat("-----------\n")


cat(
  normalizePath(
    figure_png,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nSOURCE WORKBOOK\n")
cat("---------------\n")


cat(
  normalizePath(
    source_output_workbook,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nFIGURE LEGEND\n")
cat("-------------\n")


cat(
  normalizePath(
    file.path(
      text_dir,
      "149_Figure4_caption_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nRESULTS 3.8\n")
cat("-----------\n")


cat(
  normalizePath(
    file.path(
      text_dir,
      "149_Results_3.8_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nFIGURE PLACEMENT\n")
cat("----------------\n")


cat(
  "Figure 4A-D -> Results Section 3.8\n"
)


cat(
  "Physical placement -> immediately after Section 3.8\n"
)


cat("\nINTERPRETATION GUARDRAILS\n")
cat("-------------------------\n")


cat(
  "- Primary external endpoint is NEGATIVE.\n"
)


cat(
  "- Directional external transcriptomic replication is supported.\n"
)


cat(
  "- Stronger shock-vs-uncomplicated result is SECONDARY.\n"
)


cat(
  "- Do not describe GSE154918 as validating a calibrated diagnostic assay.\n"
)


cat(
  "- Do not claim all five genes independently validated.\n"
)


cat(
  "- Follow-up samples are not part of primary analysis.\n"
)


cat("\nDone.\n")