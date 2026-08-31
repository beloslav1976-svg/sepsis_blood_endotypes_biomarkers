################################################################################
# Script 163
# FINAL v2
#
# Supplementary Table S9
#
# External evaluation of the frozen five-gene host-response score
# in GSE154918
#
# Project:
#   Sepsis_DESeq2
#
#
# PURPOSE
# -------
#
# Package the already frozen Script 141 external-evaluation results
# for GSE154918 into publication-ready Supplementary Table S9.
#
#
# IMPORTANT
# ---------
#
# This script DOES NOT:
#
#   - recompute the five-gene score
#   - rerun any Wilcoxon test
#   - rerun any ROC analysis
#   - recalculate confidence intervals
#   - rerun Spearman correlation
#   - rerun Kruskal-Wallis testing
#   - recalculate BH-adjusted P values
#   - perform feature selection
#   - refit coefficients
#   - optimize a cutoff
#   - reverse score direction
#
# All statistical values are copied from frozen Script 141 outputs.
#
#
# PRIMARY EXTERNAL COMPARISON
# ---------------------------
#
# Sepsis + septic shock versus uncomplicated infection
#
# Frozen anchors:
#
# n_case                = 39
# n_control             = 12
# case median           = 1.766097
# control median        = 0.998898
# Wilcoxon P            = 0.1074155
# BH q                  = 0.1404956
# fixed-direction AUC   = 0.6559829
# 95% CI                = 0.5055840 to 0.8063818
#
#
# SECONDARY SHOCK CONTRAST
# ------------------------
#
# Septic shock versus uncomplicated infection
#
# P                     = 0.008901336
# BH q                  = 0.02225334
# fixed-direction AUC   = 0.7850877
# 95% CI                = 0.6177149 to 0.9524605
#
#
# ORDERED BASELINE STRUCTURE
# --------------------------
#
# Spearman rho          = 0.8418590
# P                     = 1.446216e-25
#
# Kruskal-Wallis chi^2  = 68.45318
# P                     = 9.150159e-15
#
#
# COMPONENT GENES
# ---------------
#
# 5/5 directionally concordant
# 1/5 nominal P < 0.05 (CARD11)
# 0/5 BH-adjusted P < 0.05
#
################################################################################


cat("====================================================================\n")
cat("Running Script 163 FINAL v2\n")
cat("Supplementary Table S9\n")
cat("GSE154918 external evaluation\n")
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
  "readxl",
  "openxlsx"
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
  library(readxl)
  library(openxlsx)
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "163_TableS9_GSE154918_external_evaluation"
)


tables_dir <- file.path(
  output_dir,
  "tables"
)


audit_dir <- file.path(
  output_dir,
  "audit"
)


text_dir <- file.path(
  output_dir,
  "text"
)


for (
  one_dir in c(
    output_dir,
    tables_dir,
    audit_dir,
    text_dir
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
# 4. CANONICAL SCRIPT 141 SOURCE DIRECTORY
# =============================================================================

source_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "141_external_validation_GSE154918"
)


source_tables_dir <- file.path(
  source_dir,
  "tables"
)


if (!dir.exists(source_tables_dir)) {
  stop("Frozen Script 141 tables directory not found.")
}


# =============================================================================
# 5. EXACT FROZEN SOURCE FILES
# =============================================================================

status_file <- file.path(
  source_tables_dir,
  "141_GSE154918_status_count_check.csv"
)


group_summary_file <- file.path(
  source_tables_dir,
  "141_score_group_summary.csv"
)


comparisons_file <- file.path(
  source_tables_dir,
  "141_external_score_comparisons.csv"
)


ordered_file <- file.path(
  source_tables_dir,
  "141_ordered_baseline_status_tests.csv"
)


pairwise_file <- file.path(
  source_tables_dir,
  "141_pairwise_baseline_group_tests.csv"
)


gene_audit_file <- file.path(
  source_tables_dir,
  "141_external_five_gene_direction_audit.csv"
)


scaling_file <- file.path(
  source_tables_dir,
  "141_score_scaling_sensitivity.csv"
)


source_workbook <- file.path(
  source_tables_dir,
  "141_GSE154918_external_validation.xlsx"
)


source_files <- c(
  status_file,
  group_summary_file,
  comparisons_file,
  ordered_file,
  pairwise_file,
  gene_audit_file,
  scaling_file,
  source_workbook
)


missing_source_files <- source_files[
  !file.exists(
    source_files
  )
]


if (length(missing_source_files) > 0) {
  
  stop(
    paste0(
      "Missing frozen Script 141 source file(s):\n",
      paste(
        missing_source_files,
        collapse = "\n"
      )
    )
  )
}


cat("\nCANONICAL FROZEN SCRIPT 141 SOURCES\n")
cat("-----------------------------------\n")

print(
  normalizePath(
    source_files,
    winslash = "\\",
    mustWork = TRUE
  )
)


# =============================================================================
# 6. OPTIONAL FROZEN MAIN FIGURE 4 NUMERICAL AUDIT
# =============================================================================

figure4_audit_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "149_Figure4_GSE154918_external_validation",
  "tables",
  "149_Figure4_numerical_audit.csv"
)


figure4_audit_available <- file.exists(
  figure4_audit_file
)


cat("\nMain Figure 4 frozen numerical audit available: ")

cat(
  figure4_audit_available,
  "\n"
)


# =============================================================================
# 7. READ CANONICAL CSV SOURCES
# =============================================================================

status_df <- read.csv(
  status_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


group_summary_df <- read.csv(
  group_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


comparisons_df <- read.csv(
  comparisons_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


ordered_df <- read.csv(
  ordered_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


pairwise_df <- read.csv(
  pairwise_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


gene_audit_df <- read.csv(
  gene_audit_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


scaling_df <- read.csv(
  scaling_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# =============================================================================
# 8. READ MATCHED WORKBOOK COPIES
# =============================================================================

status_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "01_status_count_check"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


group_summary_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "10_score_group_summary"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


comparisons_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "05_score_comparisons"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


ordered_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "08_ordered_tests"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


pairwise_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "09_pairwise_tests"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


gene_audit_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "07_gene_direction_audit"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


scaling_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "12_scaling_sensitivity"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


run_info_df <- readxl::read_excel(
  source_workbook,
  sheet = "00_run_info"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


# =============================================================================
# 9. EXPECTED SCHEMA AUDIT
# =============================================================================

expected_status_columns <- c(
  "status",
  "expected_n",
  "observed_n",
  "count_matches_expected"
)


expected_group_columns <- c(
  "status",
  "n",
  "median",
  "q1",
  "q3",
  "mean",
  "sd"
)


expected_comparison_columns <- c(
  "comparison",
  "score",
  "case_status",
  "control_status",
  "n_case",
  "n_control",
  "case_median",
  "control_median",
  "median_difference_case_minus_control",
  "wilcoxon_W",
  "p_value",
  "auc_fixed_direction",
  "auc_ci_low",
  "auc_ci_high",
  "p_BH_across_score_comparisons"
)


expected_ordered_columns <- c(
  "analysis",
  "statistic",
  "statistic_name",
  "p_value"
)


expected_pairwise_columns <- c(
  "group_1",
  "group_2",
  "BH_adjusted_p"
)


expected_gene_columns <- c(
  "gene",
  "expected_direction",
  "observed_direction",
  "direction_concordant",
  "n_case",
  "n_control",
  "case_median",
  "control_median",
  "median_difference_case_minus_control",
  "wilcoxon_W",
  "p_value",
  "p_BH_five_genes"
)


expected_scaling_columns <- c(
  "comparison",
  "n",
  "spearman_rho",
  "p_value"
)


schema_checks <- data.frame(
  
  object = c(
    "status counts",
    "score group summary",
    "score comparisons",
    "ordered tests",
    "pairwise tests",
    "gene direction audit",
    "scaling sensitivity"
  ),
  
  schema_match = c(
    identical(
      names(status_df),
      expected_status_columns
    ),
    identical(
      names(group_summary_df),
      expected_group_columns
    ),
    identical(
      names(comparisons_df),
      expected_comparison_columns
    ),
    identical(
      names(ordered_df),
      expected_ordered_columns
    ),
    identical(
      names(pairwise_df),
      expected_pairwise_columns
    ),
    identical(
      names(gene_audit_df),
      expected_gene_columns
    ),
    identical(
      names(scaling_df),
      expected_scaling_columns
    )
  ),
  
  stringsAsFactors = FALSE
)


cat("\nEXACT FROZEN SCHEMA AUDIT\n")
cat("-------------------------\n")

print(
  schema_checks,
  row.names = FALSE
)


if (
  !all(
    schema_checks$schema_match
  )
) {
  
  stop(
    "At least one frozen Script 141 table has an unexpected schema."
  )
}


# =============================================================================
# 10. CSV-vs-XLSX EQUIVALENCE FUNCTIONS
# =============================================================================

normalize_text <- function(x) {
  
  out <- trimws(
    as.character(x)
  )
  
  out[
    out == ""
  ] <- NA_character_
  
  out
}


column_equivalent <- function(a, b) {
  
  if (
    length(a) !=
    length(b)
  ) {
    return(FALSE)
  }
  
  
  a_num <- suppressWarnings(
    as.numeric(
      as.character(a)
    )
  )
  
  
  b_num <- suppressWarnings(
    as.numeric(
      as.character(b)
    )
  )
  
  
  numeric_candidate <-
    all(
      is.na(a) |
        is.finite(a_num)
    ) &&
    all(
      is.na(b) |
        is.finite(b_num)
    )
  
  
  if (numeric_candidate) {
    
    both_missing <-
      is.na(a_num) &
      is.na(b_num)
    
    
    both_present <-
      is.finite(a_num) &
      is.finite(b_num)
    
    
    if (
      !all(
        both_missing |
        both_present
      )
    ) {
      return(FALSE)
    }
    
    
    differences <- abs(
      a_num[both_present] -
        b_num[both_present]
    )
    
    
    if (
      length(differences) ==
      0
    ) {
      return(TRUE)
    }
    
    
    return(
      all(
        differences <
          1e-10
      )
    )
  }
  
  
  a_text <- normalize_text(a)
  b_text <- normalize_text(b)
  
  
  same_missing <-
    is.na(a_text) ==
    is.na(b_text)
  
  
  same_value <-
    is.na(a_text) |
    a_text ==
    b_text
  
  
  all(
    same_missing &
      same_value
  )
}


table_equivalence_audit <- function(
    csv_table,
    xlsx_table,
    source_name
) {
  
  dimensions_match <-
    identical(
      dim(csv_table),
      dim(xlsx_table)
    )
  
  
  names_match <-
    identical(
      names(csv_table),
      names(xlsx_table)
    )
  
  
  if (
    !dimensions_match ||
    !names_match
  ) {
    
    return(
      data.frame(
        source = source_name,
        dimensions_match = dimensions_match,
        columns_match = names_match,
        content_match = FALSE,
        overall_match = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  column_matches <- vapply(
    names(csv_table),
    function(column_name) {
      
      column_equivalent(
        csv_table[[column_name]],
        xlsx_table[[column_name]]
      )
    },
    logical(1)
  )
  
  
  content_match <- all(
    column_matches
  )
  
  
  data.frame(
    source = source_name,
    dimensions_match = dimensions_match,
    columns_match = names_match,
    content_match = content_match,
    overall_match =
      dimensions_match &&
      names_match &&
      content_match,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 11. DUPLICATE-SOURCE EQUIVALENCE AUDIT
# =============================================================================

equivalence_audit <- dplyr::bind_rows(
  
  table_equivalence_audit(
    status_df,
    status_xlsx,
    "Status counts"
  ),
  
  table_equivalence_audit(
    group_summary_df,
    group_summary_xlsx,
    "Score group summary"
  ),
  
  table_equivalence_audit(
    comparisons_df,
    comparisons_xlsx,
    "Score comparisons"
  ),
  
  table_equivalence_audit(
    ordered_df,
    ordered_xlsx,
    "Ordered baseline tests"
  ),
  
  table_equivalence_audit(
    pairwise_df,
    pairwise_xlsx,
    "Pairwise baseline tests"
  ),
  
  table_equivalence_audit(
    gene_audit_df,
    gene_audit_xlsx,
    "Five-gene direction audit"
  ),
  
  table_equivalence_audit(
    scaling_df,
    scaling_xlsx,
    "Scaling sensitivity"
  )
)


cat("\nCSV-vs-WORKBOOK EQUIVALENCE AUDIT\n")
cat("---------------------------------\n")

print(
  equivalence_audit,
  row.names = FALSE
)


if (
  !all(
    equivalence_audit$overall_match
  )
) {
  
  stop(
    paste0(
      "At least one canonical Script 141 CSV differs from ",
      "its corresponding workbook copy."
    )
  )
}


cat(
  "All seven canonical Script 141 CSV tables are equivalent ",
  "to their workbook copies.\n"
)


# =============================================================================
# 12. RUN-INFO GUARDRAIL AUDIT
# =============================================================================

run_info_lookup <- setNames(
  as.character(
    run_info_df$value
  ),
  as.character(
    run_info_df$parameter
  )
)


required_run_info <- c(
  "GEO_accession",
  "primary_panel",
  "score_definition",
  "primary_external_comparison",
  "feature_selection_external",
  "coefficient_refitting_external",
  "cutoff_optimization_external",
  "direction_flipping_external",
  "followup_in_primary_analysis"
)


missing_run_info <- setdiff(
  required_run_info,
  names(
    run_info_lookup
  )
)


if (
  length(
    missing_run_info
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing required Script 141 run-info parameters: ",
      paste(
        missing_run_info,
        collapse = ", "
      )
    )
  )
}


guardrail_audit <- data.frame(
  
  item = required_run_info,
  
  value = unname(
    run_info_lookup[
      required_run_info
    ]
  ),
  
  stringsAsFactors = FALSE
)


cat("\nSCRIPT 141 RUN-INFO GUARDRAILS\n")
cat("------------------------------\n")

print(
  guardrail_audit,
  row.names = FALSE
)


if (
  run_info_lookup[["GEO_accession"]] !=
  "GSE154918"
) {
  
  stop(
    "Unexpected GEO accession in frozen Script 141 run information."
  )
}


expected_no_fields <- c(
  "feature_selection_external",
  "coefficient_refitting_external",
  "cutoff_optimization_external",
  "direction_flipping_external",
  "followup_in_primary_analysis"
)


if (
  !all(
    toupper(
      unname(
        run_info_lookup[
          expected_no_fields
        ]
      )
    ) ==
    "NO"
  )
) {
  
  stop(
    "At least one frozen external-evaluation guardrail is not NO."
  )
}


# =============================================================================
# 13. SAMPLE-COUNT AUDIT
# =============================================================================

expected_status_counts <- data.frame(
  
  status = c(
    "Hlty",
    "Inf1_P",
    "Seps_P",
    "Shock_P",
    "Seps_FU",
    "Shock_FU"
  ),
  
  expected_n_anchor = c(
    40,
    12,
    20,
    19,
    4,
    10
  ),
  
  stringsAsFactors = FALSE
)


status_anchor_audit <- expected_status_counts %>%
  
  dplyr::left_join(
    status_df,
    by = "status"
  ) %>%
  
  dplyr::mutate(
    
    anchor_match =
      observed_n ==
      expected_n_anchor
  )


cat("\nSTATUS COUNT AUDIT\n")
cat("------------------\n")

print(
  status_anchor_audit,
  row.names = FALSE
)


if (
  nrow(status_anchor_audit) !=
  6 ||
  !all(
    status_anchor_audit$anchor_match
  )
) {
  
  stop(
    "Frozen GSE154918 sample-count audit failed."
  )
}


baseline_status <- c(
  "Hlty",
  "Inf1_P",
  "Seps_P",
  "Shock_P"
)


baseline_n <- sum(
  status_df$observed_n[
    status_df$status %in%
      baseline_status
  ]
)


followup_n <- sum(
  status_df$observed_n[
    status_df$status %in%
      c(
        "Seps_FU",
        "Shock_FU"
      )
  ]
)


if (
  baseline_n !=
  91
) {
  
  stop(
    paste0(
      "Expected 91 baseline samples; observed ",
      baseline_n,
      "."
    )
  )
}


if (
  followup_n !=
  14
) {
  
  stop(
    paste0(
      "Expected 14 follow-up samples; observed ",
      followup_n,
      "."
    )
  )
}


# =============================================================================
# 14. GROUP-SUMMARY AUDIT
# =============================================================================

expected_baseline_groups <- c(
  "Hlty",
  "Inf1_P",
  "Seps_P",
  "Shock_P"
)


if (
  nrow(
    group_summary_df
  ) !=
  4 ||
  !setequal(
    group_summary_df$status,
    expected_baseline_groups
  )
) {
  
  stop(
    "Frozen score-group summary does not contain the expected four baseline groups."
  )
}


group_summary_publication <- group_summary_df %>%
  
  dplyr::mutate(
    
    Group =
      dplyr::recode(
        status,
        "Hlty" = "Healthy",
        "Inf1_P" = "Uncomplicated infection",
        "Seps_P" = "Sepsis",
        "Shock_P" = "Septic shock"
      )
  ) %>%
  
  dplyr::select(
    Group,
    n,
    median,
    q1,
    q3,
    mean,
    sd
  )


# =============================================================================
# 15. PRIMARY COMPARISON
# =============================================================================

primary_comparison <- comparisons_df %>%
  
  dplyr::filter(
    comparison ==
      "Sepsis_or_shock_vs_uncomplicated"
  )


if (
  nrow(
    primary_comparison
  ) !=
  1
) {
  
  stop(
    "Expected exactly one frozen primary GSE154918 comparison."
  )
}


# =============================================================================
# 16. PRIMARY FROZEN ANCHOR AUDIT
# =============================================================================

primary_anchor_audit <- data.frame(
  
  metric = c(
    "n_case",
    "n_control",
    "case_median",
    "control_median",
    "P",
    "BH_q",
    "AUC",
    "AUC_CI_low",
    "AUC_CI_high"
  ),
  
  observed = c(
    primary_comparison$n_case,
    primary_comparison$n_control,
    primary_comparison$case_median,
    primary_comparison$control_median,
    primary_comparison$p_value,
    primary_comparison$p_BH_across_score_comparisons,
    primary_comparison$auc_fixed_direction,
    primary_comparison$auc_ci_low,
    primary_comparison$auc_ci_high
  ),
  
  expected = c(
    39,
    12,
    1.766097,
    0.998898,
    0.1074155,
    0.1404956,
    0.6559829,
    0.5055840,
    0.8063818
  ),
  
  tolerance = c(
    0,
    0,
    1e-5,
    1e-5,
    1e-6,
    1e-6,
    1e-6,
    1e-6,
    1e-6
  ),
  
  stringsAsFactors = FALSE
) %>%
  
  dplyr::mutate(
    
    absolute_difference =
      abs(
        observed -
          expected
      ),
    
    pass =
      absolute_difference <=
      tolerance
  )


cat("\nPRIMARY EXTERNAL COMPARISON AUDIT\n")
cat("---------------------------------\n")

print(
  primary_anchor_audit,
  row.names = FALSE
)


if (
  !all(
    primary_anchor_audit$pass
  )
) {
  
  stop(
    "Frozen primary GSE154918 external-comparison anchors failed."
  )
}


# =============================================================================
# 17. ALL FIVE FROZEN SCORE COMPARISONS
# =============================================================================

if (
  nrow(
    comparisons_df
  ) !=
  5
) {
  
  stop(
    "Expected exactly five frozen score comparisons."
  )
}


comparison_labels <- c(
  
  "Sepsis_or_shock_vs_uncomplicated" =
    "Sepsis/septic shock vs uncomplicated infection",
  
  "Sepsis_or_shock_vs_healthy" =
    "Sepsis/septic shock vs healthy",
  
  "Sepsis_vs_uncomplicated" =
    "Sepsis vs uncomplicated infection",
  
  "Shock_vs_uncomplicated" =
    "Septic shock vs uncomplicated infection",
  
  "Shock_vs_sepsis" =
    "Septic shock vs sepsis"
)


comparisons_publication <- comparisons_df %>%
  
  dplyr::mutate(
    
    Comparison =
      unname(
        comparison_labels[
          comparison
        ]
      ),
    
    Analysis_role =
      dplyr::case_when(
        
        comparison ==
          "Sepsis_or_shock_vs_uncomplicated" ~
          "Prespecified primary external comparison",
        
        comparison ==
          "Shock_vs_uncomplicated" ~
          "Secondary severity-oriented comparison",
        
        comparison ==
          "Sepsis_or_shock_vs_healthy" ~
          "Contextual healthy-reference comparison",
        
        TRUE ~
          "Secondary exploratory comparison"
      )
  ) %>%
  
  dplyr::select(
    Analysis_role,
    Comparison,
    case_status,
    control_status,
    n_case,
    n_control,
    case_median,
    control_median,
    median_difference_case_minus_control,
    wilcoxon_W,
    p_value,
    p_BH_across_score_comparisons,
    auc_fixed_direction,
    auc_ci_low,
    auc_ci_high
  )


# =============================================================================
# 18. SHOCK-vs-UNCOMPLICATED FROZEN ANCHOR
# =============================================================================

shock_comparison <- comparisons_df %>%
  
  dplyr::filter(
    comparison ==
      "Shock_vs_uncomplicated"
  )


if (
  nrow(
    shock_comparison
  ) !=
  1
) {
  
  stop(
    "Expected exactly one shock-versus-uncomplicated comparison."
  )
}


shock_anchor_audit <- data.frame(
  
  metric = c(
    "n_case",
    "n_control",
    "P",
    "BH_q",
    "AUC",
    "AUC_CI_low",
    "AUC_CI_high"
  ),
  
  observed = c(
    shock_comparison$n_case,
    shock_comparison$n_control,
    shock_comparison$p_value,
    shock_comparison$p_BH_across_score_comparisons,
    shock_comparison$auc_fixed_direction,
    shock_comparison$auc_ci_low,
    shock_comparison$auc_ci_high
  ),
  
  expected = c(
    19,
    12,
    0.008901336,
    0.02225334,
    0.7850877,
    0.6177149,
    0.9524605
  ),
  
  tolerance = c(
    0,
    0,
    1e-8,
    1e-8,
    1e-6,
    1e-6,
    1e-6
  ),
  
  stringsAsFactors = FALSE
) %>%
  
  dplyr::mutate(
    
    absolute_difference =
      abs(
        observed -
          expected
      ),
    
    pass =
      absolute_difference <=
      tolerance
  )


cat("\nSHOCK-vs-UNCOMPLICATED AUDIT\n")
cat("----------------------------\n")

print(
  shock_anchor_audit,
  row.names = FALSE
)


if (
  !all(
    shock_anchor_audit$pass
  )
) {
  
  stop(
    "Frozen shock-versus-uncomplicated anchors failed."
  )
}


# =============================================================================
# 19. ORDERED BASELINE TESTS
# =============================================================================

if (
  nrow(
    ordered_df
  ) !=
  2
) {
  
  stop(
    "Expected exactly two frozen ordered-baseline tests."
  )
}


spearman_row <- ordered_df %>%
  
  dplyr::filter(
    analysis ==
      "Spearman_score_vs_ordered_status"
  )


kw_row <- ordered_df %>%
  
  dplyr::filter(
    analysis ==
      "Kruskal_Wallis_across_four_baseline_groups"
  )


if (
  nrow(
    spearman_row
  ) !=
  1 ||
  nrow(
    kw_row
  ) !=
  1
) {
  
  stop(
    "Frozen ordered-baseline tests could not be uniquely identified."
  )
}


ordered_anchor_audit <- data.frame(
  
  metric = c(
    "Ordered_Spearman_rho",
    "Ordered_Spearman_P",
    "Kruskal_Wallis_chi_square",
    "Kruskal_Wallis_P"
  ),
  
  observed = c(
    spearman_row$statistic,
    spearman_row$p_value,
    kw_row$statistic,
    kw_row$p_value
  ),
  
  expected = c(
    0.8418590,
    1.446216e-25,
    68.45318,
    9.150159e-15
  ),
  
  tolerance = c(
    1e-6,
    1e-30,
    1e-4,
    1e-20
  ),
  
  stringsAsFactors = FALSE
) %>%
  
  dplyr::mutate(
    
    absolute_difference =
      abs(
        observed -
          expected
      ),
    
    pass =
      absolute_difference <=
      tolerance
  )


cat("\nORDERED BASELINE-STATE AUDIT\n")
cat("----------------------------\n")

print(
  ordered_anchor_audit,
  row.names = FALSE
)


if (
  !all(
    ordered_anchor_audit$pass
  )
) {
  
  stop(
    "Frozen ordered-baseline anchors failed."
  )
}


ordered_publication <- ordered_df %>%
  
  dplyr::mutate(
    
    Analysis =
      dplyr::recode(
        analysis,
        
        "Spearman_score_vs_ordered_status" =
          "Five-gene score vs ordered baseline disease state",
        
        "Kruskal_Wallis_across_four_baseline_groups" =
          "Five-gene score across four baseline groups"
      )
  ) %>%
  
  dplyr::select(
    Analysis,
    statistic_name,
    statistic,
    p_value
  )


# =============================================================================
# 20. PAIRWISE BASELINE TESTS
# =============================================================================

if (
  nrow(
    pairwise_df
  ) !=
  6
) {
  
  stop(
    "Expected exactly six frozen pairwise baseline-group tests."
  )
}


status_display <- c(
  "Hlty" = "Healthy",
  "Inf1_P" = "Uncomplicated infection",
  "Seps_P" = "Sepsis",
  "Shock_P" = "Septic shock"
)


pairwise_publication <- pairwise_df %>%
  
  dplyr::mutate(
    
    Group_1 =
      unname(
        status_display[
          group_1
        ]
      ),
    
    Group_2 =
      unname(
        status_display[
          group_2
        ]
      )
  ) %>%
  
  dplyr::select(
    Group_1,
    Group_2,
    BH_adjusted_p
  )


# =============================================================================
# 21. COMPONENT-GENE DIRECTIONAL AUDIT
# =============================================================================

expected_genes <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "CARD11",
  "IKZF2"
)


if (
  nrow(
    gene_audit_df
  ) !=
  5 ||
  !setequal(
    gene_audit_df$gene,
    expected_genes
  )
) {
  
  stop(
    "Frozen five-gene direction audit does not contain the expected five genes."
  )
}


direction_concordant_n <- sum(
  gene_audit_df$direction_concordant
)


nominal_gene_n <- sum(
  gene_audit_df$p_value <
    0.05,
  na.rm = TRUE
)


bh_gene_n <- sum(
  gene_audit_df$p_BH_five_genes <
    0.05,
  na.rm = TRUE
)


cat("\nCOMPONENT-GENE DIRECTIONAL REPLICATION\n")
cat("--------------------------------------\n")

cat(
  "Directionally concordant genes = ",
  direction_concordant_n,
  "/5\n",
  sep = ""
)


cat(
  "Nominal P < 0.05 = ",
  nominal_gene_n,
  "/5\n",
  sep = ""
)


cat(
  "BH-adjusted P < 0.05 = ",
  bh_gene_n,
  "/5\n",
  sep = ""
)


if (
  direction_concordant_n !=
  5
) {
  
  stop(
    "Expected all five component genes to show concordant direction."
  )
}


if (
  nominal_gene_n !=
  1
) {
  
  stop(
    "Expected exactly one component gene with nominal P <0.05."
  )
}


if (
  bh_gene_n !=
  0
) {
  
  stop(
    "Expected zero component genes with BH-adjusted P <0.05."
  )
}


nominal_gene <- gene_audit_df$gene[
  gene_audit_df$p_value <
    0.05
]


if (
  length(
    nominal_gene
  ) !=
  1 ||
  nominal_gene !=
  "CARD11"
) {
  
  stop(
    "Expected CARD11 to be the only nominally significant component gene."
  )
}


gene_publication <- gene_audit_df %>%
  
  dplyr::select(
    gene,
    expected_direction,
    observed_direction,
    direction_concordant,
    n_case,
    n_control,
    case_median,
    control_median,
    median_difference_case_minus_control,
    wilcoxon_W,
    p_value,
    p_BH_five_genes
  )


# =============================================================================
# 22. SCALING SENSITIVITY
# =============================================================================

if (
  nrow(
    scaling_df
  ) !=
  1
) {
  
  stop(
    "Expected exactly one frozen scaling-sensitivity result."
  )
}


if (
  abs(
    scaling_df$spearman_rho -
    0.9938685
  ) >
  1e-6
) {
  
  stop(
    "Frozen scaling-sensitivity rho does not match expected anchor."
  )
}


# =============================================================================
# 23. OPTIONAL MAIN FIGURE 4 CONSISTENCY AUDIT
# =============================================================================

figure4_consistency <- data.frame(
  
  metric = character(),
  tableS9_value = numeric(),
  figure4_value = numeric(),
  absolute_difference = numeric(),
  pass = logical(),
  stringsAsFactors = FALSE
)


if (
  figure4_audit_available
) {
  
  fig4_audit <- read.csv(
    figure4_audit_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  figure4_targets <- data.frame(
    
    metric = c(
      "Primary Wilcoxon P",
      "Primary BH q",
      "Primary AUC",
      "Primary AUC CI low",
      "Primary AUC CI high",
      "Ordered Spearman rho",
      "Ordered Spearman P",
      "Ordered Kruskal-Wallis P",
      "Shock vs uncomplicated P",
      "Shock vs uncomplicated BH q"
    ),
    
    tableS9_value = c(
      primary_comparison$p_value,
      primary_comparison$p_BH_across_score_comparisons,
      primary_comparison$auc_fixed_direction,
      primary_comparison$auc_ci_low,
      primary_comparison$auc_ci_high,
      spearman_row$statistic,
      spearman_row$p_value,
      kw_row$p_value,
      shock_comparison$p_value,
      shock_comparison$p_BH_across_score_comparisons
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  figure4_consistency <- figure4_targets %>%
    
    dplyr::left_join(
      fig4_audit[
        ,
        c(
          "metric",
          "observed"
        )
      ],
      by = "metric"
    ) %>%
    
    dplyr::rename(
      figure4_value =
        observed
    ) %>%
    
    dplyr::mutate(
      
      absolute_difference =
        abs(
          tableS9_value -
            figure4_value
        ),
      
      pass =
        is.finite(
          figure4_value
        ) &
        absolute_difference <
        1e-10
    )
  
  
  cat("\nMAIN FIGURE 4 CONSISTENCY AUDIT\n")
  cat("-------------------------------\n")
  
  print(
    figure4_consistency,
    row.names = FALSE
  )
  
  
  if (
    !all(
      figure4_consistency$pass
    )
  ) {
    
    stop(
      "Table S9 frozen values are not fully consistent with frozen Main Figure 4."
    )
  }
  
  
  cat(
    "Table S9 is numerically consistent with frozen Main Figure 4.\n"
  )
}


# =============================================================================
# 24. PUBLICATION SUMMARY TABLE
# =============================================================================

summary_table <- data.frame(
  
  Domain = c(
    "Cohort",
    "Cohort",
    "Cohort",
    "Cohort",
    "Primary external comparison",
    "Primary external comparison",
    "Primary external comparison",
    "Primary external comparison",
    "Primary external comparison",
    "Ordered baseline structure",
    "Ordered baseline structure",
    "Secondary external comparison",
    "Secondary external comparison",
    "Secondary external comparison",
    "Component-gene replication",
    "Component-gene replication",
    "Scaling sensitivity"
  ),
  
  Result = c(
    "Healthy baseline samples",
    "Uncomplicated infection baseline samples",
    "Sepsis baseline samples",
    "Septic shock baseline samples",
    
    "Sepsis/septic shock vs uncomplicated infection: Wilcoxon P",
    "Sepsis/septic shock vs uncomplicated infection: BH q",
    "Sepsis/septic shock vs uncomplicated infection: fixed-direction AUC",
    "Sepsis/septic shock vs uncomplicated infection: AUC 95% CI",
    "Sepsis/septic shock vs uncomplicated infection: median score",
    
    "Score vs ordered baseline disease state",
    "Four-group Kruskal-Wallis",
    
    "Septic shock vs uncomplicated infection: Wilcoxon P",
    "Septic shock vs uncomplicated infection: BH q",
    "Septic shock vs uncomplicated infection: fixed-direction AUC",
    
    "Component genes with expected direction",
    "Component genes with BH-adjusted P <0.05",
    
    "Cohort-standardized vs healthy-reference score"
  ),
  
  Value = c(
    "n=40",
    "n=12",
    "n=20",
    "n=19",
    
    sprintf(
      "%.6g",
      primary_comparison$p_value
    ),
    
    sprintf(
      "%.6g",
      primary_comparison$p_BH_across_score_comparisons
    ),
    
    sprintf(
      "%.3f",
      primary_comparison$auc_fixed_direction
    ),
    
    paste0(
      sprintf(
        "%.3f",
        primary_comparison$auc_ci_low
      ),
      "–",
      sprintf(
        "%.3f",
        primary_comparison$auc_ci_high
      )
    ),
    
    paste0(
      sprintf(
        "%.3f",
        primary_comparison$case_median
      ),
      " vs ",
      sprintf(
        "%.3f",
        primary_comparison$control_median
      )
    ),
    
    paste0(
      "Spearman rho=",
      sprintf(
        "%.3f",
        spearman_row$statistic
      ),
      "; P=",
      format(
        spearman_row$p_value,
        scientific = TRUE,
        digits = 3
      )
    ),
    
    paste0(
      "chi-square=",
      sprintf(
        "%.3f",
        kw_row$statistic
      ),
      "; P=",
      format(
        kw_row$p_value,
        scientific = TRUE,
        digits = 3
      )
    ),
    
    sprintf(
      "%.6g",
      shock_comparison$p_value
    ),
    
    sprintf(
      "%.6g",
      shock_comparison$p_BH_across_score_comparisons
    ),
    
    paste0(
      sprintf(
        "%.3f",
        shock_comparison$auc_fixed_direction
      ),
      " (95% CI ",
      sprintf(
        "%.3f",
        shock_comparison$auc_ci_low
      ),
      "–",
      sprintf(
        "%.3f",
        shock_comparison$auc_ci_high
      ),
      ")"
    ),
    
    "5/5",
    
    "0/5",
    
    paste0(
      "Spearman rho=",
      sprintf(
        "%.3f",
        scaling_df$spearman_rho
      )
    )
  ),
  
  Interpretation = c(
    rep(
      "Baseline cohort composition",
      4
    ),
    
    "Prespecified primary comparison; not statistically significant",
    "Prespecified primary comparison; not statistically significant",
    "Modest discrimination in the prespecified primary comparison",
    "Confidence interval includes weak-to-moderate discrimination",
    "Higher median score in sepsis/septic shock, but primary test not significant",
    
    "Strong ordered molecular gradient across baseline disease states",
    "Strong overall difference across the four baseline groups",
    
    "Secondary comparison; statistically significant after frozen BH correction",
    "Secondary comparison; statistically significant after frozen BH correction",
    "Moderate discrimination in the secondary shock comparison",
    
    "Directional replication of the complete five-gene architecture",
    "No individual component gene retained significance after five-gene BH correction",
    
    "Score ordering was highly robust to the alternative scaling reference"
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 25. READ-ME TABLE
# =============================================================================

readme_table <- data.frame(
  
  Item = c(
    "Table",
    "Dataset",
    "Data type",
    "Baseline samples",
    "Follow-up samples",
    "Follow-up use",
    "Frozen five-gene signature",
    "Score definition",
    "Primary external comparison",
    "Feature selection in external cohort",
    "Coefficient refitting in external cohort",
    "Cutoff optimization",
    "Direction flipping",
    "Primary interpretation",
    "Statistical provenance"
  ),
  
  Description = c(
    "Supplementary Table S9",
    "GSE154918",
    "Whole-blood RNA-seq; GEO processed expression",
    "91: healthy 40, uncomplicated infection 12, sepsis 20, septic shock 19",
    "14: sepsis follow-up 4, septic-shock follow-up 10",
    "Excluded from the primary baseline external analysis",
    "CD177, HK3, IRAK3, CARD11, IKZF2",
    "mean z(CD177, HK3, IRAK3) - mean z(CARD11, IKZF2)",
    "Sepsis plus septic shock versus uncomplicated infection",
    "No",
    "No",
    "No",
    "No",
    paste0(
      "Directional and ordinal external transcriptomic replication, ",
      "with limited discrimination between sepsis/septic shock and ",
      "uncomplicated infection in the prespecified primary comparison"
    ),
    paste0(
      "All statistical quantities are copied from frozen Script 141 outputs; ",
      "Script 163 performs publication packaging and auditing only."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 26. INTERPRETATION TABLE
# =============================================================================

interpretation_table <- data.frame(
  
  Topic = c(
    "Primary endpoint",
    "Ordered disease-state association",
    "Secondary shock comparison",
    "Component-gene direction",
    "Individual component-gene significance",
    "Healthy controls",
    "Clinical interpretation"
  ),
  
  Interpretation = c(
    
    paste0(
      "The prespecified sepsis/septic-shock versus uncomplicated-infection ",
      "comparison was not statistically significant and must not be described ",
      "as successful diagnostic validation."
    ),
    
    paste0(
      "The score showed a strong monotonic association with ordered baseline ",
      "clinical state, supporting replication of a disease-state molecular gradient."
    ),
    
    paste0(
      "The septic-shock versus uncomplicated-infection contrast was significant ",
      "after the frozen multiple-testing correction and had a fixed-direction AUC ",
      "of approximately 0.785."
    ),
    
    paste0(
      "All five component genes changed in the direction expected from the ",
      "discovery cohort."
    ),
    
    paste0(
      "Only CARD11 was nominally significant in the primary external gene-level ",
      "comparison, and none of the five genes remained significant after BH ",
      "correction."
    ),
    
    paste0(
      "Sepsis/septic-shock samples separated strongly from healthy controls, ",
      "but healthy controls are not the clinically challenging primary comparator."
    ),
    
    paste0(
      "GSE154918 supports directional and ordinal biological replication, ",
      "but not validated clinical diagnostic performance."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 27. CREATE SUPPLEMENTARY TABLE S9 WORKBOOK
# =============================================================================

tableS9_file <- file.path(
  tables_dir,
  "163_TableS9_GSE154918_external_evaluation.xlsx"
)


wb <- openxlsx::createWorkbook()


sheet_data <- list(
  
  "S9_ReadMe" =
    readme_table,
  
  "External_summary" =
    summary_table,
  
  "Baseline_groups" =
    group_summary_publication,
  
  "Score_comparisons" =
    comparisons_publication,
  
  "Ordered_baseline" =
    ordered_publication,
  
  "Pairwise_baseline" =
    pairwise_publication,
  
  "Gene_direction_audit" =
    gene_publication,
  
  "Scaling_sensitivity" =
    scaling_df,
  
  "Interpretation" =
    interpretation_table
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE,
  valign = "top"
)


wrap_style <- openxlsx::createStyle(
  wrapText = TRUE,
  valign = "top"
)


for (
  sheet_name in names(
    sheet_data
  )
) {
  
  data_to_write <- sheet_data[[sheet_name]]
  
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  
  openxlsx::writeData(
    wb,
    sheet_name,
    data_to_write,
    withFilter = nrow(data_to_write) > 1
  )
  
  
  openxlsx::addStyle(
    wb,
    sheet_name,
    header_style,
    rows = 1,
    cols = seq_len(
      ncol(
        data_to_write
      )
    ),
    gridExpand = TRUE
  )
  
  
  if (
    nrow(
      data_to_write
    ) >
    0
  ) {
    
    openxlsx::addStyle(
      wb,
      sheet_name,
      wrap_style,
      rows = 2:(nrow(data_to_write) + 1),
      cols = seq_len(
        ncol(
          data_to_write
        )
      ),
      gridExpand = TRUE
    )
  }
  
  
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
        data_to_write
      )
    ),
    widths = "auto"
  )
}


openxlsx::saveWorkbook(
  wb,
  tableS9_file,
  overwrite = TRUE
)


# =============================================================================
# 28. INTERNAL AUDIT WORKBOOK
# =============================================================================

audit_file <- file.path(
  audit_dir,
  "163_INTERNAL_AUDIT_TableS9_GSE154918.xlsx"
)


audit_wb <- openxlsx::createWorkbook()


audit_sheets <- list(
  
  "Schema_audit" =
    schema_checks,
  
  "CSV_XLSX_equivalence" =
    equivalence_audit,
  
  "Run_info_guardrails" =
    guardrail_audit,
  
  "Status_anchor_audit" =
    status_anchor_audit,
  
  "Primary_anchor_audit" =
    primary_anchor_audit,
  
  "Shock_anchor_audit" =
    shock_anchor_audit,
  
  "Ordered_anchor_audit" =
    ordered_anchor_audit,
  
  "MainFig4_consistency" =
    figure4_consistency,
  
  "Frozen_status_source" =
    status_df,
  
  "Frozen_group_source" =
    group_summary_df,
  
  "Frozen_comparison_source" =
    comparisons_df,
  
  "Frozen_ordered_source" =
    ordered_df,
  
  "Frozen_pairwise_source" =
    pairwise_df,
  
  "Frozen_gene_source" =
    gene_audit_df,
  
  "Frozen_scaling_source" =
    scaling_df
)


for (
  sheet_name in names(
    audit_sheets
  )
) {
  
  data_to_write <- audit_sheets[[sheet_name]]
  
  
  openxlsx::addWorksheet(
    audit_wb,
    sheet_name
  )
  
  
  openxlsx::writeData(
    audit_wb,
    sheet_name,
    data_to_write
  )
  
  
  if (
    ncol(
      data_to_write
    ) >
    0
  ) {
    
    openxlsx::addStyle(
      audit_wb,
      sheet_name,
      header_style,
      rows = 1,
      cols = seq_len(
        ncol(
          data_to_write
        )
      ),
      gridExpand = TRUE
    )
    
    
    openxlsx::setColWidths(
      audit_wb,
      sheet_name,
      cols = seq_len(
        ncol(
          data_to_write
        )
      ),
      widths = "auto"
    )
  }
}


openxlsx::saveWorkbook(
  audit_wb,
  audit_file,
  overwrite = TRUE
)


# =============================================================================
# 29. TABLE TITLE AND NOTE
# =============================================================================

table_title_note <- paste0(
  
  "Supplementary Table S9. External evaluation of the frozen five-gene ",
  "host-response score in GSE154918.\n\n",
  
  "The external baseline cohort comprised 91 whole-blood RNA-seq samples: ",
  "40 healthy controls, 12 patients with uncomplicated infection, 20 with ",
  "sepsis, and 19 with septic shock. Four sepsis follow-up samples and ten ",
  "septic-shock follow-up samples were excluded from the primary baseline ",
  "analysis. The five-gene score comprised CD177, HK3, and IRAK3 in the ",
  "positive component and CARD11 and IKZF2 in the negative component. No ",
  "external feature selection, coefficient refitting, cutoff optimization, ",
  "or score-direction reversal was performed. The prespecified primary ",
  "comparison was sepsis plus septic shock versus uncomplicated infection. ",
  "All statistical quantities shown in this table were taken directly from ",
  "the frozen Script 141 external-evaluation outputs and were not recalculated ",
  "during preparation of the supplementary table. BH, Benjamini-Hochberg; ",
  "AUC, area under the receiver-operating-characteristic curve; CI, confidence ",
  "interval. The GSE154918 analysis should be interpreted as independent ",
  "external transcriptomic evaluation rather than as validation of a ",
  "calibrated clinical diagnostic assay."
)


title_note_file <- file.path(
  text_dir,
  "163_TableS9_title_and_note_EN.txt"
)


writeLines(
  table_title_note,
  title_note_file
)


# =============================================================================
# 30. PROPOSED RESULTS 3.8
# =============================================================================

results_3_8 <- paste0(
  
  "The frozen five-gene host-response score was next evaluated in the ",
  "independent whole-blood RNA-seq cohort GSE154918. The baseline analysis ",
  "included 91 samples comprising 40 healthy controls, 12 patients with ",
  "uncomplicated infection, 20 with sepsis, and 19 with septic shock; ",
  "follow-up samples were excluded from the primary analysis. Across the ",
  "four ordered baseline states, the score showed a strong monotonic ",
  "association with increasing disease state (Spearman rho=",
  sprintf(
    "%.3f",
    spearman_row$statistic
  ),
  ", P=",
  format(
    spearman_row$p_value,
    scientific = TRUE,
    digits = 3
  ),
  "), with a corresponding overall difference across groups ",
  "(Kruskal-Wallis chi-square=",
  sprintf(
    "%.2f",
    kw_row$statistic
  ),
  ", P=",
  format(
    kw_row$p_value,
    scientific = TRUE,
    digits = 3
  ),
  "; Fig. 4; Supplementary Table S9). ",
  
  "However, the prespecified primary comparison between sepsis/septic shock ",
  "and uncomplicated infection was not statistically significant: median ",
  "scores were ",
  sprintf(
    "%.3f",
    primary_comparison$case_median
  ),
  " and ",
  sprintf(
    "%.3f",
    primary_comparison$control_median
  ),
  ", respectively (Wilcoxon P=",
  sprintf(
    "%.3f",
    primary_comparison$p_value
  ),
  "; BH-adjusted q=",
  sprintf(
    "%.3f",
    primary_comparison$p_BH_across_score_comparisons
  ),
  "), with a fixed-direction AUC of ",
  sprintf(
    "%.3f",
    primary_comparison$auc_fixed_direction
  ),
  " (95% CI ",
  sprintf(
    "%.3f",
    primary_comparison$auc_ci_low
  ),
  "-",
  sprintf(
    "%.3f",
    primary_comparison$auc_ci_high
  ),
  "). ",
  
  "In a secondary comparison, septic shock was more clearly separated from ",
  "uncomplicated infection (Wilcoxon P=",
  sprintf(
    "%.4f",
    shock_comparison$p_value
  ),
  "; BH-adjusted q=",
  sprintf(
    "%.4f",
    shock_comparison$p_BH_across_score_comparisons
  ),
  "; AUC=",
  sprintf(
    "%.3f",
    shock_comparison$auc_fixed_direction
  ),
  ", 95% CI ",
  sprintf(
    "%.3f",
    shock_comparison$auc_ci_low
  ),
  "-",
  sprintf(
    "%.3f",
    shock_comparison$auc_ci_high
  ),
  "). ",
  
  "All five component genes changed in the direction expected from the ",
  "discovery cohort, although only CARD11 reached nominal significance in ",
  "the primary external gene-level comparison and none remained significant ",
  "after Benjamini-Hochberg correction across the five genes. These findings ",
  "therefore support directional and ordinal replication of the host-response ",
  "axis in GSE154918, while also demonstrating limited discrimination between ",
  "sepsis/septic shock and uncomplicated infection in the prespecified ",
  "primary comparison."
)


results_file <- file.path(
  text_dir,
  "163_proposed_Results_3.8_GSE154918_external_evaluation_EN.txt"
)


writeLines(
  results_3_8,
  results_file
)


# =============================================================================
# 31. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "163_sessionInfo.txt"
  )
)


# =============================================================================
# 32. FINAL REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 163 FINAL v2 completed successfully.\n")
cat("====================================================================\n\n")


cat("GSE154918 COHORT\n")
cat("-----------------\n")

cat(
  "Baseline samples = ",
  baseline_n,
  "\n",
  sep = ""
)

cat(
  "Healthy = 40\n"
)

cat(
  "Uncomplicated infection = 12\n"
)

cat(
  "Sepsis = 20\n"
)

cat(
  "Septic shock = 19\n"
)

cat(
  "Follow-up samples excluded from primary baseline analysis = ",
  followup_n,
  "\n",
  sep = ""
)


cat("\nPRIMARY EXTERNAL COMPARISON\n")
cat("---------------------------\n")

cat(
  "Comparison: sepsis + septic shock vs uncomplicated infection\n"
)

cat(
  "n = ",
  primary_comparison$n_case,
  " vs ",
  primary_comparison$n_control,
  "\n",
  sep = ""
)

cat(
  "Median score = ",
  primary_comparison$case_median,
  " vs ",
  primary_comparison$control_median,
  "\n",
  sep = ""
)

cat(
  "Wilcoxon P = ",
  primary_comparison$p_value,
  "\n",
  sep = ""
)

cat(
  "BH q = ",
  primary_comparison$p_BH_across_score_comparisons,
  "\n",
  sep = ""
)

cat(
  "Fixed-direction AUC = ",
  primary_comparison$auc_fixed_direction,
  "\n",
  sep = ""
)

cat(
  "95% CI = ",
  primary_comparison$auc_ci_low,
  " to ",
  primary_comparison$auc_ci_high,
  "\n",
  sep = ""
)


cat("\nORDERED BASELINE STRUCTURE\n")
cat("--------------------------\n")

cat(
  "Spearman rho = ",
  spearman_row$statistic,
  "\n",
  sep = ""
)

cat(
  "Spearman P = ",
  spearman_row$p_value,
  "\n",
  sep = ""
)

cat(
  "Kruskal-Wallis chi-square = ",
  kw_row$statistic,
  "\n",
  sep = ""
)

cat(
  "Kruskal-Wallis P = ",
  kw_row$p_value,
  "\n",
  sep = ""
)


cat("\nSECONDARY SHOCK COMPARISON\n")
cat("--------------------------\n")

cat(
  "Shock vs uncomplicated infection AUC = ",
  shock_comparison$auc_fixed_direction,
  "\n",
  sep = ""
)

cat(
  "95% CI = ",
  shock_comparison$auc_ci_low,
  " to ",
  shock_comparison$auc_ci_high,
  "\n",
  sep = ""
)

cat(
  "P = ",
  shock_comparison$p_value,
  "\n",
  sep = ""
)

cat(
  "BH q = ",
  shock_comparison$p_BH_across_score_comparisons,
  "\n",
  sep = ""
)


cat("\nCOMPONENT-GENE REPLICATION\n")
cat("--------------------------\n")

cat(
  "Expected direction = ",
  direction_concordant_n,
  "/5\n",
  sep = ""
)

cat(
  "Nominal P < 0.05 = ",
  nominal_gene_n,
  "/5 (",
  nominal_gene,
  ")\n",
  sep = ""
)

cat(
  "BH-adjusted P < 0.05 = ",
  bh_gene_n,
  "/5\n",
  sep = ""
)


cat("\nSCALING SENSITIVITY\n")
cat("-------------------\n")

cat(
  "Cohort-standardized vs healthy-reference rho = ",
  scaling_df$spearman_rho,
  "\n",
  sep = ""
)


cat("\nOUTPUT FILES\n")
cat("------------\n")

cat(
  "Supplementary Table S9:\n",
  normalizePath(
    tableS9_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Internal audit:\n",
  normalizePath(
    audit_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Table title/note:\n",
  normalizePath(
    title_note_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Proposed Results 3.8:\n",
  normalizePath(
    results_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")

cat(
  "- All statistical quantities come directly from frozen Script 141 outputs.\n"
)

cat(
  "- No statistical test, ROC analysis, confidence interval, or BH value is recalculated.\n"
)

cat(
  "- Follow-up samples are excluded from the primary baseline analysis.\n"
)

cat(
  "- The prespecified primary comparison is formally negative.\n"
)

cat(
  "- AUC 0.656 must not be presented as successful diagnostic validation.\n"
)

cat(
  "- The ordered disease-state association provides strong external molecular replication.\n"
)

cat(
  "- The shock-versus-uncomplicated result is a secondary comparison.\n"
)

cat(
  "- Five of five genes replicate the expected direction, but zero of five are BH-significant individually.\n"
)

cat(
  "- Healthy-control discrimination is contextual and is not the primary clinically relevant test.\n"
)

cat(
  "- GSE154918 supports directional and ordinal biological replication, not a calibrated clinical assay.\n"
)

cat("\nDone.\n")