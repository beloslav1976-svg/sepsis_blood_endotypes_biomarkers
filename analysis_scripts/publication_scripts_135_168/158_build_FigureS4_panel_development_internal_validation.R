################################################################################
# Script 158
# FINAL
#
# Supplementary Figure S4
#
# Five-gene panel development and internal validation
#
# Project:
#   Sepsis_DESeq2
#
# Manuscript:
#   Blood-only sepsis transcriptomic endotypes /
#   five-gene host-response signature
#
#
# FIGURE S4
# ---------
#
# A. Biology-guided 13-gene candidate pool
#
# B. Exhaustive eligible-panel search by panel size
#
# C. Composition of the primary IKZF2-containing and alternative
#    DCAF17-containing five-gene signatures
#
# D. Saturation of the five-gene exhaustive search
#
# E. 100 x stratified five-fold internal cross-validation summary
#
#
# IMPORTANT INTERPRETATION
# ------------------------
#
# The exhaustive search demonstrated PERFORMANCE SATURATION.
#
# Therefore:
#
#   - the primary five-gene signature must NOT be described as a
#     uniquely optimal diagnostic classifier;
#
#   - rank/order among tied AUC=1 panels must NOT be interpreted as
#     evidence of biological or statistical superiority;
#
#   - repeated internal CV assesses stability within the same discovery
#     cohort and is NOT independent external validation.
#
#
# INPUTS
# ------
#
# Frozen Supplementary Table S5:
#
#   results/blood_endotypes_biomarkers/
#   156_TableS5_candidate_panel_screening/tables/
#   156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx
#
#
# Frozen Supplementary Table S6:
#
#   results/blood_endotypes_biomarkers/
#   157_TableS6_internal_cross_validation/tables/
#   157_TableS6_repeated_internal_cross_validation.xlsx
#
#
# THIS SCRIPT DOES NOT:
#
#   - rerun feature selection;
#   - rerun exhaustive panel screening;
#   - rerun LOOCV;
#   - rerun repeated cross-validation;
#   - rerun ridge regression;
#   - use SRS/CTS for feature selection;
#   - use external cohorts for feature selection.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 158\n")
cat("Supplementary Figure S4\n")
cat("Five-gene panel development and internal validation\n")
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
  "158_FigureS4_panel_development_internal_validation"
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
  directory_name in c(
    output_dir,
    figures_dir,
    tables_dir,
    text_dir,
    audit_dir
  )
) {
  
  dir.create(
    directory_name,
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
# 4. INPUT FILES
# =============================================================================

tableS5_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "156_TableS5_candidate_panel_screening",
  "tables",
  "156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx"
)


tableS6_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "157_TableS6_internal_cross_validation",
  "tables",
  "157_TableS6_repeated_internal_cross_validation.xlsx"
)


if (
  !file.exists(
    tableS5_file
  )
) {
  
  stop(
    paste0(
      "Frozen Table S5 not found:\n",
      tableS5_file
    )
  )
}


if (
  !file.exists(
    tableS6_file
  )
) {
  
  stop(
    paste0(
      "Frozen Table S6 not found:\n",
      tableS6_file
    )
  )
}


cat("\nINPUT FILES\n")
cat("-----------\n")


cat(
  "Frozen Table S5:\n  ",
  normalizePath(
    tableS5_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Frozen Table S6:\n  ",
  normalizePath(
    tableS6_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 5. HELPER FUNCTIONS
# =============================================================================

normalize_gene <- function(x) {
  
  toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
}


safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      x
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
    names(
      data
    )
  )
  
  
  if (
    length(
      missing_columns
    ) >
    0
  ) {
    
    cat(
      "\nAvailable columns in ",
      object_name,
      ":\n",
      sep = ""
    )
    
    
    print(
      names(
        data
      )
    )
    
    
    stop(
      paste0(
        "Missing column(s) in ",
        object_name,
        ": ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
}


# =============================================================================
# 6. READ TABLE S5
# =============================================================================

s5_sheets <- readxl::excel_sheets(
  tableS5_file
)


cat("\nTable S5 sheets:\n")

print(
  s5_sheets
)


required_s5_sheets <- c(
  "Candidate_pool",
  "Panel_size_summary",
  "Manuscript_panels",
  "Manuscript_panel_genes",
  "All_eligible_panels"
)


missing_s5 <- setdiff(
  required_s5_sheets,
  s5_sheets
)


if (
  length(
    missing_s5
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing Table S5 sheet(s): ",
      paste(
        missing_s5,
        collapse = ", "
      )
    )
  )
}


candidate_pool <- readxl::read_excel(
  tableS5_file,
  sheet = "Candidate_pool"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


panel_size_summary <- readxl::read_excel(
  tableS5_file,
  sheet = "Panel_size_summary"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


manuscript_panels <- readxl::read_excel(
  tableS5_file,
  sheet = "Manuscript_panels"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


manuscript_panel_genes <- readxl::read_excel(
  tableS5_file,
  sheet = "Manuscript_panel_genes"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


all_eligible_panels <- readxl::read_excel(
  tableS5_file,
  sheet = "All_eligible_panels"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


# =============================================================================
# 7. READ TABLE S6
# =============================================================================

s6_sheets <- readxl::excel_sheets(
  tableS6_file
)


cat("\nTable S6 sheets:\n")

print(
  s6_sheets
)


required_s6_sheets <- c(
  "CV_summary",
  "CV_source_complete"
)


missing_s6 <- setdiff(
  required_s6_sheets,
  s6_sheets
)


if (
  length(
    missing_s6
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing Table S6 sheet(s): ",
      paste(
        missing_s6,
        collapse = ", "
      )
    )
  )
}


cv_summary <- readxl::read_excel(
  tableS6_file,
  sheet = "CV_summary"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cv_source_complete <- readxl::read_excel(
  tableS6_file,
  sheet = "CV_source_complete"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


# =============================================================================
# 8. INPUT AUDIT
# =============================================================================

check_columns(
  candidate_pool,
  c(
    "gene",
    "direction"
  ),
  "Table S5 Candidate_pool"
)


check_columns(
  panel_size_summary,
  c(
    "search",
    "n_genes",
    "n_panels",
    "n_LOOCV_AUC_ge_0_95",
    "n_LOOCV_AUC_eq_1"
  ),
  "Table S5 Panel_size_summary"
)


check_columns(
  manuscript_panels,
  c(
    "manuscript_panel",
    "genes",
    "apparent_auc_oriented",
    "loocv_auc_oriented",
    "loocv_p_wilcox"
  ),
  "Table S5 Manuscript_panels"
)


check_columns(
  manuscript_panel_genes,
  c(
    "manuscript_panel",
    "gene",
    "expected_direction"
  ),
  "Table S5 Manuscript_panel_genes"
)


check_columns(
  all_eligible_panels,
  c(
    "genes",
    "n_genes",
    "loocv_auc_oriented"
  ),
  "Table S5 All_eligible_panels"
)


check_columns(
  cv_summary,
  c(
    "panel_original",
    "panel",
    "mean_AUC",
    "median_AUC",
    "q025_AUC",
    "q975_AUC",
    "n_repeats"
  ),
  "Table S6 CV_summary"
)


# =============================================================================
# 9. HARD FROZEN ANCHORS
# =============================================================================

if (
  nrow(
    candidate_pool
  ) !=
  13
) {
  
  stop(
    paste0(
      "Expected 13 blood candidate genes; observed ",
      nrow(
        candidate_pool
      )
    )
  )
}


if (
  nrow(
    all_eligible_panels
  ) !=
  5432
) {
  
  stop(
    paste0(
      "Expected 5,432 eligible panels; observed ",
      nrow(
        all_eligible_panels
      )
    )
  )
}


if (
  nrow(
    manuscript_panels
  ) !=
  2
) {
  
  stop(
    "Expected exactly two manuscript five-gene configurations."
  )
}


if (
  nrow(
    cv_summary
  ) !=
  5
) {
  
  stop(
    paste0(
      "Expected five repeated-CV model rows; observed ",
      nrow(
        cv_summary
      )
    )
  )
}


expected_panel_sizes <- data.frame(
  
  n_genes = c(
    5,
    6,
    7,
    8
  ),
  
  expected_n = c(
    945,
    1540,
    1666,
    1281
  )
)


observed_panel_sizes <- panel_size_summary %>%
  
  dplyr::filter(
    search ==
      "All eligible blood panels"
  ) %>%
  
  dplyr::select(
    n_genes,
    observed_n = n_panels
  ) %>%
  
  dplyr::mutate(
    n_genes =
      as.numeric(
        n_genes
      ),
    observed_n =
      as.numeric(
        observed_n
      )
  )


panel_size_audit <- dplyr::left_join(
  expected_panel_sizes,
  observed_panel_sizes,
  by = "n_genes"
) %>%
  
  dplyr::mutate(
    match =
      expected_n ==
      observed_n
  )


if (
  !all(
    panel_size_audit$match
  )
) {
  
  print(
    panel_size_audit
  )
  
  stop(
    "Eligible-panel size distribution audit failed."
  )
}


# =============================================================================
# 10. PANEL A DATA — CANDIDATE POOL
# =============================================================================

candidate_A <- candidate_pool %>%
  
  dplyr::mutate(
    
    gene =
      normalize_gene(
        gene
      ),
    
    direction =
      toupper(
        trimws(
          direction
        )
      ),
    
    direction_label =
      dplyr::case_when(
        direction ==
          "UP" ~
          "UP in sepsis",
        direction ==
          "DOWN" ~
          "DOWN in sepsis",
        TRUE ~
          direction
      )
  )


primary_genes <- manuscript_panel_genes %>%
  
  dplyr::filter(
    manuscript_panel ==
      "Primary_5_gene"
  ) %>%
  
  dplyr::pull(
    gene
  ) %>%
  
  normalize_gene()


alternative_genes <- manuscript_panel_genes %>%
  
  dplyr::filter(
    manuscript_panel ==
      "DCAF17_5_gene"
  ) %>%
  
  dplyr::pull(
    gene
  ) %>%
  
  normalize_gene()


candidate_A <- candidate_A %>%
  
  dplyr::mutate(
    
    membership =
      dplyr::case_when(
        
        gene %in%
          primary_genes &
          gene %in%
          alternative_genes ~
          "Shared selected gene",
        
        gene %in%
          primary_genes ~
          "Primary-specific",
        
        gene %in%
          alternative_genes ~
          "Alternative-specific",
        
        TRUE ~
          "Candidate not retained"
      ),
    
    direction_x =
      dplyr::case_when(
        direction_label ==
          "DOWN in sepsis" ~
          -1,
        direction_label ==
          "UP in sepsis" ~
          1,
        TRUE ~
          0
      )
  )


candidate_A$gene <- factor(
  candidate_A$gene,
  levels = rev(
    c(
      "CD177",
      "HK3",
      "IRAK3",
      "PFKFB3",
      "S100A12",
      "MMP9",
      "CARD11",
      "IKZF2",
      "NR1D2",
      "P2RY10",
      "RPS6",
      "ST6GAL1",
      "DCAF17"
    )
  )
)


# =============================================================================
# 11. PANEL B DATA — PANEL-SIZE SCREENING
# =============================================================================

panel_B_data <- panel_size_summary %>%
  
  dplyr::filter(
    search ==
      "All eligible blood panels"
  ) %>%
  
  dplyr::mutate(
    
    n_genes =
      as.numeric(
        n_genes
      ),
    
    n_panels =
      as.numeric(
        n_panels
      ),
    
    n_LOOCV_AUC_ge_0_95 =
      as.numeric(
        n_LOOCV_AUC_ge_0_95
      ),
    
    n_LOOCV_AUC_eq_1 =
      as.numeric(
        n_LOOCV_AUC_eq_1
      ),
    
    percent_AUC_ge_0_95 =
      100 *
      n_LOOCV_AUC_ge_0_95 /
      n_panels,
    
    percent_AUC_eq_1 =
      100 *
      n_LOOCV_AUC_eq_1 /
      n_panels
  ) %>%
  
  dplyr::arrange(
    n_genes
  )


overall_perfect <- sum(
  panel_B_data$n_LOOCV_AUC_eq_1
)


overall_panels <- sum(
  panel_B_data$n_panels
)


overall_perfect_percent <-
  100 *
  overall_perfect /
  overall_panels


cat("\nEXHAUSTIVE SEARCH SATURATION\n")
cat("----------------------------\n")


print(
  panel_B_data[
    ,
    c(
      "n_genes",
      "n_panels",
      "n_LOOCV_AUC_ge_0_95",
      "n_LOOCV_AUC_eq_1",
      "percent_AUC_ge_0_95",
      "percent_AUC_eq_1"
    )
  ],
  row.names = FALSE
)


cat(
  "Perfect AUC=1 panels overall = ",
  overall_perfect,
  "/",
  overall_panels,
  " (",
  sprintf(
    "%.2f",
    overall_perfect_percent
  ),
  "%)\n",
  sep = ""
)


# =============================================================================
# 12. PANEL C DATA — SIGNATURE COMPOSITION
# =============================================================================

panel_C_data <- manuscript_panel_genes %>%
  
  dplyr::mutate(
    
    gene =
      normalize_gene(
        gene
      ),
    
    panel_label =
      dplyr::case_when(
        manuscript_panel ==
          "Primary_5_gene" ~
          "Primary: IKZF2",
        manuscript_panel ==
          "DCAF17_5_gene" ~
          "Alternative: DCAF17",
        TRUE ~
          manuscript_panel
      ),
    
    direction_label =
      dplyr::case_when(
        grepl(
          "UP",
          expected_direction,
          ignore.case = TRUE
        ) ~
          "UP in sepsis",
        grepl(
          "DOWN",
          expected_direction,
          ignore.case = TRUE
        ) ~
          "DOWN in sepsis",
        TRUE ~
          expected_direction
      )
  )


unique_panel_genes <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "CARD11",
  "IKZF2",
  "DCAF17"
)


panel_C_data$gene <- factor(
  panel_C_data$gene,
  levels = rev(
    unique_panel_genes
  )
)


panel_C_data$panel_label <- factor(
  panel_C_data$panel_label,
  levels = c(
    "Primary: IKZF2",
    "Alternative: DCAF17"
  )
)


# =============================================================================
# 13. PANEL D DATA — FIVE-GENE SATURATION
# =============================================================================

five_gene_panels <- all_eligible_panels %>%
  
  dplyr::mutate(
    
    n_genes =
      as.numeric(
        n_genes
      ),
    
    loocv_auc_oriented =
      safe_numeric(
        loocv_auc_oriented
      )
  ) %>%
  
  dplyr::filter(
    n_genes ==
      5
  )


if (
  nrow(
    five_gene_panels
  ) !=
  945
) {
  
  stop(
    paste0(
      "Expected 945 five-gene panels; observed ",
      nrow(
        five_gene_panels
      )
    )
  )
}


five_gene_perfect <- sum(
  abs(
    five_gene_panels$loocv_auc_oriented -
      1
  ) <
    1e-12,
  na.rm = TRUE
)


five_gene_high_nonperfect <- sum(
  five_gene_panels$loocv_auc_oriented >=
    0.95 &
    five_gene_panels$loocv_auc_oriented <
    1 -
    1e-12,
  na.rm = TRUE
)


five_gene_below_095 <- sum(
  five_gene_panels$loocv_auc_oriented <
    0.95,
  na.rm = TRUE
)


panel_D_data <- data.frame(
  
  Category = factor(
    c(
      "AUC = 1.00",
      "0.95 ≤ AUC < 1.00",
      "AUC < 0.95"
    ),
    levels = c(
      "AUC = 1.00",
      "0.95 ≤ AUC < 1.00",
      "AUC < 0.95"
    )
  ),
  
  n = c(
    five_gene_perfect,
    five_gene_high_nonperfect,
    five_gene_below_095
  ),
  
  stringsAsFactors = FALSE
)


panel_D_data$percent <-
  100 *
  panel_D_data$n /
  sum(
    panel_D_data$n
  )


# =============================================================================
# 14. PANEL E DATA — REPEATED CV
# =============================================================================

panel_E_data <- cv_summary %>%
  
  dplyr::mutate(
    
    mean_AUC =
      safe_numeric(
        mean_AUC
      ),
    
    median_AUC =
      safe_numeric(
        median_AUC
      ),
    
    q025_AUC =
      safe_numeric(
        q025_AUC
      ),
    
    q975_AUC =
      safe_numeric(
        q975_AUC
      ),
    
    n_repeats =
      safe_numeric(
        n_repeats
      ),
    
    method_type =
      dplyr::case_when(
        
        grepl(
          "signed_score",
          panel_original,
          ignore.case = TRUE
        ) ~
          "Fixed signed score",
        
        grepl(
          "ridge",
          panel_original,
          ignore.case = TRUE
        ) ~
          "Ridge sensitivity",
        
        TRUE ~
          "Other"
      ),
    
    display_label =
      dplyr::case_when(
        
        panel_original ==
          "Primary_5_gene_signed_score" ~
          "Primary 5-gene\nsigned score",
        
        panel_original ==
          "DCAF17_5_gene_signed_score" ~
          "DCAF17 5-gene\nsigned score",
        
        panel_original ==
          "Primary_5_gene_ridge" ~
          "Primary 5-gene\nridge",
        
        panel_original ==
          "DCAF17_5_gene_ridge" ~
          "DCAF17 5-gene\nridge",
        
        grepl(
          "SeptiCyte",
          panel_original,
          ignore.case = TRUE
        ) ~
          "SeptiCyte-related\n4-gene ridge",
        
        TRUE ~
          panel_original
      )
  )


expected_cv_names <- c(
  "Primary_5_gene_signed_score",
  "DCAF17_5_gene_signed_score",
  "Primary_5_gene_ridge",
  "DCAF17_5_gene_ridge",
  "SeptiCyte_related_4_gene_ridge"
)


if (
  !setequal(
    panel_E_data$panel_original,
    expected_cv_names
  )
) {
  
  cat(
    "\nObserved CV model names:\n"
  )
  
  print(
    panel_E_data$panel_original
  )
  
  
  stop(
    "Unexpected repeated-CV model set."
  )
}


panel_E_data$display_label <- factor(
  panel_E_data$display_label,
  levels = rev(
    c(
      "Primary 5-gene\nsigned score",
      "DCAF17 5-gene\nsigned score",
      "Primary 5-gene\nridge",
      "DCAF17 5-gene\nridge",
      "SeptiCyte-related\n4-gene ridge"
    )
  )
)


# =============================================================================
# 15. COLORS
# =============================================================================

direction_colors <- c(
  
  "UP in sepsis" =
    "#B2182B",
  
  "DOWN in sepsis" =
    "#2166AC"
)


membership_colors <- c(
  
  "Shared selected gene" =
    "#54278F",
  
  "Primary-specific" =
    "#1B7837",
  
  "Alternative-specific" =
    "#E08214",
  
  "Candidate not retained" =
    "#BDBDBD"
)


saturation_colors <- c(
  
  "AUC = 1.00" =
    "#54278F",
  
  "0.95 ≤ AUC < 1.00" =
    "#9E9AC8",
  
  "AUC < 0.95" =
    "#D9D9D9"
)


cv_colors <- c(
  
  "Fixed signed score" =
    "#54278F",
  
  "Ridge sensitivity" =
    "#4D4D4D",
  
  "Other" =
    "#969696"
)


# =============================================================================
# 16. PANEL A — CANDIDATE POOL
# =============================================================================

panel_A <- ggplot2::ggplot(
  
  candidate_A,
  
  ggplot2::aes(
    x = direction_x,
    y = gene
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    color = "#BDBDBD"
  ) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend = direction_x,
      yend = gene
    ),
    
    linewidth = 0.7,
    color = "#D9D9D9"
  ) +
  
  ggplot2::geom_point(
    
    ggplot2::aes(
      color = membership,
      shape = direction_label
    ),
    
    size = 3.7
  ) +
  
  ggplot2::scale_color_manual(
    values = membership_colors
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "UP in sepsis" = 16,
      "DOWN in sepsis" = 17
    )
  ) +
  
  ggplot2::scale_x_continuous(
    
    breaks = c(
      -1,
      1
    ),
    
    labels = c(
      "DOWN",
      "UP"
    ),
    
    limits = c(
      -1.35,
      1.35
    )
  ) +
  
  ggplot2::labs(
    title = "A  Biology-guided candidate pool",
    subtitle = "13 blood genes: 6 UP and 7 DOWN candidates",
    x = "Expected direction in sepsis",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  
  ggplot2::theme_bw(
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
        size = 9.5
      ),
    
    panel.grid.major.y =
      ggplot2::element_blank(),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom"
  )


# =============================================================================
# 17. PANEL B — EXHAUSTIVE SCREENING BY PANEL SIZE
# =============================================================================

panel_B <- ggplot2::ggplot(
  
  panel_B_data,
  
  ggplot2::aes(
    x = factor(
      n_genes
    ),
    y = percent_AUC_eq_1
  )
  
) +
  
  ggplot2::geom_col(
    width = 0.62,
    fill = "#54278F",
    alpha = 0.82
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        paste0(
          n_LOOCV_AUC_eq_1,
          "/",
          n_panels,
          "\n",
          sprintf(
            "%.1f",
            percent_AUC_eq_1
          ),
          "%"
        )
    ),
    
    vjust = -0.35,
    size = 3.3
  ) +
  
  ggplot2::geom_hline(
    yintercept = 100,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  ggplot2::scale_y_continuous(
    
    limits = c(
      0,
      105
    ),
    
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  ggplot2::labs(
    title = "B  Exhaustive panel-search saturation",
    subtitle = paste0(
      overall_perfect,
      "/",
      overall_panels,
      " eligible panels achieved LOOCV AUC = 1.00"
    ),
    x = "Panel size (genes)",
    y = "Panels with LOOCV AUC = 1.00 (%)"
  ) +
  
  ggplot2::theme_bw(
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
        size = 9.5
      )
  )


# =============================================================================
# 18. PANEL C — SELECTED SIGNATURE COMPOSITION
# =============================================================================

panel_C <- ggplot2::ggplot(
  
  panel_C_data,
  
  ggplot2::aes(
    x = panel_label,
    y = gene
  )
  
) +
  
  ggplot2::geom_tile(
    ggplot2::aes(
      fill = direction_label
    ),
    width = 0.72,
    height = 0.72,
    color = "white",
    linewidth = 0.8
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = as.character(
        gene
      )
    ),
    color = "white",
    fontface = "bold",
    size = 3.3
  ) +
  
  ggplot2::scale_fill_manual(
    values = direction_colors
  ) +
  
  ggplot2::labs(
    title = "C  Five-gene signature composition",
    subtitle = "Four genes are shared; the fifth DOWN gene differs",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  ggplot2::theme_bw(
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
        size = 9.5
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 15,
        hjust = 1
      ),
    
    axis.text.y =
      ggplot2::element_blank(),
    
    axis.ticks.y =
      ggplot2::element_blank(),
    
    panel.grid =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom"
  )


# =============================================================================
# 19. PANEL D — FIVE-GENE PERFORMANCE SATURATION
# =============================================================================

panel_D <- ggplot2::ggplot(
  
  panel_D_data,
  
  ggplot2::aes(
    x = Category,
    y = n,
    fill = Category
  )
  
) +
  
  ggplot2::geom_col(
    width = 0.64
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        paste0(
          n,
          "\n(",
          sprintf(
            "%.1f",
            percent
          ),
          "%)"
        )
    ),
    
    vjust = -0.35,
    size = 3.5
  ) +
  
  ggplot2::scale_fill_manual(
    values = saturation_colors
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.13
      )
    )
  ) +
  
  ggplot2::labs(
    title = "D  Performance saturation among five-gene panels",
    subtitle = "Internal discrimination did not define a unique optimum",
    x = NULL,
    y = "Number of eligible five-gene panels"
  ) +
  
  ggplot2::theme_bw(
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
        size = 9.5
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 15,
        hjust = 1
      ),
    
    legend.position =
      "none"
  )


# =============================================================================
# 20. PANEL E — 100 x 5-FOLD REPEATED CV
# =============================================================================

panel_E <- ggplot2::ggplot(
  
  panel_E_data,
  
  ggplot2::aes(
    x = mean_AUC,
    y = display_label,
    color = method_type
  )
  
) +
  
  ggplot2::geom_errorbarh(
    
    ggplot2::aes(
      xmin = q025_AUC,
      xmax = q975_AUC
    ),
    
    height = 0.18,
    linewidth = 0.8
  ) +
  
  ggplot2::geom_point(
    size = 3.6
  ) +
  
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  ggplot2::scale_color_manual(
    values = cv_colors
  ) +
  
  ggplot2::scale_x_continuous(
    
    limits = c(
      0.985,
      1.001
    ),
    
    breaks = c(
      0.985,
      0.990,
      0.995,
      1.000
    )
  ) +
  
  ggplot2::labs(
    title = "E  Repeated internal cross-validation",
    subtitle = "100 repeats of stratified five-fold CV; bars are empirical 2.5th–97.5th percentiles",
    x = "Mean repeat AUC",
    y = NULL,
    color = NULL
  ) +
  
  ggplot2::theme_bw(
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
        size = 9.5
      ),
    
    panel.grid.major.y =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom"
  )


# =============================================================================
# 21. COMBINE FIGURE
# =============================================================================
#
# Layout:
#
# A | B
# C | D
# E | E
#
# =============================================================================

layout_matrix <- matrix(
  c(
    1, 2,
    3, 4,
    5, 5
  ),
  nrow = 3,
  byrow = TRUE
)


figure_S4_grob <- gridExtra::arrangeGrob(
  
  panel_A,
  panel_B,
  panel_C,
  panel_D,
  panel_E,
  
  layout_matrix =
    layout_matrix,
  
  heights = c(
    1.15,
    1,
    0.85
  )
)


# =============================================================================
# 22. SAVE COMBINED FIGURE
# =============================================================================

figure_pdf <- file.path(
  figures_dir,
  "158_FigureS4_panel_development_internal_validation.pdf"
)


figure_png <- file.path(
  figures_dir,
  "158_FigureS4_panel_development_internal_validation.png"
)


figure_tiff <- file.path(
  figures_dir,
  "158_FigureS4_panel_development_internal_validation.tiff"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure_S4_grob,
  width = 15,
  height = 16,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure_S4_grob,
  width = 15,
  height = 16,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure_S4_grob,
  width = 15,
  height = 16,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 23. SAVE INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  
  A_candidate_pool =
    panel_A,
  
  B_exhaustive_search =
    panel_B,
  
  C_panel_composition =
    panel_C,
  
  D_five_gene_saturation =
    panel_D,
  
  E_repeated_CV =
    panel_E
)


individual_widths <- c(
  7.5,
  7,
  7,
  7,
  9
)


individual_heights <- c(
  7,
  5.5,
  6,
  5.5,
  6
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
        "158_FigureS4_",
        panel_name,
        ".pdf"
      )
    ),
    
    plot =
      individual_panels[[panel_name]],
    
    width =
      individual_widths[i],
    
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
# 24. SOURCE-DATA WORKBOOK
# =============================================================================

source_data_file <- file.path(
  tables_dir,
  "158_FigureS4_source_data.xlsx"
)


wb <- openxlsx::createWorkbook()


source_objects <- list(
  
  Candidate_pool_A =
    candidate_A,
  
  Panel_size_B =
    panel_B_data,
  
  Panel_composition_C =
    panel_C_data,
  
  Five_gene_saturation_D =
    panel_D_data,
  
  Repeated_CV_E =
    panel_E_data,
  
  Panel_size_audit =
    panel_size_audit
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
    cols = 1:ncol(
      data_object
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
    cols = 1:ncol(
      data_object
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
# 25. AUDIT SUMMARY
# =============================================================================

audit_summary <- data.frame(
  
  metric = c(
    
    "Candidate genes",
    
    "UP candidate genes",
    
    "DOWN candidate genes",
    
    "All eligible panels",
    
    "Eligible five-gene panels",
    
    "Five-gene panels with AUC=1",
    
    "Five-gene panels with AUC>=0.95",
    
    "All eligible panels with AUC=1",
    
    "All eligible panels with AUC=1 percent",
    
    "Repeated-CV rows",
    
    "Primary signed-score mean AUC",
    
    "Primary signed-score q025 AUC",
    
    "Primary signed-score q975 AUC",
    
    "Primary ridge mean AUC",
    
    "DCAF17 signed-score mean AUC",
    
    "DCAF17 ridge mean AUC",
    
    "SeptiCyte-related ridge mean AUC",
    
    "SeptiCyte-related ridge q025 AUC",
    
    "SeptiCyte-related ridge q975 AUC"
  ),
  
  value = c(
    
    nrow(
      candidate_A
    ),
    
    sum(
      candidate_A$direction_label ==
        "UP in sepsis"
    ),
    
    sum(
      candidate_A$direction_label ==
        "DOWN in sepsis"
    ),
    
    nrow(
      all_eligible_panels
    ),
    
    nrow(
      five_gene_panels
    ),
    
    five_gene_perfect,
    
    sum(
      five_gene_panels$loocv_auc_oriented >=
        0.95,
      na.rm = TRUE
    ),
    
    overall_perfect,
    
    overall_perfect_percent,
    
    nrow(
      panel_E_data
    ),
    
    panel_E_data$mean_AUC[
      panel_E_data$panel_original ==
        "Primary_5_gene_signed_score"
    ],
    
    panel_E_data$q025_AUC[
      panel_E_data$panel_original ==
        "Primary_5_gene_signed_score"
    ],
    
    panel_E_data$q975_AUC[
      panel_E_data$panel_original ==
        "Primary_5_gene_signed_score"
    ],
    
    panel_E_data$mean_AUC[
      panel_E_data$panel_original ==
        "Primary_5_gene_ridge"
    ],
    
    panel_E_data$mean_AUC[
      panel_E_data$panel_original ==
        "DCAF17_5_gene_signed_score"
    ],
    
    panel_E_data$mean_AUC[
      panel_E_data$panel_original ==
        "DCAF17_5_gene_ridge"
    ],
    
    panel_E_data$mean_AUC[
      panel_E_data$panel_original ==
        "SeptiCyte_related_4_gene_ridge"
    ],
    
    panel_E_data$q025_AUC[
      panel_E_data$panel_original ==
        "SeptiCyte_related_4_gene_ridge"
    ],
    
    panel_E_data$q975_AUC[
      panel_E_data$panel_original ==
        "SeptiCyte_related_4_gene_ridge"
    ]
  ),
  
  stringsAsFactors = FALSE
)


audit_file <- file.path(
  audit_dir,
  "158_FigureS4_audit.xlsx"
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
  cols = 1:ncol(
    audit_summary
  ),
  gridExpand = TRUE
)


openxlsx::setColWidths(
  wb_audit,
  "Audit_summary",
  cols = 1:ncol(
    audit_summary
  ),
  widths = "auto"
)


openxlsx::saveWorkbook(
  wb_audit,
  audit_file,
  overwrite = TRUE
)


# =============================================================================
# 26. FIGURE LEGEND
# =============================================================================

figure_legend <- paste0(
  
  "Supplementary Figure S4. Development and internal evaluation of the ",
  "five-gene blood host-response signature. ",
  
  "(A) Biology-guided blood candidate pool comprising six genes with ",
  "increased and seven genes with decreased expression in sepsis. Genes ",
  "retained in both five-gene configurations, genes specific to the primary ",
  "IKZF2-containing configuration, and the alternative DCAF17-containing ",
  "configuration are indicated. ",
  
  "(B) Exhaustive evaluation of eligible 5-8-gene blood panels. Bars show ",
  "the proportion of configurations achieving a leave-one-out ",
  "cross-validation (LOOCV) AUC of 1.00; labels report the corresponding ",
  "number of perfect-performing panels and the total number evaluated at ",
  "each panel size. ",
  
  "(C) Composition of the primary five-gene signature ",
  "(CD177, HK3, IRAK3, CARD11 and IKZF2) and the DCAF17-containing ",
  "alternative sensitivity signature. Both configurations contain three ",
  "genes increased in sepsis and two genes decreased in sepsis. ",
  
  "(D) Internal performance saturation among the 945 eligible five-gene ",
  "configurations. The large number of configurations with near-perfect or ",
  "perfect internal discrimination demonstrates that discovery-cohort ",
  "performance did not define a unique mathematical optimum. ",
  
  "(E) Internal stability assessed using 100 repeats of stratified five-fold ",
  "cross-validation. Points show mean repeat AUC and horizontal bars show ",
  "empirical 2.5th-97.5th percentiles across repeats. The fixed signed-score ",
  "branch used training-fold z-standardization without fitted coefficients; ",
  "ridge logistic regression was evaluated separately as a sensitivity ",
  "analysis. The SeptiCyte-related four-gene result represents an RNA-seq ",
  "implementation of the component genes and not the proprietary clinical ",
  "SeptiCyte score. Internal cross-validation was performed within the same ",
  "discovery cohort used for candidate-panel development and therefore ",
  "assesses internal stability rather than independent clinical validation."
)


legend_file <- file.path(
  text_dir,
  "158_FigureS4_legend_EN.txt"
)


writeLines(
  figure_legend,
  legend_file
)


# =============================================================================
# 27. PROPOSED RESULTS 3.4 TEXT
# =============================================================================

results_text <- paste0(
  
  "A biology-guided candidate pool of 13 blood genes was defined, comprising ",
  "six genes with increased and seven genes with decreased expression in ",
  "sepsis. Exhaustive screening evaluated 5,432 eligible 5-8-gene ",
  "configurations containing at least two genes from each directional arm ",
  "(Supplementary Fig. S4A-B and Supplementary Table S5). Internal ",
  "discrimination was highly saturated: ",
  overall_perfect,
  " of ",
  overall_panels,
  " eligible configurations (",
  sprintf(
    "%.1f",
    overall_perfect_percent
  ),
  "%) achieved a LOOCV AUC of 1.00. Among the ",
  nrow(
    five_gene_panels
  ),
  " eligible five-gene configurations, ",
  five_gene_perfect,
  " achieved AUC=1.00. Thus, discovery-cohort discrimination did not ",
  "identify a unique mathematical optimum. For subsequent endotype-focused ",
  "analyses, the five-gene configuration comprising CD177, HK3 and IRAK3 ",
  "with increased expression and CARD11 and IKZF2 with decreased expression ",
  "was designated as the primary biology-guided host-response signature; ",
  "the corresponding DCAF17-containing configuration was retained as an ",
  "alternative sensitivity signature (Supplementary Fig. S4C and ",
  "Supplementary Table S5). Internal stability was further evaluated using ",
  "100 repeats of stratified five-fold cross-validation. The primary signed ",
  "score and its ridge-logistic sensitivity implementation both showed mean ",
  "AUC=1.00 across repeated resampling, as did the corresponding DCAF17 ",
  "configurations, whereas the SeptiCyte-related four-gene ridge model ",
  "showed a mean AUC of ",
  sprintf(
    "%.3f",
    panel_E_data$mean_AUC[
      panel_E_data$panel_original ==
        "SeptiCyte_related_4_gene_ridge"
    ]
  ),
  " (empirical 2.5th-97.5th percentiles ",
  sprintf(
    "%.3f",
    panel_E_data$q025_AUC[
      panel_E_data$panel_original ==
        "SeptiCyte_related_4_gene_ridge"
    ]
  ),
  "-",
  sprintf(
    "%.3f",
    panel_E_data$q975_AUC[
      panel_E_data$panel_original ==
        "SeptiCyte_related_4_gene_ridge"
    ]
  ),
  "; Supplementary Fig. S4E and Supplementary Table S6). Because panel ",
  "development and internal resampling were performed within the same ",
  "discovery cohort, these estimates were interpreted as measures of ",
  "internal stability rather than independent validation."
)


results_file <- file.path(
  text_dir,
  "158_proposed_Results_3.4_panel_development_EN.txt"
)


writeLines(
  results_text,
  results_file
)


# =============================================================================
# 28. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "158_sessionInfo.txt"
  )
)


# =============================================================================
# 29. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 158 completed successfully.\n")
cat("====================================================================\n\n")


cat("CANDIDATE POOL\n")
cat("--------------\n")


cat(
  "Candidate genes = ",
  nrow(
    candidate_A
  ),
  "\n",
  sep = ""
)


cat(
  "UP = ",
  sum(
    candidate_A$direction_label ==
      "UP in sepsis"
  ),
  "\n",
  sep = ""
)


cat(
  "DOWN = ",
  sum(
    candidate_A$direction_label ==
      "DOWN in sepsis"
  ),
  "\n",
  sep = ""
)


cat("\nEXHAUSTIVE SEARCH\n")
cat("-----------------\n")


cat(
  "All eligible panels = ",
  overall_panels,
  "\n",
  sep = ""
)


cat(
  "Panels with LOOCV AUC = 1 = ",
  overall_perfect,
  " (",
  sprintf(
    "%.3f",
    overall_perfect_percent
  ),
  "%)\n",
  sep = ""
)


cat("\nBY PANEL SIZE\n")
cat("-------------\n")


print(
  panel_B_data[
    ,
    c(
      "n_genes",
      "n_panels",
      "n_LOOCV_AUC_ge_0_95",
      "n_LOOCV_AUC_eq_1",
      "percent_AUC_eq_1"
    )
  ],
  row.names = FALSE
)


cat("\nFIVE-GENE SEARCH\n")
cat("----------------\n")


cat(
  "Eligible five-gene panels = ",
  nrow(
    five_gene_panels
  ),
  "\n",
  sep = ""
)


cat(
  "AUC = 1.00 = ",
  five_gene_perfect,
  "\n",
  sep = ""
)


cat(
  "0.95 <= AUC < 1.00 = ",
  five_gene_high_nonperfect,
  "\n",
  sep = ""
)


cat(
  "AUC < 0.95 = ",
  five_gene_below_095,
  "\n",
  sep = ""
)


cat("\nREPEATED INTERNAL CV\n")
cat("--------------------\n")


print(
  panel_E_data[
    ,
    c(
      "panel_original",
      "mean_AUC",
      "median_AUC",
      "q025_AUC",
      "q975_AUC",
      "n_repeats"
    )
  ],
  row.names = FALSE
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Figure S4 PDF:\n  ",
  normalizePath(
    figure_pdf,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S4 PNG:\n  ",
  normalizePath(
    figure_png,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S4 TIFF:\n  ",
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
  "Proposed Results 3.4 text:\n  ",
  normalizePath(
    results_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nINTERPRETATION GUARDRAILS\n")
cat("-------------------------\n")


cat(
  "- Figure S4 uses frozen Tables S5 and S6 only.\n"
)


cat(
  "- No panel selection or cross-validation is rerun.\n"
)


cat(
  "- Exhaustive screening demonstrates performance saturation.\n"
)


cat(
  "- The primary five-gene signature is not presented as a unique mathematical optimum.\n"
)


cat(
  "- Search rank among tied AUC=1 configurations is not interpreted as superiority.\n"
)


cat(
  "- SRS and CTS were not used for original feature selection.\n"
)


cat(
  "- Repeated CV is internal validation only.\n"
)


cat(
  "- Empirical q2.5-q97.5 values are not formal confidence intervals.\n"
)


cat(
  "- External validation is evaluated separately in GSE154918 and GSE185263.\n"
)


cat("\nDone.\n")