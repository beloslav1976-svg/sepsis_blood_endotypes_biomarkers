################################################################################
# Script 164
# FINAL
#
# Supplementary Figure S7
#
# GSE154918 external-evaluation sensitivity and component-level details
#
# Project:
#   Sepsis_DESeq2
#
#
# FIGURE ARCHITECTURE
# -------------------
#
# A. Component-gene directional replication
#    Frozen primary external comparison:
#    sepsis + septic shock versus uncomplicated infection.
#
#    Displays the frozen median expression difference for:
#      CD177
#      HK3
#      IRAK3
#      CARD11
#      IKZF2
#
#    All five genes retain the expected direction.
#    1/5 nominal P < 0.05 (CARD11).
#    0/5 BH-adjusted P < 0.05.
#
#
# B. Pairwise baseline-state evidence
#
#    Displays the six frozen pairwise baseline comparisons using:
#
#      -log10(BH-adjusted P)
#
#    The transformation is for visualization only.
#    No P or BH value is recalculated.
#
#
# C. Score-scaling sensitivity
#
#    Scatterplot of:
#
#      cohort-standardized five-gene score
#
#    versus:
#
#      healthy-reference-scaled five-gene score
#
#    Frozen Script 141 sensitivity:
#
#      Spearman rho = 0.9938685
#      P = 6.486456e-87
#
#    Correlation is NOT recalculated in Script 164.
#
#
# IMPORTANT
# ---------
#
# This script:
#
#   - does not rerun any statistical test
#   - does not recompute the five-gene score
#   - does not recompute gene-level contrasts
#   - does not recompute BH-adjusted P values
#   - does not recompute Spearman correlation
#   - does not perform external feature selection
#   - does not refit coefficients
#   - does not optimize a cutoff
#
# Figure S7 is a visualization/package of frozen Script 141 / Table S9 results.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 164 FINAL\n")
cat("Supplementary Figure S7\n")
cat("GSE154918 sensitivity and component-level details\n")
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
  "ggplot2",
  "readxl",
  "openxlsx",
  "patchwork"
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
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(patchwork)
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "164_FigureS7_GSE154918_sensitivity_details"
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
# 4. FROZEN TABLE S9
# =============================================================================

tableS9_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "163_TableS9_GSE154918_external_evaluation",
  "tables",
  "163_TableS9_GSE154918_external_evaluation.xlsx"
)


if (!file.exists(tableS9_file)) {
  stop("Frozen Supplementary Table S9 not found.")
}


cat("\nFrozen Table S9:\n")

cat(
  normalizePath(
    tableS9_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


tableS9_sheets <- readxl::excel_sheets(
  tableS9_file
)


cat("\nTable S9 sheets:\n")

print(
  tableS9_sheets
)


required_S9_sheets <- c(
  "Gene_direction_audit",
  "Pairwise_baseline",
  "Scaling_sensitivity"
)


if (
  !all(
    required_S9_sheets %in%
    tableS9_sheets
  )
) {
  
  stop(
    paste0(
      "Frozen Table S9 lacks required sheet(s): ",
      paste(
        setdiff(
          required_S9_sheets,
          tableS9_sheets
        ),
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 5. READ FROZEN TABLE S9 DATA
# =============================================================================

gene_df <- readxl::read_excel(
  tableS9_file,
  sheet = "Gene_direction_audit"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


pairwise_df <- readxl::read_excel(
  tableS9_file,
  sheet = "Pairwise_baseline"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


scaling_df <- readxl::read_excel(
  tableS9_file,
  sheet = "Scaling_sensitivity"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


# =============================================================================
# 6. EXACT TABLE S9 SCHEMA AUDIT
# =============================================================================

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


expected_pairwise_columns <- c(
  "Group_1",
  "Group_2",
  "BH_adjusted_p"
)


expected_scaling_columns <- c(
  "comparison",
  "n",
  "spearman_rho",
  "p_value"
)


schema_audit <- data.frame(
  
  object = c(
    "Gene direction audit",
    "Pairwise baseline tests",
    "Scaling sensitivity"
  ),
  
  schema_match = c(
    identical(
      names(gene_df),
      expected_gene_columns
    ),
    identical(
      names(pairwise_df),
      expected_pairwise_columns
    ),
    identical(
      names(scaling_df),
      expected_scaling_columns
    )
  ),
  
  stringsAsFactors = FALSE
)


cat("\nFROZEN TABLE S9 SCHEMA AUDIT\n")
cat("----------------------------\n")

print(
  schema_audit,
  row.names = FALSE
)


if (
  !all(
    schema_audit$schema_match
  )
) {
  
  stop(
    "At least one frozen Table S9 sheet has an unexpected schema."
  )
}


# =============================================================================
# 7. FROZEN SCRIPT 141 SAMPLE-LEVEL SCORE SOURCE
# =============================================================================

score_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "141_external_validation_GSE154918",
  "tables",
  "141_GSE154918_five_gene_scores.csv"
)


source_workbook <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "141_external_validation_GSE154918",
  "tables",
  "141_GSE154918_external_validation.xlsx"
)


if (!file.exists(score_file)) {
  stop("Frozen Script 141 sample-level score CSV not found.")
}


if (!file.exists(source_workbook)) {
  stop("Frozen Script 141 workbook not found.")
}


score_df <- read.csv(
  score_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


score_xlsx <- readxl::read_excel(
  source_workbook,
  sheet = "04_sample_scores"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat("\nFrozen sample-score dimensions:\n")

cat(
  "CSV = ",
  nrow(score_df),
  " x ",
  ncol(score_df),
  "\n",
  sep = ""
)

cat(
  "XLSX = ",
  nrow(score_xlsx),
  " x ",
  ncol(score_xlsx),
  "\n",
  sep = ""
)


# =============================================================================
# 8. SAMPLE-SCORE CSV-vs-XLSX EQUIVALENCE
# =============================================================================

normalize_text <- function(x) {
  
  out <- trimws(
    as.character(x)
  )
  
  out[out == ""] <- NA_character_
  
  out
}


column_equivalent <- function(a, b) {
  
  if (length(a) != length(b)) {
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
    
    
    if (
      !any(
        both_present
      )
    ) {
      return(TRUE)
    }
    
    
    return(
      all(
        abs(
          a_num[both_present] -
            b_num[both_present]
        ) <
          1e-10
      )
    )
  }
  
  
  a_text <- normalize_text(a)
  b_text <- normalize_text(b)
  
  
  all(
    is.na(a_text) ==
      is.na(b_text) &
      (
        is.na(a_text) |
          a_text ==
          b_text
      )
  )
}


score_dimensions_match <- identical(
  dim(score_df),
  dim(score_xlsx)
)


score_names_match <- identical(
  names(score_df),
  names(score_xlsx)
)


score_content_match <- FALSE


if (
  score_dimensions_match &&
  score_names_match
) {
  
  score_column_matches <- vapply(
    names(score_df),
    function(column_name) {
      
      column_equivalent(
        score_df[[column_name]],
        score_xlsx[[column_name]]
      )
    },
    logical(1)
  )
  
  
  score_content_match <- all(
    score_column_matches
  )
}


score_equivalence_audit <- data.frame(
  
  audit = c(
    "Dimensions identical",
    "Column names/order identical",
    "Content equivalent"
  ),
  
  result = c(
    score_dimensions_match,
    score_names_match,
    score_content_match
  ),
  
  stringsAsFactors = FALSE
)


cat("\nSAMPLE-SCORE SOURCE EQUIVALENCE AUDIT\n")
cat("-------------------------------------\n")

print(
  score_equivalence_audit,
  row.names = FALSE
)


if (
  !all(
    score_equivalence_audit$result
  )
) {
  
  stop(
    "Frozen Script 141 sample-score CSV and workbook copy are not equivalent."
  )
}


# =============================================================================
# 9. SAMPLE-LEVEL STRUCTURE AUDIT
# =============================================================================

required_score_columns <- c(
  "geo_accession",
  "status",
  "five_gene_score",
  "five_gene_score_healthy_reference"
)


missing_score_columns <- setdiff(
  required_score_columns,
  names(score_df)
)


if (
  length(
    missing_score_columns
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing required sample-score columns: ",
      paste(
        missing_score_columns,
        collapse = ", "
      )
    )
  )
}


if (
  nrow(
    score_df
  ) !=
  91
) {
  
  stop(
    paste0(
      "Expected 91 frozen baseline sample scores; observed ",
      nrow(score_df),
      "."
    )
  )
}


expected_status_counts <- data.frame(
  
  status = c(
    "Hlty",
    "Inf1_P",
    "Seps_P",
    "Shock_P"
  ),
  
  expected_n = c(
    40,
    12,
    20,
    19
  ),
  
  stringsAsFactors = FALSE
)


observed_status_counts <- score_df %>%
  
  dplyr::count(
    status,
    name = "observed_n"
  )


status_audit <- expected_status_counts %>%
  
  dplyr::left_join(
    observed_status_counts,
    by = "status"
  ) %>%
  
  dplyr::mutate(
    pass =
      expected_n ==
      observed_n
  )


cat("\nBASELINE SAMPLE-SCORE STATUS AUDIT\n")
cat("----------------------------------\n")

print(
  status_audit,
  row.names = FALSE
)


if (
  !all(
    status_audit$pass
  )
) {
  
  stop(
    "Frozen baseline sample-score status counts do not match expected GSE154918 composition."
  )
}


# =============================================================================
# 10. GENE-LEVEL FROZEN ANCHOR AUDIT
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
    gene_df
  ) !=
  5 ||
  !identical(
    as.character(
      gene_df$gene
    ),
    expected_genes
  )
) {
  
  stop(
    "Frozen Table S9 component-gene order/content is not the expected five-gene signature."
  )
}


if (
  !all(
    gene_df$direction_concordant
  )
) {
  
  stop(
    "Expected all five frozen component genes to be directionally concordant."
  )
}


n_nominal <- sum(
  gene_df$p_value <
    0.05,
  na.rm = TRUE
)


n_BH <- sum(
  gene_df$p_BH_five_genes <
    0.05,
  na.rm = TRUE
)


if (
  n_nominal !=
  1
) {
  
  stop(
    "Expected exactly one nominally significant component gene."
  )
}


if (
  n_BH !=
  0
) {
  
  stop(
    "Expected zero BH-significant component genes."
  )
}


nominal_gene <- gene_df$gene[
  gene_df$p_value <
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


cat("\nCOMPONENT-GENE FROZEN AUDIT\n")
cat("---------------------------\n")

cat(
  "Directionally concordant = ",
  sum(gene_df$direction_concordant),
  "/5\n",
  sep = ""
)

cat(
  "Nominal P < 0.05 = ",
  n_nominal,
  "/5 (",
  nominal_gene,
  ")\n",
  sep = ""
)

cat(
  "BH-adjusted P < 0.05 = ",
  n_BH,
  "/5\n",
  sep = ""
)


# =============================================================================
# 11. PAIRWISE BASELINE FROZEN AUDIT
# =============================================================================

if (
  nrow(
    pairwise_df
  ) !=
  6
) {
  
  stop(
    "Expected six frozen pairwise baseline comparisons."
  )
}


if (
  any(
    !is.finite(
      pairwise_df$BH_adjusted_p
    )
  ) ||
  any(
    pairwise_df$BH_adjusted_p <=
    0
  ) ||
  any(
    pairwise_df$BH_adjusted_p >
    1
  )
) {
  
  stop(
    "Invalid frozen pairwise BH-adjusted P value."
  )
}


# =============================================================================
# 12. SCALING FROZEN ANCHOR AUDIT
# =============================================================================

if (
  nrow(
    scaling_df
  ) !=
  1
) {
  
  stop(
    "Expected one frozen scaling-sensitivity row."
  )
}


expected_scaling_rho <- 0.9938685

expected_scaling_p <- 6.486456e-87


rho_difference <- abs(
  scaling_df$spearman_rho -
    expected_scaling_rho
)


p_difference <- abs(
  scaling_df$p_value -
    expected_scaling_p
)


if (
  rho_difference >
  1e-6
) {
  
  stop(
    "Frozen scaling-sensitivity rho failed expected anchor."
  )
}


if (
  p_difference >
  1e-92
) {
  
  stop(
    "Frozen scaling-sensitivity P value failed expected anchor."
  )
}


cat("\nSCALING-SENSITIVITY FROZEN AUDIT\n")
cat("--------------------------------\n")

cat(
  "Spearman rho = ",
  scaling_df$spearman_rho,
  "\n",
  sep = ""
)

cat(
  "P = ",
  scaling_df$p_value,
  "\n",
  sep = ""
)


# =============================================================================
# 13. PREPARE PANEL A DATA
# =============================================================================

panelA_df <- gene_df %>%
  
  dplyr::mutate(
    
    gene = factor(
      gene,
      levels = rev(
        expected_genes
      )
    ),
    
    Expected_component =
      dplyr::case_when(
        
        expected_direction ==
          "UP" ~
          "UP component",
        
        expected_direction ==
          "DOWN" ~
          "DOWN component",
        
        TRUE ~
          "Other"
      ),
    
    Nominal_significance =
      dplyr::if_else(
        p_value <
          0.05,
        "Nominal P < 0.05",
        "P >= 0.05"
      ),
    
    q_label =
      paste0(
        "q=",
        format(
          p_BH_five_genes,
          digits = 3,
          scientific = FALSE,
          trim = TRUE
        )
      )
  )


# =============================================================================
# 14. PREPARE PANEL B DATA
# =============================================================================

pairwise_label <- function(
    group_1,
    group_2
) {
  
  paste0(
    group_1,
    " vs ",
    group_2
  )
}


panelB_df <- pairwise_df %>%
  
  dplyr::mutate(
    
    Comparison =
      mapply(
        pairwise_label,
        Group_1,
        Group_2,
        USE.NAMES = FALSE
      ),
    
    minus_log10_BH =
      -log10(
        BH_adjusted_p
      ),
    
    significant_BH =
      BH_adjusted_p <
      0.05,
    
    q_label =
      ifelse(
        BH_adjusted_p <
          0.001,
        format(
          BH_adjusted_p,
          scientific = TRUE,
          digits = 2
        ),
        sprintf(
          "%.3f",
          BH_adjusted_p
        )
      )
  )


panelB_df <- panelB_df %>%
  
  dplyr::arrange(
    minus_log10_BH
  ) %>%
  
  dplyr::mutate(
    
    Comparison =
      factor(
        Comparison,
        levels = Comparison
      )
  )


BH_threshold <- -log10(
  0.05
)


# =============================================================================
# 15. PREPARE PANEL C DATA
# =============================================================================

status_labels <- c(
  
  "Hlty" =
    "Healthy",
  
  "Inf1_P" =
    "Uncomplicated infection",
  
  "Seps_P" =
    "Sepsis",
  
  "Shock_P" =
    "Septic shock"
)


panelC_df <- score_df %>%
  
  dplyr::mutate(
    
    Clinical_group =
      factor(
        unname(
          status_labels[
            status
          ]
        ),
        levels = c(
          "Healthy",
          "Uncomplicated infection",
          "Sepsis",
          "Septic shock"
        )
      )
  )


if (
  any(
    is.na(
      panelC_df$Clinical_group
    )
  )
) {
  
  stop(
    "Unexpected baseline status in frozen sample-level score table."
  )
}


# =============================================================================
# 16. COMMON FIGURE THEME
# =============================================================================

theme_publication <- ggplot2::theme_bw(
  base_size = 11
) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 9.5,
        hjust = 0,
        margin = ggplot2::margin(
          b = 6
        )
      ),
    
    axis.title =
      ggplot2::element_text(
        size = 10.5
      ),
    
    axis.text =
      ggplot2::element_text(
        size = 9.5
      ),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    panel.grid.major.y =
      ggplot2::element_blank(),
    
    legend.title =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom",
    
    plot.margin =
      ggplot2::margin(
        8,
        10,
        8,
        8
      )
  )


# =============================================================================
# 17. PANEL A
# =============================================================================

panelA <- ggplot2::ggplot(
  panelA_df,
  ggplot2::aes(
    x = median_difference_case_minus_control,
    y = gene
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = median_difference_case_minus_control,
      yend = gene,
      color = Expected_component
    ),
    linewidth = 0.9,
    alpha = 0.65
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      color = Expected_component,
      shape = Nominal_significance
    ),
    size = 3.7,
    stroke = 1
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = q_label
    ),
    hjust = ifelse(
      panelA_df$median_difference_case_minus_control >=
        0,
      -0.15,
      1.15
    ),
    size = 3.1,
    show.legend = FALSE
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "UP component" = "#B2182B",
      "DOWN component" = "#2166AC"
    )
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "P >= 0.05" = 16,
      "Nominal P < 0.05" = 8
    )
  ) +
  
  ggplot2::labs(
    title = "A  Component-gene directional replication",
    subtitle = paste0(
      "Primary comparison: sepsis/septic shock vs uncomplicated infection; ",
      "5/5 concordant, 0/5 BH-adjusted P < 0.05"
    ),
    x = "Median expression difference (case minus control)",
    y = NULL
  ) +
  
  ggplot2::expand_limits(
    x = c(
      min(
        panelA_df$median_difference_case_minus_control
      ) -
        0.12,
      max(
        panelA_df$median_difference_case_minus_control
      ) +
        0.12
    )
  ) +
  
  theme_publication


# =============================================================================
# 18. PANEL B
# =============================================================================

panelB <- ggplot2::ggplot(
  panelB_df,
  ggplot2::aes(
    x = minus_log10_BH,
    y = Comparison
  )
) +
  
  ggplot2::geom_vline(
    xintercept = BH_threshold,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = minus_log10_BH,
      yend = Comparison,
      color = significant_BH
    ),
    linewidth = 0.9,
    alpha = 0.7
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      color = significant_BH
    ),
    size = 3.5
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "q=",
        q_label
      )
    ),
    hjust = -0.12,
    size = 3,
    show.legend = FALSE
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "FALSE" = "#777777",
      "TRUE" = "#7B3294"
    ),
    labels = c(
      "FALSE" = "BH q >= 0.05",
      "TRUE" = "BH q < 0.05"
    )
  ) +
  
  ggplot2::expand_limits(
    x =
      max(
        panelB_df$minus_log10_BH
      ) *
      1.15
  ) +
  
  ggplot2::labs(
    title = "B  Pairwise baseline-state contrasts",
    subtitle = paste0(
      "Frozen pairwise BH-adjusted P values; dashed line denotes q = 0.05"
    ),
    x = expression(
      -log[10](
        "BH-adjusted P"
      )
    ),
    y = NULL,
    color = "BH significance"
  ) +
  
  theme_publication


# =============================================================================
# 19. PANEL C
# =============================================================================

panelC <- ggplot2::ggplot(
  panelC_df,
  ggplot2::aes(
    x = five_gene_score,
    y = five_gene_score_healthy_reference,
    color = Clinical_group
  )
) +
  
  ggplot2::geom_point(
    size = 2.6,
    alpha = 0.82
  ) +
  
  ggplot2::annotate(
    "label",
    x = Inf,
    y = -Inf,
    hjust = 1.05,
    vjust = -0.55,
    label = paste0(
      "Frozen Spearman rho = ",
      sprintf(
        "%.3f",
        scaling_df$spearman_rho
      ),
      "\nP = ",
      format(
        scaling_df$p_value,
        scientific = TRUE,
        digits = 3
      )
    ),
    size = 3.5,
    label.size = 0.3
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "Healthy" = "#4D4D4D",
      "Uncomplicated infection" = "#2C7FB8",
      "Sepsis" = "#E08214",
      "Septic shock" = "#B2182B"
    )
  ) +
  
  ggplot2::labs(
    title = "C  Score-scaling sensitivity",
    subtitle = paste0(
      "Frozen sample-level scores; correlation annotation is copied from Script 141"
    ),
    x = "Cohort-standardized five-gene score",
    y = "Healthy-reference-scaled five-gene score",
    color = "Baseline group"
  ) +
  
  theme_publication +
  
  ggplot2::theme(
    panel.grid.major.y =
      ggplot2::element_line(
        linewidth = 0.3,
        color = "grey90"
      )
  )


# =============================================================================
# 20. COMBINE FIGURE
# =============================================================================

combined_figure <-
  (
    panelA |
      panelB
  ) /
  panelC +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1.15
    )
  ) +
  
  patchwork::plot_annotation(
    
    title =
      "Supplementary Figure S7. GSE154918 component-level and sensitivity analyses",
    
    theme =
      ggplot2::theme(
        
        plot.title =
          ggplot2::element_text(
            face = "bold",
            size = 14,
            hjust = 0
          ),
        
        plot.margin =
          ggplot2::margin(
            8,
            8,
            8,
            8
          )
      )
  )


# =============================================================================
# 21. OUTPUT FIGURE FILES
# =============================================================================

pdf_file <- file.path(
  figures_dir,
  "164_FigureS7_GSE154918_sensitivity_details.pdf"
)


png_file <- file.path(
  figures_dir,
  "164_FigureS7_GSE154918_sensitivity_details.png"
)


tiff_file <- file.path(
  figures_dir,
  "164_FigureS7_GSE154918_sensitivity_details.tiff"
)


ggplot2::ggsave(
  filename = pdf_file,
  plot = combined_figure,
  width = 14,
  height = 11,
  units = "in",
  device = cairo_pdf
)


ggplot2::ggsave(
  filename = png_file,
  plot = combined_figure,
  width = 14,
  height = 11,
  units = "in",
  dpi = 400
)


ggplot2::ggsave(
  filename = tiff_file,
  plot = combined_figure,
  width = 14,
  height = 11,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# =============================================================================
# 22. SOURCE-DATA WORKBOOK
# =============================================================================

source_data_file <- file.path(
  tables_dir,
  "164_FigureS7_source_data.xlsx"
)


source_wb <- openxlsx::createWorkbook()


source_tables <- list(
  
  "PanelA_gene_direction" =
    panelA_df,
  
  "PanelB_pairwise_BH" =
    panelB_df,
  
  "PanelC_sample_scores" =
    panelC_df,
  
  "Frozen_scaling_result" =
    scaling_df,
  
  "Schema_audit" =
    schema_audit,
  
  "Score_source_audit" =
    score_equivalence_audit,
  
  "Status_audit" =
    status_audit
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
    source_tables
  )
) {
  
  one_table <- source_tables[[sheet_name]]
  
  
  openxlsx::addWorksheet(
    source_wb,
    sheet_name
  )
  
  
  openxlsx::writeData(
    source_wb,
    sheet_name,
    one_table
  )
  
  
  if (
    ncol(
      one_table
    ) >
    0
  ) {
    
    openxlsx::addStyle(
      source_wb,
      sheet_name,
      header_style,
      rows = 1,
      cols = seq_len(
        ncol(
          one_table
        )
      ),
      gridExpand = TRUE
    )
    
    
    openxlsx::setColWidths(
      source_wb,
      sheet_name,
      cols = seq_len(
        ncol(
          one_table
        )
      ),
      widths = "auto"
    )
  }
  
  
  openxlsx::freezePane(
    source_wb,
    sheet_name,
    firstActiveRow = 2
  )
}


openxlsx::saveWorkbook(
  source_wb,
  source_data_file,
  overwrite = TRUE
)


# =============================================================================
# 23. FIGURE LEGEND
# =============================================================================

figure_legend <- paste0(
  
  "Supplementary Figure S7. Component-level replication and sensitivity ",
  "analyses of the frozen five-gene host-response score in GSE154918. ",
  
  "(A) Median component-gene expression differences for the prespecified ",
  "primary external comparison of sepsis/septic shock versus uncomplicated ",
  "infection. Positive values indicate higher expression in sepsis/septic ",
  "shock and negative values indicate lower expression. CD177, HK3, and ",
  "IRAK3 showed positive differences, whereas CARD11 and IKZF2 showed ",
  "negative differences, yielding directional concordance for all five ",
  "genes. CARD11 was nominally significant, but none of the five genes ",
  "remained significant after Benjamini-Hochberg correction. ",
  
  "(B) Pairwise comparisons among the four baseline clinical groups, ",
  "displayed as -log10 of the frozen BH-adjusted P value. The dashed line ",
  "indicates BH-adjusted P=0.05. These pairwise adjusted P values belong ",
  "to the frozen six-comparison baseline pairwise family and are distinct ",
  "from the five-comparison score-analysis multiplicity family reported ",
  "for the prespecified primary and secondary external comparisons. ",
  
  "(C) Sensitivity of sample ordering to score standardization. Frozen ",
  "cohort-standardized scores are plotted against the alternative ",
  "healthy-reference-scaled scores for the 91 baseline samples. The ",
  "correlation annotation (Spearman rho=",
  sprintf(
    "%.3f",
    scaling_df$spearman_rho
  ),
  ", P=",
  format(
    scaling_df$p_value,
    scientific = TRUE,
    digits = 3
  ),
  ") is copied directly from the frozen Script 141 sensitivity analysis ",
  "and was not recalculated for this figure. ",
  
  "All panels visualize frozen external-evaluation results; no feature ",
  "selection, coefficient refitting, cutoff optimization, score-direction ",
  "reversal, or new statistical testing was performed."
)


legend_file <- file.path(
  text_dir,
  "164_FigureS7_legend_EN.txt"
)


writeLines(
  figure_legend,
  legend_file
)


# =============================================================================
# 24. PROPOSED RESULTS 3.8 SUPPLEMENTARY SENTENCE
# =============================================================================

results_text <- paste0(
  
  "Component-level analyses further supported directional replication of ",
  "the frozen signature: CD177, HK3, and IRAK3 were higher and CARD11 and ",
  "IKZF2 were lower in sepsis/septic shock than in uncomplicated infection, ",
  "yielding concordant direction for all five genes, although none remained ",
  "significant after Benjamini-Hochberg correction across the five ",
  "components (Supplementary Fig. S7A and Supplementary Table S9). ",
  "Pairwise baseline-state analyses provided additional context for the ",
  "ordered disease-state gradient (Supplementary Fig. S7B). The relative ",
  "sample ordering was also highly insensitive to the choice of score ",
  "standardization reference, with a frozen Spearman correlation of ",
  sprintf(
    "%.3f",
    scaling_df$spearman_rho
  ),
  " between cohort-standardized and healthy-reference-scaled scores ",
  "(Supplementary Fig. S7C)."
)


results_file <- file.path(
  text_dir,
  "164_proposed_Results_3.8_supplementary_GSE154918_EN.txt"
)


writeLines(
  results_text,
  results_file
)


# =============================================================================
# 25. INTERNAL AUDIT WORKBOOK
# =============================================================================

audit_file <- file.path(
  audit_dir,
  "164_INTERNAL_AUDIT_FigureS7_GSE154918.xlsx"
)


audit_wb <- openxlsx::createWorkbook()


audit_tables <- list(
  
  "Schema_audit" =
    schema_audit,
  
  "Score_source_equivalence" =
    score_equivalence_audit,
  
  "Status_audit" =
    status_audit,
  
  "Frozen_gene_data" =
    gene_df,
  
  "Frozen_pairwise_data" =
    pairwise_df,
  
  "Frozen_scaling_data" =
    scaling_df
)


for (
  sheet_name in names(
    audit_tables
  )
) {
  
  one_table <- audit_tables[[sheet_name]]
  
  
  openxlsx::addWorksheet(
    audit_wb,
    sheet_name
  )
  
  
  openxlsx::writeData(
    audit_wb,
    sheet_name,
    one_table
  )
  
  
  if (
    ncol(
      one_table
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
          one_table
        )
      ),
      gridExpand = TRUE
    )
    
    
    openxlsx::setColWidths(
      audit_wb,
      sheet_name,
      cols = seq_len(
        ncol(
          one_table
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
# 26. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "164_sessionInfo.txt"
  )
)


# =============================================================================
# 27. FINAL REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 164 FINAL completed successfully.\n")
cat("====================================================================\n\n")


cat("COMPONENT-GENE REPLICATION\n")
cat("--------------------------\n")

cat(
  "Expected-direction concordance = ",
  sum(gene_df$direction_concordant),
  "/5\n",
  sep = ""
)

cat(
  "Nominal P < 0.05 = ",
  n_nominal,
  "/5 (",
  nominal_gene,
  ")\n",
  sep = ""
)

cat(
  "BH-adjusted P < 0.05 = ",
  n_BH,
  "/5\n",
  sep = ""
)


cat("\nPAIRWISE BASELINE CONTRASTS\n")
cat("---------------------------\n")

print(
  panelB_df[
    ,
    c(
      "Comparison",
      "BH_adjusted_p",
      "minus_log10_BH",
      "significant_BH"
    )
  ],
  row.names = FALSE
)


cat("\nSCALING SENSITIVITY\n")
cat("-------------------\n")

cat(
  "Baseline samples = ",
  nrow(panelC_df),
  "\n",
  sep = ""
)

cat(
  "Frozen Spearman rho = ",
  scaling_df$spearman_rho,
  "\n",
  sep = ""
)

cat(
  "Frozen P = ",
  scaling_df$p_value,
  "\n",
  sep = ""
)


cat("\nOUTPUT FILES\n")
cat("------------\n")

cat(
  "Figure S7 PDF:\n",
  normalizePath(
    pdf_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Figure S7 PNG:\n",
  normalizePath(
    png_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Figure S7 TIFF:\n",
  normalizePath(
    tiff_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Figure source data:\n",
  normalizePath(
    source_data_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Figure legend:\n",
  normalizePath(
    legend_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n"
)

cat(
  "Proposed Results 3.8 supplementary text:\n",
  normalizePath(
    results_file,
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
  "\n"
)


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")

cat(
  "- Figure S7 is supplementary to frozen Main Figure 4 and does not replace it.\n"
)

cat(
  "- Panel A uses frozen component-gene statistics from Table S9.\n"
)

cat(
  "- Five of five genes reproduce the expected expression direction.\n"
)

cat(
  "- Zero of five individual genes are BH-significant in the primary external comparison.\n"
)

cat(
  "- Panel B uses frozen six-comparison pairwise baseline BH-adjusted P values.\n"
)

cat(
  "- Pairwise baseline BH values must not be confused with BH correction across the five score comparisons.\n"
)

cat(
  "- Panel C uses frozen sample-level scores and the frozen scaling-sensitivity correlation.\n"
)

cat(
  "- Spearman correlation is not recalculated by Script 164.\n"
)

cat(
  "- The primary sepsis/septic-shock versus uncomplicated-infection comparison remains formally negative.\n"
)

cat(
  "- GSE154918 provides directional and ordinal molecular replication, not calibrated diagnostic validation.\n"
)

cat("\nDone.\n")