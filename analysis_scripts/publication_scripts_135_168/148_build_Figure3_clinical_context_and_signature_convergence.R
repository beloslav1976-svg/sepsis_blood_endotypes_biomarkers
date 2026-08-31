################################################################################
# Script 148
# FINAL Main Figure 3
# Clinical context and convergence of blood transcriptomic signatures
#
# Project: Sepsis_DESeq2
#
# FIGURE STRUCTURE
# ----------------
# A. Primary five-gene host-response score vs CRP
# B. SRSq vs CRP
# C. Convergence of transcriptomic signatures with SRSq
# D. Association of transcriptomic signatures with CTS
#
#
# PURPOSE
# -------
# Publication-packaging of FINALIZED outputs from:
#
#   Script 136
#   Script 136b
#   Script 137
#
#
# IMPORTANT
# ---------
# This script does NOT:
#   - rerun clinical association testing;
#   - rerun multiplicity correction;
#   - recalculate SRS/SRSq;
#   - recalculate CTS;
#   - recalculate biomarker scores;
#   - reconstruct published signatures;
#   - rerun Kruskal-Wallis tests;
#   - rerun Spearman tests;
#   - perform feature selection;
#   - perform diagnostic benchmarking;
#   - perform new hypothesis testing.
#
# Scatter-plot regression lines in Panels A and B are descriptive
# visual guides only. All reported rho, P, and BH-adjusted P values
# are read directly from finalized upstream results.
#
#
# PRIMARY SOURCE FILES
# --------------------
# Script 136b:
#
# results/blood_endotypes_biomarkers/
#   136b_demographic_sensitivity/tables/
#     136b_BP_molecular_clinical_demographics.csv
#
# results/blood_endotypes_biomarkers/
#   136b_demographic_sensitivity/tables/
#     136b_all_clinical_tests_updated_FDR.csv
#
#
# Script 137:
#
# results/blood_endotypes_biomarkers/
#   137_benchmarking/tables/
#     137_signature_correlations_primary_SRSq.csv
#
# results/blood_endotypes_biomarkers/
#   137_benchmarking/tables/
#     137_signature_endotype_associations.csv
#
#
# EXPECTED KEY RESULTS
# --------------------
# Global evaluable clinical family:
#   60 tests
#
# Global BH significant:
#   exactly 2 tests
#
# Primary 5-gene score vs CRP:
#   n = 35
#   rho = 0.5743503664
#   P = 0.0003085409
#   global BH = 0.0185124559
#
# SRSq vs CRP:
#   n = 35
#   rho = 0.5260208843
#   P = 0.0011724381
#   global BH = 0.0351731433
#
#
# SRSq correlations:
#
#   LIFTS-like                 0.851821
#   DCAF17 alternative         0.794398
#   Primary 5-gene            0.764986
#   FAIM3:PLAC8-like          0.716246
#   Sepsis MetaScore-like     0.684874
#   PLAC8-PLA2G7 contrast     0.494678
#   SeptiCyte LAB-like        0.462465
#
#
# CTS epsilon-squared:
#
#   LIFTS-like                 0.697868
#   Primary 5-gene            0.660677
#   FAIM3:PLAC8-like          0.649247
#   DCAF17 alternative         0.609646
#   Sepsis MetaScore-like     0.529236
#   PLAC8-PLA2G7 contrast     0.436113
#   SeptiCyte LAB-like        0.328125
#
#
# INTERPRETATION
# --------------
# These analyses assess biological convergence on a common host-response axis.
#
# They do NOT constitute:
#   - superiority testing;
#   - head-to-head diagnostic validation;
#   - sepsis-versus-SIRS validation;
#   - external validation.
#
#
# OUTPUT
# ------
# results/blood_endotypes_biomarkers/
#   148_Figure3_clinical_signature_convergence/
#
################################################################################


cat("====================================================================\n")
cat("Running Script 148\n")
cat("FINAL Main Figure 3\n")
cat("Clinical context and transcriptomic-signature convergence\n")
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
  "tidyr",
  "ggplot2",
  "patchwork",
  "openxlsx",
  "scales",
  "stringr",
  "forcats"
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
  library(scales)
  library(stringr)
  library(forcats)
  
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


read_source_table <- function(
    file,
    xlsx_sheet = NULL
) {
  
  ext <- tolower(
    tools::file_ext(
      file
    )
  )
  
  
  if (ext == "csv") {
    
    return(
      read.csv(
        file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }
  
  
  if (ext == "xlsx") {
    
    if (is.null(xlsx_sheet)) {
      
      stop(
        paste0(
          "xlsx_sheet must be supplied for:\n",
          file
        )
      )
    }
    
    
    return(
      openxlsx::read.xlsx(
        file,
        sheet = xlsx_sheet
      )
    )
  }
  
  
  stop(
    paste0(
      "Unsupported source format: ",
      ext
    )
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
  
  x <- as.numeric(x)
  
  
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
  
  
  return(
    sprintf(
      "%.4f",
      x
    )
  )
}


format_effect <- function(
    x,
    digits = 3
) {
  
  sprintf(
    paste0(
      "%.",
      digits,
      "f"
    ),
    as.numeric(x)
  )
}


# =============================================================================
# 4. INPUT FILES
# =============================================================================


# -----------------------------------------------------------------------------
# Script 136b — BP molecular/clinical table
# -----------------------------------------------------------------------------

bp_clinical_file <- find_project_file(
  
  candidates = c(
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "136b_demographic_sensitivity",
      "tables",
      "136b_BP_molecular_clinical_demographics.csv"
    ),
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "136b_demographic_sensitivity",
      "tables",
      "136b_demographic_sensitivity.xlsx"
    )
  ),
  
  recursive_pattern =
    "136b_BP_molecular_clinical_demographics\\.csv$",
  
  description =
    "Script 136b BP molecular-clinical table"
)


if (
  is.na(bp_clinical_file) ||
  !file.exists(bp_clinical_file)
) {
  
  stop(
    "Could not locate Script 136b BP clinical table."
  )
}


# -----------------------------------------------------------------------------
# Script 136b — final clinical tests after demographic integration
# -----------------------------------------------------------------------------

clinical_tests_file <- find_project_file(
  
  candidates = c(
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "136b_demographic_sensitivity",
      "tables",
      "136b_all_clinical_tests_updated_FDR.csv"
    ),
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "136b_demographic_sensitivity",
      "tables",
      "136b_demographic_sensitivity.xlsx"
    )
  ),
  
  recursive_pattern =
    "136b_all_clinical_tests_updated_FDR\\.csv$",
  
  description =
    "Script 136b final clinical association tests"
)


if (
  is.na(clinical_tests_file) ||
  !file.exists(clinical_tests_file)
) {
  
  stop(
    "Could not locate Script 136b final clinical test table."
  )
}


# -----------------------------------------------------------------------------
# Script 137 — correlations with primary score / SRSq
# -----------------------------------------------------------------------------

signature_corr_file <- find_project_file(
  
  candidates = c(
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "137_benchmarking",
      "tables",
      "137_signature_correlations_primary_SRSq.csv"
    ),
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "137_benchmarking",
      "tables",
      "137_blood_biomarker_benchmarking.xlsx"
    )
  ),
  
  recursive_pattern =
    "137_signature_correlations_primary_SRSq\\.csv$",
  
  description =
    "Script 137 signature-SRSq correlations"
)


if (
  is.na(signature_corr_file) ||
  !file.exists(signature_corr_file)
) {
  
  stop(
    "Could not locate Script 137 signature correlation table."
  )
}


# -----------------------------------------------------------------------------
# Script 137 — endotype associations
# -----------------------------------------------------------------------------

endotype_file <- find_project_file(
  
  candidates = c(
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "137_benchmarking",
      "tables",
      "137_signature_endotype_associations.csv"
    ),
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "137_benchmarking",
      "tables",
      "137_blood_biomarker_benchmarking.xlsx"
    )
  ),
  
  recursive_pattern =
    "137_signature_endotype_associations\\.csv$",
  
  description =
    "Script 137 signature-endotype associations"
)


if (
  is.na(endotype_file) ||
  !file.exists(endotype_file)
) {
  
  stop(
    "Could not locate Script 137 endotype-association table."
  )
}


cat("\n====================================================================\n")
cat("FIGURE 3 SOURCE FILES\n")
cat("====================================================================\n")


cat(
  "136b BP clinical table:\n  ",
  normalizePath(
    bp_clinical_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "136b clinical tests:\n  ",
  normalizePath(
    clinical_tests_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "137 correlations:\n  ",
  normalizePath(
    signature_corr_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "137 endotype associations:\n  ",
  normalizePath(
    endotype_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 5. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "148_Figure3_clinical_signature_convergence"
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
# 6. LOAD SOURCE TABLES
# =============================================================================


# -----------------------------------------------------------------------------
# BP molecular / clinical table
# -----------------------------------------------------------------------------

if (
  tolower(
    tools::file_ext(
      bp_clinical_file
    )
  ) == "xlsx"
) {
  
  bp_clinical <- read_source_table(
    bp_clinical_file,
    xlsx_sheet = "02_BP_demographics"
  )
  
} else {
  
  bp_clinical <- read_source_table(
    bp_clinical_file
  )
}


# -----------------------------------------------------------------------------
# Clinical test table
# -----------------------------------------------------------------------------

if (
  tolower(
    tools::file_ext(
      clinical_tests_file
    )
  ) == "xlsx"
) {
  
  clinical_tests <- read_source_table(
    clinical_tests_file,
    xlsx_sheet = "04_updated_all_tests"
  )
  
} else {
  
  clinical_tests <- read_source_table(
    clinical_tests_file
  )
}


# -----------------------------------------------------------------------------
# Signature correlation table
# -----------------------------------------------------------------------------

if (
  tolower(
    tools::file_ext(
      signature_corr_file
    )
  ) == "xlsx"
) {
  
  signature_corr <- read_source_table(
    signature_corr_file,
    xlsx_sheet = "07_corr_primary_SRSq"
  )
  
} else {
  
  signature_corr <- read_source_table(
    signature_corr_file
  )
}


# -----------------------------------------------------------------------------
# Signature endotype table
# -----------------------------------------------------------------------------

if (
  tolower(
    tools::file_ext(
      endotype_file
    )
  ) == "xlsx"
) {
  
  endotype_tests <- read_source_table(
    endotype_file,
    xlsx_sheet = "08_endotype_tests"
  )
  
} else {
  
  endotype_tests <- read_source_table(
    endotype_file
  )
}


# =============================================================================
# 7. REQUIRED COLUMN AUDIT
# =============================================================================

required_bp_columns <- c(
  "sample_id",
  "primary_5gene_score",
  "SRSq",
  "crp_numeric"
)


missing_bp_columns <- setdiff(
  required_bp_columns,
  names(
    bp_clinical
  )
)


if (length(missing_bp_columns) > 0) {
  
  stop(
    paste0(
      "Missing BP clinical columns:\n",
      paste(
        missing_bp_columns,
        collapse = ", "
      )
    )
  )
}


required_test_columns <- c(
  "framework",
  "clinical_variable",
  "test",
  "n",
  "effect",
  "p_value",
  "BH_global"
)


missing_test_columns <- setdiff(
  required_test_columns,
  names(
    clinical_tests
  )
)


if (length(missing_test_columns) > 0) {
  
  stop(
    paste0(
      "Missing clinical-test columns:\n",
      paste(
        missing_test_columns,
        collapse = ", "
      )
    )
  )
}


required_corr_columns <- c(
  "score_column",
  "target",
  "n",
  "Spearman_rho",
  "p_value",
  "BH_correlation"
)


missing_corr_columns <- setdiff(
  required_corr_columns,
  names(
    signature_corr
  )
)


if (length(missing_corr_columns) > 0) {
  
  stop(
    paste0(
      "Missing signature-correlation columns:\n",
      paste(
        missing_corr_columns,
        collapse = ", "
      )
    )
  )
}


required_endotype_columns <- c(
  "score_column",
  "framework",
  "n",
  "effect",
  "effect_name",
  "p_value",
  "BH_endotype"
)


missing_endotype_columns <- setdiff(
  required_endotype_columns,
  names(
    endotype_tests
  )
)


if (length(missing_endotype_columns) > 0) {
  
  stop(
    paste0(
      "Missing endotype-test columns:\n",
      paste(
        missing_endotype_columns,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 8. NORMALIZE NUMERIC COLUMNS
# =============================================================================

bp_clinical <- bp_clinical %>%
  
  dplyr::mutate(
    
    primary_5gene_score =
      as.numeric(
        primary_5gene_score
      ),
    
    SRSq =
      as.numeric(
        SRSq
      ),
    
    crp_numeric =
      as.numeric(
        crp_numeric
      )
  )


clinical_tests <- clinical_tests %>%
  
  dplyr::mutate(
    
    n =
      as.numeric(
        n
      ),
    
    effect =
      as.numeric(
        effect
      ),
    
    p_value =
      as.numeric(
        p_value
      ),
    
    BH_global =
      as.numeric(
        BH_global
      )
  )


signature_corr <- signature_corr %>%
  
  dplyr::mutate(
    
    n =
      as.numeric(
        n
      ),
    
    Spearman_rho =
      as.numeric(
        Spearman_rho
      ),
    
    p_value =
      as.numeric(
        p_value
      ),
    
    BH_correlation =
      as.numeric(
        BH_correlation
      )
  )


endotype_tests <- endotype_tests %>%
  
  dplyr::mutate(
    
    n =
      as.numeric(
        n
      ),
    
    effect =
      as.numeric(
        effect
      ),
    
    p_value =
      as.numeric(
        p_value
      ),
    
    BH_endotype =
      as.numeric(
        BH_endotype
      )
  )


# =============================================================================
# 9. EXTRACT FINAL CRP RESULTS
# =============================================================================


score_crp_test <- clinical_tests %>%
  
  dplyr::filter(
    framework ==
      "Primary_5gene_score",
    clinical_variable ==
      "crp_numeric",
    test ==
      "Spearman"
  )


srsq_crp_test <- clinical_tests %>%
  
  dplyr::filter(
    framework ==
      "SRSq",
    clinical_variable ==
      "crp_numeric",
    test ==
      "Spearman"
  )


if (nrow(score_crp_test) != 1) {
  
  stop(
    paste0(
      "Expected exactly one Primary_5gene_score vs CRP row; observed ",
      nrow(
        score_crp_test
      )
    )
  )
}


if (nrow(srsq_crp_test) != 1) {
  
  stop(
    paste0(
      "Expected exactly one SRSq vs CRP row; observed ",
      nrow(
        srsq_crp_test
      )
    )
  )
}


# =============================================================================
# 10. EXTRACT PANEL C — SRSq CONVERGENCE
# =============================================================================

panel_C <- signature_corr %>%
  
  dplyr::filter(
    target ==
      "SRSq"
  )


if (nrow(panel_C) != 7) {
  
  stop(
    paste0(
      "Expected 7 SRSq signature correlations; observed ",
      nrow(
        panel_C
      )
    )
  )
}


# =============================================================================
# 11. EXTRACT PANEL D — CTS ASSOCIATIONS
# =============================================================================

panel_D <- endotype_tests %>%
  
  dplyr::filter(
    framework ==
      "CTS",
    effect_name ==
      "epsilon_squared"
  )


if (nrow(panel_D) != 7) {
  
  stop(
    paste0(
      "Expected 7 CTS epsilon-squared results; observed ",
      nrow(
        panel_D
      )
    )
  )
}


# =============================================================================
# 12. SIGNATURE REGISTRY FOR FIGURE DISPLAY
# =============================================================================

signature_registry <- data.frame(
  
  score_column = c(
    "LIFTS_like",
    "DCAF17_alt_score",
    "primary_5gene_score",
    "FAIM3_PLAC8_like",
    "Sepsis_MetaScore_like",
    "SeptiCyte_RAPID_like",
    "SeptiCyte_LAB_like"
  ),
  
  figure_label = c(
    "LIFTS-like",
    "DCAF17 alternative",
    "Primary 5-gene",
    "FAIM3:PLAC8-like",
    "Sepsis MetaScore-like",
    "PLAC8\u2013PLA2G7 contrast",
    "SeptiCyte LAB-like"
  ),
  
  origin = c(
    "Published comparator",
    "Current-study sensitivity",
    "Current-study primary",
    "Published comparator",
    "Published comparator",
    "Published comparator",
    "Published comparator"
  ),
  
  stringsAsFactors = FALSE
)


signature_order <- signature_registry$figure_label


panel_C <- panel_C %>%
  
  dplyr::left_join(
    signature_registry,
    by = "score_column"
  )


panel_D <- panel_D %>%
  
  dplyr::left_join(
    signature_registry,
    by = "score_column"
  )


if (
  any(
    is.na(
      panel_C$figure_label
    )
  )
) {
  
  stop(
    "Unmapped signature detected in Panel C."
  )
}


if (
  any(
    is.na(
      panel_D$figure_label
    )
  )
) {
  
  stop(
    "Unmapped signature detected in Panel D."
  )
}


panel_C$figure_label <- factor(
  panel_C$figure_label,
  levels =
    rev(
      signature_order
    )
)


panel_D$figure_label <- factor(
  panel_D$figure_label,
  levels =
    rev(
      signature_order
    )
)


# =============================================================================
# 13. NUMERICAL AUDIT
# =============================================================================


expected_score_crp_rho <- 0.574350366395657
expected_score_crp_p <- 0.00030854093164314
expected_score_crp_q <- 0.0185124558985884

expected_srsq_crp_rho <- 0.526020884345291
expected_srsq_crp_p <- 0.00117243811059945
expected_srsq_crp_q <- 0.0351731433179835


expected_rho <- c(
  
  LIFTS_like =
    0.851821,
  
  DCAF17_alt_score =
    0.794398,
  
  primary_5gene_score =
    0.764986,
  
  FAIM3_PLAC8_like =
    0.716246,
  
  Sepsis_MetaScore_like =
    0.684874,
  
  SeptiCyte_RAPID_like =
    0.494678,
  
  SeptiCyte_LAB_like =
    0.462465
)


expected_cts_epsilon <- c(
  
  LIFTS_like =
    0.697868,
  
  primary_5gene_score =
    0.660677,
  
  FAIM3_PLAC8_like =
    0.649247,
  
  DCAF17_alt_score =
    0.609646,
  
  Sepsis_MetaScore_like =
    0.529236,
  
  SeptiCyte_RAPID_like =
    0.436113,
  
  SeptiCyte_LAB_like =
    0.328125
)


global_sig <- clinical_tests %>%
  
  dplyr::filter(
    is.finite(
      BH_global
    ),
    BH_global < 0.05
  )


audit_table <- data.frame(
  
  metric = c(
    "BP samples",
    "Clinical tests after demographic integration",
    "Global BH significant clinical tests",
    "Primary score vs CRP n",
    "Primary score vs CRP rho",
    "Primary score vs CRP P",
    "Primary score vs CRP global BH",
    "SRSq vs CRP n",
    "SRSq vs CRP rho",
    "SRSq vs CRP P",
    "SRSq vs CRP global BH",
    "Panel C signatures",
    "Panel D signatures"
  ),
  
  observed = c(
    nrow(
      bp_clinical
    ),
    nrow(
      clinical_tests
    ),
    nrow(
      global_sig
    ),
    score_crp_test$n[1],
    score_crp_test$effect[1],
    score_crp_test$p_value[1],
    score_crp_test$BH_global[1],
    srsq_crp_test$n[1],
    srsq_crp_test$effect[1],
    srsq_crp_test$p_value[1],
    srsq_crp_test$BH_global[1],
    nrow(
      panel_C
    ),
    nrow(
      panel_D
    )
  ),
  
  expected = c(
    35,
    60,
    2,
    35,
    expected_score_crp_rho,
    expected_score_crp_p,
    expected_score_crp_q,
    35,
    expected_srsq_crp_rho,
    expected_srsq_crp_p,
    expected_srsq_crp_q,
    7,
    7
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


# -----------------------------------------------------------------------------
# Hard central checks
# -----------------------------------------------------------------------------

if (
  nrow(
    bp_clinical
  ) != 35
) {
  
  stop(
    "Expected 35 BP samples."
  )
}


if (
  nrow(
    clinical_tests
  ) != 60
) {
  
  stop(
    paste0(
      "Expected 60 final clinical tests after demographic integration; observed ",
      nrow(
        clinical_tests
      )
    )
  )
}


if (
  nrow(
    global_sig
  ) != 2
) {
  
  stop(
    paste0(
      "Expected exactly 2 global-BH-significant clinical tests; observed ",
      nrow(
        global_sig
      )
    )
  )
}


expected_global_frameworks <- sort(
  c(
    "Primary_5gene_score",
    "SRSq"
  )
)


observed_global_frameworks <- sort(
  global_sig$framework
)


if (
  !identical(
    expected_global_frameworks,
    observed_global_frameworks
  )
) {
  
  stop(
    paste0(
      "Unexpected globally significant framework(s): ",
      paste(
        observed_global_frameworks,
        collapse = ", "
      )
    )
  )
}


if (
  any(
    global_sig$clinical_variable !=
    "crp_numeric"
  )
) {
  
  stop(
    "A globally significant final clinical test was not CRP."
  )
}


# -----------------------------------------------------------------------------
# Exact CRP result checks
# -----------------------------------------------------------------------------

if (
  abs(
    score_crp_test$effect[1] -
    expected_score_crp_rho
  ) > 1e-8
) {
  
  stop(
    "Primary five-gene score vs CRP rho mismatch."
  )
}


if (
  abs(
    score_crp_test$p_value[1] -
    expected_score_crp_p
  ) > 1e-10
) {
  
  stop(
    "Primary five-gene score vs CRP P mismatch."
  )
}


if (
  abs(
    score_crp_test$BH_global[1] -
    expected_score_crp_q
  ) > 1e-10
) {
  
  stop(
    "Primary five-gene score vs CRP global BH mismatch."
  )
}


if (
  abs(
    srsq_crp_test$effect[1] -
    expected_srsq_crp_rho
  ) > 1e-8
) {
  
  stop(
    "SRSq vs CRP rho mismatch."
  )
}


if (
  abs(
    srsq_crp_test$p_value[1] -
    expected_srsq_crp_p
  ) > 1e-10
) {
  
  stop(
    "SRSq vs CRP P mismatch."
  )
}


if (
  abs(
    srsq_crp_test$BH_global[1] -
    expected_srsq_crp_q
  ) > 1e-10
) {
  
  stop(
    "SRSq vs CRP global BH mismatch."
  )
}


# -----------------------------------------------------------------------------
# Panel C expected-value audit
# -----------------------------------------------------------------------------

for (score_name in names(
  expected_rho
)) {
  
  observed_row <- panel_C %>%
    
    dplyr::filter(
      score_column ==
        score_name
    )
  
  
  if (
    nrow(
      observed_row
    ) != 1
  ) {
    
    stop(
      paste0(
        "Missing or duplicated Panel C score: ",
        score_name
      )
    )
  }
  
  
  if (
    abs(
      observed_row$Spearman_rho[1] -
      expected_rho[score_name]
    ) > 1e-5
  ) {
    
    stop(
      paste0(
        "Panel C rho mismatch for ",
        score_name,
        ". Observed=",
        observed_row$Spearman_rho[1],
        "; expected~",
        expected_rho[score_name]
      )
    )
  }
}


# -----------------------------------------------------------------------------
# Panel D expected-value audit
# -----------------------------------------------------------------------------

for (score_name in names(
  expected_cts_epsilon
)) {
  
  observed_row <- panel_D %>%
    
    dplyr::filter(
      score_column ==
        score_name
    )
  
  
  if (
    nrow(
      observed_row
    ) != 1
  ) {
    
    stop(
      paste0(
        "Missing or duplicated Panel D score: ",
        score_name
      )
    )
  }
  
  
  if (
    abs(
      observed_row$effect[1] -
      expected_cts_epsilon[score_name]
    ) > 1e-5
  ) {
    
    stop(
      paste0(
        "Panel D epsilon-squared mismatch for ",
        score_name,
        ". Observed=",
        observed_row$effect[1],
        "; expected~",
        expected_cts_epsilon[score_name]
      )
    )
  }
}


cat("\nNumerical audit passed successfully.\n")


# =============================================================================
# 14. PUBLICATION COLORS
# =============================================================================

col_primary <- "#D55E00"
col_sensitivity <- "#E69F00"
col_published <- "#0072B2"

col_score_scatter <- "#D55E00"
col_srsq_scatter <- "#0072B2"

col_trend <- "#4D4D4D"
col_segment <- "#B8B8B8"


origin_colors <- c(
  
  "Current-study primary" =
    col_primary,
  
  "Current-study sensitivity" =
    col_sensitivity,
  
  "Published comparator" =
    col_published
)


# =============================================================================
# 15. PANEL A — FIVE-GENE SCORE vs CRP
# =============================================================================

panel_A_data <- bp_clinical %>%
  
  dplyr::filter(
    is.finite(
      primary_5gene_score
    ),
    is.finite(
      crp_numeric
    )
  )


if (
  nrow(
    panel_A_data
  ) !=
  score_crp_test$n[1]
) {
  
  stop(
    "Panel A complete-case n does not match finalized clinical test."
  )
}


panel_A_annotation <- paste0(
  
  "Spearman \u03c1 = ",
  sprintf(
    "%.3f",
    score_crp_test$effect[1]
  ),
  
  "\nP = ",
  format_p(
    score_crp_test$p_value[1]
  ),
  
  "\nGlobal BH q = ",
  sprintf(
    "%.4f",
    score_crp_test$BH_global[1]
  )
)


p_A <- ggplot2::ggplot(
  
  panel_A_data,
  
  ggplot2::aes(
    x = crp_numeric,
    y = primary_5gene_score
  )
  
) +
  
  ggplot2::geom_point(
    size = 2.8,
    alpha = 0.82,
    color = col_score_scatter
  ) +
  
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.75,
    linetype = "solid",
    color = col_trend
  ) +
  
  ggplot2::annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = panel_A_annotation,
    hjust = -0.05,
    vjust = 1.08,
    size = 3.3,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    
    tag =
      "A",
    
    title =
      "Five-gene host-response score and CRP",
    
    subtitle =
      "Discovery sepsis cohort; n=35",
    
    x =
      "C-reactive protein (mg/L)",
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 16. PANEL B — SRSq vs CRP
# =============================================================================

panel_B_data <- bp_clinical %>%
  
  dplyr::filter(
    is.finite(
      SRSq
    ),
    is.finite(
      crp_numeric
    )
  )


if (
  nrow(
    panel_B_data
  ) !=
  srsq_crp_test$n[1]
) {
  
  stop(
    "Panel B complete-case n does not match finalized clinical test."
  )
}


panel_B_annotation <- paste0(
  
  "Spearman \u03c1 = ",
  sprintf(
    "%.3f",
    srsq_crp_test$effect[1]
  ),
  
  "\nP = ",
  format_p(
    srsq_crp_test$p_value[1]
  ),
  
  "\nGlobal BH q = ",
  sprintf(
    "%.4f",
    srsq_crp_test$BH_global[1]
  )
)


p_B <- ggplot2::ggplot(
  
  panel_B_data,
  
  ggplot2::aes(
    x = crp_numeric,
    y = SRSq
  )
  
) +
  
  ggplot2::geom_point(
    size = 2.8,
    alpha = 0.82,
    color = col_srsq_scatter
  ) +
  
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.75,
    linetype = "solid",
    color = col_trend
  ) +
  
  ggplot2::annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = panel_B_annotation,
    hjust = -0.05,
    vjust = 1.08,
    size = 3.3,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    
    tag =
      "B",
    
    title =
      "Continuous SRSq host-response state and CRP",
    
    subtitle =
      "Discovery sepsis cohort; n=35",
    
    x =
      "C-reactive protein (mg/L)",
    
    y =
      "SRSq"
  )


# =============================================================================
# 17. PANEL C — CONVERGENCE WITH SRSq
# =============================================================================

panel_C <- panel_C %>%
  
  dplyr::mutate(
    
    rho_label =
      sprintf(
        "%.3f",
        Spearman_rho
      )
  )


p_C <- ggplot2::ggplot(
  
  panel_C,
  
  ggplot2::aes(
    x = Spearman_rho,
    y = figure_label
  )
  
) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend = Spearman_rho,
      yend = figure_label
    ),
    
    linewidth = 0.8,
    color = col_segment
  ) +
  
  ggplot2::geom_point(
    
    ggplot2::aes(
      color = origin
    ),
    
    size = 4.0,
    alpha = 0.95
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = rho_label
    ),
    
    hjust = -0.35,
    size = 3.1,
    color = "black"
  ) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    color = "grey45"
  ) +
  
  ggplot2::scale_color_manual(
    values = origin_colors
  ) +
  
  ggplot2::scale_x_continuous(
    
    limits = c(
      0,
      0.95
    ),
    
    breaks = seq(
      0,
      0.8,
      by = 0.2
    ),
    
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    axis.title.y =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom"
  ) +
  
  ggplot2::labs(
    
    tag =
      "C",
    
    title =
      "Convergence with continuous SRSq",
    
    subtitle =
      "Spearman correlations within sepsis blood; n=35",
    
    x =
      "Spearman \u03c1",
    
    color =
      NULL
  )


# =============================================================================
# 18. PANEL D — ASSOCIATION WITH CTS
# =============================================================================

panel_D <- panel_D %>%
  
  dplyr::mutate(
    
    epsilon_label =
      sprintf(
        "%.3f",
        effect
      )
  )


p_D <- ggplot2::ggplot(
  
  panel_D,
  
  ggplot2::aes(
    x = effect,
    y = figure_label
  )
  
) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend = effect,
      yend = figure_label
    ),
    
    linewidth = 0.8,
    color = col_segment
  ) +
  
  ggplot2::geom_point(
    
    ggplot2::aes(
      color = origin
    ),
    
    size = 4.0,
    alpha = 0.95
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = epsilon_label
    ),
    
    hjust = -0.35,
    size = 3.1,
    color = "black"
  ) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    color = "grey45"
  ) +
  
  ggplot2::scale_color_manual(
    values = origin_colors
  ) +
  
  ggplot2::scale_x_continuous(
    
    limits = c(
      0,
      0.80
    ),
    
    breaks = seq(
      0,
      0.8,
      by = 0.2
    ),
    
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    axis.title.y =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom"
  ) +
  
  ggplot2::labs(
    
    tag =
      "D",
    
    title =
      "Association with Consensus Transcriptomic Subtypes",
    
    subtitle =
      "Kruskal\u2013Wallis effect size within sepsis blood; n=35",
    
    x =
      expression(
        epsilon^2
      ),
    
    color =
      NULL
  )


# =============================================================================
# 19. ASSEMBLE FINAL FIGURE 3
# =============================================================================
#
# No global title inside the figure.
#
# Scientific figure title is retained in the figure legend.
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
      1,
      1
    )
  )


figure3 <- (
  top_row /
    bottom_row
) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1.04
    ),
    guides = "collect"
  ) &
  
  ggplot2::theme(
    legend.position = "bottom"
  )


# =============================================================================
# 20. EXPORT MAIN FIGURE
# =============================================================================

figure_png <- file.path(
  figures_dir,
  "148_Figure3_clinical_context_signature_convergence.png"
)


figure_pdf <- file.path(
  figures_dir,
  "148_Figure3_clinical_context_signature_convergence.pdf"
)


figure_tiff <- file.path(
  figures_dir,
  "148_Figure3_clinical_context_signature_convergence.tiff"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure3,
  width = 13.5,
  height = 9.3,
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure3,
  width = 13.5,
  height = 9.3,
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
  plot = figure3,
  width = 13.5,
  height = 9.3,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 21. EXPORT INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  
  A_score_vs_CRP =
    p_A,
  
  B_SRSq_vs_CRP =
    p_B,
  
  C_signature_vs_SRSq =
    p_C,
  
  D_signature_vs_CTS =
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
          "148_panel_",
          nm,
          ".png"
        )
      ),
    
    plot =
      individual_panels[[nm]],
    
    width =
      6.7,
    
    height =
      5.3,
    
    dpi =
      600,
    
    bg =
      "white"
  )
}


# =============================================================================
# 22. FIGURE SOURCE WORKBOOK
# =============================================================================

run_info <- data.frame(
  
  item = c(
    "script",
    "status",
    "analysis_mode",
    "new_hypothesis_tests",
    "clinical_multiplicity_source",
    "clinical_test_family_size",
    "benchmark_source",
    "benchmark_interpretation",
    "Panel_A",
    "Panel_B",
    "Panel_C",
    "Panel_D"
  ),
  
  value = c(
    "148_build_Figure3_clinical_context_and_signature_convergence.R",
    "publication packaging",
    "frozen upstream results only",
    "NO",
    "Script 136b",
    "60 evaluable tests",
    "Script 137",
    "biological convergence; not superiority or diagnostic validation",
    "Primary five-gene score vs CRP",
    "SRSq vs CRP",
    "Signature correlations with SRSq",
    "Signature associations with CTS"
  ),
  
  stringsAsFactors = FALSE
)


crp_test_source <- dplyr::bind_rows(
  score_crp_test,
  srsq_crp_test
)


source_workbook <- file.path(
  tables_dir,
  "148_Figure3_source_data.xlsx"
)


openxlsx::write.xlsx(
  
  list(
    
    Run_info =
      run_info,
    
    Numerical_audit =
      audit_table,
    
    Global_BH_significant =
      global_sig,
    
    Clinical_BP_data =
      bp_clinical,
    
    CRP_tests =
      crp_test_source,
    
    Panel_A_score_CRP =
      panel_A_data,
    
    Panel_B_SRSq_CRP =
      panel_B_data,
    
    Panel_C_SRSq_convergence =
      panel_C,
    
    Panel_D_CTS_convergence =
      panel_D
    
  ),
  
  source_workbook,
  
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
    "Primary five-gene score vs CRP",
    "SRSq vs CRP",
    "Transcriptomic-signature correlations with SRSq",
    "Transcriptomic-signature association with CTS"
  ),
  
  Upstream_script = c(
    "136 + 136b",
    "136 + 136b",
    "137",
    "137"
  ),
  
  Source_file = c(
    basename(
      clinical_tests_file
    ),
    basename(
      clinical_tests_file
    ),
    basename(
      signature_corr_file
    ),
    basename(
      endotype_file
    )
  ),
  
  New_statistical_analysis = rep(
    "NO",
    4
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  provenance,
  file.path(
    tables_dir,
    "148_Figure3_provenance.csv"
  ),
  row.names = FALSE
)


write.csv(
  audit_table,
  file.path(
    tables_dir,
    "148_Figure3_numerical_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 24. FINAL FIGURE LEGEND — ENGLISH
# =============================================================================

caption_en <- paste0(
  
  "Figure 3. The five-gene host-response axis is associated with systemic ",
  "inflammation and converges with previously published sepsis transcriptomic ",
  "signatures. ",
  
  "(A) Association between the primary five-gene host-response score and ",
  "circulating C-reactive protein (CRP) concentration in the discovery sepsis ",
  "cohort (n=35). The score was positively associated with CRP (Spearman ",
  "rho=",
  sprintf(
    "%.3f",
    score_crp_test$effect[1]
  ),
  ", P=",
  format(
    score_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; global Benjamini-Hochberg-adjusted P=",
  sprintf(
    "%.4f",
    score_crp_test$BH_global[1]
  ),
  "). ",
  
  "(B) Association between the continuous SRSq host-response measure and CRP ",
  "(Spearman rho=",
  sprintf(
    "%.3f",
    srsq_crp_test$effect[1]
  ),
  ", P=",
  format(
    srsq_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; global Benjamini-Hochberg-adjusted P=",
  sprintf(
    "%.4f",
    srsq_crp_test$BH_global[1]
  ),
  "). These were the only clinical associations that remained significant ",
  "after global multiplicity correction across the 60 evaluable exploratory ",
  "clinical-association tests after demographic integration. Lines in panels ",
  "A and B are descriptive linear trends for visualization; inferential ",
  "statistics are the prespecified Spearman analyses shown in the panels. ",
  
  "(C) Spearman correlations between SRSq and the primary five-gene score, ",
  "the DCAF17-containing alternative score, and platform-adapted implementations ",
  "of previously published sepsis transcriptomic signatures. ",
  
  "(D) Association of the same signatures with Consensus Transcriptomic ",
  "Subtypes (CTS), summarized using epsilon-squared effect sizes from the ",
  "corresponding Kruskal-Wallis analyses. All benchmarking analyses were ",
  "performed within the same 35-sample sepsis blood cohort. SeptiCyte-derived ",
  "RNA-seq contrasts are research implementations and are not official ",
  "SeptiScores. These comparisons assess convergence on a shared molecular ",
  "host-response axis and should not be interpreted as head-to-head diagnostic ",
  "performance or evidence of superiority of any individual signature."
)


writeLines(
  caption_en,
  file.path(
    text_dir,
    "148_Figure3_caption_EN.txt"
  )
)


# =============================================================================
# 25. FINAL FIGURE LEGEND — RUSSIAN
# =============================================================================

caption_ru <- paste0(
  
  "Рисунок 3. Пятигенная ось host response связана с системной воспалительной ",
  "активностью и сходится с ранее опубликованными транскриптомными сигнатурами ",
  "сепсиса. ",
  
  "(A) Связь первичного пятигенного host-response score с концентрацией ",
  "C-реактивного белка (CRP) в discovery-когорте пациентов с сепсисом (n=35). ",
  "Score положительно коррелировал с CRP (Spearman rho=",
  sprintf(
    "%.3f",
    score_crp_test$effect[1]
  ),
  ", P=",
  format(
    score_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; global BH-adjusted P=",
  sprintf(
    "%.4f",
    score_crp_test$BH_global[1]
  ),
  "). ",
  
  "(B) Связь непрерывного показателя SRSq с CRP (Spearman rho=",
  sprintf(
    "%.3f",
    srsq_crp_test$effect[1]
  ),
  ", P=",
  format(
    srsq_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; global BH-adjusted P=",
  sprintf(
    "%.4f",
    srsq_crp_test$BH_global[1]
  ),
  "). Эти две связи были единственными клиническими ассоциациями, сохранившими ",
  "значимость после глобальной коррекции множественных сравнений среди 60 ",
  "доступных exploratory clinical tests после интеграции демографических ",
  "переменных. Линии в панелях A и B являются только описательными визуальными ",
  "трендами; статистический вывод основан на Spearman rank correlation. ",
  
  "(C) Корреляции Spearman между SRSq и первичной пятигенной сигнатурой, ",
  "альтернативной DCAF17-сигнатурой и platform-adapted реализациями ранее ",
  "опубликованных транскриптомных сигнатур сепсиса. ",
  
  "(D) Ассоциации тех же сигнатур с Consensus Transcriptomic Subtypes (CTS), ",
  "представленные как epsilon-squared effect sizes соответствующих ",
  "Kruskal-Wallis analyses. Все benchmark analyses выполнены в одной и той же ",
  "когорте из 35 образцов крови пациентов с сепсисом. SeptiCyte-derived ",
  "RNA-seq contrasts являются исследовательскими реализациями и не являются ",
  "официальными SeptiScores. Сравнение отражает биологическую конвергенцию ",
  "сигнатур на общей molecular host-response axis, а не их диагностическое ",
  "превосходство."
)


writeLines(
  caption_ru,
  file.path(
    text_dir,
    "148_Figure3_caption_RU.txt"
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
    "148_Figure3_caption_EN_RU.txt"
  )
)


# =============================================================================
# 26. RESULTS 3.6 — READY-TO-USE ENGLISH TEXT
# =============================================================================

results_36_en <- paste0(
  
  "Among the 60 evaluable exploratory clinical-association tests retained ",
  "after demographic integration, only associations with circulating CRP ",
  "remained significant after global Benjamini-Hochberg correction. The ",
  "primary five-gene host-response score correlated positively with CRP ",
  "(Spearman rho=",
  sprintf(
    "%.3f",
    score_crp_test$effect[1]
  ),
  ", P=",
  format(
    score_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  ", global BH-adjusted P=",
  sprintf(
    "%.4f",
    score_crp_test$BH_global[1]
  ),
  "), and a concordant association was observed for SRSq (rho=",
  sprintf(
    "%.3f",
    srsq_crp_test$effect[1]
  ),
  ", P=",
  format(
    srsq_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  ", global BH-adjusted P=",
  sprintf(
    "%.4f",
    srsq_crp_test$BH_global[1]
  ),
  ") (Fig. 3A,B). No other evaluated clinical association remained significant ",
  "after correction across the complete exploratory family. These findings ",
  "link the molecular host-response axis primarily to systemic inflammatory ",
  "activity within the discovery cohort rather than establishing a prognostic ",
  "relationship."
)


writeLines(
  results_36_en,
  file.path(
    text_dir,
    "148_Results_3.6_EN.txt"
  )
)


# =============================================================================
# 27. RESULTS 3.7 — READY-TO-USE ENGLISH TEXT
# =============================================================================

primary_corr <- panel_C %>%
  dplyr::filter(
    score_column ==
      "primary_5gene_score"
  )


primary_cts <- panel_D %>%
  dplyr::filter(
    score_column ==
      "primary_5gene_score"
  )


min_rho <- min(
  panel_C$Spearman_rho,
  na.rm = TRUE
)


max_rho <- max(
  panel_C$Spearman_rho,
  na.rm = TRUE
)


min_epsilon <- min(
  panel_D$effect,
  na.rm = TRUE
)


max_epsilon <- max(
  panel_D$effect,
  na.rm = TRUE
)


results_37_en <- paste0(
  
  "To determine whether the molecular structure captured by the five-gene ",
  "score was specific to the newly derived signature or reflected a broader ",
  "sepsis host-response program, we compared it with platform-adapted ",
  "implementations of previously published blood transcriptomic signatures. ",
  "All seven evaluated scores were positively associated with SRSq, with ",
  "Spearman correlations ranging from ",
  sprintf(
    "%.3f",
    min_rho
  ),
  " to ",
  sprintf(
    "%.3f",
    max_rho
  ),
  " (Fig. 3C). The primary five-gene score showed rho=",
  sprintf(
    "%.3f",
    primary_corr$Spearman_rho[1]
  ),
  ". The same signatures also showed substantial associations with CTS, with ",
  "epsilon-squared effect sizes ranging from ",
  sprintf(
    "%.3f",
    min_epsilon
  ),
  " to ",
  sprintf(
    "%.3f",
    max_epsilon
  ),
  "; the primary five-gene score yielded epsilon-squared=",
  sprintf(
    "%.3f",
    primary_cts$effect[1]
  ),
  " (Fig. 3D). Thus, independently developed transcriptomic signatures ",
  "converged on a shared host-response axis rather than defining mutually ",
  "exclusive molecular constructs. These within-cohort comparisons were ",
  "intended to assess biological convergence and not comparative diagnostic ",
  "superiority."
)


writeLines(
  results_37_en,
  file.path(
    text_dir,
    "148_Results_3.7_EN.txt"
  )
)


# =============================================================================
# 28. RESULTS 3.6 / 3.7 — RUSSIAN WORKING TEXT
# =============================================================================

results_36_ru <- paste0(
  
  "После интеграции демографических переменных в окончательную exploratory ",
  "family вошло 60 доступных клинических ассоциационных тестов. После глобальной ",
  "коррекции Benjamini-Hochberg статистическую значимость сохранили только ",
  "ассоциации с CRP. Первичный пятигенный host-response score положительно ",
  "коррелировал с CRP (Spearman rho=",
  sprintf(
    "%.3f",
    score_crp_test$effect[1]
  ),
  ", P=",
  format(
    score_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  ", global BH-adjusted P=",
  sprintf(
    "%.4f",
    score_crp_test$BH_global[1]
  ),
  "), аналогичная связь наблюдалась для SRSq (rho=",
  sprintf(
    "%.3f",
    srsq_crp_test$effect[1]
  ),
  ", P=",
  format(
    srsq_crp_test$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  ", global BH-adjusted P=",
  sprintf(
    "%.4f",
    srsq_crp_test$BH_global[1]
  ),
  ") (рис. 3A,B). Другие клинические ассоциации не сохраняли значимость после ",
  "коррекции по всей exploratory family. Эти результаты связывают молекулярную ",
  "host-response axis преимущественно с системной воспалительной активностью, ",
  "но не устанавливают прогностическую роль сигнатуры."
)


results_37_ru <- paste0(
  
  "Для оценки того, является ли молекулярная структура, отражаемая пятигенным ",
  "score, специфичной для новой сигнатуры или представляет более общую ",
  "host-response программу сепсиса, были сопоставлены platform-adapted ",
  "реализации ранее опубликованных транскриптомных сигнатур крови. Все семь ",
  "оцененных scores положительно коррелировали с SRSq; значения Spearman rho ",
  "находились в диапазоне ",
  sprintf(
    "%.3f",
    min_rho
  ),
  "\u2013",
  sprintf(
    "%.3f",
    max_rho
  ),
  " (рис. 3C), при rho=",
  sprintf(
    "%.3f",
    primary_corr$Spearman_rho[1]
  ),
  " для первичной пятигенной сигнатуры. Те же scores демонстрировали выраженную ",
  "связь с CTS: epsilon-squared находилось в диапазоне ",
  sprintf(
    "%.3f",
    min_epsilon
  ),
  "\u2013",
  sprintf(
    "%.3f",
    max_epsilon
  ),
  ", а для первичной пятигенной сигнатуры составляло ",
  sprintf(
    "%.3f",
    primary_cts$effect[1]
  ),
  " (рис. 3D). Таким образом, независимо разработанные транскриптомные ",
  "сигнатуры сходятся на общей molecular host-response axis. Эти сравнения ",
  "внутри одной когорты оценивают биологическую конвергенцию, а не ",
  "диагностическое превосходство."
)


writeLines(
  results_36_ru,
  file.path(
    text_dir,
    "148_Results_3.6_RU.txt"
  )
)


writeLines(
  results_37_ru,
  file.path(
    text_dir,
    "148_Results_3.7_RU.txt"
  )
)


# =============================================================================
# 29. FIGURE-TO-RESULTS PLACEMENT RECORD
# =============================================================================

placement <- data.frame(
  
  Results_section = c(
    "3.6",
    "3.6",
    "3.7",
    "3.7"
  ),
  
  Panel = c(
    "3A",
    "3B",
    "3C",
    "3D"
  ),
  
  Content = c(
    "Primary five-gene score vs CRP",
    "SRSq vs CRP",
    "Published/current signatures vs SRSq",
    "Published/current signatures vs CTS"
  ),
  
  Recommended_first_citation = c(
    "End of sentence reporting five-gene score-CRP association",
    "End of sentence reporting SRSq-CRP association",
    "Sentence reporting range of signature-SRSq correlations",
    "Sentence reporting CTS epsilon-squared range"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  placement,
  file.path(
    tables_dir,
    "148_Figure3_Results_placement.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 30. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "148_sessionInfo.txt"
  )
)


# =============================================================================
# 31. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 148 completed successfully.\n")
cat("====================================================================\n\n")


cat("CLINICAL MULTIPLICITY FRAMEWORK\n")
cat("-------------------------------\n")


cat(
  "Final evaluable clinical tests: ",
  nrow(
    clinical_tests
  ),
  "\n",
  sep = ""
)


cat(
  "Global BH significant tests: ",
  nrow(
    global_sig
  ),
  "\n",
  sep = ""
)


cat("\nCRP ASSOCIATIONS\n")
cat("----------------\n")


cat(
  "Primary 5-gene vs CRP:\n",
  "  n = ",
  score_crp_test$n[1],
  "\n",
  "  rho = ",
  sprintf(
    "%.6f",
    score_crp_test$effect[1]
  ),
  "\n",
  "  P = ",
  format(
    score_crp_test$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  "  global BH = ",
  sprintf(
    "%.8f",
    score_crp_test$BH_global[1]
  ),
  "\n",
  sep = ""
)


cat(
  "SRSq vs CRP:\n",
  "  n = ",
  srsq_crp_test$n[1],
  "\n",
  "  rho = ",
  sprintf(
    "%.6f",
    srsq_crp_test$effect[1]
  ),
  "\n",
  "  P = ",
  format(
    srsq_crp_test$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  "  global BH = ",
  sprintf(
    "%.8f",
    srsq_crp_test$BH_global[1]
  ),
  "\n",
  sep = ""
)


cat("\nSRSq SIGNATURE CONVERGENCE\n")
cat("--------------------------\n")


print(
  panel_C %>%
    
    dplyr::arrange(
      dplyr::desc(
        Spearman_rho
      )
    ) %>%
    
    dplyr::select(
      figure_label,
      n,
      Spearman_rho,
      p_value,
      BH_correlation
    ),
  row.names = FALSE
)


cat("\nCTS SIGNATURE ASSOCIATIONS\n")
cat("--------------------------\n")


print(
  panel_D %>%
    
    dplyr::arrange(
      dplyr::desc(
        effect
      )
    ) %>%
    
    dplyr::select(
      figure_label,
      n,
      effect,
      p_value,
      BH_endotype
    ),
  row.names = FALSE
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
    source_workbook,
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
      "148_Figure3_caption_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nRESULTS TEXT\n")
cat("------------\n")


cat(
  normalizePath(
    file.path(
      text_dir,
      "148_Results_3.6_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat(
  normalizePath(
    file.path(
      text_dir,
      "148_Results_3.7_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nFIGURE PLACEMENT\n")
cat("----------------\n")

cat(
  "Figure 3A-B -> Results 3.6\n"
)

cat(
  "Figure 3C-D -> Results 3.7\n"
)

cat(
  "Physical placement of complete Figure 3 -> after Section 3.7\n"
)


cat("\nIMPORTANT INTERPRETATION\n")
cat("------------------------\n")

cat(
  "- Biological convergence, NOT signature superiority.\n"
)

cat(
  "- Same discovery BP cohort, n=35.\n"
)

cat(
  "- SeptiCyte-derived contrasts are NOT official SeptiScores.\n"
)

cat(
  "- No external diagnostic validation claim.\n"
)

cat(
  "- No prognostic claim from CRP association.\n"
)


cat("\nDone.\n")