################################################################################
# Script 157
# FINAL
#
# Supplementary Table S6
#
# Repeated internal cross-validation and sensitivity analyses of the
# five-gene blood host-response score
#
# Project:
#   Sepsis_DESeq2
#
# Manuscript:
#   Blood-only sepsis transcriptomic endotypes /
#   five-gene host-response signature
#
#
# PURPOSE
# -------
#
# Assemble the frozen repeated-CV results generated previously by Script 135.
#
# This script DOES NOT rerun cross-validation.
#
#
# PRIMARY FIVE-GENE SIGNATURE
# ---------------------------
#
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
# ALTERNATIVE FIVE-GENE SIGNATURE
# -------------------------------
#
# UP:
#   CD177
#   HK3
#   IRAK3
#
# DOWN:
#   CARD11
#   DCAF17
#
#
# COMPARATOR
# ----------
#
# SeptiCyte-derived four-gene RNA-expression implementation:
#
#   CEACAM4
#   LAMP1
#   PLA2G7
#   PLAC8
#
# IMPORTANT:
# This is an RNA-seq implementation of the component genes and is NOT
# the proprietary clinical SeptiCyte score.
#
#
# FROZEN REPEATED-CV DESIGN
# -------------------------
#
#   100 repeats
#   stratified 5-fold cross-validation
#
# Fixed signed-score branch:
#
#   - training-fold gene means and SDs were estimated;
#   - test samples were standardized using TRAINING-fold parameters;
#   - score = mean(z_UP) - mean(z_DOWN);
#   - no coefficients were fitted;
#   - no cutoff was optimized.
#
# Ridge-logistic sensitivity branch:
#
#   glmnet alpha = 0
#   lambda = lambda.1se selected by inner stratified CV
#   standardize = FALSE
#
#
# IMPORTANT INTERPRETATION
# ------------------------
#
# Candidate genes/panels were developed in the same discovery cohort.
# Repeated CV therefore assesses INTERNAL stability and is NOT
# independent external validation.
#
#
# EXPECTED FROZEN ANCHORS
# -----------------------
#
# Fixed signed-score repeated CV:
#
#   Primary_5_gene:
#       mean AUC approximately 1.000
#
#   DCAF17_5_gene:
#       mean AUC approximately 1.000
#
#   SeptiCyte_4_gene:
#       mean AUC approximately 0.995
#       empirical 2.5th percentile approximately 0.991
#
################################################################################


cat("====================================================================\n")
cat("Running Script 157\n")
cat("Supplementary Table S6\n")
cat("Repeated internal cross-validation and sensitivity analyses\n")
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
  library(stringr)
  library(openxlsx)
  
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "157_TableS6_internal_cross_validation"
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
  directory_name in c(
    output_dir,
    tables_dir,
    audit_dir,
    text_dir
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
# 4. HELPER FUNCTIONS
# =============================================================================

find_file_by_basename <- function(
    search_root,
    target_basename
) {
  
  hits <- list.files(
    search_root,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  
  
  hits <- hits[
    basename(
      hits
    ) ==
      target_basename
  ]
  
  
  hits <- sort(
    unique(
      hits
    )
  )
  
  
  if (
    length(
      hits
    ) ==
    0
  ) {
    
    stop(
      paste0(
        "Required file not found: ",
        target_basename
      )
    )
  }
  
  
  if (
    length(
      hits
    ) >
    1
  ) {
    
    cat(
      "\nMultiple files found for ",
      target_basename,
      ":\n",
      sep = ""
    )
    
    
    print(
      hits
    )
    
    
    stop(
      paste0(
        "Ambiguous provenance for ",
        target_basename,
        "."
      )
    )
  }
  
  
  hits[1]
}


find_column <- function(
    data,
    exact_candidates,
    regex = NULL,
    label = "column",
    required = TRUE
) {
  
  nm <- names(
    data
  )
  
  
  nm_lower <- tolower(
    nm
  )
  
  
  for (
    candidate in exact_candidates
  ) {
    
    idx <- which(
      nm_lower ==
        tolower(
          candidate
        )
    )
    
    
    if (
      length(
        idx
      ) ==
      1
    ) {
      
      return(
        nm[idx]
      )
    }
  }
  
  
  if (!is.null(regex)) {
    
    idx <- grep(
      regex,
      nm,
      ignore.case = TRUE,
      perl = TRUE
    )
    
    
    if (
      length(
        idx
      ) ==
      1
    ) {
      
      return(
        nm[idx]
      )
    }
    
    
    if (
      length(
        idx
      ) >
      1
    ) {
      
      cat(
        "\nMultiple candidates for ",
        label,
        ":\n",
        sep = ""
      )
      
      
      print(
        nm[idx]
      )
    }
  }
  
  
  if (required) {
    
    cat(
      "\nAvailable columns while searching for ",
      label,
      ":\n",
      sep = ""
    )
    
    
    print(
      nm
    )
    
    
    stop(
      paste0(
        "Could not uniquely identify ",
        label,
        "."
      )
    )
  }
  
  
  NA_character_
}


as_numeric_safe <- function(x) {
  
  suppressWarnings(
    as.numeric(
      x
    )
  )
}


standardize_panel_name <- function(x) {
  
  x2 <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  
  dplyr::case_when(
    
    grepl(
      "primary",
      x2
    ) &
      grepl(
        "5",
        x2
      ) ~
      "Primary_5_gene",
    
    grepl(
      "dcaf17",
      x2
    ) &
      grepl(
        "5",
        x2
      ) ~
      "DCAF17_5_gene",
    
    grepl(
      "septicyte",
      x2
    ) ~
      "SeptiCyte_4_gene",
    
    grepl(
      "dcaf17",
      x2
    ) &
      grepl(
        "single|1_gene|one_gene",
        x2
      ) ~
      "DCAF17_single",
    
    TRUE ~
      as.character(
        x
      )
  )
}


# =============================================================================
# 5. LOCATE SCRIPT 135 OUTPUTS
# =============================================================================

results_root <- file.path(
  project_dir,
  "results"
)


cv_summary_file <- find_file_by_basename(
  results_root,
  "135_cross_validation_summary.csv"
)


validation_workbook_hits <- list.files(
  results_root,
  pattern = "^135_blood_endotype_biomarker_validation\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)


validation_workbook_hits <- sort(
  unique(
    validation_workbook_hits
  )
)


validation_workbook_file <- if (
  length(
    validation_workbook_hits
  ) ==
  1
) {
  
  validation_workbook_hits[1]
  
} else {
  
  NA_character_
}


score_file_hits <- list.files(
  results_root,
  pattern = "^135_sepsis_blood_scores_with_SRS_CTS\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)


score_file_hits <- sort(
  unique(
    score_file_hits
  )
)


score_file <- if (
  length(
    score_file_hits
  ) ==
  1
) {
  
  score_file_hits[1]
  
} else {
  
  NA_character_
}


cat("\nSCRIPT 135 SOURCE FILES\n")
cat("-----------------------\n")


cat(
  "Cross-validation summary:\n  ",
  normalizePath(
    cv_summary_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Validation workbook:\n  ",
  if (
    !is.na(
      validation_workbook_file
    )
  ) {
    normalizePath(
      validation_workbook_file,
      winslash = "\\",
      mustWork = TRUE
    )
  } else {
    "not uniquely detected"
  },
  "\n\n",
  sep = ""
)


cat(
  "Sepsis score file:\n  ",
  if (
    !is.na(
      score_file
    )
  ) {
    normalizePath(
      score_file,
      winslash = "\\",
      mustWork = TRUE
    )
  } else {
    "not uniquely detected"
  },
  "\n",
  sep = ""
)


# =============================================================================
# 6. READ FROZEN CROSS-VALIDATION SUMMARY
# =============================================================================

cv_raw <- read.csv(
  cv_summary_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("\nCV SUMMARY DIMENSIONS\n")
cat("---------------------\n")


cat(
  nrow(
    cv_raw
  ),
  " rows x ",
  ncol(
    cv_raw
  ),
  " columns\n",
  sep = ""
)


cat("\nCV SUMMARY COLUMNS\n")
cat("------------------\n")


print(
  names(
    cv_raw
  )
)


if (
  nrow(
    cv_raw
  ) <
  3
) {
  
  stop(
    "Cross-validation summary contains unexpectedly few rows."
  )
}


# =============================================================================
# 7. IDENTIFY KEY COLUMNS
# =============================================================================

panel_col <- find_column(
  
  cv_raw,
  
  exact_candidates = c(
    "panel",
    "panel_name",
    "signature",
    "score",
    "model"
  ),
  
  regex =
    "panel|signature|score.*name|model",
  
  label =
    "panel/signature name"
)


method_col <- find_column(
  
  cv_raw,
  
  exact_candidates = c(
    "method",
    "cv_method",
    "model_type",
    "evaluation_method"
  ),
  
  regex =
    "method|model.*type|evaluation",
  
  label =
    "CV method",
  
  required =
    FALSE
)


mean_auc_col <- find_column(
  
  cv_raw,
  
  exact_candidates = c(
    "mean_auc",
    "AUC_mean",
    "mean_AUC",
    "mean_cv_auc"
  ),
  
  regex =
    "mean.*auc|auc.*mean",
  
  label =
    "mean AUC"
)


median_auc_col <- find_column(
  
  cv_raw,
  
  exact_candidates = c(
    "median_auc",
    "AUC_median",
    "median_AUC"
  ),
  
  regex =
    "median.*auc|auc.*median",
  
  label =
    "median AUC",
  
  required =
    FALSE
)


q025_auc_col <- find_column(
  
  cv_raw,
  
  exact_candidates = c(
    "q025_auc",
    "auc_q025",
    "AUC_q025",
    "p025_auc",
    "q2.5_auc"
  ),
  
  regex =
    "(q025|2\\.5|025).*auc|auc.*(q025|2\\.5|025)",
  
  label =
    "2.5th percentile AUC",
  
  required =
    FALSE
)


q975_auc_col <- find_column(
  
  cv_raw,
  
  exact_candidates = c(
    "q975_auc",
    "auc_q975",
    "AUC_q975",
    "p975_auc",
    "q97.5_auc"
  ),
  
  regex =
    "(q975|97\\.5|975).*auc|auc.*(q975|97\\.5|975)",
  
  label =
    "97.5th percentile AUC",
  
  required =
    FALSE
)


n_repeats_col <- find_column(
  
  cv_raw,
  
  exact_candidates = c(
    "n_repeats",
    "repeats",
    "n_repeat"
  ),
  
  regex =
    "repeat",
  
  label =
    "number of CV repeats",
  
  required =
    FALSE
)


cat("\nDETECTED CV COLUMN MAPPING\n")
cat("--------------------------\n")


column_mapping <- data.frame(
  
  role = c(
    "Panel/signature",
    "Method",
    "Mean AUC",
    "Median AUC",
    "AUC q2.5",
    "AUC q97.5",
    "Number of repeats"
  ),
  
  source_column = c(
    panel_col,
    ifelse(
      is.na(
        method_col
      ),
      "<not available>",
      method_col
    ),
    mean_auc_col,
    ifelse(
      is.na(
        median_auc_col
      ),
      "<not available>",
      median_auc_col
    ),
    ifelse(
      is.na(
        q025_auc_col
      ),
      "<not available>",
      q025_auc_col
    ),
    ifelse(
      is.na(
        q975_auc_col
      ),
      "<not available>",
      q975_auc_col
    ),
    ifelse(
      is.na(
        n_repeats_col
      ),
      "<not available>",
      n_repeats_col
    )
  ),
  
  stringsAsFactors = FALSE
)


print(
  column_mapping,
  row.names = FALSE
)


# =============================================================================
# 8. STANDARDIZED CV SUMMARY
# =============================================================================

cv_standardized <- data.frame(
  
  panel_original =
    as.character(
      cv_raw[[panel_col]]
    ),
  
  panel =
    standardize_panel_name(
      cv_raw[[panel_col]]
    ),
  
  method =
    if (
      !is.na(
        method_col
      )
    ) {
      as.character(
        cv_raw[[method_col]]
      )
    } else {
      "Not separately labelled in source"
    },
  
  mean_AUC =
    as_numeric_safe(
      cv_raw[[mean_auc_col]]
    ),
  
  median_AUC =
    if (
      !is.na(
        median_auc_col
      )
    ) {
      as_numeric_safe(
        cv_raw[[median_auc_col]]
      )
    } else {
      NA_real_
    },
  
  q025_AUC =
    if (
      !is.na(
        q025_auc_col
      )
    ) {
      as_numeric_safe(
        cv_raw[[q025_auc_col]]
      )
    } else {
      NA_real_
    },
  
  q975_AUC =
    if (
      !is.na(
        q975_auc_col
      )
    ) {
      as_numeric_safe(
        cv_raw[[q975_auc_col]]
      )
    } else {
      NA_real_
    },
  
  n_repeats =
    if (
      !is.na(
        n_repeats_col
      )
    ) {
      as_numeric_safe(
        cv_raw[[n_repeats_col]]
      )
    } else {
      NA_real_
    },
  
  source_row =
    seq_len(
      nrow(
        cv_raw
      )
    ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 9. PANEL ROW AUDIT
# =============================================================================

target_panels <- c(
  "Primary_5_gene",
  "DCAF17_5_gene",
  "SeptiCyte_4_gene"
)


panel_row_audit <- data.frame()

for (
  target_panel in target_panels
) {
  
  hits <- cv_standardized %>%
    
    dplyr::filter(
      panel ==
        target_panel
    )
  
  
  if (
    nrow(
      hits
    ) ==
    0
  ) {
    
    cat(
      "\nCould not find standardized row for ",
      target_panel,
      ".\n",
      sep = ""
    )
    
    
    cat(
      "Available standardized panel names:\n"
    )
    
    
    print(
      unique(
        cv_standardized$panel
      )
    )
    
    
    stop(
      paste0(
        "Missing CV result for ",
        target_panel,
        "."
      )
    )
  }
  
  
  panel_row_audit <- dplyr::bind_rows(
    panel_row_audit,
    hits
  )
}


cat("\nTARGET PANEL CV ROWS\n")
cat("--------------------\n")


print(
  panel_row_audit,
  row.names = FALSE
)


# =============================================================================
# 10. FROZEN PERFORMANCE ANCHOR AUDIT
# =============================================================================
#
# Multiple rows per panel may exist if fixed-score and ridge-logistic
# sensitivity branches are both represented.
#
# We therefore audit whether at least one row for each target panel
# reproduces the expected fixed-score result.
#
# =============================================================================

primary_hits <- panel_row_audit %>%
  
  dplyr::filter(
    panel ==
      "Primary_5_gene"
  )


dcaf_hits <- panel_row_audit %>%
  
  dplyr::filter(
    panel ==
      "DCAF17_5_gene"
  )


septi_hits <- panel_row_audit %>%
  
  dplyr::filter(
    panel ==
      "SeptiCyte_4_gene"
  )


primary_anchor_ok <- any(
  abs(
    primary_hits$mean_AUC -
      1
  ) <
    1e-8,
  na.rm = TRUE
)


dcaf_anchor_ok <- any(
  abs(
    dcaf_hits$mean_AUC -
      1
  ) <
    1e-8,
  na.rm = TRUE
)


septi_anchor_ok <- any(
  septi_hits$mean_AUC >
    0.98 &
    septi_hits$mean_AUC <=
    1,
  na.rm = TRUE
)


if (
  !primary_anchor_ok
) {
  
  stop(
    "Primary_5_gene repeated-CV mean AUC anchor was not reproduced."
  )
}


if (
  !dcaf_anchor_ok
) {
  
  stop(
    "DCAF17_5_gene repeated-CV mean AUC anchor was not reproduced."
  )
}


if (
  !septi_anchor_ok
) {
  
  stop(
    "SeptiCyte_4_gene repeated-CV mean AUC is outside expected range."
  )
}


if (
  any(
    is.finite(
      septi_hits$q025_AUC
    )
  )
) {
  
  septi_q025_ok <- any(
    septi_hits$q025_AUC >
      0.97 &
      septi_hits$q025_AUC <=
      1,
    na.rm = TRUE
  )
  
  
  if (
    !septi_q025_ok
  ) {
    
    stop(
      "SeptiCyte_4_gene empirical 2.5th-percentile AUC is outside expected range."
    )
  }
}


cat(
  "\nRepeated-CV frozen performance anchors passed.\n"
)


# =============================================================================
# 11. SCORE DEFINITIONS
# =============================================================================

score_definitions <- data.frame(
  
  signature = c(
    rep(
      "Primary_5_gene",
      5
    ),
    rep(
      "DCAF17_5_gene",
      5
    ),
    rep(
      "SeptiCyte_4_gene",
      4
    )
  ),
  
  gene = c(
    
    "CD177",
    "HK3",
    "IRAK3",
    "CARD11",
    "IKZF2",
    
    "CD177",
    "HK3",
    "IRAK3",
    "CARD11",
    "DCAF17",
    
    "CEACAM4",
    "LAMP1",
    "PLA2G7",
    "PLAC8"
  ),
  
  direction_in_score = c(
    
    "+",
    "+",
    "+",
    "-",
    "-",
    
    "+",
    "+",
    "+",
    "-",
    "-",
    
    "-",
    "+",
    "-",
    "+"
  ),
  
  manuscript_role = c(
    
    rep(
      "Primary biology-guided host-response signature",
      5
    ),
    
    rep(
      "Alternative sensitivity signature",
      5
    ),
    
    rep(
      "Published-gene comparator implementation",
      4
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 12. CROSS-VALIDATION METHOD TABLE
# =============================================================================

cv_method <- data.frame(
  
  Component = c(
    
    "Population",
    "Resampling",
    "Fixed signed-score scaling",
    "Fixed signed-score formula",
    "Coefficient fitting",
    "Threshold optimization",
    "Direction optimization",
    "Ridge sensitivity",
    "Ridge penalty selection",
    "Ridge standardization",
    "Summary interval",
    "Interpretation"
  ),
  
  Description = c(
    
    "Discovery blood cohort: 35 sepsis BP and 10 healthy-control BC samples.",
    
    "100 repeats of stratified five-fold cross-validation.",
    
    paste0(
      "Within each outer training fold, gene-specific means and standard ",
      "deviations were estimated; held-out samples were standardized using ",
      "training-fold parameters only."
    ),
    
    paste0(
      "Primary and DCAF17 signatures: mean z-score of the three UP genes ",
      "minus mean z-score of the two DOWN genes."
    ),
    
    paste0(
      "The fixed signed-score branch used no fitted regression coefficients."
    ),
    
    "No classification cutoff was optimized for the fixed-score AUC analysis.",
    
    "The prespecified score direction was retained; no test-fold direction reversal was used.",
    
    paste0(
      "Ridge logistic regression was evaluated as a separate sensitivity ",
      "analysis."
    ),
    
    "glmnet alpha=0; lambda.1se selected by inner stratified cross-validation.",
    
    "glmnet standardize=FALSE.",
    
    paste0(
      "Reported 2.5th and 97.5th percentiles are empirical percentiles across ",
      "repeated resampling runs and are not formal confidence intervals."
    ),
    
    paste0(
      "Because candidate panels were developed in the same discovery cohort, ",
      "these analyses assess internal stability rather than independent ",
      "validation."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 13. INTERNAL-VALIDATION INTERPRETIVE GUARDRAILS
# =============================================================================

interpretation_guardrails <- data.frame(
  
  Item = c(
    
    "Discovery-cohort reuse",
    "Performance saturation",
    "Clinical validation",
    "Absolute cutoff",
    "External validation"
  ),
  
  Statement = c(
    
    paste0(
      "Candidate-gene and panel development occurred in the same discovery ",
      "cohort used for internal resampling."
    ),
    
    paste0(
      "The original exhaustive search showed near-complete discrimination ",
      "saturation across eligible panels; repeated CV therefore demonstrates ",
      "stability of the strong discovery-cohort signal rather than uniqueness ",
      "of the selected five-gene configuration."
    ),
    
    paste0(
      "Internal AUC estimates should not be described as clinical diagnostic ",
      "validation."
    ),
    
    paste0(
      "The cohort-standardized z-score is not a pre-calibrated transferable ",
      "clinical threshold."
    ),
    
    paste0(
      "Independent external evaluation is reported separately for GSE154918 ",
      "and GSE185263."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 14. OPTIONAL SOURCE INVENTORY
# =============================================================================

source_inventory <- data.frame(
  
  component = c(
    "Repeated-CV summary",
    "Script 135 workbook",
    "Sepsis score table"
  ),
  
  path = c(
    
    normalizePath(
      cv_summary_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    if (
      !is.na(
        validation_workbook_file
      )
    ) {
      normalizePath(
        validation_workbook_file,
        winslash = "\\",
        mustWork = TRUE
      )
    } else {
      "Not uniquely detected"
    },
    
    if (
      !is.na(
        score_file
      )
    ) {
      normalizePath(
        score_file,
        winslash = "\\",
        mustWork = TRUE
      )
    } else {
      "Not uniquely detected"
    }
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 15. OPTIONAL WORKBOOK SHEET INVENTORY
# =============================================================================

workbook_sheet_inventory <- data.frame(
  sheet = character(),
  n_rows = integer(),
  n_columns = integer(),
  stringsAsFactors = FALSE
)


if (
  !is.na(
    validation_workbook_file
  ) &&
  requireNamespace(
    "readxl",
    quietly = TRUE
  )
) {
  
  workbook_sheets <- readxl::excel_sheets(
    validation_workbook_file
  )
  
  
  for (
    one_sheet in workbook_sheets
  ) {
    
    temp <- tryCatch(
      
      readxl::read_excel(
        validation_workbook_file,
        sheet = one_sheet
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    if (!is.null(temp)) {
      
      workbook_sheet_inventory <- dplyr::bind_rows(
        
        workbook_sheet_inventory,
        
        data.frame(
          
          sheet =
            one_sheet,
          
          n_rows =
            nrow(
              temp
            ),
          
          n_columns =
            ncol(
              temp
            ),
          
          stringsAsFactors =
            FALSE
        )
      )
    }
  }
}


# =============================================================================
# 16. TABLE S6 README
# =============================================================================

s6_readme <- data.frame(
  
  Item = c(
    
    "Title",
    "Source analysis",
    "Primary signature",
    "Alternative signature",
    "Comparator",
    "Resampling",
    "Outer-fold preprocessing",
    "Primary metric",
    "Ridge sensitivity",
    "Reported percentile range",
    "Interpretive boundary"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S6. Repeated internal cross-validation and ",
      "sensitivity analyses of the five-gene blood host-response score."
    ),
    
    "Frozen outputs from Script 135; no cross-validation is rerun here.",
    
    "CD177, HK3 and IRAK3 UP; CARD11 and IKZF2 DOWN.",
    
    "CD177, HK3 and IRAK3 UP; CARD11 and DCAF17 DOWN.",
    
    paste0(
      "SeptiCyte-related four-gene RNA-expression implementation using ",
      "CEACAM4, LAMP1, PLA2G7 and PLAC8; this is not the proprietary ",
      "clinical SeptiCyte score."
    ),
    
    "100 repeats of stratified five-fold cross-validation.",
    
    paste0(
      "Gene-level z-standardization for held-out samples used parameters ",
      "estimated from the corresponding training fold."
    ),
    
    "Area under the ROC curve.",
    
    paste0(
      "Ridge logistic regression was evaluated as a separate sensitivity ",
      "analysis using alpha=0 and lambda.1se from inner CV."
    ),
    
    paste0(
      "Empirical 2.5th and 97.5th percentiles across repeated resampling ",
      "runs; these are not formal confidence intervals."
    ),
    
    paste0(
      "Internal cross-validation assesses stability within the discovery ",
      "cohort and is not independent external or clinical validation."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 17. OUTPUT PATHS
# =============================================================================

submission_file <- file.path(
  tables_dir,
  "157_TableS6_repeated_internal_cross_validation.xlsx"
)


audit_file <- file.path(
  audit_dir,
  "157_INTERNAL_AUDIT_TableS6_cross_validation.xlsx"
)


note_file <- file.path(
  text_dir,
  "157_TableS6_title_and_note_EN.txt"
)


# =============================================================================
# 18. EXCEL STYLES
# =============================================================================

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE,
  valign = "center"
)


readme_header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#EEF3F7",
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE
)


# =============================================================================
# 19. WRITE SUBMISSION WORKBOOK
# =============================================================================

wb <- openxlsx::createWorkbook()


submission_objects <- list(
  
  S6_ReadMe =
    s6_readme,
  
  CV_summary =
    cv_standardized,
  
  CV_source_complete =
    cv_raw,
  
  Score_definitions =
    score_definitions,
  
  CV_method =
    cv_method,
  
  Interpretation =
    interpretation_guardrails
)


for (
  sheet_name in names(
    submission_objects
  )
) {
  
  data_object <- submission_objects[[sheet_name]]
  
  
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
    if (
      sheet_name ==
      "S6_ReadMe"
    ) {
      readme_header_style
    } else {
      header_style
    },
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


openxlsx::setColWidths(
  wb,
  "S6_ReadMe",
  cols = 1,
  widths = 30
)


openxlsx::setColWidths(
  wb,
  "S6_ReadMe",
  cols = 2,
  widths = 95
)


openxlsx::setColWidths(
  wb,
  "CV_method",
  cols = 2,
  widths = 95
)


openxlsx::setColWidths(
  wb,
  "Interpretation",
  cols = 2,
  widths = 95
)


openxlsx::saveWorkbook(
  wb,
  submission_file,
  overwrite = TRUE
)


# =============================================================================
# 20. WRITE INTERNAL AUDIT WORKBOOK
# =============================================================================

wb_audit <- openxlsx::createWorkbook()


audit_objects <- list(
  
  Column_mapping =
    column_mapping,
  
  Target_panel_rows =
    panel_row_audit,
  
  Source_inventory =
    source_inventory,
  
  Workbook_inventory =
    workbook_sheet_inventory
)


for (
  sheet_name in names(
    audit_objects
  )
) {
  
  data_object <- audit_objects[[sheet_name]]
  
  
  openxlsx::addWorksheet(
    wb_audit,
    sheet_name
  )
  
  
  openxlsx::writeData(
    wb_audit,
    sheet_name,
    data_object,
    withFilter = TRUE
  )
  
  
  if (
    ncol(
      data_object
    ) >
    0
  ) {
    
    openxlsx::addStyle(
      wb_audit,
      sheet_name,
      header_style,
      rows = 1,
      cols = 1:ncol(
        data_object
      ),
      gridExpand = TRUE
    )
    
    
    openxlsx::setColWidths(
      wb_audit,
      sheet_name,
      cols = 1:ncol(
        data_object
      ),
      widths = "auto"
    )
  }
  
  
  openxlsx::freezePane(
    wb_audit,
    sheet_name,
    firstActiveRow = 2
  )
}


openxlsx::saveWorkbook(
  wb_audit,
  audit_file,
  overwrite = TRUE
)


# =============================================================================
# 21. TABLE NOTE
# =============================================================================

table_note <- c(
  
  paste0(
    "Supplementary Table S6. Repeated internal cross-validation and ",
    "sensitivity analyses of the five-gene blood host-response score."
  ),
  
  "",
  
  paste0(
    "Internal stability was assessed using 100 repeats of stratified ",
    "five-fold cross-validation. For the fixed signed-score branch, ",
    "gene-specific scaling parameters were estimated within each training ",
    "fold and applied to the held-out samples. No regression coefficient or ",
    "classification cutoff was fitted for this branch. Ridge logistic ",
    "regression was evaluated separately as a sensitivity analysis. Because ",
    "candidate panels were developed in the same discovery cohort, repeated ",
    "cross-validation estimates should be interpreted as internal stability ",
    "rather than independent clinical validation."
  )
)


writeLines(
  table_note,
  note_file
)


# =============================================================================
# 22. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "157_sessionInfo.txt"
  )
)


# =============================================================================
# 23. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 157 completed successfully.\n")
cat("====================================================================\n\n")


cat("CROSS-VALIDATION SOURCE\n")
cat("-----------------------\n")


cat(
  "Rows = ",
  nrow(
    cv_raw
  ),
  "\n",
  sep = ""
)


cat(
  "Columns = ",
  ncol(
    cv_raw
  ),
  "\n",
  sep = ""
)


cat("\nSTANDARDIZED CV SUMMARY\n")
cat("-----------------------\n")


print(
  cv_standardized,
  row.names = FALSE
)


cat("\nTARGET PANELS\n")
cat("-------------\n")


print(
  panel_row_audit,
  row.names = FALSE
)


cat("\nFROZEN ANCHOR STATUS\n")
cat("--------------------\n")


cat(
  "Primary_5_gene mean AUC = 1 reproduced: ",
  primary_anchor_ok,
  "\n",
  sep = ""
)


cat(
  "DCAF17_5_gene mean AUC = 1 reproduced: ",
  dcaf_anchor_ok,
  "\n",
  sep = ""
)


cat(
  "SeptiCyte_4_gene mean AUC in expected range: ",
  septi_anchor_ok,
  "\n",
  sep = ""
)


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")


cat(
  "- 100 x stratified 5-fold CV is INTERNAL validation only.\n"
)


cat(
  "- Test-fold z-scaling uses training-fold parameters.\n"
)


cat(
  "- Fixed signed score has no fitted coefficients.\n"
)


cat(
  "- No clinical cutoff is optimized.\n"
)


cat(
  "- Ridge logistic regression is a sensitivity analysis.\n"
)


cat(
  "- Empirical q2.5/q97.5 are not formal confidence intervals.\n"
)


cat(
  "- Near-perfect internal discrimination does not establish unique panel optimality.\n"
)


cat(
  "- Independent external evaluation is reported separately.\n"
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Supplementary Table S6:\n  ",
  normalizePath(
    submission_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Internal audit:\n  ",
  normalizePath(
    audit_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Table title/note:\n  ",
  normalizePath(
    note_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nDone.\n")