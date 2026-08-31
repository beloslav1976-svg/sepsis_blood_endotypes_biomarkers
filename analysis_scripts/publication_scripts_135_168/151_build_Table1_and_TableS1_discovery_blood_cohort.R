################################################################################
# Script 151
# FINAL
#
# Main Table 1 + Supplementary Table S1
#
# Project:
# Sepsis_DESeq2
#
# Manuscript:
# Blood-only sepsis transcriptomic endotypes /
# five-gene host-response signature
#
#
# INPUT STRATEGY
# --------------
#
# Input Excel files:
#   readxl
#
# Output Excel files:
#   openxlsx
#
#
# INPUT FILES
# -----------
#
# 1. Table_S1_Deidentified_RNAseq_Sample_Annotation_REVISED.xlsx
#
#    Primary de-identified source for:
#      - participant_id
#      - cohort
#      - sex
#      - corrected age
#      - blood RNA-seq sample ID
#      - outcome
#      - mechanical ventilation
#      - culture status
#      - serum creatinine
#      - blood sequencing batch
#
#
# 2. Table_S1_Deidentified_RNAseq_Sample_Annotation.xlsx
#
#    Additional de-identified clinical source:
#      - CRP
#      - lactate
#      - albumin
#      - procalcitonin
#      - ALT
#      - AST
#      - WBC
#      - platelets
#
#
# 3. Выборка Сепсис.xlsx
#
#    Optional raw-source integrity audit only.
#
#    Direct identifiers from this file are NEVER exported.
#
#
# DISCOVERY BLOOD COHORT
# ----------------------
#
# Sepsis:
#   n = 35
#
# Healthy controls:
#   n = 10
#
# Blood RNA-seq:
#   BP = 35
#   BC = 10
#
#
# CROSS-AUDIT LOGIC
# -----------------
#
# Across all 45 participants:
#   - cohort
#   - sex
#   - blood_rnaseq_sample_id
#   - serum_creatinine_umol_l
#
# Across sepsis participants only (n = 35):
#   - outcome_status
#   - mechanical_ventilation_status
#   - culture_status
#
# This avoids treating blank vs NA acute-care fields among healthy controls
# as source-data discrepancies.
#
#
# IMPORTANT REPORTING RULES
# -------------------------
#
# - blood-only outputs
# - no urine metadata exported
# - no direct identifiers exported
# - no missing values imputed
# - SOFA is not reported because source SOFA values are unavailable (0/35)
#
################################################################################


cat("====================================================================\n")
cat("Running Script 151\n")
cat("Main Table 1 + Supplementary Table S1\n")
cat("Discovery blood transcriptomic cohort\n")
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
  library(stringr)
  library(readxl)
  library(openxlsx)
  
})


# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

find_input_file <- function(
    filename,
    required = TRUE
) {
  
  candidates <- c(
    
    file.path(
      project_dir,
      "data",
      "clinical",
      filename
    ),
    
    file.path(
      project_dir,
      "data",
      "metadata",
      filename
    ),
    
    file.path(
      project_dir,
      "data",
      filename
    ),
    
    file.path(
      project_dir,
      filename
    ),
    
    file.path(
      path.expand("~/Downloads"),
      filename
    )
  )
  
  
  existing <- candidates[
    file.exists(candidates)
  ]
  
  
  if (length(existing) > 0) {
    
    return(
      existing[1]
    )
  }
  
  
  all_project_files <- list.files(
    project_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  
  
  hits <- all_project_files[
    basename(all_project_files) ==
      filename
  ]
  
  
  if (length(hits) > 0) {
    
    hits <- sort(
      unique(hits)
    )
    
    
    if (length(hits) > 1) {
      
      cat(
        "\nMultiple copies found for ",
        filename,
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
        "Using:\n  ",
        hits[1],
        "\n",
        sep = ""
      )
    }
    
    
    return(
      hits[1]
    )
  }
  
  
  if (required) {
    
    stop(
      paste0(
        "Required input file not found:\n",
        filename,
        "\n\nRecommended location:\n",
        file.path(
          project_dir,
          "data",
          "clinical",
          filename
        )
      )
    )
  }
  
  
  return(
    NA_character_
  )
}


check_required_columns <- function(
    data,
    required_columns,
    table_name
) {
  
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        "Missing required column(s) in ",
        table_name,
        ":\n",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
}


compare_character <- function(
    x,
    y
) {
  
  x2 <- ifelse(
    is.na(x),
    "<NA>",
    stringr::str_trim(
      as.character(x)
    )
  )
  
  
  y2 <- ifelse(
    is.na(y),
    "<NA>",
    stringr::str_trim(
      as.character(y)
    )
  )
  
  
  x2 == y2
}


compare_numeric <- function(
    x,
    y,
    tolerance = 1e-10
) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  
  y <- suppressWarnings(
    as.numeric(y)
  )
  
  
  result <- rep(
    FALSE,
    length(x)
  )
  
  
  both_na <-
    is.na(x) &
    is.na(y)
  
  
  both_finite <-
    is.finite(x) &
    is.finite(y)
  
  
  result[both_na] <- TRUE
  
  
  result[both_finite] <-
    abs(
      x[both_finite] -
        y[both_finite]
    ) <=
    tolerance
  
  
  result
}


format_number <- function(
    x,
    digits = 1
) {
  
  if (!is.finite(x)) {
    return("NA")
  }
  
  
  formatC(
    x,
    format = "f",
    digits = digits,
    big.mark = ","
  )
}


median_iqr_string <- function(
    x,
    digits = 1
) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  
  x <- x[
    is.finite(x)
  ]
  
  
  if (length(x) == 0) {
    
    return(
      "NA"
    )
  }
  
  
  qq <- stats::quantile(
    x,
    probs = c(
      0.25,
      0.50,
      0.75
    ),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  
  
  paste0(
    format_number(
      qq[2],
      digits
    ),
    " [",
    format_number(
      qq[1],
      digits
    ),
    "\u2013",
    format_number(
      qq[3],
      digits
    ),
    "]"
  )
}


n_percent_string <- function(
    n,
    denominator,
    digits = 1
) {
  
  paste0(
    n,
    " (",
    formatC(
      100 * n / denominator,
      format = "f",
      digits = digits
    ),
    "%)"
  )
}


format_p_table <- function(p) {
  
  if (!is.finite(p)) {
    
    return(
      "\u2014"
    )
  }
  
  
  if (p < 0.001) {
    
    return(
      format(
        p,
        scientific = TRUE,
        digits = 3
      )
    )
  }
  
  
  sprintf(
    "%.4f",
    p
  )
}


# =============================================================================
# 4. INPUT FILES
# =============================================================================

revised_file <- find_input_file(
  "Table_S1_Deidentified_RNAseq_Sample_Annotation_REVISED.xlsx",
  required = TRUE
)


extended_file <- find_input_file(
  "Table_S1_Deidentified_RNAseq_Sample_Annotation.xlsx",
  required = TRUE
)


raw_file <- find_input_file(
  "Выборка Сепсис.xlsx",
  required = FALSE
)


cat("\n====================================================================\n")
cat("INPUT FILES\n")
cat("====================================================================\n")


cat(
  "REVISED metadata:\n  ",
  normalizePath(
    revised_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Extended de-identified metadata:\n  ",
  normalizePath(
    extended_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


if (
  !is.na(raw_file) &&
  file.exists(raw_file)
) {
  
  cat(
    "Raw clinical workbook:\n  ",
    normalizePath(
      raw_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    "\n",
    sep = ""
  )
  
} else {
  
  cat(
    "Raw clinical workbook not found.\n",
    "Optional raw-source audit will be skipped.\n",
    sep = ""
  )
}


# =============================================================================
# 5. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "151_Table1_TableS1_discovery_cohort"
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


for (
  directory_name in c(
    output_dir,
    tables_dir,
    text_dir,
    logs_dir
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
# 6. SHEET AUDIT — READXL
# =============================================================================

revised_sheets <- readxl::excel_sheets(
  revised_file
)


extended_sheets <- readxl::excel_sheets(
  extended_file
)


cat("\nREVISED workbook sheets:\n")

print(
  revised_sheets
)


cat("\nExtended workbook sheets:\n")

print(
  extended_sheets
)


required_revised_sheets <- c(
  "Participant_Metadata",
  "Sample_Metadata"
)


missing_revised_sheets <- setdiff(
  required_revised_sheets,
  revised_sheets
)


if (length(missing_revised_sheets) > 0) {
  
  stop(
    paste0(
      "Missing required REVISED sheet(s): ",
      paste(
        missing_revised_sheets,
        collapse = ", "
      )
    )
  )
}


if (
  !("Participant_Metadata" %in%
    extended_sheets)
) {
  
  stop(
    "Participant_Metadata sheet missing from extended workbook."
  )
}


# =============================================================================
# 7. READ INPUT TABLES — READXL
# =============================================================================

participant_revised <- readxl::read_excel(
  revised_file,
  sheet = "Participant_Metadata"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )


sample_revised <- readxl::read_excel(
  revised_file,
  sheet = "Sample_Metadata"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )


participant_extended <- readxl::read_excel(
  extended_file,
  sheet = "Participant_Metadata"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )


cat("\nInput dimensions:\n")


cat(
  "REVISED Participant_Metadata: ",
  nrow(participant_revised),
  " x ",
  ncol(participant_revised),
  "\n",
  sep = ""
)


cat(
  "REVISED Sample_Metadata: ",
  nrow(sample_revised),
  " x ",
  ncol(sample_revised),
  "\n",
  sep = ""
)


cat(
  "Extended Participant_Metadata: ",
  nrow(participant_extended),
  " x ",
  ncol(participant_extended),
  "\n",
  sep = ""
)


# =============================================================================
# 8. REQUIRED COLUMN AUDIT
# =============================================================================

required_revised_participant_columns <- c(
  "participant_id",
  "cohort",
  "sex",
  "age_years",
  "blood_rnaseq_sample_id",
  "outcome_status",
  "mechanical_ventilation_status",
  "culture_status",
  "serum_creatinine_umol_l",
  "creatinine_group_110_umol_l"
)


required_revised_sample_columns <- c(
  "sample_id",
  "participant_id",
  "cohort",
  "biofluid",
  "sample_group",
  "sequencing_batch"
)


required_extended_columns <- c(
  "participant_id",
  "original_study_code",
  "cohort",
  "sex",
  "age_years",
  "blood_rnaseq_sample_id",
  "outcome_status",
  "mechanical_ventilation_status",
  "culture_status",
  "serum_creatinine_umol_l",
  "crp_mg_l",
  "lactate_mmol_l",
  "albumin_g_l",
  "procalcitonin_ng_ml",
  "alt_u_l",
  "ast_u_l",
  "wbc_10e9_l",
  "platelets_10e9_l",
  "sofa_score"
)


check_required_columns(
  participant_revised,
  required_revised_participant_columns,
  "REVISED Participant_Metadata"
)


check_required_columns(
  sample_revised,
  required_revised_sample_columns,
  "REVISED Sample_Metadata"
)


check_required_columns(
  participant_extended,
  required_extended_columns,
  "Extended Participant_Metadata"
)


# =============================================================================
# 9. CLEAN AND AUDIT PARTICIPANTS
# =============================================================================

participant_revised <- participant_revised %>%
  
  dplyr::filter(
    !is.na(participant_id),
    participant_id != ""
  )


participant_extended <- participant_extended %>%
  
  dplyr::filter(
    !is.na(participant_id),
    participant_id != ""
  )


if (
  nrow(participant_revised) != 45
) {
  
  stop(
    paste0(
      "Expected 45 participants in REVISED metadata; observed ",
      nrow(participant_revised)
    )
  )
}


if (
  nrow(participant_extended) != 45
) {
  
  stop(
    paste0(
      "Expected 45 participants in extended metadata; observed ",
      nrow(participant_extended)
    )
  )
}


if (
  dplyr::n_distinct(
    participant_revised$participant_id
  ) != 45
) {
  
  stop(
    "Duplicate participant IDs detected in REVISED metadata."
  )
}


if (
  dplyr::n_distinct(
    participant_extended$participant_id
  ) != 45
) {
  
  stop(
    "Duplicate participant IDs detected in extended metadata."
  )
}


if (
  !setequal(
    participant_revised$participant_id,
    participant_extended$participant_id
  )
) {
  
  stop(
    "Participant IDs differ between REVISED and extended metadata."
  )
}


# =============================================================================
# 10. CROSS-AUDIT SHARED DE-IDENTIFIED FIELDS
# =============================================================================
#
# ALL 45 PARTICIPANTS:
#   cohort
#   sex
#   blood_rnaseq_sample_id
#   serum_creatinine_umol_l
#
# SEPSIS ONLY:
#   outcome_status
#   mechanical_ventilation_status
#   culture_status
#
# Acute-care variables are not applicable to healthy controls.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 10.1 Fields audited across all 45 participants
# -----------------------------------------------------------------------------

shared_fields_all <- c(
  "cohort",
  "sex",
  "blood_rnaseq_sample_id",
  "serum_creatinine_umol_l"
)


audit_revised_all <- participant_revised %>%
  
  dplyr::select(
    participant_id,
    dplyr::all_of(
      shared_fields_all
    )
  ) %>%
  
  dplyr::rename_with(
    ~ paste0(
      .x,
      "_revised"
    ),
    -participant_id
  )


audit_extended_all <- participant_extended %>%
  
  dplyr::select(
    participant_id,
    dplyr::all_of(
      shared_fields_all
    )
  ) %>%
  
  dplyr::rename_with(
    ~ paste0(
      .x,
      "_extended"
    ),
    -participant_id
  )


shared_audit_all <- dplyr::inner_join(
  audit_revised_all,
  audit_extended_all,
  by = "participant_id"
)


shared_field_audit_all <- data.frame(
  
  field =
    shared_fields_all,
  
  audit_scope =
    "All participants",
  
  matching_records =
    NA_integer_,
  
  total_records =
    nrow(
      shared_audit_all
    ),
  
  stringsAsFactors = FALSE
)


for (
  i in seq_len(
    nrow(
      shared_field_audit_all
    )
  )
) {
  
  field_name <-
    shared_field_audit_all$field[i]
  
  
  revised_column_name <- paste0(
    field_name,
    "_revised"
  )
  
  
  extended_column_name <- paste0(
    field_name,
    "_extended"
  )
  
  
  revised_values <-
    shared_audit_all[[revised_column_name]]
  
  
  extended_values <-
    shared_audit_all[[extended_column_name]]
  
  
  if (
    field_name ==
    "serum_creatinine_umol_l"
  ) {
    
    matches <- compare_numeric(
      revised_values,
      extended_values
    )
    
  } else {
    
    matches <- compare_character(
      revised_values,
      extended_values
    )
  }
  
  
  shared_field_audit_all$matching_records[i] <-
    sum(
      matches
    )
}


shared_field_audit_all$all_match <-
  shared_field_audit_all$matching_records ==
  shared_field_audit_all$total_records


# -----------------------------------------------------------------------------
# 10.2 Acute-care fields audited in sepsis only
# -----------------------------------------------------------------------------

shared_fields_sepsis <- c(
  "outcome_status",
  "mechanical_ventilation_status",
  "culture_status"
)


audit_revised_sepsis <- participant_revised %>%
  
  dplyr::filter(
    cohort ==
      "Sepsis"
  ) %>%
  
  dplyr::select(
    participant_id,
    dplyr::all_of(
      shared_fields_sepsis
    )
  ) %>%
  
  dplyr::rename_with(
    ~ paste0(
      .x,
      "_revised"
    ),
    -participant_id
  )


audit_extended_sepsis <- participant_extended %>%
  
  dplyr::filter(
    cohort ==
      "Sepsis"
  ) %>%
  
  dplyr::select(
    participant_id,
    dplyr::all_of(
      shared_fields_sepsis
    )
  ) %>%
  
  dplyr::rename_with(
    ~ paste0(
      .x,
      "_extended"
    ),
    -participant_id
  )


shared_audit_sepsis <- dplyr::inner_join(
  audit_revised_sepsis,
  audit_extended_sepsis,
  by = "participant_id"
)


if (
  nrow(shared_audit_sepsis) != 35
) {
  
  stop(
    paste0(
      "Expected 35 sepsis participants for acute-care audit; observed ",
      nrow(
        shared_audit_sepsis
      )
    )
  )
}


shared_field_audit_sepsis <- data.frame(
  
  field =
    shared_fields_sepsis,
  
  audit_scope =
    "Sepsis only",
  
  matching_records =
    NA_integer_,
  
  total_records =
    nrow(
      shared_audit_sepsis
    ),
  
  stringsAsFactors = FALSE
)


for (
  i in seq_len(
    nrow(
      shared_field_audit_sepsis
    )
  )
) {
  
  field_name <-
    shared_field_audit_sepsis$field[i]
  
  
  revised_column_name <- paste0(
    field_name,
    "_revised"
  )
  
  
  extended_column_name <- paste0(
    field_name,
    "_extended"
  )
  
  
  revised_values <-
    shared_audit_sepsis[[revised_column_name]]
  
  
  extended_values <-
    shared_audit_sepsis[[extended_column_name]]
  
  
  matches <- compare_character(
    revised_values,
    extended_values
  )
  
  
  shared_field_audit_sepsis$matching_records[i] <-
    sum(
      matches
    )
}


shared_field_audit_sepsis$all_match <-
  shared_field_audit_sepsis$matching_records ==
  shared_field_audit_sepsis$total_records


# -----------------------------------------------------------------------------
# 10.3 Combined source audit
# -----------------------------------------------------------------------------

shared_field_audit <- dplyr::bind_rows(
  shared_field_audit_all,
  shared_field_audit_sepsis
)


cat("\nDe-identified source cross-audit:\n")


print(
  shared_field_audit,
  row.names = FALSE
)


if (
  !all(
    shared_field_audit$all_match
  )
) {
  
  stop(
    "REVISED vs extended metadata cross-audit failed."
  )
}


cat(
  "\nREVISED vs extended metadata cross-audit passed successfully.\n"
)


# =============================================================================
# 11. SEPSIS AGE CROSS-AUDIT
# =============================================================================

age_revised <- participant_revised %>%
  
  dplyr::filter(
    cohort ==
      "Sepsis"
  ) %>%
  
  dplyr::select(
    participant_id,
    age_years
  ) %>%
  
  dplyr::rename(
    age_revised =
      age_years
  )


age_extended <- participant_extended %>%
  
  dplyr::filter(
    cohort ==
      "Sepsis"
  ) %>%
  
  dplyr::select(
    participant_id,
    age_years
  ) %>%
  
  dplyr::rename(
    age_extended =
      age_years
  )


age_audit <- dplyr::left_join(
  age_revised,
  age_extended,
  by = "participant_id"
)


age_matches <- compare_numeric(
  age_audit$age_revised,
  age_audit$age_extended
)


if (
  !all(
    age_matches
  )
) {
  
  stop(
    "Sepsis age mismatch between REVISED and extended metadata."
  )
}


cat(
  "Sepsis age cross-audit passed successfully.\n"
)


# =============================================================================
# 12. BUILD FINAL BLOOD PARTICIPANT DATASET
# =============================================================================

additional_clinical <- participant_extended %>%
  
  dplyr::select(
    participant_id,
    crp_mg_l,
    lactate_mmol_l,
    albumin_g_l,
    procalcitonin_ng_ml,
    alt_u_l,
    ast_u_l,
    wbc_10e9_l,
    platelets_10e9_l
  )


participant_blood <- participant_revised %>%
  
  dplyr::select(
    participant_id,
    cohort,
    sex,
    age_years,
    blood_rnaseq_sample_id,
    outcome_status,
    mechanical_ventilation_status,
    culture_status,
    serum_creatinine_umol_l,
    creatinine_group_110_umol_l
  ) %>%
  
  dplyr::left_join(
    additional_clinical,
    by = "participant_id"
  ) %>%
  
  dplyr::mutate(
    
    age_years =
      suppressWarnings(
        as.numeric(
          age_years
        )
      ),
    
    serum_creatinine_umol_l =
      suppressWarnings(
        as.numeric(
          serum_creatinine_umol_l
        )
      ),
    
    crp_mg_l =
      suppressWarnings(
        as.numeric(
          crp_mg_l
        )
      ),
    
    lactate_mmol_l =
      suppressWarnings(
        as.numeric(
          lactate_mmol_l
        )
      ),
    
    albumin_g_l =
      suppressWarnings(
        as.numeric(
          albumin_g_l
        )
      ),
    
    procalcitonin_ng_ml =
      suppressWarnings(
        as.numeric(
          procalcitonin_ng_ml
        )
      ),
    
    alt_u_l =
      suppressWarnings(
        as.numeric(
          alt_u_l
        )
      ),
    
    ast_u_l =
      suppressWarnings(
        as.numeric(
          ast_u_l
        )
      ),
    
    wbc_10e9_l =
      suppressWarnings(
        as.numeric(
          wbc_10e9_l
        )
      ),
    
    platelets_10e9_l =
      suppressWarnings(
        as.numeric(
          platelets_10e9_l
        )
      )
  )


# =============================================================================
# 13. BLOOD SAMPLE METADATA
# =============================================================================

blood_samples <- sample_revised %>%
  
  dplyr::filter(
    biofluid ==
      "Whole blood"
  ) %>%
  
  dplyr::select(
    sample_id,
    participant_id,
    cohort,
    biofluid,
    sample_group,
    sequencing_batch
  )


if (
  nrow(
    blood_samples
  ) != 45
) {
  
  stop(
    paste0(
      "Expected 45 blood RNA-seq samples; observed ",
      nrow(
        blood_samples
      )
    )
  )
}


if (
  dplyr::n_distinct(
    blood_samples$sample_id
  ) != 45
) {
  
  stop(
    "Duplicate blood RNA-seq sample IDs detected."
  )
}


n_BP <- sum(
  blood_samples$sample_group ==
    "BP",
  na.rm = TRUE
)


n_BC <- sum(
  blood_samples$sample_group ==
    "BC",
  na.rm = TRUE
)


if (
  n_BP != 35
) {
  
  stop(
    "Expected BP n=35."
  )
}


if (
  n_BC != 10
) {
  
  stop(
    "Expected BC n=10."
  )
}


if (
  !setequal(
    blood_samples$sample_id,
    participant_blood$blood_rnaseq_sample_id
  )
) {
  
  stop(
    "Participant blood RNA-seq IDs and Sample_Metadata IDs do not match."
  )
}


# =============================================================================
# 14. DEFINE DISCOVERY COHORTS
# =============================================================================

sepsis <- participant_blood %>%
  
  dplyr::filter(
    cohort ==
      "Sepsis"
  )


controls <- participant_blood %>%
  
  dplyr::filter(
    cohort ==
      "Healthy control"
  )


if (
  nrow(
    sepsis
  ) != 35
) {
  
  stop(
    paste0(
      "Expected sepsis n=35; observed ",
      nrow(
        sepsis
      )
    )
  )
}


if (
  nrow(
    controls
  ) != 10
) {
  
  stop(
    paste0(
      "Expected healthy-control n=10; observed ",
      nrow(
        controls
      )
    )
  )
}


# =============================================================================
# 15. OPTIONAL RAW-WORKBOOK SOURCE AUDIT
# =============================================================================
#
# Raw source:
#   Выборка Сепсис.xlsx
#
# Sheet verified in source workbook:
#   Пациенты Сепсис
#
# Raw/direct identifiers are used only transiently for source-integrity
# matching and are NEVER exported.
#
# =============================================================================

raw_audit_summary <- data.frame(
  
  metric = c(
    "Raw source workbook found",
    "Raw sepsis records evaluated",
    "Study-code mapping complete",
    "Numeric cells checked",
    "Numeric mismatches",
    "Categorical cells checked",
    "Categorical mismatches"
  ),
  
  value = c(
    "No",
    NA,
    NA,
    NA,
    NA,
    NA,
    NA
  ),
  
  stringsAsFactors = FALSE
)


if (
  !is.na(
    raw_file
  ) &&
  file.exists(
    raw_file
  )
) {
  
  raw_sheets <- readxl::excel_sheets(
    raw_file
  )
  
  
  cat("\nRaw workbook sheets:\n")
  
  print(
    raw_sheets
  )
  
  
  if (
    "Пациенты Сепсис" %in%
    raw_sheets
  ) {
    
    raw_sepsis <- readxl::read_excel(
      raw_file,
      sheet = "Пациенты Сепсис"
    ) %>%
      as.data.frame(
        stringsAsFactors = FALSE
      )
    
    
    required_raw_columns <- c(
      "#P",
      "Исходы/дата",
      "Обнаруженные патогены",
      "Креатинин",
      "ИВЛ",
      "crp_mg_l",
      "lactate_mmol_l",
      "albumin_g_l",
      "procalcitonin_ng_ml",
      "alt_u_l",
      "ast_u_l",
      "wbc_10e9_l",
      "platelets_10e9_l"
    )
    
    
    missing_raw_columns <- setdiff(
      required_raw_columns,
      names(
        raw_sepsis
      )
    )
    
    
    if (
      length(
        missing_raw_columns
      ) ==
      0
    ) {
      
      raw_clean <- raw_sepsis %>%
        
        dplyr::filter(
          !is.na(
            .data[["#P"]]
          ),
          .data[["#P"]] != ""
        ) %>%
        
        dplyr::transmute(
          
          original_study_code =
            as.character(
              .data[["#P"]]
            ),
          
          serum_creatinine_umol_l =
            suppressWarnings(
              as.numeric(
                .data[["Креатинин"]]
              )
            ),
          
          crp_mg_l =
            suppressWarnings(
              as.numeric(
                .data[["crp_mg_l"]]
              )
            ),
          
          lactate_mmol_l =
            suppressWarnings(
              as.numeric(
                .data[["lactate_mmol_l"]]
              )
            ),
          
          albumin_g_l =
            suppressWarnings(
              as.numeric(
                .data[["albumin_g_l"]]
              )
            ),
          
          procalcitonin_ng_ml =
            suppressWarnings(
              as.numeric(
                .data[["procalcitonin_ng_ml"]]
              )
            ),
          
          alt_u_l =
            suppressWarnings(
              as.numeric(
                .data[["alt_u_l"]]
              )
            ),
          
          ast_u_l =
            suppressWarnings(
              as.numeric(
                .data[["ast_u_l"]]
              )
            ),
          
          wbc_10e9_l =
            suppressWarnings(
              as.numeric(
                .data[["wbc_10e9_l"]]
              )
            ),
          
          platelets_10e9_l =
            suppressWarnings(
              as.numeric(
                .data[["platelets_10e9_l"]]
              )
            ),
          
          raw_outcome_status =
            dplyr::if_else(
              stringr::str_detect(
                stringr::str_to_lower(
                  stringr::str_trim(
                    as.character(
                      .data[["Исходы/дата"]]
                    )
                  )
                ),
                "^смерть"
              ),
              "Death",
              "Discharged"
            ),
          
          raw_ventilation_status =
            dplyr::if_else(
              suppressWarnings(
                as.numeric(
                  .data[["ИВЛ"]]
                )
              ) > 0,
              "Yes",
              "No"
            ),
          
          raw_culture_status =
            dplyr::if_else(
              stringr::str_detect(
                stringr::str_to_lower(
                  stringr::str_trim(
                    as.character(
                      .data[["Обнаруженные патогены"]]
                    )
                  )
                ),
                "^роста нет"
              ),
              "Negative",
              "Positive"
            )
        )
      
      
      if (
        nrow(
          raw_clean
        ) != 35
      ) {
        
        stop(
          paste0(
            "Expected 35 non-empty raw sepsis records; observed ",
            nrow(
              raw_clean
            )
          )
        )
      }
      
      
      extended_raw_key <- participant_extended %>%
        
        dplyr::filter(
          cohort ==
            "Sepsis"
        ) %>%
        
        dplyr::select(
          original_study_code,
          serum_creatinine_umol_l,
          crp_mg_l,
          lactate_mmol_l,
          albumin_g_l,
          procalcitonin_ng_ml,
          alt_u_l,
          ast_u_l,
          wbc_10e9_l,
          platelets_10e9_l,
          outcome_status,
          mechanical_ventilation_status,
          culture_status
        )
      
      
      raw_compare <- dplyr::left_join(
        extended_raw_key,
        raw_clean,
        by = "original_study_code",
        suffix = c(
          "_deid",
          "_raw"
        )
      )
      
      
      raw_source_codes <-
        extended_raw_key$original_study_code
      
      
      mapping_complete <-
        nrow(
          raw_compare
        ) ==
        35 &&
        all(
          raw_source_codes %in%
            raw_clean$original_study_code
        )
      
      
      numeric_fields <- c(
        "serum_creatinine_umol_l",
        "crp_mg_l",
        "lactate_mmol_l",
        "albumin_g_l",
        "procalcitonin_ng_ml",
        "alt_u_l",
        "ast_u_l",
        "wbc_10e9_l",
        "platelets_10e9_l"
      )
      
      
      numeric_mismatches <- 0L
      numeric_cells_checked <- 0L
      
      
      for (
        field_name in numeric_fields
      ) {
        
        deid_column_name <- paste0(
          field_name,
          "_deid"
        )
        
        
        raw_column_name <- paste0(
          field_name,
          "_raw"
        )
        
        
        deid_values <-
          raw_compare[[deid_column_name]]
        
        
        raw_values <-
          raw_compare[[raw_column_name]]
        
        
        comparison <- compare_numeric(
          deid_values,
          raw_values
        )
        
        
        numeric_mismatches <-
          numeric_mismatches +
          sum(
            !comparison
          )
        
        
        numeric_cells_checked <-
          numeric_cells_checked +
          length(
            comparison
          )
      }
      
      
      outcome_match <- compare_character(
        raw_compare$outcome_status,
        raw_compare$raw_outcome_status
      )
      
      
      ventilation_match <- compare_character(
        raw_compare$mechanical_ventilation_status,
        raw_compare$raw_ventilation_status
      )
      
      
      culture_match <- compare_character(
        raw_compare$culture_status,
        raw_compare$raw_culture_status
      )
      
      
      categorical_checks <- c(
        outcome_match,
        ventilation_match,
        culture_match
      )
      
      
      categorical_mismatches <-
        sum(
          !categorical_checks
        )
      
      
      raw_audit_summary$value <- c(
        "Yes",
        nrow(
          raw_clean
        ),
        ifelse(
          mapping_complete,
          "Yes",
          "No"
        ),
        numeric_cells_checked,
        numeric_mismatches,
        length(
          categorical_checks
        ),
        categorical_mismatches
      )
      
      
      if (
        !mapping_complete
      ) {
        
        stop(
          "Raw-source study-code mapping audit failed."
        )
      }
      
      
      if (
        numeric_mismatches !=
        0
      ) {
        
        stop(
          paste0(
            "Raw-source numeric audit detected ",
            numeric_mismatches,
            " mismatch(es)."
          )
        )
      }
      
      
      if (
        categorical_mismatches !=
        0
      ) {
        
        stop(
          paste0(
            "Raw-source categorical audit detected ",
            categorical_mismatches,
            " mismatch(es)."
          )
        )
      }
      
      
      cat(
        "\nRaw-source integrity audit passed successfully.\n"
      )
      
    } else {
      
      cat(
        "\nRaw workbook found, but expected raw columns are missing:\n"
      )
      
      
      print(
        missing_raw_columns
      )
      
      
      cat(
        "Optional raw-source integrity audit skipped.\n"
      )
    }
    
  } else {
    
    cat(
      "\nSheet 'Пациенты Сепсис' was not found in raw workbook.\n"
    )
    
    
    cat(
      "Optional raw-source integrity audit skipped.\n"
    )
  }
}


# =============================================================================
# 16. MISSINGNESS AUDIT
# =============================================================================

clinical_variables <- c(
  "age_years",
  "serum_creatinine_umol_l",
  "crp_mg_l",
  "lactate_mmol_l",
  "albumin_g_l",
  "procalcitonin_ng_ml",
  "alt_u_l",
  "ast_u_l",
  "wbc_10e9_l",
  "platelets_10e9_l"
)


missingness_sepsis <- data.frame(
  
  variable =
    clinical_variables,
  
  n_total =
    nrow(
      sepsis
    ),
  
  n_available =
    vapply(
      clinical_variables,
      function(variable_name) {
        
        values <-
          sepsis[[variable_name]]
        
        
        values <- suppressWarnings(
          as.numeric(
            values
          )
        )
        
        
        sum(
          is.finite(
            values
          )
        )
      },
      numeric(1)
    ),
  
  stringsAsFactors = FALSE
)


missingness_sepsis$n_missing <-
  missingness_sepsis$n_total -
  missingness_sepsis$n_available


missingness_sepsis$percent_missing <-
  100 *
  missingness_sepsis$n_missing /
  missingness_sepsis$n_total


# =============================================================================
# 17. AGE AND SEX TESTS
# =============================================================================

age_sepsis <- sepsis$age_years[
  is.finite(
    sepsis$age_years
  )
]


age_controls <- controls$age_years[
  is.finite(
    controls$age_years
  )
]


age_test <- stats::wilcox.test(
  age_sepsis,
  age_controls,
  alternative = "two.sided",
  exact = FALSE,
  correct = TRUE
)


age_p <- as.numeric(
  age_test$p.value
)


sex_table <- matrix(
  
  c(
    sum(
      sepsis$sex ==
        "Male",
      na.rm = TRUE
    ),
    
    sum(
      sepsis$sex ==
        "Female",
      na.rm = TRUE
    ),
    
    sum(
      controls$sex ==
        "Male",
      na.rm = TRUE
    ),
    
    sum(
      controls$sex ==
        "Female",
      na.rm = TRUE
    )
  ),
  
  nrow = 2,
  
  byrow = TRUE,
  
  dimnames = list(
    
    Cohort = c(
      "Sepsis",
      "Healthy control"
    ),
    
    Sex = c(
      "Male",
      "Female"
    )
  )
)


sex_test <- stats::fisher.test(
  sex_table,
  alternative = "two.sided"
)


sex_p <- as.numeric(
  sex_test$p.value
)


# =============================================================================
# 18. DEMOGRAPHIC HARD AUDIT
# =============================================================================
#
# These anchors reproduce the previously audited blood-cohort demographics.
# Tolerance is intentionally modest to avoid machine/R-version rounding issues.
#
# =============================================================================

expected_age_p <-
  0.0074308


expected_sex_p <-
  0.0101468


demographic_audit <- data.frame(
  
  metric = c(
    "Sepsis n",
    "Control n",
    "Sepsis age available n",
    "Control age available n",
    "Sepsis male n",
    "Sepsis female n",
    "Control male n",
    "Control female n",
    "Age P",
    "Sex P"
  ),
  
  observed = c(
    nrow(
      sepsis
    ),
    nrow(
      controls
    ),
    length(
      age_sepsis
    ),
    length(
      age_controls
    ),
    
    sum(
      sepsis$sex ==
        "Male"
    ),
    
    sum(
      sepsis$sex ==
        "Female"
    ),
    
    sum(
      controls$sex ==
        "Male"
    ),
    
    sum(
      controls$sex ==
        "Female"
    ),
    
    age_p,
    sex_p
  ),
  
  expected = c(
    35,
    10,
    34,
    10,
    24,
    11,
    2,
    8,
    expected_age_p,
    expected_sex_p
  ),
  
  stringsAsFactors = FALSE
)


demographic_audit$difference <-
  demographic_audit$observed -
  demographic_audit$expected


count_rows <- seq_len(
  8
)


if (
  any(
    demographic_audit$difference[count_rows] !=
    0
  )
) {
  
  print(
    demographic_audit
  )
  
  
  stop(
    "Demographic count audit failed."
  )
}


if (
  abs(
    age_p -
    expected_age_p
  ) >
  1e-6
) {
  
  stop(
    "Age P-value audit failed."
  )
}


if (
  abs(
    sex_p -
    expected_sex_p
  ) >
  1e-6
) {
  
  stop(
    "Sex P-value audit failed."
  )
}


cat(
  "\nDemographic hard audit passed successfully.\n"
)


# =============================================================================
# 19. CLINICAL SUMMARY
# =============================================================================

clinical_summary_labels <- c(
  "Age, years",
  "Serum creatinine, \u00b5mol/L",
  "C-reactive protein, mg/L",
  "Lactate, mmol/L",
  "Albumin, g/L",
  "Procalcitonin, ng/mL",
  "ALT, U/L",
  "AST, U/L",
  "WBC, \u00d710^9/L",
  "Platelets, \u00d710^9/L"
)


clinical_summary <- data.frame(
  
  variable =
    clinical_summary_labels,
  
  source_variable =
    clinical_variables,
  
  n_available =
    NA_integer_,
  
  median =
    NA_real_,
  
  q1 =
    NA_real_,
  
  q3 =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (
  i in seq_along(
    clinical_variables
  )
) {
  
  variable_name <-
    clinical_variables[i]
  
  
  values <-
    sepsis[[variable_name]]
  
  
  values <- suppressWarnings(
    as.numeric(
      values
    )
  )
  
  
  finite_values <- values[
    is.finite(
      values
    )
  ]
  
  
  clinical_summary$n_available[i] <-
    length(
      finite_values
    )
  
  
  clinical_summary$median[i] <-
    stats::median(
      finite_values
    )
  
  
  clinical_summary$q1[i] <-
    as.numeric(
      stats::quantile(
        finite_values,
        0.25,
        type = 7
      )
    )
  
  
  clinical_summary$q3[i] <-
    as.numeric(
      stats::quantile(
        finite_values,
        0.75,
        type = 7
      )
    )
}


# =============================================================================
# 20. CATEGORICAL CLINICAL CHARACTERISTICS
# =============================================================================

n_sepsis <- nrow(
  sepsis
)


n_controls <- nrow(
  controls
)


n_death <- sum(
  sepsis$outcome_status ==
    "Death",
  na.rm = TRUE
)


n_vent <- sum(
  sepsis$mechanical_ventilation_status ==
    "Yes",
  na.rm = TRUE
)


n_culture_positive <- sum(
  sepsis$culture_status ==
    "Positive",
  na.rm = TRUE
)


if (
  n_death != 18
) {
  
  stop(
    paste0(
      "Expected 18 deaths; observed ",
      n_death
    )
  )
}


if (
  n_vent != 18
) {
  
  stop(
    paste0(
      "Expected 18 mechanically ventilated patients; observed ",
      n_vent
    )
  )
}


if (
  n_culture_positive != 13
) {
  
  stop(
    paste0(
      "Expected 13 culture-positive patients; observed ",
      n_culture_positive
    )
  )
}


# =============================================================================
# 21. BUILD MAIN TABLE 1
# =============================================================================

table1 <- data.frame(
  
  Characteristic = c(
    
    "Demographic characteristics",
    
    "Age, years",
    
    "Male sex, n (%)",
    
    "Clinical characteristics of patients with sepsis",
    
    "In-hospital death, n (%)",
    
    "Invasive mechanical ventilation, n (%)",
    
    "Positive microbiological culture, n (%)",
    
    "Serum creatinine, \u00b5mol/L",
    
    "C-reactive protein, mg/L",
    
    "Lactate, mmol/L",
    
    "Albumin, g/L",
    
    "Procalcitonin, ng/mL",
    
    "ALT, U/L",
    
    "AST, U/L",
    
    "WBC, \u00d710^9/L",
    
    "Platelets, \u00d710^9/L"
  ),
  
  Sepsis = c(
    
    "",
    
    median_iqr_string(
      sepsis$age_years,
      digits = 1
    ),
    
    n_percent_string(
      sum(
        sepsis$sex ==
          "Male"
      ),
      n_sepsis
    ),
    
    "",
    
    n_percent_string(
      n_death,
      n_sepsis
    ),
    
    n_percent_string(
      n_vent,
      n_sepsis
    ),
    
    n_percent_string(
      n_culture_positive,
      n_sepsis
    ),
    
    median_iqr_string(
      sepsis$serum_creatinine_umol_l,
      digits = 1
    ),
    
    median_iqr_string(
      sepsis$crp_mg_l,
      digits = 0
    ),
    
    median_iqr_string(
      sepsis$lactate_mmol_l,
      digits = 2
    ),
    
    median_iqr_string(
      sepsis$albumin_g_l,
      digits = 1
    ),
    
    median_iqr_string(
      sepsis$procalcitonin_ng_ml,
      digits = 1
    ),
    
    median_iqr_string(
      sepsis$alt_u_l,
      digits = 1
    ),
    
    median_iqr_string(
      sepsis$ast_u_l,
      digits = 1
    ),
    
    median_iqr_string(
      sepsis$wbc_10e9_l,
      digits = 1
    ),
    
    median_iqr_string(
      sepsis$platelets_10e9_l,
      digits = 1
    )
  ),
  
  `Healthy controls` = c(
    
    "",
    
    median_iqr_string(
      controls$age_years,
      digits = 1
    ),
    
    n_percent_string(
      sum(
        controls$sex ==
          "Male"
      ),
      n_controls
    ),
    
    "",
    
    rep(
      "\u2014",
      12
    )
  ),
  
  `P value` = c(
    
    "",
    
    format_p_table(
      age_p
    ),
    
    format_p_table(
      sex_p
    ),
    
    "",
    
    rep(
      "\u2014",
      12
    )
  ),
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
)


table1_export <-
  table1


names(
  table1_export
) <- c(
  "Characteristic",
  paste0(
    "Sepsis (n=",
    n_sepsis,
    ")"
  ),
  paste0(
    "Healthy controls (n=",
    n_controls,
    ")"
  ),
  "P value"
)


table1_title <- paste0(
  "Table 1. Demographic and clinical characteristics of the ",
  "discovery blood transcriptomic cohort"
)


table1_footnotes <- c(
  
  paste0(
    "Data are presented as median [interquartile range] for continuous ",
    "variables and n (%) for categorical variables."
  ),
  
  paste0(
    "Age was available for 34 of 35 participants with sepsis and for ",
    "all 10 healthy controls."
  ),
  
  paste0(
    "Procalcitonin was available for 30 of 35 participants with sepsis; ",
    "all other sepsis clinical variables shown were available for all ",
    "35 participants."
  ),
  
  paste0(
    "P values are reported only for variables applicable to both cohorts. ",
    "Age was compared using a two-sided Wilcoxon rank-sum test and sex ",
    "using Fisher's exact test."
  ),
  
  paste0(
    "Clinical variables specific to acute illness were not applicable ",
    "to healthy controls."
  ),
  
  paste0(
    "SOFA is not shown because SOFA values were unavailable in the ",
    "discovery-cohort clinical source used for this table."
  ),
  
  paste0(
    "ALT, alanine aminotransferase; AST, aspartate aminotransferase; ",
    "CRP, C-reactive protein; WBC, white blood cell count."
  )
)


# =============================================================================
# 22. BUILD SUPPLEMENTARY TABLE S1 — PARTICIPANTS
# =============================================================================

tableS1_participant <- participant_blood %>%
  
  dplyr::select(
    participant_id,
    cohort,
    sex,
    age_years,
    blood_rnaseq_sample_id,
    outcome_status,
    mechanical_ventilation_status,
    culture_status,
    serum_creatinine_umol_l,
    creatinine_group_110_umol_l,
    crp_mg_l,
    lactate_mmol_l,
    albumin_g_l,
    procalcitonin_ng_ml,
    alt_u_l,
    ast_u_l,
    wbc_10e9_l,
    platelets_10e9_l
  ) %>%
  
  dplyr::arrange(
    factor(
      cohort,
      levels = c(
        "Sepsis",
        "Healthy control"
      )
    ),
    participant_id
  )


# =============================================================================
# 23. BUILD SUPPLEMENTARY TABLE S1 — BLOOD SAMPLES
# =============================================================================

tableS1_samples <- blood_samples %>%
  
  dplyr::arrange(
    factor(
      sample_group,
      levels = c(
        "BP",
        "BC"
      )
    ),
    participant_id
  )


# =============================================================================
# 24. DATA DICTIONARY
# =============================================================================

data_dictionary <- data.frame(
  
  variable_name = c(
    "participant_id",
    "cohort",
    "sex",
    "age_years",
    "blood_rnaseq_sample_id",
    "outcome_status",
    "mechanical_ventilation_status",
    "culture_status",
    "serum_creatinine_umol_l",
    "creatinine_group_110_umol_l",
    "crp_mg_l",
    "lactate_mmol_l",
    "albumin_g_l",
    "procalcitonin_ng_ml",
    "alt_u_l",
    "ast_u_l",
    "wbc_10e9_l",
    "platelets_10e9_l",
    "sample_id",
    "biofluid",
    "sample_group",
    "sequencing_batch"
  ),
  
  description = c(
    "De-identified participant identifier.",
    "Clinical cohort.",
    "Biological sex.",
    "Age in completed years.",
    "Blood RNA-seq sample identifier.",
    "Hospital outcome.",
    "Invasive mechanical ventilation status.",
    "Microbiological culture status.",
    "Serum creatinine concentration.",
    "Exploratory creatinine category.",
    "C-reactive protein concentration.",
    "Blood lactate concentration.",
    "Serum albumin concentration.",
    "Procalcitonin concentration.",
    "Alanine aminotransferase activity.",
    "Aspartate aminotransferase activity.",
    "White blood cell count.",
    "Platelet count.",
    "Blood RNA-seq sample identifier.",
    "RNA source material.",
    "Analytical blood sample group.",
    "Sequencing batch/chip variable."
  ),
  
  coding_or_units = c(
    "De-identified code",
    "Sepsis; Healthy control",
    "Male; Female",
    "years",
    "BP# or BC#",
    "Death; Discharged; blank = not applicable",
    "Yes; No; blank = not applicable",
    "Positive; Negative; blank = not applicable",
    "\u00b5mol/L",
    ">110; \u2264110 \u00b5mol/L",
    "mg/L",
    "mmol/L",
    "g/L",
    "ng/mL",
    "U/L",
    "U/L",
    "\u00d710^9/L",
    "\u00d710^9/L",
    "BP# or BC#",
    "Whole blood",
    "BP; BC",
    "chip/batch code"
  ),
  
  worksheet = c(
    rep(
      "Participant_Metadata",
      18
    ),
    rep(
      "Blood_Sample_Metadata",
      4
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 25. TABLE S1 README
# =============================================================================

tableS1_readme <- data.frame(
  
  Item = c(
    "Title",
    "Scope",
    "Participants",
    "Blood samples",
    "Data minimization",
    "Age",
    "Clinical variables",
    "Urine data",
    "Direct identifiers",
    "Missing values"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S1. De-identified participant- and blood ",
      "RNA-seq sample-level metadata for the discovery cohort"
    ),
    
    "Blood-only manuscript metadata.",
    
    "35 participants with sepsis and 10 healthy controls.",
    
    "35 BP samples and 10 BC samples.",
    
    paste0(
      "Only variables needed for the blood manuscript and reported ",
      "analyses are retained."
    ),
    
    paste0(
      "Age was available for 34/35 participants with sepsis and ",
      "10/10 healthy controls."
    ),
    
    paste0(
      "Acute clinical variables are applicable to the sepsis cohort. ",
      "Procalcitonin was available for 30/35 sepsis participants."
    ),
    
    paste0(
      "Urine sample IDs, paired urine metadata, and urine-specific ",
      "variables are intentionally excluded."
    ),
    
    paste0(
      "Names, raw patient IDs, dates, birth years, original study codes, ",
      "hospital identifiers, pathogen names, and free-text infection sites ",
      "are not included."
    ),
    
    "Blank cells indicate unavailable or non-applicable values."
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 26. SOFA AUDIT
# =============================================================================

sofa_values <- suppressWarnings(
  as.numeric(
    participant_extended$sofa_score[
      participant_extended$cohort ==
        "Sepsis"
    ]
  )
)


sofa_available_n <- sum(
  is.finite(
    sofa_values
  )
)


if (
  sofa_available_n != 0
) {
  
  stop(
    paste0(
      "Expected 0 discovery SOFA values; observed ",
      sofa_available_n,
      ". Re-audit source before Table 1 is generated."
    )
  )
}


# =============================================================================
# 27. DATA-MINIMIZATION AUDIT
# =============================================================================

forbidden_columns <- c(
  "original_study_code",
  "ID",
  "#P",
  "name",
  "birth_year",
  "detected_pathogen",
  "infection_focus_category",
  "urine_rnaseq_sample_id",
  "paired_sample_id"
)


exported_columns <- unique(
  c(
    names(
      tableS1_participant
    ),
    names(
      tableS1_samples
    )
  )
)


forbidden_present <- intersect(
  forbidden_columns,
  exported_columns
)


if (
  length(
    forbidden_present
  ) >
  0
) {
  
  stop(
    paste0(
      "Data-minimization audit failed. Forbidden columns present: ",
      paste(
        forbidden_present,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 28. VALIDATION AUDIT
# =============================================================================

validation_audit <- data.frame(
  
  metric = c(
    "Total participants",
    "Sepsis participants",
    "Healthy controls",
    "Blood RNA-seq samples",
    "BP samples",
    "BC samples",
    "Unique participant IDs",
    "Unique blood sample IDs",
    "Sepsis age available n",
    "Control age available n",
    "Sepsis male",
    "Sepsis female",
    "Control male",
    "Control female",
    "Age Wilcoxon P",
    "Sex Fisher exact P",
    "Deaths",
    "Mechanical ventilation yes",
    "Positive culture",
    "Procalcitonin available n",
    "SOFA available n",
    "REVISED vs extended shared fields all match",
    "Forbidden export columns",
    "Urine sample records exported"
  ),
  
  observed = c(
    nrow(
      participant_blood
    ),
    nrow(
      sepsis
    ),
    nrow(
      controls
    ),
    nrow(
      blood_samples
    ),
    n_BP,
    n_BC,
    
    dplyr::n_distinct(
      participant_blood$participant_id
    ),
    
    dplyr::n_distinct(
      blood_samples$sample_id
    ),
    
    length(
      age_sepsis
    ),
    
    length(
      age_controls
    ),
    
    sum(
      sepsis$sex ==
        "Male"
    ),
    
    sum(
      sepsis$sex ==
        "Female"
    ),
    
    sum(
      controls$sex ==
        "Male"
    ),
    
    sum(
      controls$sex ==
        "Female"
    ),
    
    age_p,
    sex_p,
    n_death,
    n_vent,
    n_culture_positive,
    
    sum(
      is.finite(
        sepsis$procalcitonin_ng_ml
      )
    ),
    
    sofa_available_n,
    
    all(
      shared_field_audit$all_match
    ),
    
    length(
      forbidden_present
    ),
    
    0
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 29. OUTPUT PATHS
# =============================================================================

table1_xlsx <- file.path(
  tables_dir,
  "151_Table1_discovery_blood_cohort.xlsx"
)


table1_csv <- file.path(
  tables_dir,
  "151_Table1_discovery_blood_cohort.csv"
)


tableS1_xlsx <- file.path(
  tables_dir,
  "151_TableS1_deidentified_blood_metadata.xlsx"
)


tableS1_participant_csv <- file.path(
  tables_dir,
  "151_TableS1_participant_metadata.csv"
)


tableS1_sample_csv <- file.path(
  tables_dir,
  "151_TableS1_blood_sample_metadata.csv"
)


# =============================================================================
# 30. WRITE MAIN TABLE 1 — OPENXLSX
# =============================================================================

wb1 <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wb1,
  "Table1"
)


openxlsx::addWorksheet(
  wb1,
  "Clinical_summary_source"
)


openxlsx::addWorksheet(
  wb1,
  "Missingness"
)


openxlsx::addWorksheet(
  wb1,
  "Validation_audit"
)


openxlsx::addWorksheet(
  wb1,
  "Source_cross_audit"
)


openxlsx::addWorksheet(
  wb1,
  "Raw_source_audit"
)


title_style <- openxlsx::createStyle(
  fontSize = 14,
  textDecoration = "bold",
  halign = "left"
)


header_style <- openxlsx::createStyle(
  fontSize = 11,
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  halign = "center",
  valign = "center"
)


section_style <- openxlsx::createStyle(
  fontSize = 10.5,
  textDecoration = "bold",
  fgFill = "#EEF3F7",
  border = "Bottom",
  borderStyle = "thin"
)


body_left_style <- openxlsx::createStyle(
  fontSize = 10,
  halign = "left",
  valign = "center"
)


body_center_style <- openxlsx::createStyle(
  fontSize = 10,
  halign = "center",
  valign = "center"
)


footnote_style <- openxlsx::createStyle(
  fontSize = 9,
  fontColour = "#404040",
  wrapText = TRUE
)


openxlsx::writeData(
  wb1,
  "Table1",
  table1_title,
  startRow = 1,
  startCol = 1
)


openxlsx::mergeCells(
  wb1,
  "Table1",
  cols = 1:4,
  rows = 1
)


openxlsx::addStyle(
  wb1,
  "Table1",
  title_style,
  rows = 1,
  cols = 1:4,
  gridExpand = TRUE
)


openxlsx::writeData(
  wb1,
  "Table1",
  table1_export,
  startRow = 3,
  startCol = 1,
  headerStyle = header_style
)


first_data_row <- 4


last_data_row <-
  first_data_row +
  nrow(
    table1_export
  ) -
  1


section_indices <- which(
  table1$Characteristic %in%
    c(
      "Demographic characteristics",
      "Clinical characteristics of patients with sepsis"
    )
)


section_rows <-
  first_data_row +
  section_indices -
  1


for (
  rr in section_rows
) {
  
  openxlsx::addStyle(
    wb1,
    "Table1",
    section_style,
    rows = rr,
    cols = 1:4,
    gridExpand = TRUE,
    stack = TRUE
  )
}


openxlsx::addStyle(
  wb1,
  "Table1",
  body_left_style,
  rows = first_data_row:last_data_row,
  cols = 1,
  gridExpand = TRUE,
  stack = TRUE
)


openxlsx::addStyle(
  wb1,
  "Table1",
  body_center_style,
  rows = first_data_row:last_data_row,
  cols = 2:4,
  gridExpand = TRUE,
  stack = TRUE
)


footnote_start_row <-
  last_data_row +
  2


for (
  i in seq_along(
    table1_footnotes
  )
) {
  
  rr <-
    footnote_start_row +
    i -
    1
  
  
  openxlsx::writeData(
    wb1,
    "Table1",
    paste0(
      i,
      ". ",
      table1_footnotes[i]
    ),
    startRow = rr,
    startCol = 1
  )
  
  
  openxlsx::mergeCells(
    wb1,
    "Table1",
    cols = 1:4,
    rows = rr
  )
  
  
  openxlsx::addStyle(
    wb1,
    "Table1",
    footnote_style,
    rows = rr,
    cols = 1:4,
    gridExpand = TRUE
  )
}


openxlsx::setColWidths(
  wb1,
  "Table1",
  cols = 1,
  widths = 45
)


openxlsx::setColWidths(
  wb1,
  "Table1",
  cols = 2:3,
  widths = 25
)


openxlsx::setColWidths(
  wb1,
  "Table1",
  cols = 4,
  widths = 14
)


openxlsx::freezePane(
  wb1,
  "Table1",
  firstActiveRow = 4
)


openxlsx::writeData(
  wb1,
  "Clinical_summary_source",
  clinical_summary,
  withFilter = TRUE
)


openxlsx::writeData(
  wb1,
  "Missingness",
  missingness_sepsis,
  withFilter = TRUE
)


openxlsx::writeData(
  wb1,
  "Validation_audit",
  validation_audit,
  withFilter = TRUE
)


openxlsx::writeData(
  wb1,
  "Source_cross_audit",
  shared_field_audit,
  withFilter = TRUE
)


openxlsx::writeData(
  wb1,
  "Raw_source_audit",
  raw_audit_summary,
  withFilter = TRUE
)


for (
  sheet_name in c(
    "Clinical_summary_source",
    "Missingness",
    "Validation_audit",
    "Source_cross_audit",
    "Raw_source_audit"
  )
) {
  
  openxlsx::setColWidths(
    wb1,
    sheet_name,
    cols = 1:20,
    widths = "auto"
  )
  
  
  openxlsx::freezePane(
    wb1,
    sheet_name,
    firstActiveRow = 2
  )
}


openxlsx::saveWorkbook(
  wb1,
  table1_xlsx,
  overwrite = TRUE
)


write.csv(
  table1_export,
  table1_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 31. WRITE TABLE S1 — OPENXLSX
# =============================================================================

wbS1 <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wbS1,
  "S1_ReadMe"
)


openxlsx::addWorksheet(
  wbS1,
  "Participant_Metadata"
)


openxlsx::addWorksheet(
  wbS1,
  "Blood_Sample_Metadata"
)


openxlsx::addWorksheet(
  wbS1,
  "Data_Dictionary"
)


openxlsx::addWorksheet(
  wbS1,
  "Validation_Audit"
)


openxlsx::addWorksheet(
  wbS1,
  "Source_Cross_Audit"
)


openxlsx::writeData(
  wbS1,
  "S1_ReadMe",
  tableS1_readme
)


openxlsx::writeData(
  wbS1,
  "Participant_Metadata",
  tableS1_participant,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS1,
  "Blood_Sample_Metadata",
  tableS1_samples,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS1,
  "Data_Dictionary",
  data_dictionary,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS1,
  "Validation_Audit",
  validation_audit,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS1,
  "Source_Cross_Audit",
  shared_field_audit,
  withFilter = TRUE
)


supp_header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  valign = "center",
  wrapText = TRUE
)


openxlsx::addStyle(
  wbS1,
  "S1_ReadMe",
  supp_header_style,
  rows = 1,
  cols = 1:ncol(
    tableS1_readme
  ),
  gridExpand = TRUE
)


openxlsx::addStyle(
  wbS1,
  "Participant_Metadata",
  supp_header_style,
  rows = 1,
  cols = 1:ncol(
    tableS1_participant
  ),
  gridExpand = TRUE
)


openxlsx::addStyle(
  wbS1,
  "Blood_Sample_Metadata",
  supp_header_style,
  rows = 1,
  cols = 1:ncol(
    tableS1_samples
  ),
  gridExpand = TRUE
)


openxlsx::addStyle(
  wbS1,
  "Data_Dictionary",
  supp_header_style,
  rows = 1,
  cols = 1:ncol(
    data_dictionary
  ),
  gridExpand = TRUE
)


openxlsx::addStyle(
  wbS1,
  "Validation_Audit",
  supp_header_style,
  rows = 1,
  cols = 1:ncol(
    validation_audit
  ),
  gridExpand = TRUE
)


openxlsx::addStyle(
  wbS1,
  "Source_Cross_Audit",
  supp_header_style,
  rows = 1,
  cols = 1:ncol(
    shared_field_audit
  ),
  gridExpand = TRUE
)


for (
  sheet_name in c(
    "S1_ReadMe",
    "Participant_Metadata",
    "Blood_Sample_Metadata",
    "Data_Dictionary",
    "Validation_Audit",
    "Source_Cross_Audit"
  )
) {
  
  openxlsx::freezePane(
    wbS1,
    sheet_name,
    firstActiveRow = 2
  )
}


openxlsx::setColWidths(
  wbS1,
  "S1_ReadMe",
  cols = 1,
  widths = 28
)


openxlsx::setColWidths(
  wbS1,
  "S1_ReadMe",
  cols = 2,
  widths = 80
)


openxlsx::setColWidths(
  wbS1,
  "Participant_Metadata",
  cols = 1:ncol(
    tableS1_participant
  ),
  widths = "auto"
)


openxlsx::setColWidths(
  wbS1,
  "Blood_Sample_Metadata",
  cols = 1:ncol(
    tableS1_samples
  ),
  widths = "auto"
)


openxlsx::setColWidths(
  wbS1,
  "Data_Dictionary",
  cols = 1,
  widths = 30
)


openxlsx::setColWidths(
  wbS1,
  "Data_Dictionary",
  cols = 2,
  widths = 65
)


openxlsx::setColWidths(
  wbS1,
  "Data_Dictionary",
  cols = 3,
  widths = 32
)


openxlsx::setColWidths(
  wbS1,
  "Data_Dictionary",
  cols = 4,
  widths = 24
)


openxlsx::setColWidths(
  wbS1,
  "Validation_Audit",
  cols = 1,
  widths = 48
)


openxlsx::setColWidths(
  wbS1,
  "Validation_Audit",
  cols = 2,
  widths = 25
)


openxlsx::setColWidths(
  wbS1,
  "Source_Cross_Audit",
  cols = 1:5,
  widths = "auto"
)


openxlsx::saveWorkbook(
  wbS1,
  tableS1_xlsx,
  overwrite = TRUE
)


write.csv(
  tableS1_participant,
  tableS1_participant_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  tableS1_samples,
  tableS1_sample_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 32. TEXT OUTPUTS
# =============================================================================

writeLines(
  table1_title,
  file.path(
    text_dir,
    "151_Table1_title_EN.txt"
  )
)


writeLines(
  table1_footnotes,
  file.path(
    text_dir,
    "151_Table1_footnotes_EN.txt"
  )
)


results_31_text <- paste0(
  
  "The discovery blood transcriptomic cohort comprised 35 patients with ",
  "sepsis and 10 healthy controls. Age was available for 34 of 35 patients ",
  "with sepsis and for all controls. Patients with sepsis were older than ",
  "healthy controls (median ",
  median_iqr_string(
    sepsis$age_years,
    digits = 1
  ),
  " versus ",
  median_iqr_string(
    controls$age_years,
    digits = 1
  ),
  " years; P=",
  sprintf(
    "%.4f",
    age_p
  ),
  "), and the sex distribution differed between groups (male sex, ",
  n_percent_string(
    sum(
      sepsis$sex ==
        "Male"
    ),
    n_sepsis
  ),
  " versus ",
  n_percent_string(
    sum(
      controls$sex ==
        "Male"
    ),
    n_controls
  ),
  "; Fisher exact P=",
  sprintf(
    "%.4f",
    sex_p
  ),
  "). Among patients with sepsis, ",
  n_percent_string(
    n_death,
    n_sepsis
  ),
  " died during hospitalization, ",
  n_percent_string(
    n_vent,
    n_sepsis
  ),
  " received invasive mechanical ventilation, and ",
  n_percent_string(
    n_culture_positive,
    n_sepsis
  ),
  " had a positive microbiological culture (Table 1)."
)


writeLines(
  results_31_text,
  file.path(
    text_dir,
    "151_Results_3.1_cohort_description_EN.txt"
  )
)


# =============================================================================
# 33. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "151_sessionInfo.txt"
  )
)


# =============================================================================
# 34. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 151 completed successfully.\n")
cat("====================================================================\n\n")


cat("DISCOVERY BLOOD COHORT\n")
cat("----------------------\n")


cat(
  "Sepsis participants = ",
  n_sepsis,
  "\n",
  sep = ""
)


cat(
  "Healthy controls = ",
  n_controls,
  "\n",
  sep = ""
)


cat(
  "Whole-blood RNA-seq samples = ",
  nrow(
    blood_samples
  ),
  "\n",
  sep = ""
)


cat(
  "BP = ",
  n_BP,
  "\n",
  sep = ""
)


cat(
  "BC = ",
  n_BC,
  "\n",
  sep = ""
)


cat("\nDEMOGRAPHICS\n")
cat("------------\n")


cat(
  "Sepsis age available n = ",
  length(
    age_sepsis
  ),
  "\n",
  sep = ""
)


cat(
  "Sepsis age = ",
  median_iqr_string(
    sepsis$age_years,
    digits = 1
  ),
  "\n",
  sep = ""
)


cat(
  "Control age available n = ",
  length(
    age_controls
  ),
  "\n",
  sep = ""
)


cat(
  "Control age = ",
  median_iqr_string(
    controls$age_years,
    digits = 1
  ),
  "\n",
  sep = ""
)


cat(
  "Age Wilcoxon P = ",
  format(
    age_p,
    scientific = TRUE,
    digits = 8
  ),
  "\n",
  sep = ""
)


cat(
  "Sepsis sex: Male ",
  sum(
    sepsis$sex ==
      "Male"
  ),
  "; Female ",
  sum(
    sepsis$sex ==
      "Female"
  ),
  "\n",
  sep = ""
)


cat(
  "Control sex: Male ",
  sum(
    controls$sex ==
      "Male"
  ),
  "; Female ",
  sum(
    controls$sex ==
      "Female"
  ),
  "\n",
  sep = ""
)


cat(
  "Sex Fisher exact P = ",
  format(
    sex_p,
    scientific = TRUE,
    digits = 8
  ),
  "\n",
  sep = ""
)


cat("\nSEPSIS CLINICAL CHARACTERISTICS\n")
cat("-------------------------------\n")


cat(
  "Deaths = ",
  n_death,
  "/35\n",
  sep = ""
)


cat(
  "Mechanical ventilation = ",
  n_vent,
  "/35\n",
  sep = ""
)


cat(
  "Positive culture = ",
  n_culture_positive,
  "/35\n",
  sep = ""
)


cat(
  "Creatinine = ",
  median_iqr_string(
    sepsis$serum_creatinine_umol_l,
    digits = 1
  ),
  " umol/L\n",
  sep = ""
)


cat(
  "CRP = ",
  median_iqr_string(
    sepsis$crp_mg_l,
    digits = 0
  ),
  " mg/L\n",
  sep = ""
)


cat(
  "Lactate = ",
  median_iqr_string(
    sepsis$lactate_mmol_l,
    digits = 2
  ),
  " mmol/L\n",
  sep = ""
)


cat(
  "Albumin = ",
  median_iqr_string(
    sepsis$albumin_g_l,
    digits = 1
  ),
  " g/L\n",
  sep = ""
)


cat(
  "Procalcitonin available n = ",
  sum(
    is.finite(
      sepsis$procalcitonin_ng_ml
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Procalcitonin = ",
  median_iqr_string(
    sepsis$procalcitonin_ng_ml,
    digits = 1
  ),
  " ng/mL\n",
  sep = ""
)


cat(
  "ALT = ",
  median_iqr_string(
    sepsis$alt_u_l,
    digits = 1
  ),
  " U/L\n",
  sep = ""
)


cat(
  "AST = ",
  median_iqr_string(
    sepsis$ast_u_l,
    digits = 1
  ),
  " U/L\n",
  sep = ""
)


cat(
  "WBC = ",
  median_iqr_string(
    sepsis$wbc_10e9_l,
    digits = 1
  ),
  " x10^9/L\n",
  sep = ""
)


cat(
  "Platelets = ",
  median_iqr_string(
    sepsis$platelets_10e9_l,
    digits = 1
  ),
  " x10^9/L\n",
  sep = ""
)


cat("\nDEMOGRAPHIC AUDIT\n")
cat("-----------------\n")


print(
  demographic_audit,
  row.names = FALSE
)


cat("\nDE-IDENTIFIED SOURCE CROSS-AUDIT\n")
cat("--------------------------------\n")


print(
  shared_field_audit,
  row.names = FALSE
)


cat("\nRAW SOURCE AUDIT\n")
cat("----------------\n")


print(
  raw_audit_summary,
  row.names = FALSE
)


cat("\nMISSINGNESS\n")
cat("-----------\n")


print(
  missingness_sepsis,
  row.names = FALSE
)


cat("\nDATA-MINIMIZATION CHECK\n")
cat("-----------------------\n")


cat(
  "Forbidden/direct-identifier columns in Table S1 outputs = ",
  length(
    forbidden_present
  ),
  "\n",
  sep = ""
)


cat(
  "Urine sample records exported = 0\n"
)


cat(
  "SOFA values available in discovery source = ",
  sofa_available_n,
  "/35\n",
  sep = ""
)


cat("\nMAIN TABLE 1\n")
cat("------------\n")


print(
  table1_export,
  row.names = FALSE
)


cat("\nFILES CREATED\n")
cat("-------------\n")


cat(
  "Main Table 1:\n  ",
  normalizePath(
    table1_xlsx,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Supplementary Table S1:\n  ",
  normalizePath(
    tableS1_xlsx,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Table 1 CSV:\n  ",
  normalizePath(
    table1_csv,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Table S1 participant CSV:\n  ",
  normalizePath(
    tableS1_participant_csv,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Table S1 blood-sample CSV:\n  ",
  normalizePath(
    tableS1_sample_csv,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nTABLE PLACEMENT\n")
cat("---------------\n")


cat(
  "Main Table 1 -> Results Section 3.1\n"
)


cat(
  "Supplementary Table S1 -> Supplementary Information\n"
)


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")


cat(
  "- Input workbooks read with readxl.\n"
)


cat(
  "- Output workbooks written with openxlsx.\n"
)


cat(
  "- Blood-only outputs.\n"
)


cat(
  "- No urine metadata exported.\n"
)


cat(
  "- No direct identifiers exported.\n"
)


cat(
  "- No missing values imputed.\n"
)


cat(
  "- SOFA unavailable in discovery source and therefore not reported.\n"
)


cat(
  "- Age comparison: two-sided Wilcoxon rank-sum test.\n"
)


cat(
  "- Sex comparison: Fisher exact test.\n"
)


cat(
  "- Acute-care source cross-audit restricted to sepsis participants.\n"
)


cat("\nDone.\n")
