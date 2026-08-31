################################################################################
# Script 159
# FINAL v6
#
# Supplementary Table S7
#
# Complete exploratory clinical-association landscape
#
# Project:
#   Sepsis_DESeq2
#
# PURPOSE
# -------
#
# Build Supplementary Table S7 from frozen Script 136b outputs.
#
# Canonical source:
#
#   results/blood_endotypes_biomarkers/
#   136b_demographic_sensitivity/tables/
#   136b_all_clinical_tests_updated_FDR.csv
#
# Duplicate provenance source:
#
#   results/blood_endotypes_biomarkers/
#   136b_demographic_sensitivity/tables/
#   136b_demographic_sensitivity.xlsx
#
#   sheet:
#   04_updated_all_tests
#
#
# IMPORTANT
# ---------
#
# This script DOES NOT:
#
#   - rerun clinical association analyses
#   - recalculate P values
#   - recalculate BH/FDR
#   - add tests
#   - remove tests
#   - redefine the multiplicity family
#
#
# SOURCE-EQUIVALENCE RULES
# ------------------------
#
# Numeric representations such as:
#
#   35
#   35.0
#   3.5e+01
#
# are equivalent.
#
# Textual missing values represented as:
#
#   NA
#   ""
#   whitespace-only strings
#
# are equivalent.
#
# Any disagreement between two NON-MISSING textual values remains
# a hard provenance error.
#
#
# FROZEN SCHEMA
# -------------
#
# framework
# clinical_variable
# clinical_label
# test_family
# test
# n
# statistic
# effect
# effect_name
# p_value
# group_summary
# note
# test_id
# BH_global
# BH_within_framework
# BH_within_test_family
# significance_global
#
#
# FROZEN ANALYTIC ANCHORS
# -----------------------
#
# Complete exploratory family:
#
#   60 evaluable tests
#
#
# Global BH significant:
#
# 1. Primary five-gene score vs CRP
#
#      Spearman rho ~ 0.574
#      P             ~ 3.09e-4
#      global BH q   ~ 0.0185
#
#
# 2. SRSq vs CRP
#
#      Spearman rho ~ 0.526
#      P             ~ 0.00117
#      global BH q   ~ 0.0352
#
#
# Exactly two associations should remain significant after
# global BH correction across all 60 evaluable tests.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 159 FINAL v6\n")
cat("Supplementary Table S7\n")
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
      paste(missing_packages, collapse = ", ")
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
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "159_TableS7_complete_clinical_associations"
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


for (one_dir in c(
  output_dir,
  tables_dir,
  audit_dir,
  text_dir
)) {
  
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
# 4. HELPER FUNCTIONS
# =============================================================================

safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


is_missing_text <- function(x) {
  
  x_chr <- as.character(x)
  
  is.na(x) |
    is.na(x_chr) |
    trimws(x_chr) == ""
}


normalize_text_for_comparison <- function(x) {
  
  x_chr <- as.character(x)
  
  missing_flag <- is_missing_text(x)
  
  x_chr <- stringr::str_squish(x_chr)
  
  x_chr[missing_flag] <- NA_character_
  
  x_chr
}


is_numeric_like_pair <- function(x, y) {
  
  xy <- c(
    as.character(x),
    as.character(y)
  )
  
  missing_flag <-
    is.na(xy) |
    trimws(xy) == ""
  
  observed <- xy[!missing_flag]
  
  if (length(observed) == 0) {
    return(FALSE)
  }
  
  parsed <- suppressWarnings(
    as.numeric(observed)
  )
  
  all(is.finite(parsed))
}


compare_column_pair <- function(x, y) {
  
  if (length(x) != length(y)) {
    stop("Column lengths differ during source-equivalence audit.")
  }
  
  
  numeric_like <- is_numeric_like_pair(x, y)
  
  n_values <- length(x)
  
  equal <- rep(FALSE, n_values)
  
  
  if (numeric_like) {
    
    x_num <- safe_numeric(x)
    y_num <- safe_numeric(y)
    
    x_missing <- is.na(x_num)
    y_missing <- is.na(y_num)
    
    both_missing <- x_missing & y_missing
    
    equal[both_missing] <- TRUE
    
    comparable <- !x_missing & !y_missing
    
    
    if (any(comparable)) {
      
      tolerance <-
        1e-12 +
        1e-10 *
        pmax(
          1,
          abs(x_num[comparable]),
          abs(y_num[comparable])
        )
      
      
      equal[comparable] <-
        abs(
          x_num[comparable] -
            y_num[comparable]
        ) <= tolerance
    }
    
    
    x_canonical <- ifelse(
      is.na(x_num),
      NA_character_,
      sprintf("%.15g", x_num)
    )
    
    
    y_canonical <- ifelse(
      is.na(y_num),
      NA_character_,
      sprintf("%.15g", y_num)
    )
    
    
    comparison_mode <- "numeric"
    
    
  } else {
    
    x_text <- normalize_text_for_comparison(x)
    y_text <- normalize_text_for_comparison(y)
    
    x_missing <- is.na(x_text)
    y_missing <- is.na(y_text)
    
    both_missing <- x_missing & y_missing
    
    equal[both_missing] <- TRUE
    
    comparable <- !x_missing & !y_missing
    
    
    if (any(comparable)) {
      
      equal[comparable] <-
        x_text[comparable] ==
        y_text[comparable]
    }
    
    
    x_canonical <- x_text
    y_canonical <- y_text
    
    comparison_mode <- "text"
  }
  
  
  list(
    equal = equal,
    comparison_mode = comparison_mode,
    x_canonical = x_canonical,
    y_canonical = y_canonical
  )
}


compare_tables_equivalent <- function(
    table_csv,
    table_xlsx
) {
  
  dimensions_match <- identical(
    dim(table_csv),
    dim(table_xlsx)
  )
  
  
  columns_match <- identical(
    names(table_csv),
    names(table_xlsx)
  )
  
  
  if (!dimensions_match || !columns_match) {
    
    return(
      list(
        dimensions_match = dimensions_match,
        columns_match = columns_match,
        rowwise_match = FALSE,
        reordered_match = FALSE,
        equivalent = FALSE,
        column_audit = data.frame(),
        differences = data.frame(),
        canonical_csv = NULL,
        canonical_xlsx = NULL
      )
    )
  }
  
  
  n_rows <- nrow(table_csv)
  n_cols <- ncol(table_csv)
  
  
  # Create zero-column data frames with the correct row count.
  
  canonical_csv <- data.frame(
    row.names = seq_len(n_rows)
  )
  
  
  canonical_xlsx <- data.frame(
    row.names = seq_len(n_rows)
  )
  
  
  column_audit_list <- vector(
    "list",
    n_cols
  )
  
  
  difference_list <- list()
  
  difference_counter <- 0L
  
  
  for (j in seq_len(n_cols)) {
    
    column_name <- names(table_csv)[j]
    
    
    comparison <- compare_column_pair(
      table_csv[[column_name]],
      table_xlsx[[column_name]]
    )
    
    
    canonical_csv[[column_name]] <-
      comparison$x_canonical
    
    
    canonical_xlsx[[column_name]] <-
      comparison$y_canonical
    
    
    n_different <- sum(
      !comparison$equal
    )
    
    
    column_audit_list[[j]] <- data.frame(
      
      column =
        column_name,
      
      comparison_mode =
        comparison$comparison_mode,
      
      n_rows =
        n_rows,
      
      n_different =
        n_different,
      
      match =
        n_different == 0,
      
      stringsAsFactors =
        FALSE
    )
    
    
    if (n_different > 0) {
      
      bad_rows <- which(
        !comparison$equal
      )
      
      
      for (one_row in bad_rows) {
        
        difference_counter <-
          difference_counter + 1L
        
        
        difference_list[[difference_counter]] <- data.frame(
          
          row =
            one_row,
          
          column =
            column_name,
          
          comparison_mode =
            comparison$comparison_mode,
          
          value_CSV =
            if (
              is.na(
                table_csv[[column_name]][one_row]
              )
            ) {
              "<NA>"
            } else {
              as.character(
                table_csv[[column_name]][one_row]
              )
            },
          
          value_XLSX =
            if (
              is.na(
                table_xlsx[[column_name]][one_row]
              )
            ) {
              "<NA>"
            } else {
              as.character(
                table_xlsx[[column_name]][one_row]
              )
            },
          
          stringsAsFactors =
            FALSE
        )
      }
    }
  }
  
  
  column_audit <- dplyr::bind_rows(
    column_audit_list
  )
  
  
  differences <- if (
    length(difference_list) >
    0
  ) {
    
    dplyr::bind_rows(
      difference_list
    )
    
  } else {
    
    data.frame(
      
      row =
        integer(),
      
      column =
        character(),
      
      comparison_mode =
        character(),
      
      value_CSV =
        character(),
      
      value_XLSX =
        character(),
      
      stringsAsFactors =
        FALSE
    )
  }
  
  
  rowwise_match <- all(
    column_audit$match
  )
  
  
  # ---------------------------------------------------------------------------
  # Secondary comparison allowing different row order
  # ---------------------------------------------------------------------------
  
  row_signature <- function(data) {
    
    apply(
      data,
      1,
      function(values) {
        
        values[
          is.na(values)
        ] <- "<MISSING>"
        
        paste(
          values,
          collapse = "|||"
        )
      }
    )
  }
  
  
  csv_signatures <- sort(
    row_signature(
      canonical_csv
    )
  )
  
  
  xlsx_signatures <- sort(
    row_signature(
      canonical_xlsx
    )
  )
  
  
  reordered_match <- identical(
    csv_signatures,
    xlsx_signatures
  )
  
  
  equivalent <-
    rowwise_match ||
    reordered_match
  
  
  list(
    dimensions_match = dimensions_match,
    columns_match = columns_match,
    rowwise_match = rowwise_match,
    reordered_match = reordered_match,
    equivalent = equivalent,
    column_audit = column_audit,
    differences = differences,
    canonical_csv = canonical_csv,
    canonical_xlsx = canonical_xlsx
  )
}


# =============================================================================
# 5. INPUT FILES
# =============================================================================

source_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity",
  "tables"
)


canonical_csv_file <- file.path(
  source_dir,
  "136b_all_clinical_tests_updated_FDR.csv"
)


workbook_file <- file.path(
  source_dir,
  "136b_demographic_sensitivity.xlsx"
)


workbook_sheet <- "04_updated_all_tests"


if (!file.exists(canonical_csv_file)) {
  
  stop(
    paste0(
      "Canonical CSV not found:\n",
      canonical_csv_file
    )
  )
}


if (!file.exists(workbook_file)) {
  
  stop(
    paste0(
      "Duplicate workbook not found:\n",
      workbook_file
    )
  )
}


if (
  !(workbook_sheet %in%
    readxl::excel_sheets(workbook_file))
) {
  
  stop(
    paste0(
      "Required workbook sheet not found: ",
      workbook_sheet
    )
  )
}


cat("\nCANONICAL FINAL SOURCES\n")
cat("-----------------------\n")


cat(
  "Standalone CSV:\n  ",
  normalizePath(
    canonical_csv_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Workbook:\n  ",
  normalizePath(
    workbook_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Workbook sheet: ",
  workbook_sheet,
  "\n",
  sep = ""
)


# =============================================================================
# 6. READ BOTH COPIES
# =============================================================================

clinical_csv <- read.csv(
  canonical_csv_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


clinical_xlsx <- readxl::read_excel(
  workbook_file,
  sheet = workbook_sheet
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat("\nSOURCE DIMENSIONS\n")
cat("-----------------\n")


cat(
  "CSV = ",
  nrow(clinical_csv),
  " x ",
  ncol(clinical_csv),
  "\n",
  sep = ""
)


cat(
  "XLSX sheet = ",
  nrow(clinical_xlsx),
  " x ",
  ncol(clinical_xlsx),
  "\n",
  sep = ""
)


# =============================================================================
# 7. FROZEN SCHEMA AUDIT
# =============================================================================

expected_columns <- c(
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
  "significance_global"
)


if (
  !identical(
    names(clinical_csv),
    expected_columns
  )
) {
  
  cat("\nObserved CSV columns:\n")
  
  print(
    names(clinical_csv)
  )
  
  
  cat("\nExpected frozen schema:\n")
  
  print(
    expected_columns
  )
  
  
  stop(
    "Canonical clinical-association schema has changed."
  )
}


if (
  !identical(
    names(clinical_xlsx),
    expected_columns
  )
) {
  
  stop(
    "Workbook duplicate schema differs from canonical frozen schema."
  )
}


cat(
  "\nFrozen 17-column schema audit passed.\n"
)


# =============================================================================
# 8. SOURCE-EQUIVALENCE AUDIT
# =============================================================================

equivalence <- compare_tables_equivalent(
  clinical_csv,
  clinical_xlsx
)


equivalence_summary <- data.frame(
  
  audit = c(
    "Dimensions identical",
    "Column names/order identical",
    "Rowwise content equivalent",
    "Content equivalent allowing row reordering",
    "Overall source equivalence"
  ),
  
  result = c(
    equivalence$dimensions_match,
    equivalence$columns_match,
    equivalence$rowwise_match,
    equivalence$reordered_match,
    equivalence$equivalent
  ),
  
  stringsAsFactors = FALSE
)


cat("\nDUPLICATE-SOURCE EQUIVALENCE AUDIT\n")
cat("----------------------------------\n")


print(
  equivalence_summary,
  row.names = FALSE
)


cat("\nCOLUMN-LEVEL EQUIVALENCE AUDIT\n")
cat("------------------------------\n")


print(
  equivalence$column_audit,
  row.names = FALSE
)


if (!equivalence$equivalent) {
  
  cat(
    "\nGenuine content differences remain after:\n"
  )
  
  cat(
    "- numeric normalization\n"
  )
  
  cat(
    "- blank / NA textual normalization\n"
  )
  
  
  if (
    nrow(equivalence$differences) >
    0
  ) {
    
    cat(
      "\nFirst 50 remaining differences:\n"
    )
    
    
    print(
      utils::head(
        equivalence$differences,
        50
      ),
      row.names = FALSE
    )
  }
  
  
  stop(
    "CSV and XLSX copies are genuinely non-equivalent."
  )
}


cat(
  "\nDuplicate-source equivalence audit PASSED.\n"
)


if (equivalence$rowwise_match) {
  
  cat(
    "CSV and XLSX copies are equivalent in the same row order.\n"
  )
  
} else {
  
  cat(
    "CSV and XLSX copies are equivalent after row reordering.\n"
  )
}


cat(
  "Standalone CSV selected as canonical frozen source.\n"
)


# =============================================================================
# 9. CANONICAL FINAL TABLE
# =============================================================================

clinical_raw <- clinical_csv


if (
  nrow(clinical_raw) !=
  60
) {
  
  stop(
    paste0(
      "Expected 60 frozen tests; observed ",
      nrow(clinical_raw),
      "."
    )
  )
}


# =============================================================================
# 10. NUMERIC FIELDS
# =============================================================================

nominal_p <- safe_numeric(
  clinical_raw$p_value
)


global_bh <- safe_numeric(
  clinical_raw$BH_global
)


within_framework_bh <- safe_numeric(
  clinical_raw$BH_within_framework
)


within_family_bh <- safe_numeric(
  clinical_raw$BH_within_test_family
)


if (
  any(
    !is.finite(nominal_p)
  )
) {
  
  stop(
    "Non-finite nominal P value detected."
  )
}


if (
  any(
    !is.finite(global_bh)
  )
) {
  
  stop(
    "Non-finite global BH value detected."
  )
}


if (
  any(
    nominal_p < 0 |
    nominal_p > 1
  )
) {
  
  stop(
    "Nominal P value outside [0,1]."
  )
}


if (
  any(
    global_bh < 0 |
    global_bh > 1
  )
) {
  
  stop(
    "Global BH value outside [0,1]."
  )
}


# =============================================================================
# 11. GLOBAL MULTIPLICITY AUDIT
# =============================================================================

global_significant <- global_bh <
  0.05


n_global_significant <- sum(
  global_significant
)


cat("\nGLOBAL MULTIPLICITY AUDIT\n")
cat("-------------------------\n")


cat(
  "Evaluable tests = ",
  nrow(clinical_raw),
  "\n",
  sep = ""
)


cat(
  "Global-BH significant tests (q<0.05) = ",
  n_global_significant,
  "\n",
  sep = ""
)


if (
  n_global_significant !=
  2
) {
  
  cat(
    "\nUnexpected globally significant rows:\n"
  )
  
  
  print(
    clinical_raw[
      global_significant,
      ,
      drop = FALSE
    ]
  )
  
  
  stop(
    paste0(
      "Expected exactly two globally significant tests; observed ",
      n_global_significant,
      "."
    )
  )
}


# =============================================================================
# 12. IDENTIFY SIGNIFICANT ASSOCIATIONS
# =============================================================================

framework_text <- tolower(
  as.character(
    clinical_raw$framework
  )
)


clinical_text <- tolower(
  paste(
    clinical_raw$clinical_variable,
    clinical_raw$clinical_label,
    sep = " | "
  )
)


primary_flag <-
  grepl(
    "primary",
    framework_text
  ) &
  grepl(
    "score|5gene|5_gene|five",
    framework_text
  )


srsq_flag <- grepl(
  "srsq",
  framework_text
)


crp_flag <- grepl(
  "(^|[^a-z])crp([^a-z]|$)|c-reactive|c reactive",
  clinical_text,
  perl = TRUE
)


significant_idx <- which(
  global_significant
)


cat("\nGLOBAL-BH SIGNIFICANT ASSOCIATIONS\n")
cat("----------------------------------\n")


print(
  clinical_raw[
    significant_idx,
    ,
    drop = FALSE
  ],
  row.names = FALSE
)


if (
  !all(
    crp_flag[significant_idx]
  )
) {
  
  stop(
    "At least one globally significant association is not CRP-related."
  )
}


if (
  sum(
    primary_flag[significant_idx]
  ) !=
  1
) {
  
  stop(
    "Expected exactly one Primary-score CRP association."
  )
}


if (
  sum(
    srsq_flag[significant_idx]
  ) !=
  1
) {
  
  stop(
    "Expected exactly one SRSq CRP association."
  )
}


cat(
  "\nSignificant-association identity audit passed.\n"
)


# =============================================================================
# 13. CRP NUMERIC ANCHORS
# =============================================================================

effect_value <- safe_numeric(
  clinical_raw$effect
)


effect_name_text <- tolower(
  as.character(
    clinical_raw$effect_name
  )
)


primary_idx <- significant_idx[
  primary_flag[significant_idx]
]


srsq_idx <- significant_idx[
  srsq_flag[significant_idx]
]


if (
  length(primary_idx) !=
  1 ||
  length(srsq_idx) !=
  1
) {
  
  stop(
    "Could not uniquely identify frozen CRP anchor rows."
  )
}


if (
  !grepl(
    "rho|spearman|correlation",
    effect_name_text[primary_idx]
  )
) {
  
  stop(
    "Primary-score CRP effect is not labelled as correlation/rho."
  )
}


if (
  !grepl(
    "rho|spearman|correlation",
    effect_name_text[srsq_idx]
  )
) {
  
  stop(
    "SRSq CRP effect is not labelled as correlation/rho."
  )
}


primary_rho <- effect_value[primary_idx]
primary_p <- nominal_p[primary_idx]
primary_q <- global_bh[primary_idx]

srsq_rho <- effect_value[srsq_idx]
srsq_p <- nominal_p[srsq_idx]
srsq_q <- global_bh[srsq_idx]


cat("\nCRP NUMERIC ANCHORS\n")
cat("-------------------\n")


cat(
  "Primary score vs CRP: rho = ",
  primary_rho,
  "; P = ",
  primary_p,
  "; global BH = ",
  primary_q,
  "\n",
  sep = ""
)


cat(
  "SRSq vs CRP: rho = ",
  srsq_rho,
  "; P = ",
  srsq_p,
  "; global BH = ",
  srsq_q,
  "\n",
  sep = ""
)


if (
  abs(
    primary_rho -
    0.574
  ) >
  0.02
) {
  
  stop(
    "Primary-score vs CRP rho anchor failed."
  )
}


if (
  abs(
    primary_p -
    3.09e-4
  ) >
  1e-4
) {
  
  stop(
    "Primary-score vs CRP nominal-P anchor failed."
  )
}


if (
  abs(
    primary_q -
    0.0185
  ) >
  0.005
) {
  
  stop(
    "Primary-score vs CRP global-BH anchor failed."
  )
}


if (
  abs(
    srsq_rho -
    0.526
  ) >
  0.02
) {
  
  stop(
    "SRSq vs CRP rho anchor failed."
  )
}


if (
  abs(
    srsq_p -
    0.00117
  ) >
  4e-4
) {
  
  stop(
    "SRSq vs CRP nominal-P anchor failed."
  )
}


if (
  abs(
    srsq_q -
    0.0352
  ) >
  0.008
) {
  
  stop(
    "SRSq vs CRP global-BH anchor failed."
  )
}


cat(
  "Frozen CRP numeric-anchor audit PASSED.\n"
)


# =============================================================================
# 14. TEST-FAMILY STRUCTURE
# =============================================================================

test_family_summary <- clinical_raw %>%
  
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
  test_family_summary,
  row.names = FALSE
)


# =============================================================================
# 15. BUILD SUBMISSION TABLE
# =============================================================================

clinical_submission <- clinical_raw


clinical_submission$Global_BH_significant <-
  global_significant


clinical_submission$Multiplicity_family <-
  "Global exploratory family: 60 evaluable tests"


clinical_submission$Interpretation_level <-
  ifelse(
    global_significant,
    "Global-BH significant exploratory association",
    "Not significant after global BH correction"
  )


# =============================================================================
# 16. CONTINUOUS / CATEGORICAL SPLIT
# =============================================================================

test_family_lower <- tolower(
  as.character(
    clinical_submission$test_family
  )
)


test_lower <- tolower(
  as.character(
    clinical_submission$test
  )
)


continuous_flag <-
  grepl(
    "continuous|numeric|spearman|correlation",
    test_family_lower
  ) |
  grepl(
    "spearman|correlation",
    test_lower
  )


categorical_flag <-
  grepl(
    "categorical|binary|group",
    test_family_lower
  ) |
  grepl(
    "wilcoxon|kruskal|fisher",
    test_lower
  )


continuous_table <- clinical_submission[
  continuous_flag,
  ,
  drop = FALSE
]


categorical_table <- clinical_submission[
  categorical_flag,
  ,
  drop = FALSE
]


# =============================================================================
# 17. TABLE S7 README
# =============================================================================

s7_readme <- data.frame(
  
  Item = c(
    "Title",
    "Scope",
    "Canonical source",
    "Duplicate source",
    "Duplicate-source equivalence",
    "Multiplicity family",
    "Primary multiplicity interpretation",
    "Global-BH significant findings",
    "Within-framework BH",
    "Within-test-family BH",
    "Effect-size reporting",
    "Clinical interpretation",
    "Important limitation"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S7. Complete exploratory clinical associations ",
      "of blood transcriptomic host-response measures."
    ),
    
    paste0(
      "Complete frozen family of 60 evaluable exploratory clinical-association ",
      "tests from the final Script 136b analysis."
    ),
    
    basename(
      canonical_csv_file
    ),
    
    paste0(
      basename(workbook_file),
      "::",
      workbook_sheet
    ),
    
    paste0(
      "The CSV and workbook copies were compared column by column. Numeric ",
      "values were compared numerically, and blank textual cells were treated ",
      "as equivalent to missing values. No substantive source differences ",
      "were permitted."
    ),
    
    "60 evaluable tests.",
    
    paste0(
      "Benjamini-Hochberg correction across all 60 evaluable tests is the ",
      "primary multiplicity interpretation."
    ),
    
    paste0(
      "Two associations remained significant after global BH correction: ",
      "the primary five-gene host-response score versus CRP and SRSq versus CRP."
    ),
    
    paste0(
      "BH_within_framework values are retained as supplementary multiplicity ",
      "information and do not replace the global 60-test correction."
    ),
    
    paste0(
      "BH_within_test_family values are retained as supplementary multiplicity ",
      "information and do not replace the global 60-test correction."
    ),
    
    paste0(
      "The effect and effect_name fields retain the effect metric generated ",
      "by the frozen source analysis; continuous correlations use Spearman rho."
    ),
    
    paste0(
      "The globally significant CRP associations support a relationship ",
      "between the transcriptomic host-response axis and systemic inflammatory ",
      "activity."
    ),
    
    paste0(
      "These exploratory analyses do not constitute independent diagnostic, ",
      "prognostic, causal, or treatment-response validation."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 18. MULTIPLICITY SUMMARY
# =============================================================================

multiplicity_summary <- data.frame(
  
  Metric = c(
    "Evaluable tests",
    "Global-BH significant tests",
    "Global-BH nonsignificant tests",
    "Primary score vs CRP rho",
    "Primary score vs CRP P",
    "Primary score vs CRP global BH",
    "SRSq vs CRP rho",
    "SRSq vs CRP P",
    "SRSq vs CRP global BH"
  ),
  
  Value = c(
    60,
    2,
    58,
    primary_rho,
    primary_p,
    primary_q,
    srsq_rho,
    srsq_p,
    srsq_q
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 19. SOURCE MANIFEST
# =============================================================================

source_manifest <- data.frame(
  
  Component = c(
    "Canonical final 60-test source",
    "Equivalent workbook duplicate",
    "Workbook sheet"
  ),
  
  Source = c(
    
    normalizePath(
      canonical_csv_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      workbook_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    workbook_sheet
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 20. OUTPUT PATHS
# =============================================================================

submission_file <- file.path(
  tables_dir,
  "159_TableS7_complete_exploratory_clinical_associations.xlsx"
)


audit_file <- file.path(
  audit_dir,
  "159_INTERNAL_AUDIT_TableS7_clinical_associations.xlsx"
)


note_file <- file.path(
  text_dir,
  "159_TableS7_title_and_note_EN.txt"
)


# =============================================================================
# 21. EXCEL STYLES
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


significant_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#FFF2CC"
)


# =============================================================================
# 22. WRITE SUBMISSION WORKBOOK
# =============================================================================

wb <- openxlsx::createWorkbook()


submission_objects <- list(
  
  S7_ReadMe =
    s7_readme,
  
  Complete_60_tests =
    clinical_submission,
  
  Global_BH_significant =
    clinical_submission[
      global_significant,
      ,
      drop = FALSE
    ],
  
  Continuous_associations =
    continuous_table,
  
  Categorical_associations =
    categorical_table,
  
  Test_family_summary =
    test_family_summary,
  
  Multiplicity_summary =
    multiplicity_summary
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
      "S7_ReadMe"
    ) {
      readme_header_style
    } else {
      header_style
    },
    rows = 1,
    cols = seq_len(
      ncol(data_object)
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
      ncol(data_object)
    ),
    widths = "auto"
  )
}


openxlsx::setColWidths(
  wb,
  "S7_ReadMe",
  cols = 1,
  widths = 33
)


openxlsx::setColWidths(
  wb,
  "S7_ReadMe",
  cols = 2,
  widths = 95
)


significant_excel_rows <- which(
  global_significant
) +
  1


openxlsx::addStyle(
  wb,
  "Complete_60_tests",
  significant_style,
  rows = significant_excel_rows,
  cols = seq_len(
    ncol(clinical_submission)
  ),
  gridExpand = TRUE
)


openxlsx::saveWorkbook(
  wb,
  submission_file,
  overwrite = TRUE
)


# =============================================================================
# 23. INTERNAL AUDIT WORKBOOK
# =============================================================================

wb_audit <- openxlsx::createWorkbook()


audit_objects <- list(
  
  Source_equivalence =
    equivalence_summary,
  
  Column_equivalence =
    equivalence$column_audit,
  
  Source_manifest =
    source_manifest,
  
  Frozen_schema =
    data.frame(
      position =
        seq_along(
          expected_columns
        ),
      column =
        expected_columns,
      stringsAsFactors =
        FALSE
    ),
  
  Significant_rows =
    clinical_submission[
      global_significant,
      ,
      drop = FALSE
    ],
  
  Canonical_CSV =
    clinical_csv,
  
  Workbook_copy =
    clinical_xlsx
)


if (
  nrow(equivalence$differences) >
  0
) {
  
  audit_objects$Remaining_differences <-
    equivalence$differences
}


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
  
  
  openxlsx::addStyle(
    wb_audit,
    sheet_name,
    header_style,
    rows = 1,
    cols = seq_len(
      ncol(data_object)
    ),
    gridExpand = TRUE
  )
  
  
  openxlsx::freezePane(
    wb_audit,
    sheet_name,
    firstActiveRow = 2
  )
  
  
  openxlsx::setColWidths(
    wb_audit,
    sheet_name,
    cols = seq_len(
      ncol(data_object)
    ),
    widths = "auto"
  )
}


openxlsx::saveWorkbook(
  wb_audit,
  audit_file,
  overwrite = TRUE
)


# =============================================================================
# 24. TABLE NOTE
# =============================================================================

table_note <- c(
  
  paste0(
    "Supplementary Table S7. Complete exploratory clinical associations ",
    "of blood transcriptomic host-response measures."
  ),
  
  "",
  
  paste0(
    "The table reports the complete frozen family of 60 evaluable exploratory ",
    "clinical-association tests. The standalone CSV source was verified against ",
    "the corresponding workbook copy using numeric-aware comparison and ",
    "normalization of blank versus missing textual cells. Benjamini-Hochberg ",
    "correction across the full 60-test family is used as the primary ",
    "multiplicity interpretation. Two associations remained significant after ",
    "global correction: the primary five-gene host-response score versus ",
    "C-reactive protein and SRSq versus C-reactive protein. Within-framework ",
    "and within-test-family adjusted values are retained as supplementary ",
    "information. These analyses provide exploratory molecular-clinical context ",
    "and do not constitute independent diagnostic or prognostic validation."
  )
)


writeLines(
  table_note,
  note_file
)


# =============================================================================
# 25. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "159_sessionInfo.txt"
  )
)


# =============================================================================
# 26. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 159 FINAL v6 completed successfully.\n")
cat("====================================================================\n\n")


cat("SOURCE-EQUIVALENCE STATUS\n")
cat("-------------------------\n")


print(
  equivalence_summary,
  row.names = FALSE
)


cat("\nMULTIPLICITY\n")
cat("------------\n")


cat(
  "Evaluable tests = 60\n"
)


cat(
  "Global BH q<0.05 = 2\n"
)


cat("\nCRP ANCHORS\n")
cat("-----------\n")


cat(
  "Primary score vs CRP: rho = ",
  primary_rho,
  "; P = ",
  primary_p,
  "; global BH = ",
  primary_q,
  "\n",
  sep = ""
)


cat(
  "SRSq vs CRP: rho = ",
  srsq_rho,
  "; P = ",
  srsq_p,
  "; global BH = ",
  srsq_q,
  "\n",
  sep = ""
)


cat("\nGLOBAL-BH SIGNIFICANT ASSOCIATIONS\n")
cat("----------------------------------\n")


print(
  clinical_submission[
    global_significant,
    ,
    drop = FALSE
  ],
  row.names = FALSE
)


cat("\nTEST-FAMILY STRUCTURE\n")
cat("---------------------\n")


print(
  test_family_summary,
  row.names = FALSE
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Supplementary Table S7:\n  ",
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


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")


cat(
  "- Blank and NA textual cells are equivalent missing values.\n"
)


cat(
  "- Non-missing textual disagreements remain hard errors.\n"
)


cat(
  "- Numeric fields are compared numerically.\n"
)


cat(
  "- No clinical association or FDR value is recalculated.\n"
)


cat(
  "- Complete frozen multiplicity family = 60 evaluable tests.\n"
)


cat(
  "- Primary multiplicity interpretation = global BH across all 60 tests.\n"
)


cat(
  "- Exactly two globally significant results are expected.\n"
)


cat(
  "- These must be Primary score vs CRP and SRSq vs CRP.\n"
)


cat(
  "- Clinical associations remain exploratory.\n"
)


cat("\nDone.\n")