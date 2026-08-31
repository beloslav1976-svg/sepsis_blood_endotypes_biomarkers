################################################################################
# Script 161
# FINAL v3
#
# Supplementary Table S8
#
# Published-signature benchmarking
#
# Project:
#   Sepsis_DESeq2
#
#
# PURPOSE
# -------
#
# Package the FROZEN Script 137 blood-only benchmarking results into
# Supplementary Table S8.
#
#
# CANONICAL FROZEN SOURCES
# ------------------------
#
# Correlations:
#
#   137_signature_correlations_primary_SRSq.csv
#
# Frozen schema:
#
#   score_column
#   target
#   n
#   Spearman_rho
#   p_value
#   display_name
#   panel_origin
#   same_cohort_selected
#   available
#   BH_correlation
#
#
# Endotype associations:
#
#   137_signature_endotype_associations.csv
#
# Frozen schema:
#
#   score_column
#   framework
#   test
#   n
#   effect
#   effect_name
#   p_value
#   display_name
#   panel_origin
#   same_cohort_selected
#   available
#   BH_endotype
#
#
# WORKBOOK DUPLICATES
# -------------------
#
#   137_blood_biomarker_benchmarking.xlsx
#
#   07_corr_primary_SRSq
#   08_endotype_tests
#
#
# MANUSCRIPT BENCHMARKS
# ---------------------
#
# For each of seven signatures:
#
#   1. Spearman rho with SRSq
#   2. Nominal P for SRSq association
#   3. BH-adjusted P for SRSq association
#   4. CTS epsilon-squared
#   5. Nominal P for CTS association
#   6. BH-adjusted P for CTS association
#
#
# SEVEN SIGNATURES
# ----------------
#
#   LIFTS-like
#   DCAF17 five-gene alternative
#   Primary five-gene
#   FAIM3:PLAC8-related
#   MetaScore-like
#   RAPID-related / PLAC8-PLA2G7 contrast
#   SeptiCyte LAB-like
#
#
# IMPORTANT INTERPRETATION
# ------------------------
#
# This analysis evaluates convergence on a shared host-response/endotype axis.
#
# It DOES NOT test whether one signature is statistically superior to another.
#
# Descriptive ranks are included only to summarize frozen effect magnitudes.
#
#
# THIS SCRIPT DOES NOT
# --------------------
#
#   - recalculate signature scores
#   - recalculate SRSq
#   - recalculate CTS
#   - rerun Spearman correlations
#   - rerun Kruskal-Wallis tests
#   - recalculate P values
#   - recalculate BH-adjusted P values
#   - optimize score directions
#   - refit published models
#
################################################################################


cat("====================================================================\n")
cat("Running Script 161 FINAL v3\n")
cat("Supplementary Table S8\n")
cat("Published-signature benchmarking\n")
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
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "161_TableS8_published_signature_benchmarking"
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
# 4. HELPER FUNCTIONS
# =============================================================================

safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


normalize_text <- function(x) {
  
  x_chr <- as.character(x)
  
  missing_flag <-
    is.na(x_chr) |
    trimws(x_chr) == ""
  
  x_chr <- stringr::str_squish(x_chr)
  
  x_chr[missing_flag] <- NA_character_
  
  x_chr
}


standardize_signature <- function(x) {
  
  x <- as.character(x)
  
  
  dplyr::case_when(
    
    x == "LIFTS_like" ~
      "LIFTS_like",
    
    x == "DCAF17_alt_score" ~
      "DCAF17_5_gene",
    
    x == "primary_5gene_score" ~
      "Primary_5_gene",
    
    x == "FAIM3_PLAC8_like" ~
      "FAIM3_PLAC8_related",
    
    x == "Sepsis_MetaScore_like" ~
      "MetaScore_like",
    
    x == "SeptiCyte_RAPID_like" ~
      "RAPID_related",
    
    x == "SeptiCyte_LAB_like" ~
      "SeptiCyte_LAB_like",
    
    TRUE ~
      NA_character_
  )
}


manuscript_signature_name <- function(x) {
  
  dplyr::recode(
    
    x,
    
    "LIFTS_like" =
      "LIFTS-like",
    
    "DCAF17_5_gene" =
      "DCAF17 five-gene alternative",
    
    "Primary_5_gene" =
      "Primary five-gene",
    
    "FAIM3_PLAC8_related" =
      "FAIM3:PLAC8-related",
    
    "MetaScore_like" =
      "MetaScore-like",
    
    "RAPID_related" =
      "RAPID-related (PLAC8-PLA2G7 contrast)",
    
    "SeptiCyte_LAB_like" =
      "SeptiCyte LAB-like"
  )
}


compare_column_pair <- function(x, y) {
  
  x_text <- normalize_text(x)
  y_text <- normalize_text(y)
  
  
  observed <- c(
    x_text,
    y_text
  )
  
  observed <- observed[
    !is.na(observed)
  ]
  
  
  numeric_like <- FALSE
  
  
  if (length(observed) > 0) {
    
    numeric_test <- suppressWarnings(
      as.numeric(observed)
    )
    
    numeric_like <- all(
      is.finite(
        numeric_test
      )
    )
  }
  
  
  if (numeric_like) {
    
    x_num <- safe_numeric(x)
    y_num <- safe_numeric(y)
    
    x_missing <- is.na(x_num)
    y_missing <- is.na(y_num)
    
    equal <- rep(
      FALSE,
      length(x_num)
    )
    
    equal[
      x_missing &
        y_missing
    ] <- TRUE
    
    
    comparable <-
      !x_missing &
      !y_missing
    
    
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
    
    
    mode <- "numeric"
    
    
  } else {
    
    x_missing <- is.na(x_text)
    y_missing <- is.na(y_text)
    
    equal <- rep(
      FALSE,
      length(x_text)
    )
    
    equal[
      x_missing &
        y_missing
    ] <- TRUE
    
    
    comparable <-
      !x_missing &
      !y_missing
    
    
    equal[comparable] <-
      x_text[comparable] ==
      y_text[comparable]
    
    
    mode <- "text"
  }
  
  
  list(
    equal = equal,
    mode = mode
  )
}


compare_tables <- function(
    csv_data,
    xlsx_data,
    source_label
) {
  
  if (
    !identical(
      dim(csv_data),
      dim(xlsx_data)
    )
  ) {
    
    stop(
      paste0(
        source_label,
        ": CSV/XLSX dimensions differ."
      )
    )
  }
  
  
  if (
    !identical(
      names(csv_data),
      names(xlsx_data)
    )
  ) {
    
    stop(
      paste0(
        source_label,
        ": CSV/XLSX columns differ."
      )
    )
  }
  
  
  audit_list <- vector(
    "list",
    ncol(csv_data)
  )
  
  
  difference_list <- list()
  
  difference_counter <- 0L
  
  
  for (
    j in seq_len(
      ncol(csv_data)
    )
  ) {
    
    column_name <- names(csv_data)[j]
    
    
    comparison <- compare_column_pair(
      csv_data[[column_name]],
      xlsx_data[[column_name]]
    )
    
    
    n_different <- sum(
      !comparison$equal
    )
    
    
    audit_list[[j]] <- data.frame(
      
      source =
        source_label,
      
      column =
        column_name,
      
      comparison_mode =
        comparison$mode,
      
      n_rows =
        nrow(csv_data),
      
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
      
      
      for (
        one_row in bad_rows
      ) {
        
        difference_counter <-
          difference_counter + 1L
        
        
        difference_list[[difference_counter]] <- data.frame(
          
          source =
            source_label,
          
          row =
            one_row,
          
          column =
            column_name,
          
          value_CSV =
            ifelse(
              is.na(
                csv_data[[column_name]][one_row]
              ),
              "<NA>",
              as.character(
                csv_data[[column_name]][one_row]
              )
            ),
          
          value_XLSX =
            ifelse(
              is.na(
                xlsx_data[[column_name]][one_row]
              ),
              "<NA>",
              as.character(
                xlsx_data[[column_name]][one_row]
              )
            ),
          
          stringsAsFactors =
            FALSE
        )
      }
    }
  }
  
  
  audit <- dplyr::bind_rows(
    audit_list
  )
  
  
  differences <- if (
    length(
      difference_list
    ) >
    0
  ) {
    
    dplyr::bind_rows(
      difference_list
    )
    
  } else {
    
    data.frame(
      source = character(),
      row = integer(),
      column = character(),
      value_CSV = character(),
      value_XLSX = character(),
      stringsAsFactors = FALSE
    )
  }
  
  
  if (
    any(
      !audit$match
    )
  ) {
    
    cat(
      "\nSOURCE EQUIVALENCE FAILURE\n"
    )
    
    
    print(
      audit[
        !audit$match,
        ,
        drop = FALSE
      ],
      row.names = FALSE
    )
    
    
    if (nrow(differences) > 0) {
      
      print(
        utils::head(
          differences,
          40
        ),
        row.names = FALSE
      )
    }
    
    
    stop(
      paste0(
        source_label,
        ": CSV and workbook copy are not equivalent."
      )
    )
  }
  
  
  list(
    audit = audit,
    differences = differences
  )
}


# =============================================================================
# 5. CANONICAL SCRIPT 137 SOURCES
# =============================================================================

source_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "137_benchmarking",
  "tables"
)


corr_csv_file <- file.path(
  source_dir,
  "137_signature_correlations_primary_SRSq.csv"
)


endotype_csv_file <- file.path(
  source_dir,
  "137_signature_endotype_associations.csv"
)


workbook_file <- file.path(
  source_dir,
  "137_blood_biomarker_benchmarking.xlsx"
)


corr_sheet <- "07_corr_primary_SRSq"

endotype_sheet <- "08_endotype_tests"


required_files <- c(
  corr_csv_file,
  endotype_csv_file,
  workbook_file
)


if (
  any(
    !file.exists(
      required_files
    )
  )
) {
  
  cat("\nRequired files:\n")
  
  print(
    required_files
  )
  
  stop(
    "At least one canonical Script 137 source is missing."
  )
}


cat("\nCANONICAL SCRIPT 137 SOURCES\n")
cat("----------------------------\n")


cat(
  "SRSq correlations:\n  ",
  normalizePath(
    corr_csv_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Endotype associations:\n  ",
  normalizePath(
    endotype_csv_file,
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


# =============================================================================
# 6. READ SOURCES
# =============================================================================

corr_csv <- read.csv(
  corr_csv_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


endotype_csv <- read.csv(
  endotype_csv_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


corr_xlsx <- readxl::read_excel(
  workbook_file,
  sheet = corr_sheet
) %>%
  
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


endotype_xlsx <- readxl::read_excel(
  workbook_file,
  sheet = endotype_sheet
) %>%
  
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat("\nSOURCE DIMENSIONS\n")
cat("-----------------\n")


cat(
  "Correlation CSV = ",
  nrow(corr_csv),
  " x ",
  ncol(corr_csv),
  "\n",
  sep = ""
)


cat(
  "Correlation XLSX = ",
  nrow(corr_xlsx),
  " x ",
  ncol(corr_xlsx),
  "\n",
  sep = ""
)


cat(
  "Endotype CSV = ",
  nrow(endotype_csv),
  " x ",
  ncol(endotype_csv),
  "\n",
  sep = ""
)


cat(
  "Endotype XLSX = ",
  nrow(endotype_xlsx),
  " x ",
  ncol(endotype_xlsx),
  "\n",
  sep = ""
)


# =============================================================================
# 7. EXACT FROZEN SCHEMA AUDIT
# =============================================================================

expected_corr_columns <- c(
  "score_column",
  "target",
  "n",
  "Spearman_rho",
  "p_value",
  "display_name",
  "panel_origin",
  "same_cohort_selected",
  "available",
  "BH_correlation"
)


expected_endotype_columns <- c(
  "score_column",
  "framework",
  "test",
  "n",
  "effect",
  "effect_name",
  "p_value",
  "display_name",
  "panel_origin",
  "same_cohort_selected",
  "available",
  "BH_endotype"
)


if (
  !identical(
    names(corr_csv),
    expected_corr_columns
  )
) {
  
  cat("\nObserved correlation columns:\n")
  
  print(
    names(corr_csv)
  )
  
  stop(
    "Correlation source schema differs from frozen Script 137 schema."
  )
}


if (
  !identical(
    names(endotype_csv),
    expected_endotype_columns
  )
) {
  
  cat("\nObserved endotype columns:\n")
  
  print(
    names(endotype_csv)
  )
  
  stop(
    "Endotype source schema differs from frozen Script 137 schema."
  )
}


if (nrow(corr_csv) != 14) {
  
  stop(
    paste0(
      "Expected 14 frozen correlation rows; observed ",
      nrow(corr_csv),
      "."
    )
  )
}


if (nrow(endotype_csv) != 14) {
  
  stop(
    paste0(
      "Expected 14 frozen endotype rows; observed ",
      nrow(endotype_csv),
      "."
    )
  )
}


cat(
  "\nExact frozen Script 137 schemas passed.\n"
)


# =============================================================================
# 8. CSV / WORKBOOK EQUIVALENCE
# =============================================================================

corr_equivalence <- compare_tables(
  corr_csv,
  corr_xlsx,
  "SRSq correlations"
)


endotype_equivalence <- compare_tables(
  endotype_csv,
  endotype_xlsx,
  "Endotype associations"
)


source_equivalence <- dplyr::bind_rows(
  corr_equivalence$audit,
  endotype_equivalence$audit
)


cat("\nSOURCE-EQUIVALENCE AUDIT\n")
cat("------------------------\n")


print(
  source_equivalence,
  row.names = FALSE
)


cat(
  "\nBoth canonical Script 137 CSV tables are equivalent to their workbook copies.\n"
)


# =============================================================================
# 9. SELECT EXACT SEVEN SRSq ROWS
# =============================================================================

corr_srsq <- corr_csv %>%
  
  dplyr::filter(
    target ==
      "SRSq",
    available
  ) %>%
  
  dplyr::mutate(
    
    signature =
      standardize_signature(
        score_column
      ),
    
    SRSq_rho =
      safe_numeric(
        Spearman_rho
      ),
    
    SRSq_p =
      safe_numeric(
        p_value
      ),
    
    SRSq_BH =
      safe_numeric(
        BH_correlation
      )
  )


cat("\nSELECTED SRSq BENCHMARK ROWS\n")
cat("----------------------------\n")


print(
  corr_srsq[
    ,
    c(
      "score_column",
      "display_name",
      "signature",
      "n",
      "SRSq_rho",
      "SRSq_p",
      "SRSq_BH",
      "panel_origin",
      "same_cohort_selected"
    )
  ],
  row.names = FALSE
)


if (nrow(corr_srsq) != 7) {
  
  stop(
    paste0(
      "Expected seven SRSq benchmark rows; observed ",
      nrow(corr_srsq),
      "."
    )
  )
}


if (
  any(
    is.na(
      corr_srsq$signature
    )
  )
) {
  
  stop(
    "At least one SRSq score_column could not be standardized."
  )
}


if (
  length(
    unique(
      corr_srsq$signature
    )
  ) !=
  7
) {
  
  stop(
    "Selected SRSq table does not contain seven unique signatures."
  )
}


if (
  any(
    !is.finite(
      corr_srsq$SRSq_rho
    )
  )
) {
  
  stop(
    "At least one SRSq benchmark lacks finite rho."
  )
}


# =============================================================================
# 10. SELECT EXACT SEVEN CTS ROWS
# =============================================================================

endotype_cts <- endotype_csv %>%
  
  dplyr::filter(
    framework ==
      "CTS",
    test ==
      "Kruskal_Wallis",
    effect_name ==
      "epsilon_squared",
    available
  ) %>%
  
  dplyr::mutate(
    
    signature =
      standardize_signature(
        score_column
      ),
    
    CTS_epsilon_squared =
      safe_numeric(
        effect
      ),
    
    CTS_p =
      safe_numeric(
        p_value
      ),
    
    CTS_BH =
      safe_numeric(
        BH_endotype
      )
  )


cat("\nSELECTED CTS BENCHMARK ROWS\n")
cat("---------------------------\n")


print(
  endotype_cts[
    ,
    c(
      "score_column",
      "display_name",
      "signature",
      "n",
      "CTS_epsilon_squared",
      "CTS_p",
      "CTS_BH",
      "panel_origin",
      "same_cohort_selected"
    )
  ],
  row.names = FALSE
)


if (nrow(endotype_cts) != 7) {
  
  stop(
    paste0(
      "Expected seven CTS benchmark rows; observed ",
      nrow(endotype_cts),
      "."
    )
  )
}


if (
  any(
    is.na(
      endotype_cts$signature
    )
  )
) {
  
  stop(
    "At least one CTS score_column could not be standardized."
  )
}


if (
  length(
    unique(
      endotype_cts$signature
    )
  ) !=
  7
) {
  
  stop(
    "Selected CTS table does not contain seven unique signatures."
  )
}


if (
  any(
    !is.finite(
      endotype_cts$CTS_epsilon_squared
    )
  )
) {
  
  stop(
    "At least one CTS benchmark lacks finite epsilon-squared."
  )
}


# =============================================================================
# 11. CROSS-SOURCE SIGNATURE CONSISTENCY
# =============================================================================

if (
  !setequal(
    corr_srsq$signature,
    endotype_cts$signature
  )
) {
  
  cat("\nSRSq signatures:\n")
  
  print(
    sort(
      corr_srsq$signature
    )
  )
  
  
  cat("\nCTS signatures:\n")
  
  print(
    sort(
      endotype_cts$signature
    )
  )
  
  
  stop(
    "SRSq and CTS benchmark signature sets differ."
  )
}


cat(
  "\nSRSq and CTS signature sets are identical: 7/7.\n"
)


# =============================================================================
# 12. MERGE FROZEN BENCHMARKS
# =============================================================================

srsq_summary <- corr_srsq %>%
  
  dplyr::transmute(
    
    signature,
    
    source_score_column =
      score_column,
    
    source_display_name =
      display_name,
    
    panel_origin,
    
    same_cohort_selected,
    
    n_SRSq =
      safe_numeric(n),
    
    SRSq_rho,
    
    SRSq_p,
    
    SRSq_BH
  )


cts_summary <- endotype_cts %>%
  
  dplyr::transmute(
    
    signature,
    
    n_CTS =
      safe_numeric(n),
    
    CTS_epsilon_squared,
    
    CTS_p,
    
    CTS_BH
  )


benchmark_summary <- dplyr::inner_join(
  srsq_summary,
  cts_summary,
  by = "signature"
)


if (nrow(benchmark_summary) != 7) {
  
  stop(
    paste0(
      "Merged benchmark table should contain seven rows; observed ",
      nrow(benchmark_summary),
      "."
    )
  )
}


# =============================================================================
# 13. FROZEN ANCHOR AUDIT
# =============================================================================

expected_anchor <- data.frame(
  
  signature = c(
    "LIFTS_like",
    "DCAF17_5_gene",
    "Primary_5_gene",
    "FAIM3_PLAC8_related",
    "MetaScore_like",
    "RAPID_related",
    "SeptiCyte_LAB_like"
  ),
  
  expected_SRSq_rho = c(
    0.8518207,
    0.7943978,
    0.7649860,
    0.7162465,
    0.6848739,
    0.4946779,
    0.4624650
  ),
  
  expected_CTS_epsilon_squared = c(
    0.6978685,
    0.6096457,
    0.6606774,
    0.6492474,
    0.5292361,
    0.4361125,
    0.3281250
  ),
  
  stringsAsFactors = FALSE
)


anchor_audit <- expected_anchor %>%
  
  dplyr::left_join(
    benchmark_summary,
    by = "signature"
  ) %>%
  
  dplyr::mutate(
    
    SRSq_absolute_difference =
      abs(
        SRSq_rho -
          expected_SRSq_rho
      ),
    
    CTS_absolute_difference =
      abs(
        CTS_epsilon_squared -
          expected_CTS_epsilon_squared
      ),
    
    SRSq_pass =
      SRSq_absolute_difference <
      1e-6,
    
    CTS_pass =
      CTS_absolute_difference <
      1e-6,
    
    overall_pass =
      SRSq_pass &
      CTS_pass
  )


cat("\nFROZEN BENCHMARKING ANCHOR AUDIT\n")
cat("--------------------------------\n")


print(
  anchor_audit[
    ,
    c(
      "signature",
      "SRSq_rho",
      "expected_SRSq_rho",
      "SRSq_absolute_difference",
      "CTS_epsilon_squared",
      "expected_CTS_epsilon_squared",
      "CTS_absolute_difference",
      "overall_pass"
    )
  ],
  row.names = FALSE
)


if (
  !all(
    anchor_audit$overall_pass
  )
) {
  
  stop(
    "At least one frozen Script 137 anchor failed."
  )
}


cat(
  "All seven frozen benchmarking anchors passed.\n"
)


# =============================================================================
# 14. P-VALUE / BH AUDIT
# =============================================================================

numeric_probability_fields <- c(
  benchmark_summary$SRSq_p,
  benchmark_summary$SRSq_BH,
  benchmark_summary$CTS_p,
  benchmark_summary$CTS_BH
)


if (
  any(
    !is.finite(
      numeric_probability_fields
    )
  )
) {
  
  stop(
    "At least one benchmark P/BH value is non-finite."
  )
}


if (
  any(
    numeric_probability_fields < 0 |
    numeric_probability_fields > 1
  )
) {
  
  stop(
    "At least one benchmark P/BH value lies outside [0,1]."
  )
}


# =============================================================================
# 15. ADD MANUSCRIPT-FACING METADATA
# =============================================================================

benchmark_summary <- benchmark_summary %>%
  
  dplyr::mutate(
    
    Signature =
      manuscript_signature_name(
        signature
      ),
    
    Signature_type =
      dplyr::case_when(
        
        signature ==
          "Primary_5_gene" ~
          "Primary study-derived signature",
        
        signature ==
          "DCAF17_5_gene" ~
          "Alternative study-derived sensitivity signature",
        
        TRUE ~
          "Published-signature RNA-seq analogue"
      ),
    
    SRSq_rank_descriptive =
      rank(
        -SRSq_rho,
        ties.method = "min"
      ),
    
    CTS_rank_descriptive =
      rank(
        -CTS_epsilon_squared,
        ties.method = "min"
      ),
    
    Interpretation =
      dplyr::case_when(
        
        signature ==
          "SeptiCyte_LAB_like" ~
          paste0(
            "RNA-seq gene-based analogue; not the proprietary ",
            "SeptiCyte LAB clinical assay output."
          ),
        
        signature ==
          "RAPID_related" ~
          paste0(
            "RNA-seq PLAC8-PLA2G7 contrast used as the ",
            "RAPID-related comparator implementation."
          ),
        
        TRUE ~
          paste0(
            "Frozen Script 137 implementation used to assess ",
            "alignment with the shared host-response axis."
          )
      )
  )


# =============================================================================
# 16. FINAL SUBMISSION TABLE
# =============================================================================

submission_table <- benchmark_summary %>%
  
  dplyr::select(
    
    Signature,
    
    Signature_type,
    
    Source_score_column =
      source_score_column,
    
    Panel_origin =
      panel_origin,
    
    Same_cohort_selected =
      same_cohort_selected,
    
    n =
      n_SRSq,
    
    SRSq_Spearman_rho =
      SRSq_rho,
    
    SRSq_nominal_P =
      SRSq_p,
    
    SRSq_BH_adjusted_P =
      SRSq_BH,
    
    CTS_epsilon_squared,
    
    CTS_nominal_P =
      CTS_p,
    
    CTS_BH_adjusted_P =
      CTS_BH,
    
    SRSq_rank_descriptive,
    
    CTS_rank_descriptive,
    
    Interpretation
  ) %>%
  
  dplyr::arrange(
    SRSq_rank_descriptive
  )


cat("\nFINAL TABLE S8 BENCHMARK SUMMARY\n")
cat("--------------------------------\n")


print(
  submission_table,
  row.names = FALSE
)


# =============================================================================
# 17. CORE MANUSCRIPT AUDIT
# =============================================================================

if (
  !all(
    submission_table$n ==
    35
  )
) {
  
  stop(
    "Expected n=35 for all benchmark signatures."
  )
}


if (
  any(
    abs(
      submission_table$SRSq_Spearman_rho
    ) >
    1
  )
) {
  
  stop(
    "Invalid Spearman rho outside [-1,1]."
  )
}


if (
  any(
    submission_table$CTS_epsilon_squared < 0 |
    submission_table$CTS_epsilon_squared > 1
  )
) {
  
  stop(
    "Invalid CTS epsilon-squared outside [0,1]."
  )
}


# =============================================================================
# 18. MULTIPLE-TESTING SUMMARY
# =============================================================================

significance_summary <- data.frame(
  
  Benchmark = c(
    "SRSq alignment",
    "CTS alignment"
  ),
  
  n_tests = c(
    7,
    7
  ),
  
  n_BH_lt_0_05 = c(
    
    sum(
      submission_table$SRSq_BH_adjusted_P <
        0.05
    ),
    
    sum(
      submission_table$CTS_BH_adjusted_P <
        0.05
    )
  ),
  
  stringsAsFactors = FALSE
)


cat("\nBENCHMARK SIGNIFICANCE SUMMARY\n")
cat("------------------------------\n")


print(
  significance_summary,
  row.names = FALSE
)


# =============================================================================
# 19. OPTIONAL REFERENCE REGISTRY
# =============================================================================

reference_file <- file.path(
  source_dir,
  "137_signature_reference_registry.csv"
)


if (
  file.exists(
    reference_file
  )
) {
  
  reference_registry <- read.csv(
    reference_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
} else {
  
  reference_registry <- data.frame(
    Note =
      "137_signature_reference_registry.csv not found.",
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 20. OPTIONAL GENE COMPLETENESS
# =============================================================================

gene_completeness_file <- file.path(
  source_dir,
  "137_signature_gene_completeness.csv"
)


if (
  file.exists(
    gene_completeness_file
  )
) {
  
  gene_completeness <- read.csv(
    gene_completeness_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
} else {
  
  gene_completeness <- data.frame(
    Note =
      "137_signature_gene_completeness.csv not found.",
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 21. README
# =============================================================================

s8_readme <- data.frame(
  
  Item = c(
    
    "Title",
    
    "Scope",
    
    "Canonical correlation source",
    
    "Canonical endotype source",
    
    "Source verification",
    
    "SRSq benchmark",
    
    "CTS benchmark",
    
    "Multiplicity",
    
    "Descriptive ranks",
    
    "Primary purpose",
    
    "Published comparator terminology",
    
    "SeptiCyte LAB-like terminology",
    
    "RAPID-related terminology",
    
    "Important limitation"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S8. Alignment of the primary five-gene ",
      "host-response signature and comparator transcriptomic signatures ",
      "with SRS and CTS endotype structure."
    ),
    
    paste0(
      "Blood-only frozen benchmarking results generated by Script 137."
    ),
    
    basename(
      corr_csv_file
    ),
    
    basename(
      endotype_csv_file
    ),
    
    paste0(
      "Each canonical CSV source was verified column-by-column against ",
      "the corresponding worksheet in ",
      "137_blood_biomarker_benchmarking.xlsx."
    ),
    
    "Spearman correlation with the continuous SRSq output.",
    
    paste0(
      "Kruskal-Wallis epsilon-squared effect size across CTS classes."
    ),
    
    paste0(
      "BH-adjusted values generated in frozen Script 137 are retained. ",
      "No multiple-testing correction is recalculated here."
    ),
    
    paste0(
      "Ranks summarize frozen effect magnitudes only and are not formal ",
      "pairwise tests of superiority."
    ),
    
    paste0(
      "The analysis evaluates convergence of study-derived and published ",
      "transcriptomic signatures on a shared host-response/endotype axis."
    ),
    
    paste0(
      "Published comparators are represented by the RNA-seq gene-based ",
      "implementations used in Script 137."
    ),
    
    paste0(
      "SeptiCyte LAB-like denotes the study RNA-seq analogue and not the ",
      "proprietary clinical assay output."
    ),
    
    paste0(
      "The RAPID-related comparator is represented by the frozen ",
      "PLAC8-PLA2G7 contrast used in Script 137."
    ),
    
    paste0(
      "All benchmarking was performed within the same discovery cohort ",
      "(n=35 sepsis samples); these comparisons do not constitute ",
      "independent external validation."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 22. INTERPRETIVE SUMMARY
# =============================================================================

top_srsq_idx <- which.max(
  submission_table$SRSq_Spearman_rho
)


top_cts_idx <- which.max(
  submission_table$CTS_epsilon_squared
)


interpretation_table <- data.frame(
  
  Finding = c(
    
    "Strongest SRSq alignment",
    
    "Strongest CTS alignment",
    
    "Primary five-gene SRSq alignment",
    
    "Primary five-gene CTS alignment",
    
    "Published-signature convergence",
    
    "Interpretive boundary"
  ),
  
  Result = c(
    
    paste0(
      submission_table$Signature[top_srsq_idx],
      ": rho = ",
      sprintf(
        "%.3f",
        submission_table$SRSq_Spearman_rho[top_srsq_idx]
      )
    ),
    
    paste0(
      submission_table$Signature[top_cts_idx],
      ": epsilon-squared = ",
      sprintf(
        "%.3f",
        submission_table$CTS_epsilon_squared[top_cts_idx]
      )
    ),
    
    paste0(
      "Primary five-gene: rho = ",
      sprintf(
        "%.3f",
        submission_table$SRSq_Spearman_rho[
          submission_table$Signature ==
            "Primary five-gene"
        ]
      )
    ),
    
    paste0(
      "Primary five-gene: epsilon-squared = ",
      sprintf(
        "%.3f",
        submission_table$CTS_epsilon_squared[
          submission_table$Signature ==
            "Primary five-gene"
        ]
      )
    ),
    
    paste0(
      "All seven benchmark implementations showed positive SRSq ",
      "correlations and substantial CTS-associated effect sizes, supporting ",
      "convergence on a shared host-response structure."
    ),
    
    paste0(
      "Differences in rho or epsilon-squared are descriptive and do not ",
      "constitute formal evidence that one signature is superior to another."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 23. OUTPUT PATHS
# =============================================================================

submission_file <- file.path(
  tables_dir,
  "161_TableS8_published_signature_benchmarking.xlsx"
)


audit_file <- file.path(
  audit_dir,
  "161_INTERNAL_AUDIT_TableS8_benchmarking.xlsx"
)


note_file <- file.path(
  text_dir,
  "161_TableS8_title_and_note_EN.txt"
)


results_file <- file.path(
  text_dir,
  "161_proposed_Results_3.7_signature_benchmarking_EN.txt"
)


# =============================================================================
# 24. EXCEL STYLES
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
  wrapText = TRUE,
  valign = "center"
)


# =============================================================================
# 25. WRITE TABLE S8
# =============================================================================

wb <- openxlsx::createWorkbook()


submission_objects <- list(
  
  S8_ReadMe =
    s8_readme,
  
  Benchmark_summary =
    submission_table,
  
  Significance_summary =
    significance_summary,
  
  Interpretation =
    interpretation_table,
  
  Reference_registry =
    reference_registry,
  
  Gene_completeness =
    gene_completeness
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
      "S8_ReadMe"
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
  "S8_ReadMe",
  cols = 1,
  widths = 34
)


openxlsx::setColWidths(
  wb,
  "S8_ReadMe",
  cols = 2,
  widths = 95
)


openxlsx::setColWidths(
  wb,
  "Benchmark_summary",
  cols = 15,
  widths = 65
)


openxlsx::saveWorkbook(
  wb,
  submission_file,
  overwrite = TRUE
)


# =============================================================================
# 26. WRITE INTERNAL AUDIT WORKBOOK
# =============================================================================

wb_audit <- openxlsx::createWorkbook()


audit_objects <- list(
  
  Source_equivalence =
    source_equivalence,
  
  Correlation_source =
    corr_csv,
  
  Endotype_source =
    endotype_csv,
  
  Selected_SRSq =
    corr_srsq,
  
  Selected_CTS =
    endotype_cts,
  
  Benchmark_merged =
    benchmark_summary,
  
  Frozen_anchor_audit =
    anchor_audit,
  
  Significance_summary =
    significance_summary
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
# 27. TABLE TITLE / NOTE
# =============================================================================

table_note <- c(
  
  paste0(
    "Supplementary Table S8. Alignment of the primary five-gene ",
    "host-response signature and comparator transcriptomic signatures ",
    "with SRS and CTS endotype structure."
  ),
  
  "",
  
  paste0(
    "Benchmarking values were reproduced directly from the frozen Script 137 ",
    "analysis without recalculation. Continuous alignment with SRSq is ",
    "summarized using Spearman rho, and categorical CTS alignment is ",
    "summarized using epsilon-squared from the Kruskal-Wallis analysis. ",
    "Nominal and Benjamini-Hochberg-adjusted P values are retained from the ",
    "original benchmarking analysis. Descriptive ranks are included only to ",
    "facilitate comparison of effect magnitudes and do not represent formal ",
    "pairwise tests of superiority. Published comparator signatures are ",
    "represented by the RNA-seq gene-based implementations used in this study. ",
    "The analysis evaluates biological convergence on a shared host-response ",
    "axis rather than diagnostic superiority."
  )
)


writeLines(
  table_note,
  note_file
)


# =============================================================================
# 28. PROPOSED RESULTS 3.7
# =============================================================================

get_value <- function(
    signature_name,
    metric_name
) {
  
  benchmark_summary[
    benchmark_summary$signature ==
      signature_name,
    metric_name,
    drop = TRUE
  ]
}


results_3_7 <- paste0(
  
  "We next examined whether the molecular continuum captured by the primary ",
  "five-gene score was specific to this gene configuration or was also ",
  "represented by other study-derived and previously published blood ",
  "transcriptomic signatures. All seven benchmark implementations showed ",
  "positive associations with continuous SRSq and significant differences ",
  "across CTS classes in the frozen discovery-cohort analysis ",
  "(Supplementary Fig. S6 and Supplementary Table S8). The strongest ",
  "alignment with SRSq was observed for the LIFTS-like score ",
  "(Spearman rho=",
  sprintf(
    "%.3f",
    get_value(
      "LIFTS_like",
      "SRSq_rho"
    )
  ),
  "), followed by the DCAF17 alternative five-gene score (rho=",
  sprintf(
    "%.3f",
    get_value(
      "DCAF17_5_gene",
      "SRSq_rho"
    )
  ),
  "), the primary five-gene score (rho=",
  sprintf(
    "%.3f",
    get_value(
      "Primary_5_gene",
      "SRSq_rho"
    )
  ),
  "), and the FAIM3:PLAC8-related score (rho=",
  sprintf(
    "%.3f",
    get_value(
      "FAIM3_PLAC8_related",
      "SRSq_rho"
    )
  ),
  "). CTS-associated effect sizes showed a broadly similar pattern, with ",
  "epsilon-squared=",
  sprintf(
    "%.3f",
    get_value(
      "LIFTS_like",
      "CTS_epsilon_squared"
    )
  ),
  " for LIFTS-like, ",
  sprintf(
    "%.3f",
    get_value(
      "Primary_5_gene",
      "CTS_epsilon_squared"
    )
  ),
  " for the primary five-gene score, and ",
  sprintf(
    "%.3f",
    get_value(
      "FAIM3_PLAC8_related",
      "CTS_epsilon_squared"
    )
  ),
  " for the FAIM3:PLAC8-related score. The primary five-gene signature ",
  "therefore did not uniquely define the observed molecular continuum; ",
  "rather, multiple transcriptomic signatures converged on the same ",
  "underlying host-response structure."
)


writeLines(
  results_3_7,
  results_file
)


# =============================================================================
# 29. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "161_sessionInfo.txt"
  )
)


# =============================================================================
# 30. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 161 FINAL v3 completed successfully.\n")
cat("====================================================================\n\n")


cat("FROZEN BENCHMARK SUMMARY\n")
cat("------------------------\n")


print(
  submission_table[
    ,
    c(
      "Signature",
      "SRSq_Spearman_rho",
      "SRSq_nominal_P",
      "SRSq_BH_adjusted_P",
      "CTS_epsilon_squared",
      "CTS_nominal_P",
      "CTS_BH_adjusted_P",
      "SRSq_rank_descriptive",
      "CTS_rank_descriptive"
    )
  ],
  row.names = FALSE
)


cat("\nBENCHMARK SIGNIFICANCE SUMMARY\n")
cat("------------------------------\n")


print(
  significance_summary,
  row.names = FALSE
)


cat("\nKEY INTERPRETATION\n")
cat("------------------\n")


cat(
  "Highest SRSq alignment: ",
  submission_table$Signature[
    which.max(
      submission_table$SRSq_Spearman_rho
    )
  ],
  " (rho = ",
  max(
    submission_table$SRSq_Spearman_rho
  ),
  ")\n",
  sep = ""
)


cat(
  "Highest CTS alignment: ",
  submission_table$Signature[
    which.max(
      submission_table$CTS_epsilon_squared
    )
  ],
  " (epsilon-squared = ",
  max(
    submission_table$CTS_epsilon_squared
  ),
  ")\n",
  sep = ""
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Supplementary Table S8:\n  ",
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
  "\n\n",
  sep = ""
)


cat(
  "Proposed Results 3.7:\n  ",
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
  "- Frozen Script 137 statistical sources are used directly.\n"
)


cat(
  "- Both CSV sources are verified against workbook duplicates.\n"
)


cat(
  "- Exactly seven SRSq and seven CTS benchmark rows are required.\n"
)


cat(
  "- No association, effect estimate, P value or BH value is recalculated.\n"
)


cat(
  "- Descriptive ranks are not formal superiority tests.\n"
)


cat(
  "- Benchmarking evaluates biological convergence, not diagnostic superiority.\n"
)


cat(
  "- Published comparators are RNA-seq gene-based implementations.\n"
)


cat(
  "- SeptiCyte LAB-like is not the proprietary clinical assay output.\n"
)


cat(
  "- RAPID-related refers to the frozen PLAC8-PLA2G7 contrast implementation.\n"
)


cat(
  "- All benchmarking remains internal to the discovery cohort.\n"
)


cat("\nDone.\n")