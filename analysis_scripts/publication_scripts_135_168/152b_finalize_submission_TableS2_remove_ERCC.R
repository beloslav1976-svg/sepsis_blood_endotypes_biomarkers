################################################################################
# Script 152b
# FINAL SUBMISSION CLEANUP
#
# Supplementary Table S2
#
# Purpose:
#   Create the final submission-ready version of Supplementary Table S2
#   containing biological/human targets only.
#
# Project:
#   Sepsis_DESeq2
#
# Manuscript:
#   Blood-only sepsis transcriptomic endotypes /
#   five-gene host-response signature
#
#
# IMPORTANT PROVENANCE
# --------------------
#
# The original blood DESeq2 workflow applied the expression prefilter:
#
#       >=10 raw counts in >=3 blood samples
#
# BEFORE exclusion of ERCC-prefixed technical-control features.
#
# Therefore:
#
#   Original frozen DESeq2 universe:
#       12,400 retained features
#
# comprising:
#
#       12,393 biological/human targets
#            7 ERCC-prefixed technical-control features
#
#
# Script 153 demonstrated that all seven retained ERCC-prefixed features were:
#
#   primary DEG          = FALSE
#   batch-adjusted DEG   = FALSE
#   robust-core DEG      = FALSE
#
#
# Therefore removing them from the SUBMISSION table:
#
#   DOES NOT change:
#
#       primary DEG total          = 2,659
#       batch-adjusted DEG total   = 4,125
#       robust-core DEG total      = 1,796
#
#
# THIS SCRIPT DOES NOT:
#
#   - re-run DESeq2
#   - recalculate fold changes
#   - recalculate P values
#   - recalculate adjusted P values
#   - redefine DEG thresholds
#   - redefine robust core
#
# It only creates a cleaner supplementary submission table.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 152b\n")
cat("Final submission cleanup of Supplementary Table S2\n")
cat("Remove ERCC-prefixed technical features from submission table only\n")
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
# 3. INPUT FILES
# =============================================================================

script152_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "152_Tables_S2_S3_DE_enrichment"
)


tableS2_source_file <- file.path(
  script152_dir,
  "tables",
  "152_TableS2_complete_blood_differential_expression.xlsx"
)


model_comparison_file <- file.path(
  project_dir,
  "results",
  "blood_BP_vs_BC_model_comparison",
  "blood_model_comparison_all_genes.csv"
)


if (
  !file.exists(
    tableS2_source_file
  )
) {
  
  stop(
    paste0(
      "Script 152 Table S2 workbook not found:\n",
      tableS2_source_file
    )
  )
}


if (
  !file.exists(
    model_comparison_file
  )
) {
  
  stop(
    paste0(
      "Frozen model-comparison source not found:\n",
      model_comparison_file
    )
  )
}


cat("\nINPUT FILES\n")
cat("-----------\n")


cat(
  "Script 152 Table S2:\n  ",
  normalizePath(
    tableS2_source_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Frozen model comparison:\n  ",
  normalizePath(
    model_comparison_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 4. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "152b_final_submission_TableS2"
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
# 5. HELPER FUNCTIONS
# =============================================================================

normalize_gene <- function(x) {
  
  toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
}


as_logical_flag <- function(x) {
  
  if (is.logical(x)) {
    
    x[is.na(x)] <- FALSE
    
    return(x)
  }
  
  
  if (is.numeric(x)) {
    
    x[is.na(x)] <- 0
    
    return(
      x != 0
    )
  }
  
  
  x2 <- toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
  
  
  x2[is.na(x2)] <- ""
  
  
  x2 %in% c(
    "TRUE",
    "T",
    "YES",
    "Y",
    "1",
    "DEG",
    "SIGNIFICANT",
    "SIG"
  )
}


find_gene_column <- function(data) {
  
  candidates <- c(
    "Gene",
    "gene",
    "gene_symbol",
    "GeneSymbol",
    "symbol"
  )
  
  
  hits <- candidates[
    candidates %in%
      names(data)
  ]
  
  
  if (
    length(
      hits
    ) ==
    0
  ) {
    
    stop(
      paste0(
        "Could not identify gene column. Available columns:\n",
        paste(
          names(data),
          collapse = ", "
        )
      )
    )
  }
  
  
  hits[1]
}


# =============================================================================
# 6. READ SCRIPT 152 TABLE S2
# =============================================================================

source_sheets <- readxl::excel_sheets(
  tableS2_source_file
)


cat("\nSource workbook sheets:\n")

print(
  source_sheets
)


if (
  !("Complete_DE" %in%
    source_sheets)
) {
  
  stop(
    "Sheet 'Complete_DE' not found in Script 152 Table S2 workbook."
  )
}


tableS2_all <- readxl::read_excel(
  tableS2_source_file,
  sheet = "Complete_DE"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "\nScript 152 Complete_DE dimensions: ",
  nrow(
    tableS2_all
  ),
  " rows x ",
  ncol(
    tableS2_all
  ),
  " columns\n",
  sep = ""
)


if (
  nrow(
    tableS2_all
  ) !=
  12400
) {
  
  stop(
    paste0(
      "Expected 12,400 rows in Script 152 Complete_DE; observed ",
      nrow(
        tableS2_all
      )
    )
  )
}


# =============================================================================
# 7. IDENTIFY ERCC-PREFIXED FEATURES
# =============================================================================

gene_col <- find_gene_column(
  tableS2_all
)


genes <- as.character(
  tableS2_all[[gene_col]]
)


ercc_flag <- grepl(
  "^ERCC[-_]",
  genes,
  ignore.case = TRUE
)


n_ercc <- sum(
  ercc_flag
)


n_biological <- sum(
  !ercc_flag
)


cat("\nFEATURE CLASSIFICATION\n")
cat("----------------------\n")


cat(
  "All retained features = ",
  nrow(
    tableS2_all
  ),
  "\n",
  sep = ""
)


cat(
  "ERCC-prefixed technical features = ",
  n_ercc,
  "\n",
  sep = ""
)


cat(
  "Biological/human targets = ",
  n_biological,
  "\n",
  sep = ""
)


if (
  n_ercc !=
  7
) {
  
  stop(
    paste0(
      "Expected 7 retained ERCC-prefixed features; observed ",
      n_ercc
    )
  )
}


if (
  n_biological !=
  12393
) {
  
  stop(
    paste0(
      "Expected 12,393 biological targets; observed ",
      n_biological
    )
  )
}


# =============================================================================
# 8. VERIFY AGAINST ORIGINAL FROZEN MODEL COMPARISON
# =============================================================================

model_df <- read.csv(
  model_comparison_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


if (
  nrow(
    model_df
  ) !=
  12400
) {
  
  stop(
    paste0(
      "Frozen model-comparison table should contain 12,400 rows; observed ",
      nrow(
        model_df
      )
    )
  )
}


model_gene_col <- find_gene_column(
  model_df
)


tableS2_gene_set <- sort(
  normalize_gene(
    tableS2_all[[gene_col]]
  )
)


model_gene_set <- sort(
  normalize_gene(
    model_df[[model_gene_col]]
  )
)


if (
  !identical(
    tableS2_gene_set,
    model_gene_set
  )
) {
  
  stop(
    paste0(
      "Script 152 Complete_DE gene set does not exactly match the ",
      "frozen model-comparison source."
    )
  )
}


cat(
  "\nGene-set audit passed: Script 152 Complete_DE exactly matches frozen model source.\n"
)


# =============================================================================
# 9. VERIFY REQUIRED FROZEN FLAGS
# =============================================================================

required_flag_columns <- c(
  "primary_DEG",
  "batch_adjusted_DEG",
  "robust_core_DEG"
)


missing_flag_columns <- setdiff(
  required_flag_columns,
  names(
    tableS2_all
  )
)


if (
  length(
    missing_flag_columns
  ) >
  0
) {
  
  stop(
    paste0(
      "Required Script 152 flag column(s) missing: ",
      paste(
        missing_flag_columns,
        collapse = ", "
      )
    )
  )
}


primary_deg <- as_logical_flag(
  tableS2_all$primary_DEG
)


batch_deg <- as_logical_flag(
  tableS2_all$batch_adjusted_DEG
)


core_deg <- as_logical_flag(
  tableS2_all$robust_core_DEG
)


# =============================================================================
# 10. AUDIT ERCC FLAGS
# =============================================================================

ercc_primary_n <- sum(
  primary_deg[
    ercc_flag
  ]
)


ercc_batch_n <- sum(
  batch_deg[
    ercc_flag
  ]
)


ercc_core_n <- sum(
  core_deg[
    ercc_flag
  ]
)


cat("\nERCC DEG AUDIT\n")
cat("--------------\n")


cat(
  "ERCC primary DEG = ",
  ercc_primary_n,
  "\n",
  sep = ""
)


cat(
  "ERCC batch-adjusted DEG = ",
  ercc_batch_n,
  "\n",
  sep = ""
)


cat(
  "ERCC robust-core DEG = ",
  ercc_core_n,
  "\n",
  sep = ""
)


if (
  ercc_primary_n !=
  0 ||
  ercc_batch_n !=
  0 ||
  ercc_core_n !=
  0
) {
  
  stop(
    "At least one ERCC-prefixed feature is classified as DEG."
  )
}


# =============================================================================
# 11. FULL VS BIOLOGICAL DEG COUNTS
# =============================================================================

full_counts <- data.frame(
  
  dataset = "Original frozen 12,400-feature DE universe",
  
  n_features =
    nrow(
      tableS2_all
    ),
  
  primary_DEG =
    sum(
      primary_deg
    ),
  
  batch_adjusted_DEG =
    sum(
      batch_deg
    ),
  
  robust_core_DEG =
    sum(
      core_deg
    ),
  
  stringsAsFactors = FALSE
)


tableS2_submission <- tableS2_all[
  !ercc_flag,
  ,
  drop = FALSE
]


submission_primary_deg <- as_logical_flag(
  tableS2_submission$primary_DEG
)


submission_batch_deg <- as_logical_flag(
  tableS2_submission$batch_adjusted_DEG
)


submission_core_deg <- as_logical_flag(
  tableS2_submission$robust_core_DEG
)


submission_counts <- data.frame(
  
  dataset = "Submission biological-target table",
  
  n_features =
    nrow(
      tableS2_submission
    ),
  
  primary_DEG =
    sum(
      submission_primary_deg
    ),
  
  batch_adjusted_DEG =
    sum(
      submission_batch_deg
    ),
  
  robust_core_DEG =
    sum(
      submission_core_deg
    ),
  
  stringsAsFactors = FALSE
)


count_audit <- dplyr::bind_rows(
  full_counts,
  submission_counts
)


cat("\nDEG COUNT AUDIT BEFORE/AFTER ERCC EXCLUSION\n")
cat("-------------------------------------------\n")


print(
  count_audit,
  row.names = FALSE
)


# =============================================================================
# 12. HARD EXPECTED RESULT AUDIT
# =============================================================================

expected_submission <- data.frame(
  
  n_features =
    12393,
  
  primary_DEG =
    2659,
  
  batch_adjusted_DEG =
    4125,
  
  robust_core_DEG =
    1796
)


observed_submission <- submission_counts %>%
  
  dplyr::select(
    n_features,
    primary_DEG,
    batch_adjusted_DEG,
    robust_core_DEG
  )


if (
  !identical(
    as.numeric(
      observed_submission[1, ]
    ),
    as.numeric(
      expected_submission[1, ]
    )
  )
) {
  
  cat("\nObserved submission counts:\n")
  
  print(
    observed_submission
  )
  
  
  cat("\nExpected submission counts:\n")
  
  print(
    expected_submission
  )
  
  
  stop(
    "Submission Table S2 hard audit failed."
  )
}


cat(
  "\nSubmission Table S2 hard audit passed successfully.\n"
)


# =============================================================================
# 13. DIRECTION COUNTS
# =============================================================================

if (
  !("primary_log2FoldChange" %in%
    names(
      tableS2_submission
    ))
) {
  
  stop(
    "primary_log2FoldChange column missing."
  )
}


if (
  !("batch_adjusted_log2FoldChange" %in%
    names(
      tableS2_submission
    ))
) {
  
  stop(
    "batch_adjusted_log2FoldChange column missing."
  )
}


primary_lfc <- suppressWarnings(
  as.numeric(
    tableS2_submission$primary_log2FoldChange
  )
)


batch_lfc <- suppressWarnings(
  as.numeric(
    tableS2_submission$batch_adjusted_log2FoldChange
  )
)


primary_up <- sum(
  submission_primary_deg &
    primary_lfc >
    0
)


primary_down <- sum(
  submission_primary_deg &
    primary_lfc <
    0
)


batch_up <- sum(
  submission_batch_deg &
    batch_lfc >
    0
)


batch_down <- sum(
  submission_batch_deg &
    batch_lfc <
    0
)


core_up <- sum(
  submission_core_deg &
    primary_lfc >
    0
)


core_down <- sum(
  submission_core_deg &
    primary_lfc <
    0
)


direction_audit <- data.frame(
  
  analysis = c(
    "Primary/simple",
    "Batch-adjusted",
    "Robust core"
  ),
  
  total = c(
    sum(
      submission_primary_deg
    ),
    sum(
      submission_batch_deg
    ),
    sum(
      submission_core_deg
    )
  ),
  
  UP = c(
    primary_up,
    batch_up,
    core_up
  ),
  
  DOWN = c(
    primary_down,
    batch_down,
    core_down
  ),
  
  expected_total = c(
    2659,
    4125,
    1796
  ),
  
  expected_UP = c(
    1660,
    2093,
    1133
  ),
  
  expected_DOWN = c(
    999,
    2032,
    663
  ),
  
  stringsAsFactors = FALSE
)


direction_audit$match <-
  direction_audit$total ==
  direction_audit$expected_total &
  direction_audit$UP ==
  direction_audit$expected_UP &
  direction_audit$DOWN ==
  direction_audit$expected_DOWN


cat("\nDIRECTION AUDIT\n")
cat("---------------\n")


print(
  direction_audit,
  row.names = FALSE
)


if (
  !all(
    direction_audit$match
  )
) {
  
  stop(
    "UP/DOWN direction audit failed after ERCC exclusion."
  )
}


# =============================================================================
# 14. ERCC INTERNAL AUDIT TABLE
# =============================================================================

ercc_audit <- tableS2_all[
  ercc_flag,
  ,
  drop = FALSE
]


ercc_audit <- ercc_audit %>%
  
  dplyr::arrange(
    .data[[gene_col]]
  )


# =============================================================================
# 15. FINAL SUBMISSION TABLE ORDER
# =============================================================================
#
# Keep the Script 152 column structure intact.
#
# Sort:
#   robust core first
#   then primary DEGs
#   then remaining tested biological targets
#
# Within these classes:
#   decreasing absolute primary log2FC
#
# This sorting is for usability only and changes no statistics.
#
# =============================================================================

tableS2_submission <- tableS2_submission %>%
  
  dplyr::mutate(
    
    .sort_class =
      dplyr::case_when(
        
        as_logical_flag(
          robust_core_DEG
        ) ~
          1L,
        
        as_logical_flag(
          primary_DEG
        ) ~
          2L,
        
        TRUE ~
          3L
      ),
    
    .sort_abs_lfc =
      abs(
        suppressWarnings(
          as.numeric(
            primary_log2FoldChange
          )
        )
      )
  ) %>%
  
  dplyr::arrange(
    .sort_class,
    dplyr::desc(
      .sort_abs_lfc
    ),
    .data[[gene_col]]
  ) %>%
  
  dplyr::select(
    -.sort_class,
    -.sort_abs_lfc
  )


# =============================================================================
# 16. FINAL README
# =============================================================================

tableS2_readme <- data.frame(
  
  Item = c(
    "Title",
    "Scope",
    "Original targeted matrix",
    "Original DESeq2 prefilter",
    "Original frozen DESeq2 analysis universe",
    "Submission table universe",
    "ERCC-prefixed technical features",
    "Primary model",
    "Batch-adjusted sensitivity model",
    "DEG threshold",
    "Primary DEG result",
    "Batch-adjusted DEG result",
    "Robust-core result",
    "Multiple testing",
    "Important provenance note"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S2. Complete biological-target differential-",
      "expression results for sepsis versus healthy-control blood."
    ),
    
    paste0(
      "Gene-level differential-expression results for all biological targets ",
      "retained in the blood discovery analysis after the prespecified ",
      "expression prefilter."
    ),
    
    paste0(
      "The complete targeted input matrix contained 20,812 rows, including ",
      "20,802 biological targets and 10 ERCC-prefixed technical-control ",
      "features."
    ),
    
    paste0(
      "Features with at least 10 raw counts in at least three blood samples ",
      "were retained."
    ),
    
    paste0(
      "The historical DESeq2 workflow retained 12,400 features: 12,393 ",
      "biological targets and seven ERCC-prefixed technical-control features."
    ),
    
    paste0(
      "This submission table contains the 12,393 retained biological targets. ",
      "The seven retained ERCC-prefixed technical-control features are omitted ",
      "from the submission table and retained only in the internal provenance ",
      "audit."
    ),
    
    paste0(
      "None of the seven retained ERCC-prefixed features was classified as ",
      "differentially expressed in the primary or batch-adjusted model, and ",
      "none entered the robust-core DEG set."
    ),
    
    "Primary DESeq2 design: ~ condition.",
    
    "Sensitivity DESeq2 design: ~ batch + condition.",
    
    paste0(
      "Differentially expressed genes were defined by Benjamini-Hochberg ",
      "adjusted P <0.05 and absolute log2 fold change >=1."
    ),
    
    "2,659 DEGs: 1,660 increased and 999 decreased in sepsis.",
    
    "4,125 DEGs: 2,093 increased and 2,032 decreased in sepsis.",
    
    paste0(
      "1,796 concordant robust-core DEGs: 1,133 increased and 663 decreased ",
      "in sepsis."
    ),
    
    "Benjamini-Hochberg false-discovery-rate correction.",
    
    paste0(
      "Removal of ERCC-prefixed rows from the submission table is a reporting ",
      "cleanup only. Differential-expression statistics were not recalculated."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 17. INTERNAL PROVENANCE SUMMARY
# =============================================================================

provenance_summary <- data.frame(
  
  metric = c(
    "Complete targeted input rows",
    "Biological targets in complete input",
    "ERCC-prefixed rows in complete input",
    "Original DESeq2 retained features",
    "Original DESeq2 retained biological targets",
    "Original DESeq2 retained ERCC-prefixed features",
    "Submission Table S2 rows",
    "ERCC rows removed from submission table",
    "Primary DEG before ERCC cleanup",
    "Primary DEG after ERCC cleanup",
    "Batch DEG before ERCC cleanup",
    "Batch DEG after ERCC cleanup",
    "Robust core before ERCC cleanup",
    "Robust core after ERCC cleanup"
  ),
  
  value = c(
    20812,
    20802,
    10,
    12400,
    12393,
    7,
    nrow(
      tableS2_submission
    ),
    nrow(
      ercc_audit
    ),
    full_counts$primary_DEG,
    submission_counts$primary_DEG,
    full_counts$batch_adjusted_DEG,
    submission_counts$batch_adjusted_DEG,
    full_counts$robust_core_DEG,
    submission_counts$robust_core_DEG
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 18. OUTPUT PATHS
# =============================================================================

submission_xlsx <- file.path(
  tables_dir,
  "152b_TableS2_FINAL_SUBMISSION_biological_targets.xlsx"
)


submission_csv <- file.path(
  tables_dir,
  "152b_TableS2_FINAL_SUBMISSION_biological_targets.csv"
)


internal_audit_xlsx <- file.path(
  audit_dir,
  "152b_INTERNAL_AUDIT_TableS2_ERCC_cleanup.xlsx"
)


readme_txt <- file.path(
  text_dir,
  "152b_TableS2_submission_note_EN.txt"
)


# =============================================================================
# 19. EXCEL STYLES
# =============================================================================

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  valign = "center",
  halign = "center",
  wrapText = TRUE
)


readme_header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#EEF3F7",
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE
)


# =============================================================================
# 20. WRITE FINAL SUBMISSION TABLE S2
# =============================================================================
#
# Submission workbook intentionally contains ONLY:
#
#   S2_ReadMe
#   Complete_DE
#
# Internal ERCC audit is NOT included in the journal-facing workbook.
#
# =============================================================================

wb_submission <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wb_submission,
  "S2_ReadMe"
)


openxlsx::addWorksheet(
  wb_submission,
  "Complete_DE"
)


openxlsx::writeData(
  wb_submission,
  "S2_ReadMe",
  tableS2_readme
)


openxlsx::writeData(
  wb_submission,
  "Complete_DE",
  tableS2_submission,
  withFilter = TRUE
)


openxlsx::addStyle(
  wb_submission,
  "S2_ReadMe",
  readme_header_style,
  rows = 1,
  cols = 1:ncol(
    tableS2_readme
  ),
  gridExpand = TRUE
)


openxlsx::addStyle(
  wb_submission,
  "Complete_DE",
  header_style,
  rows = 1,
  cols = 1:ncol(
    tableS2_submission
  ),
  gridExpand = TRUE
)


openxlsx::freezePane(
  wb_submission,
  "S2_ReadMe",
  firstActiveRow = 2
)


openxlsx::freezePane(
  wb_submission,
  "Complete_DE",
  firstActiveRow = 2,
  firstActiveCol = 2
)


openxlsx::setColWidths(
  wb_submission,
  "S2_ReadMe",
  cols = 1,
  widths = 32
)


openxlsx::setColWidths(
  wb_submission,
  "S2_ReadMe",
  cols = 2,
  widths = 90
)


for (
  col_idx in seq_len(
    ncol(
      tableS2_submission
    )
  )
) {
  
  desired_width <- max(
    11,
    nchar(
      names(
        tableS2_submission
      )[col_idx]
    ) +
      2
  )
  
  
  desired_width <- min(
    desired_width,
    26
  )
  
  
  openxlsx::setColWidths(
    wb_submission,
    "Complete_DE",
    cols = col_idx,
    widths = desired_width
  )
}


openxlsx::saveWorkbook(
  wb_submission,
  submission_xlsx,
  overwrite = TRUE
)


write.csv(
  tableS2_submission,
  submission_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 21. WRITE INTERNAL AUDIT WORKBOOK
# =============================================================================

wb_audit <- openxlsx::createWorkbook()


audit_objects <- list(
  
  Provenance_summary =
    provenance_summary,
  
  Count_audit =
    count_audit,
  
  Direction_audit =
    direction_audit,
  
  Excluded_ERCC =
    ercc_audit
)


for (
  sheet_name in names(
    audit_objects
  )
) {
  
  data_object <-
    audit_objects[[sheet_name]]
  
  
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
    cols = 1:ncol(
      data_object
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
    cols = 1:ncol(
      data_object
    ),
    widths = "auto"
  )
}


openxlsx::saveWorkbook(
  wb_audit,
  internal_audit_xlsx,
  overwrite = TRUE
)


# =============================================================================
# 22. WRITE MANUSCRIPT NOTE
# =============================================================================

submission_note <- paste0(
  
  "Supplementary Table S2 contains complete differential-expression results ",
  "for the 12,393 biological targets retained after the prespecified blood ",
  "expression prefilter. The original DESeq2 workflow retained 12,400 ",
  "features because seven ERCC-prefixed technical-control features also ",
  "passed the expression threshold. None of these seven technical features ",
  "was differentially expressed in the primary or batch-adjusted model, and ",
  "none entered the robust-core DEG set. They were therefore omitted from ",
  "the journal-facing biological supplementary table without recalculation ",
  "of any differential-expression statistic."
)


writeLines(
  submission_note,
  readme_txt
)


# =============================================================================
# 23. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "152b_sessionInfo.txt"
  )
)


# =============================================================================
# 24. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 152b completed successfully.\n")
cat("====================================================================\n\n")


cat("ORIGINAL FROZEN DE UNIVERSE\n")
cat("---------------------------\n")


cat(
  "Total retained features = ",
  nrow(
    tableS2_all
  ),
  "\n",
  sep = ""
)


cat(
  "Biological targets = ",
  n_biological,
  "\n",
  sep = ""
)


cat(
  "ERCC-prefixed technical features = ",
  n_ercc,
  "\n",
  sep = ""
)


cat("\nERCC AUDIT\n")
cat("----------\n")


cat(
  "Primary DEG among ERCC = ",
  ercc_primary_n,
  "\n",
  sep = ""
)


cat(
  "Batch-adjusted DEG among ERCC = ",
  ercc_batch_n,
  "\n",
  sep = ""
)


cat(
  "Robust-core DEG among ERCC = ",
  ercc_core_n,
  "\n",
  sep = ""
)


cat("\nFINAL SUBMISSION TABLE S2\n")
cat("-------------------------\n")


cat(
  "Biological targets = ",
  nrow(
    tableS2_submission
  ),
  "\n",
  sep = ""
)


cat(
  "Primary/simple DEG = ",
  sum(
    submission_primary_deg
  ),
  " (UP ",
  primary_up,
  "; DOWN ",
  primary_down,
  ")\n",
  sep = ""
)


cat(
  "Batch-adjusted DEG = ",
  sum(
    submission_batch_deg
  ),
  " (UP ",
  batch_up,
  "; DOWN ",
  batch_down,
  ")\n",
  sep = ""
)


cat(
  "Robust core = ",
  sum(
    submission_core_deg
  ),
  " (UP ",
  core_up,
  "; DOWN ",
  core_down,
  ")\n",
  sep = ""
)


cat("\nBEFORE/AFTER AUDIT\n")
cat("------------------\n")


print(
  count_audit,
  row.names = FALSE
)


cat("\nDIRECTION AUDIT\n")
cat("---------------\n")


print(
  direction_audit,
  row.names = FALSE
)


cat("\nEXCLUDED ERCC FEATURES\n")
cat("----------------------\n")


print(
  ercc_audit[
    ,
    intersect(
      c(
        "Gene",
        "primary_log2FoldChange",
        "primary_padj",
        "primary_DEG",
        "batch_adjusted_log2FoldChange",
        "batch_adjusted_padj",
        "batch_adjusted_DEG",
        "robust_core_DEG"
      ),
      names(
        ercc_audit
      )
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "FINAL submission Table S2:\n  ",
  normalizePath(
    submission_xlsx,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "FINAL submission Table S2 CSV:\n  ",
  normalizePath(
    submission_csv,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Internal ERCC audit:\n  ",
  normalizePath(
    internal_audit_xlsx,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Manuscript note:\n  ",
  normalizePath(
    readme_txt,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")


cat(
  "- Original historical DESeq2 universe remains 12,400 retained features.\n"
)


cat(
  "- Submission Table S2 contains 12,393 biological targets.\n"
)


cat(
  "- Seven retained ERCC-prefixed features are excluded from submission only.\n"
)


cat(
  "- No DESeq2 statistic has been recalculated.\n"
)


cat(
  "- Primary DEG count remains 2,659.\n"
)


cat(
  "- Batch-adjusted DEG count remains 4,125.\n"
)


cat(
  "- Robust-core DEG count remains 1,796.\n"
)


cat(
  "- Table S3 requires no ERCC correction because no ERCC feature entered the robust core.\n"
)


cat("\nDone.\n")