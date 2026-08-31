################################################################################
# Script 155
# FINAL v3
#
# Supplementary Table S4
#
# Blood SRS and CTS molecular-endotype assignments and robustness summary
#
# Project:
#   Sepsis_DESeq2
#
# Manuscript:
#   Blood-only sepsis transcriptomic endotypes /
#   five-gene host-response signature
#
#
# TABLE S4 CONTENT
# ----------------
#
# 1. Primary SepstratifieR SRS assignments for all 45 blood samples
# 2. Integrated SRS + CTS assignments for the 35 sepsis blood samples
# 3. SRS distribution summary
# 4. CTS distribution summary
# 5. CTS x SRS cross-table
# 6. Frozen SRS/CTS robustness summary
# 7. Source/provenance manifest
#
#
# IMPORTANT
# ---------
#
# CTS source MUST be:
#
#   results/cts_consensus/tables/
#   107_CTS_predictions_BP_only.xlsx
#
# NOT:
#
#   107_CTS_predictions_BP_BC.xlsx
#
#
# Primary SRS source:
#
#   results/sepstratifier/publication_summary/
#   43_SRS_by_each_sample.xlsx
#
#
# EXPECTED PRIMARY RESULTS
# ------------------------
#
# SRS:
#
# Sepsis BP:
#   SRS1 = 28
#   SRS2 =  7
#   SRS3 =  0
#
# Controls BC:
#   SRS1 = 0
#   SRS2 = 3
#   SRS3 = 7
#
# Primary mNN outliers:
#
#   BP10
#   BP26
#   BP27
#   BP31
#
#
# CTS BP only:
#
#   CTS1 = 14
#   CTS2 =  6
#   CTS3 = 15
#
# NOTE:
#   The frozen source may encode these classes as:
#
#     1, 2, 3
#
#   rather than:
#
#     CTS1, CTS2, CTS3
#
#   Script 155 standardizes both forms to CTS1/CTS2/CTS3
#   without changing the original assignments.
#
#
# CTS x SRS:
#
#            SRS1  SRS2  SRS3
#   CTS1       14     0     0
#   CTS2        6     0     0
#   CTS3        8     7     0
#
#
# THIS SCRIPT DOES NOT:
#
#   - rerun SepstratifieR
#   - rerun CTS classifier
#   - rerun robustness analyses
#   - redefine any endotype
#
#
# CODING RULE
# -----------
#
# Never split [[ across lines.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 155 FINAL v3\n")
cat("Supplementary Table S4\n")
cat("Blood SRS and CTS assignments + robustness summary\n")
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
  "155_TableS4_SRS_CTS_endotypes"
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

normalize_sample_id <- function(x) {
  
  toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
}


as_logical_safe <- function(x) {
  
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
    "1"
  )
}


find_column <- function(
    data,
    exact_candidates,
    regex = NULL,
    label = "column",
    required = TRUE
) {
  
  nm <- names(data)
  
  nm_lower <- tolower(
    nm
  )
  
  
  for (
    candidate in exact_candidates
  ) {
    
    hit <- which(
      nm_lower ==
        tolower(
          candidate
        )
    )
    
    
    if (
      length(
        hit
      ) ==
      1
    ) {
      
      return(
        nm[hit]
      )
    }
  }
  
  
  if (!is.null(regex)) {
    
    hit <- grep(
      regex,
      nm,
      ignore.case = TRUE,
      perl = TRUE
    )
    
    
    if (
      length(
        hit
      ) ==
      1
    ) {
      
      return(
        nm[hit]
      )
    }
    
    
    if (
      length(
        hit
      ) >
      1
    ) {
      
      cat(
        "\nAmbiguous candidates for ",
        label,
        ":\n",
        sep = ""
      )
      
      print(
        nm[hit]
      )
    }
  }
  
  
  if (required) {
    
    cat(
      "\nAvailable columns for ",
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
  
  
  return(
    NA_character_
  )
}


find_best_sheet <- function(
    path,
    required_patterns
) {
  
  sheets <- readxl::excel_sheets(
    path
  )
  
  
  if (
    length(
      sheets
    ) ==
    0
  ) {
    
    stop(
      paste0(
        "Workbook contains no readable sheets:\n",
        path
      )
    )
  }
  
  
  scores <- rep(
    0L,
    length(
      sheets
    )
  )
  
  
  for (
    i in seq_along(
      sheets
    )
  ) {
    
    temp <- tryCatch(
      
      readxl::read_excel(
        path,
        sheet = sheets[i],
        n_max = 10
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    if (is.null(temp)) {
      next
    }
    
    
    nm <- names(
      temp
    )
    
    
    for (
      pattern in required_patterns
    ) {
      
      if (
        any(
          grepl(
            pattern,
            nm,
            ignore.case = TRUE,
            perl = TRUE
          )
        )
      ) {
        
        scores[i] <-
          scores[i] +
          1L
      }
    }
  }
  
  
  max_score <- max(
    scores
  )
  
  
  if (
    max_score ==
    0
  ) {
    
    cat(
      "\nWorkbook sheets:\n"
    )
    
    print(
      sheets
    )
    
    stop(
      paste0(
        "Could not identify usable sheet in:\n",
        path
      )
    )
  }
  
  
  best <- which(
    scores ==
      max_score
  )
  
  
  sheets[best[1]]
}


find_file_by_basename <- function(
    preferred_path,
    target_basename
) {
  
  if (
    file.exists(
      preferred_path
    )
  ) {
    
    return(
      preferred_path
    )
  }
  
  
  hits <- list.files(
    file.path(
      project_dir,
      "results"
    ),
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
    1
  ) {
    
    return(
      hits[1]
    )
  }
  
  
  if (
    length(
      hits
    ) ==
    0
  ) {
    
    stop(
      paste0(
        "Required source not found: ",
        target_basename
      )
    )
  }
  
  
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
      "Ambiguous source provenance for ",
      target_basename,
      "."
    )
  )
}


choose_bp_sample_column <- function(data) {
  
  candidate_names <- intersect(
    c(
      "sample_id",
      "Sample",
      "sample",
      "Sample_ID"
    ),
    names(
      data
    )
  )
  
  
  if (
    length(
      candidate_names
    ) ==
    0
  ) {
    
    candidate_names <- names(
      data
    )[
      grepl(
        "sample",
        names(
          data
        ),
        ignore.case = TRUE
      )
    ]
  }
  
  
  if (
    length(
      candidate_names
    ) ==
    0
  ) {
    
    stop(
      "No candidate CTS sample-ID column found."
    )
  }
  
  
  candidate_audit <- data.frame(
    
    column =
      candidate_names,
    
    n_nonmissing =
      NA_integer_,
    
    n_BP_format =
      NA_integer_,
    
    n_unique_BP =
      NA_integer_,
    
    stringsAsFactors = FALSE
  )
  
  
  for (
    i in seq_along(
      candidate_names
    )
  ) {
    
    x <- normalize_sample_id(
      data[[candidate_names[i]]]
    )
    
    
    valid <- grepl(
      "^BP[0-9]+$",
      x
    )
    
    
    candidate_audit$n_nonmissing[i] <-
      sum(
        !is.na(
          x
        ) &
          x !=
          ""
      )
    
    
    candidate_audit$n_BP_format[i] <-
      sum(
        valid
      )
    
    
    candidate_audit$n_unique_BP[i] <-
      length(
        unique(
          x[valid]
        )
      )
  }
  
  
  cat(
    "\nCTS sample-ID candidate audit:\n"
  )
  
  print(
    candidate_audit,
    row.names = FALSE
  )
  
  
  perfect <- candidate_audit$column[
    candidate_audit$n_BP_format ==
      35 &
      candidate_audit$n_unique_BP ==
      35
  ]
  
  
  if (
    length(
      perfect
    ) ==
    1
  ) {
    
    return(
      list(
        column = perfect[1],
        audit = candidate_audit
      )
    )
  }
  
  
  if (
    length(
      perfect
    ) >
    1
  ) {
    
    reference <- normalize_sample_id(
      data[[perfect[1]]]
    )
    
    
    all_identical <- all(
      vapply(
        perfect[-1],
        function(col_name) {
          
          identical(
            reference,
            normalize_sample_id(
              data[[col_name]]
            )
          )
        },
        logical(1)
      )
    )
    
    
    if (all_identical) {
      
      preferred_order <- c(
        "sample_id",
        "Sample",
        "sample",
        "Sample_ID"
      )
      
      
      chosen <- preferred_order[
        preferred_order %in%
          perfect
      ][1]
      
      
      cat(
        "Multiple equivalent BP sample-ID columns found; using ",
        chosen,
        ".\n",
        sep = ""
      )
      
      
      return(
        list(
          column = chosen,
          audit = candidate_audit
        )
      )
    }
    
    
    stop(
      paste0(
        "Multiple non-identical CTS columns contain 35 valid BP sample IDs: ",
        paste(
          perfect,
          collapse = ", "
        )
      )
    )
  }
  
  
  best_index <- which.max(
    candidate_audit$n_unique_BP
  )
  
  
  stop(
    paste0(
      "No CTS sample-ID column contained all 35 BP identifiers. ",
      "Best candidate: ",
      candidate_audit$column[best_index],
      " with ",
      candidate_audit$n_unique_BP[best_index],
      " unique BP IDs."
    )
  )
}


standardize_cts_label <- function(x) {
  
  x2 <- toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
  
  
  x2 <- gsub(
    "\\s+",
    "",
    x2
  )
  
  
  dplyr::case_when(
    
    x2 %in%
      c(
        "1",
        "CTS1",
        "CTS_1",
        "CTS-1"
      ) ~
      "CTS1",
    
    x2 %in%
      c(
        "2",
        "CTS2",
        "CTS_2",
        "CTS-2"
      ) ~
      "CTS2",
    
    x2 %in%
      c(
        "3",
        "CTS3",
        "CTS_3",
        "CTS-3"
      ) ~
      "CTS3",
    
    TRUE ~
      NA_character_
  )
}


# =============================================================================
# 5. PRIMARY INPUT FILES
# =============================================================================

srs_file <- find_file_by_basename(
  
  preferred_path = file.path(
    project_dir,
    "results",
    "sepstratifier",
    "publication_summary",
    "43_SRS_by_each_sample.xlsx"
  ),
  
  target_basename =
    "43_SRS_by_each_sample.xlsx"
)


cts_file <- find_file_by_basename(
  
  preferred_path = file.path(
    project_dir,
    "results",
    "cts_consensus",
    "tables",
    "107_CTS_predictions_BP_only.xlsx"
  ),
  
  target_basename =
    "107_CTS_predictions_BP_only.xlsx"
)


sample_metadata_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "151_Table1_TableS1_discovery_cohort",
  "tables",
  "151_TableS1_blood_sample_metadata.csv"
)


participant_metadata_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "151_Table1_TableS1_discovery_cohort",
  "tables",
  "151_TableS1_participant_metadata.csv"
)


required_files <- c(
  srs_file,
  cts_file,
  sample_metadata_file,
  participant_metadata_file
)


missing_files <- required_files[
  !file.exists(
    required_files
  )
]


if (
  length(
    missing_files
  ) >
  0
) {
  
  cat(
    "\nMissing required file(s):\n"
  )
  
  print(
    missing_files
  )
  
  stop(
    "Required frozen source file(s) missing."
  )
}


cat("\nPRIMARY SOURCE FILES\n")
cat("--------------------\n")


cat(
  "SRS:\n  ",
  normalizePath(
    srs_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "CTS BP-only:\n  ",
  normalizePath(
    cts_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


# =============================================================================
# 6. READ PRIMARY SRS ASSIGNMENTS
# =============================================================================

cat("SRS workbook sheets:\n")

print(
  readxl::excel_sheets(
    srs_file
  )
)


srs_sheet <- find_best_sheet(
  srs_file,
  required_patterns = c(
    "sample|barcode",
    "SRS",
    "SRSq"
  )
)


cat(
  "Selected SRS sheet = ",
  srs_sheet,
  "\n",
  sep = ""
)


srs_raw <- readxl::read_excel(
  srs_file,
  sheet = srs_sheet
) %>%
  
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "SRS source dimensions = ",
  nrow(
    srs_raw
  ),
  " x ",
  ncol(
    srs_raw
  ),
  "\n",
  sep = ""
)


cat("SRS source columns:\n")

print(
  names(
    srs_raw
  )
)


srs_sample_col <- find_column(
  
  srs_raw,
  
  exact_candidates = c(
    "sample_id",
    "Sample",
    "sample",
    "barcode_id",
    "barcode"
  ),
  
  regex = "sample.*id|^sample$|barcode",
  
  label = "SRS sample ID"
)


srs_class_col <- find_column(
  
  srs_raw,
  
  exact_candidates = c(
    "SRS",
    "srs"
  ),
  
  regex = "^SRS$",
  
  label = "SRS class"
)


srsq_col <- find_column(
  
  srs_raw,
  
  exact_candidates = c(
    "SRSq",
    "srsq"
  ),
  
  regex = "^SRSq$",
  
  label = "SRSq"
)


mnn_col <- find_column(
  
  srs_raw,
  
  exact_candidates = c(
    "mNN_outlier",
    "mnn_outlier",
    "mNN",
    "outlier"
  ),
  
  regex = "mNN.*outlier|outlier.*mNN|^mNN$|^outlier$",
  
  label = "mNN outlier",
  
  required = FALSE
)


cat("\nDetected SRS columns:\n")

cat(
  "Sample ID = ",
  srs_sample_col,
  "\n",
  sep = ""
)

cat(
  "SRS = ",
  srs_class_col,
  "\n",
  sep = ""
)

cat(
  "SRSq = ",
  srsq_col,
  "\n",
  sep = ""
)

cat(
  "mNN outlier = ",
  ifelse(
    is.na(
      mnn_col
    ),
    "not found in source",
    mnn_col
  ),
  "\n",
  sep = ""
)


srs <- data.frame(
  
  sample_id =
    normalize_sample_id(
      srs_raw[[srs_sample_col]]
    ),
  
  SRS =
    as.character(
      srs_raw[[srs_class_col]]
    ),
  
  SRSq =
    suppressWarnings(
      as.numeric(
        srs_raw[[srsq_col]]
      )
    ),
  
  stringsAsFactors = FALSE
)


if (
  !is.na(
    mnn_col
  )
) {
  
  srs$mNN_outlier <-
    as_logical_safe(
      srs_raw[[mnn_col]]
    )
  
} else {
  
  srs$mNN_outlier <-
    FALSE
}


srs <- srs %>%
  
  dplyr::filter(
    grepl(
      "^(BP|BC)[0-9]+$",
      sample_id
    )
  ) %>%
  
  dplyr::distinct(
    sample_id,
    .keep_all = TRUE
  )


if (
  nrow(
    srs
  ) !=
  45
) {
  
  stop(
    paste0(
      "Expected 45 unique blood SRS assignments; observed ",
      nrow(
        srs
      )
    )
  )
}


if (
  any(
    is.na(
      srs$SRS
    )
  )
) {
  
  stop(
    "Missing SRS categorical assignments detected."
  )
}


if (
  any(
    !is.finite(
      srs$SRSq
    )
  )
) {
  
  stop(
    "Missing/non-finite SRSq values detected."
  )
}


# =============================================================================
# 7. READ FROZEN BLOOD SAMPLE METADATA
# =============================================================================

sample_meta <- read.csv(
  sample_metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


sample_meta$sample_id <- normalize_sample_id(
  sample_meta$sample_id
)


sample_meta <- sample_meta %>%
  
  dplyr::filter(
    biofluid ==
      "Whole blood"
  )


if (
  nrow(
    sample_meta
  ) !=
  45
) {
  
  stop(
    paste0(
      "Expected 45 frozen whole-blood metadata rows; observed ",
      nrow(
        sample_meta
      )
    )
  )
}


if (
  dplyr::n_distinct(
    sample_meta$sample_id
  ) !=
  45
) {
  
  stop(
    "Duplicate blood sample IDs in Script 151 metadata."
  )
}


participant_meta <- read.csv(
  participant_metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# =============================================================================
# 8. MERGE PRIMARY SRS WITH SAMPLE METADATA
# =============================================================================

srs_all_blood <- sample_meta %>%
  
  dplyr::select(
    sample_id,
    participant_id,
    cohort,
    sample_group,
    sequencing_batch
  ) %>%
  
  dplyr::left_join(
    srs,
    by = "sample_id"
  )


if (
  any(
    is.na(
      srs_all_blood$SRS
    )
  )
) {
  
  cat(
    "\nSamples without SRS assignment:\n"
  )
  
  print(
    srs_all_blood$sample_id[
      is.na(
        srs_all_blood$SRS
      )
    ]
  )
  
  
  stop(
    "Blood metadata/SRS merge failed."
  )
}


# =============================================================================
# 9. SRS PRIMARY AUDIT
# =============================================================================

srs_summary <- srs_all_blood %>%
  
  dplyr::count(
    sample_group,
    SRS,
    name = "n"
  ) %>%
  
  dplyr::group_by(
    sample_group
  ) %>%
  
  dplyr::mutate(
    percent =
      100 *
      n /
      sum(
        n
      )
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::arrange(
    sample_group,
    SRS
  )


cat("\nSRS DISTRIBUTION\n")
cat("----------------\n")


print(
  srs_summary,
  row.names = FALSE
)


expected_srs <- data.frame(
  
  sample_group = c(
    "BC",
    "BC",
    "BP",
    "BP"
  ),
  
  SRS = c(
    "SRS2",
    "SRS3",
    "SRS1",
    "SRS2"
  ),
  
  n = c(
    3,
    7,
    28,
    7
  ),
  
  stringsAsFactors = FALSE
)


observed_srs <- srs_summary %>%
  
  dplyr::select(
    sample_group,
    SRS,
    n
  ) %>%
  
  dplyr::mutate(
    n =
      as.numeric(
        n
      )
  )


srs_comparison <- dplyr::full_join(
  
  expected_srs %>%
    dplyr::rename(
      expected_n = n
    ),
  
  observed_srs %>%
    dplyr::rename(
      observed_n = n
    ),
  
  by = c(
    "sample_group",
    "SRS"
  )
) %>%
  
  dplyr::mutate(
    
    expected_n =
      dplyr::coalesce(
        as.numeric(
          expected_n
        ),
        0
      ),
    
    observed_n =
      dplyr::coalesce(
        as.numeric(
          observed_n
        ),
        0
      ),
    
    match =
      expected_n ==
      observed_n
  ) %>%
  
  dplyr::arrange(
    sample_group,
    SRS
  )


if (
  !all(
    srs_comparison$match
  )
) {
  
  print(
    srs_comparison
  )
  
  stop(
    "Primary SRS distribution audit failed."
  )
}


# =============================================================================
# 10. mNN OUTLIER AUDIT
# =============================================================================

expected_outliers <- sort(
  c(
    "BP10",
    "BP26",
    "BP27",
    "BP31"
  )
)


if (
  !is.na(
    mnn_col
  )
) {
  
  observed_outliers <- sort(
    srs_all_blood$sample_id[
      srs_all_blood$mNN_outlier
    ]
  )
  
  
  if (
    !identical(
      observed_outliers,
      expected_outliers
    )
  ) {
    
    cat(
      "\nObserved mNN outliers:\n"
    )
    
    print(
      observed_outliers
    )
    
    
    cat(
      "\nExpected frozen mNN outliers:\n"
    )
    
    print(
      expected_outliers
    )
    
    
    stop(
      "mNN outlier audit failed."
    )
  }
  
} else {
  
  srs_all_blood$mNN_outlier <-
    srs_all_blood$sample_id %in%
    expected_outliers
  
  
  observed_outliers <-
    expected_outliers
}


cat(
  "\nPrimary SRS distribution and mNN audit passed.\n"
)


# =============================================================================
# 11. READ PRIMARY CTS BP-ONLY ASSIGNMENTS
# =============================================================================

cat("\nCTS workbook sheets:\n")

print(
  readxl::excel_sheets(
    cts_file
  )
)


cts_sheet <- find_best_sheet(
  cts_file,
  required_patterns = c(
    "sample",
    "CTS|subtype"
  )
)


cat(
  "Selected CTS sheet = ",
  cts_sheet,
  "\n",
  sep = ""
)


cts_raw <- readxl::read_excel(
  cts_file,
  sheet = cts_sheet
) %>%
  
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "CTS source dimensions = ",
  nrow(
    cts_raw
  ),
  " x ",
  ncol(
    cts_raw
  ),
  "\n",
  sep = ""
)


cat("CTS source columns:\n")

print(
  names(
    cts_raw
  )
)


# =============================================================================
# 12. CHOOSE CTS SAMPLE-ID COLUMN BY CONTENT
# =============================================================================

cts_sample_selection <- choose_bp_sample_column(
  cts_raw
)


cts_sample_col <-
  cts_sample_selection$column


cts_sample_column_audit <-
  cts_sample_selection$audit


cts_class_col <- find_column(
  
  cts_raw,
  
  exact_candidates = c(
    "CTS",
    "cts",
    "Subtype",
    "subtype",
    "Prediction",
    "prediction"
  ),
  
  regex = "^CTS$|subtype|prediction",
  
  label = "CTS class"
)


cat("\nDetected CTS columns:\n")

cat(
  "Sample ID = ",
  cts_sample_col,
  "\n",
  sep = ""
)

cat(
  "CTS raw class = ",
  cts_class_col,
  "\n",
  sep = ""
)


# =============================================================================
# 13. STANDARDIZE CTS LABELS
# =============================================================================

cts_raw_labels <- as.character(
  cts_raw[[cts_class_col]]
)


cts_standardized_labels <- standardize_cts_label(
  cts_raw_labels
)


cts_label_audit <- data.frame(
  
  raw_label =
    sort(
      unique(
        as.character(
          cts_raw_labels
        )
      )
    ),
  
  stringsAsFactors = FALSE
)


cts_label_audit$standardized_label <-
  standardize_cts_label(
    cts_label_audit$raw_label
  )


cat("\nCTS label standardization audit:\n")

print(
  cts_label_audit,
  row.names = FALSE
)


if (
  any(
    is.na(
      cts_standardized_labels
    )
  )
) {
  
  bad_labels <- sort(
    unique(
      cts_raw_labels[
        is.na(
          cts_standardized_labels
        )
      ]
    )
  )
  
  
  cat(
    "\nUnrecognized raw CTS label(s):\n"
  )
  
  print(
    bad_labels
  )
  
  
  stop(
    "CTS label standardization failed."
  )
}


cts <- data.frame(
  
  sample_id =
    normalize_sample_id(
      cts_raw[[cts_sample_col]]
    ),
  
  CTS_raw =
    as.character(
      cts_raw[[cts_class_col]]
    ),
  
  CTS =
    cts_standardized_labels,
  
  stringsAsFactors = FALSE
) %>%
  
  dplyr::filter(
    grepl(
      "^BP[0-9]+$",
      sample_id
    )
  ) %>%
  
  dplyr::distinct(
    sample_id,
    .keep_all = TRUE
  )


if (
  nrow(
    cts
  ) !=
  35
) {
  
  cat(
    "\nCTS BP sample IDs found:\n"
  )
  
  print(
    cts$sample_id
  )
  
  
  stop(
    paste0(
      "Expected 35 BP-only CTS assignments; observed ",
      nrow(
        cts
      )
    )
  )
}


if (
  any(
    !cts$CTS %in%
    c(
      "CTS1",
      "CTS2",
      "CTS3"
    )
  )
) {
  
  stop(
    "Unexpected standardized CTS label(s)."
  )
}


# =============================================================================
# 14. CTS DISTRIBUTION AUDIT
# =============================================================================

cts_summary <- cts %>%
  
  dplyr::count(
    CTS,
    name = "n"
  ) %>%
  
  dplyr::mutate(
    n =
      as.numeric(
        n
      ),
    percent =
      100 *
      n /
      sum(
        n
      )
  ) %>%
  
  dplyr::arrange(
    CTS
  )


cat("\nCTS DISTRIBUTION\n")
cat("----------------\n")


print(
  cts_summary,
  row.names = FALSE
)


expected_cts <- data.frame(
  
  CTS = c(
    "CTS1",
    "CTS2",
    "CTS3"
  ),
  
  expected_n = c(
    14,
    6,
    15
  ),
  
  stringsAsFactors = FALSE
)


cts_comparison <- dplyr::full_join(
  
  expected_cts,
  
  cts_summary %>%
    dplyr::select(
      CTS,
      observed_n = n
    ),
  
  by = "CTS"
) %>%
  
  dplyr::mutate(
    
    expected_n =
      dplyr::coalesce(
        as.numeric(
          expected_n
        ),
        0
      ),
    
    observed_n =
      dplyr::coalesce(
        as.numeric(
          observed_n
        ),
        0
      ),
    
    match =
      expected_n ==
      observed_n
  ) %>%
  
  dplyr::arrange(
    CTS
  )


if (
  !all(
    cts_comparison$match
  )
) {
  
  cat(
    "\nCTS expected-versus-observed comparison:\n"
  )
  
  print(
    cts_comparison
  )
  
  
  stop(
    "Primary BP-only CTS distribution audit failed."
  )
}


cat(
  "Primary CTS distribution audit passed.\n"
)


# =============================================================================
# 15. BUILD INTEGRATED SEPSIS SRS + CTS TABLE
# =============================================================================

sepsis_srs_cts <- srs_all_blood %>%
  
  dplyr::filter(
    sample_group ==
      "BP"
  ) %>%
  
  dplyr::left_join(
    cts %>%
      dplyr::select(
        sample_id,
        CTS
      ),
    by = "sample_id"
  )


if (
  nrow(
    sepsis_srs_cts
  ) !=
  35
) {
  
  stop(
    "Integrated BP SRS/CTS table does not contain 35 samples."
  )
}


if (
  any(
    is.na(
      sepsis_srs_cts$CTS
    )
  )
) {
  
  cat(
    "\nBP samples without CTS assignment:\n"
  )
  
  print(
    sepsis_srs_cts$sample_id[
      is.na(
        sepsis_srs_cts$CTS
      )
    ]
  )
  
  
  stop(
    "Missing CTS assignment after joining BP SRS and CTS."
  )
}


# =============================================================================
# 16. CTS x SRS CROSS-TABLE
# =============================================================================

cts_srs_matrix <- table(
  
  factor(
    sepsis_srs_cts$CTS,
    levels = c(
      "CTS1",
      "CTS2",
      "CTS3"
    )
  ),
  
  factor(
    sepsis_srs_cts$SRS,
    levels = c(
      "SRS1",
      "SRS2",
      "SRS3"
    )
  )
)


cts_srs_cross <- data.frame(
  
  CTS =
    rownames(
      cts_srs_matrix
    ),
  
  as.data.frame.matrix(
    cts_srs_matrix
  ),
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
)


cat("\nCTS x SRS\n")
cat("---------\n")


print(
  cts_srs_cross,
  row.names = FALSE
)


expected_cross <- matrix(
  
  c(
    14, 0, 0,
    6, 0, 0,
    8, 7, 0
  ),
  
  nrow = 3,
  
  byrow = TRUE
)


observed_cross <- unname(
  as.matrix(
    cts_srs_matrix
  )
)


if (
  !all(
    dim(
      observed_cross
    ) ==
    dim(
      expected_cross
    )
  ) ||
  !all(
    as.numeric(
      observed_cross
    ) ==
    as.numeric(
      expected_cross
    )
  )
) {
  
  cat(
    "\nObserved cross-table matrix:\n"
  )
  
  print(
    observed_cross
  )
  
  
  cat(
    "\nExpected cross-table matrix:\n"
  )
  
  print(
    expected_cross
  )
  
  
  stop(
    "CTS x SRS cross-table audit failed."
  )
}


cat(
  "CTS x SRS audit passed successfully.\n"
)


# =============================================================================
# 17. PRIMARY SAMPLE-LEVEL FINAL TABLES
# =============================================================================

sepsis_srs_cts <- sepsis_srs_cts %>%
  
  dplyr::select(
    participant_id,
    sample_id,
    sequencing_batch,
    SRS,
    SRSq,
    mNN_outlier,
    CTS
  ) %>%
  
  dplyr::arrange(
    CTS,
    dplyr::desc(
      SRSq
    )
  )


srs_all_blood <- srs_all_blood %>%
  
  dplyr::select(
    participant_id,
    sample_id,
    sample_group,
    sequencing_batch,
    SRS,
    SRSq,
    mNN_outlier
  ) %>%
  
  dplyr::arrange(
    sample_group,
    SRS,
    dplyr::desc(
      SRSq
    )
  )


# =============================================================================
# 18. FROZEN ROBUSTNESS SUMMARY
# =============================================================================

robustness_summary <- data.frame(
  
  Framework = c(
    
    "SRS",
    "SRS",
    "SRS",
    "SRS",
    "SRS",
    
    "CTS",
    "CTS",
    "CTS",
    "CTS"
  ),
  
  Sensitivity_analysis = c(
    
    "Exact primary rerun",
    
    "Alternative mNN k values",
    
    "BP-only rerun",
    
    "BP-only recommended k range",
    
    "BP-only batch-adjusted sensitivity",
    
    "50-seed reproducibility",
    
    "BP-only VST preprocessing",
    
    "BP-only batch-adjusted preprocessing",
    
    "CTS x sequencing batch"
  ),
  
  Result = c(
    
    paste0(
      "45/45 categorical assignments reproduced; maximum absolute ",
      "SRSq difference = 4.44e-16; mNN outliers identical."
    ),
    
    paste0(
      "Across tested k values 5, 9, 10, 11, 12, 14, 15, 20 and 25, ",
      "categorical agreement was 97.8-100%; minimum SRSq Spearman ",
      "correlation across tested k values was 0.938."
    ),
    
    paste0(
      "Original all-blood versus BP-only k=20: 33/35 assignments ",
      "unchanged (94.3%); SRSq Spearman rho = 0.948."
    ),
    
    paste0(
      "Across the recommended BP-only k range 7-10, 34/35 samples ",
      "retained the primary SRS assignment."
    ),
    
    paste0(
      "Original BP versus BP-only batch-adjusted k=20: 30/35 ",
      "assignments unchanged (85.7%); SRSq Spearman rho = 0.804."
    ),
    
    paste0(
      "Primary CTS modal assignment matched the original classifier ",
      "for 35/35 samples; all 35 samples had >=95% seed stability, ",
      "34/35 were 100% stable, and minimum individual stability was 98%."
    ),
    
    paste0(
      "BP-only VST sensitivity: modal CTS assignment matched the ",
      "primary assignment for 35/35 samples."
    ),
    
    paste0(
      "BP-only batch-adjusted sensitivity: modal CTS assignment ",
      "matched the primary assignment for 29/35 samples (82.9%). ",
      "BP30 and BP7 showed the lowest modal support (54% and 56%)."
    ),
    
    paste0(
      "Monte Carlo Fisher exact P = 0.102049; ",
      "Cramer's V = 0.4536; n = 35."
    )
  ),
  
  Interpretation = c(
    
    "Exact computational reproducibility of the primary SRS analysis.",
    
    paste0(
      "Primary SRS assignments and continuous SRSq were highly stable ",
      "to reasonable mNN neighborhood size."
    ),
    
    paste0(
      "Primary SRS structure was largely preserved when healthy controls ",
      "were omitted from the projection."
    ),
    
    paste0(
      "Categorical SRS assignment was highly stable in the recommended ",
      "BP-only neighborhood range."
    ),
    
    paste0(
      "Batch-adjusted projection produced greater categorical movement ",
      "and was treated as a conservative sensitivity analysis."
    ),
    
    "CTS classifier assignments were highly reproducible to random seed.",
    
    paste0(
      "CTS assignments were invariant to the BP-only VST ",
      "sensitivity preprocessing."
    ),
    
    paste0(
      "Batch-adjusted preprocessing altered a subset of CTS assignments; ",
      "primary CTS labels were therefore retained."
    ),
    
    paste0(
      "No statistically significant CTS-by-batch association was detected, ",
      "but moderate structure was present and was not interpreted as ",
      "evidence of batch independence."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 19. ROBUSTNESS PROVENANCE FILES
# =============================================================================

all_result_files <- list.files(
  file.path(
    project_dir,
    "results"
  ),
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)


all_script_files <- list.files(
  file.path(
    project_dir,
    "scripts"
  ),
  recursive = FALSE,
  full.names = TRUE,
  include.dirs = FALSE
)


robustness_candidates <- c(
  all_result_files,
  all_script_files
)


robustness_candidates <- robustness_candidates[
  grepl(
    "(40b|40c|107d|146)",
    basename(
      robustness_candidates
    ),
    ignore.case = TRUE
  )
]


robustness_candidates <- sort(
  unique(
    robustness_candidates
  )
)


if (
  length(
    robustness_candidates
  ) >
  0
) {
  
  robustness_file_manifest <- data.frame(
    
    path =
      normalizePath(
        robustness_candidates,
        winslash = "\\",
        mustWork = TRUE
      ),
    
    stringsAsFactors = FALSE
  )
  
} else {
  
  robustness_file_manifest <- data.frame(
    
    path =
      "No matching 40b/40c/107d/146 filenames detected automatically",
    
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 20. SOURCE MANIFEST
# =============================================================================

source_manifest <- data.frame(
  
  Component = c(
    "Primary SRS assignments",
    "Primary BP-only CTS assignments",
    "Frozen blood sample metadata",
    "Frozen participant metadata",
    "Robustness summary"
  ),
  
  Source = c(
    
    normalizePath(
      srs_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      cts_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      sample_metadata_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      participant_metadata_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    paste0(
      "Previously completed Scripts 40b, 40c, 107d and ",
      "Supplementary Figure S3 / Script 146; no robustness ",
      "analysis rerun in Script 155."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 21. TABLE S4 README
# =============================================================================

s4_readme <- data.frame(
  
  Item = c(
    "Title",
    "Primary SRS source",
    "Primary CTS source",
    "CTS source label encoding",
    "SRS population",
    "CTS population",
    "Primary SRS result",
    "Primary CTS result",
    "CTS-SRS relationship",
    "SRSq",
    "mNN outliers",
    "Robustness",
    "Important limitation"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S4. Blood SRS and Consensus Transcriptomic ",
      "Subtype assignments and robustness analyses."
    ),
    
    paste0(
      "Final SepstratifieR blood publication summary: ",
      basename(
        srs_file
      ),
      "."
    ),
    
    paste0(
      "ConsensusTranscriptomicSubtype BP-only predictions: ",
      basename(
        cts_file
      ),
      ". The BP+BC CTS output is not used."
    ),
    
    paste0(
      "The frozen BP-only source encodes CTS classes numerically as ",
      "1, 2 and 3; these labels are standardized for reporting as ",
      "CTS1, CTS2 and CTS3 without reclassification."
    ),
    
    paste0(
      "All 45 blood transcriptomes: 35 sepsis BP and ",
      "10 healthy-control BC samples."
    ),
    
    "CTS classification is reported only for the 35 sepsis blood transcriptomes.",
    
    paste0(
      "Sepsis: SRS1 28/35, SRS2 7/35, SRS3 0/35. ",
      "Controls: SRS1 0/10, SRS2 3/10, SRS3 7/10."
    ),
    
    "CTS1 14/35, CTS2 6/35, CTS3 15/35.",
    
    paste0(
      "CTS1 and CTS2 were entirely SRS1; CTS3 comprised ",
      "8 SRS1 and 7 SRS2 samples."
    ),
    
    paste0(
      "SRSq is reported as the continuous quantitative output of ",
      "SepstratifieR and is not interpreted as a probability."
    ),
    
    "Primary mNN outliers: BP10, BP26, BP27 and BP31.",
    
    paste0(
      "The Robustness_summary sheet reports frozen sensitivity results ",
      "from prior analyses; Script 155 does not rerun classifiers."
    ),
    
    paste0(
      "Batch-adjusted SRS/CTS analyses are sensitivity analyses because ",
      "sequencing batch is partially structured by biological group; ",
      "they are not treated as superior replacement classifications."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 22. OUTPUT FILES
# =============================================================================

submission_file <- file.path(
  tables_dir,
  "155_TableS4_SRS_CTS_assignments_and_robustness.xlsx"
)


audit_file <- file.path(
  audit_dir,
  "155_INTERNAL_AUDIT_TableS4_SRS_CTS.xlsx"
)


note_file <- file.path(
  text_dir,
  "155_TableS4_title_and_note_EN.txt"
)


# =============================================================================
# 23. EXCEL STYLES
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
# 24. WRITE SUBMISSION WORKBOOK
# =============================================================================

wb <- openxlsx::createWorkbook()


submission_objects <- list(
  
  S4_ReadMe =
    s4_readme,
  
  SRS_all_blood =
    srs_all_blood,
  
  Sepsis_SRS_CTS =
    sepsis_srs_cts,
  
  SRS_summary =
    srs_summary,
  
  CTS_summary =
    cts_summary,
  
  CTSxSRS =
    cts_srs_cross,
  
  Robustness_summary =
    robustness_summary
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
      "S4_ReadMe"
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
  "S4_ReadMe",
  cols = 1,
  widths = 30
)


openxlsx::setColWidths(
  wb,
  "S4_ReadMe",
  cols = 2,
  widths = 90
)


openxlsx::setColWidths(
  wb,
  "Robustness_summary",
  cols = 3:4,
  widths = 65
)


openxlsx::saveWorkbook(
  wb,
  submission_file,
  overwrite = TRUE
)


# =============================================================================
# 25. INTERNAL AUDIT WORKBOOK
# =============================================================================

wb_audit <- openxlsx::createWorkbook()


audit_objects <- list(
  
  Source_manifest =
    source_manifest,
  
  Robustness_files =
    robustness_file_manifest,
  
  SRS_comparison =
    srs_comparison,
  
  CTS_sample_columns =
    cts_sample_column_audit,
  
  CTS_label_mapping =
    cts_label_audit,
  
  CTS_comparison =
    cts_comparison,
  
  CTSxSRS =
    cts_srs_cross,
  
  SRS_raw =
    srs_raw,
  
  CTS_raw =
    cts_raw
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
  audit_file,
  overwrite = TRUE
)


# =============================================================================
# 26. TABLE TITLE / NOTE
# =============================================================================

table_note <- c(
  
  paste0(
    "Supplementary Table S4. Blood SRS and Consensus Transcriptomic ",
    "Subtype assignments and robustness analyses."
  ),
  
  "",
  
  paste0(
    "Primary SepstratifieR SRS assignments are reported for all 45 blood ",
    "transcriptomes, whereas CTS assignments are restricted to the 35 ",
    "sepsis blood transcriptomes. The frozen BP-only CTS source encoded ",
    "the subtype labels numerically as 1, 2 and 3; these were standardized ",
    "for reporting as CTS1, CTS2 and CTS3 without changing any assignment. ",
    "SRSq denotes the continuous SepstratifieR output and is not interpreted ",
    "as a probability. mNN outlier status is retained as a quality-control ",
    "annotation. Robustness analyses evaluate sensitivity to neighborhood ",
    "size, sample composition, preprocessing and sequencing-batch adjustment; ",
    "primary endotype assignments were not replaced by sensitivity ",
    "classifications."
  )
)


writeLines(
  table_note,
  note_file
)


# =============================================================================
# 27. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "155_sessionInfo.txt"
  )
)


# =============================================================================
# 28. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 155 FINAL v3 completed successfully.\n")
cat("====================================================================\n\n")


cat("PRIMARY SRS\n")
cat("-----------\n")


cat(
  "Blood SRS samples = ",
  nrow(
    srs_all_blood
  ),
  "\n",
  sep = ""
)


cat(
  "Sepsis BP: SRS1 = ",
  sum(
    srs_all_blood$sample_group ==
      "BP" &
      srs_all_blood$SRS ==
      "SRS1"
  ),
  "; SRS2 = ",
  sum(
    srs_all_blood$sample_group ==
      "BP" &
      srs_all_blood$SRS ==
      "SRS2"
  ),
  "; SRS3 = ",
  sum(
    srs_all_blood$sample_group ==
      "BP" &
      srs_all_blood$SRS ==
      "SRS3"
  ),
  "\n",
  sep = ""
)


cat(
  "Controls BC: SRS1 = ",
  sum(
    srs_all_blood$sample_group ==
      "BC" &
      srs_all_blood$SRS ==
      "SRS1"
  ),
  "; SRS2 = ",
  sum(
    srs_all_blood$sample_group ==
      "BC" &
      srs_all_blood$SRS ==
      "SRS2"
  ),
  "; SRS3 = ",
  sum(
    srs_all_blood$sample_group ==
      "BC" &
      srs_all_blood$SRS ==
      "SRS3"
  ),
  "\n",
  sep = ""
)


cat(
  "mNN outliers = ",
  paste(
    sort(
      srs_all_blood$sample_id[
        srs_all_blood$mNN_outlier
      ]
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat("\nPRIMARY CTS — BP ONLY\n")
cat("---------------------\n")


cat(
  "CTS sample-ID source column = ",
  cts_sample_col,
  "\n",
  sep = ""
)


cat(
  "Raw CTS labels = ",
  paste(
    sort(
      unique(
        cts$CTS_raw
      )
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "Standardized CTS labels = ",
  paste(
    sort(
      unique(
        cts$CTS
      )
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "CTS samples = ",
  nrow(
    cts
  ),
  "\n",
  sep = ""
)


cat(
  "CTS1 = ",
  sum(
    cts$CTS ==
      "CTS1"
  ),
  "\n",
  sep = ""
)


cat(
  "CTS2 = ",
  sum(
    cts$CTS ==
      "CTS2"
  ),
  "\n",
  sep = ""
)


cat(
  "CTS3 = ",
  sum(
    cts$CTS ==
      "CTS3"
  ),
  "\n",
  sep = ""
)


cat("\nCTS x SRS\n")
cat("---------\n")


print(
  cts_srs_cross,
  row.names = FALSE
)


cat("\nROBUSTNESS SUMMARY\n")
cat("------------------\n")


print(
  robustness_summary[
    ,
    c(
      "Framework",
      "Sensitivity_analysis",
      "Result"
    )
  ],
  row.names = FALSE
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Supplementary Table S4:\n  ",
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
  "- Primary SRS source is the final Script 43 blood publication summary.\n"
)


cat(
  "- Primary CTS source is 107_CTS_predictions_BP_only.xlsx.\n"
)


cat(
  "- Frozen numeric CTS labels 1/2/3 are standardized only for reporting as CTS1/CTS2/CTS3.\n"
)


cat(
  "- No CTS sample is reclassified by Script 155.\n"
)


cat(
  "- Do not use 107_CTS_predictions_BP_BC.xlsx for manuscript CTS assignments.\n"
)


cat(
  "- SRSq is a continuous classifier output, not a probability.\n"
)


cat(
  "- CTS is reported for sepsis blood only.\n"
)


cat(
  "- Batch-adjusted SRS/CTS classifications remain sensitivity analyses.\n"
)


cat(
  "- Script 155 does not rerun SepstratifieR or CTS classifiers.\n"
)


cat("\nDone.\n")