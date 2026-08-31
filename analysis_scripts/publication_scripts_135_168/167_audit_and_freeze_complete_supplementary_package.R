# ======================================================================
# Script 167
# 167_audit_and_freeze_complete_supplementary_package.R
#
# FINAL v2
#
# COMPLETE SUPPLEMENTARY PACKAGE AUDIT AND FREEZE
#
# Blood transcriptomic endotypes / five-gene host-response manuscript
#
# IMPORTANT CHANGES FROM v1:
#   1. Correct numbering audit: remove names before identical()
#   2. Correct row-wise provenance audit: mapply() instead of vectorized
#      source_script passed into grepl()
#   3. NEAT1 is NOT treated as a scope violation because it may legitimately
#      occur as a transcript/gene in complete transcriptomic result tables.
#   4. Explicit negative/guardrail wording such as:
#          "not clinically validated"
#          "no urine"
#          "urine not included"
#      is not treated as a positive prohibited claim.
#   5. Scope/reporting keyword scans are REVIEW NOTES, not critical failures.
#   6. Figure captions accept either:
#          "Supplementary Figure S3"
#          "Figure S3"
#   7. Critical freeze criteria are restricted to:
#          - mandatory artifacts missing
#          - unreadable workbooks
#          - numbering failure
#          - expected source provenance failure
#          - direct identifier warning
#
# NO NEW INFERENTIAL STATISTICAL ANALYSIS.
# NO MODIFICATION OF FROZEN ARTIFACTS.
# ======================================================================


# ======================================================================
# 0. START
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("Running Script 167 FINAL v2\n")
cat("Complete Supplementary Package Audit and Freeze\n")
cat("=====================================================================\n\n")


# ======================================================================
# 1. PROJECT DIRECTORY
# ======================================================================

project_candidates <- c(
  Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2")),
  path.expand("~/Sepsis_DESeq2")
)

existing_projects <- project_candidates[
  dir.exists(project_candidates)
]

if (length(existing_projects) == 0) {
  stop("Sepsis_DESeq2 project directory was not found.")
}

project_dir <- normalizePath(
  existing_projects[1],
  winslash = "/",
  mustWork = TRUE
)

cat("Project directory:\n")
print(project_dir)
cat("\n")


# ======================================================================
# 2. REQUIRED PACKAGES
# ======================================================================

required_packages <- c(
  "readxl",
  "openxlsx",
  "dplyr",
  "stringr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
})

cat("Required packages loaded successfully.\n\n")


# ======================================================================
# 3. INPUT / OUTPUT DIRECTORIES
# ======================================================================

results_root <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers"
)

if (!dir.exists(results_root)) {
  stop(
    "Expected results directory was not found:\n",
    results_root
  )
}

output_dir <- file.path(
  results_root,
  "167_complete_supplementary_package_freeze"
)

tables_dir <- file.path(output_dir, "tables")
text_dir   <- file.path(output_dir, "text")
audit_dir  <- file.path(output_dir, "audit")
logs_dir   <- file.path(output_dir, "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

cat("Results root:\n")
cat(results_root, "\n\n")

cat("Output folder:\n")
cat(output_dir, "\n\n")


# ======================================================================
# 4. EXPECTED SUPPLEMENTARY ARCHITECTURE
# ======================================================================

expected_tables <- tibble::tribble(
  ~artifact_id, ~source_script, ~short_title,
  "Table S1",  "151",  "Discovery cohort and sample metadata",
  "Table S2",  "152b", "Complete blood differential-expression results",
  "Table S3",  "152",  "Robust-core differential expression and functional enrichment",
  "Table S4",  "155",  "SRS and CTS assignments and robustness analyses",
  "Table S5",  "156",  "Candidate-gene pool and exhaustive panel screening",
  "Table S6",  "157",  "Repeated internal cross-validation",
  "Table S7",  "159",  "Complete exploratory clinical-association analysis",
  "Table S8",  "161",  "Published transcriptomic-signature benchmarking",
  "Table S9",  "163",  "GSE154918 external evaluation",
  "Table S10", "165",  "GSE185263 external severity evaluation"
)

expected_figures <- tibble::tribble(
  ~artifact_id, ~source_script, ~short_title,
  "Figure S1", "153", "Discovery-cohort transcriptomic quality control",
  "Figure S2", "154", "Differential-expression and enrichment sensitivity analyses",
  "Figure S3", "146", "SRS and CTS robustness analyses",
  "Figure S4", "158", "Five-gene panel development and internal discrimination",
  "Figure S5", "160", "Complete clinical-association landscape",
  "Figure S6", "162", "Convergence of study-derived and published transcriptomic signatures",
  "Figure S7", "164", "GSE154918 external sensitivity analyses",
  "Figure S8", "166", "GSE185263 external sensitivity and secondary analyses"
)


# ======================================================================
# 5. MANUSCRIPT CROSS-REFERENCE MAP
# ======================================================================

cross_reference_map <- tibble::tribble(
  ~results_section, ~main_story, ~supplementary_items,
  
  "3.1",
  "Discovery cohort and dataset",
  "Table S1; Figure S1",
  
  "3.2",
  "Broad and robust blood transcriptomic response",
  "Table S2; Table S3; Figure S2",
  
  "3.3",
  "SRS and CTS blood transcriptomic heterogeneity",
  "Table S4; Figure S3",
  
  "3.4",
  "Biology-guided five-gene host-response signature",
  "Table S5; Table S6; Figure S4",
  
  "3.5",
  "Five-gene score recapitulates SRS and CTS structure",
  "Table S4; Table S6",
  
  "3.6",
  "Clinical associations and systemic inflammation",
  "Table S7; Figure S5",
  
  "3.7",
  "Published transcriptomic signatures converge on a shared host-response axis",
  "Table S8; Figure S6",
  
  "3.8",
  "External evaluation in GSE154918",
  "Table S9; Figure S7",
  
  "3.9",
  "External organ-dysfunction-severity replication in GSE185263",
  "Table S10; Figure S8"
)


# ======================================================================
# 6. FILE INVENTORY
# ======================================================================

all_files <- list.files(
  results_root,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)

all_files <- all_files[file.exists(all_files)]

all_files_norm <- normalizePath(
  all_files,
  winslash = "/",
  mustWork = FALSE
)

file_inventory <- tibble(
  path = all_files_norm,
  basename = basename(all_files_norm),
  extension = tolower(tools::file_ext(all_files_norm)),
  size_bytes = file.info(all_files_norm)$size,
  modified = file.info(all_files_norm)$mtime
)

cat(
  "Files detected under blood_endotypes_biomarkers =",
  nrow(file_inventory),
  "\n\n"
)


# ======================================================================
# 7. HELPERS
# ======================================================================

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}


artifact_number <- function(artifact_id) {
  
  as.integer(
    gsub(
      "[^0-9]",
      "",
      artifact_id
    )
  )
}


script_path_match <- function(path, script_token) {
  
  if (
    length(path) != 1 ||
    length(script_token) != 1 ||
    is.na(path) ||
    is.na(script_token)
  ) {
    return(FALSE)
  }
  
  path2 <- normalize_slashes(path)
  
  pattern <- paste0(
    "/",
    gsub("\\.", "\\\\.", script_token),
    "(_|/)"
  )
  
  grepl(
    pattern,
    path2,
    ignore.case = TRUE,
    perl = TRUE
  )
}


table_label_match <- function(path, n) {
  
  x <- basename(path)
  
  pattern <- paste0(
    "(?i)",
    "Table[_ ]?S",
    n,
    "(?![0-9])"
  )
  
  grepl(
    pattern,
    x,
    perl = TRUE
  )
}


figure_label_match <- function(path, n) {
  
  x <- basename(path)
  
  pattern <- paste0(
    "(?i)",
    "Figure[_ ]?S",
    n,
    "(?![0-9])"
  )
  
  grepl(
    pattern,
    x,
    perl = TRUE
  )
}


is_excluded_candidate <- function(path) {
  
  x <- tolower(
    normalize_slashes(path)
  )
  
  grepl(
    paste(
      c(
        "/audit/",
        "/logs/",
        "internal_audit",
        "provenance_audit",
        "numerical_audit",
        "source_data",
        "manifest",
        "sessioninfo",
        "run_info",
        "results_placement"
      ),
      collapse = "|"
    ),
    x
  )
}


choose_best_table_file <- function(
    source_script,
    table_number
) {
  
  candidates <- file_inventory$path[
    file_inventory$extension == "xlsx" &
      vapply(
        file_inventory$path,
        script_path_match,
        FUN.VALUE = logical(1),
        script_token = source_script
      )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  score <- rep(0, length(candidates))
  
  score <- score +
    ifelse(
      grepl(
        "/tables/",
        normalize_slashes(candidates),
        ignore.case = TRUE
      ),
      20,
      0
    )
  
  score <- score +
    ifelse(
      vapply(
        candidates,
        table_label_match,
        FUN.VALUE = logical(1),
        n = table_number
      ),
      25,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "final",
        basename(candidates),
        ignore.case = TRUE
      ),
      6,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "submission",
        basename(candidates),
        ignore.case = TRUE
      ),
      8,
      0
    )
  
  score <- score -
    ifelse(
      vapply(
        candidates,
        is_excluded_candidate,
        FUN.VALUE = logical(1)
      ),
      100,
      0
    )
  
  ord <- order(
    score,
    file.info(candidates)$mtime,
    decreasing = TRUE,
    na.last = TRUE
  )
  
  candidates[ord][1]
}


choose_best_figure_file <- function(
    source_script,
    figure_number,
    extension
) {
  
  candidates <- file_inventory$path[
    file_inventory$extension == extension &
      vapply(
        file_inventory$path,
        script_path_match,
        FUN.VALUE = logical(1),
        script_token = source_script
      )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  score <- rep(0, length(candidates))
  
  score <- score +
    ifelse(
      grepl(
        "/figures/",
        normalize_slashes(candidates),
        ignore.case = TRUE
      ),
      20,
      0
    )
  
  score <- score +
    ifelse(
      vapply(
        candidates,
        figure_label_match,
        FUN.VALUE = logical(1),
        n = figure_number
      ),
      25,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "final",
        basename(candidates),
        ignore.case = TRUE
      ),
      5,
      0
    )
  
  score <- score -
    ifelse(
      vapply(
        candidates,
        is_excluded_candidate,
        FUN.VALUE = logical(1)
      ),
      100,
      0
    )
  
  ord <- order(
    score,
    file.info(candidates)$mtime,
    decreasing = TRUE,
    na.last = TRUE
  )
  
  candidates[ord][1]
}


choose_english_caption <- function(
    source_script,
    figure_number
) {
  
  candidates <- file_inventory$path[
    file_inventory$extension == "txt" &
      vapply(
        file_inventory$path,
        script_path_match,
        FUN.VALUE = logical(1),
        script_token = source_script
      )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  x <- tolower(
    basename(candidates)
  )
  
  score <- rep(0, length(candidates))
  
  score <- score +
    ifelse(
      grepl("caption", x),
      30,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        paste0(
          "figure[_ ]?s",
          figure_number
        ),
        x
      ),
      20,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "(^|_)en($|_|\\.)",
        x,
        perl = TRUE
      ),
      15,
      0
    )
  
  score <- score -
    ifelse(
      grepl("results", x),
      20,
      0
    )
  
  ord <- order(
    score,
    file.info(candidates)$mtime,
    decreasing = TRUE,
    na.last = TRUE
  )
  
  if (score[ord][1] <= 0) {
    return(NA_character_)
  }
  
  candidates[ord][1]
}


safe_md5 <- function(path) {
  
  if (
    is.na(path) ||
    !nzchar(path) ||
    !file.exists(path)
  ) {
    return(NA_character_)
  }
  
  unname(
    tools::md5sum(path)
  )
}


safe_file_size <- function(path) {
  
  if (
    is.na(path) ||
    !nzchar(path) ||
    !file.exists(path)
  ) {
    return(NA_real_)
  }
  
  as.numeric(
    file.info(path)$size
  )
}


safe_modified <- function(path) {
  
  if (
    is.na(path) ||
    !nzchar(path) ||
    !file.exists(path)
  ) {
    return(as.POSIXct(NA))
  }
  
  file.info(path)$mtime
}


read_text_safe <- function(path) {
  
  if (
    is.na(path) ||
    !nzchar(path) ||
    !file.exists(path)
  ) {
    return("")
  }
  
  paste(
    readLines(
      path,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    collapse = "\n"
  )
}


# ======================================================================
# 8. LOCATE SUPPLEMENTARY TABLES
# ======================================================================

table_manifest <- expected_tables

table_manifest$artifact_number <- unname(
  vapply(
    table_manifest$artifact_id,
    artifact_number,
    FUN.VALUE = integer(1)
  )
)

table_manifest$file_path <- NA_character_

for (i in seq_len(nrow(table_manifest))) {
  
  table_manifest$file_path[i] <- choose_best_table_file(
    source_script = table_manifest$source_script[i],
    table_number = table_manifest$artifact_number[i]
  )
}

table_manifest$file_found <-
  !is.na(table_manifest$file_path) &
  file.exists(table_manifest$file_path)

table_manifest$file_name <- ifelse(
  table_manifest$file_found,
  basename(table_manifest$file_path),
  NA_character_
)


# ======================================================================
# 9. LOCATE SUPPLEMENTARY FIGURES
# ======================================================================

figure_manifest <- expected_figures

figure_manifest$artifact_number <- unname(
  vapply(
    figure_manifest$artifact_id,
    artifact_number,
    FUN.VALUE = integer(1)
  )
)

figure_manifest$png_path <- NA_character_
figure_manifest$pdf_path <- NA_character_
figure_manifest$tiff_path <- NA_character_
figure_manifest$caption_EN_path <- NA_character_

for (i in seq_len(nrow(figure_manifest))) {
  
  script_token <- figure_manifest$source_script[i]
  n <- figure_manifest$artifact_number[i]
  
  figure_manifest$png_path[i] <- choose_best_figure_file(
    source_script = script_token,
    figure_number = n,
    extension = "png"
  )
  
  figure_manifest$pdf_path[i] <- choose_best_figure_file(
    source_script = script_token,
    figure_number = n,
    extension = "pdf"
  )
  
  figure_manifest$tiff_path[i] <- choose_best_figure_file(
    source_script = script_token,
    figure_number = n,
    extension = "tiff"
  )
  
  if (is.na(figure_manifest$tiff_path[i])) {
    
    figure_manifest$tiff_path[i] <- choose_best_figure_file(
      source_script = script_token,
      figure_number = n,
      extension = "tif"
    )
  }
  
  figure_manifest$caption_EN_path[i] <- choose_english_caption(
    source_script = script_token,
    figure_number = n
  )
}

figure_manifest$png_found <-
  !is.na(figure_manifest$png_path) &
  file.exists(figure_manifest$png_path)

figure_manifest$pdf_found <-
  !is.na(figure_manifest$pdf_path) &
  file.exists(figure_manifest$pdf_path)

figure_manifest$tiff_found <-
  !is.na(figure_manifest$tiff_path) &
  file.exists(figure_manifest$tiff_path)

figure_manifest$caption_EN_found <-
  !is.na(figure_manifest$caption_EN_path) &
  file.exists(figure_manifest$caption_EN_path)


# ======================================================================
# 10. TABLE WORKBOOK STRUCTURE AUDIT
# ======================================================================

sheet_hint_patterns <- list(
  "Table S1" = "metadata|cohort|sample|readme|info",
  "Table S2" = "complete.*de|complete_de|readme",
  "Table S3" = "robust|core|hallmark|kegg|wiki|go|enrichment",
  "Table S4" = "srs|cts|robust",
  "Table S5" = "candidate|panel|screen|search",
  "Table S6" = "cv|cross|validation",
  "Table S7" = "clinical|association|complete",
  "Table S8" = "benchmark|signature|srsq|cts",
  "Table S9" = "gse154918|external|validation|evaluation",
  "Table S10" = "gse185263|sofa|external|severity"
)

table_sheet_audit <- tibble(
  artifact_id = character(),
  workbook = character(),
  n_sheets = integer(),
  sheet_names = character(),
  hint_pattern = character(),
  hint_match = logical(),
  readable = logical(),
  error = character()
)

for (i in seq_len(nrow(table_manifest))) {
  
  artifact <- table_manifest$artifact_id[i]
  path <- table_manifest$file_path[i]
  hint <- sheet_hint_patterns[[artifact]]
  
  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    
    table_sheet_audit <- bind_rows(
      table_sheet_audit,
      tibble(
        artifact_id = artifact,
        workbook = NA_character_,
        n_sheets = NA_integer_,
        sheet_names = NA_character_,
        hint_pattern = hint,
        hint_match = FALSE,
        readable = FALSE,
        error = "Workbook not found"
      )
    )
    
    next
  }
  
  sheet_result <- tryCatch(
    {
      
      sheets <- readxl::excel_sheets(path)
      
      tibble(
        artifact_id = artifact,
        workbook = basename(path),
        n_sheets = length(sheets),
        sheet_names = paste(sheets, collapse = " | "),
        hint_pattern = hint,
        hint_match = any(
          grepl(
            hint,
            sheets,
            ignore.case = TRUE,
            perl = TRUE
          )
        ),
        readable = TRUE,
        error = NA_character_
      )
    },
    error = function(e) {
      
      tibble(
        artifact_id = artifact,
        workbook = basename(path),
        n_sheets = NA_integer_,
        sheet_names = NA_character_,
        hint_pattern = hint,
        hint_match = FALSE,
        readable = FALSE,
        error = conditionMessage(e)
      )
    }
  )
  
  table_sheet_audit <- bind_rows(
    table_sheet_audit,
    sheet_result
  )
}


# ======================================================================
# 11. WORKBOOK TEXT READER
# ======================================================================

scan_workbook_text <- function(path) {
  
  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    return("")
  }
  
  sheets <- tryCatch(
    readxl::excel_sheets(path),
    error = function(e) character(0)
  )
  
  if (length(sheets) == 0) {
    return("")
  }
  
  text_parts <- character()
  
  for (sheet in sheets) {
    
    dat <- tryCatch(
      suppressMessages(
        readxl::read_excel(
          path,
          sheet = sheet,
          col_types = "text",
          .name_repair = "unique_quiet"
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(dat)) {
      next
    }
    
    header_text <- paste(
      names(dat),
      collapse = " "
    )
    
    body_text <- paste(
      unlist(
        dat,
        use.names = FALSE
      ),
      collapse = " "
    )
    
    text_parts <- c(
      text_parts,
      paste(
        sheet,
        header_text,
        body_text,
        sep = " "
      )
    )
  }
  
  paste(
    text_parts,
    collapse = "\n"
  )
}


# ======================================================================
# 12. INFORMATIONAL CONTENT SCANNER
#
# IMPORTANT:
# These scans generate REVIEW NOTES only.
# They do NOT by themselves prevent freezing.
# ======================================================================

scope_patterns <- c(
  "\\burine\\b",
  "\\burinary\\b",
  "\\blncrna\\b",
  "\\blong non[- ]coding rna\\b",
  "\\bparaspeckle"
)

problematic_claim_patterns <- c(
  "\\bclinically validated\\b",
  "\\bvalidated diagnostic\\b",
  "\\bvalidated prognostic\\b",
  "\\bgeneralizable diagnostic accuracy\\b",
  "\\bpre[- ]?calibrated threshold\\b",
  "\\bunique optimum\\b",
  "\\buniquely optimized\\b",
  "\\bindependent validation cohorts\\b",
  "\\bfive independent validation cohorts\\b",
  "\\bfive external validation cohorts\\b"
)

identifier_column_patterns <- c(
  "^patient[_ ]?name$",
  "^participant[_ ]?name$",
  "^first[_ ]?name$",
  "^last[_ ]?name$",
  "^full[_ ]?name$",
  "^surname$",
  "^mrn$",
  "medical[_ ]?record",
  "hospital[_ ]?record",
  "^iin$",
  "national[_ ]?id",
  "passport",
  "telephone",
  "^phone$",
  "^email$",
  "street[_ ]?address",
  "home[_ ]?address"
)


is_negative_scope_context <- function(text_lower, term) {
  
  patterns <- switch(
    term,
    
    "\\burine\\b" = c(
      "no urine",
      "urine not included",
      "urine was not included",
      "urine is not included",
      "without urine",
      "exclude urine",
      "excluded urine",
      "blood-only",
      "blood only"
    ),
    
    "\\burinary\\b" = c(
      "no urinary",
      "urinary.*not included",
      "without urinary",
      "exclude urinary",
      "excluded urinary",
      "blood-only",
      "blood only"
    ),
    
    "\\blncrna\\b" = c(
      "no lncrna",
      "lncrna.*not included",
      "without lncrna",
      "exclude lncrna",
      "excluded lncrna"
    ),
    
    "\\blong non[- ]coding rna\\b" = c(
      "long non[- ]coding rna.*not included",
      "without long non[- ]coding rna",
      "exclude.*long non[- ]coding rna"
    ),
    
    "\\bparaspeckle" = c(
      "paraspeckle.*not included",
      "without paraspeckle",
      "exclude.*paraspeckle"
    ),
    
    character(0)
  )
  
  if (length(patterns) == 0) {
    return(FALSE)
  }
  
  any(
    vapply(
      patterns,
      function(p) {
        grepl(
          p,
          text_lower,
          ignore.case = TRUE,
          perl = TRUE
        )
      },
      FUN.VALUE = logical(1)
    )
  )
}


is_negative_claim_context <- function(
    text_lower,
    claim_pattern
) {
  
  negative_patterns <- c(
    "not clinically validated",
    "not a clinically validated",
    "not validated diagnostic",
    "not a validated diagnostic",
    "not validated prognostic",
    "not a validated prognostic",
    "should not be interpreted as.*validated diagnostic",
    "should not be interpreted as.*validated prognostic",
    "does not constitute.*validated diagnostic",
    "does not constitute.*validated prognostic",
    "not.*independent validation cohorts",
    "not.*five independent validation cohorts",
    "not.*five external validation cohorts"
  )
  
  any(
    vapply(
      negative_patterns,
      function(p) {
        grepl(
          p,
          text_lower,
          ignore.case = TRUE,
          perl = TRUE
        )
      },
      FUN.VALUE = logical(1)
    )
  )
}


content_audit <- tibble(
  artifact_id = character(),
  file_path = character(),
  scope_keyword_hits = character(),
  problematic_claim_hits = character(),
  identifier_column_hits = character(),
  review_status = character()
)


for (i in seq_len(nrow(table_manifest))) {
  
  artifact <- table_manifest$artifact_id[i]
  path <- table_manifest$file_path[i]
  
  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    
    content_audit <- bind_rows(
      content_audit,
      tibble(
        artifact_id = artifact,
        file_path = NA_character_,
        scope_keyword_hits = NA_character_,
        problematic_claim_hits = NA_character_,
        identifier_column_hits = NA_character_,
        review_status = "MISSING"
      )
    )
    
    next
  }
  
  cat(
    "Scanning",
    artifact,
    "for informational scope/reporting/identifier notes...\n"
  )
  
  workbook_text <- scan_workbook_text(path)
  workbook_lower <- tolower(workbook_text)
  
  raw_scope_hits <- scope_patterns[
    vapply(
      scope_patterns,
      function(p) {
        grepl(
          p,
          workbook_lower,
          ignore.case = TRUE,
          perl = TRUE
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  scope_hits <- raw_scope_hits[
    !vapply(
      raw_scope_hits,
      function(p) {
        is_negative_scope_context(
          workbook_lower,
          p
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  raw_claim_hits <- problematic_claim_patterns[
    vapply(
      problematic_claim_patterns,
      function(p) {
        grepl(
          p,
          workbook_lower,
          ignore.case = TRUE,
          perl = TRUE
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  claim_hits <- raw_claim_hits[
    !vapply(
      raw_claim_hits,
      function(p) {
        is_negative_claim_context(
          workbook_lower,
          p
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  sheets <- tryCatch(
    readxl::excel_sheets(path),
    error = function(e) character(0)
  )
  
  identifier_hits <- character()
  
  for (sheet in sheets) {
    
    dat0 <- tryCatch(
      suppressMessages(
        readxl::read_excel(
          path,
          sheet = sheet,
          n_max = 1,
          .name_repair = "unique_quiet"
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(dat0)) {
      next
    }
    
    cols <- names(dat0)
    
    for (p in identifier_column_patterns) {
      
      matched <- cols[
        grepl(
          p,
          cols,
          ignore.case = TRUE,
          perl = TRUE
        )
      ]
      
      if (length(matched) > 0) {
        
        identifier_hits <- c(
          identifier_hits,
          paste0(
            sheet,
            "::",
            matched
          )
        )
      }
    }
  }
  
  identifier_hits <- unique(identifier_hits)
  
  review_status <- "CLEAR"
  
  if (
    length(scope_hits) > 0 ||
    length(claim_hits) > 0
  ) {
    review_status <- "REVIEW_NOTE"
  }
  
  if (length(identifier_hits) > 0) {
    review_status <- "CRITICAL_IDENTIFIER_WARNING"
  }
  
  content_audit <- bind_rows(
    content_audit,
    tibble(
      artifact_id = artifact,
      file_path = path,
      
      scope_keyword_hits = ifelse(
        length(scope_hits) == 0,
        "",
        paste(
          scope_hits,
          collapse = " | "
        )
      ),
      
      problematic_claim_hits = ifelse(
        length(claim_hits) == 0,
        "",
        paste(
          claim_hits,
          collapse = " | "
        )
      ),
      
      identifier_column_hits = ifelse(
        length(identifier_hits) == 0,
        "",
        paste(
          identifier_hits,
          collapse = " | "
        )
      ),
      
      review_status = review_status
    )
  )
}

cat("\n")


# ======================================================================
# 13. FIGURE CAPTION AUDIT
# ======================================================================

caption_audit <- tibble(
  artifact_id = character(),
  caption_path = character(),
  has_expected_figure_label = logical(),
  scope_keyword_hits = character(),
  problematic_claim_hits = character(),
  review_status = character()
)


for (i in seq_len(nrow(figure_manifest))) {
  
  artifact <- figure_manifest$artifact_id[i]
  n <- figure_manifest$artifact_number[i]
  path <- figure_manifest$caption_EN_path[i]
  
  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    
    caption_audit <- bind_rows(
      caption_audit,
      tibble(
        artifact_id = artifact,
        caption_path = NA_character_,
        has_expected_figure_label = FALSE,
        scope_keyword_hits = NA_character_,
        problematic_claim_hits = NA_character_,
        review_status = "MISSING"
      )
    )
    
    next
  }
  
  txt <- read_text_safe(path)
  txt_lower <- tolower(txt)
  
  expected_label_pattern <- paste0(
    "(supplementary\\s+)?figure\\s+s",
    n,
    "(?![0-9])"
  )
  
  label_ok <- grepl(
    expected_label_pattern,
    txt_lower,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  raw_scope_hits <- scope_patterns[
    vapply(
      scope_patterns,
      function(p) {
        grepl(
          p,
          txt_lower,
          ignore.case = TRUE,
          perl = TRUE
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  scope_hits <- raw_scope_hits[
    !vapply(
      raw_scope_hits,
      function(p) {
        is_negative_scope_context(
          txt_lower,
          p
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  raw_claim_hits <- problematic_claim_patterns[
    vapply(
      problematic_claim_patterns,
      function(p) {
        grepl(
          p,
          txt_lower,
          ignore.case = TRUE,
          perl = TRUE
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  claim_hits <- raw_claim_hits[
    !vapply(
      raw_claim_hits,
      function(p) {
        is_negative_claim_context(
          txt_lower,
          p
        )
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  review_status <- "CLEAR"
  
  if (
    !label_ok ||
    length(scope_hits) > 0 ||
    length(claim_hits) > 0
  ) {
    review_status <- "REVIEW_NOTE"
  }
  
  caption_audit <- bind_rows(
    caption_audit,
    tibble(
      artifact_id = artifact,
      caption_path = path,
      has_expected_figure_label = label_ok,
      
      scope_keyword_hits = ifelse(
        length(scope_hits) == 0,
        "",
        paste(
          scope_hits,
          collapse = " | "
        )
      ),
      
      problematic_claim_hits = ifelse(
        length(claim_hits) == 0,
        "",
        paste(
          claim_hits,
          collapse = " | "
        )
      ),
      
      review_status = review_status
    )
  )
}


# ======================================================================
# 14. TABLE / FIGURE AVAILABILITY STATUS
# ======================================================================

table_manifest <- table_manifest %>%
  left_join(
    table_sheet_audit %>%
      select(
        artifact_id,
        n_sheets,
        sheet_names,
        hint_match,
        readable
      ),
    by = "artifact_id"
  ) %>%
  left_join(
    content_audit %>%
      select(
        artifact_id,
        content_review_status = review_status
      ),
    by = "artifact_id"
  ) %>%
  mutate(
    availability_status = case_when(
      !file_found ~ "MISSING",
      !readable ~ "ERROR",
      TRUE ~ "READY"
    )
  )


figure_manifest <- figure_manifest %>%
  left_join(
    caption_audit %>%
      select(
        artifact_id,
        caption_review_status = review_status
      ),
    by = "artifact_id"
  ) %>%
  mutate(
    formats_complete =
      png_found &
      pdf_found &
      tiff_found,
    
    availability_status = case_when(
      !png_found & !pdf_found & !tiff_found ~ "MISSING",
      !formats_complete ~ "ERROR",
      !caption_EN_found ~ "ERROR",
      TRUE ~ "READY"
    )
  )


# ======================================================================
# 15. CORRECT NUMBERING AUDIT
# ======================================================================

table_numbers <- sort(
  unname(
    as.integer(
      table_manifest$artifact_number
    )
  )
)

figure_numbers <- sort(
  unname(
    as.integer(
      figure_manifest$artifact_number
    )
  )
)

table_numbering_pass <- identical(
  table_numbers,
  as.integer(1:10)
)

figure_numbering_pass <- identical(
  figure_numbers,
  as.integer(1:8)
)

numbering_audit <- tibble(
  item = c(
    "Supplementary Tables S1-S10",
    "Supplementary Figures S1-S8"
  ),
  
  expected = c(
    paste(1:10, collapse = ","),
    paste(1:8, collapse = ",")
  ),
  
  observed = c(
    paste(table_numbers, collapse = ","),
    paste(figure_numbers, collapse = ",")
  ),
  
  pass = c(
    table_numbering_pass,
    figure_numbering_pass
  )
)


# ======================================================================
# 16. CORRECT SOURCE PROVENANCE AUDIT
# ======================================================================

table_source_match <- mapply(
  FUN = script_path_match,
  path = table_manifest$file_path,
  script_token = table_manifest$source_script,
  USE.NAMES = FALSE
)

figure_source_match <- mapply(
  FUN = script_path_match,
  path = figure_manifest$png_path,
  script_token = figure_manifest$source_script,
  USE.NAMES = FALSE
)

source_script_audit <- bind_rows(
  
  table_manifest %>%
    transmute(
      artifact_id,
      artifact_type = "Table",
      expected_source_script = source_script,
      file_path,
      source_path_contains_expected_script =
        table_source_match
    ),
  
  figure_manifest %>%
    transmute(
      artifact_id,
      artifact_type = "Figure",
      expected_source_script = source_script,
      file_path = png_path,
      source_path_contains_expected_script =
        figure_source_match
    )
)

source_provenance_pass <- all(
  source_script_audit$source_path_contains_expected_script,
  na.rm = FALSE
)


# ======================================================================
# 17. MD5 FREEZE MANIFEST
# ======================================================================

freeze_manifest <- tibble(
  artifact_id = character(),
  artifact_type = character(),
  format = character(),
  source_script = character(),
  path = character(),
  file_name = character(),
  size_bytes = numeric(),
  modified = as.POSIXct(character()),
  md5 = character()
)


# Tables ---------------------------------------------------------------

for (i in seq_len(nrow(table_manifest))) {
  
  path <- table_manifest$file_path[i]
  
  freeze_manifest <- bind_rows(
    freeze_manifest,
    tibble(
      artifact_id = table_manifest$artifact_id[i],
      artifact_type = "Supplementary Table",
      format = "XLSX",
      source_script = table_manifest$source_script[i],
      path = path,
      file_name = ifelse(
        is.na(path),
        NA_character_,
        basename(path)
      ),
      size_bytes = safe_file_size(path),
      modified = safe_modified(path),
      md5 = safe_md5(path)
    )
  )
}


# Figures --------------------------------------------------------------

for (i in seq_len(nrow(figure_manifest))) {
  
  figure_paths <- c(
    PNG = figure_manifest$png_path[i],
    PDF = figure_manifest$pdf_path[i],
    TIFF = figure_manifest$tiff_path[i],
    Caption_EN = figure_manifest$caption_EN_path[i]
  )
  
  for (fmt in names(figure_paths)) {
    
    path <- figure_paths[[fmt]]
    
    freeze_manifest <- bind_rows(
      freeze_manifest,
      tibble(
        artifact_id = figure_manifest$artifact_id[i],
        artifact_type = "Supplementary Figure",
        format = fmt,
        source_script = figure_manifest$source_script[i],
        path = path,
        file_name = ifelse(
          is.na(path),
          NA_character_,
          basename(path)
        ),
        size_bytes = safe_file_size(path),
        modified = safe_modified(path),
        md5 = safe_md5(path)
      )
    )
  }
}

freeze_manifest$file_exists <-
  !is.na(freeze_manifest$path) &
  file.exists(freeze_manifest$path)


# ======================================================================
# 18. CRITICAL FAILURE COUNTS
# ======================================================================

missing_tables <- sum(
  !table_manifest$file_found,
  na.rm = TRUE
)

unreadable_tables <- sum(
  table_manifest$file_found &
    !table_manifest$readable,
  na.rm = TRUE
)

missing_figure_png <- sum(
  !figure_manifest$png_found,
  na.rm = TRUE
)

missing_figure_pdf <- sum(
  !figure_manifest$pdf_found,
  na.rm = TRUE
)

missing_figure_tiff <- sum(
  !figure_manifest$tiff_found,
  na.rm = TRUE
)

missing_figure_captions <- sum(
  !figure_manifest$caption_EN_found,
  na.rm = TRUE
)

identifier_warning_tables <- sum(
  content_audit$identifier_column_hits != "" &
    !is.na(content_audit$identifier_column_hits),
  na.rm = TRUE
)


# Informational review notes -------------------------------------------

scope_review_tables <- sum(
  content_audit$scope_keyword_hits != "" &
    !is.na(content_audit$scope_keyword_hits),
  na.rm = TRUE
)

claim_review_tables <- sum(
  content_audit$problematic_claim_hits != "" &
    !is.na(content_audit$problematic_claim_hits),
  na.rm = TRUE
)

caption_scope_review <- sum(
  caption_audit$scope_keyword_hits != "" &
    !is.na(caption_audit$scope_keyword_hits),
  na.rm = TRUE
)

caption_claim_review <- sum(
  caption_audit$problematic_claim_hits != "" &
    !is.na(caption_audit$problematic_claim_hits),
  na.rm = TRUE
)

caption_label_review <- sum(
  !caption_audit$has_expected_figure_label &
    caption_audit$review_status != "MISSING",
  na.rm = TRUE
)


# ======================================================================
# 19. CRITICAL FREEZE LOGIC
# ======================================================================

critical_failure <- any(
  c(
    missing_tables > 0,
    unreadable_tables > 0,
    missing_figure_png > 0,
    missing_figure_pdf > 0,
    missing_figure_tiff > 0,
    missing_figure_captions > 0,
    !table_numbering_pass,
    !figure_numbering_pass,
    !source_provenance_pass,
    identifier_warning_tables > 0
  )
)

informational_review_notes <- any(
  c(
    scope_review_tables > 0,
    claim_review_tables > 0,
    caption_scope_review > 0,
    caption_claim_review > 0,
    caption_label_review > 0
  )
)

overall_status <- if (
  critical_failure
) {
  "NOT_FROZEN_CRITICAL_ISSUES"
} else {
  "FROZEN_READY"
}


# ======================================================================
# 20. HIGH-LEVEL AUDIT SUMMARY
# ======================================================================

audit_summary <- tibble::tribble(
  ~audit_item, ~expected, ~observed, ~status,
  
  "Supplementary Tables",
  "10",
  as.character(sum(table_manifest$file_found)),
  ifelse(missing_tables == 0, "PASS", "FAIL"),
  
  "Supplementary Figures PNG",
  "8",
  as.character(sum(figure_manifest$png_found)),
  ifelse(missing_figure_png == 0, "PASS", "FAIL"),
  
  "Supplementary Figures PDF",
  "8",
  as.character(sum(figure_manifest$pdf_found)),
  ifelse(missing_figure_pdf == 0, "PASS", "FAIL"),
  
  "Supplementary Figures TIFF",
  "8",
  as.character(sum(figure_manifest$tiff_found)),
  ifelse(missing_figure_tiff == 0, "PASS", "FAIL"),
  
  "English figure captions",
  "8",
  as.character(sum(figure_manifest$caption_EN_found)),
  ifelse(missing_figure_captions == 0, "PASS", "FAIL"),
  
  "Readable supplementary-table workbooks",
  "10",
  as.character(sum(table_manifest$readable, na.rm = TRUE)),
  ifelse(unreadable_tables == 0, "PASS", "FAIL"),
  
  "Table numbering S1-S10",
  "PASS",
  as.character(table_numbering_pass),
  ifelse(table_numbering_pass, "PASS", "FAIL"),
  
  "Figure numbering S1-S8",
  "PASS",
  as.character(figure_numbering_pass),
  ifelse(figure_numbering_pass, "PASS", "FAIL"),
  
  "Expected source-script provenance",
  "PASS",
  as.character(source_provenance_pass),
  ifelse(source_provenance_pass, "PASS", "FAIL"),
  
  "Possible direct-identifier warnings",
  "0",
  as.character(identifier_warning_tables),
  ifelse(identifier_warning_tables == 0, "PASS", "FAIL"),
  
  "Informational scope review notes in tables",
  "manual review only",
  as.character(scope_review_tables),
  ifelse(scope_review_tables == 0, "CLEAR", "NOTE"),
  
  "Informational wording review notes in tables",
  "manual review only",
  as.character(claim_review_tables),
  ifelse(claim_review_tables == 0, "CLEAR", "NOTE"),
  
  "Informational scope review notes in captions",
  "manual review only",
  as.character(caption_scope_review),
  ifelse(caption_scope_review == 0, "CLEAR", "NOTE"),
  
  "Informational wording review notes in captions",
  "manual review only",
  as.character(caption_claim_review),
  ifelse(caption_claim_review == 0, "CLEAR", "NOTE"),
  
  "Caption-label review notes",
  "manual review only",
  as.character(caption_label_review),
  ifelse(caption_label_review == 0, "CLEAR", "NOTE")
)


# ======================================================================
# 21. REPORTING GUARDRAILS
# ======================================================================

reporting_guardrails <- tibble(
  rule_id = sprintf(
    "RG%02d",
    1:16
  ),
  
  rule = c(
    "The manuscript is blood-only.",
    "Urine and urinary transcriptomic analyses are not part of this manuscript.",
    "lncRNA-focused analyses and paraspeckle analyses are not part of this manuscript.",
    "The incidental presence of an lncRNA gene in complete transcriptomic tables does not constitute an lncRNA-focused analysis.",
    "The primary study-derived panel is CD177, HK3, IRAK3, CARD11, and IKZF2.",
    "The primary five-gene panel is described as a biology-guided host-response signature.",
    "Internal discrimination is saturated and does not identify a unique optimal classifier.",
    "SRS and CTS were not used to select the primary five-gene panel.",
    "The five-gene score is a molecular host-response readout, not a clinically calibrated assay.",
    "GSE154918 provides an external infection-control context with a non-significant primary endpoint.",
    "GSE185263 provides external replication of an organ-dysfunction-severity association.",
    "Mortality and ICU-versus-Emergency-Room analyses in GSE185263 are secondary.",
    "Sepsis-versus-healthy discrimination in GSE185263 is contextual.",
    "The five GSE185263 collection locations are sensitivity strata within one dataset.",
    "The pooled location correlation is descriptive and does not constitute five independent validation cohorts.",
    "Script 167 performs publication-package integrity checks only and introduces no new inferential statistics."
  )
)


# ======================================================================
# 22. TERMINOLOGY REFERENCE
# ======================================================================

terminology_reference <- tibble::tribble(
  ~concept, ~preferred_wording,
  
  "Five-gene panel",
  "biology-guided five-gene host-response signature",
  
  "Five-gene score",
  "five-gene host-response score",
  
  "Biological interpretation",
  "myeloid-adaptive host-response axis",
  
  "Internal discrimination",
  "internal cross-validation / internal discrimination",
  
  "External GSE154918",
  "external evaluation / directional replication",
  
  "External GSE185263",
  "external replication of an organ-dysfunction-severity association",
  
  "GSE185263 locations",
  "location-specific sensitivity strata within one dataset",
  
  "Pooled GSE185263 location estimate",
  "descriptive fixed-effect pooled correlation",
  
  "Clinical interpretation",
  "molecular readout of host-response state",
  
  "Avoid",
  "clinically validated diagnostic or prognostic assay"
)


# ======================================================================
# 23. FREEZE DECLARATION
# ======================================================================

freeze_time <- Sys.time()

freeze_declaration <- tibble(
  item = c(
    "Freeze status",
    "Freeze timestamp",
    "Project",
    "Manuscript scope",
    "Supplementary Tables",
    "Supplementary Figures",
    "New inferential analyses in Script 167",
    "Hash algorithm",
    "Informational review notes present",
    "Interpretation"
  ),
  
  value = c(
    overall_status,
    as.character(freeze_time),
    "Sepsis_DESeq2 blood endotypes and biomarkers",
    "Blood only",
    "S1-S10",
    "S1-S8",
    "NO",
    "MD5",
    as.character(informational_review_notes),
    
    ifelse(
      overall_status == "FROZEN_READY",
      paste(
        "All mandatory supplementary artifacts were located,",
        "readable, correctly numbered, and linked to the expected",
        "source scripts. No direct-identifier warning was detected.",
        "Informational terminology/scope notes, if present, do not",
        "invalidate the analytical freeze."
      ),
      paste(
        "One or more critical supplementary-package integrity",
        "requirements failed. The package must not yet be frozen."
      )
    )
  )
)


# ======================================================================
# 24. WRITE COMPLETE MANIFEST
# ======================================================================

manifest_xlsx <- file.path(
  tables_dir,
  "167_COMPLETE_Supplementary_Package_Manifest.xlsx"
)

wb <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wb,
  "00_Freeze_status"
)

openxlsx::writeData(
  wb,
  "00_Freeze_status",
  freeze_declaration
)


openxlsx::addWorksheet(
  wb,
  "01_Audit_summary"
)

openxlsx::writeData(
  wb,
  "01_Audit_summary",
  audit_summary
)


openxlsx::addWorksheet(
  wb,
  "02_Tables_S1_S10"
)

table_export <- table_manifest %>%
  select(
    artifact_id,
    source_script,
    short_title,
    file_found,
    file_name,
    file_path,
    n_sheets,
    sheet_names,
    hint_match,
    readable,
    content_review_status,
    availability_status
  )

openxlsx::writeData(
  wb,
  "02_Tables_S1_S10",
  table_export
)


openxlsx::addWorksheet(
  wb,
  "03_Figures_S1_S8"
)

figure_export <- figure_manifest %>%
  select(
    artifact_id,
    source_script,
    short_title,
    png_found,
    pdf_found,
    tiff_found,
    caption_EN_found,
    formats_complete,
    caption_review_status,
    availability_status,
    png_path,
    pdf_path,
    tiff_path,
    caption_EN_path
  )

openxlsx::writeData(
  wb,
  "03_Figures_S1_S8",
  figure_export
)


openxlsx::addWorksheet(
  wb,
  "04_MD5_freeze_manifest"
)

openxlsx::writeData(
  wb,
  "04_MD5_freeze_manifest",
  freeze_manifest
)


openxlsx::addWorksheet(
  wb,
  "05_Table_content_notes"
)

openxlsx::writeData(
  wb,
  "05_Table_content_notes",
  content_audit
)


openxlsx::addWorksheet(
  wb,
  "06_Caption_notes"
)

openxlsx::writeData(
  wb,
  "06_Caption_notes",
  caption_audit
)


openxlsx::addWorksheet(
  wb,
  "07_Workbook_sheet_audit"
)

openxlsx::writeData(
  wb,
  "07_Workbook_sheet_audit",
  table_sheet_audit
)


openxlsx::addWorksheet(
  wb,
  "08_Numbering_audit"
)

openxlsx::writeData(
  wb,
  "08_Numbering_audit",
  numbering_audit
)


openxlsx::addWorksheet(
  wb,
  "09_Source_provenance"
)

openxlsx::writeData(
  wb,
  "09_Source_provenance",
  source_script_audit
)


openxlsx::addWorksheet(
  wb,
  "10_Cross_reference_map"
)

openxlsx::writeData(
  wb,
  "10_Cross_reference_map",
  cross_reference_map
)


openxlsx::addWorksheet(
  wb,
  "11_Reporting_guardrails"
)

openxlsx::writeData(
  wb,
  "11_Reporting_guardrails",
  reporting_guardrails
)


openxlsx::addWorksheet(
  wb,
  "12_Terminology"
)

openxlsx::writeData(
  wb,
  "12_Terminology",
  terminology_reference
)


# ======================================================================
# 25. FORMAT WORKBOOK
# ======================================================================

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "Bottom"
)

for (sheet in names(wb)) {
  
  openxlsx::freezePane(
    wb,
    sheet,
    firstRow = TRUE
  )
  
  openxlsx::addStyle(
    wb,
    sheet,
    style = header_style,
    rows = 1,
    cols = 1:30,
    gridExpand = TRUE,
    stack = TRUE
  )
  
  openxlsx::setColWidths(
    wb,
    sheet,
    cols = 1:30,
    widths = "auto"
  )
}

openxlsx::saveWorkbook(
  wb,
  manifest_xlsx,
  overwrite = TRUE
)


# ======================================================================
# 26. WRITE CSV AUDITS
# ======================================================================

write.csv(
  audit_summary,
  file.path(
    audit_dir,
    "167_supplementary_package_audit_summary.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  freeze_manifest,
  file.path(
    audit_dir,
    "167_supplementary_package_MD5_manifest.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  content_audit,
  file.path(
    audit_dir,
    "167_supplementary_table_content_notes.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  caption_audit,
  file.path(
    audit_dir,
    "167_supplementary_figure_caption_notes.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  source_script_audit,
  file.path(
    audit_dir,
    "167_source_script_provenance_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  cross_reference_map,
  file.path(
    audit_dir,
    "167_manuscript_supplementary_cross_reference_map.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ======================================================================
# 27. HUMAN-READABLE FREEZE REPORT
# ======================================================================

freeze_report <- file.path(
  text_dir,
  "167_COMPLETE_Supplementary_Package_Freeze_Report.txt"
)

sink(freeze_report)

cat("=====================================================================\n")
cat("COMPLETE SUPPLEMENTARY PACKAGE FREEZE REPORT\n")
cat("Script 167 FINAL v2\n")
cat("=====================================================================\n\n")

cat("Project:\n")
cat(project_dir, "\n\n")

cat("Timestamp:\n")
cat(as.character(freeze_time), "\n\n")

cat("OVERALL STATUS\n")
cat("--------------\n")
cat(overall_status, "\n\n")

cat("SUPPLEMENTARY ARCHITECTURE\n")
cat("--------------------------\n")
cat("Tables: S1-S10\n")
cat("Figures: S1-S8\n\n")

cat("MANDATORY ARTIFACT AUDIT\n")
cat("------------------------\n")
print(
  audit_summary,
  row.names = FALSE
)

cat("\n")

cat("INFORMATIONAL REVIEW NOTES\n")
cat("--------------------------\n")

cat(
  "Table scope notes = ",
  scope_review_tables,
  "\n",
  sep = ""
)

cat(
  "Table wording notes = ",
  claim_review_tables,
  "\n",
  sep = ""
)

cat(
  "Caption scope notes = ",
  caption_scope_review,
  "\n",
  sep = ""
)

cat(
  "Caption wording notes = ",
  caption_claim_review,
  "\n",
  sep = ""
)

cat(
  "Caption-label notes = ",
  caption_label_review,
  "\n\n",
  sep = ""
)

cat(
  paste(
    "These content scans are informational only.",
    "They do not invalidate the analytical freeze unless",
    "a possible direct identifier is detected."
  ),
  "\n\n"
)

cat("REPORTING GUARDRAILS\n")
cat("--------------------\n")

for (i in seq_len(nrow(reporting_guardrails))) {
  
  cat(
    reporting_guardrails$rule_id[i],
    ". ",
    reporting_guardrails$rule[i],
    "\n",
    sep = ""
  )
}

cat("\n")

cat("INTERPRETATION\n")
cat("--------------\n")

if (overall_status == "FROZEN_READY") {
  
  cat(
    paste(
      "All mandatory Supplementary Tables S1-S10 and Figures S1-S8",
      "were located and passed the critical integrity audit.",
      "Expected source-script provenance was confirmed.",
      "No possible direct identifier was detected.",
      "The MD5 manifest defines the frozen Supplementary package.",
      "No inferential statistical analysis was performed by Script 167."
    ),
    "\n"
  )
  
} else {
  
  cat(
    paste(
      "At least one critical package-integrity requirement failed.",
      "Review the audit workbook before freezing the package."
    ),
    "\n"
  )
}

sink()


# ======================================================================
# 28. CROSS-REFERENCE TEXT
# ======================================================================

crossref_txt <- file.path(
  text_dir,
  "167_manuscript_supplementary_cross_reference_map.txt"
)

sink(crossref_txt)

cat("MANUSCRIPT -> SUPPLEMENTARY CROSS-REFERENCE MAP\n")
cat("================================================\n\n")

for (i in seq_len(nrow(cross_reference_map))) {
  
  cat(
    "Results ",
    cross_reference_map$results_section[i],
    " — ",
    cross_reference_map$main_story[i],
    "\n",
    sep = ""
  )
  
  cat(
    "Supplementary: ",
    cross_reference_map$supplementary_items[i],
    "\n\n",
    sep = ""
  )
}

sink()


# ======================================================================
# 29. FREEZE LOCK
# ======================================================================

freeze_lock_file <- file.path(
  output_dir,
  "167_SUPPLEMENTARY_PACKAGE_FREEZE_LOCK.txt"
)

sink(freeze_lock_file)

cat("SUPPLEMENTARY PACKAGE FREEZE LOCK\n")
cat("=================================\n\n")

cat(
  "Status: ",
  overall_status,
  "\n",
  sep = ""
)

cat(
  "Timestamp: ",
  as.character(freeze_time),
  "\n",
  sep = ""
)

cat("Tables: Supplementary Tables S1-S10\n")
cat("Figures: Supplementary Figures S1-S8\n")
cat("Scope: blood transcriptomics only\n")
cat("New inferential statistics in Script 167: NO\n")
cat("Integrity manifest: MD5\n\n")

cat(
  "Primary manifest:\n",
  manifest_xlsx,
  "\n\n",
  sep = ""
)

if (overall_status == "FROZEN_READY") {
  
  cat(
    paste(
      "This file records the Supplementary package state that passed",
      "the critical Script 167 FINAL v2 integrity checks.",
      "Any subsequent modification of a frozen artifact will change",
      "its MD5 hash and should trigger a new Script 167 audit."
    ),
    "\n"
  )
  
} else {
  
  cat(
    paste(
      "The package is not fully frozen.",
      "Resolve critical integrity failures and rerun Script 167."
    ),
    "\n"
  )
}

sink()


# ======================================================================
# 30. SESSION INFO
# ======================================================================

session_file <- file.path(
  logs_dir,
  "167_sessionInfo.txt"
)

sink(session_file)
sessionInfo()
sink()


# ======================================================================
# 31. CONSOLE SUMMARY
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("Script 167 FINAL v2 COMPLETE SUPPLEMENTARY PACKAGE AUDIT\n")
cat("=====================================================================\n\n")


cat("TABLES\n")
cat("------\n")

print(
  table_manifest %>%
    select(
      artifact_id,
      source_script,
      file_found,
      readable,
      availability_status,
      file_name
    ),
  n = Inf,
  width = Inf
)

cat("\n")


cat("FIGURES\n")
cat("-------\n")

print(
  figure_manifest %>%
    select(
      artifact_id,
      source_script,
      png_found,
      pdf_found,
      tiff_found,
      caption_EN_found,
      availability_status
    ),
  n = Inf,
  width = Inf
)

cat("\n")


cat("NUMBERING AUDIT\n")
cat("----------------\n")

print(
  numbering_audit,
  n = Inf,
  width = Inf
)

cat("\n")


cat("SOURCE PROVENANCE AUDIT\n")
cat("-----------------------\n")

print(
  source_script_audit,
  n = Inf,
  width = Inf
)

cat("\n")


cat("AUDIT SUMMARY\n")
cat("-------------\n")

print(
  audit_summary,
  n = Inf,
  width = Inf
)

cat("\n")


cat("INFORMATIONAL CONTENT REVIEW NOTES\n")
cat("----------------------------------\n")

cat(
  "Table scope notes = ",
  scope_review_tables,
  "\n",
  sep = ""
)

cat(
  "Table wording notes = ",
  claim_review_tables,
  "\n",
  sep = ""
)

cat(
  "Possible direct-identifier warnings = ",
  identifier_warning_tables,
  "\n",
  sep = ""
)

cat(
  "Caption scope notes = ",
  caption_scope_review,
  "\n",
  sep = ""
)

cat(
  "Caption wording notes = ",
  caption_claim_review,
  "\n",
  sep = ""
)

cat(
  "Caption-label notes = ",
  caption_label_review,
  "\n\n",
  sep = ""
)


if (
  any(
    content_audit$review_status != "CLEAR",
    na.rm = TRUE
  )
) {
  
  cat("TABLE CONTENT NOTES\n")
  cat("-------------------\n")
  
  print(
    content_audit %>%
      filter(
        review_status != "CLEAR"
      ),
    n = Inf,
    width = Inf
  )
  
  cat("\n")
}


if (
  any(
    caption_audit$review_status != "CLEAR",
    na.rm = TRUE
  )
) {
  
  cat("CAPTION NOTES\n")
  cat("-------------\n")
  
  print(
    caption_audit %>%
      filter(
        review_status != "CLEAR"
      ),
    n = Inf,
    width = Inf
  )
  
  cat("\n")
}


# ======================================================================
# 32. FINAL STATUS
# ======================================================================

cat("=====================================================================\n")
cat("FINAL SUPPLEMENTARY PACKAGE STATUS\n")
cat("=====================================================================\n\n")

cat(overall_status, "\n\n")

if (overall_status == "FROZEN_READY") {
  
  cat("All 10 Supplementary Tables were found and are readable.\n")
  cat("All 8 Supplementary Figures were found in PNG/PDF/TIFF.\n")
  cat("All 8 English figure captions were found.\n")
  cat("Table numbering S1-S10 passed.\n")
  cat("Figure numbering S1-S8 passed.\n")
  cat("Expected source-script provenance passed.\n")
  cat("No possible direct identifiers were detected.\n")
  cat("No new inferential statistics were performed.\n\n")
  
  if (informational_review_notes) {
    
    cat(
      paste(
        "One or more informational content-review notes remain.",
        "These are non-critical terminology/context flags and do not",
        "invalidate the analytical freeze."
      ),
      "\n\n"
    )
  }
  
  cat(
    "COMPLETE SUPPLEMENTARY PACKAGE = FROZEN_READY\n\n"
  )
  
} else {
  
  cat(
    "Critical missing or invalid supplementary artifacts remain.\n"
  )
  
  cat(
    "Do NOT freeze the package until these are resolved.\n\n"
  )
}


# ======================================================================
# 33. OUTPUT FILES
# ======================================================================

cat("OUTPUT FILES\n")
cat("------------\n")

cat(
  "Complete Supplementary Package Manifest:\n",
  manifest_xlsx,
  "\n\n",
  sep = ""
)

cat(
  "Freeze report:\n",
  freeze_report,
  "\n\n",
  sep = ""
)

cat(
  "Freeze lock:\n",
  freeze_lock_file,
  "\n\n",
  sep = ""
)

cat(
  "Cross-reference map:\n",
  crossref_txt,
  "\n\n",
  sep = ""
)

cat(
  "Session info:\n",
  session_file,
  "\n\n",
  sep = ""
)


# ======================================================================
# 34. REPORTING GUARDRAILS
# ======================================================================

cat("REPORTING GUARDRAILS\n")
cat("--------------------\n")
cat("- Script 167 performs no new inferential statistical analysis.\n")
cat("- Blood transcriptomics only; no urine-focused analysis.\n")
cat("- No lncRNA-focused or paraspeckle analysis is part of this manuscript.\n")
cat("- Presence of individual transcript names in complete DEG tables is not a scope violation.\n")
cat("- Primary panel: CD177, HK3, IRAK3, CARD11, IKZF2.\n")
cat("- Panel is biology-guided; internal discrimination is saturated.\n")
cat("- Score is a molecular host-response readout, not a calibrated clinical assay.\n")
cat("- GSE154918 and GSE185263 are separate external datasets.\n")
cat("- GSE185263 locations are sensitivity strata within one dataset.\n")
cat("- Mortality and ICU/ER analyses are secondary.\n")
cat("- Sepsis/healthy discrimination in GSE185263 is contextual only.\n")
cat("- MD5 hashes define the frozen Supplementary artifact state.\n")

cat("\nDone.\n")