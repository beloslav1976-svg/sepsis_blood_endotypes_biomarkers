# ======================================================================
# Script 168
# 168_audit_and_freeze_complete_main_manuscript_package.R
#
# FINAL
#
# COMPLETE MAIN MANUSCRIPT PACKAGE AUDIT AND FREEZE
#
# Blood transcriptomic endotypes / five-gene host-response manuscript
#
# PURPOSE
# -------
# Freeze the publication-facing MAIN manuscript assets:
#
#   Main Table 1
#   Main Figure 1
#   Main Figure 2
#   Main Figure 3
#   Main Figure 4
#   Main Figure 5
#
# and verify:
#
#   - required files are present
#   - figures exist in PNG / PDF / TIFF
#   - English captions exist
#   - expected source-script provenance is preserved
#   - Main Table 1 remains de-identified
#   - frozen Supplementary Package from Script 167 is present
#   - MD5 hashes define the frozen state
#
# NO NEW INFERENTIAL STATISTICAL ANALYSIS.
# NO MODIFICATION OF ANY FROZEN MANUSCRIPT ARTIFACT.
#
# Expected frozen sources:
#
#   Main Table 1 : Script 151
#   Figure 1     : Script 147
#   Figure 2     : Script 145b preferred; Script 145 allowed fallback
#   Figure 3     : Script 148
#   Figure 4     : Script 149
#   Figure 5     : Script 150
#
# ======================================================================


# ======================================================================
# 0. START
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("Running Script 168 FINAL\n")
cat("Complete Main Manuscript Package Audit and Freeze\n")
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

cat("Required packages loaded successfully.\n\n")


# ======================================================================
# 3. DIRECTORIES
# ======================================================================

results_root <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers"
)

if (!dir.exists(results_root)) {
  stop(
    "Expected results root not found:\n",
    results_root
  )
}

output_dir <- file.path(
  results_root,
  "168_complete_main_manuscript_package_freeze"
)

tables_dir <- file.path(output_dir, "tables")
audit_dir  <- file.path(output_dir, "audit")
text_dir   <- file.path(output_dir, "text")
logs_dir   <- file.path(output_dir, "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

cat("Results root:\n")
cat(results_root, "\n\n")

cat("Output folder:\n")
cat(output_dir, "\n\n")


# ======================================================================
# 4. FROZEN SUPPLEMENTARY PACKAGE REQUIREMENT
# ======================================================================

supplementary_freeze_lock <- file.path(
  results_root,
  "167_complete_supplementary_package_freeze",
  "167_SUPPLEMENTARY_PACKAGE_FREEZE_LOCK.txt"
)

supplementary_manifest <- file.path(
  results_root,
  "167_complete_supplementary_package_freeze",
  "tables",
  "167_COMPLETE_Supplementary_Package_Manifest.xlsx"
)

supplementary_lock_found <- file.exists(
  supplementary_freeze_lock
)

supplementary_manifest_found <- file.exists(
  supplementary_manifest
)

supplementary_lock_text <- ""

if (supplementary_lock_found) {
  
  supplementary_lock_text <- paste(
    readLines(
      supplementary_freeze_lock,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    collapse = "\n"
  )
}

supplementary_frozen_ready <- supplementary_lock_found &&
  grepl(
    "FROZEN_READY",
    supplementary_lock_text,
    fixed = TRUE
  )

cat("SUPPLEMENTARY PACKAGE DEPENDENCY\n")
cat("--------------------------------\n")

cat(
  "Freeze lock found =",
  supplementary_lock_found,
  "\n"
)

cat(
  "Manifest found =",
  supplementary_manifest_found,
  "\n"
)

cat(
  "Supplementary status FROZEN_READY =",
  supplementary_frozen_ready,
  "\n\n"
)


# ======================================================================
# 5. EXPECTED MAIN MANUSCRIPT ARCHITECTURE
# ======================================================================

expected_table <- tibble::tibble(
  artifact_id = "Main Table 1",
  source_scripts = "151",
  short_title = "Characteristics of the discovery blood transcriptomic cohort"
)

expected_figures <- tibble::tribble(
  ~artifact_id, ~figure_number, ~source_scripts, ~short_title,
  
  "Main Figure 1",
  1L,
  "147",
  "Robust blood transcriptomic response to sepsis",
  
  "Main Figure 2",
  2L,
  "145b;145",
  "SRS/CTS heterogeneity and five-gene host-response continuum",
  
  "Main Figure 3",
  3L,
  "148",
  "Clinical inflammatory context and published-signature convergence",
  
  "Main Figure 4",
  4L,
  "149",
  "External evaluation in GSE154918",
  
  "Main Figure 5",
  5L,
  "150",
  "External organ-dysfunction-severity replication in GSE185263"
)


# ======================================================================
# 6. MANUSCRIPT STORY MAP
# ======================================================================

story_map <- tibble::tribble(
  ~results_section,
  ~main_artifact,
  ~main_message,
  
  "3.1",
  "Main Table 1",
  "Discovery cohort and available clinical characteristics",
  
  "3.2",
  "Main Figure 1",
  "Sepsis produces a broad and robust blood transcriptomic response",
  
  "3.3–3.5",
  "Main Figure 2",
  "SRS and CTS define structured molecular heterogeneity captured by the five-gene score",
  
  "3.6–3.7",
  "Main Figure 3",
  "The host-response axis relates primarily to inflammation and converges with published signatures",
  
  "3.8",
  "Main Figure 4",
  "GSE154918 provides directional external replication but limited primary discrimination",
  
  "3.9",
  "Main Figure 5",
  "GSE185263 externally replicates the association between the five-gene score and organ-dysfunction severity"
)


# ======================================================================
# 7. FILE INVENTORY
# ======================================================================

all_files <- list.files(
  results_root,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)

all_files <- all_files[
  file.exists(all_files)
]

all_files <- normalizePath(
  all_files,
  winslash = "/",
  mustWork = FALSE
)

file_inventory <- tibble::tibble(
  path = all_files,
  basename = basename(all_files),
  extension = tolower(
    tools::file_ext(all_files)
  ),
  size_bytes = as.numeric(
    file.info(all_files)$size
  ),
  modified = file.info(all_files)$mtime
)

cat(
  "Files detected under blood_endotypes_biomarkers =",
  nrow(file_inventory),
  "\n\n"
)


# ======================================================================
# 8. HELPER FUNCTIONS
# ======================================================================

normalize_slashes <- function(x) {
  
  gsub(
    "\\\\",
    "/",
    x
  )
}


split_script_tokens <- function(x) {
  
  tokens <- unlist(
    strsplit(
      x,
      ";",
      fixed = TRUE
    )
  )
  
  tokens <- trimws(tokens)
  
  tokens[
    nzchar(tokens)
  ]
}


script_path_match <- function(
    path,
    script_token
) {
  
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
    gsub(
      "\\.",
      "\\\\.",
      script_token
    ),
    "(_|/)"
  )
  
  grepl(
    pattern,
    path2,
    ignore.case = TRUE,
    perl = TRUE
  )
}


script_path_match_any <- function(
    path,
    script_tokens
) {
  
  tokens <- split_script_tokens(
    script_tokens
  )
  
  if (length(tokens) == 0) {
    return(FALSE)
  }
  
  any(
    vapply(
      tokens,
      function(token) {
        
        script_path_match(
          path,
          token
        )
      },
      FUN.VALUE = logical(1)
    )
  )
}


script_priority_score <- function(
    path,
    script_tokens
) {
  
  tokens <- split_script_tokens(
    script_tokens
  )
  
  if (length(tokens) == 0) {
    return(0)
  }
  
  for (i in seq_along(tokens)) {
    
    if (
      script_path_match(
        path,
        tokens[i]
      )
    ) {
      
      return(
        50 - ((i - 1) * 5)
      )
    }
  }
  
  0
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
        "provenance",
        "numerical_audit",
        "source_data",
        "sessioninfo",
        "manifest",
        "results_placement",
        "/167_",
        "/168_"
      ),
      collapse = "|"
    ),
    x,
    perl = TRUE
  )
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


safe_size <- function(path) {
  
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


# ======================================================================
# 9. FIND MAIN TABLE 1
# ======================================================================

find_main_table1 <- function() {
  
  candidates <- file_inventory$path[
    file_inventory$extension == "xlsx" &
      vapply(
        file_inventory$path,
        script_path_match_any,
        FUN.VALUE = logical(1),
        script_tokens = "151"
      )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  file_names <- basename(candidates)
  
  score <- rep(
    0,
    length(candidates)
  )
  
  score <- score +
    vapply(
      candidates,
      script_priority_score,
      FUN.VALUE = numeric(1),
      script_tokens = "151"
    )
  
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
      grepl(
        "main.*table.?1|table.?1.*main",
        file_names,
        ignore.case = TRUE,
        perl = TRUE
      ),
      50,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "table.?1",
        file_names,
        ignore.case = TRUE,
        perl = TRUE
      ),
      20,
      0
    )
  
  score <- score -
    ifelse(
      grepl(
        "table.?s1",
        file_names,
        ignore.case = TRUE,
        perl = TRUE
      ),
      70,
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
  
  cat("MAIN TABLE 1 CANDIDATES\n")
  cat("-----------------------\n")
  
  print(
    tibble::tibble(
      path = candidates,
      score = score
    )[
      ord,
      ,
      drop = FALSE
    ],
    n = Inf,
    width = Inf
  )
  
  cat("\n")
  
  candidates[
    ord[1]
  ]
}


main_table_path <- find_main_table1()

main_table_found <- !is.na(main_table_path) &&
  file.exists(main_table_path)

cat("Selected Main Table 1:\n")
print(main_table_path)
cat("\n")


# ======================================================================
# 10. FIGURE FILE FINDER
# ======================================================================

find_figure_file <- function(
    figure_number,
    source_scripts,
    extension
) {
  
  candidates <- file_inventory$path[
    file_inventory$extension == extension &
      vapply(
        file_inventory$path,
        script_path_match_any,
        FUN.VALUE = logical(1),
        script_tokens = source_scripts
      )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  file_names <- basename(
    candidates
  )
  
  score <- vapply(
    candidates,
    script_priority_score,
    FUN.VALUE = numeric(1),
    script_tokens = source_scripts
  )
  
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
  
  figure_pattern <- paste0(
    "figure[_ ]?",
    figure_number,
    "(?![0-9])"
  )
  
  score <- score +
    ifelse(
      grepl(
        figure_pattern,
        file_names,
        ignore.case = TRUE,
        perl = TRUE
      ),
      35,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "final",
        file_names,
        ignore.case = TRUE
      ),
      5,
      0
    )
  
  score <- score -
    ifelse(
      grepl(
        "supplement",
        file_names,
        ignore.case = TRUE
      ),
      50,
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
  
  candidates[
    ord[1]
  ]
}


# ======================================================================
# 11. CAPTION FINDER
# ======================================================================

find_figure_caption <- function(
    figure_number,
    source_scripts
) {
  
  candidates <- file_inventory$path[
    file_inventory$extension == "txt" &
      vapply(
        file_inventory$path,
        script_path_match_any,
        FUN.VALUE = logical(1),
        script_tokens = source_scripts
      )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  file_names <- basename(
    candidates
  )
  
  score <- vapply(
    candidates,
    script_priority_score,
    FUN.VALUE = numeric(1),
    script_tokens = source_scripts
  )
  
  score <- score +
    ifelse(
      grepl(
        "caption",
        file_names,
        ignore.case = TRUE
      ),
      40,
      0
    )
  
  figure_pattern <- paste0(
    "figure[_ ]?",
    figure_number,
    "(?![0-9])"
  )
  
  score <- score +
    ifelse(
      grepl(
        figure_pattern,
        file_names,
        ignore.case = TRUE,
        perl = TRUE
      ),
      25,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "caption[_ ]?en\\.txt$",
        file_names,
        ignore.case = TRUE,
        perl = TRUE
      ),
      20,
      0
    )
  
  score <- score +
    ifelse(
      grepl(
        "_en\\.txt$",
        file_names,
        ignore.case = TRUE,
        perl = TRUE
      ),
      10,
      0
    )
  
  score <- score -
    ifelse(
      grepl(
        "results",
        file_names,
        ignore.case = TRUE
      ),
      30,
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
  
  if (score[ord[1]] <= 0) {
    return(NA_character_)
  }
  
  candidates[
    ord[1]
  ]
}


# ======================================================================
# 12. LOCATE MAIN FIGURES
# ======================================================================

figure_manifest <- expected_figures

figure_manifest$png_path <- NA_character_
figure_manifest$pdf_path <- NA_character_
figure_manifest$tiff_path <- NA_character_
figure_manifest$caption_EN_path <- NA_character_

for (i in seq_len(nrow(figure_manifest))) {
  
  n <- figure_manifest$figure_number[i]
  scripts <- figure_manifest$source_scripts[i]
  
  figure_manifest$png_path[i] <- find_figure_file(
    figure_number = n,
    source_scripts = scripts,
    extension = "png"
  )
  
  figure_manifest$pdf_path[i] <- find_figure_file(
    figure_number = n,
    source_scripts = scripts,
    extension = "pdf"
  )
  
  figure_manifest$tiff_path[i] <- find_figure_file(
    figure_number = n,
    source_scripts = scripts,
    extension = "tiff"
  )
  
  if (is.na(figure_manifest$tiff_path[i])) {
    
    figure_manifest$tiff_path[i] <- find_figure_file(
      figure_number = n,
      source_scripts = scripts,
      extension = "tif"
    )
  }
  
  figure_manifest$caption_EN_path[i] <- find_figure_caption(
    figure_number = n,
    source_scripts = scripts
  )
}


figure_manifest <- figure_manifest %>%
  dplyr::mutate(
    png_found =
      !is.na(png_path) &
      file.exists(png_path),
    
    pdf_found =
      !is.na(pdf_path) &
      file.exists(pdf_path),
    
    tiff_found =
      !is.na(tiff_path) &
      file.exists(tiff_path),
    
    caption_EN_found =
      !is.na(caption_EN_path) &
      file.exists(caption_EN_path)
  )


# ======================================================================
# 13. MAIN TABLE READABILITY AUDIT
# ======================================================================

main_table_readable <- FALSE
main_table_sheets <- character()
main_table_error <- NA_character_

if (main_table_found) {
  
  table_test <- tryCatch(
    {
      
      main_table_sheets <- readxl::excel_sheets(
        main_table_path
      )
      
      if (length(main_table_sheets) > 0) {
        
        invisible(
          readxl::read_excel(
            main_table_path,
            sheet = main_table_sheets[1],
            n_max = 3
          )
        )
      }
      
      TRUE
    },
    error = function(e) {
      
      main_table_error <<- conditionMessage(e)
      
      FALSE
    }
  )
  
  main_table_readable <- table_test
}


table_manifest <- expected_table %>%
  dplyr::mutate(
    file_path = main_table_path,
    file_found = main_table_found,
    readable = main_table_readable,
    n_sheets = length(main_table_sheets),
    sheet_names = paste(
      main_table_sheets,
      collapse = " | "
    ),
    error = main_table_error
  )


# ======================================================================
# 14. DIRECT IDENTIFIER AUDIT FOR MAIN TABLE 1
# ======================================================================

identifier_patterns <- c(
  "^patient[_ ]?name$",
  "^participant[_ ]?name$",
  "^full[_ ]?name$",
  "^first[_ ]?name$",
  "^last[_ ]?name$",
  "^surname$",
  "^mrn$",
  "medical[_ ]?record",
  "hospital[_ ]?record",
  "^iin$",
  "national[_ ]?id",
  "passport",
  "^telephone$",
  "^phone$",
  "^email$",
  "street[_ ]?address",
  "home[_ ]?address"
)

identifier_hits <- character()

if (
  main_table_found &&
  main_table_readable
) {
  
  for (sheet in main_table_sheets) {
    
    dat <- tryCatch(
      suppressMessages(
        readxl::read_excel(
          main_table_path,
          sheet = sheet,
          n_max = 2,
          .name_repair = "unique_quiet"
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(dat)) {
      next
    }
    
    cols <- names(dat)
    
    for (pattern in identifier_patterns) {
      
      matched <- cols[
        grepl(
          pattern,
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
}

identifier_hits <- unique(
  identifier_hits
)

identifier_warning_count <- length(
  identifier_hits
)


# ======================================================================
# 15. SOURCE PROVENANCE AUDIT
# ======================================================================

main_table_source_match <- FALSE

if (main_table_found) {
  
  main_table_source_match <- script_path_match_any(
    main_table_path,
    expected_table$source_scripts
  )
}


figure_source_match <- vapply(
  seq_len(nrow(figure_manifest)),
  function(i) {
    
    path <- figure_manifest$png_path[i]
    
    if (
      is.na(path) ||
      !file.exists(path)
    ) {
      return(FALSE)
    }
    
    script_path_match_any(
      path,
      figure_manifest$source_scripts[i]
    )
  },
  FUN.VALUE = logical(1)
)


source_provenance <- dplyr::bind_rows(
  
  tibble::tibble(
    artifact_id = "Main Table 1",
    artifact_type = "Table",
    expected_source_scripts = "151",
    selected_path = main_table_path,
    provenance_match = main_table_source_match
  ),
  
  figure_manifest %>%
    dplyr::transmute(
      artifact_id,
      artifact_type = "Figure",
      expected_source_scripts = source_scripts,
      selected_path = png_path,
      provenance_match = figure_source_match
    )
)

source_provenance_pass <- all(
  source_provenance$provenance_match
)


# ======================================================================
# 16. FIGURE CAPTION LABEL AUDIT
# ======================================================================

caption_audit <- tibble::tibble(
  artifact_id = character(),
  caption_path = character(),
  expected_label = character(),
  label_found = logical(),
  status = character()
)


for (i in seq_len(nrow(figure_manifest))) {
  
  artifact <- figure_manifest$artifact_id[i]
  n <- figure_manifest$figure_number[i]
  path <- figure_manifest$caption_EN_path[i]
  
  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    
    caption_audit <- dplyr::bind_rows(
      caption_audit,
      tibble::tibble(
        artifact_id = artifact,
        caption_path = path,
        expected_label = paste0(
          "Figure ",
          n
        ),
        label_found = FALSE,
        status = "MISSING"
      )
    )
    
    next
  }
  
  caption_text <- paste(
    readLines(
      path,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    collapse = "\n"
  )
  
  pattern <- paste0(
    "(^|\\b)",
    "figure\\s*",
    n,
    "(?![0-9])"
  )
  
  label_found <- grepl(
    pattern,
    caption_text,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  caption_audit <- dplyr::bind_rows(
    caption_audit,
    tibble::tibble(
      artifact_id = artifact,
      caption_path = path,
      expected_label = paste0(
        "Figure ",
        n
      ),
      label_found = label_found,
      status = ifelse(
        label_found,
        "PASS",
        "REVIEW_NOTE"
      )
    )
  )
}


# ======================================================================
# 17. FIGURE COMPLETENESS STATUS
# ======================================================================

figure_manifest <- figure_manifest %>%
  dplyr::mutate(
    formats_complete =
      png_found &
      pdf_found &
      tiff_found,
    
    availability_status = dplyr::case_when(
      !png_found &
        !pdf_found &
        !tiff_found ~ "MISSING",
      
      !formats_complete ~ "ERROR",
      
      !caption_EN_found ~ "ERROR",
      
      TRUE ~ "READY"
    )
  )


# ======================================================================
# 18. MD5 FREEZE MANIFEST
# ======================================================================

freeze_manifest <- tibble::tibble(
  artifact_id = character(),
  artifact_type = character(),
  format = character(),
  source_scripts = character(),
  path = character(),
  file_name = character(),
  size_bytes = numeric(),
  modified = as.POSIXct(character()),
  md5 = character()
)


# Main Table 1 ---------------------------------------------------------

freeze_manifest <- dplyr::bind_rows(
  freeze_manifest,
  tibble::tibble(
    artifact_id = "Main Table 1",
    artifact_type = "Main Table",
    format = "XLSX",
    source_scripts = "151",
    path = main_table_path,
    file_name = ifelse(
      is.na(main_table_path),
      NA_character_,
      basename(main_table_path)
    ),
    size_bytes = safe_size(
      main_table_path
    ),
    modified = safe_modified(
      main_table_path
    ),
    md5 = safe_md5(
      main_table_path
    )
  )
)


# Main Figures ---------------------------------------------------------

for (i in seq_len(nrow(figure_manifest))) {
  
  artifact <- figure_manifest$artifact_id[i]
  scripts <- figure_manifest$source_scripts[i]
  
  paths <- c(
    PNG = figure_manifest$png_path[i],
    PDF = figure_manifest$pdf_path[i],
    TIFF = figure_manifest$tiff_path[i],
    Caption_EN = figure_manifest$caption_EN_path[i]
  )
  
  for (fmt in names(paths)) {
    
    path <- paths[[fmt]]
    
    freeze_manifest <- dplyr::bind_rows(
      freeze_manifest,
      tibble::tibble(
        artifact_id = artifact,
        artifact_type = "Main Figure",
        format = fmt,
        source_scripts = scripts,
        path = path,
        file_name = ifelse(
          is.na(path),
          NA_character_,
          basename(path)
        ),
        size_bytes = safe_size(
          path
        ),
        modified = safe_modified(
          path
        ),
        md5 = safe_md5(
          path
        )
      )
    )
  }
}

freeze_manifest <- freeze_manifest %>%
  dplyr::mutate(
    file_exists =
      !is.na(path) &
      file.exists(path)
  )


# ======================================================================
# 19. CRITICAL AUDIT COUNTS
# ======================================================================

missing_main_table <- as.integer(
  !main_table_found
)

unreadable_main_table <- as.integer(
  main_table_found &&
    !main_table_readable
)

missing_png <- sum(
  !figure_manifest$png_found
)

missing_pdf <- sum(
  !figure_manifest$pdf_found
)

missing_tiff <- sum(
  !figure_manifest$tiff_found
)

missing_captions <- sum(
  !figure_manifest$caption_EN_found
)

caption_label_notes <- sum(
  caption_audit$status == "REVIEW_NOTE"
)

provenance_failures <- sum(
  !source_provenance$provenance_match
)


# ======================================================================
# 20. CRITICAL FREEZE LOGIC
# ======================================================================

critical_failure <- any(
  c(
    !supplementary_lock_found,
    !supplementary_manifest_found,
    !supplementary_frozen_ready,
    missing_main_table > 0,
    unreadable_main_table > 0,
    missing_png > 0,
    missing_pdf > 0,
    missing_tiff > 0,
    missing_captions > 0,
    provenance_failures > 0,
    identifier_warning_count > 0
  )
)

overall_status <- ifelse(
  critical_failure,
  "NOT_FROZEN_CRITICAL_ISSUES",
  "FROZEN_READY"
)


# ======================================================================
# 21. AUDIT SUMMARY
# ======================================================================

audit_summary <- tibble::tribble(
  ~audit_item, ~expected, ~observed, ~status,
  
  "Supplementary package freeze lock",
  "FOUND",
  as.character(supplementary_lock_found),
  ifelse(
    supplementary_lock_found,
    "PASS",
    "FAIL"
  ),
  
  "Supplementary package FROZEN_READY",
  "TRUE",
  as.character(supplementary_frozen_ready),
  ifelse(
    supplementary_frozen_ready,
    "PASS",
    "FAIL"
  ),
  
  "Main Table 1",
  "1",
  as.character(
    as.integer(main_table_found)
  ),
  ifelse(
    main_table_found,
    "PASS",
    "FAIL"
  ),
  
  "Readable Main Table 1",
  "1",
  as.character(
    as.integer(main_table_readable)
  ),
  ifelse(
    main_table_readable,
    "PASS",
    "FAIL"
  ),
  
  "Main Figures PNG",
  "5",
  as.character(
    sum(figure_manifest$png_found)
  ),
  ifelse(
    missing_png == 0,
    "PASS",
    "FAIL"
  ),
  
  "Main Figures PDF",
  "5",
  as.character(
    sum(figure_manifest$pdf_found)
  ),
  ifelse(
    missing_pdf == 0,
    "PASS",
    "FAIL"
  ),
  
  "Main Figures TIFF",
  "5",
  as.character(
    sum(figure_manifest$tiff_found)
  ),
  ifelse(
    missing_tiff == 0,
    "PASS",
    "FAIL"
  ),
  
  "English main-figure captions",
  "5",
  as.character(
    sum(figure_manifest$caption_EN_found)
  ),
  ifelse(
    missing_captions == 0,
    "PASS",
    "FAIL"
  ),
  
  "Expected source-script provenance",
  "PASS",
  as.character(
    source_provenance_pass
  ),
  ifelse(
    source_provenance_pass,
    "PASS",
    "FAIL"
  ),
  
  "Possible direct identifiers in Main Table 1",
  "0",
  as.character(
    identifier_warning_count
  ),
  ifelse(
    identifier_warning_count == 0,
    "PASS",
    "FAIL"
  ),
  
  "Caption-label review notes",
  "manual review only",
  as.character(
    caption_label_notes
  ),
  ifelse(
    caption_label_notes == 0,
    "CLEAR",
    "NOTE"
  )
)


# ======================================================================
# 22. MANUSCRIPT NUMERICAL GUARDRAILS
#
# These are reporting anchors only.
# No statistical calculation is performed here.
# ======================================================================

numerical_guardrails <- tibble::tribble(
  ~domain,
  ~metric,
  ~frozen_value,
  ~reporting_note,
  
  "Discovery cohort",
  "Sepsis blood samples",
  "35",
  "BP",
  
  "Discovery cohort",
  "Healthy blood controls",
  "10",
  "BC",
  
  "Primary DE",
  "Differentially expressed genes",
  "2659",
  "1660 up; 999 down",
  
  "Batch-adjusted DE",
  "Differentially expressed genes",
  "4125",
  "2093 up; 2032 down",
  
  "Robust core",
  "Genes",
  "1796",
  "1133 up; 663 down",
  
  "SRS",
  "BP SRS distribution",
  "SRS1=28; SRS2=7; SRS3=0",
  "Primary SepstratifieR assignment",
  
  "CTS",
  "BP CTS distribution",
  "CTS1=14; CTS2=6; CTS3=15",
  "Primary BP-only CTS assignment",
  
  "Five-gene score",
  "Primary panel",
  "CD177; HK3; IRAK3; CARD11; IKZF2",
  "Biology-guided primary configuration",
  
  "Five-gene score",
  "Score vs SRSq",
  "rho=0.764986",
  "Discovery blood cohort",
  
  "Five-gene score",
  "Score across CTS",
  "Kruskal-Wallis P=9.437e-06",
  "Discovery blood cohort",
  
  "Clinical context",
  "Score vs CRP",
  "rho=0.574350; q=0.018512",
  "Only globally BH-significant score-clinical association",
  
  "GSE154918",
  "Primary infection-control AUC",
  "0.656",
  "Primary endpoint non-significant; P=0.1074",
  
  "GSE185263",
  "Score vs continuous 24-h SOFA",
  "rho=0.311496; P=3.369e-09",
  "Prespecified primary external severity endpoint",
  
  "GSE185263",
  "SOFA >=2 vs 0-1 AUC",
  "0.669887",
  "Secondary binary severity context",
  
  "GSE185263",
  "Adjusted SOFA beta",
  "0.116609",
  "Adjusted for age, sex, and collection location",
  
  "GSE185263",
  "Component genes with expected SOFA direction",
  "5/5",
  "All five BH-significant",
  
  "GSE185263",
  "Mortality AUC",
  "0.627264",
  "Secondary analysis",
  
  "GSE185263",
  "ICU vs Emergency Room AUC",
  "0.643866",
  "Secondary analysis",
  
  "GSE185263",
  "Sepsis vs healthy AUC",
  "0.949451",
  "Contextual only"
)


# ======================================================================
# 23. REPORTING GUARDRAILS
# ======================================================================

reporting_guardrails <- tibble::tibble(
  rule_id = sprintf(
    "RG%02d",
    1:18
  ),
  
  rule = c(
    "The manuscript concerns blood transcriptomics only.",
    "The primary discovery cohort contains 35 sepsis blood samples and 10 healthy blood controls.",
    "The own-cohort platform should be described as targeted whole-blood transcriptomic profiling or targeted RNA sequencing.",
    "The five-gene panel comprises CD177, HK3, IRAK3, CARD11, and IKZF2.",
    "The five-gene panel is biology-guided rather than uniquely optimized by internal discrimination.",
    "Internal discrimination is highly saturated across eligible candidate panels.",
    "SRS and CTS were not used as feature-selection criteria for the primary panel.",
    "The five-gene score represents a molecular host-response state rather than a calibrated clinical assay.",
    "SRSq should be described as a continuous SepstratifieR output, not as a probability.",
    "GSE154918 is an independent external dataset.",
    "The primary GSE154918 infection-control endpoint was not statistically significant.",
    "GSE154918 therefore supports directional replication rather than strong primary clinical discrimination.",
    "GSE185263 is a separate independent external dataset.",
    "The primary GSE185263 endpoint is association with continuous 24-h SOFA among sepsis samples.",
    "GSE185263 collection locations are sensitivity strata within one dataset and not independent validation cohorts.",
    "Mortality and ICU-versus-Emergency-Room analyses in GSE185263 are secondary.",
    "GSE185263 sepsis-versus-healthy discrimination is contextual only.",
    "No main manuscript artifact should be modified after Script 168 freeze without rerunning Script 168."
  )
)


# ======================================================================
# 24. FREEZE DECLARATION
# ======================================================================

freeze_time <- Sys.time()

freeze_declaration <- tibble::tibble(
  item = c(
    "Freeze status",
    "Freeze timestamp",
    "Project",
    "Manuscript scope",
    "Main Table",
    "Main Figures",
    "Supplementary dependency",
    "New inferential analyses in Script 168",
    "Hash algorithm"
  ),
  
  value = c(
    overall_status,
    as.character(freeze_time),
    "Sepsis_DESeq2 blood endotypes and biomarkers",
    "Blood transcriptomics only",
    "Main Table 1",
    "Main Figures 1-5",
    ifelse(
      supplementary_frozen_ready,
      "Script 167 FROZEN_READY",
      "NOT CONFIRMED"
    ),
    "NO",
    "MD5"
  )
)


# ======================================================================
# 25. SAVE MAIN PACKAGE MANIFEST
# ======================================================================

manifest_xlsx <- file.path(
  tables_dir,
  "168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx"
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
  "02_Main_Table1"
)

openxlsx::writeData(
  wb,
  "02_Main_Table1",
  table_manifest
)


openxlsx::addWorksheet(
  wb,
  "03_Main_Figures1_5"
)

openxlsx::writeData(
  wb,
  "03_Main_Figures1_5",
  figure_manifest
)


openxlsx::addWorksheet(
  wb,
  "04_Source_provenance"
)

openxlsx::writeData(
  wb,
  "04_Source_provenance",
  source_provenance
)


openxlsx::addWorksheet(
  wb,
  "05_Caption_audit"
)

openxlsx::writeData(
  wb,
  "05_Caption_audit",
  caption_audit
)


openxlsx::addWorksheet(
  wb,
  "06_MD5_manifest"
)

openxlsx::writeData(
  wb,
  "06_MD5_manifest",
  freeze_manifest
)


openxlsx::addWorksheet(
  wb,
  "07_Numerical_guardrails"
)

openxlsx::writeData(
  wb,
  "07_Numerical_guardrails",
  numerical_guardrails
)


openxlsx::addWorksheet(
  wb,
  "08_Story_map"
)

openxlsx::writeData(
  wb,
  "08_Story_map",
  story_map
)


openxlsx::addWorksheet(
  wb,
  "09_Reporting_guardrails"
)

openxlsx::writeData(
  wb,
  "09_Reporting_guardrails",
  reporting_guardrails
)


# ======================================================================
# 26. FORMAT WORKBOOK
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
# 27. WRITE AUDIT CSV FILES
# ======================================================================

write.csv(
  audit_summary,
  file.path(
    audit_dir,
    "168_main_manuscript_package_audit_summary.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  source_provenance,
  file.path(
    audit_dir,
    "168_main_manuscript_source_provenance.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  freeze_manifest,
  file.path(
    audit_dir,
    "168_main_manuscript_MD5_manifest.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  numerical_guardrails,
  file.path(
    audit_dir,
    "168_main_manuscript_numerical_guardrails.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ======================================================================
# 28. FREEZE REPORT
# ======================================================================

freeze_report <- file.path(
  text_dir,
  "168_COMPLETE_Main_Manuscript_Package_Freeze_Report.txt"
)

sink(
  freeze_report
)

cat("=====================================================================\n")
cat("COMPLETE MAIN MANUSCRIPT PACKAGE FREEZE REPORT\n")
cat("Script 168 FINAL\n")
cat("=====================================================================\n\n")

cat("Project:\n")
cat(project_dir, "\n\n")

cat("Timestamp:\n")
cat(as.character(freeze_time), "\n\n")

cat("STATUS\n")
cat("------\n")
cat(overall_status, "\n\n")

cat("MAIN MANUSCRIPT ARTIFACTS\n")
cat("-------------------------\n")
cat("Main Table 1\n")
cat("Main Figures 1-5\n\n")

cat("SUPPLEMENTARY PACKAGE\n")
cat("---------------------\n")
cat(
  "Script 167 FROZEN_READY =",
  supplementary_frozen_ready,
  "\n\n"
)

cat("AUDIT SUMMARY\n")
cat("-------------\n")
print(
  audit_summary,
  row.names = FALSE
)

cat("\n")

cat("NUMERICAL REPORTING GUARDRAILS\n")
cat("------------------------------\n")
print(
  numerical_guardrails,
  row.names = FALSE
)

cat("\n")

cat("INTERPRETATION\n")
cat("--------------\n")

if (overall_status == "FROZEN_READY") {
  
  cat(
    paste(
      "The complete publication-facing main manuscript package passed",
      "the critical Script 168 integrity audit.",
      "Main Table 1 and Main Figures 1-5 are linked to the expected",
      "frozen source scripts.",
      "Required figure formats and English captions are present.",
      "The previously frozen Supplementary Package was confirmed.",
      "MD5 hashes define the frozen state.",
      "No inferential statistical analysis was performed."
    ),
    "\n"
  )
  
} else {
  
  cat(
    paste(
      "At least one critical manuscript-package integrity requirement",
      "failed. Review the audit outputs before final manuscript assembly."
    ),
    "\n"
  )
}

sink()


# ======================================================================
# 29. FREEZE LOCK
# ======================================================================

freeze_lock_file <- file.path(
  output_dir,
  "168_MAIN_MANUSCRIPT_PACKAGE_FREEZE_LOCK.txt"
)

sink(
  freeze_lock_file
)

cat("MAIN MANUSCRIPT PACKAGE FREEZE LOCK\n")
cat("===================================\n\n")

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

cat("Main Table: Table 1\n")
cat("Main Figures: Figures 1-5\n")
cat("Supplementary Package: Tables S1-S10; Figures S1-S8\n")
cat("Scope: blood transcriptomics only\n")
cat("New inferential statistics in Script 168: NO\n")
cat("Integrity hash: MD5\n\n")

cat(
  "Manifest:\n",
  manifest_xlsx,
  "\n\n",
  sep = ""
)

if (overall_status == "FROZEN_READY") {
  
  cat(
    paste(
      "This file records the main manuscript artifact state",
      "that passed Script 168 FINAL.",
      "Any subsequent modification of a frozen main manuscript artifact",
      "will change its MD5 hash and requires rerunning Script 168."
    ),
    "\n"
  )
  
} else {
  
  cat(
    "The main manuscript package is NOT frozen.\n"
  )
}

sink()


# ======================================================================
# 30. STORY MAP TEXT
# ======================================================================

story_map_txt <- file.path(
  text_dir,
  "168_main_manuscript_story_and_figure_map.txt"
)

sink(
  story_map_txt
)

cat("MAIN MANUSCRIPT STORY / ARTIFACT MAP\n")
cat("====================================\n\n")

for (i in seq_len(nrow(story_map))) {
  
  cat(
    "Results ",
    story_map$results_section[i],
    "\n",
    sep = ""
  )
  
  cat(
    "Artifact: ",
    story_map$main_artifact[i],
    "\n",
    sep = ""
  )
  
  cat(
    "Message: ",
    story_map$main_message[i],
    "\n\n",
    sep = ""
  )
}

sink()


# ======================================================================
# 31. SESSION INFO
# ======================================================================

session_file <- file.path(
  logs_dir,
  "168_sessionInfo.txt"
)

sink(
  session_file
)

sessionInfo()

sink()


# ======================================================================
# 32. CONSOLE REPORT
# ======================================================================

cat("\n")
cat("=====================================================================\n")
cat("Script 168 FINAL COMPLETE MAIN MANUSCRIPT PACKAGE AUDIT\n")
cat("=====================================================================\n\n")


cat("MAIN TABLE 1\n")
cat("------------\n")

print(
  table_manifest,
  n = Inf,
  width = Inf
)

cat("\n")


cat("MAIN FIGURES\n")
cat("------------\n")

print(
  figure_manifest %>%
    dplyr::select(
      artifact_id,
      source_scripts,
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


cat("SOURCE PROVENANCE\n")
cat("-----------------\n")

print(
  source_provenance,
  n = Inf,
  width = Inf
)

cat("\n")


cat("CAPTION LABEL AUDIT\n")
cat("-------------------\n")

print(
  caption_audit,
  n = Inf,
  width = Inf
)

cat("\n")


cat("DIRECT IDENTIFIER AUDIT\n")
cat("-----------------------\n")

cat(
  "Possible identifier-column hits =",
  identifier_warning_count,
  "\n"
)

if (identifier_warning_count > 0) {
  
  print(
    identifier_hits
  )
}

cat("\n")


cat("AUDIT SUMMARY\n")
cat("-------------\n")

print(
  audit_summary,
  n = Inf,
  width = Inf
)

cat("\n")


cat("=====================================================================\n")
cat("FINAL MAIN MANUSCRIPT PACKAGE STATUS\n")
cat("=====================================================================\n\n")

cat(
  overall_status,
  "\n\n"
)

if (overall_status == "FROZEN_READY") {
  
  cat("Main Table 1 found and readable.\n")
  cat("Main Figures 1-5 found in PNG/PDF/TIFF.\n")
  cat("English captions found for Main Figures 1-5.\n")
  cat("Expected source-script provenance passed.\n")
  cat("No possible direct identifiers detected in Main Table 1.\n")
  cat("Supplementary package confirmed FROZEN_READY.\n")
  cat("No new inferential statistical analysis performed.\n\n")
  
  cat(
    "COMPLETE MAIN MANUSCRIPT PACKAGE = FROZEN_READY\n\n"
  )
  
} else {
  
  cat(
    "Critical manuscript-package integrity issues remain.\n"
  )
  
  cat(
    "Do not proceed to final manuscript assembly until resolved.\n\n"
  )
}


cat("OUTPUT FILES\n")
cat("------------\n")

cat(
  "Main manuscript manifest:\n",
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
  "Story map:\n",
  story_map_txt,
  "\n\n",
  sep = ""
)

cat(
  "Session info:\n",
  session_file,
  "\n\n",
  sep = ""
)


cat("REPORTING GUARDRAILS\n")
cat("--------------------\n")
cat("- Script 168 performs no new inferential statistical analysis.\n")
cat("- Blood transcriptomics only.\n")
cat("- Main Table 1 comes from frozen Script 151 outputs.\n")
cat("- Main Figure 1 comes from Script 147.\n")
cat("- Main Figure 2 comes from Script 145b, with Script 145 allowed only as legacy fallback.\n")
cat("- Main Figure 3 comes from Script 148.\n")
cat("- Main Figure 4 comes from Script 149.\n")
cat("- Main Figure 5 comes from Script 150.\n")
cat("- Supplementary Tables S1-S10 and Figures S1-S8 remain frozen by Script 167.\n")
cat("- Primary five-gene panel: CD177, HK3, IRAK3, CARD11, IKZF2.\n")
cat("- Internal discrimination is saturated; the panel is biology-guided.\n")
cat("- GSE154918 primary infection-control endpoint is not significant.\n")
cat("- GSE185263 primary external endpoint is continuous 24-h SOFA.\n")
cat("- GSE185263 locations are sensitivity strata, not independent cohorts.\n")
cat("- The score is not a calibrated diagnostic or prognostic clinical assay.\n")
cat("- MD5 hashes define the frozen main-manuscript artifact state.\n")

cat("\nDone.\n")