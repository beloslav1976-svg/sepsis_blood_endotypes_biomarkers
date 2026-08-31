################################################################################
# Script 156
# FINAL
#
# Supplementary Table S5
#
# Candidate-gene pool and exhaustive blood-panel screening underlying
# the five-gene host-response signature
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
# Assemble a publication-ready supplementary table documenting:
#
#   1. the biology-guided 13-gene blood candidate pool;
#   2. exhaustive screening of eligible 5-8 gene blood panels;
#   3. the DCAF17-forced sensitivity search;
#   4. recommended configurations from the original Script 126;
#   5. the manuscript-specific designation of:
#
#        PRIMARY:
#        CD177, HK3, IRAK3, CARD11, IKZF2
#
#        ALTERNATIVE SENSITIVITY:
#        CD177, HK3, IRAK3, CARD11, DCAF17
#
#
# IMPORTANT
# ---------
#
# THIS SCRIPT DOES NOT:
#
#   - rerun panel selection;
#   - rerun LOOCV;
#   - rerun ROC analysis;
#   - use SRS or CTS for feature selection;
#   - use external-validation cohorts for feature selection;
#   - modify any previously computed AUC or P value.
#
#
# FROZEN SCRIPT 126 SEARCH
# ------------------------
#
# Blood candidate pool:
#
# UP:
#   CD177
#   HK3
#   IRAK3
#   PFKFB3
#   S100A12
#   MMP9
#
# DOWN:
#   CARD11
#   IKZF2
#   NR1D2
#   P2RY10
#   RPS6
#   ST6GAL1
#   DCAF17
#
#
# Eligibility:
#
#   panel size 5-8 genes
#   >=2 UP genes
#   >=2 DOWN genes
#
#
# Complete eligible blood search:
#
#   5,432 panels
#
# DCAF17-forced search:
#
#   2,707 panels
#
#
# Original Script 126 recommended blood configurations:
#
# IKZF2 configuration:
#   CD177;HK3;IRAK3;CARD11;IKZF2
#
#   apparent oriented AUC = 1.0
#   LOOCV oriented AUC    = 1.0
#   LOOCV Wilcoxon P      = 1.897784e-06
#
#
# DCAF17 configuration:
#   CD177;HK3;IRAK3;CARD11;DCAF17
#
#   apparent oriented AUC = 1.0
#   LOOCV oriented AUC    = 1.0
#   LOOCV Wilcoxon P      = 1.897784e-06
#
#
# MANUSCRIPT ROLE
# ---------------
#
# For the current endotype-focused blood manuscript:
#
#   IKZF2 configuration = PRIMARY five-gene signature
#
#   DCAF17 configuration = ALTERNATIVE sensitivity signature
#
# This manuscript designation does not imply that SRS/CTS were used
# to generate the original candidate panel search.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 156\n")
cat("Supplementary Table S5\n")
cat("Candidate pool + exhaustive blood-panel screening\n")
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
  "156_TableS5_candidate_panel_screening"
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

normalize_gene <- function(x) {
  
  toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
}


normalize_gene_set <- function(x) {
  
  x <- as.character(x)
  
  
  vapply(
    x,
    function(one_set) {
      
      if (
        is.na(one_set) ||
        trimws(one_set) ==
        ""
      ) {
        
        return(
          NA_character_
        )
      }
      
      
      genes <- unlist(
        strsplit(
          one_set,
          ";",
          fixed = TRUE
        )
      )
      
      
      genes <- normalize_gene(
        genes
      )
      
      
      genes <- genes[
        genes !=
          ""
      ]
      
      
      paste(
        sort(
          unique(
            genes
          )
        ),
        collapse = ";"
      )
    },
    character(1)
  )
}


count_genes_in_set <- function(x) {
  
  x <- as.character(x)
  
  
  vapply(
    x,
    function(one_set) {
      
      if (
        is.na(one_set) ||
        trimws(one_set) ==
        ""
      ) {
        
        return(
          NA_integer_
        )
      }
      
      
      genes <- unlist(
        strsplit(
          one_set,
          ";",
          fixed = TRUE
        )
      )
      
      
      genes <- trimws(
        genes
      )
      
      
      genes <- genes[
        genes !=
          ""
      ]
      
      
      length(
        unique(
          genes
        )
      )
    },
    integer(1)
  )
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
        "\nAmbiguous candidates for ",
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
      "\nAvailable columns in object while searching for ",
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


clean_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      x
    )
  )
}


# =============================================================================
# 5. INPUT — ORIGINAL SCRIPT 126 WORKBOOK
# =============================================================================

source_file <- file.path(
  project_dir,
  "results",
  "minimal_qPCR_ddPCR_panel_126",
  "tables",
  "126_minimal_qPCR_ddPCR_panel_selection.xlsx"
)


if (
  !file.exists(
    source_file
  )
) {
  
  stop(
    paste0(
      "Original Script 126 workbook not found:\n",
      source_file
    )
  )
}


cat("\nOriginal Script 126 source:\n")

cat(
  normalizePath(
    source_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


source_sheets <- readxl::excel_sheets(
  source_file
)


cat("\nWorkbook sheets:\n")

print(
  source_sheets
)


required_sheets <- c(
  "01_candidate_pool_table",
  "05_blood_panel_search",
  "06_blood_DCAF17_panel_search",
  "09_recommended_panels",
  "10_final_comparison_summary"
)


missing_sheets <- setdiff(
  required_sheets,
  source_sheets
)


if (
  length(
    missing_sheets
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing required Script 126 sheet(s): ",
      paste(
        missing_sheets,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 6. READ SOURCE TABLES
# =============================================================================

candidate_pool_raw <- readxl::read_excel(
  source_file,
  sheet = "01_candidate_pool_table"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


blood_search <- readxl::read_excel(
  source_file,
  sheet = "05_blood_panel_search"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


dcaf17_search <- readxl::read_excel(
  source_file,
  sheet = "06_blood_DCAF17_panel_search"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


recommended_raw <- readxl::read_excel(
  source_file,
  sheet = "09_recommended_panels"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


comparison_raw <- readxl::read_excel(
  source_file,
  sheet = "10_final_comparison_summary"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat("\nSOURCE DIMENSIONS\n")
cat("-----------------\n")


cat(
  "01_candidate_pool_table = ",
  nrow(
    candidate_pool_raw
  ),
  " x ",
  ncol(
    candidate_pool_raw
  ),
  "\n",
  sep = ""
)


cat(
  "05_blood_panel_search = ",
  nrow(
    blood_search
  ),
  " x ",
  ncol(
    blood_search
  ),
  "\n",
  sep = ""
)


cat(
  "06_blood_DCAF17_panel_search = ",
  nrow(
    dcaf17_search
  ),
  " x ",
  ncol(
    dcaf17_search
  ),
  "\n",
  sep = ""
)


cat(
  "09_recommended_panels = ",
  nrow(
    recommended_raw
  ),
  " x ",
  ncol(
    recommended_raw
  ),
  "\n",
  sep = ""
)


cat(
  "10_final_comparison_summary = ",
  nrow(
    comparison_raw
  ),
  " x ",
  ncol(
    comparison_raw
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 7. HARD SEARCH-SIZE AUDIT
# =============================================================================

if (
  nrow(
    blood_search
  ) !=
  5432
) {
  
  stop(
    paste0(
      "Expected 5,432 exhaustive blood panels; observed ",
      nrow(
        blood_search
      ),
      "."
    )
  )
}


if (
  nrow(
    dcaf17_search
  ) !=
  2707
) {
  
  stop(
    paste0(
      "Expected 2,707 DCAF17-forced panels; observed ",
      nrow(
        dcaf17_search
      ),
      "."
    )
  )
}


cat(
  "\nSearch-size audit passed: 5,432 total eligible panels; ",
  "2,707 DCAF17-forced panels.\n",
  sep = ""
)


# =============================================================================
# 8. BLOOD CANDIDATE POOL
# =============================================================================

expected_candidate_genes <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "PFKFB3",
  "S100A12",
  "MMP9",
  "CARD11",
  "IKZF2",
  "NR1D2",
  "P2RY10",
  "RPS6",
  "ST6GAL1",
  "DCAF17"
)


expected_up <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "PFKFB3",
  "S100A12",
  "MMP9"
)


expected_down <- c(
  "CARD11",
  "IKZF2",
  "NR1D2",
  "P2RY10",
  "RPS6",
  "ST6GAL1",
  "DCAF17"
)


candidate_gene_col <- find_column(
  
  candidate_pool_raw,
  
  exact_candidates = c(
    "gene",
    "Gene"
  ),
  
  regex =
    "^gene$",
  
  label =
    "candidate-pool gene"
)


candidate_direction_col <- find_column(
  
  candidate_pool_raw,
  
  exact_candidates = c(
    "direction",
    "Direction"
  ),
  
  regex =
    "direction",
  
  label =
    "candidate-pool direction"
)


candidate_material_col <- find_column(
  
  candidate_pool_raw,
  
  exact_candidates = c(
    "material",
    "Material"
  ),
  
  regex =
    "^material$",
  
  label =
    "candidate-pool material",
  
  required =
    FALSE
)


candidate_source_col <- find_column(
  
  candidate_pool_raw,
  
  exact_candidates = c(
    "source",
    "Source"
  ),
  
  regex =
    "^source$",
  
  label =
    "candidate-pool source",
  
  required =
    FALSE
)


candidate_pool <- candidate_pool_raw


if (
  !is.na(
    candidate_material_col
  )
) {
  
  candidate_pool <- candidate_pool %>%
    
    dplyr::filter(
      tolower(
        trimws(
          as.character(
            .data[[candidate_material_col]]
          )
        )
      ) ==
        "blood"
    )
  
} else {
  
  candidate_pool <- candidate_pool %>%
    
    dplyr::filter(
      normalize_gene(
        .data[[candidate_gene_col]]
      ) %in%
        expected_candidate_genes
    )
}


candidate_pool <- candidate_pool %>%
  
  dplyr::mutate(
    
    gene =
      normalize_gene(
        .data[[candidate_gene_col]]
      ),
    
    direction =
      tolower(
        trimws(
          as.character(
            .data[[candidate_direction_col]]
          )
        )
      )
  )


candidate_genes_observed <- sort(
  unique(
    candidate_pool$gene
  )
)


candidate_genes_expected <- sort(
  expected_candidate_genes
)


if (
  !identical(
    candidate_genes_observed,
    candidate_genes_expected
  )
) {
  
  cat(
    "\nObserved candidate genes:\n"
  )
  
  print(
    candidate_genes_observed
  )
  
  
  cat(
    "\nExpected candidate genes:\n"
  )
  
  print(
    candidate_genes_expected
  )
  
  
  stop(
    "Blood candidate-pool audit failed."
  )
}


if (
  sum(
    candidate_pool$gene %in%
    expected_up &
    candidate_pool$direction ==
    "up"
  ) !=
  6
) {
  
  stop(
    "Expected six UP blood candidate genes."
  )
}


if (
  sum(
    candidate_pool$gene %in%
    expected_down &
    candidate_pool$direction ==
    "down"
  ) !=
  7
) {
  
  stop(
    "Expected seven DOWN blood candidate genes."
  )
}


cat(
  "\nBlood candidate pool audit passed: 13 genes = 6 UP + 7 DOWN.\n"
)


# =============================================================================
# 9. STANDARDIZE BLOOD PANEL-SEARCH COLUMNS
# =============================================================================

standardize_panel_search <- function(
    data,
    search_name
) {
  
  genes_col <- find_column(
    
    data,
    
    exact_candidates = c(
      "genes",
      "Genes"
    ),
    
    regex =
      "^genes$|gene.*set|panel.*genes",
    
    label =
      paste0(
        search_name,
        " genes"
      )
  )
  
  
  n_genes_col <- find_column(
    
    data,
    
    exact_candidates = c(
      "n_genes",
      "panel_size",
      "n_genes_defined"
    ),
    
    regex =
      "n.*genes|panel.*size",
    
    label =
      paste0(
        search_name,
        " panel size"
      ),
    
    required =
      FALSE
  )
  
  
  n_up_col <- find_column(
    
    data,
    
    exact_candidates = c(
      "n_up",
      "n_up_defined"
    ),
    
    regex =
      "^n_up$|n.*up",
    
    label =
      paste0(
        search_name,
        " n_up"
      ),
    
    required =
      FALSE
  )
  
  
  n_down_col <- find_column(
    
    data,
    
    exact_candidates = c(
      "n_down",
      "n_down_defined"
    ),
    
    regex =
      "^n_down$|n.*down",
    
    label =
      paste0(
        search_name,
        " n_down"
      ),
    
    required =
      FALSE
  )
  
  
  apparent_auc_col <- find_column(
    
    data,
    
    exact_candidates = c(
      "apparent_auc_oriented",
      "auc_oriented",
      "apparent_auc"
    ),
    
    regex =
      "apparent.*auc|auc.*apparent",
    
    label =
      paste0(
        search_name,
        " apparent AUC"
      ),
    
    required =
      FALSE
  )
  
  
  loocv_auc_col <- find_column(
    
    data,
    
    exact_candidates = c(
      "loocv_auc_oriented",
      "LOOCV_AUC",
      "loocv_auc"
    ),
    
    regex =
      "loocv.*auc|auc.*loocv",
    
    label =
      paste0(
        search_name,
        " LOOCV AUC"
      )
  )
  
  
  loocv_p_col <- find_column(
    
    data,
    
    exact_candidates = c(
      "loocv_p_wilcox",
      "LOOCV_p",
      "loocv_p"
    ),
    
    regex =
      "loocv.*p|p.*loocv",
    
    label =
      paste0(
        search_name,
        " LOOCV P"
      )
  )
  
  
  result <- data.frame(
    
    search =
      search_name,
    
    genes =
      as.character(
        data[[genes_col]]
      ),
    
    gene_set_normalized =
      normalize_gene_set(
        data[[genes_col]]
      ),
    
    n_genes =
      if (
        !is.na(
          n_genes_col
        )
      ) {
        
        as.integer(
          clean_numeric(
            data[[n_genes_col]]
          )
        )
        
      } else {
        
        count_genes_in_set(
          data[[genes_col]]
        )
      },
    
    n_up =
      if (
        !is.na(
          n_up_col
        )
      ) {
        
        as.integer(
          clean_numeric(
            data[[n_up_col]]
          )
        )
        
      } else {
        
        NA_integer_
      },
    
    n_down =
      if (
        !is.na(
          n_down_col
        )
      ) {
        
        as.integer(
          clean_numeric(
            data[[n_down_col]]
          )
        )
        
      } else {
        
        NA_integer_
      },
    
    apparent_auc_oriented =
      if (
        !is.na(
          apparent_auc_col
        )
      ) {
        
        clean_numeric(
          data[[apparent_auc_col]]
        )
        
      } else {
        
        NA_real_
      },
    
    loocv_auc_oriented =
      clean_numeric(
        data[[loocv_auc_col]]
      ),
    
    loocv_p_wilcox =
      clean_numeric(
        data[[loocv_p_col]]
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  remaining_columns <- setdiff(
    names(
      data
    ),
    unique(
      na.omit(
        c(
          genes_col,
          n_genes_col,
          n_up_col,
          n_down_col,
          apparent_auc_col,
          loocv_auc_col,
          loocv_p_col
        )
      )
    )
  )
  
  
  if (
    length(
      remaining_columns
    ) >
    0
  ) {
    
    result <- cbind(
      result,
      data[
        remaining_columns
      ]
    )
  }
  
  
  result
}


blood_search_std <- standardize_panel_search(
  blood_search,
  "All eligible blood panels"
)


dcaf17_search_std <- standardize_panel_search(
  dcaf17_search,
  "DCAF17-forced blood panels"
)


# =============================================================================
# 10. PANEL ELIGIBILITY AUDIT
# =============================================================================

if (
  any(
    !blood_search_std$n_genes %in%
    5:8
  )
) {
  
  stop(
    "At least one exhaustive blood panel has size outside 5-8 genes."
  )
}


if (
  any(
    !dcaf17_search_std$n_genes %in%
    5:8
  )
) {
  
  stop(
    "At least one DCAF17-forced blood panel has size outside 5-8 genes."
  )
}


if (
  all(
    !is.na(
      blood_search_std$n_up
    )
  ) &&
  any(
    blood_search_std$n_up <
    2
  )
) {
  
  stop(
    "At least one eligible blood panel contains fewer than two UP genes."
  )
}


if (
  all(
    !is.na(
      blood_search_std$n_down
    )
  ) &&
  any(
    blood_search_std$n_down <
    2
  )
) {
  
  stop(
    "At least one eligible blood panel contains fewer than two DOWN genes."
  )
}


# =============================================================================
# 11. DCAF17-FORCED SEARCH AUDIT
# =============================================================================

dcaf17_presence <- vapply(
  strsplit(
    dcaf17_search_std$gene_set_normalized,
    ";",
    fixed = TRUE
  ),
  function(x) {
    "DCAF17" %in% x
  },
  logical(1)
)


if (
  !all(
    dcaf17_presence
  )
) {
  
  stop(
    "At least one panel in DCAF17-forced search does not contain DCAF17."
  )
}


cat(
  "DCAF17-forced search audit passed: 2,707/2,707 panels contain DCAF17.\n"
)


# =============================================================================
# 12. PANEL-SIZE SUMMARY
# =============================================================================

summarise_search_by_size <- function(
    data
) {
  
  data %>%
    
    dplyr::group_by(
      search,
      n_genes
    ) %>%
    
    dplyr::summarise(
      
      n_panels =
        dplyr::n(),
      
      median_LOOCV_AUC =
        stats::median(
          loocv_auc_oriented,
          na.rm = TRUE
        ),
      
      q1_LOOCV_AUC =
        as.numeric(
          stats::quantile(
            loocv_auc_oriented,
            0.25,
            na.rm = TRUE,
            type = 7
          )
        ),
      
      q3_LOOCV_AUC =
        as.numeric(
          stats::quantile(
            loocv_auc_oriented,
            0.75,
            na.rm = TRUE,
            type = 7
          )
        ),
      
      min_LOOCV_AUC =
        min(
          loocv_auc_oriented,
          na.rm = TRUE
        ),
      
      max_LOOCV_AUC =
        max(
          loocv_auc_oriented,
          na.rm = TRUE
        ),
      
      n_LOOCV_AUC_ge_0_95 =
        sum(
          loocv_auc_oriented >=
            0.95,
          na.rm = TRUE
        ),
      
      n_LOOCV_AUC_eq_1 =
        sum(
          abs(
            loocv_auc_oriented -
              1
          ) <
            1e-12,
          na.rm = TRUE
        ),
      
      .groups =
        "drop"
    )
}


panel_size_summary <- dplyr::bind_rows(
  
  summarise_search_by_size(
    blood_search_std
  ),
  
  summarise_search_by_size(
    dcaf17_search_std
  )
)


# =============================================================================
# 13. READ AND FILTER ORIGINAL RECOMMENDED PANELS
# =============================================================================

recommended_type_col <- find_column(
  
  recommended_raw,
  
  exact_candidates = c(
    "recommendation_type"
  ),
  
  regex =
    "recommendation.*type",
  
  label =
    "recommendation type"
)


recommended_genes_col <- find_column(
  
  recommended_raw,
  
  exact_candidates = c(
    "genes"
  ),
  
  regex =
    "^genes$|panel.*genes",
  
  label =
    "recommended-panel genes"
)


recommended_auc_col <- find_column(
  
  recommended_raw,
  
  exact_candidates = c(
    "loocv_auc_oriented"
  ),
  
  regex =
    "loocv.*auc",
  
  label =
    "recommended-panel LOOCV AUC"
)


recommended_p_col <- find_column(
  
  recommended_raw,
  
  exact_candidates = c(
    "loocv_p_wilcox"
  ),
  
  regex =
    "loocv.*p",
  
  label =
    "recommended-panel LOOCV P"
)


recommended_apparent_col <- find_column(
  
  recommended_raw,
  
  exact_candidates = c(
    "apparent_auc_oriented"
  ),
  
  regex =
    "apparent.*auc",
  
  label =
    "recommended-panel apparent AUC",
  
  required =
    FALSE
)


recommended_blood <- recommended_raw %>%
  
  dplyr::mutate(
    
    gene_set_normalized =
      normalize_gene_set(
        .data[[recommended_genes_col]]
      )
  ) %>%
  
  dplyr::filter(
    grepl(
      "blood",
      .data[[recommended_type_col]],
      ignore.case = TRUE
    )
  )


if (
  nrow(
    recommended_blood
  ) !=
  2
) {
  
  stop(
    paste0(
      "Expected exactly two original recommended blood panels; observed ",
      nrow(
        recommended_blood
      ),
      "."
    )
  )
}


# =============================================================================
# 14. DEFINE CURRENT MANUSCRIPT PANELS
# =============================================================================

primary_genes <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "CARD11",
  "IKZF2"
)


alternative_genes <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "CARD11",
  "DCAF17"
)


primary_gene_set <- paste(
  sort(
    primary_genes
  ),
  collapse = ";"
)


alternative_gene_set <- paste(
  sort(
    alternative_genes
  ),
  collapse = ";"
)


primary_match <- recommended_blood %>%
  
  dplyr::filter(
    gene_set_normalized ==
      primary_gene_set
  )


alternative_match <- recommended_blood %>%
  
  dplyr::filter(
    gene_set_normalized ==
      alternative_gene_set
  )


if (
  nrow(
    primary_match
  ) !=
  1
) {
  
  stop(
    "Could not uniquely recover the IKZF2-containing five-gene panel."
  )
}


if (
  nrow(
    alternative_match
  ) !=
  1
) {
  
  stop(
    "Could not uniquely recover the DCAF17-containing five-gene panel."
  )
}


primary_auc <- clean_numeric(
  primary_match[[recommended_auc_col]]
)


primary_p <- clean_numeric(
  primary_match[[recommended_p_col]]
)


primary_apparent <- if (
  !is.na(
    recommended_apparent_col
  )
) {
  clean_numeric(
    primary_match[[recommended_apparent_col]]
  )
} else {
  NA_real_
}


alternative_auc <- clean_numeric(
  alternative_match[[recommended_auc_col]]
)


alternative_p <- clean_numeric(
  alternative_match[[recommended_p_col]]
)


alternative_apparent <- if (
  !is.na(
    recommended_apparent_col
  )
) {
  clean_numeric(
    alternative_match[[recommended_apparent_col]]
  )
} else {
  NA_real_
}


# =============================================================================
# 15. HARD PANEL-PERFORMANCE AUDIT
# =============================================================================

if (
  abs(
    primary_auc -
    1
  ) >
  1e-12
) {
  
  stop(
    "Primary IKZF2 panel LOOCV AUC is not 1.0."
  )
}


if (
  abs(
    alternative_auc -
    1
  ) >
  1e-12
) {
  
  stop(
    "Alternative DCAF17 panel LOOCV AUC is not 1.0."
  )
}


if (
  abs(
    primary_p -
    1.897784e-06
  ) >
  1e-10
) {
  
  stop(
    paste0(
      "Unexpected primary panel LOOCV Wilcoxon P: ",
      primary_p
    )
  )
}


if (
  abs(
    alternative_p -
    1.897784e-06
  ) >
  1e-10
) {
  
  stop(
    paste0(
      "Unexpected alternative panel LOOCV Wilcoxon P: ",
      alternative_p
    )
  )
}


cat(
  "\nFive-gene panel performance audit passed.\n"
)


# =============================================================================
# 16. CURRENT MANUSCRIPT PANEL DEFINITION
# =============================================================================

manuscript_panel_definition <- dplyr::bind_rows(
  
  data.frame(
    
    manuscript_panel =
      "Primary_5_gene",
    
    manuscript_role =
      "Primary biology-guided host-response signature",
    
    gene =
      c(
        "CD177",
        "HK3",
        "IRAK3"
      ),
    
    expected_direction =
      "UP in sepsis",
    
    biological_arm =
      "Myeloid/neutrophil-associated",
    
    stringsAsFactors =
      FALSE
  ),
  
  data.frame(
    
    manuscript_panel =
      "Primary_5_gene",
    
    manuscript_role =
      "Primary biology-guided host-response signature",
    
    gene =
      c(
        "CARD11",
        "IKZF2"
      ),
    
    expected_direction =
      "DOWN in sepsis",
    
    biological_arm =
      "Adaptive/T-cell-associated",
    
    stringsAsFactors =
      FALSE
  ),
  
  data.frame(
    
    manuscript_panel =
      "DCAF17_5_gene",
    
    manuscript_role =
      "Alternative sensitivity signature",
    
    gene =
      c(
        "CD177",
        "HK3",
        "IRAK3"
      ),
    
    expected_direction =
      "UP in sepsis",
    
    biological_arm =
      "Myeloid/neutrophil-associated",
    
    stringsAsFactors =
      FALSE
  ),
  
  data.frame(
    
    manuscript_panel =
      "DCAF17_5_gene",
    
    manuscript_role =
      "Alternative sensitivity signature",
    
    gene =
      c(
        "CARD11",
        "DCAF17"
      ),
    
    expected_direction =
      "DOWN in sepsis",
    
    biological_arm =
      "Downregulated/supporting arm",
    
    stringsAsFactors =
      FALSE
  )
)


manuscript_panel_summary <- data.frame(
  
  manuscript_panel = c(
    "Primary_5_gene",
    "DCAF17_5_gene"
  ),
  
  manuscript_role = c(
    "Primary biology-guided host-response signature",
    "Alternative sensitivity signature"
  ),
  
  genes = c(
    "CD177;HK3;IRAK3;CARD11;IKZF2",
    "CD177;HK3;IRAK3;CARD11;DCAF17"
  ),
  
  n_genes = c(
    5,
    5
  ),
  
  n_up = c(
    3,
    3
  ),
  
  n_down = c(
    2,
    2
  ),
  
  apparent_auc_oriented = c(
    primary_apparent,
    alternative_apparent
  ),
  
  loocv_auc_oriented = c(
    primary_auc,
    alternative_auc
  ),
  
  loocv_p_wilcox = c(
    primary_p,
    alternative_p
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 17. PANEL SEARCH RANK INFORMATION
# =============================================================================

rank_search <- blood_search_std %>%
  
  dplyr::arrange(
    n_genes,
    dplyr::desc(
      loocv_auc_oriented
    ),
    loocv_p_wilcox
  ) %>%
  
  dplyr::mutate(
    ordered_search_position =
      dplyr::row_number()
  )


primary_search_hit <- rank_search %>%
  
  dplyr::filter(
    gene_set_normalized ==
      primary_gene_set
  )


alternative_search_hit <- rank_search %>%
  
  dplyr::filter(
    gene_set_normalized ==
      alternative_gene_set
  )


if (
  nrow(
    primary_search_hit
  ) !=
  1
) {
  
  stop(
    "Primary five-gene panel not uniquely found in exhaustive 5,432-panel search."
  )
}


if (
  nrow(
    alternative_search_hit
  ) !=
  1
) {
  
  stop(
    "Alternative five-gene panel not uniquely found in exhaustive search."
  )
}


selected_panel_search_audit <- data.frame(
  
  panel = c(
    "Primary_5_gene",
    "DCAF17_5_gene"
  ),
  
  gene_set_normalized = c(
    primary_gene_set,
    alternative_gene_set
  ),
  
  exhaustive_search_position = c(
    primary_search_hit$ordered_search_position[1],
    alternative_search_hit$ordered_search_position[1]
  ),
  
  n_genes = c(
    primary_search_hit$n_genes[1],
    alternative_search_hit$n_genes[1]
  ),
  
  loocv_auc_oriented = c(
    primary_search_hit$loocv_auc_oriented[1],
    alternative_search_hit$loocv_auc_oriented[1]
  ),
  
  loocv_p_wilcox = c(
    primary_search_hit$loocv_p_wilcox[1],
    alternative_search_hit$loocv_p_wilcox[1]
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 18. SELECTION PROVENANCE
# =============================================================================

selection_provenance <- data.frame(
  
  Step = c(
    "Candidate-gene pool",
    "Direction constraint",
    "Panel-size constraint",
    "Exhaustive blood search",
    "DCAF17-forced sensitivity search",
    "Internal performance metric",
    "Original practical recommendation threshold",
    "Original ranking logic",
    "Current manuscript primary configuration",
    "Current manuscript alternative configuration",
    "Endotype independence",
    "External-validation independence",
    "Interpretive boundary"
  ),
  
  Description = c(
    
    paste0(
      "Thirteen biology-guided blood candidates: six genes expected to be ",
      "increased in sepsis and seven expected to be decreased."
    ),
    
    "Eligible panels required at least two UP and two DOWN genes.",
    
    "Eligible blood panels contained 5-8 genes.",
    
    "All 5,432 eligible blood configurations were evaluated.",
    
    paste0(
      "A separate search evaluated 2,707 eligible configurations containing ",
      "DCAF17."
    ),
    
    paste0(
      "The signed multi-gene score was evaluated using apparent discrimination ",
      "and leave-one-out cross-validation within the discovery blood cohort."
    ),
    
    "LOOCV AUC >=0.95 was used as the practical recommendation threshold.",
    
    paste0(
      "Among panels meeting the recommendation criterion, ranking prioritized ",
      "smaller panel size, then higher LOOCV AUC, then lower LOOCV Wilcoxon P."
    ),
    
    paste0(
      "For the present endotype-focused manuscript, the IKZF2-containing ",
      "configuration CD177/HK3/IRAK3/CARD11/IKZF2 was designated as the ",
      "primary five-gene host-response signature."
    ),
    
    paste0(
      "The DCAF17-containing configuration ",
      "CD177/HK3/IRAK3/CARD11/DCAF17 was retained as an alternative ",
      "sensitivity signature."
    ),
    
    paste0(
      "Neither SRS nor CTS assignment was used to generate the original ",
      "candidate-gene pool or exhaustive panel search."
    ),
    
    paste0(
      "External cohorts were not used for gene substitution, feature selection, ",
      "coefficient refitting, cutoff optimization, or direction reversal."
    ),
    
    paste0(
      "Internal apparent and cross-validated discrimination within the discovery ",
      "cohort is not equivalent to independent clinical validation."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 19. TABLE S5 README
# =============================================================================

s5_readme <- data.frame(
  
  Item = c(
    "Title",
    "Scope",
    "Candidate pool",
    "Panel eligibility",
    "Exhaustive search",
    "DCAF17-forced search",
    "Primary manuscript signature",
    "Alternative sensitivity signature",
    "Score orientation",
    "Internal validation",
    "Endotype independence",
    "External-data independence",
    "Important limitation"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S5. Candidate-gene pool and exhaustive blood-panel ",
      "screening underlying the five-gene host-response signature."
    ),
    
    paste0(
      "Blood-only candidate-panel development. Urine candidate-panel results ",
      "from the broader Script 126 project are intentionally excluded."
    ),
    
    "13 blood candidate genes: 6 UP and 7 DOWN candidates.",
    
    "Panel size 5-8 genes with at least 2 UP and 2 DOWN genes.",
    
    "5,432 eligible blood panels were evaluated.",
    
    "2,707 eligible DCAF17-containing panels were evaluated separately.",
    
    "CD177, HK3, IRAK3, CARD11 and IKZF2.",
    
    "CD177, HK3, IRAK3, CARD11 and DCAF17.",
    
    paste0(
      "Signed scores were constructed so that greater values represent the ",
      "myeloid-high/adaptive-low sepsis-associated direction."
    ),
    
    paste0(
      "Apparent and LOOCV metrics are internal discovery-cohort estimates; ",
      "they should not be interpreted as independent validation."
    ),
    
    "SRS and CTS assignments were not used for the original panel search.",
    
    paste0(
      "External cohorts were not used for feature selection, coefficient ",
      "refitting, cutoff optimization or gene substitution."
    ),
    
    paste0(
      "The primary five-gene signature is a molecular host-response index, ",
      "not a clinically validated diagnostic assay."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 20. SOURCE MANIFEST
# =============================================================================

source_manifest <- data.frame(
  
  Component = c(
    "Candidate pool",
    "Complete blood panel search",
    "DCAF17-forced panel search",
    "Original recommended configurations",
    "Original comparison summary"
  ),
  
  Workbook_sheet = c(
    "01_candidate_pool_table",
    "05_blood_panel_search",
    "06_blood_DCAF17_panel_search",
    "09_recommended_panels",
    "10_final_comparison_summary"
  ),
  
  Source_file =
    normalizePath(
      source_file,
      winslash = "\\",
      mustWork = TRUE
    ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 21. AUDIT SUMMARY
# =============================================================================

audit_summary <- data.frame(
  
  metric = c(
    "Blood candidate genes",
    "Candidate UP genes",
    "Candidate DOWN genes",
    "Complete eligible blood panels",
    "DCAF17-forced eligible panels",
    "Primary manuscript genes",
    "Primary panel LOOCV AUC",
    "Primary panel LOOCV P",
    "Alternative panel genes",
    "Alternative panel LOOCV AUC",
    "Alternative panel LOOCV P"
  ),
  
  value = c(
    13,
    6,
    7,
    nrow(
      blood_search_std
    ),
    nrow(
      dcaf17_search_std
    ),
    "CD177;HK3;IRAK3;CARD11;IKZF2",
    primary_auc,
    primary_p,
    "CD177;HK3;IRAK3;CARD11;DCAF17",
    alternative_auc,
    alternative_p
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 22. OUTPUT PATHS
# =============================================================================

submission_file <- file.path(
  tables_dir,
  "156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx"
)


audit_file <- file.path(
  audit_dir,
  "156_INTERNAL_AUDIT_TableS5_panel_selection.xlsx"
)


note_file <- file.path(
  text_dir,
  "156_TableS5_title_and_note_EN.txt"
)


# =============================================================================
# 23. STYLES
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
  
  S5_ReadMe =
    s5_readme,
  
  Candidate_pool =
    candidate_pool,
  
  Panel_size_summary =
    panel_size_summary,
  
  Manuscript_panels =
    manuscript_panel_summary,
  
  Manuscript_panel_genes =
    manuscript_panel_definition,
  
  Recommended_blood =
    recommended_blood,
  
  All_eligible_panels =
    blood_search_std,
  
  DCAF17_forced_panels =
    dcaf17_search_std,
  
  Selection_provenance =
    selection_provenance
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
      "S5_ReadMe"
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
  "S5_ReadMe",
  cols = 1,
  widths = 31
)


openxlsx::setColWidths(
  wb,
  "S5_ReadMe",
  cols = 2,
  widths = 90
)


openxlsx::setColWidths(
  wb,
  "Selection_provenance",
  cols = 2,
  widths = 100
)


openxlsx::setColWidths(
  wb,
  "All_eligible_panels",
  cols = 2:3,
  widths = 55
)


openxlsx::setColWidths(
  wb,
  "DCAF17_forced_panels",
  cols = 2:3,
  widths = 55
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
  
  Audit_summary =
    audit_summary,
  
  Source_manifest =
    source_manifest,
  
  Selected_panel_search =
    selected_panel_search_audit,
  
  Candidate_pool_raw =
    candidate_pool_raw,
  
  Recommended_raw =
    recommended_raw,
  
  Comparison_raw =
    comparison_raw
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
# 26. TABLE NOTE
# =============================================================================

table_note <- c(
  
  paste0(
    "Supplementary Table S5. Candidate-gene pool and exhaustive blood-panel ",
    "screening underlying the five-gene host-response signature."
  ),
  
  "",
  
  paste0(
    "Thirteen biology-guided blood candidate genes were evaluated in eligible ",
    "5-8 gene configurations containing at least two genes expected to be ",
    "increased and at least two expected to be decreased in sepsis. The ",
    "original exhaustive search evaluated 5,432 configurations, with a separate ",
    "2,707-panel search requiring DCAF17. For the present endotype-focused ",
    "manuscript, the CD177/HK3/IRAK3/CARD11/IKZF2 configuration was designated ",
    "as the primary five-gene signature, whereas the corresponding DCAF17 ",
    "configuration was retained as an alternative sensitivity signature. ",
    "SRS and CTS assignments were not used in the original feature-selection ",
    "process, and internal LOOCV performance should not be interpreted as ",
    "independent clinical validation."
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
    "156_sessionInfo.txt"
  )
)


# =============================================================================
# 28. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 156 completed successfully.\n")
cat("====================================================================\n\n")


cat("BLOOD CANDIDATE POOL\n")
cat("--------------------\n")


cat(
  "Candidate genes = ",
  length(
    candidate_genes_observed
  ),
  "\n",
  sep = ""
)


cat(
  "UP candidates = ",
  sum(
    candidate_pool$direction ==
      "up"
  ),
  "\n",
  sep = ""
)


cat(
  "DOWN candidates = ",
  sum(
    candidate_pool$direction ==
      "down"
  ),
  "\n",
  sep = ""
)


cat("\nEXHAUSTIVE PANEL SEARCH\n")
cat("-----------------------\n")


cat(
  "All eligible blood panels = ",
  nrow(
    blood_search_std
  ),
  "\n",
  sep = ""
)


cat(
  "DCAF17-forced blood panels = ",
  nrow(
    dcaf17_search_std
  ),
  "\n",
  sep = ""
)


cat("\nPANEL SIZE SUMMARY\n")
cat("------------------\n")


print(
  panel_size_summary,
  row.names = FALSE
)


cat("\nCURRENT MANUSCRIPT PANELS\n")
cat("-------------------------\n")


print(
  manuscript_panel_summary,
  row.names = FALSE
)


cat("\nSELECTED PANEL SEARCH AUDIT\n")
cat("---------------------------\n")


print(
  selected_panel_search_audit,
  row.names = FALSE
)


cat("\nSELECTION GUARDRAILS\n")
cat("--------------------\n")


cat(
  "- Candidate pool = 13 biology-guided blood genes.\n"
)


cat(
  "- Eligible panel sizes = 5-8 genes.\n"
)


cat(
  "- Eligible panels required >=2 UP and >=2 DOWN genes.\n"
)


cat(
  "- Complete search = 5,432 panels.\n"
)


cat(
  "- DCAF17-forced search = 2,707 panels.\n"
)


cat(
  "- Primary manuscript panel = CD177/HK3/IRAK3/CARD11/IKZF2.\n"
)


cat(
  "- DCAF17 configuration is retained as an alternative sensitivity signature.\n"
)


cat(
  "- SRS and CTS were not used for original feature selection.\n"
)


cat(
  "- External cohorts were not used for feature selection or coefficient refitting.\n"
)


cat(
  "- Internal LOOCV is not independent validation.\n"
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Supplementary Table S5:\n  ",
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