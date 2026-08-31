################################################################################
# Script 152
# FINAL v2 — EXPLICIT FROZEN ENRICHMENT SOURCES
#
# Supplementary Tables S2 and S3
#
# Project:
#   Sepsis_DESeq2
#
# Manuscript:
#   Blood-only sepsis transcriptomic endotypes /
#   five-gene host-response signature
#
#
# TABLE S2
# --------
# Complete differential-expression results for all genes retained after
# the prespecified expression prefilter.
#
#
# TABLE S3
# --------
# Robust-core differential expression and functional enrichment:
#
#   - complete robust core
#   - robust core after predefined sex-linked QC filtering
#   - complete GO Biological Process ORA
#   - complete KEGG ORA
#   - complete WikiPathways ORA
#   - complete Hallmark GSEA
#
#
# IMPORTANT
# ---------
#
# THIS SCRIPT DOES NOT RE-RUN:
#
#   - DESeq2
#   - ORA
#   - GSEA
#
# It only imports previously generated frozen outputs, audits them,
# and assembles manuscript-ready supplementary workbooks.
#
#
# FROZEN DE RESULTS
# -----------------
#
# Original targeted matrix:
#   approximately 20,800 targets
#
# Prespecified expression prefilter:
#   >=10 raw counts in >=3 samples
#
# Frozen DE model-comparison universe:
#   12,400 genes
#
# Primary/simple DESeq2:
#   2,659 DEG
#   1,660 UP
#     999 DOWN
#
# Batch-adjusted DESeq2:
#   4,125 DEG
#   2,093 UP
#   2,032 DOWN
#
# Robust core:
#   1,796 DEG
#   1,133 UP
#     663 DOWN
#
# Expected effect-size concordance:
#   Pearson r    ~ 0.815
#   Spearman rho ~ 0.859
#
#
# FROZEN ORA SOURCES
# ------------------
#
# results/blood_enrichment_core/
#
#   GO_BP_core_UP_blood.csv
#   GO_BP_core_DOWN_blood.csv
#   KEGG_core_UP_blood.csv
#   KEGG_core_DOWN_blood.csv
#   WikiPathways_core_UP_blood.csv
#   WikiPathways_core_DOWN_blood.csv
#
#
# FROZEN HALLMARK SOURCE
# ----------------------
#
# 20b_FINAL_Hallmark_GSEA_blood_all_pathways.csv
#
#
# SEX-LINKED QC RULE
# ------------------
#
# Sex-linked genes were NOT removed before DESeq2.
#
# Table S2:
#   retains them and flags them.
#
# Table S3:
#   contains both:
#
#      Robust_core_full
#      Robust_core_no_sex
#
# Selected robust-core ORA/QC branch:
#   sex-linked genes excluded.
#
# Hallmark GSEA:
#   based on the complete ranked gene list and NOT on the
#   sex-linked-filtered robust core.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 152 FINAL v2\n")
cat("Supplementary Tables S2 and S3\n")
cat("Blood differential expression + functional enrichment\n")
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
  "152_Tables_S2_S3_DE_enrichment"
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
  
  x <- as.character(x)
  
  x <- stringr::str_trim(x)
  
  toupper(x)
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
    "SIG",
    "UP",
    "DOWN"
  )
}


find_exact_file <- function(
    preferred_path,
    basename_target,
    required = TRUE
) {
  
  if (
    !is.null(preferred_path) &&
    file.exists(preferred_path)
  ) {
    
    return(
      preferred_path
    )
  }
  
  
  results_root <- file.path(
    project_dir,
    "results"
  )
  
  
  all_files <- list.files(
    results_root,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  
  
  hits <- all_files[
    basename(all_files) ==
      basename_target
  ]
  
  
  hits <- sort(
    unique(hits)
  )
  
  
  if (length(hits) == 1) {
    
    return(
      hits[1]
    )
  }
  
  
  if (length(hits) > 1) {
    
    cat(
      "\nMultiple copies found for:\n",
      basename_target,
      "\n",
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
    
    
    if (
      !is.null(preferred_path)
    ) {
      
      preferred_normalized <- normalizePath(
        dirname(preferred_path),
        winslash = "/",
        mustWork = FALSE
      )
      
      
      hits_normalized <- normalizePath(
        dirname(hits),
        winslash = "/",
        mustWork = FALSE
      )
      
      
      preferred_hits <- hits[
        hits_normalized ==
          preferred_normalized
      ]
      
      
      if (length(preferred_hits) == 1) {
        
        cat(
          "Using preferred frozen source:\n  ",
          preferred_hits[1],
          "\n",
          sep = ""
        )
        
        
        return(
          preferred_hits[1]
        )
      }
    }
    
    
    stop(
      paste0(
        "Multiple candidate files found for ",
        basename_target,
        ". Resolve provenance before continuing."
      )
    )
  }
  
  
  if (required) {
    
    stop(
      paste0(
        "Required source file not found:\n",
        basename_target,
        "\nPreferred path:\n",
        preferred_path
      )
    )
  }
  
  
  return(
    NA_character_
  )
}


find_column <- function(
    data,
    exact_candidates = character(0),
    regex = NULL,
    label = "column",
    required = TRUE
) {
  
  nm <- names(data)
  
  nm_lower <- tolower(nm)
  
  
  if (
    length(
      exact_candidates
    ) >
    0
  ) {
    
    for (
      candidate in exact_candidates
    ) {
      
      idx <- which(
        nm_lower ==
          tolower(candidate)
      )
      
      
      if (length(idx) == 1) {
        
        return(
          nm[idx]
        )
      }
    }
  }
  
  
  if (!is.null(regex)) {
    
    idx <- grep(
      regex,
      nm,
      ignore.case = TRUE,
      perl = TRUE
    )
    
    
    if (length(idx) == 1) {
      
      return(
        nm[idx]
      )
    }
    
    
    if (length(idx) > 1) {
      
      cat(
        "\nAmbiguous ",
        label,
        " candidates:\n",
        sep = ""
      )
      
      
      print(
        nm[idx]
      )
      
      
      stop(
        paste0(
          "Could not uniquely identify ",
          label,
          "."
        )
      )
    }
  }
  
  
  if (required) {
    
    cat(
      "\nAvailable columns:\n"
    )
    
    print(
      nm
    )
    
    
    stop(
      paste0(
        "Could not identify required ",
        label,
        "."
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
  
  
  if (
    length(
      missing_columns
    ) >
    0
  ) {
    
    cat(
      "\nColumns available in ",
      table_name,
      ":\n",
      sep = ""
    )
    
    
    print(
      names(data)
    )
    
    
    stop(
      paste0(
        "Missing required column(s) in ",
        table_name,
        ": ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
}


sanitize_for_excel <- function(data) {
  
  data <- as.data.frame(
    data,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  
  for (
    column_name in names(data)
  ) {
    
    if (
      is.list(
        data[[column_name]]
      )
    ) {
      
      data[[column_name]] <- vapply(
        data[[column_name]],
        function(x) {
          
          if (length(x) == 0) {
            return("")
          }
          
          
          paste(
            as.character(x),
            collapse = ";"
          )
        },
        character(1)
      )
    }
  }
  
  
  data
}


# =============================================================================
# 5. PRIMARY MODEL-COMPARISON SOURCE
# =============================================================================

model_comparison_file <- find_exact_file(
  
  preferred_path = file.path(
    project_dir,
    "results",
    "blood_BP_vs_BC_model_comparison",
    "blood_model_comparison_all_genes.csv"
  ),
  
  basename_target =
    "blood_model_comparison_all_genes.csv",
  
  required = TRUE
)


cat("\nPrimary model-comparison source:\n")

cat(
  normalizePath(
    model_comparison_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


model_df <- read.csv(
  model_comparison_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat(
  "Model-comparison dimensions: ",
  nrow(model_df),
  " rows x ",
  ncol(model_df),
  " columns\n",
  sep = ""
)


# =============================================================================
# 6. ANALYSIS-UNIVERSE AUDIT
# =============================================================================

n_model_genes <- nrow(
  model_df
)


cat(
  "Genes retained in frozen DE model-comparison universe = ",
  n_model_genes,
  "\n",
  sep = ""
)


if (
  n_model_genes !=
  12400
) {
  
  stop(
    paste0(
      "Frozen model-comparison universe differs from expected 12,400 genes. ",
      "Observed: ",
      n_model_genes,
      ". Re-audit provenance before continuing."
    )
  )
}


cat(
  "Analysis-universe audit passed: 12,400 genes.\n"
)


# =============================================================================
# 7. DETECT MODEL-COMPARISON COLUMNS
# =============================================================================

gene_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "Gene",
    "gene",
    "gene_symbol",
    "GeneSymbol",
    "symbol"
  ),
  
  regex =
    "^(gene|gene_symbol|genesymbol|symbol)$",
  
  label =
    "gene-symbol column"
)


primary_lfc_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "log2FoldChange_simple",
    "log2FC_simple",
    "simple_log2FoldChange",
    "primary_log2FoldChange",
    "log2FoldChange_primary"
  ),
  
  regex =
    "(log2foldchange|log2fc).*(simple|primary)|(simple|primary).*(log2foldchange|log2fc)",
  
  label =
    "primary/simple log2FoldChange"
)


batch_lfc_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "log2FoldChange_batch",
    "log2FC_batch",
    "batch_log2FoldChange",
    "log2FoldChange_batch_adjusted",
    "batch_adjusted_log2FoldChange"
  ),
  
  regex =
    "(log2foldchange|log2fc).*(batch)|(batch).*(log2foldchange|log2fc)",
  
  label =
    "batch-adjusted log2FoldChange"
)


primary_padj_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "padj_simple",
    "simple_padj",
    "padj_primary",
    "primary_padj"
  ),
  
  regex =
    "(padj|fdr|adj.*p).*(simple|primary)|(simple|primary).*(padj|fdr|adj.*p)",
  
  label =
    "primary/simple adjusted P value"
)


batch_padj_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "padj_batch",
    "batch_padj",
    "padj_batch_adjusted",
    "batch_adjusted_padj"
  ),
  
  regex =
    "(padj|fdr|adj.*p).*(batch)|(batch).*(padj|fdr|adj.*p)",
  
  label =
    "batch-adjusted adjusted P value"
)


primary_deg_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "DEG_simple",
    "deg_simple",
    "simple_DEG",
    "DEG_primary",
    "primary_DEG"
  ),
  
  regex =
    "deg.*(simple|primary)|(simple|primary).*deg",
  
  label =
    "primary DEG flag",
  
  required = FALSE
)


batch_deg_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "DEG_batch",
    "deg_batch",
    "batch_DEG",
    "DEG_batch_adjusted",
    "batch_adjusted_DEG"
  ),
  
  regex =
    "deg.*batch|batch.*deg",
  
  label =
    "batch DEG flag",
  
  required = FALSE
)


core_deg_col <- find_column(
  
  model_df,
  
  exact_candidates = c(
    "core_DEG",
    "Core_DEG",
    "core_deg",
    "robust_core_DEG",
    "robust_core"
  ),
  
  regex =
    "core.*deg|robust.*core",
  
  label =
    "robust-core DEG flag",
  
  required = FALSE
)


column_mapping <- data.frame(
  
  role = c(
    "Gene symbol",
    "Primary log2FoldChange",
    "Primary adjusted P",
    "Primary DEG flag",
    "Batch-adjusted log2FoldChange",
    "Batch-adjusted adjusted P",
    "Batch-adjusted DEG flag",
    "Robust-core DEG flag"
  ),
  
  source_column = c(
    gene_col,
    primary_lfc_col,
    primary_padj_col,
    primary_deg_col,
    batch_lfc_col,
    batch_padj_col,
    batch_deg_col,
    core_deg_col
  ),
  
  stringsAsFactors = FALSE
)


cat("\nDetected model-comparison column mapping:\n")

print(
  column_mapping,
  row.names = FALSE
)


# =============================================================================
# 8. STANDARDIZE DE VARIABLES
# =============================================================================

gene_values <- as.character(
  model_df[[gene_col]]
)


primary_lfc <- suppressWarnings(
  as.numeric(
    model_df[[primary_lfc_col]]
  )
)


batch_lfc <- suppressWarnings(
  as.numeric(
    model_df[[batch_lfc_col]]
  )
)


primary_padj <- suppressWarnings(
  as.numeric(
    model_df[[primary_padj_col]]
  )
)


batch_padj <- suppressWarnings(
  as.numeric(
    model_df[[batch_padj_col]]
  )
)


# =============================================================================
# 9. DEG FLAGS
# =============================================================================

primary_deg_calculated <-
  is.finite(primary_padj) &
  primary_padj <
  0.05 &
  is.finite(primary_lfc) &
  abs(primary_lfc) >=
  1


batch_deg_calculated <-
  is.finite(batch_padj) &
  batch_padj <
  0.05 &
  is.finite(batch_lfc) &
  abs(batch_lfc) >=
  1


if (!is.na(primary_deg_col)) {
  
  primary_deg <- as_logical_flag(
    model_df[[primary_deg_col]]
  )
  
  
  primary_flag_mismatch <- sum(
    primary_deg !=
      primary_deg_calculated
  )
  
} else {
  
  primary_deg <-
    primary_deg_calculated
  
  
  primary_flag_mismatch <-
    NA_integer_
}


if (!is.na(batch_deg_col)) {
  
  batch_deg <- as_logical_flag(
    model_df[[batch_deg_col]]
  )
  
  
  batch_flag_mismatch <- sum(
    batch_deg !=
      batch_deg_calculated
  )
  
} else {
  
  batch_deg <-
    batch_deg_calculated
  
  
  batch_flag_mismatch <-
    NA_integer_
}


core_calculated <-
  primary_deg &
  batch_deg &
  is.finite(primary_lfc) &
  is.finite(batch_lfc) &
  sign(primary_lfc) ==
  sign(batch_lfc)


if (!is.na(core_deg_col)) {
  
  robust_core <- as_logical_flag(
    model_df[[core_deg_col]]
  )
  
  
  core_flag_vs_reconstructed_mismatch <- sum(
    robust_core !=
      core_calculated
  )
  
} else {
  
  robust_core <-
    core_calculated
  
  
  core_flag_vs_reconstructed_mismatch <-
    NA_integer_
}


# =============================================================================
# 10. DEG COUNT AUDIT
# =============================================================================

primary_up <-
  primary_deg &
  primary_lfc >
  0


primary_down <-
  primary_deg &
  primary_lfc <
  0


batch_up <-
  batch_deg &
  batch_lfc >
  0


batch_down <-
  batch_deg &
  batch_lfc <
  0


core_up <-
  robust_core &
  primary_lfc >
  0


core_down <-
  robust_core &
  primary_lfc <
  0


de_summary <- data.frame(
  
  analysis = c(
    "Primary/simple",
    "Batch-adjusted",
    "Robust core"
  ),
  
  total_DEG = c(
    sum(primary_deg),
    sum(batch_deg),
    sum(robust_core)
  ),
  
  UP = c(
    sum(primary_up),
    sum(batch_up),
    sum(core_up)
  ),
  
  DOWN = c(
    sum(primary_down),
    sum(batch_down),
    sum(core_down)
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


de_summary$total_match <-
  de_summary$total_DEG ==
  de_summary$expected_total


de_summary$UP_match <-
  de_summary$UP ==
  de_summary$expected_UP


de_summary$DOWN_match <-
  de_summary$DOWN ==
  de_summary$expected_DOWN


cat("\nDifferential-expression count audit:\n")

print(
  de_summary,
  row.names = FALSE
)


if (
  !all(de_summary$total_match) ||
  !all(de_summary$UP_match) ||
  !all(de_summary$DOWN_match)
) {
  
  stop(
    "Differential-expression count audit failed."
  )
}


cat(
  "\nDEG count audit passed successfully.\n"
)


# =============================================================================
# 11. EFFECT-SIZE CONCORDANCE AUDIT
# =============================================================================

valid_lfc <-
  is.finite(primary_lfc) &
  is.finite(batch_lfc)


pearson_lfc <- stats::cor(
  primary_lfc[valid_lfc],
  batch_lfc[valid_lfc],
  method = "pearson"
)


spearman_lfc <- stats::cor(
  primary_lfc[valid_lfc],
  batch_lfc[valid_lfc],
  method = "spearman"
)


concordance_audit <- data.frame(
  
  metric = c(
    "Pearson r",
    "Spearman rho"
  ),
  
  observed = c(
    pearson_lfc,
    spearman_lfc
  ),
  
  expected_approx = c(
    0.815,
    0.859
  ),
  
  difference = c(
    pearson_lfc - 0.815,
    spearman_lfc - 0.859
  ),
  
  stringsAsFactors = FALSE
)


cat("\nEffect-size concordance audit:\n")

print(
  concordance_audit,
  row.names = FALSE
)


if (
  abs(
    pearson_lfc -
    0.815
  ) >
  0.02
) {
  
  stop(
    "Primary vs batch Pearson correlation audit failed."
  )
}


if (
  abs(
    spearman_lfc -
    0.859
  ) >
  0.02
) {
  
  stop(
    "Primary vs batch Spearman correlation audit failed."
  )
}


cat(
  "\nEffect-size concordance audit passed successfully.\n"
)


# =============================================================================
# 12. SEX-LINKED QC SET
# =============================================================================

sex_linked_genes <- c(
  "EIF1AY",
  "KDM5D",
  "DDX3Y",
  "ZFY",
  "TTTY15",
  "USP9Y",
  "RPS4Y1",
  "UTY",
  "SRY",
  "TMSB4Y",
  "NLGN4Y",
  "TXLNGY",
  "XIST",
  "TSIX"
)


gene_upper <- normalize_gene(
  gene_values
)


sex_linked_flag <-
  gene_upper %in%
  sex_linked_genes


sex_linked_audit <- data.frame(
  
  gene =
    sex_linked_genes,
  
  present_in_analysis_universe =
    sex_linked_genes %in%
    gene_upper,
  
  present_in_primary_DEG =
    sex_linked_genes %in%
    gene_upper[
      primary_deg
    ],
  
  present_in_batch_DEG =
    sex_linked_genes %in%
    gene_upper[
      batch_deg
    ],
  
  present_in_robust_core =
    sex_linked_genes %in%
    gene_upper[
      robust_core
    ],
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 13. BUILD TABLE S2
# =============================================================================

standardized_columns <- data.frame(
  
  Gene =
    gene_values,
  
  primary_log2FoldChange =
    primary_lfc,
  
  primary_padj =
    primary_padj,
  
  primary_DEG =
    primary_deg,
  
  batch_adjusted_log2FoldChange =
    batch_lfc,
  
  batch_adjusted_padj =
    batch_padj,
  
  batch_adjusted_DEG =
    batch_deg,
  
  robust_core_DEG =
    robust_core,
  
  robust_core_direction =
    ifelse(
      robust_core &
        primary_lfc >
        0,
      "UP in sepsis",
      ifelse(
        robust_core &
          primary_lfc <
          0,
        "DOWN in sepsis",
        ""
      )
    ),
  
  sex_linked_QC_gene =
    sex_linked_flag,
  
  stringsAsFactors = FALSE
)


mapped_source_columns <- unique(
  na.omit(
    c(
      gene_col,
      primary_lfc_col,
      batch_lfc_col,
      primary_padj_col,
      batch_padj_col,
      primary_deg_col,
      batch_deg_col,
      core_deg_col
    )
  )
)


remaining_source_columns <- setdiff(
  names(model_df),
  mapped_source_columns
)


if (
  length(
    remaining_source_columns
  ) >
  0
) {
  
  tableS2 <- cbind(
    standardized_columns,
    model_df[
      remaining_source_columns
    ]
  )
  
} else {
  
  tableS2 <-
    standardized_columns
}


tableS2 <- sanitize_for_excel(
  tableS2
)


if (
  nrow(tableS2) !=
  12400
) {
  
  stop(
    "Table S2 row-count audit failed."
  )
}


# =============================================================================
# 14. BUILD ROBUST-CORE TABLES
# =============================================================================

tableS3_core_full <- tableS2[
  robust_core,
  ,
  drop = FALSE
]


tableS3_core_no_sex <- tableS2[
  robust_core &
    !sex_linked_flag,
  ,
  drop = FALSE
]


if (
  nrow(
    tableS3_core_full
  ) !=
  1796
) {
  
  stop(
    "Expected Robust_core_full n=1796."
  )
}


tableS3_core_full <- tableS3_core_full %>%
  
  dplyr::arrange(
    factor(
      robust_core_direction,
      levels = c(
        "UP in sepsis",
        "DOWN in sepsis"
      )
    ),
    dplyr::desc(
      abs(
        primary_log2FoldChange
      )
    )
  )


tableS3_core_no_sex <- tableS3_core_no_sex %>%
  
  dplyr::arrange(
    factor(
      robust_core_direction,
      levels = c(
        "UP in sepsis",
        "DOWN in sepsis"
      )
    ),
    dplyr::desc(
      abs(
        primary_log2FoldChange
      )
    )
  )


# =============================================================================
# 15. OPTIONAL AUDIT AGAINST EXISTING NO-SEX ROBUST-CORE CSV FILES
# =============================================================================

sex_filter_dir <- file.path(
  project_dir,
  "results",
  "blood_QC_publication",
  "sex_linked_filter"
)


existing_no_sex_up <- file.path(
  sex_filter_dir,
  "blood_core_UP_no_sex_linked.csv"
)


existing_no_sex_down <- file.path(
  sex_filter_dir,
  "blood_core_DOWN_no_sex_linked.csv"
)


sex_filter_source_audit <- data.frame(
  
  source = c(
    "blood_core_UP_no_sex_linked.csv",
    "blood_core_DOWN_no_sex_linked.csv"
  ),
  
  exists = c(
    file.exists(
      existing_no_sex_up
    ),
    file.exists(
      existing_no_sex_down
    )
  ),
  
  set_match =
    NA,
  
  stringsAsFactors = FALSE
)


extract_gene_set <- function(path) {
  
  tmp <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  
  gcol <- find_column(
    
    tmp,
    
    exact_candidates = c(
      "Gene",
      "gene",
      "gene_symbol",
      "symbol"
    ),
    
    regex =
      "^(gene|gene_symbol|symbol)$",
    
    label =
      paste0(
        "gene column in ",
        basename(path)
      )
  )
  
  
  unique(
    normalize_gene(
      tmp[[gcol]]
    )
  )
}


derived_no_sex_up_genes <- sort(
  unique(
    gene_upper[
      robust_core &
        !sex_linked_flag &
        primary_lfc >
        0
    ]
  )
)


derived_no_sex_down_genes <- sort(
  unique(
    gene_upper[
      robust_core &
        !sex_linked_flag &
        primary_lfc <
        0
    ]
  )
)


if (
  file.exists(
    existing_no_sex_up
  )
) {
  
  existing_up_genes <- sort(
    extract_gene_set(
      existing_no_sex_up
    )
  )
  
  
  sex_filter_source_audit$set_match[1] <-
    setequal(
      existing_up_genes,
      derived_no_sex_up_genes
    )
}


if (
  file.exists(
    existing_no_sex_down
  )
) {
  
  existing_down_genes <- sort(
    extract_gene_set(
      existing_no_sex_down
    )
  )
  
  
  sex_filter_source_audit$set_match[2] <-
    setequal(
      existing_down_genes,
      derived_no_sex_down_genes
    )
}


if (
  any(
    sex_filter_source_audit$exists &
    !sex_filter_source_audit$set_match,
    na.rm = TRUE
  )
) {
  
  print(
    sex_filter_source_audit
  )
  
  
  stop(
    "Existing no-sex robust-core lists do not match reconstructed QC branch."
  )
}


# =============================================================================
# 16. EXPLICIT FROZEN ORA SOURCES
# =============================================================================
#
# IMPORTANT:
#
# These are the canonical source files used by Scripts 16a and 16b.
#
# Do NOT use files from blood_enrichment_publication as the source of
# complete ORA results, because that directory contains publication-selected
# and visualization-oriented derivative files.
#
# =============================================================================

enrichment_core_dir <- file.path(
  project_dir,
  "results",
  "blood_enrichment_core"
)


if (
  !dir.exists(
    enrichment_core_dir
  )
) {
  
  stop(
    paste0(
      "Required frozen enrichment directory not found:\n",
      enrichment_core_dir
    )
  )
}


ora_files <- c(
  
  GO_BP_UP = file.path(
    enrichment_core_dir,
    "GO_BP_core_UP_blood.csv"
  ),
  
  GO_BP_DOWN = file.path(
    enrichment_core_dir,
    "GO_BP_core_DOWN_blood.csv"
  ),
  
  KEGG_UP = file.path(
    enrichment_core_dir,
    "KEGG_core_UP_blood.csv"
  ),
  
  KEGG_DOWN = file.path(
    enrichment_core_dir,
    "KEGG_core_DOWN_blood.csv"
  ),
  
  WikiPathways_UP = file.path(
    enrichment_core_dir,
    "WikiPathways_core_UP_blood.csv"
  ),
  
  WikiPathways_DOWN = file.path(
    enrichment_core_dir,
    "WikiPathways_core_DOWN_blood.csv"
  )
)


missing_ora_files <- ora_files[
  !file.exists(
    ora_files
  )
]


if (
  length(
    missing_ora_files
  ) >
  0
) {
  
  cat(
    "\nMissing frozen ORA file(s):\n"
  )
  
  
  print(
    missing_ora_files
  )
  
  
  stop(
    paste0(
      "One or more canonical blood_enrichment_core files are missing.\n",
      "Do not substitute publication-selected files."
    )
  )
}


cat("\nFrozen ORA source manifest:\n")


ora_manifest <- data.frame(
  
  Analysis =
    names(
      ora_files
    ),
  
  File =
    normalizePath(
      unname(
        ora_files
      ),
      winslash = "\\",
      mustWork = TRUE
    ),
  
  Exists =
    file.exists(
      ora_files
    ),
  
  stringsAsFactors = FALSE
)


print(
  ora_manifest,
  row.names = FALSE
)


# =============================================================================
# 17. READ COMPLETE ORA TABLES
# =============================================================================

go_up <- read.csv(
  ora_files[["GO_BP_UP"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


go_down <- read.csv(
  ora_files[["GO_BP_DOWN"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


kegg_up <- read.csv(
  ora_files[["KEGG_UP"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


kegg_down <- read.csv(
  ora_files[["KEGG_DOWN"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


wiki_up <- read.csv(
  ora_files[["WikiPathways_UP"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


wiki_down <- read.csv(
  ora_files[["WikiPathways_DOWN"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("\nFrozen ORA table dimensions:\n")


cat(
  "GO BP UP = ",
  nrow(go_up),
  " rows\n",
  sep = ""
)


cat(
  "GO BP DOWN = ",
  nrow(go_down),
  " rows\n",
  sep = ""
)


cat(
  "KEGG UP = ",
  nrow(kegg_up),
  " rows\n",
  sep = ""
)


cat(
  "KEGG DOWN = ",
  nrow(kegg_down),
  " rows\n",
  sep = ""
)


cat(
  "WikiPathways UP = ",
  nrow(wiki_up),
  " rows\n",
  sep = ""
)


cat(
  "WikiPathways DOWN = ",
  nrow(wiki_down),
  " rows\n",
  sep = ""
)


# =============================================================================
# 18. ORA STRUCTURE AUDIT
# =============================================================================

minimum_ora_columns <- c(
  "ID",
  "Description",
  "p.adjust",
  "Count"
)


check_required_columns(
  go_up,
  minimum_ora_columns,
  "GO BP UP"
)


check_required_columns(
  go_down,
  minimum_ora_columns,
  "GO BP DOWN"
)


check_required_columns(
  kegg_up,
  minimum_ora_columns,
  "KEGG UP"
)


check_required_columns(
  kegg_down,
  minimum_ora_columns,
  "KEGG DOWN"
)


check_required_columns(
  wiki_up,
  minimum_ora_columns,
  "WikiPathways UP"
)


check_required_columns(
  wiki_down,
  minimum_ora_columns,
  "WikiPathways DOWN"
)


if (
  any(
    c(
      nrow(go_up),
      nrow(go_down),
      nrow(kegg_up),
      nrow(kegg_down),
      nrow(wiki_up),
      nrow(wiki_down)
    ) ==
    0
  )
) {
  
  stop(
    "At least one frozen ORA table contains zero rows."
  )
}


# =============================================================================
# 19. COMBINE ORA TABLES
# =============================================================================

add_ora_metadata <- function(
    data,
    database,
    direction
) {
  
  data.frame(
    
    Database =
      database,
    
    Direction =
      direction,
    
    data,
    
    check.names = FALSE,
    
    stringsAsFactors = FALSE
  )
}


go_table <- dplyr::bind_rows(
  
  add_ora_metadata(
    go_up,
    "GO Biological Process",
    "UP in sepsis"
  ),
  
  add_ora_metadata(
    go_down,
    "GO Biological Process",
    "DOWN in sepsis"
  )
)


kegg_table <- dplyr::bind_rows(
  
  add_ora_metadata(
    kegg_up,
    "KEGG",
    "UP in sepsis"
  ),
  
  add_ora_metadata(
    kegg_down,
    "KEGG",
    "DOWN in sepsis"
  )
)


wiki_table <- dplyr::bind_rows(
  
  add_ora_metadata(
    wiki_up,
    "WikiPathways",
    "UP in sepsis"
  ),
  
  add_ora_metadata(
    wiki_down,
    "WikiPathways",
    "DOWN in sepsis"
  )
)


go_table <- sanitize_for_excel(
  go_table
)


kegg_table <- sanitize_for_excel(
  kegg_table
)


wiki_table <- sanitize_for_excel(
  wiki_table
)


# =============================================================================
# 20. ORA SUMMARY AUDIT
# =============================================================================

ora_summary <- data.frame(
  
  database = c(
    "GO Biological Process",
    "GO Biological Process",
    "KEGG",
    "KEGG",
    "WikiPathways",
    "WikiPathways"
  ),
  
  direction = c(
    "UP in sepsis",
    "DOWN in sepsis",
    "UP in sepsis",
    "DOWN in sepsis",
    "UP in sepsis",
    "DOWN in sepsis"
  ),
  
  n_rows = c(
    nrow(go_up),
    nrow(go_down),
    nrow(kegg_up),
    nrow(kegg_down),
    nrow(wiki_up),
    nrow(wiki_down)
  ),
  
  n_FDR_lt_0_05 = c(
    sum(
      go_up$p.adjust <
        0.05,
      na.rm = TRUE
    ),
    sum(
      go_down$p.adjust <
        0.05,
      na.rm = TRUE
    ),
    sum(
      kegg_up$p.adjust <
        0.05,
      na.rm = TRUE
    ),
    sum(
      kegg_down$p.adjust <
        0.05,
      na.rm = TRUE
    ),
    sum(
      wiki_up$p.adjust <
        0.05,
      na.rm = TRUE
    ),
    sum(
      wiki_down$p.adjust <
        0.05,
      na.rm = TRUE
    )
  ),
  
  min_FDR = c(
    min(
      go_up$p.adjust,
      na.rm = TRUE
    ),
    min(
      go_down$p.adjust,
      na.rm = TRUE
    ),
    min(
      kegg_up$p.adjust,
      na.rm = TRUE
    ),
    min(
      kegg_down$p.adjust,
      na.rm = TRUE
    ),
    min(
      wiki_up$p.adjust,
      na.rm = TRUE
    ),
    min(
      wiki_down$p.adjust,
      na.rm = TRUE
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 21. FROZEN HALLMARK GSEA SOURCE
# =============================================================================

hallmark_file <- find_exact_file(
  
  preferred_path = file.path(
    project_dir,
    "results",
    "blood_GSEA_Hallmark_publication_20b_FINAL",
    "20b_FINAL_Hallmark_GSEA_blood_all_pathways.csv"
  ),
  
  basename_target =
    "20b_FINAL_Hallmark_GSEA_blood_all_pathways.csv",
  
  required = TRUE
)


cat("\nFrozen Hallmark GSEA source:\n")

cat(
  normalizePath(
    hallmark_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


hallmark_table <- read.csv(
  hallmark_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


if (
  nrow(
    hallmark_table
  ) !=
  50
) {
  
  stop(
    paste0(
      "Expected 50 Hallmark pathways; observed ",
      nrow(
        hallmark_table
      ),
      "."
    )
  )
}


hallmark_nes_col <- find_column(
  
  hallmark_table,
  
  exact_candidates = c(
    "NES",
    "nes"
  ),
  
  regex =
    "^nes$",
  
  label =
    "Hallmark NES"
)


hallmark_padj_col <- find_column(
  
  hallmark_table,
  
  exact_candidates = c(
    "p.adjust",
    "padj",
    "FDR",
    "qvalue"
  ),
  
  regex =
    "^p\\.adjust$|^padj$|^FDR$|^qvalue$",
  
  label =
    "Hallmark adjusted P value"
)


hallmark_nes <- suppressWarnings(
  as.numeric(
    hallmark_table[[hallmark_nes_col]]
  )
)


hallmark_padj <- suppressWarnings(
  as.numeric(
    hallmark_table[[hallmark_padj_col]]
  )
)


hallmark_table <- data.frame(
  
  Enrichment_direction =
    ifelse(
      hallmark_nes >
        0,
      "Higher in sepsis",
      ifelse(
        hallmark_nes <
          0,
        "Higher in healthy controls",
        "Neutral"
      )
    ),
  
  hallmark_table,
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
)


hallmark_table <- sanitize_for_excel(
  hallmark_table
)


hallmark_audit <- data.frame(
  
  metric = c(
    "Hallmark pathways tested",
    "Positive NES pathways",
    "Negative NES pathways",
    "Positive NES FDR <0.25",
    "Negative NES FDR <0.25",
    "Positive NES FDR <0.05",
    "Negative NES FDR <0.05"
  ),
  
  value = c(
    length(
      hallmark_nes
    ),
    sum(
      hallmark_nes >
        0,
      na.rm = TRUE
    ),
    sum(
      hallmark_nes <
        0,
      na.rm = TRUE
    ),
    sum(
      hallmark_nes >
        0 &
        hallmark_padj <
        0.25,
      na.rm = TRUE
    ),
    sum(
      hallmark_nes <
        0 &
        hallmark_padj <
        0.25,
      na.rm = TRUE
    ),
    sum(
      hallmark_nes >
        0 &
        hallmark_padj <
        0.05,
      na.rm = TRUE
    ),
    sum(
      hallmark_nes <
        0 &
        hallmark_padj <
        0.05,
      na.rm = TRUE
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 22. TABLE S2 README
# =============================================================================

tableS2_readme <- data.frame(
  
  Item = c(
    "Title",
    "Scope",
    "Analysis universe",
    "Primary model",
    "Sensitivity model",
    "DEG threshold",
    "Robust core",
    "Sex-linked genes",
    "Multiple testing",
    "Important note"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S2. Complete differential-expression results ",
      "for sepsis versus healthy-control blood."
    ),
    
    paste0(
      "Complete gene-level differential-expression results for all genes ",
      "retained after the prespecified expression prefilter in the targeted ",
      "whole-blood transcriptomic discovery cohort."
    ),
    
    paste0(
      "The original targeted expression matrix contained approximately ",
      "20,800 targets. Genes with >=10 raw counts in >=3 samples were ",
      "retained. The frozen model-comparison DESeq2 universe contains ",
      n_model_genes,
      " genes."
    ),
    
    "Primary DESeq2 design: ~ condition.",
    
    "Batch/chip-adjusted sensitivity design: ~ batch + condition.",
    
    paste0(
      "Differentially expressed genes were defined by adjusted P <0.05 ",
      "and absolute log2 fold change >=1."
    ),
    
    paste0(
      "Robust-core genes were significant in both models with concordant ",
      "effect direction."
    ),
    
    paste0(
      "Sex-linked QC genes are retained in Table S2 and marked by the ",
      "sex_linked_QC_gene variable. They were not removed before DESeq2."
    ),
    
    "Benjamini-Hochberg false-discovery-rate correction.",
    
    paste0(
      "Script 152 imports frozen DESeq2 outputs and does not re-run ",
      "differential-expression analysis."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 23. TABLE S3 README
# =============================================================================

tableS3_readme <- data.frame(
  
  Item = c(
    "Title",
    "Robust_core_full",
    "Robust_core_no_sex",
    "GO_BP_ORA",
    "KEGG_ORA",
    "WikiPathways_ORA",
    "Hallmark_GSEA",
    "Sex-linked QC",
    "Hallmark ranking",
    "Interpretation limitation",
    "Important note"
  ),
  
  Description = c(
    
    paste0(
      "Supplementary Table S3. Robust-core blood transcriptional response ",
      "and functional enrichment."
    ),
    
    paste0(
      "Complete robust-core DEG set defined by concordant significance ",
      "and effect direction in the primary and batch-adjusted DESeq2 models."
    ),
    
    paste0(
      "Robust-core branch after removal of the predefined sex-linked QC ",
      "gene set used for selected ORA and visualizations."
    ),
    
    "Complete GO Biological Process over-representation results.",
    
    "Complete KEGG over-representation results.",
    
    "Complete WikiPathways over-representation results.",
    
    paste0(
      "Complete Hallmark gene-set enrichment analysis imported from the ",
      "frozen 20b_FINAL workflow."
    ),
    
    paste0(
      "Sex-linked genes were not removed before DESeq2 and were not removed ",
      "from the ranked-list Hallmark GSEA."
    ),
    
    paste0(
      "The frozen Hallmark GSEA workflow used the DESeq2 Wald statistic ",
      "as the preferred gene-ranking metric."
    ),
    
    paste0(
      "ORA enrichment of directional DEG subsets indicates statistical ",
      "over-representation and should not by itself be interpreted as direct ",
      "functional activation or inhibition of a pathway."
    ),
    
    paste0(
      "Script 152 performs table assembly and provenance auditing only; ",
      "ORA and GSEA are not re-run."
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 24. SOURCE / PROVENANCE TABLE
# =============================================================================

source_files <- data.frame(
  
  component = c(
    "Complete DE model comparison",
    "GO BP UP ORA",
    "GO BP DOWN ORA",
    "KEGG UP ORA",
    "KEGG DOWN ORA",
    "WikiPathways UP ORA",
    "WikiPathways DOWN ORA",
    "Hallmark GSEA"
  ),
  
  source = c(
    
    normalizePath(
      model_comparison_file,
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      ora_files[["GO_BP_UP"]],
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      ora_files[["GO_BP_DOWN"]],
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      ora_files[["KEGG_UP"]],
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      ora_files[["KEGG_DOWN"]],
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      ora_files[["WikiPathways_UP"]],
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      ora_files[["WikiPathways_DOWN"]],
      winslash = "\\",
      mustWork = TRUE
    ),
    
    normalizePath(
      hallmark_file,
      winslash = "\\",
      mustWork = TRUE
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 25. FULL AUDIT SUMMARY
# =============================================================================

audit_summary <- data.frame(
  
  metric = c(
    "Genes in frozen DE analysis universe",
    "Primary DEG total",
    "Primary UP",
    "Primary DOWN",
    "Batch-adjusted DEG total",
    "Batch-adjusted UP",
    "Batch-adjusted DOWN",
    "Robust-core total",
    "Robust-core UP",
    "Robust-core DOWN",
    "Sex-linked QC genes in robust core",
    "Robust-core no-sex total",
    "Robust-core no-sex UP",
    "Robust-core no-sex DOWN",
    "Pearson primary vs batch log2FC",
    "Spearman primary vs batch log2FC",
    "Primary stored-vs-derived flag mismatches",
    "Batch stored-vs-derived flag mismatches",
    "Core stored-vs-reconstructed mismatches",
    "GO BP total ORA rows",
    "KEGG total ORA rows",
    "WikiPathways total ORA rows",
    "Hallmark pathways"
  ),
  
  value = c(
    n_model_genes,
    sum(primary_deg),
    sum(primary_up),
    sum(primary_down),
    sum(batch_deg),
    sum(batch_up),
    sum(batch_down),
    sum(robust_core),
    sum(core_up),
    sum(core_down),
    sum(
      sex_linked_audit$present_in_robust_core
    ),
    nrow(
      tableS3_core_no_sex
    ),
    sum(
      tableS3_core_no_sex$robust_core_direction ==
        "UP in sepsis"
    ),
    sum(
      tableS3_core_no_sex$robust_core_direction ==
        "DOWN in sepsis"
    ),
    pearson_lfc,
    spearman_lfc,
    ifelse(
      is.na(
        primary_flag_mismatch
      ),
      "Source flag absent; reconstructed",
      primary_flag_mismatch
    ),
    ifelse(
      is.na(
        batch_flag_mismatch
      ),
      "Source flag absent; reconstructed",
      batch_flag_mismatch
    ),
    ifelse(
      is.na(
        core_flag_vs_reconstructed_mismatch
      ),
      "Source flag absent; reconstructed",
      core_flag_vs_reconstructed_mismatch
    ),
    nrow(
      go_table
    ),
    nrow(
      kegg_table
    ),
    nrow(
      wiki_table
    ),
    nrow(
      hallmark_table
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 26. OUTPUT PATHS
# =============================================================================

tableS2_file <- file.path(
  tables_dir,
  "152_TableS2_complete_blood_differential_expression.xlsx"
)


tableS3_file <- file.path(
  tables_dir,
  "152_TableS3_robust_core_and_functional_enrichment.xlsx"
)


audit_file <- file.path(
  audit_dir,
  "152_INTERNAL_AUDIT_Tables_S2_S3.xlsx"
)


audit_text_file <- file.path(
  text_dir,
  "152_Tables_S2_S3_audit_summary.txt"
)


# =============================================================================
# 27. STYLES
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
# 28. WRITE TABLE S2
# =============================================================================

wbS2 <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wbS2,
  "S2_ReadMe"
)


openxlsx::addWorksheet(
  wbS2,
  "Complete_DE"
)


openxlsx::writeData(
  wbS2,
  "S2_ReadMe",
  tableS2_readme
)


openxlsx::writeData(
  wbS2,
  "Complete_DE",
  tableS2,
  withFilter = TRUE
)


openxlsx::addStyle(
  wbS2,
  "S2_ReadMe",
  readme_header_style,
  rows = 1,
  cols = 1:ncol(
    tableS2_readme
  ),
  gridExpand = TRUE
)


openxlsx::addStyle(
  wbS2,
  "Complete_DE",
  header_style,
  rows = 1,
  cols = 1:ncol(
    tableS2
  ),
  gridExpand = TRUE
)


openxlsx::freezePane(
  wbS2,
  "Complete_DE",
  firstActiveRow = 2,
  firstActiveCol = 2
)


openxlsx::setColWidths(
  wbS2,
  "S2_ReadMe",
  cols = 1,
  widths = 28
)


openxlsx::setColWidths(
  wbS2,
  "S2_ReadMe",
  cols = 2,
  widths = 85
)


for (
  col_idx in seq_len(
    ncol(
      tableS2
    )
  )
) {
  
  desired_width <- max(
    11,
    nchar(
      names(
        tableS2
      )[col_idx]
    ) +
      2
  )
  
  
  desired_width <- min(
    desired_width,
    26
  )
  
  
  openxlsx::setColWidths(
    wbS2,
    "Complete_DE",
    cols = col_idx,
    widths = desired_width
  )
}


openxlsx::saveWorkbook(
  wbS2,
  tableS2_file,
  overwrite = TRUE
)


# =============================================================================
# 29. WRITE TABLE S3
# =============================================================================

wbS3 <- openxlsx::createWorkbook()


s3_sheets <- c(
  "S3_ReadMe",
  "Robust_core_full",
  "Robust_core_no_sex",
  "GO_BP_ORA",
  "KEGG_ORA",
  "WikiPathways_ORA",
  "Hallmark_GSEA"
)


for (
  sheet_name in s3_sheets
) {
  
  openxlsx::addWorksheet(
    wbS3,
    sheet_name
  )
}


openxlsx::writeData(
  wbS3,
  "S3_ReadMe",
  tableS3_readme
)


openxlsx::writeData(
  wbS3,
  "Robust_core_full",
  tableS3_core_full,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS3,
  "Robust_core_no_sex",
  tableS3_core_no_sex,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS3,
  "GO_BP_ORA",
  go_table,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS3,
  "KEGG_ORA",
  kegg_table,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS3,
  "WikiPathways_ORA",
  wiki_table,
  withFilter = TRUE
)


openxlsx::writeData(
  wbS3,
  "Hallmark_GSEA",
  hallmark_table,
  withFilter = TRUE
)


openxlsx::addStyle(
  wbS3,
  "S3_ReadMe",
  readme_header_style,
  rows = 1,
  cols = 1:ncol(
    tableS3_readme
  ),
  gridExpand = TRUE
)


data_sheet_objects <- list(
  Robust_core_full =
    tableS3_core_full,
  Robust_core_no_sex =
    tableS3_core_no_sex,
  GO_BP_ORA =
    go_table,
  KEGG_ORA =
    kegg_table,
  WikiPathways_ORA =
    wiki_table,
  Hallmark_GSEA =
    hallmark_table
)


for (
  sheet_name in names(
    data_sheet_objects
  )
) {
  
  data_object <- data_sheet_objects[[sheet_name]]
  
  
  openxlsx::addStyle(
    wbS3,
    sheet_name,
    header_style,
    rows = 1,
    cols = 1:ncol(
      data_object
    ),
    gridExpand = TRUE
  )
  
  
  openxlsx::freezePane(
    wbS3,
    sheet_name,
    firstActiveRow = 2,
    firstActiveCol = 2
  )
  
  
  openxlsx::setColWidths(
    wbS3,
    sheet_name,
    cols = 1:ncol(
      data_object
    ),
    widths = "auto"
  )
}


openxlsx::setColWidths(
  wbS3,
  "S3_ReadMe",
  cols = 1,
  widths = 30
)


openxlsx::setColWidths(
  wbS3,
  "S3_ReadMe",
  cols = 2,
  widths = 90
)


openxlsx::saveWorkbook(
  wbS3,
  tableS3_file,
  overwrite = TRUE
)


# =============================================================================
# 30. INTERNAL AUDIT WORKBOOK
# =============================================================================

wbAudit <- openxlsx::createWorkbook()


audit_objects <- list(
  
  Audit_summary =
    audit_summary,
  
  DE_summary =
    de_summary,
  
  Column_mapping =
    column_mapping,
  
  Effect_concordance =
    concordance_audit,
  
  Sex_linked_QC =
    sex_linked_audit,
  
  Sex_filter_sources =
    sex_filter_source_audit,
  
  ORA_manifest =
    ora_manifest,
  
  ORA_summary =
    ora_summary,
  
  Hallmark_audit =
    hallmark_audit,
  
  Source_files =
    source_files
)


for (
  sheet_name in names(
    audit_objects
  )
) {
  
  openxlsx::addWorksheet(
    wbAudit,
    sheet_name
  )
  
  
  audit_data <- audit_objects[[sheet_name]]
  
  
  openxlsx::writeData(
    wbAudit,
    sheet_name,
    audit_data
  )
  
  
  openxlsx::addStyle(
    wbAudit,
    sheet_name,
    header_style,
    rows = 1,
    cols = 1:ncol(
      audit_data
    ),
    gridExpand = TRUE
  )
  
  
  openxlsx::freezePane(
    wbAudit,
    sheet_name,
    firstActiveRow = 2
  )
  
  
  openxlsx::setColWidths(
    wbAudit,
    sheet_name,
    cols = 1:ncol(
      audit_data
    ),
    widths = "auto"
  )
}


openxlsx::saveWorkbook(
  wbAudit,
  audit_file,
  overwrite = TRUE
)


# =============================================================================
# 31. TEXT AUDIT SUMMARY
# =============================================================================

audit_text <- c(
  
  "Script 152 FINAL v2 - Tables S2 and S3 audit summary",
  
  "",
  
  paste0(
    "Frozen DE analysis universe: ",
    n_model_genes,
    " genes."
  ),
  
  paste0(
    "Primary DEG: ",
    sum(primary_deg),
    " total; ",
    sum(primary_up),
    " UP; ",
    sum(primary_down),
    " DOWN."
  ),
  
  paste0(
    "Batch-adjusted DEG: ",
    sum(batch_deg),
    " total; ",
    sum(batch_up),
    " UP; ",
    sum(batch_down),
    " DOWN."
  ),
  
  paste0(
    "Robust core: ",
    sum(robust_core),
    " total; ",
    sum(core_up),
    " UP; ",
    sum(core_down),
    " DOWN."
  ),
  
  paste0(
    "Robust core after sex-linked QC removal: ",
    nrow(
      tableS3_core_no_sex
    ),
    "."
  ),
  
  paste0(
    "Pearson primary-vs-batch log2FC: ",
    signif(
      pearson_lfc,
      7
    ),
    "."
  ),
  
  paste0(
    "Spearman primary-vs-batch log2FC: ",
    signif(
      spearman_lfc,
      7
    ),
    "."
  ),
  
  "",
  
  paste0(
    "GO BP ORA rows: ",
    nrow(
      go_table
    ),
    "."
  ),
  
  paste0(
    "KEGG ORA rows: ",
    nrow(
      kegg_table
    ),
    "."
  ),
  
  paste0(
    "WikiPathways ORA rows: ",
    nrow(
      wiki_table
    ),
    "."
  ),
  
  paste0(
    "Hallmark pathways: ",
    nrow(
      hallmark_table
    ),
    "."
  ),
  
  "",
  
  paste0(
    "Table S2: ",
    tableS2_file
  ),
  
  paste0(
    "Table S3: ",
    tableS3_file
  ),
  
  paste0(
    "Internal audit: ",
    audit_file
  )
)


writeLines(
  audit_text,
  audit_text_file
)


# =============================================================================
# 32. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "152_sessionInfo.txt"
  )
)


# =============================================================================
# 33. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 152 FINAL v2 completed successfully.\n")
cat("====================================================================\n\n")


cat("ANALYSIS UNIVERSE\n")
cat("-----------------\n")


cat(
  "Original targeted matrix = approximately 20,800 targets\n"
)


cat(
  "Expression prefilter = >=10 raw counts in >=3 samples\n"
)


cat(
  "Frozen DE universe = ",
  n_model_genes,
  " genes\n",
  sep = ""
)


cat("\nDIFFERENTIAL EXPRESSION\n")
cat("-----------------------\n")


cat(
  "Primary/simple DEG = ",
  sum(primary_deg),
  " (UP ",
  sum(primary_up),
  "; DOWN ",
  sum(primary_down),
  ")\n",
  sep = ""
)


cat(
  "Batch-adjusted DEG = ",
  sum(batch_deg),
  " (UP ",
  sum(batch_up),
  "; DOWN ",
  sum(batch_down),
  ")\n",
  sep = ""
)


cat(
  "Robust core = ",
  sum(robust_core),
  " (UP ",
  sum(core_up),
  "; DOWN ",
  sum(core_down),
  ")\n",
  sep = ""
)


cat(
  "Robust core after sex-linked QC removal = ",
  nrow(
    tableS3_core_no_sex
  ),
  "\n",
  sep = ""
)


cat("\nEFFECT-SIZE CONCORDANCE\n")
cat("-----------------------\n")


cat(
  "Pearson r = ",
  signif(
    pearson_lfc,
    7
  ),
  "\n",
  sep = ""
)


cat(
  "Spearman rho = ",
  signif(
    spearman_lfc,
    7
  ),
  "\n",
  sep = ""
)


cat("\nSTORED-FLAG AUDIT\n")
cat("-----------------\n")


cat(
  "Primary stored-vs-derived mismatches = ",
  ifelse(
    is.na(
      primary_flag_mismatch
    ),
    "source flag absent",
    primary_flag_mismatch
  ),
  "\n",
  sep = ""
)


cat(
  "Batch stored-vs-derived mismatches = ",
  ifelse(
    is.na(
      batch_flag_mismatch
    ),
    "source flag absent",
    batch_flag_mismatch
  ),
  "\n",
  sep = ""
)


cat(
  "Core stored-vs-reconstructed mismatches = ",
  ifelse(
    is.na(
      core_flag_vs_reconstructed_mismatch
    ),
    "source flag absent",
    core_flag_vs_reconstructed_mismatch
  ),
  "\n",
  sep = ""
)


cat("\nSEX-LINKED QC\n")
cat("-------------\n")


cat(
  "Predefined sex-linked QC genes = ",
  length(
    sex_linked_genes
  ),
  "\n",
  sep = ""
)


cat(
  "Sex-linked QC genes present in robust core = ",
  sum(
    sex_linked_audit$present_in_robust_core
  ),
  "\n",
  sep = ""
)


cat(
  "Robust_core_full = ",
  nrow(
    tableS3_core_full
  ),
  "\n",
  sep = ""
)


cat(
  "Robust_core_no_sex = ",
  nrow(
    tableS3_core_no_sex
  ),
  "\n",
  sep = ""
)


cat("\nFUNCTIONAL ENRICHMENT\n")
cat("---------------------\n")


print(
  ora_summary,
  row.names = FALSE
)


cat(
  "\nGO BP total rows = ",
  nrow(
    go_table
  ),
  "\n",
  sep = ""
)


cat(
  "KEGG total rows = ",
  nrow(
    kegg_table
  ),
  "\n",
  sep = ""
)


cat(
  "WikiPathways total rows = ",
  nrow(
    wiki_table
  ),
  "\n",
  sep = ""
)


cat(
  "Hallmark pathways = ",
  nrow(
    hallmark_table
  ),
  "\n",
  sep = ""
)


cat("\nFROZEN ORA SOURCES\n")
cat("------------------\n")


print(
  ora_manifest,
  row.names = FALSE
)


cat("\nHALLMARK AUDIT\n")
cat("--------------\n")


print(
  hallmark_audit,
  row.names = FALSE
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Table S2:\n  ",
  normalizePath(
    tableS2_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Table S3:\n  ",
  normalizePath(
    tableS3_file,
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
  "\n",
  sep = ""
)


cat("\nMANUSCRIPT GUARDRAILS\n")
cat("---------------------\n")


cat(
  "- Table S2 contains all 12,400 genes retained for DE testing.\n"
)


cat(
  "- Low-expression targets removed by the prespecified prefilter are not restored.\n"
)


cat(
  "- Sex-linked genes were not removed before DESeq2.\n"
)


cat(
  "- Sex-linked genes remain in Complete_DE and Robust_core_full.\n"
)


cat(
  "- Robust_core_no_sex represents the selected ORA/QC branch only.\n"
)


cat(
  "- GO/KEGG/WikiPathways are imported from frozen blood_enrichment_core files.\n"
)


cat(
  "- Publication-selected enrichment files are not used as complete ORA sources.\n"
)


cat(
  "- Hallmark GSEA remains a threshold-free ranked-list analysis.\n"
)


cat(
  "- Script 152 does not re-run DESeq2, ORA, or GSEA.\n"
)


cat("\nDone.\n")