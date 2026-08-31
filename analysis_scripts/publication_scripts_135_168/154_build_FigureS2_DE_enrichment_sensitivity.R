################################################################################
# Script 154
# FINAL v2
#
# Supplementary Figure S2
#
# Differential-expression robustness and functional enrichment
#
# Project:
#   Sepsis_DESeq2
#
# Manuscript:
#   Blood-only sepsis transcriptomic endotypes /
#   five-gene host-response signature
#
#
# FIGURE S2 PANELS
# ----------------
#
# A. Number of UP and DOWN DEGs in:
#      - primary/simple model
#      - batch-adjusted sensitivity model
#
# B. DEG overlap:
#      - primary only
#      - shared concordant / robust core
#      - batch-adjusted only
#
# C. Gene-level log2FC concordance between primary and batch-adjusted models
#
# D. Top significant GO Biological Process ORA terms
#
# E. Top significant KEGG ORA terms
#
# F. Top significant WikiPathways ORA terms
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
# It imports:
#
#   FINAL submission Table S2:
#     12,393 biological targets
#
#   Frozen Table S3:
#     robust core + GO/KEGG/WikiPathways enrichment
#
#
# EXPECTED FROZEN RESULTS
# -----------------------
#
# Primary/simple:
#   2,659 DEG
#   1,660 UP
#     999 DOWN
#
# Batch-adjusted:
#   4,125 DEG
#   2,093 UP
#   2,032 DOWN
#
# Shared concordant robust core:
#   1,796
#   1,133 UP
#     663 DOWN
#
# Expected model overlap:
#
#   Primary only       =   863
#   Shared robust core = 1,796
#   Batch only         = 2,329
#   Union              = 4,988
#
#
# EFFECT-SIZE CONCORDANCE
# -----------------------
#
# IMPORTANT provenance distinction:
#
# Historical 12,400-feature model-comparison universe:
#
#   Pearson r    = 0.8149663
#   Spearman rho = 0.8594785
#
# Final biological-only submission universe (12,393 targets):
#
#   Pearson r    ~ 0.815176
#   Spearman rho ~ 0.859938
#
# Figure S2 uses the FINAL biological-only universe.
#
#
# ORA TABLES
# ----------
#
# Complete significant ORA results:
#
#   GO BP:
#     411 UP
#     129 DOWN
#
#   KEGG:
#      41 UP
#      22 DOWN
#
#   WikiPathways:
#      35 UP
#      18 DOWN
#
#
# HALLMARK
# --------
#
# Hallmark GSEA is intentionally NOT repeated in Figure S2 because:
#
#   - Hallmark is already shown in Main Figure 1E;
#   - all 50 Hallmark results are reported in Supplementary Table S3.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 154 FINAL v2\n")
cat("Supplementary Figure S2\n")
cat("Differential-expression robustness and functional enrichment\n")
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
  "ggplot2",
  "readxl",
  "openxlsx",
  "gridExtra"
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
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(gridExtra)
  
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "154_FigureS2_DE_enrichment_sensitivity"
)


figures_dir <- file.path(
  output_dir,
  "figures"
)


tables_dir <- file.path(
  output_dir,
  "tables"
)


text_dir <- file.path(
  output_dir,
  "text"
)


audit_dir <- file.path(
  output_dir,
  "audit"
)


for (
  directory_name in c(
    output_dir,
    figures_dir,
    tables_dir,
    text_dir,
    audit_dir
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


normalize_gene <- function(x) {
  
  toupper(
    stringr::str_trim(
      as.character(x)
    )
  )
}


check_columns <- function(
    data,
    required_columns,
    object_name
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
      "\nAvailable columns in ",
      object_name,
      ":\n",
      sep = ""
    )
    
    
    print(
      names(
        data
      )
    )
    
    
    stop(
      paste0(
        "Missing required column(s) in ",
        object_name,
        ": ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
}


prepare_top_ora <- function(
    data,
    database_name,
    n_per_direction = 5
) {
  
  check_columns(
    data,
    c(
      "Direction",
      "Description",
      "p.adjust",
      "Count"
    ),
    database_name
  )
  
  
  temp <- data %>%
    
    dplyr::mutate(
      
      p.adjust =
        suppressWarnings(
          as.numeric(
            p.adjust
          )
        ),
      
      Count =
        suppressWarnings(
          as.numeric(
            Count
          )
        )
    ) %>%
    
    dplyr::filter(
      is.finite(
        p.adjust
      ),
      p.adjust <
        0.05,
      is.finite(
        Count
      )
    ) %>%
    
    dplyr::group_by(
      Direction
    ) %>%
    
    dplyr::slice_min(
      order_by = p.adjust,
      n = n_per_direction,
      with_ties = FALSE
    ) %>%
    
    dplyr::ungroup() %>%
    
    dplyr::mutate(
      
      neg_log10_FDR =
        -log10(
          pmax(
            p.adjust,
            .Machine$double.xmin
          )
        ),
      
      Description_wrapped =
        stringr::str_wrap(
          as.character(
            Description
          ),
          width = 40
        ),
      
      term_key =
        paste0(
          Description_wrapped,
          "|||",
          Direction
        )
    )
  
  
  level_order <- temp %>%
    
    dplyr::arrange(
      Direction,
      p.adjust
    ) %>%
    
    dplyr::pull(
      term_key
    )
  
  
  temp$term_key <- factor(
    temp$term_key,
    levels = rev(
      unique(
        level_order
      )
    )
  )
  
  
  temp
}


clean_term_labels <- function(x) {
  
  sub(
    "\\|\\|\\|.*$",
    "",
    x
  )
}


# =============================================================================
# 5. INPUT FILES
# =============================================================================

tableS2_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "152b_final_submission_TableS2",
  "tables",
  "152b_TableS2_FINAL_SUBMISSION_biological_targets.xlsx"
)


tableS3_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "152_Tables_S2_S3_DE_enrichment",
  "tables",
  "152_TableS3_robust_core_and_functional_enrichment.xlsx"
)


if (
  !file.exists(
    tableS2_file
  )
) {
  
  stop(
    paste0(
      "FINAL Table S2 not found:\n",
      tableS2_file
    )
  )
}


if (
  !file.exists(
    tableS3_file
  )
) {
  
  stop(
    paste0(
      "Frozen Table S3 not found:\n",
      tableS3_file
    )
  )
}


cat("\nINPUT FILES\n")
cat("-----------\n")


cat(
  "FINAL Table S2:\n  ",
  normalizePath(
    tableS2_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Frozen Table S3:\n  ",
  normalizePath(
    tableS3_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 6. READ FINAL TABLE S2
# =============================================================================

s2_sheets <- readxl::excel_sheets(
  tableS2_file
)


if (
  !("Complete_DE" %in%
    s2_sheets)
) {
  
  stop(
    "Sheet Complete_DE missing from FINAL Table S2."
  )
}


de <- readxl::read_excel(
  tableS2_file,
  sheet = "Complete_DE"
) %>%
  
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "\nFINAL Table S2 dimensions = ",
  nrow(
    de
  ),
  " x ",
  ncol(
    de
  ),
  "\n",
  sep = ""
)


required_de_columns <- c(
  "Gene",
  "primary_log2FoldChange",
  "primary_padj",
  "primary_DEG",
  "batch_adjusted_log2FoldChange",
  "batch_adjusted_padj",
  "batch_adjusted_DEG",
  "robust_core_DEG",
  "robust_core_direction",
  "sex_linked_QC_gene"
)


check_columns(
  de,
  required_de_columns,
  "FINAL Table S2 Complete_DE"
)


if (
  nrow(
    de
  ) !=
  12393
) {
  
  stop(
    paste0(
      "Expected 12,393 biological targets in FINAL Table S2; observed ",
      nrow(
        de
      )
    )
  )
}


if (
  any(
    grepl(
      "^ERCC[-_]",
      de$Gene,
      ignore.case = TRUE
    )
  )
) {
  
  stop(
    "ERCC-prefixed features unexpectedly present in FINAL submission Table S2."
  )
}


# =============================================================================
# 7. STANDARDIZE DE VARIABLES
# =============================================================================

de <- de %>%
  
  dplyr::mutate(
    
    primary_log2FoldChange =
      suppressWarnings(
        as.numeric(
          primary_log2FoldChange
        )
      ),
    
    batch_adjusted_log2FoldChange =
      suppressWarnings(
        as.numeric(
          batch_adjusted_log2FoldChange
        )
      ),
    
    primary_padj =
      suppressWarnings(
        as.numeric(
          primary_padj
        )
      ),
    
    batch_adjusted_padj =
      suppressWarnings(
        as.numeric(
          batch_adjusted_padj
        )
      ),
    
    primary_DEG =
      as_logical_flag(
        primary_DEG
      ),
    
    batch_adjusted_DEG =
      as_logical_flag(
        batch_adjusted_DEG
      ),
    
    robust_core_DEG =
      as_logical_flag(
        robust_core_DEG
      ),
    
    sex_linked_QC_gene =
      as_logical_flag(
        sex_linked_QC_gene
      )
  )


# =============================================================================
# 8. HARD DE AUDIT
# =============================================================================

primary_up <- sum(
  de$primary_DEG &
    de$primary_log2FoldChange >
    0
)


primary_down <- sum(
  de$primary_DEG &
    de$primary_log2FoldChange <
    0
)


batch_up <- sum(
  de$batch_adjusted_DEG &
    de$batch_adjusted_log2FoldChange >
    0
)


batch_down <- sum(
  de$batch_adjusted_DEG &
    de$batch_adjusted_log2FoldChange <
    0
)


core_up <- sum(
  de$robust_core_DEG &
    de$primary_log2FoldChange >
    0
)


core_down <- sum(
  de$robust_core_DEG &
    de$primary_log2FoldChange <
    0
)


de_audit <- data.frame(
  
  analysis = c(
    "Primary/simple",
    "Batch-adjusted",
    "Robust core"
  ),
  
  total = c(
    sum(
      de$primary_DEG
    ),
    sum(
      de$batch_adjusted_DEG
    ),
    sum(
      de$robust_core_DEG
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


de_audit$match <-
  de_audit$total ==
  de_audit$expected_total &
  de_audit$UP ==
  de_audit$expected_UP &
  de_audit$DOWN ==
  de_audit$expected_DOWN


cat("\nDE AUDIT\n")
cat("--------\n")


print(
  de_audit,
  row.names = FALSE
)


if (
  !all(
    de_audit$match
  )
) {
  
  stop(
    "Frozen DE anchor audit failed."
  )
}


# =============================================================================
# 9. MODEL OVERLAP AUDIT
# =============================================================================

primary_only <- sum(
  de$primary_DEG &
    !de$batch_adjusted_DEG
)


batch_only <- sum(
  !de$primary_DEG &
    de$batch_adjusted_DEG
)


shared_sig <-
  de$primary_DEG &
  de$batch_adjusted_DEG


shared_total <- sum(
  shared_sig
)


shared_concordant <- sum(
  shared_sig &
    sign(
      de$primary_log2FoldChange
    ) ==
    sign(
      de$batch_adjusted_log2FoldChange
    )
)


shared_discordant <- sum(
  shared_sig &
    sign(
      de$primary_log2FoldChange
    ) !=
    sign(
      de$batch_adjusted_log2FoldChange
    )
)


union_deg <- sum(
  de$primary_DEG |
    de$batch_adjusted_DEG
)


overlap_audit <- data.frame(
  
  category = c(
    "Primary only",
    "Significant in both models",
    "Shared concordant / robust core",
    "Shared discordant",
    "Batch-adjusted only",
    "Union of either model"
  ),
  
  n = c(
    primary_only,
    shared_total,
    shared_concordant,
    shared_discordant,
    batch_only,
    union_deg
  ),
  
  stringsAsFactors = FALSE
)


cat("\nMODEL OVERLAP AUDIT\n")
cat("-------------------\n")


print(
  overlap_audit,
  row.names = FALSE
)


if (
  primary_only !=
  863
) {
  
  stop(
    paste0(
      "Expected primary-only n=863; observed ",
      primary_only
    )
  )
}


if (
  batch_only !=
  2329
) {
  
  stop(
    paste0(
      "Expected batch-only n=2329; observed ",
      batch_only
    )
  )
}


if (
  shared_total !=
  1796
) {
  
  stop(
    paste0(
      "Expected shared significant n=1796; observed ",
      shared_total
    )
  )
}


if (
  shared_concordant !=
  1796 ||
  shared_discordant !=
  0
) {
  
  stop(
    "Shared DEG concordance audit failed."
  )
}


if (
  union_deg !=
  4988
) {
  
  stop(
    paste0(
      "Expected DEG union n=4988; observed ",
      union_deg
    )
  )
}


# =============================================================================
# 10. EFFECT-SIZE CONCORDANCE — BIOLOGICAL TARGETS ONLY
# =============================================================================

valid_lfc <-
  is.finite(
    de$primary_log2FoldChange
  ) &
  is.finite(
    de$batch_adjusted_log2FoldChange
  )


pearson_r <- stats::cor(
  de$primary_log2FoldChange[valid_lfc],
  de$batch_adjusted_log2FoldChange[valid_lfc],
  method = "pearson"
)


spearman_rho <- stats::cor(
  de$primary_log2FoldChange[valid_lfc],
  de$batch_adjusted_log2FoldChange[valid_lfc],
  method = "spearman"
)


concordance_audit <- data.frame(
  
  metric = c(
    "Pearson r",
    "Spearman rho"
  ),
  
  value = c(
    pearson_r,
    spearman_rho
  ),
  
  expected_approx = c(
    0.81518,
    0.85994
  ),
  
  difference = c(
    pearson_r - 0.81518,
    spearman_rho - 0.85994
  ),
  
  stringsAsFactors = FALSE
)


cat("\nEFFECT-SIZE CONCORDANCE — BIOLOGICAL TARGETS ONLY\n")
cat("-------------------------------------------------\n")


print(
  concordance_audit,
  row.names = FALSE
)


# Use a scientifically meaningful tolerance.
# The purpose is source-provenance auditing, not floating-point identity.

if (
  abs(
    pearson_r -
    0.81518
  ) >
  0.002
) {
  
  stop(
    paste0(
      "Unexpected biological-only Pearson r: ",
      pearson_r
    )
  )
}


if (
  abs(
    spearman_rho -
    0.85994
  ) >
  0.002
) {
  
  stop(
    paste0(
      "Unexpected biological-only Spearman rho: ",
      spearman_rho
    )
  )
}


cat(
  "\nBiological-only effect-size concordance audit passed.\n"
)


# =============================================================================
# 11. SEX-LINKED QC AUDIT
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


sex_linked_present <- de %>%
  
  dplyr::filter(
    sex_linked_QC_gene
  )


sex_linked_audit <- data.frame(
  
  metric = c(
    "Predefined sex-linked QC genes",
    "Present after biological expression prefilter",
    "Primary DEGs among present sex-linked genes",
    "Batch-adjusted DEGs among present sex-linked genes",
    "Robust-core DEGs among present sex-linked genes"
  ),
  
  value = c(
    length(
      sex_linked_genes
    ),
    nrow(
      sex_linked_present
    ),
    sum(
      sex_linked_present$primary_DEG
    ),
    sum(
      sex_linked_present$batch_adjusted_DEG
    ),
    sum(
      sex_linked_present$robust_core_DEG
    )
  ),
  
  stringsAsFactors = FALSE
)


cat("\nSEX-LINKED QC AUDIT\n")
cat("-------------------\n")


print(
  sex_linked_audit,
  row.names = FALSE
)


if (
  sum(
    sex_linked_present$robust_core_DEG
  ) !=
  0
) {
  
  stop(
    "A predefined sex-linked QC gene unexpectedly entered the robust core."
  )
}


# =============================================================================
# 12. READ FROZEN TABLE S3
# =============================================================================

s3_sheets <- readxl::excel_sheets(
  tableS3_file
)


required_s3_sheets <- c(
  "Robust_core_full",
  "GO_BP_ORA",
  "KEGG_ORA",
  "WikiPathways_ORA"
)


missing_s3_sheets <- setdiff(
  required_s3_sheets,
  s3_sheets
)


if (
  length(
    missing_s3_sheets
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing Table S3 sheet(s): ",
      paste(
        missing_s3_sheets,
        collapse = ", "
      )
    )
  )
}


robust_core_s3 <- readxl::read_excel(
  tableS3_file,
  sheet = "Robust_core_full"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


go_table <- readxl::read_excel(
  tableS3_file,
  sheet = "GO_BP_ORA"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


kegg_table <- readxl::read_excel(
  tableS3_file,
  sheet = "KEGG_ORA"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


wiki_table <- readxl::read_excel(
  tableS3_file,
  sheet = "WikiPathways_ORA"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat("\nTABLE S3 SOURCE AUDIT\n")
cat("---------------------\n")


cat(
  "Robust_core_full rows = ",
  nrow(
    robust_core_s3
  ),
  "\n",
  sep = ""
)


cat(
  "GO BP ORA rows = ",
  nrow(
    go_table
  ),
  "\n",
  sep = ""
)


cat(
  "KEGG ORA rows = ",
  nrow(
    kegg_table
  ),
  "\n",
  sep = ""
)


cat(
  "WikiPathways ORA rows = ",
  nrow(
    wiki_table
  ),
  "\n",
  sep = ""
)


if (
  nrow(
    robust_core_s3
  ) !=
  1796
) {
  
  stop(
    "Table S3 robust-core row count mismatch."
  )
}


if (
  nrow(
    go_table
  ) !=
  540
) {
  
  stop(
    "Expected 540 GO BP ORA rows."
  )
}


if (
  nrow(
    kegg_table
  ) !=
  63
) {
  
  stop(
    "Expected 63 KEGG ORA rows."
  )
}


if (
  nrow(
    wiki_table
  ) !=
  53
) {
  
  stop(
    "Expected 53 WikiPathways ORA rows."
  )
}


# =============================================================================
# 13. VERIFY ROBUST CORE BETWEEN FINAL S2 AND S3
# =============================================================================

core_genes_s2 <- sort(
  normalize_gene(
    de$Gene[
      de$robust_core_DEG
    ]
  )
)


core_gene_col_s3 <- if (
  "Gene" %in%
  names(
    robust_core_s3
  )
) {
  
  "Gene"
  
} else {
  
  stop(
    "Gene column missing from Table S3 Robust_core_full."
  )
}


core_genes_s3 <- sort(
  normalize_gene(
    robust_core_s3[[core_gene_col_s3]]
  )
)


if (
  !identical(
    core_genes_s2,
    core_genes_s3
  )
) {
  
  stop(
    "Robust-core gene set differs between FINAL Table S2 and Table S3."
  )
}


cat(
  "Robust-core S2/S3 gene-set audit passed: 1,796/1,796 identical.\n"
)


# =============================================================================
# 14. PREPARE TOP ORA TERMS
# =============================================================================

go_top <- prepare_top_ora(
  go_table,
  database_name = "GO Biological Process",
  n_per_direction = 5
)


kegg_top <- prepare_top_ora(
  kegg_table,
  database_name = "KEGG",
  n_per_direction = 5
)


wiki_top <- prepare_top_ora(
  wiki_table,
  database_name = "WikiPathways",
  n_per_direction = 5
)


# =============================================================================
# 15. COLOR PALETTES
# =============================================================================

direction_colors <- c(
  
  "UP in sepsis" =
    "#B2182B",
  
  "DOWN in sepsis" =
    "#2166AC"
)


overlap_colors <- c(
  
  "Primary only" =
    "#636363",
  
  "Shared robust core" =
    "#54278F",
  
  "Batch-adjusted only" =
    "#9E9AC8"
)


# =============================================================================
# 16. PANEL A — DEG COUNTS
# =============================================================================

panel_A_data <- data.frame(
  
  Model = factor(
    c(
      "Primary/simple",
      "Primary/simple",
      "Batch-adjusted",
      "Batch-adjusted"
    ),
    levels = c(
      "Primary/simple",
      "Batch-adjusted"
    )
  ),
  
  Direction = factor(
    c(
      "UP in sepsis",
      "DOWN in sepsis",
      "UP in sepsis",
      "DOWN in sepsis"
    ),
    levels = c(
      "UP in sepsis",
      "DOWN in sepsis"
    )
  ),
  
  n = c(
    primary_up,
    primary_down,
    batch_up,
    batch_down
  )
)


panel_A_totals <- panel_A_data %>%
  
  dplyr::group_by(
    Model
  ) %>%
  
  dplyr::summarise(
    total =
      sum(
        n
      ),
    .groups = "drop"
  )


panel_A <- ggplot2::ggplot(
  
  panel_A_data,
  
  ggplot2::aes(
    x = Model,
    y = n,
    fill = Direction
  )
  
) +
  
  ggplot2::geom_col(
    width = 0.62
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = n
    ),
    
    position =
      ggplot2::position_stack(
        vjust = 0.5
      ),
    
    size =
      4,
    
    color =
      "white",
    
    fontface =
      "bold"
  ) +
  
  ggplot2::geom_text(
    
    data =
      panel_A_totals,
    
    ggplot2::aes(
      x = Model,
      y = total,
      label = paste0(
        "Total ",
        total
      )
    ),
    
    inherit.aes =
      FALSE,
    
    vjust =
      -0.6,
    
    size =
      3.8
  ) +
  
  ggplot2::scale_fill_manual(
    values = direction_colors
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  ggplot2::labs(
    title = "A  Differentially expressed genes",
    x = NULL,
    y = "Number of DEGs",
    fill = NULL
  ) +
  
  ggplot2::theme_bw(
    base_size = 11
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    
    legend.position =
      "bottom"
  )


# =============================================================================
# 17. PANEL B — MODEL OVERLAP
# =============================================================================

panel_B_data <- data.frame(
  
  Category = factor(
    c(
      "Primary only",
      "Shared robust core",
      "Batch-adjusted only"
    ),
    levels = c(
      "Primary only",
      "Shared robust core",
      "Batch-adjusted only"
    )
  ),
  
  n = c(
    primary_only,
    shared_concordant,
    batch_only
  )
)


panel_B <- ggplot2::ggplot(
  
  panel_B_data,
  
  ggplot2::aes(
    x = Category,
    y = n,
    fill = Category
  )
  
) +
  
  ggplot2::geom_col(
    width = 0.64
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = n
    ),
    
    vjust =
      -0.45,
    
    size =
      4
  ) +
  
  ggplot2::scale_fill_manual(
    values = overlap_colors
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  ggplot2::labs(
    title = "B  DEG overlap across models",
    subtitle = "All 1,796 shared DEGs had concordant direction",
    x = NULL,
    y = "Number of genes"
  ) +
  
  ggplot2::theme_bw(
    base_size = 11
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 9.5
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 20,
        hjust = 1
      )
  )


# =============================================================================
# 18. PANEL C — BIOLOGICAL-TARGET EFFECT-SIZE CONCORDANCE
# =============================================================================

panel_C_data <- de %>%
  
  dplyr::mutate(
    
    Plot_class =
      dplyr::case_when(
        
        robust_core_DEG &
          primary_log2FoldChange >
          0 ~
          "Robust core UP",
        
        robust_core_DEG &
          primary_log2FoldChange <
          0 ~
          "Robust core DOWN",
        
        TRUE ~
          "Other tested gene"
      )
  )


panel_C_colors <- c(
  
  "Other tested gene" =
    "#BDBDBD",
  
  "Robust core UP" =
    "#B2182B",
  
  "Robust core DOWN" =
    "#2166AC"
)


panel_C <- ggplot2::ggplot(
  
  panel_C_data,
  
  ggplot2::aes(
    x = primary_log2FoldChange,
    y = batch_adjusted_log2FoldChange
  )
  
) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    color = "#BDBDBD"
  ) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    color = "#BDBDBD"
  ) +
  
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "#4D4D4D"
  ) +
  
  ggplot2::geom_point(
    
    data = panel_C_data %>%
      dplyr::filter(
        Plot_class ==
          "Other tested gene"
      ),
    
    ggplot2::aes(
      color = Plot_class
    ),
    
    alpha = 0.24,
    size = 1
  ) +
  
  ggplot2::geom_point(
    
    data = panel_C_data %>%
      dplyr::filter(
        Plot_class !=
          "Other tested gene"
      ),
    
    ggplot2::aes(
      color = Plot_class
    ),
    
    alpha = 0.72,
    size = 1.5
  ) +
  
  ggplot2::scale_color_manual(
    values = panel_C_colors
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x =
      -Inf,
    
    y =
      Inf,
    
    hjust =
      -0.05,
    
    vjust =
      1.3,
    
    label =
      paste0(
        "Pearson r = ",
        sprintf(
          "%.3f",
          pearson_r
        ),
        "\nSpearman \u03c1 = ",
        sprintf(
          "%.3f",
          spearman_rho
        )
      ),
    
    size =
      3.7
  ) +
  
  ggplot2::labs(
    title = "C  Effect-size concordance",
    x = "Primary model log2 fold change",
    y = "Batch-adjusted model log2 fold change",
    color = NULL
  ) +
  
  ggplot2::theme_bw(
    base_size = 11
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    
    legend.position =
      "bottom"
  )


# =============================================================================
# 19. ORA PANEL FUNCTION
# =============================================================================

make_ora_panel <- function(
    data,
    panel_title
) {
  
  ggplot2::ggplot(
    
    data,
    
    ggplot2::aes(
      x = neg_log10_FDR,
      y = term_key,
      size = Count,
      color = Direction
    )
    
  ) +
    
    ggplot2::geom_point(
      alpha = 0.85
    ) +
    
    ggplot2::facet_wrap(
      ~ Direction,
      scales = "free_y",
      ncol = 1
    ) +
    
    ggplot2::scale_color_manual(
      values = direction_colors
    ) +
    
    ggplot2::scale_y_discrete(
      labels = clean_term_labels
    ) +
    
    ggplot2::labs(
      title = panel_title,
      x = expression(
        -log[10](
          "BH-adjusted P"
        )
      ),
      y = NULL,
      size = "Gene count",
      color = NULL
    ) +
    
    ggplot2::theme_bw(
      base_size = 9
    ) +
    
    ggplot2::theme(
      
      plot.title =
        ggplot2::element_text(
          face = "bold",
          size = 12
        ),
      
      strip.text =
        ggplot2::element_text(
          face = "bold",
          size = 8.5
        ),
      
      legend.position =
        "bottom",
      
      axis.text.y =
        ggplot2::element_text(
          size = 7.5
        ),
      
      panel.grid.major.y =
        ggplot2::element_blank()
    )
}


# =============================================================================
# 20. PANELS D-F
# =============================================================================

panel_D <- make_ora_panel(
  go_top,
  "D  GO Biological Process"
)


panel_E <- make_ora_panel(
  kegg_top,
  "E  KEGG pathways"
)


panel_F <- make_ora_panel(
  wiki_top,
  "F  WikiPathways"
)


# =============================================================================
# 21. COMBINE FIGURE S2
# =============================================================================

figure_S2_grob <- gridExtra::arrangeGrob(
  
  panel_A,
  panel_B,
  
  panel_C,
  panel_D,
  
  panel_E,
  panel_F,
  
  ncol = 2,
  
  widths = c(
    1,
    1
  ),
  
  heights = c(
    0.85,
    1.25,
    1.25
  )
)


# =============================================================================
# 22. SAVE FIGURE
# =============================================================================

figure_pdf <- file.path(
  figures_dir,
  "154_FigureS2_DE_enrichment_sensitivity.pdf"
)


figure_png <- file.path(
  figures_dir,
  "154_FigureS2_DE_enrichment_sensitivity.png"
)


figure_tiff <- file.path(
  figures_dir,
  "154_FigureS2_DE_enrichment_sensitivity.tiff"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure_S2_grob,
  width = 15,
  height = 17,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure_S2_grob,
  width = 15,
  height = 17,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure_S2_grob,
  width = 15,
  height = 17,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 23. SOURCE-DATA WORKBOOK
# =============================================================================

source_data_file <- file.path(
  tables_dir,
  "154_FigureS2_source_data.xlsx"
)


wb <- openxlsx::createWorkbook()


source_objects <- list(
  
  DE_audit =
    de_audit,
  
  DEG_counts_panel_A =
    panel_A_data,
  
  Model_overlap =
    overlap_audit,
  
  Overlap_panel_B =
    panel_B_data,
  
  Effect_concordance =
    concordance_audit,
  
  Effect_scatter =
    panel_C_data %>%
    dplyr::select(
      Gene,
      primary_log2FoldChange,
      batch_adjusted_log2FoldChange,
      primary_DEG,
      batch_adjusted_DEG,
      robust_core_DEG,
      robust_core_direction,
      Plot_class
    ),
  
  Sex_linked_QC =
    sex_linked_audit,
  
  Sex_linked_genes =
    sex_linked_present,
  
  GO_top =
    go_top,
  
  KEGG_top =
    kegg_top,
  
  WikiPathways_top =
    wiki_top
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE
)


for (
  sheet_name in names(
    source_objects
  )
) {
  
  data_object <-
    source_objects[[sheet_name]]
  
  
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
    header_style,
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


openxlsx::saveWorkbook(
  wb,
  source_data_file,
  overwrite = TRUE
)


# =============================================================================
# 24. FIGURE LEGEND
# =============================================================================

figure_legend <- paste0(
  
  "Supplementary Figure S2. Robustness of blood differential-expression ",
  "results and functional interpretation. ",
  
  "(A) Numbers of genes with increased and decreased expression in sepsis ",
  "in the primary DESeq2 model and in the sequencing-batch-adjusted ",
  "sensitivity model. ",
  
  "(B) Partition of differentially expressed genes between the two models. ",
  "Of the genes significant in both models, all 1,796 showed concordant ",
  "effect direction and constituted the robust-core transcriptional response. ",
  
  "(C) Gene-level comparison of log2 fold changes between the primary and ",
  "batch-adjusted models across the 12,393 biological targets reported in ",
  "Supplementary Table S2. Robust-core genes are highlighted according to ",
  "direction of change. Pearson r=",
  sprintf(
    "%.3f",
    pearson_r
  ),
  " and Spearman rho=",
  sprintf(
    "%.3f",
    spearman_rho
  ),
  ". The dashed line denotes identity. ",
  
  "(D-F) Most statistically significant over-representation results for ",
  "GO Biological Process, KEGG, and WikiPathways, respectively, shown ",
  "separately for genes with increased and decreased expression in sepsis. ",
  "Point size indicates the number of genes contributing to each term and ",
  "the x-axis represents -log10 of the Benjamini-Hochberg-adjusted P value. ",
  
  "Complete significant ORA results are provided in Supplementary Table S3. ",
  "Hallmark ranked-list GSEA is shown in Main Fig. 1E and reported in full ",
  "in Supplementary Table S3. Sequencing batch was partially structured by ",
  "study group; accordingly, the batch-adjusted model is interpreted as a ",
  "conservative sensitivity analysis rather than as a replacement for the ",
  "prespecified primary model."
)


legend_file <- file.path(
  text_dir,
  "154_FigureS2_legend_EN.txt"
)


writeLines(
  figure_legend,
  legend_file
)


# =============================================================================
# 25. PROPOSED RESULTS 3.2 TEXT
# =============================================================================

results_text <- paste0(
  
  "The prespecified primary analysis identified 2,659 differentially ",
  "expressed biological genes, including 1,660 with increased and 999 with ",
  "decreased expression in sepsis. Because sequencing batch was partially ",
  "structured by study group (Supplementary Fig. S1D), differential ",
  "expression was additionally examined using a batch-adjusted sensitivity ",
  "model, which identified 4,125 differentially expressed genes (2,093 ",
  "increased and 2,032 decreased). Across the 12,393 retained biological ",
  "targets, gene-level effect estimates remained strongly concordant between ",
  "models (Pearson r=",
  sprintf(
    "%.3f",
    pearson_r
  ),
  "; Spearman rho=",
  sprintf(
    "%.3f",
    spearman_rho
  ),
  "). A total of 1,796 genes were significant in both models, and all showed ",
  "concordant direction of change, defining a robust transcriptional core of ",
  "1,133 genes increased and 663 decreased in sepsis (Supplementary Fig. S2; ",
  "Supplementary Tables S2-S3). Functional over-representation analysis of ",
  "this robust response demonstrated broad enrichment of inflammatory and ",
  "myeloid-associated processes among genes increased in sepsis, together ",
  "with enrichment of adaptive immune-associated processes among genes with ",
  "decreased expression (Supplementary Fig. S2D-F and Supplementary Table S3)."
)


results_file <- file.path(
  text_dir,
  "154_proposed_Results_3.2_DE_enrichment_EN.txt"
)


writeLines(
  results_text,
  results_file
)


# =============================================================================
# 26. AUDIT SUMMARY
# =============================================================================

audit_summary <- data.frame(
  
  metric = c(
    "Biological targets in FINAL Table S2",
    "Primary DEG",
    "Primary UP",
    "Primary DOWN",
    "Batch-adjusted DEG",
    "Batch UP",
    "Batch DOWN",
    "Primary-only DEG",
    "Shared significant DEG",
    "Shared concordant DEG",
    "Shared discordant DEG",
    "Batch-only DEG",
    "Union DEG",
    "Robust core",
    "Robust-core UP",
    "Robust-core DOWN",
    "Pearson biological-only log2FC correlation",
    "Spearman biological-only log2FC correlation",
    "Predefined sex-linked QC genes",
    "Sex-linked genes present after prefilter",
    "Primary sex-linked DEG",
    "Batch-adjusted sex-linked DEG",
    "Robust-core sex-linked DEG",
    "GO BP ORA rows",
    "KEGG ORA rows",
    "WikiPathways ORA rows"
  ),
  
  value = c(
    nrow(
      de
    ),
    sum(
      de$primary_DEG
    ),
    primary_up,
    primary_down,
    sum(
      de$batch_adjusted_DEG
    ),
    batch_up,
    batch_down,
    primary_only,
    shared_total,
    shared_concordant,
    shared_discordant,
    batch_only,
    union_deg,
    sum(
      de$robust_core_DEG
    ),
    core_up,
    core_down,
    pearson_r,
    spearman_rho,
    length(
      sex_linked_genes
    ),
    nrow(
      sex_linked_present
    ),
    sum(
      sex_linked_present$primary_DEG
    ),
    sum(
      sex_linked_present$batch_adjusted_DEG
    ),
    sum(
      sex_linked_present$robust_core_DEG
    ),
    nrow(
      go_table
    ),
    nrow(
      kegg_table
    ),
    nrow(
      wiki_table
    )
  ),
  
  stringsAsFactors = FALSE
)


audit_file <- file.path(
  audit_dir,
  "154_FigureS2_audit.xlsx"
)


wb_audit <- openxlsx::createWorkbook()


audit_objects <- list(
  
  Audit_summary =
    audit_summary,
  
  DE_audit =
    de_audit,
  
  Model_overlap =
    overlap_audit,
  
  Effect_concordance =
    concordance_audit,
  
  Sex_linked_QC =
    sex_linked_audit
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
    data_object
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
# 27. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "154_sessionInfo.txt"
  )
)


# =============================================================================
# 28. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 154 FINAL v2 completed successfully.\n")
cat("====================================================================\n\n")


cat("FINAL TABLE S2\n")
cat("--------------\n")


cat(
  "Biological targets = ",
  nrow(
    de
  ),
  "\n",
  sep = ""
)


cat(
  "ERCC-prefixed rows = 0\n"
)


cat("\nDIFFERENTIAL EXPRESSION\n")
cat("-----------------------\n")


cat(
  "Primary DEG = ",
  sum(
    de$primary_DEG
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
    de$batch_adjusted_DEG
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
    de$robust_core_DEG
  ),
  " (UP ",
  core_up,
  "; DOWN ",
  core_down,
  ")\n",
  sep = ""
)


cat("\nMODEL OVERLAP\n")
cat("-------------\n")


cat(
  "Primary only = ",
  primary_only,
  "\n",
  sep = ""
)


cat(
  "Shared significant = ",
  shared_total,
  "\n",
  sep = ""
)


cat(
  "Shared concordant = ",
  shared_concordant,
  "\n",
  sep = ""
)


cat(
  "Shared discordant = ",
  shared_discordant,
  "\n",
  sep = ""
)


cat(
  "Batch-adjusted only = ",
  batch_only,
  "\n",
  sep = ""
)


cat(
  "Union = ",
  union_deg,
  "\n",
  sep = ""
)


cat("\nEFFECT-SIZE CONCORDANCE — FINAL BIOLOGICAL TARGETS\n")
cat("--------------------------------------------------\n")


cat(
  "Pearson r = ",
  sprintf(
    "%.7f",
    pearson_r
  ),
  "\n",
  sep = ""
)


cat(
  "Spearman rho = ",
  sprintf(
    "%.7f",
    spearman_rho
  ),
  "\n",
  sep = ""
)


cat(
  "Manuscript rounded values: r = ",
  sprintf(
    "%.3f",
    pearson_r
  ),
  "; rho = ",
  sprintf(
    "%.3f",
    spearman_rho
  ),
  "\n",
  sep = ""
)


cat("\nSEX-LINKED QC\n")
cat("-------------\n")


print(
  sex_linked_audit,
  row.names = FALSE
)


cat("\nFUNCTIONAL ENRICHMENT SOURCE\n")
cat("----------------------------\n")


cat(
  "GO BP significant ORA rows = ",
  nrow(
    go_table
  ),
  "\n",
  sep = ""
)


cat(
  "KEGG significant ORA rows = ",
  nrow(
    kegg_table
  ),
  "\n",
  sep = ""
)


cat(
  "WikiPathways significant ORA rows = ",
  nrow(
    wiki_table
  ),
  "\n",
  sep = ""
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Figure S2 PDF:\n  ",
  normalizePath(
    figure_pdf,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S2 PNG:\n  ",
  normalizePath(
    figure_png,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S2 TIFF:\n  ",
  normalizePath(
    figure_tiff,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure source data:\n  ",
  normalizePath(
    source_data_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure legend:\n  ",
  normalizePath(
    legend_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Proposed Results 3.2 text:\n  ",
  normalizePath(
    results_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nINTERPRETATION GUARDRAILS\n")
cat("-------------------------\n")


cat(
  "- Figure S2 uses the final 12,393 biological-target Table S2.\n"
)


cat(
  "- Biological-only concordance is r ~0.815 and rho ~0.860.\n"
)


cat(
  "- Historical 12,400-feature concordance differed only because seven ERCC features were included.\n"
)


cat(
  "- No DESeq2, ORA, or GSEA analysis is re-run.\n"
)


cat(
  "- Batch-adjusted analysis remains a sensitivity analysis.\n"
)


cat(
  "- All 1,796 genes significant in both models have concordant direction.\n"
)


cat(
  "- None of the predefined sex-linked QC genes enters the robust core.\n"
)


cat(
  "- ORA panels display selected top terms; all significant ORA terms remain in Table S3.\n"
)


cat(
  "- Hallmark GSEA is not duplicated because it is already shown in Main Figure 1E.\n"
)


cat("\nDone.\n")