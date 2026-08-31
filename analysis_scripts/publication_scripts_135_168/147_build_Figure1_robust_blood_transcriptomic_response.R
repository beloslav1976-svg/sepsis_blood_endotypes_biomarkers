################################################################################
# Script 147
# FINAL Main Figure 1
# Robust blood transcriptional response to sepsis
#
# Project: Sepsis_DESeq2
#
# FINAL FIGURE
# ------------
# A. PCA: sepsis vs healthy controls
# B. Primary BP vs BC volcano plot
# C. Primary vs batch-adjusted log2FC concordance
# D. Curated robust-core expression heatmap
# E. Hallmark GSEA of the complete ranked gene list
#
#
# STATUS
# ------
# Publication-packaging / visual-polish script.
#
# IMPORTANT
# ---------
# This script does NOT:
#   - rerun DESeq2;
#   - redefine DEG;
#   - rerun the batch-adjusted DESeq2 model;
#   - redefine the robust-core DEG set;
#   - rerun enrichment;
#   - rerun Hallmark GSEA;
#   - perform feature selection;
#   - perform new hypothesis testing.
#
#
# FINAL VISUAL CHANGES
# --------------------
# Panel B:
#   - reduced to 10 biologically informative gene labels;
#   - includes all five genes of the primary host-response signature.
#
# Panel D:
#   - removed redundant Program annotation/legend;
#   - retained 15 curated UP + 15 curated DOWN genes;
#   - added horizontal gap between UP and DOWN blocks.
#
# Panel E:
#   - wording changed from "complete ranked transcriptome"
#     to "complete ranked gene list";
#   - Gene set size and FDR legends placed on separate rows.
#
# Main figure:
#   - no global title/subtitle inside the image;
#   - scientific title retained in the figure caption.
#
#
# PRIMARY EXPECTED NUMBERS
# ------------------------
# Primary DEG:
#   total  = 2659
#   UP     = 1660
#   DOWN   = 999
#
# Batch-adjusted DEG:
#   total  = 4125
#   UP     = 2093
#   DOWN   = 2032
#
# Robust core:
#   total  = 1796
#   UP     = 1133
#   DOWN   = 663
#
# Effect-size concordance:
#   Pearson r    ~ 0.815
#   Spearman rho ~ 0.859
#
# PCA:
#   PC1 ~ 33.3%
#   PC2 ~ 12.2%
#
# Hallmark:
#   50 gene sets
#   43 positive NES
#   7 negative NES
#   33 positive NES with FDR < 0.25
#   5 negative NES with FDR < 0.25
#
#
# OUTPUT
# ------
# results/blood_endotypes_biomarkers/
#   147_Figure1_robust_blood_response/
#
################################################################################


cat("====================================================================\n")
cat("Running Script 147\n")
cat("FINAL Main Figure 1: robust blood transcriptional response to sepsis\n")
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
  "DESeq2",
  "SummarizedExperiment",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "patchwork",
  "openxlsx",
  "stringr",
  "forcats",
  "grid",
  "scales"
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
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


suppressPackageStartupMessages({
  
  library(DESeq2)
  library(SummarizedExperiment)
  
  library(dplyr)
  library(tidyr)
  
  library(ggplot2)
  library(ggrepel)
  
  library(pheatmap)
  library(patchwork)
  
  library(openxlsx)
  
  library(stringr)
  library(forcats)
  
  library(grid)
  library(scales)
  
})


# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================


to_logical_safe <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  x_chr <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  x_chr %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}


safe_p <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  x[
    is.na(x) |
      !is.finite(x)
  ] <- 1
  
  pmax(
    x,
    .Machine$double.xmin
  )
}


capitalize_condition <- function(x) {
  
  x <- as.character(x)
  
  dplyr::recode(
    x,
    "control" = "Healthy control",
    "sepsis" = "Sepsis",
    "BC" = "Healthy control",
    "BP" = "Sepsis",
    .default = x
  )
}


theme_publication <- function(
    base_size = 10
) {
  
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    
    ggplot2::theme(
      
      plot.title =
        ggplot2::element_text(
          face = "bold",
          size = base_size + 1
        ),
      
      plot.subtitle =
        ggplot2::element_text(
          size = base_size - 0.5,
          color = "grey25"
        ),
      
      plot.tag =
        ggplot2::element_text(
          face = "bold",
          size = base_size + 2
        ),
      
      axis.title =
        ggplot2::element_text(
          face = "bold"
        ),
      
      legend.title =
        ggplot2::element_text(
          face = "bold"
        ),
      
      plot.margin =
        ggplot2::margin(
          7,
          7,
          7,
          7
        )
    )
}


# -----------------------------------------------------------------------------
# Find finalized upstream file
# -----------------------------------------------------------------------------

find_project_file <- function(
    candidates,
    recursive_pattern = NULL,
    description = "file"
) {
  
  candidates <- unique(
    candidates
  )
  
  existing <- candidates[
    file.exists(
      candidates
    )
  ]
  
  
  if (length(existing) > 0) {
    
    return(
      existing[1]
    )
  }
  
  
  if (!is.null(recursive_pattern)) {
    
    hits <- list.files(
      path = "results",
      pattern = recursive_pattern,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    hits <- sort(
      unique(
        hits
      )
    )
    
    
    if (length(hits) > 0) {
      
      if (length(hits) > 1) {
        
        cat(
          "\nMultiple candidate files found for ",
          description,
          ":\n",
          sep = ""
        )
        
        for (h in hits) {
          cat("  ", h, "\n", sep = "")
        }
        
        cat(
          "Using first candidate:\n  ",
          hits[1],
          "\n",
          sep = ""
        )
      }
      
      return(
        hits[1]
      )
    }
  }
  
  
  NA_character_
}


# -----------------------------------------------------------------------------
# Flexible Excel worksheet finder
# -----------------------------------------------------------------------------

find_sheet <- function(
    xlsx_file,
    exact = NULL,
    pattern = NULL,
    required = TRUE
) {
  
  sheets <- openxlsx::getSheetNames(
    xlsx_file
  )
  
  
  if (!is.null(exact)) {
    
    idx_exact <- which(
      tolower(sheets) ==
        tolower(exact)
    )
    
    if (length(idx_exact) > 0) {
      
      return(
        sheets[idx_exact[1]]
      )
    }
  }
  
  
  if (!is.null(pattern)) {
    
    idx_pattern <- grep(
      pattern,
      sheets,
      ignore.case = TRUE
    )
    
    if (length(idx_pattern) > 0) {
      
      return(
        sheets[idx_pattern[1]]
      )
    }
  }
  
  
  if (required) {
    
    stop(
      paste0(
        "Required worksheet not found in:\n",
        xlsx_file,
        "\n\nAvailable sheets:\n",
        paste(
          sheets,
          collapse = ", "
        )
      )
    )
  }
  
  
  NA_character_
}


# -----------------------------------------------------------------------------
# Hallmark pathway labels
# -----------------------------------------------------------------------------

clean_hallmark_label <- function(x) {
  
  label_map <- c(
    
    "HALLMARK_IL6_JAK_STAT3_SIGNALING" =
      "IL6/JAK/STAT3 signaling",
    
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB" =
      "TNFA signaling via NF-\u03baB",
    
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" =
      "Reactive oxygen species pathway",
    
    "HALLMARK_INFLAMMATORY_RESPONSE" =
      "Inflammatory response",
    
    "HALLMARK_COAGULATION" =
      "Coagulation",
    
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS" =
      "Cholesterol homeostasis",
    
    "HALLMARK_INTERFERON_ALPHA_RESPONSE" =
      "Interferon alpha response",
    
    "HALLMARK_COMPLEMENT" =
      "Complement",
    
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING" =
      "PI3K/AKT/mTOR signaling",
    
    "HALLMARK_HYPOXIA" =
      "Hypoxia",
    
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" =
      "Unfolded protein response",
    
    "HALLMARK_E2F_TARGETS" =
      "E2F targets",
    
    "HALLMARK_ALLOGRAFT_REJECTION" =
      "Allograft rejection",
    
    "HALLMARK_MYC_TARGETS_V1" =
      "MYC targets V1",
    
    "HALLMARK_MYC_TARGETS_V2" =
      "MYC targets V2"
  )
  
  
  x_chr <- as.character(x)
  
  mapped <- unname(
    label_map[x_chr]
  )
  
  
  fallback <- x_chr %>%
    stringr::str_replace(
      "^HALLMARK_",
      ""
    ) %>%
    stringr::str_replace_all(
      "_",
      " "
    ) %>%
    stringr::str_to_lower() %>%
    stringr::str_to_sentence()
  
  
  use_map <- !is.na(mapped)
  
  fallback[use_map] <- mapped[use_map]
  
  fallback
}


# =============================================================================
# 4. LOCATE INPUT FILES
# =============================================================================

source_table_dir <- file.path(
  "results",
  "manuscript_source_tables"
)


# -----------------------------------------------------------------------------
# S2: DESeq2 model comparison
# -----------------------------------------------------------------------------

table_s2 <- find_project_file(
  
  candidates = c(
    
    file.path(
      source_table_dir,
      "Table_S2_blood_DESeq2_model_comparison.xlsx"
    ),
    
    file.path(
      "results",
      "blood_BP_vs_BC_model_comparison",
      "blood_model_comparison_all_genes.xlsx"
    ),
    
    file.path(
      "results",
      "blood_BP_vs_BC_model_comparison",
      "blood_model_comparison_all_genes.csv"
    )
  ),
  
  recursive_pattern =
    "blood.*model.*comparison.*all.*genes.*\\.(csv|xlsx)$",
  
  description =
    "Table S2 / blood DESeq2 model comparison"
)


if (
  is.na(table_s2) ||
  !file.exists(table_s2)
) {
  
  stop(
    paste0(
      "\nCould not locate the finalized blood model-comparison source.\n\n",
      "Expected e.g.:\n",
      "results/blood_BP_vs_BC_model_comparison/",
      "blood_model_comparison_all_genes.csv\n"
    )
  )
}


table_s2_type <- tolower(
  tools::file_ext(
    table_s2
  )
)


# -----------------------------------------------------------------------------
# S3: curated robust-core gene table
# -----------------------------------------------------------------------------

table_s3 <- find_project_file(
  
  candidates = c(
    
    file.path(
      source_table_dir,
      "Table_S3_blood_top_core_gene_table_for_article.xlsx"
    ),
    
    file.path(
      "results",
      "blood_core_gene_table",
      "blood_top_core_gene_table_for_article.xlsx"
    )
  ),
  
  recursive_pattern =
    "blood.*top.*core.*gene.*table.*article.*\\.xlsx$",
  
  description =
    "Table S3 / curated robust-core gene table"
)


if (
  is.na(table_s3) ||
  !file.exists(table_s3)
) {
  
  stop(
    "Could not locate Table S3 / curated robust-core gene table."
  )
}


# -----------------------------------------------------------------------------
# S4: sex-linked robust-core audit
# Optional for figure generation
# -----------------------------------------------------------------------------

table_s4 <- find_project_file(
  
  candidates = c(
    
    file.path(
      source_table_dir,
      "Table_S4_blood_core_DEG_sex_linked_filter.xlsx"
    ),
    
    file.path(
      "results",
      "blood_QC_publication",
      "sex_linked_filter",
      "blood_core_DEG_sex_linked_filter.xlsx"
    )
  ),
  
  recursive_pattern =
    "blood.*core.*DEG.*sex.*linked.*filter.*\\.xlsx$",
  
  description =
    "Table S4 / sex-linked robust-core audit"
)


# -----------------------------------------------------------------------------
# S6: Hallmark GSEA
# -----------------------------------------------------------------------------

table_s6 <- find_project_file(
  
  candidates = c(
    
    file.path(
      source_table_dir,
      "Table_S6_Hallmark_GSEA_blood_results.xlsx"
    ),
    
    file.path(
      "results",
      "blood_GSEA_Hallmark_publication_20b_FINAL",
      "20b_FINAL_Hallmark_GSEA_blood_results.xlsx"
    )
  ),
  
  recursive_pattern =
    "Hallmark.*GSEA.*blood.*results.*\\.xlsx$",
  
  description =
    "Table S6 / Hallmark GSEA"
)


if (
  is.na(table_s6) ||
  !file.exists(table_s6)
) {
  
  stop(
    "Could not locate Table S6 / Hallmark GSEA workbook."
  )
}


# -----------------------------------------------------------------------------
# VST object
# -----------------------------------------------------------------------------

vsd_file <- find_project_file(
  
  candidates = c(
    
    file.path(
      "results",
      "blood_QC_publication",
      "blood_BP_vs_BC_vsd.rds"
    ),
    
    file.path(
      "results",
      "blood_BP_vs_BC",
      "blood_BP_vs_BC_vsd.rds"
    )
  ),
  
  recursive_pattern =
    "blood_BP_vs_BC_vsd\\.rds$",
  
  description =
    "blood VST object"
)


if (
  is.na(vsd_file) ||
  !file.exists(vsd_file)
) {
  
  stop(
    "Could not locate blood_BP_vs_BC_vsd.rds."
  )
}


cat("\n====================================================================\n")
cat("FIGURE 1 SOURCE FILES\n")
cat("====================================================================\n")


cat(
  "S2 source type: ",
  table_s2_type,
  "\n",
  sep = ""
)


cat(
  "S2: ",
  normalizePath(
    table_s2,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "S3: ",
  normalizePath(
    table_s3,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


if (
  !is.na(table_s4) &&
  file.exists(table_s4)
) {
  
  cat(
    "S4: ",
    normalizePath(
      table_s4,
      winslash = "\\",
      mustWork = TRUE
    ),
    "\n",
    sep = ""
  )
  
} else {
  
  cat(
    "S4: not found; sex-linked robust-core audit will be reported as NA.\n"
  )
}


cat(
  "S6: ",
  normalizePath(
    table_s6,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "VST: ",
  normalizePath(
    vsd_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 5. OUTPUT DIRECTORIES
# =============================================================================

out_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "147_Figure1_robust_blood_response"
)


figure_dir <- file.path(
  out_dir,
  "figures"
)


table_dir <- file.path(
  out_dir,
  "tables"
)


text_dir <- file.path(
  out_dir,
  "text"
)


for (d in c(
  out_dir,
  figure_dir,
  table_dir,
  text_dir
)) {
  
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


cat("\nOutput directory:\n")

print(
  normalizePath(
    out_dir,
    winslash = "\\",
    mustWork = FALSE
  )
)


# =============================================================================
# 6. EXPECTED VALUES
# =============================================================================

expected_primary_total <- 2659
expected_primary_up <- 1660
expected_primary_down <- 999

expected_batch_total <- 4125
expected_batch_up <- 2093
expected_batch_down <- 2032

expected_core_total <- 1796
expected_core_up <- 1133
expected_core_down <- 663

expected_pearson <- 0.815
expected_spearman <- 0.859

expected_hallmark_total <- 50


# =============================================================================
# 7. PUBLICATION COLORS
# =============================================================================

col_sepsis <- "#D55E00"
col_control <- "#0072B2"

col_core_up <- "#D55E00"
col_core_down <- "#0072B2"

col_sig_noncore <- "#8C8C8C"
col_nonsig <- "#D9D9D9"

col_fdr_strict <- "#0072B2"
col_fdr_relaxed <- "#E69F00"


# =============================================================================
# 8. LOAD TABLE S2
# =============================================================================

cat("\n====================================================================\n")
cat("LOADING TABLE S2 / MODEL COMPARISON\n")
cat("====================================================================\n")


if (table_s2_type == "xlsx") {
  
  s2_sheets <- openxlsx::getSheetNames(
    table_s2
  )
  
  cat("Available S2 sheets:\n")
  print(s2_sheets)
  
  
  s2_main_sheet <- find_sheet(
    table_s2,
    exact = "all_genes_comparison",
    pattern = "all.*genes",
    required = FALSE
  )
  
  
  if (is.na(s2_main_sheet)) {
    
    s2_main_sheet <- s2_sheets[1]
    
    warning(
      "Could not identify all_genes_comparison sheet; using first sheet: ",
      s2_main_sheet
    )
  }
  
  
  s2_all <- openxlsx::read.xlsx(
    table_s2,
    sheet = s2_main_sheet
  )
  
  
} else if (table_s2_type == "csv") {
  
  s2_all <- read.csv(
    table_s2,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
} else {
  
  stop(
    "Unsupported Table S2 file type: ",
    table_s2_type
  )
}


required_s2_columns <- c(
  "Gene",
  "log2FC_simple",
  "padj_simple",
  "log2FC_batch",
  "padj_batch",
  "DEG_simple",
  "DEG_batch",
  "core_DEG"
)


missing_s2_columns <- setdiff(
  required_s2_columns,
  colnames(
    s2_all
  )
)


if (length(missing_s2_columns) > 0) {
  
  stop(
    paste0(
      "Table S2 is missing required columns:\n",
      paste(
        missing_s2_columns,
        collapse = ", "
      ),
      "\n\nObserved columns:\n",
      paste(
        colnames(
          s2_all
        ),
        collapse = ", "
      )
    )
  )
}


s2_all <- s2_all %>%
  
  dplyr::mutate(
    
    Gene =
      as.character(
        Gene
      ),
    
    log2FC_simple =
      suppressWarnings(
        as.numeric(
          log2FC_simple
        )
      ),
    
    padj_simple =
      suppressWarnings(
        as.numeric(
          padj_simple
        )
      ),
    
    log2FC_batch =
      suppressWarnings(
        as.numeric(
          log2FC_batch
        )
      ),
    
    padj_batch =
      suppressWarnings(
        as.numeric(
          padj_batch
        )
      ),
    
    DEG_simple =
      to_logical_safe(
        DEG_simple
      ),
    
    DEG_batch =
      to_logical_safe(
        DEG_batch
      ),
    
    core_DEG =
      to_logical_safe(
        core_DEG
      )
  )


# Frozen robust-core subsets.
# No DEG criteria are recalculated.

s2_core <- s2_all %>%
  dplyr::filter(
    core_DEG
  )


s2_core_up <- s2_core %>%
  dplyr::filter(
    log2FC_simple > 0
  )


s2_core_down <- s2_core %>%
  dplyr::filter(
    log2FC_simple < 0
  )


cat("\nTable S2 loaded successfully:\n")

cat(
  "  All genes: ",
  nrow(s2_all),
  "\n",
  sep = ""
)

cat(
  "  Frozen robust core: ",
  nrow(s2_core),
  "\n",
  sep = ""
)

cat(
  "  Frozen robust UP: ",
  nrow(s2_core_up),
  "\n",
  sep = ""
)

cat(
  "  Frozen robust DOWN: ",
  nrow(s2_core_down),
  "\n",
  sep = ""
)


# =============================================================================
# 9. LOAD TABLE S3
# =============================================================================

cat("\n====================================================================\n")
cat("LOADING TABLE S3 / CURATED ROBUST-CORE GENES\n")
cat("====================================================================\n")


s3_sheets <- openxlsx::getSheetNames(
  table_s3
)

cat("Available S3 sheets:\n")
print(s3_sheets)


s3_up_sheet <- find_sheet(
  table_s3,
  exact = "top50_core_UP",
  pattern = "top.*core.*up"
)


s3_down_sheet <- find_sheet(
  table_s3,
  exact = "top50_core_DOWN",
  pattern = "top.*core.*down"
)


s3_top_up <- openxlsx::read.xlsx(
  table_s3,
  sheet = s3_up_sheet
)


s3_top_down <- openxlsx::read.xlsx(
  table_s3,
  sheet = s3_down_sheet
)


if (
  !"Gene" %in%
  colnames(s3_top_up)
) {
  
  stop(
    "Gene column missing from S3 UP sheet."
  )
}


if (
  !"Gene" %in%
  colnames(s3_top_down)
) {
  
  stop(
    "Gene column missing from S3 DOWN sheet."
  )
}


# =============================================================================
# 10. LOAD TABLE S4 — OPTIONAL SEX-LINKED AUDIT
# =============================================================================

sex_linked_core_n <- NA_integer_

s4_summary <- data.frame()
s4_sex_core <- data.frame()


if (
  !is.na(table_s4) &&
  file.exists(table_s4)
) {
  
  cat("\n====================================================================\n")
  cat("LOADING TABLE S4 / SEX-LINKED ROBUST-CORE AUDIT\n")
  cat("====================================================================\n")
  
  
  s4_sheets <- openxlsx::getSheetNames(
    table_s4
  )
  
  cat("Available S4 sheets:\n")
  print(s4_sheets)
  
  
  s4_summary_sheet <- find_sheet(
    table_s4,
    exact = "summary",
    pattern = "^summary$",
    required = FALSE
  )
  
  
  s4_sex_core_sheet <- find_sheet(
    table_s4,
    exact = "sex_linked_genes_present_in_c",
    pattern = "sex.*linked.*present",
    required = FALSE
  )
  
  
  if (!is.na(s4_summary_sheet)) {
    
    s4_summary <- openxlsx::read.xlsx(
      table_s4,
      sheet = s4_summary_sheet
    )
  }
  
  
  if (!is.na(s4_sex_core_sheet)) {
    
    s4_sex_core <- openxlsx::read.xlsx(
      table_s4,
      sheet = s4_sex_core_sheet
    )
    
    sex_linked_core_n <- nrow(
      s4_sex_core
    )
    
  } else {
    
    warning(
      "Explicit sex-linked robust-core sheet not found in Table S4."
    )
  }
}


# =============================================================================
# 11. LOAD TABLE S6 — HALLMARK GSEA
# =============================================================================

cat("\n====================================================================\n")
cat("LOADING TABLE S6 / HALLMARK GSEA\n")
cat("====================================================================\n")


s6_sheets <- openxlsx::getSheetNames(
  table_s6
)

cat("Available S6 sheets:\n")
print(s6_sheets)


s6_all_sheet <- find_sheet(
  table_s6,
  exact = "all_Hallmark_GSEA",
  pattern = "all.*hallmark.*gsea",
  required = FALSE
)


if (is.na(s6_all_sheet)) {
  
  s6_all_sheet <- find_sheet(
    table_s6,
    pattern = "hallmark",
    required = FALSE
  )
}


if (is.na(s6_all_sheet)) {
  
  stop(
    paste0(
      "Could not identify Hallmark results sheet in Table S6.\n",
      "Available sheets: ",
      paste(
        s6_sheets,
        collapse = ", "
      )
    )
  )
}


s6_all <- openxlsx::read.xlsx(
  table_s6,
  sheet = s6_all_sheet
)


# -----------------------------------------------------------------------------
# Resolve Hallmark pathway ID column
# -----------------------------------------------------------------------------

if (
  !"ID" %in%
  colnames(s6_all)
) {
  
  possible_id <- grep(
    "^(ID|pathway|gene.?set|Description)$",
    colnames(s6_all),
    ignore.case = TRUE,
    value = TRUE
  )
  
  
  if (length(possible_id) >= 1) {
    
    s6_all$ID <- s6_all[[possible_id[1]]]
    
    cat(
      "Hallmark ID column resolved from: ",
      possible_id[1],
      "\n",
      sep = ""
    )
    
  } else {
    
    stop(
      paste0(
        "Could not identify Hallmark pathway ID column.\n",
        "Observed columns:\n",
        paste(
          colnames(s6_all),
          collapse = ", "
        )
      )
    )
  }
}


# -----------------------------------------------------------------------------
# Resolve adjusted p-value / FDR column
# -----------------------------------------------------------------------------

if (
  !"p.adjust" %in%
  colnames(s6_all)
) {
  
  possible_padj <- grep(
    "^(p\\.adjust|padj|FDR|qvalue)$",
    colnames(s6_all),
    ignore.case = TRUE,
    value = TRUE
  )
  
  
  if (length(possible_padj) >= 1) {
    
    s6_all$p.adjust <- s6_all[[possible_padj[1]]]
    
    cat(
      "Hallmark adjusted-p column resolved from: ",
      possible_padj[1],
      "\n",
      sep = ""
    )
    
  } else {
    
    stop(
      paste0(
        "Could not identify Hallmark adjusted p-value/FDR column.\n",
        "Observed columns:\n",
        paste(
          colnames(s6_all),
          collapse = ", "
        )
      )
    )
  }
}


# -----------------------------------------------------------------------------
# Resolve gene-set size column
# -----------------------------------------------------------------------------

if (
  !"setSize" %in%
  colnames(s6_all)
) {
  
  possible_size <- grep(
    "set.?size|geneset.?size|size",
    colnames(s6_all),
    ignore.case = TRUE,
    value = TRUE
  )
  
  
  if (length(possible_size) >= 1) {
    
    s6_all$setSize <- s6_all[[possible_size[1]]]
    
    cat(
      "Hallmark gene-set-size column resolved from: ",
      possible_size[1],
      "\n",
      sep = ""
    )
    
  } else {
    
    s6_all$setSize <- 1
    
    warning(
      paste0(
        "No Hallmark setSize column found. ",
        "A constant point size will be used for visualization only."
      )
    )
  }
}


if (
  !"NES" %in%
  colnames(s6_all)
) {
  
  stop(
    "NES column missing from Hallmark results table."
  )
}


s6_all <- s6_all %>%
  
  dplyr::mutate(
    
    ID =
      as.character(
        ID
      ),
    
    NES =
      suppressWarnings(
        as.numeric(
          NES
        )
      ),
    
    p.adjust =
      suppressWarnings(
        as.numeric(
          p.adjust
        )
      ),
    
    setSize =
      suppressWarnings(
        as.numeric(
          setSize
        )
      )
  )


# =============================================================================
# 12. LOAD VST OBJECT
# =============================================================================

vsd <- readRDS(
  vsd_file
)


vsd_mat <- SummarizedExperiment::assay(
  vsd
)


col_data <- as.data.frame(
  SummarizedExperiment::colData(
    vsd
  )
)


if (
  !"condition" %in%
  colnames(col_data)
) {
  
  stop(
    "VST object does not contain 'condition' metadata."
  )
}


# =============================================================================
# 13. NUMERICAL AUDIT — DIFFERENTIAL EXPRESSION
# =============================================================================

primary_up_n <- s2_all %>%
  dplyr::filter(
    DEG_simple,
    log2FC_simple > 0
  ) %>%
  nrow()


primary_down_n <- s2_all %>%
  dplyr::filter(
    DEG_simple,
    log2FC_simple < 0
  ) %>%
  nrow()


primary_total_n <-
  primary_up_n +
  primary_down_n


batch_up_n <- s2_all %>%
  dplyr::filter(
    DEG_batch,
    log2FC_batch > 0
  ) %>%
  nrow()


batch_down_n <- s2_all %>%
  dplyr::filter(
    DEG_batch,
    log2FC_batch < 0
  ) %>%
  nrow()


batch_total_n <-
  batch_up_n +
  batch_down_n


core_total_n <- nrow(
  s2_core
)


core_up_n <- nrow(
  s2_core_up
)


core_down_n <- nrow(
  s2_core_down
)


fc_complete <- s2_all %>%
  dplyr::filter(
    is.finite(log2FC_simple),
    is.finite(log2FC_batch)
  )


pearson_r <- stats::cor(
  fc_complete$log2FC_simple,
  fc_complete$log2FC_batch,
  method = "pearson"
)


spearman_rho <- stats::cor(
  fc_complete$log2FC_simple,
  fc_complete$log2FC_batch,
  method = "spearman"
)


# =============================================================================
# 14. NUMERICAL AUDIT — HALLMARK
# =============================================================================

hallmark_total_n <- nrow(
  s6_all
)


hallmark_positive_n <- sum(
  s6_all$NES > 0,
  na.rm = TRUE
)


hallmark_negative_n <- sum(
  s6_all$NES < 0,
  na.rm = TRUE
)


hallmark_positive_fdr025_n <- sum(
  s6_all$NES > 0 &
    s6_all$p.adjust < 0.25,
  na.rm = TRUE
)


hallmark_negative_fdr025_n <- sum(
  s6_all$NES < 0 &
    s6_all$p.adjust < 0.25,
  na.rm = TRUE
)


audit_table <- data.frame(
  
  metric = c(
    "Primary DEG total",
    "Primary DEG UP",
    "Primary DEG DOWN",
    "Batch-adjusted DEG total",
    "Batch-adjusted DEG UP",
    "Batch-adjusted DEG DOWN",
    "Robust core total",
    "Robust core UP",
    "Robust core DOWN",
    "Sex-linked genes in robust core",
    "Pearson log2FC correlation",
    "Spearman log2FC correlation",
    "Hallmark pathways tested",
    "Hallmark positive NES",
    "Hallmark negative NES",
    "Hallmark positive NES FDR<0.25",
    "Hallmark negative NES FDR<0.25"
  ),
  
  observed = c(
    primary_total_n,
    primary_up_n,
    primary_down_n,
    batch_total_n,
    batch_up_n,
    batch_down_n,
    core_total_n,
    core_up_n,
    core_down_n,
    sex_linked_core_n,
    pearson_r,
    spearman_rho,
    hallmark_total_n,
    hallmark_positive_n,
    hallmark_negative_n,
    hallmark_positive_fdr025_n,
    hallmark_negative_fdr025_n
  ),
  
  expected = c(
    expected_primary_total,
    expected_primary_up,
    expected_primary_down,
    expected_batch_total,
    expected_batch_up,
    expected_batch_down,
    expected_core_total,
    expected_core_up,
    expected_core_down,
    0,
    expected_pearson,
    expected_spearman,
    expected_hallmark_total,
    NA,
    NA,
    NA,
    NA
  ),
  
  stringsAsFactors = FALSE
)


audit_table$difference <-
  audit_table$observed -
  audit_table$expected


cat("\n====================================================================\n")
cat("NUMERICAL AUDIT\n")
cat("====================================================================\n")

print(
  audit_table
)


# -----------------------------------------------------------------------------
# Hard checks
# -----------------------------------------------------------------------------

hard_checks <- c(
  
  primary_total_n ==
    expected_primary_total,
  
  primary_up_n ==
    expected_primary_up,
  
  primary_down_n ==
    expected_primary_down,
  
  batch_total_n ==
    expected_batch_total,
  
  batch_up_n ==
    expected_batch_up,
  
  batch_down_n ==
    expected_batch_down,
  
  core_total_n ==
    expected_core_total,
  
  core_up_n ==
    expected_core_up,
  
  core_down_n ==
    expected_core_down
)


if (!all(hard_checks)) {
  
  stop(
    paste0(
      "\nCRITICAL DESeq2/core numerical audit failed.\n",
      "Figure 1 was NOT generated.\n",
      "Resolve source-table provenance before proceeding."
    )
  )
}


if (
  abs(
    pearson_r -
    expected_pearson
  ) > 0.01
) {
  
  warning(
    "Pearson correlation differs from expected ~0.815. Observed: ",
    sprintf(
      "%.4f",
      pearson_r
    )
  )
}


if (
  abs(
    spearman_rho -
    expected_spearman
  ) > 0.01
) {
  
  warning(
    "Spearman correlation differs from expected ~0.859. Observed: ",
    sprintf(
      "%.4f",
      spearman_rho
    )
  )
}


if (
  hallmark_total_n !=
  expected_hallmark_total
) {
  
  warning(
    "Hallmark pathway count differs from expected 50. Observed: ",
    hallmark_total_n,
    ". Final S6 table will be treated as source of truth."
  )
}


# =============================================================================
# 15. PANEL A — PCA
# =============================================================================

pca_data <- DESeq2::plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)


percent_var <- round(
  100 *
    attr(
      pca_data,
      "percentVar"
    ),
  1
)


pca_data <- pca_data %>%
  
  dplyr::mutate(
    
    Condition =
      capitalize_condition(
        condition
      ),
    
    Condition =
      factor(
        Condition,
        levels = c(
          "Healthy control",
          "Sepsis"
        )
      )
  )


p_A <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    color = Condition
  )
) +
  
  ggplot2::geom_point(
    size = 3.1,
    alpha = 0.92
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "Healthy control" =
        col_control,
      "Sepsis" =
        col_sepsis
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  
  ggplot2::labs(
    
    tag = "A",
    
    title =
      "Global blood transcriptomic structure",
    
    subtitle =
      "Variance-stabilized expression; n=45",
    
    x = paste0(
      "PC1 (",
      percent_var[1],
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      percent_var[2],
      "%)"
    ),
    
    color = NULL
  )


# =============================================================================
# 16. PANEL B — PRIMARY VOLCANO
# =============================================================================

volcano_df <- s2_all %>%
  
  dplyr::mutate(
    
    padj_plot =
      safe_p(
        padj_simple
      ),
    
    neg_log10_padj =
      -log10(
        padj_plot
      ),
    
    # Visualization cap only.
    neg_log10_padj_plot =
      pmin(
        neg_log10_padj,
        42
      ),
    
    volcano_group =
      dplyr::case_when(
        
        core_DEG &
          log2FC_simple > 0 ~
          "Core up",
        
        core_DEG &
          log2FC_simple < 0 ~
          "Core down",
        
        DEG_simple &
          !core_DEG ~
          "Significant non-core",
        
        TRUE ~
          "Not significant"
      ),
    
    volcano_group =
      factor(
        volcano_group,
        levels = c(
          "Not significant",
          "Significant non-core",
          "Core down",
          "Core up"
        )
      )
  )


# Final reduced label set.
# Contains all five genes of the primary host-response signature.

label_genes <- c(
  "GRB10",
  "CD177",
  "ANXA3",
  "IL1R2",
  "S100A12",
  "MMP9",
  "HK3",
  "IRAK3",
  "CARD11",
  "IKZF2"
)


label_df <- volcano_df %>%
  
  dplyr::filter(
    Gene %in% label_genes
  ) %>%
  
  dplyr::filter(
    DEG_simple
  )


p_B <- ggplot2::ggplot(
  volcano_df,
  ggplot2::aes(
    x = log2FC_simple,
    y = neg_log10_padj_plot
  )
) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      color = volcano_group
    ),
    size = 1.0,
    alpha = 0.68
  ) +
  
  ggplot2::geom_vline(
    xintercept = c(
      -1,
      1
    ),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey45"
  ) +
  
  ggplot2::geom_hline(
    yintercept =
      -log10(
        0.05
      ),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey45"
  ) +
  
  ggrepel::geom_text_repel(
    data = label_df,
    ggplot2::aes(
      label = Gene
    ),
    size = 2.9,
    max.overlaps = Inf,
    box.padding = 0.32,
    point.padding = 0.16,
    min.segment.length = 0,
    segment.size = 0.25,
    seed = 123
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      
      "Not significant" =
        col_nonsig,
      
      "Significant non-core" =
        col_sig_noncore,
      
      "Core down" =
        col_core_down,
      
      "Core up" =
        col_core_up
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "bottom",
    
    legend.text =
      ggplot2::element_text(
        size = 8
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "B",
    
    title =
      "Differential-expression landscape",
    
    subtitle =
      paste0(
        primary_total_n,
        " DEG: ",
        primary_up_n,
        " increased, ",
        primary_down_n,
        " decreased"
      ),
    
    x =
      expression(
        log[2] *
          " fold change (sepsis vs control)"
      ),
    
    y =
      expression(
        -log[10] *
          " adjusted p-value"
      ),
    
    color =
      NULL
  )


# =============================================================================
# 17. PANEL C — PRIMARY vs BATCH-ADJUSTED EFFECT SIZES
# =============================================================================

concordance_df <- s2_all %>%
  
  dplyr::filter(
    is.finite(log2FC_simple),
    is.finite(log2FC_batch)
  ) %>%
  
  dplyr::mutate(
    
    category =
      dplyr::case_when(
        
        core_DEG &
          log2FC_simple > 0 ~
          "Core up",
        
        core_DEG &
          log2FC_simple < 0 ~
          "Core down",
        
        TRUE ~
          "Other genes"
      ),
    
    category =
      factor(
        category,
        levels = c(
          "Other genes",
          "Core down",
          "Core up"
        )
      )
  )


# Correlations use all finite genes.
# Axis clipping is visualization only.

x_quant <- stats::quantile(
  concordance_df$log2FC_simple,
  probs = c(
    0.002,
    0.998
  ),
  na.rm = TRUE
)


y_quant <- stats::quantile(
  concordance_df$log2FC_batch,
  probs = c(
    0.002,
    0.995
  ),
  na.rm = TRUE
)


common_lower <- min(
  as.numeric(
    x_quant[1]
  ),
  as.numeric(
    y_quant[1]
  ),
  -4
)


common_upper <- max(
  as.numeric(
    x_quant[2]
  ),
  as.numeric(
    y_quant[2]
  ),
  8
)


p_C <- ggplot2::ggplot(
  concordance_df,
  ggplot2::aes(
    x = log2FC_simple,
    y = log2FC_batch
  )
) +
  
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.65,
    color = "black"
  ) +
  
  ggplot2::geom_point(
    
    data =
      concordance_df %>%
      dplyr::filter(
        category ==
          "Other genes"
      ),
    
    color = "#D3D3D3",
    
    size = 0.75,
    
    alpha = 0.38
  ) +
  
  ggplot2::geom_point(
    
    data =
      concordance_df %>%
      dplyr::filter(
        category !=
          "Other genes"
      ),
    
    ggplot2::aes(
      color = category
    ),
    
    size = 1.15,
    
    alpha = 0.80
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "Core down" =
        col_core_down,
      "Core up" =
        col_core_up
    )
  ) +
  
  ggplot2::coord_cartesian(
    
    xlim = c(
      common_lower,
      common_upper
    ),
    
    ylim = c(
      common_lower,
      common_upper
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  
  ggplot2::annotate(
    
    "label",
    
    x =
      common_lower +
      0.03 *
      (
        common_upper -
          common_lower
      ),
    
    y =
      common_upper -
      0.04 *
      (
        common_upper -
          common_lower
      ),
    
    hjust = 0,
    
    vjust = 1,
    
    label = paste0(
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
    
    size = 3.2,
    
    linewidth = 0.25,
    
    fill = "white"
  ) +
  
  ggplot2::labs(
    
    tag =
      "C",
    
    title =
      "Effect-size robustness to batch adjustment",
    
    subtitle =
      paste0(
        core_total_n,
        " genes retained concordant differential expression"
      ),
    
    x =
      expression(
        "Primary model " *
          log[2] *
          "FC"
      ),
    
    y =
      expression(
        "Batch-adjusted model " *
          log[2] *
          "FC"
      ),
    
    color =
      NULL
  )


# =============================================================================
# 18. PANEL D — FINAL CURATED ROBUST-CORE HEATMAP
# =============================================================================
#
# 15 curated UP + 15 curated DOWN genes from Table S3.
#
# No Program annotation legend in the final main figure.
#
# A horizontal gap separates the UP and DOWN blocks.
#
# =============================================================================

n_heatmap_up <- 15
n_heatmap_down <- 15


if (
  nrow(s3_top_up) <
  n_heatmap_up
) {
  
  stop(
    "S3 contains fewer than 15 curated UP genes."
  )
}


if (
  nrow(s3_top_down) <
  n_heatmap_down
) {
  
  stop(
    "S3 contains fewer than 15 curated DOWN genes."
  )
}


heatmap_up <- s3_top_up %>%
  dplyr::slice_head(
    n = n_heatmap_up
  )


heatmap_down <- s3_top_down %>%
  dplyr::slice_head(
    n = n_heatmap_down
  )


heatmap_genes <- c(
  as.character(
    heatmap_up$Gene
  ),
  as.character(
    heatmap_down$Gene
  )
)


if (
  anyDuplicated(
    heatmap_genes
  ) > 0
) {
  
  stop(
    "Duplicate genes detected in curated heatmap subset."
  )
}


missing_heatmap_genes <- setdiff(
  heatmap_genes,
  rownames(vsd_mat)
)


if (
  length(missing_heatmap_genes) > 0
) {
  
  stop(
    paste0(
      "Curated heatmap genes missing from VST matrix:\n",
      paste(
        missing_heatmap_genes,
        collapse = ", "
      )
    )
  )
}


heatmap_mat <- vsd_mat[
  heatmap_genes,
  ,
  drop = FALSE
]


# Row-wise z-score for visualization.

heatmap_z <- t(
  scale(
    t(
      heatmap_mat
    )
  )
)


heatmap_z[
  !is.finite(
    heatmap_z
  )
] <- 0


annotation_col <- as.data.frame(
  SummarizedExperiment::colData(
    vsd
  )
)


annotation_col$condition <-
  capitalize_condition(
    annotation_col$condition
  )


annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c(
    "Healthy control",
    "Sepsis"
  )
)


annotation_col_plot <- data.frame(
  Condition =
    annotation_col$condition
)


rownames(
  annotation_col_plot
) <- rownames(
  annotation_col
)


annotation_colors <- list(
  
  Condition = c(
    
    "Healthy control" =
      col_control,
    
    "Sepsis" =
      col_sepsis
  )
)


heatmap_colors <- grDevices::colorRampPalette(
  c(
    "#2166AC",
    "#F7F7F7",
    "#B2182B"
  )
)(
  101
)


heatmap_object <- pheatmap::pheatmap(
  
  heatmap_z,
  
  cluster_rows = FALSE,
  
  cluster_cols = TRUE,
  
  annotation_col =
    annotation_col_plot,
  
  annotation_colors =
    annotation_colors,
  
  gaps_row =
    n_heatmap_up,
  
  show_colnames = FALSE,
  
  show_rownames = TRUE,
  
  fontsize = 7,
  
  fontsize_row = 7,
  
  border_color = NA,
  
  color =
    heatmap_colors,
  
  breaks =
    seq(
      -2.5,
      2.5,
      length.out = 102
    ),
  
  main = "",
  
  silent = TRUE
)


p_D_core <- patchwork::wrap_elements(
  full =
    heatmap_object$gtable
)


p_D_title <- ggplot2::ggplot() +
  
  ggplot2::annotate(
    "text",
    x = 0,
    y = 1,
    label =
      "D   Robust-core expression program",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  ggplot2::xlim(
    0,
    1
  ) +
  
  ggplot2::ylim(
    0,
    1
  ) +
  
  ggplot2::theme_void()


p_D <- (
  p_D_title /
    p_D_core
) +
  
  patchwork::plot_layout(
    heights = c(
      0.06,
      0.94
    )
  )


# =============================================================================
# 19. PANEL E — FINAL HALLMARK GSEA
# =============================================================================
#
# FDR < 0.25.
#
# Display:
#   - top 10 positive NES pathways;
#   - up to 10 negative NES pathways.
#
# Ranking uses only existing frozen Table S6 results.
# GSEA is NOT rerun.
#
# =============================================================================

hallmark_positive <- s6_all %>%
  
  dplyr::filter(
    NES > 0,
    p.adjust < 0.25
  ) %>%
  
  dplyr::arrange(
    p.adjust,
    dplyr::desc(
      abs(
        NES
      )
    )
  ) %>%
  
  dplyr::slice_head(
    n = 10
  )


hallmark_negative <- s6_all %>%
  
  dplyr::filter(
    NES < 0,
    p.adjust < 0.25
  ) %>%
  
  dplyr::arrange(
    p.adjust,
    dplyr::desc(
      abs(
        NES
      )
    )
  ) %>%
  
  dplyr::slice_head(
    n = 10
  )


hallmark_plot_df <- dplyr::bind_rows(
  hallmark_negative,
  hallmark_positive
) %>%
  
  dplyr::mutate(
    
    pathway_clean =
      clean_hallmark_label(
        ID
      ),
    
    pathway_label =
      stringr::str_wrap(
        pathway_clean,
        width = 32
      ),
    
    Enrichment =
      dplyr::if_else(
        NES > 0,
        "Sepsis",
        "Control"
      ),
    
    FDR_group =
      dplyr::case_when(
        
        p.adjust < 0.05 ~
          "FDR < 0.05",
        
        p.adjust < 0.25 ~
          "0.05 \u2264 FDR < 0.25",
        
        TRUE ~
          "FDR \u2265 0.25"
      ),
    
    FDR_group =
      factor(
        FDR_group,
        levels = c(
          "FDR < 0.05",
          "0.05 \u2264 FDR < 0.25"
        )
      ),
    
    pathway_label =
      forcats::fct_reorder(
        pathway_label,
        NES
      )
  )


if (
  nrow(hallmark_plot_df) == 0
) {
  
  stop(
    "No Hallmark pathways with FDR < 0.25 were found."
  )
}


size_is_informative <-
  length(
    unique(
      hallmark_plot_df$setSize[
        is.finite(
          hallmark_plot_df$setSize
        )
      ]
    )
  ) > 1


if (size_is_informative) {
  
  p_E <- ggplot2::ggplot(
    hallmark_plot_df,
    ggplot2::aes(
      x = NES,
      y = pathway_label
    )
  ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.55,
      color = "grey35"
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        size = setSize,
        color = FDR_group
      ),
      alpha = 0.95
    ) +
    
    ggplot2::scale_size_continuous(
      range = c(
        2.7,
        7.2
      )
    )
  
  
} else {
  
  p_E <- ggplot2::ggplot(
    hallmark_plot_df,
    ggplot2::aes(
      x = NES,
      y = pathway_label
    )
  ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.55,
      color = "grey35"
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        color = FDR_group
      ),
      size = 4,
      alpha = 0.95
    )
}


p_E <- p_E +
  
  ggplot2::scale_color_manual(
    values = c(
      
      "FDR < 0.05" =
        col_fdr_strict,
      
      "0.05 \u2264 FDR < 0.25" =
        col_fdr_relaxed
    ),
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = seq(
      -3,
      3,
      by = 1
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.04,
        0.04
      )
    )
  ) +
  
  ggplot2::guides(
    
    size =
      ggplot2::guide_legend(
        order = 1,
        nrow = 1,
        title.position = "left"
      ),
    
    color =
      ggplot2::guide_legend(
        order = 2,
        nrow = 1,
        title.position = "left"
      )
  ) +
  
  ggplot2::theme_bw(
    base_size = 9
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 11
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 8.7,
        color = "grey25"
      ),
    
    axis.title.x =
      ggplot2::element_text(
        face = "bold"
      ),
    
    axis.title.y =
      ggplot2::element_blank(),
    
    axis.text.y =
      ggplot2::element_text(
        size = 7.3,
        color = "grey10"
      ),
    
    axis.text.x =
      ggplot2::element_text(
        size = 8.2
      ),
    
    panel.grid.major.y =
      ggplot2::element_blank(),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    panel.grid.major.x =
      ggplot2::element_line(
        color = "grey90",
        linewidth = 0.3
      ),
    
    legend.position =
      "bottom",
    
    legend.box =
      "vertical",
    
    legend.direction =
      "horizontal",
    
    legend.box.just =
      "left",
    
    legend.spacing.y =
      grid::unit(
        0.05,
        "cm"
      ),
    
    legend.title =
      ggplot2::element_text(
        face = "bold",
        size = 8
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 7.5
      ),
    
    plot.tag =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    
    plot.margin =
      ggplot2::margin(
        6,
        6,
        6,
        6
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "E",
    
    title =
      "Coordinated pathway-level host-response state",
    
    subtitle =
      "Hallmark GSEA of the complete ranked gene list",
    
    x =
      "Normalized enrichment score (NES)",
    
    color =
      "FDR",
    
    size =
      "Gene set size"
  )


# =============================================================================
# 20. ASSEMBLE FINAL MAIN FIGURE 1
# =============================================================================
#
# No global title or subtitle inside the figure.
#
# Scientific figure title remains in the caption.
#
# =============================================================================

top_row <- (
  p_A |
    p_B |
    p_C
) +
  
  patchwork::plot_layout(
    widths = c(
      0.85,
      1.15,
      1.00
    )
  )


bottom_row <- (
  p_D |
    p_E
) +
  
  patchwork::plot_layout(
    widths = c(
      1.38,
      1.00
    )
  )


figure1 <- (
  top_row /
    bottom_row
) +
  
  patchwork::plot_layout(
    heights = c(
      0.90,
      1.30
    )
  )


# =============================================================================
# 21. EXPORT FINAL MAIN FIGURE
# =============================================================================

figure_png <- file.path(
  figure_dir,
  "147_Figure1_robust_blood_transcriptomic_response.png"
)


figure_pdf <- file.path(
  figure_dir,
  "147_Figure1_robust_blood_transcriptomic_response.pdf"
)


figure_tiff <- file.path(
  figure_dir,
  "147_Figure1_robust_blood_transcriptomic_response.tiff"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure1,
  width = 15.8,
  height = 10.5,
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure1,
  width = 15.8,
  height = 10.5,
  device =
    if (
      capabilities(
        "cairo"
      )
    ) {
      
      grDevices::cairo_pdf
      
    } else {
      
      grDevices::pdf
    },
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure1,
  width = 15.8,
  height = 10.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 22. EXPORT INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  
  A_PCA =
    p_A,
  
  B_volcano =
    p_B,
  
  C_effect_size_concordance =
    p_C,
  
  E_Hallmark_GSEA =
    p_E
)


for (nm in names(
  individual_panels
)) {
  
  ggplot2::ggsave(
    
    filename =
      file.path(
        figure_dir,
        paste0(
          "147_panel_",
          nm,
          ".png"
        )
      ),
    
    plot =
      individual_panels[[nm]],
    
    width = 6.4,
    
    height = 5.3,
    
    dpi = 600,
    
    bg = "white"
  )
}


# Heatmap exported separately.

grDevices::png(
  
  filename =
    file.path(
      figure_dir,
      "147_panel_D_robust_core_heatmap.png"
    ),
  
  width = 3200,
  
  height = 2800,
  
  res = 450
)


grid::grid.newpage()

grid::grid.draw(
  heatmap_object$gtable
)

grDevices::dev.off()


# =============================================================================
# 23. FIGURE SOURCE DATA
# =============================================================================

pca_source <- pca_data %>%
  
  dplyr::select(
    name,
    PC1,
    PC2,
    condition,
    Condition
  )


volcano_source <- volcano_df %>%
  
  dplyr::select(
    Gene,
    log2FC_simple,
    padj_simple,
    log2FC_batch,
    padj_batch,
    DEG_simple,
    DEG_batch,
    core_DEG,
    volcano_group
  )


concordance_source <- concordance_df %>%
  
  dplyr::select(
    Gene,
    log2FC_simple,
    log2FC_batch,
    DEG_simple,
    DEG_batch,
    core_DEG,
    category
  )


heatmap_gene_source <- dplyr::bind_rows(
  
  heatmap_up %>%
    dplyr::mutate(
      Figure_program =
        "Increased in sepsis"
    ),
  
  heatmap_down %>%
    dplyr::mutate(
      Figure_program =
        "Decreased in sepsis"
    )
)


heatmap_z_source <- data.frame(
  
  Gene =
    rownames(
      heatmap_z
    ),
  
  heatmap_z,
  
  check.names = FALSE
)


source_workbook <- file.path(
  table_dir,
  "147_Figure1_source_data.xlsx"
)


s4_source_string <- if (
  !is.na(table_s4) &&
  file.exists(table_s4)
) {
  
  normalizePath(
    table_s4,
    winslash = "\\",
    mustWork = TRUE
  )
  
} else {
  
  "NOT FOUND / OPTIONAL"
}


openxlsx::write.xlsx(
  
  list(
    
    Run_info =
      data.frame(
        
        item = c(
          "script",
          "status",
          "analysis_mode",
          "DESeq2_rerun",
          "robust_core_redefined",
          "Hallmark_GSEA_rerun",
          "S2_source",
          "S3_source",
          "S4_source",
          "S6_source",
          "VST_source",
          "heatmap_display",
          "heatmap_program_annotation",
          "Hallmark_subtitle",
          "global_figure_title_inside_image",
          "Figure1_architecture"
        ),
        
        value = c(
          "147_build_Figure1_robust_blood_transcriptomic_response.R",
          "FINAL visual-polish version",
          "publication packaging only",
          "NO",
          "NO",
          "NO",
          
          normalizePath(
            table_s2,
            winslash = "\\",
            mustWork = TRUE
          ),
          
          normalizePath(
            table_s3,
            winslash = "\\",
            mustWork = TRUE
          ),
          
          s4_source_string,
          
          normalizePath(
            table_s6,
            winslash = "\\",
            mustWork = TRUE
          ),
          
          normalizePath(
            vsd_file,
            winslash = "\\",
            mustWork = TRUE
          ),
          
          "First 15 UP + first 15 DOWN from curated Table S3",
          "Removed from main figure; UP/DOWN separated by row gap",
          "Hallmark GSEA of the complete ranked gene list",
          "NO",
          paste0(
            "A PCA | B volcano | C model concordance | ",
            "D curated heatmap | E Hallmark GSEA"
          )
        ),
        
        stringsAsFactors = FALSE
      ),
    
    
    Numerical_audit =
      audit_table,
    
    
    PCA =
      pca_source,
    
    
    Volcano =
      volcano_source,
    
    
    Model_concordance =
      concordance_source,
    
    
    Heatmap_genes =
      heatmap_gene_source,
    
    
    Heatmap_z =
      heatmap_z_source,
    
    
    Hallmark_main =
      hallmark_plot_df,
    
    
    Hallmark_all =
      s6_all
  ),
  
  source_workbook,
  
  overwrite = TRUE
)


# =============================================================================
# 24. FINAL FIGURE CAPTION — ENGLISH
# =============================================================================

caption_en <- paste0(
  
  "Figure 1. A robust blood transcriptional response defines an opposing ",
  "myeloid-adaptive host-response axis in sepsis. ",
  
  "(A) Principal-component analysis of variance-stabilized targeted blood ",
  "transcriptomic profiles from patients with sepsis and healthy controls. ",
  
  "(B) Differential-expression landscape from the primary DESeq2 ",
  "sepsis-versus-control model. Differential expression was defined by ",
  "adjusted p<0.05 and |log2 fold change|>=1. The primary analysis identified ",
  primary_total_n,
  " differentially expressed genes, including ",
  primary_up_n,
  " genes increased and ",
  primary_down_n,
  " genes decreased in sepsis. Robust-core genes that also met these criteria ",
  "with concordant direction in the batch-adjusted model are highlighted. ",
  
  "(C) Concordance of gene-level log2 fold-change estimates between the primary ",
  "and sequencing-batch-adjusted DESeq2 models. Correlation coefficients were ",
  "calculated using all genes with finite estimates; the displayed axes show ",
  "the central effect-size range for readability. The dashed diagonal denotes ",
  "identity. Effect estimates were strongly concordant (Pearson r=",
  sprintf(
    "%.3f",
    pearson_r
  ),
  "; Spearman rho=",
  sprintf(
    "%.3f",
    spearman_rho
  ),
  "), yielding a robust core of ",
  core_total_n,
  " genes (",
  core_up_n,
  " increased and ",
  core_down_n,
  " decreased in sepsis). ",
  
  "(D) Row-standardized variance-stabilized expression of 30 representative ",
  "genes drawn from the previously curated robust-core publication table, ",
  "comprising 15 genes increased and 15 genes decreased in sepsis. The upper ",
  "and lower heatmap blocks correspond to genes increased and decreased in ",
  "sepsis, respectively, and are separated by a horizontal gap. Samples were ",
  "clustered without use of clinical group labels; the annotation bar indicates ",
  "clinical condition. ",
  
  "(E) Hallmark gene-set enrichment analysis of the complete ranked gene list ",
  "derived from the targeted blood transcriptomic dataset. Positive normalized ",
  "enrichment scores indicate relative enrichment in sepsis blood and negative ",
  "scores indicate relative enrichment in control blood. The display contains ",
  "the highest-ranked positively enriched pathways and the available negatively ",
  "enriched pathways meeting FDR<0.25, selected from the frozen Hallmark GSEA ",
  "results without rerunning the analysis. Together, these findings demonstrate ",
  "a coordinated blood host-response state characterized by broad innate and ",
  "inflammatory activation accompanied by relative suppression of selected ",
  "adaptive, proliferative, and translational programs."
)


writeLines(
  caption_en,
  file.path(
    text_dir,
    "147_Figure1_caption_EN.txt"
  )
)


# =============================================================================
# 25. FINAL FIGURE CAPTION — RUSSIAN
# =============================================================================

caption_ru <- paste0(
  
  "Рисунок 1. Устойчивый транскриптомный ответ крови формирует ",
  "противоположно направленную миелоидно-адаптивную ось host response ",
  "при сепсисе. ",
  
  "(A) Анализ главных компонент variance-stabilized профилей таргетного ",
  "транскриптомного анализа крови пациентов с сепсисом и здоровых контролей. ",
  
  "(B) Картина дифференциальной экспрессии в основной модели DESeq2 ",
  "сепсис против контроля. DEG определяли как гены с adjusted p<0,05 и ",
  "|log2 fold change|>=1. Основная модель выявила ",
  primary_total_n,
  " DEG, из которых ",
  primary_up_n,
  " имели повышенную и ",
  primary_down_n,
  " сниженную экспрессию при сепсисе. Цветом выделены robust-core гены, ",
  "сохранявшие критерии DEG и направление эффекта в batch-adjusted модели. ",
  
  "(C) Согласованность gene-level log2 fold-change между основной и ",
  "batch-adjusted моделями DESeq2. Корреляции рассчитаны по всем генам с ",
  "конечными оценками эффекта; для наглядности оси ограничены центральным ",
  "диапазоном effect sizes. Диагональная пунктирная линия соответствует y=x. ",
  "Pearson r=",
  sprintf(
    "%.3f",
    pearson_r
  ),
  ", Spearman rho=",
  sprintf(
    "%.3f",
    spearman_rho
  ),
  ". Устойчивое ядро включало ",
  core_total_n,
  " генов: ",
  core_up_n,
  " с повышенной и ",
  core_down_n,
  " со сниженной экспрессией при сепсисе. ",
  
  "(D) Row-standardized VST-экспрессия 30 репрезентативных генов из ранее ",
  "сформированной curated robust-core publication table: 15 повышенных и ",
  "15 сниженных при сепсисе. Верхний и нижний блоки heatmap соответствуют ",
  "генам с повышенной и сниженной экспрессией соответственно и разделены ",
  "горизонтальным промежутком. Кластеризация образцов выполнялась без ",
  "использования клинической принадлежности. ",
  
  "(E) Hallmark GSEA полного ранжированного списка генов, полученного на основе ",
  "таргетного транскриптомного профилирования крови. Положительные NES ",
  "соответствуют относительному обогащению при сепсисе, отрицательные NES — ",
  "относительному обогащению в контрольной группе. Панель содержит наиболее ",
  "значимые положительно и отрицательно обогащенные Hallmark-программы с ",
  "FDR<0,25 из зафиксированных результатов GSEA; повторный анализ не ",
  "выполнялся. Совокупность результатов показывает координированный ",
  "host-response state с выраженной innate/inflammatory активацией и ",
  "относительным снижением отдельных adaptive, proliferative и translational ",
  "программ."
)


writeLines(
  caption_ru,
  file.path(
    text_dir,
    "147_Figure1_caption_RU.txt"
  )
)


writeLines(
  c(
    caption_en,
    "",
    caption_ru
  ),
  file.path(
    text_dir,
    "147_Figure1_caption_EN_RU.txt"
  )
)


# =============================================================================
# 26. PROVENANCE TABLE
# =============================================================================

provenance <- data.frame(
  
  Panel = c(
    "A",
    "B",
    "C",
    "D",
    "E"
  ),
  
  Content = c(
    "PCA",
    "Primary differential-expression volcano",
    "Primary vs batch-adjusted effect-size concordance",
    "Curated robust-core expression heatmap",
    "Hallmark GSEA"
  ),
  
  Primary_source = c(
    
    basename(
      vsd_file
    ),
    
    basename(
      table_s2
    ),
    
    basename(
      table_s2
    ),
    
    paste0(
      basename(
        table_s3
      ),
      " + ",
      basename(
        vsd_file
      )
    ),
    
    basename(
      table_s6
    )
  ),
  
  New_statistical_analysis = rep(
    "NO",
    5
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  provenance,
  file.path(
    table_dir,
    "147_Figure1_provenance.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 27. SAVE NUMERICAL AUDIT
# =============================================================================

write.csv(
  audit_table,
  file.path(
    table_dir,
    "147_Figure1_numerical_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 28. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    text_dir,
    "147_sessionInfo.txt"
  )
)


# =============================================================================
# 29. FINAL CONSOLE OUTPUT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 147 completed successfully.\n")
cat("FINAL Figure 1 visual-polish version generated.\n")
cat("====================================================================\n\n")


cat("DIFFERENTIAL EXPRESSION\n")
cat("-----------------------\n")


cat(
  "Primary DEG: ",
  primary_total_n,
  " (UP ",
  primary_up_n,
  "; DOWN ",
  primary_down_n,
  ")\n",
  sep = ""
)


cat(
  "Batch-adjusted DEG: ",
  batch_total_n,
  " (UP ",
  batch_up_n,
  "; DOWN ",
  batch_down_n,
  ")\n",
  sep = ""
)


cat(
  "Robust core: ",
  core_total_n,
  " (UP ",
  core_up_n,
  "; DOWN ",
  core_down_n,
  ")\n",
  sep = ""
)


cat(
  "Sex-linked genes in robust-core audit: ",
  ifelse(
    is.na(sex_linked_core_n),
    "NA",
    as.character(
      sex_linked_core_n
    )
  ),
  "\n",
  sep = ""
)


cat("\nMODEL CONCORDANCE\n")
cat("-----------------\n")


cat(
  "Pearson r = ",
  sprintf(
    "%.6f",
    pearson_r
  ),
  "\n",
  sep = ""
)


cat(
  "Spearman rho = ",
  sprintf(
    "%.6f",
    spearman_rho
  ),
  "\n",
  sep = ""
)


cat("\nPCA\n")
cat("---\n")


cat(
  "PC1 variance = ",
  percent_var[1],
  "%\n",
  sep = ""
)


cat(
  "PC2 variance = ",
  percent_var[2],
  "%\n",
  sep = ""
)


cat("\nVOLCANO\n")
cat("-------\n")


cat(
  "Gene labels displayed: ",
  length(
    label_genes
  ),
  "\n",
  sep = ""
)


cat(
  "Labels: ",
  paste(
    label_genes,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat("\nHEATMAP\n")
cat("-------\n")


cat(
  "Curated UP genes shown: ",
  n_heatmap_up,
  "\n",
  sep = ""
)


cat(
  "Curated DOWN genes shown: ",
  n_heatmap_down,
  "\n",
  sep = ""
)


cat(
  "Total heatmap genes: ",
  length(
    heatmap_genes
  ),
  "\n",
  sep = ""
)


cat(
  "Program legend: removed from final main figure\n"
)


cat(
  "UP/DOWN row gap: YES\n"
)


cat("\nHALLMARK GSEA\n")
cat("-------------\n")


cat(
  "Hallmark pathways tested: ",
  hallmark_total_n,
  "\n",
  sep = ""
)


cat(
  "Positive NES: ",
  hallmark_positive_n,
  "\n",
  sep = ""
)


cat(
  "Negative NES: ",
  hallmark_negative_n,
  "\n",
  sep = ""
)


cat(
  "Positive NES, FDR < 0.25: ",
  hallmark_positive_fdr025_n,
  "\n",
  sep = ""
)


cat(
  "Negative NES, FDR < 0.25: ",
  hallmark_negative_fdr025_n,
  "\n",
  sep = ""
)


cat(
  "Pathways displayed in Panel E: ",
  nrow(
    hallmark_plot_df
  ),
  "\n",
  sep = ""
)


myc_v1_row <- s6_all %>%
  dplyr::filter(
    ID ==
      "HALLMARK_MYC_TARGETS_V1"
  )


if (
  nrow(
    myc_v1_row
  ) == 1
) {
  
  cat(
    "MYC Targets V1: NES = ",
    sprintf(
      "%.6f",
      myc_v1_row$NES
    ),
    "; FDR = ",
    format(
      myc_v1_row$p.adjust,
      scientific = TRUE,
      digits = 5
    ),
    "\n",
    sep = ""
  )
}


cat("\nFINAL VISUAL STATUS\n")
cat("-------------------\n")

cat(
  "Global figure title inside image: NO\n"
)

cat(
  "Panel D Program legend: NO\n"
)

cat(
  "Panel D UP/DOWN gap: YES\n"
)

cat(
  "Panel E wording: complete ranked gene list\n"
)

cat(
  "Panel E legends: vertical legend box / separate rows\n"
)


cat("\nMAIN FIGURE\n")
cat("-----------\n")


cat(
  normalizePath(
    figure_png,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nSOURCE WORKBOOK\n")
cat("---------------\n")


cat(
  normalizePath(
    source_workbook,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nCAPTION\n")
cat("-------\n")


cat(
  normalizePath(
    file.path(
      text_dir,
      "147_Figure1_caption_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\n")
cat("If visual inspection passes, freeze as:\n")
cat("FIGURE 1 — FROZEN\n")
cat("Analytics frozen | source data frozen | composition frozen\n")
cat("\nDone.\n")