# ==============================================================================
# Script 137
# Blood biomarker benchmarking
#
# Project: Sepsis_DESeq2
#
# PURPOSE
# Benchmark the frozen five-gene blood host-response signature against
# predefined published transcriptomic signatures.
#
# CURRENT-STUDY SIGNATURES
#
# Primary five-gene:
#   UP   = CD177, HK3, IRAK3
#   DOWN = CARD11, IKZF2
#
# DCAF17 alternative:
#   UP   = CD177, HK3, IRAK3
#   DOWN = CARD11, DCAF17
#
# EXTERNAL / PUBLISHED CONTEXT
#
# SeptiCyte LAB-like:
#   PLAC8 + LAMP1 - PLA2G7 - CEACAM4
#
# SeptiCyte RAPID-like research contrast:
#   PLAC8 - PLA2G7
#
# FAIM3:PLAC8-like:
#   PLAC8 - FAIM3
#   FAIM3 may be represented as FCMR in current annotations
#
# Sepsis MetaScore-like:
#   UP:
#     CEACAM1, ZDHHC19, NMRK1/C9orf95,
#     GNA15, BATF, C3AR1
#
#   DOWN:
#     KIAA1370/FAM214A,
#     TGFBI, MTCH1, RPGRIP1, HLA-DPB1
#
# LIFTS-like:
#   -0.9305 * LRRN3
#   -0.9692 * IL2RB
#   -0.7378 * FCER1A
#   +0.8460 * TLR5
#   +0.8905 * S100A12
#
# IMPORTANT
#
# - Primary five-gene panel is frozen from Script 135.
# - No new feature selection is performed.
# - BP = sepsis blood RNA-seq.
# - BC = healthy-control blood RNA-seq.
# - Healthy controls are NOT noninfectious ICU/SIRS controls.
# - External scores are platform-adapted RNA-seq implementations.
# - SeptiCyte-derived contrasts are NOT official SeptiScores.
# - AUCs are apparent/contextual discrimination only.
# - No claim of independent diagnostic validation.
# - Blood only.
# - No urine.
# - No lncRNA.
# ==============================================================================


# ==============================================================================
# 0. SETTINGS
# ==============================================================================

options(stringsAsFactors = FALSE)

project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))

if (!dir.exists(project_dir)) {
  stop(
    paste0(
      "Project directory does not exist: ",
      project_dir
    )
  )
}

setwd(project_dir)

script_name <- "137_blood_biomarker_benchmarking.R"
run_date <- Sys.time()

cat("\n")
cat("====================================================================\n")
cat("Running Script 137\n")
cat("Blood biomarker benchmarking\n")
cat("====================================================================\n\n")

cat("Project directory:\n")
cat(
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)

cat("Run date:\n")
cat(
  as.character(run_date),
  "\n\n"
)


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_packages <- c(
  "edgeR",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "pROC",
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
      "Missing required packages:\n",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pROC)
  library(openxlsx)
})

cat("Required packages loaded successfully.\n\n")


# ==============================================================================
# 2. GENERAL HELPERS
# ==============================================================================

clean_sample_id <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  x <- gsub(
    "[^A-Z0-9]",
    "",
    x
  )
  
  return(x)
}


clean_gene_symbol <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  return(x)
}


z_safe <- function(x) {
  
  x <- as.numeric(x)
  
  sd_x <- stats::sd(
    x,
    na.rm = TRUE
  )
  
  if (
    !is.finite(sd_x) ||
    sd_x == 0
  ) {
    
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }
  
  return(
    as.numeric(
      scale(x)
    )
  )
}


safe_htest_statistic <- function(x) {
  
  if (
    is.list(x) &&
    "statistic" %in% names(x) &&
    length(x[["statistic"]]) >= 1
  ) {
    
    return(
      as.numeric(
        x[["statistic"]][1]
      )
    )
  }
  
  return(NA_real_)
}


safe_htest_pvalue <- function(x) {
  
  if (
    is.list(x) &&
    "p.value" %in% names(x) &&
    length(x[["p.value"]]) >= 1
  ) {
    
    return(
      as.numeric(
        x[["p.value"]][1]
      )
    )
  }
  
  return(NA_real_)
}


safe_htest_estimate <- function(x) {
  
  if (
    is.list(x) &&
    "estimate" %in% names(x) &&
    length(x[["estimate"]]) >= 1
  ) {
    
    return(
      as.numeric(
        x[["estimate"]][1]
      )
    )
  }
  
  return(NA_real_)
}


epsilon_squared_kw <- function(
    x,
    group
) {
  
  keep <- complete.cases(
    x,
    group
  )
  
  x2 <- x[keep]
  
  g2 <- droplevels(
    factor(
      group[keep]
    )
  )
  
  if (
    length(x2) < 6 ||
    nlevels(g2) < 2
  ) {
    
    return(NA_real_)
  }
  
  ht <- stats::kruskal.test(
    x2 ~ g2
  )
  
  H <- safe_htest_statistic(ht)
  
  n_total <- length(x2)
  k <- nlevels(g2)
  
  epsilon2 <- (
    H - k + 1
  ) / (
    n_total - k
  )
  
  return(
    max(
      0,
      epsilon2
    )
  )
}


# ==============================================================================
# 3. INPUT FILES
# ==============================================================================

counts_file <- file.path(
  "data",
  "counts_all.csv"
)


script135_scores_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "135_validation",
  "tables",
  "135_blood_scores_with_final_endotypes.csv"
)


demographics_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity",
  "tables",
  "136b_blood_demographics.csv"
)


required_files <- c(
  counts_file,
  script135_scores_file,
  demographics_file
)


if (any(!file.exists(required_files))) {
  
  missing_files <- required_files[
    !file.exists(
      required_files
    )
  ]
  
  stop(
    paste0(
      "Missing required input files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


cat("Input files:\n")

cat(
  "Raw counts: ",
  counts_file,
  "\n",
  sep = ""
)

cat(
  "Script 135 scores: ",
  script135_scores_file,
  "\n",
  sep = ""
)

cat(
  "Demographics: ",
  demographics_file,
  "\n\n",
  sep = ""
)


# ==============================================================================
# 4. OUTPUT DIRECTORIES
# ==============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "137_benchmarking"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

figures_dir <- file.path(
  output_dir,
  "figures"
)

text_dir <- file.path(
  output_dir,
  "text"
)

logs_dir <- file.path(
  output_dir,
  "logs"
)


dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  text_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("Output folder:\n")

cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


# ==============================================================================
# 5. SIGNATURE REFERENCE REGISTRY
# ==============================================================================

reference_table <- tibble(
  
  signature = c(
    "Primary five-gene",
    "DCAF17 alternative",
    "SeptiCyte LAB-like",
    "SeptiCyte RAPID-like",
    "FAIM3:PLAC8-like",
    "Sepsis MetaScore-like",
    "LIFTS-like"
  ),
  
  category = c(
    "Current study",
    "Current study sensitivity",
    "Published comparator",
    "Published comparator",
    "Published comparator",
    "Published comparator",
    "Published comparator"
  ),
  
  implementation = c(
    
    paste0(
      "Frozen Script 135 score: mean z(CD177, HK3, IRAK3) ",
      "- mean z(CARD11, IKZF2)"
    ),
    
    paste0(
      "Sensitivity score: mean z(CD177, HK3, IRAK3) ",
      "- mean z(CARD11, DCAF17)"
    ),
    
    "PLAC8 + LAMP1 - PLA2G7 - CEACAM4 on TMM logCPM",
    
    paste0(
      "PLAC8 - PLA2G7 research contrast on TMM logCPM; ",
      "not an official SeptiScore"
    ),
    
    paste0(
      "PLAC8 - FAIM3/FCMR infection-oriented contrast ",
      "on TMM logCPM"
    ),
    
    paste0(
      "Mean logCPM of six infection-up genes minus mean logCPM ",
      "of five infection-down genes"
    ),
    
    paste0(
      "-0.9305*LRRN3 -0.9692*IL2RB -0.7378*FCER1A ",
      "+0.8460*TLR5 +0.8905*S100A12"
    )
  ),
  
  interpretation_note = c(
    "Internal candidate signature selected in current cohort",
    "Internal sensitivity panel selected in current cohort",
    "Platform-adapted contextual reconstruction",
    "Two-gene research contrast only",
    "Platform-adapted contextual reconstruction",
    "Platform-adapted contextual reconstruction",
    "Published coefficients applied to current RNA-seq scale"
  )
)


# ==============================================================================
# 6. GENE ALIASES
# ==============================================================================

gene_aliases <- list(
  
  CD177 = c(
    "CD177"
  ),
  
  HK3 = c(
    "HK3"
  ),
  
  IRAK3 = c(
    "IRAK3"
  ),
  
  CARD11 = c(
    "CARD11"
  ),
  
  IKZF2 = c(
    "IKZF2"
  ),
  
  DCAF17 = c(
    "DCAF17"
  ),
  
  PLAC8 = c(
    "PLAC8"
  ),
  
  LAMP1 = c(
    "LAMP1"
  ),
  
  PLA2G7 = c(
    "PLA2G7"
  ),
  
  CEACAM4 = c(
    "CEACAM4"
  ),
  
  FAIM3 = c(
    "FAIM3",
    "FCMR"
  ),
  
  CEACAM1 = c(
    "CEACAM1"
  ),
  
  ZDHHC19 = c(
    "ZDHHC19"
  ),
  
  NMRK1 = c(
    "NMRK1",
    "C9ORF95"
  ),
  
  GNA15 = c(
    "GNA15"
  ),
  
  BATF = c(
    "BATF"
  ),
  
  C3AR1 = c(
    "C3AR1"
  ),
  
  KIAA1370 = c(
    "KIAA1370",
    "FAM214A"
  ),
  
  TGFBI = c(
    "TGFBI"
  ),
  
  MTCH1 = c(
    "MTCH1"
  ),
  
  RPGRIP1 = c(
    "RPGRIP1"
  ),
  
  `HLA-DPB1` = c(
    "HLA-DPB1",
    "HLADPB1"
  ),
  
  LRRN3 = c(
    "LRRN3"
  ),
  
  IL2RB = c(
    "IL2RB"
  ),
  
  FCER1A = c(
    "FCER1A"
  ),
  
  TLR5 = c(
    "TLR5"
  ),
  
  S100A12 = c(
    "S100A12"
  )
)


# ==============================================================================
# 7. EXPECTED DIRECTIONS
# ==============================================================================

expected_direction <- c(
  
  CD177 = "UP",
  HK3 = "UP",
  IRAK3 = "UP",
  CARD11 = "DOWN",
  IKZF2 = "DOWN",
  DCAF17 = "DOWN",
  
  PLAC8 = "UP",
  LAMP1 = "UP",
  PLA2G7 = "DOWN",
  CEACAM4 = "DOWN",
  
  FAIM3 = "DOWN",
  
  CEACAM1 = "UP",
  ZDHHC19 = "UP",
  NMRK1 = "UP",
  GNA15 = "UP",
  BATF = "UP",
  C3AR1 = "UP",
  
  KIAA1370 = "DOWN",
  TGFBI = "DOWN",
  MTCH1 = "DOWN",
  RPGRIP1 = "DOWN",
  `HLA-DPB1` = "DOWN",
  
  LRRN3 = "DOWN",
  IL2RB = "DOWN",
  FCER1A = "DOWN",
  TLR5 = "UP",
  S100A12 = "UP"
)


# ==============================================================================
# 8. SIGNATURE GENE SETS
# ==============================================================================

signature_gene_sets <- list(
  
  Primary_5gene = c(
    "CD177",
    "HK3",
    "IRAK3",
    "CARD11",
    "IKZF2"
  ),
  
  DCAF17_alternative = c(
    "CD177",
    "HK3",
    "IRAK3",
    "CARD11",
    "DCAF17"
  ),
  
  SeptiCyte_LAB_like = c(
    "PLAC8",
    "LAMP1",
    "PLA2G7",
    "CEACAM4"
  ),
  
  SeptiCyte_RAPID_like = c(
    "PLAC8",
    "PLA2G7"
  ),
  
  FAIM3_PLAC8_like = c(
    "FAIM3",
    "PLAC8"
  ),
  
  Sepsis_MetaScore_like = c(
    "CEACAM1",
    "ZDHHC19",
    "NMRK1",
    "GNA15",
    "BATF",
    "C3AR1",
    "KIAA1370",
    "TGFBI",
    "MTCH1",
    "RPGRIP1",
    "HLA-DPB1"
  ),
  
  LIFTS_like = c(
    "LRRN3",
    "IL2RB",
    "FCER1A",
    "TLR5",
    "S100A12"
  )
)


# ==============================================================================
# 9. READ SCRIPT 135 SCORES
# ==============================================================================

blood_scores <- read.csv(
  script135_scores_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (!"sample_id" %in% names(blood_scores)) {
  stop(
    "sample_id column missing from Script 135 score table."
  )
}


blood_scores$sample_id <- clean_sample_id(
  blood_scores$sample_id
)


blood_scores <- blood_scores %>%
  filter(
    grepl(
      "^(BP|BC)[0-9]+$",
      sample_id
    )
  )


if (nrow(blood_scores) != 45) {
  
  stop(
    paste0(
      "Expected 45 blood samples in Script 135 score table; observed ",
      nrow(blood_scores)
    )
  )
}


required_score_columns <- c(
  "sample_id",
  "primary_5gene_score",
  "SRS",
  "SRSq",
  "CTS"
)


missing_score_columns <- setdiff(
  required_score_columns,
  names(blood_scores)
)


if (length(missing_score_columns) > 0) {
  
  stop(
    paste0(
      "Missing Script 135 score columns: ",
      paste(
        missing_score_columns,
        collapse = ", "
      )
    )
  )
}


blood_scores$condition <- dplyr::case_when(
  
  grepl(
    "^BP[0-9]+$",
    blood_scores$sample_id
  ) ~ "BP",
  
  grepl(
    "^BC[0-9]+$",
    blood_scores$sample_id
  ) ~ "BC",
  
  TRUE ~ NA_character_
)


blood_scores$condition <- factor(
  blood_scores$condition,
  levels = c(
    "BC",
    "BP"
  )
)


if (any(is.na(blood_scores$condition))) {
  stop(
    "Unable to determine BP/BC status from sample_id."
  )
}


cat("Blood score table:\n")

print(
  table(
    blood_scores$condition
  )
)

cat("\n")


# ==============================================================================
# 10. READ DEMOGRAPHICS
# ==============================================================================

demographics <- read.csv(
  demographics_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_demo_columns <- c(
  "sample_id",
  "age_years",
  "sex"
)


missing_demo_columns <- setdiff(
  required_demo_columns,
  names(demographics)
)


if (length(missing_demo_columns) > 0) {
  
  stop(
    paste0(
      "Missing demographic columns: ",
      paste(
        missing_demo_columns,
        collapse = ", "
      )
    )
  )
}


demographics$sample_id <- clean_sample_id(
  demographics$sample_id
)


# ==============================================================================
# 11. READ RAW COUNTS
# ==============================================================================

cat("Reading raw count matrix:\n")
cat(
  counts_file,
  "\n\n"
)


counts_raw <- read.csv(
  counts_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Raw count dimensions: ",
  nrow(counts_raw),
  " rows x ",
  ncol(counts_raw),
  " columns\n\n",
  sep = ""
)


# ==============================================================================
# 12. DETECT BLOOD SAMPLE COLUMNS
# ==============================================================================

raw_column_names <- names(
  counts_raw
)


clean_column_names <- clean_sample_id(
  raw_column_names
)


blood_column_index <- which(
  grepl(
    "^(BP|BC)[0-9]+$",
    clean_column_names
  )
)


if (length(blood_column_index) != 45) {
  
  stop(
    paste0(
      "Expected exactly 45 BP/BC columns in counts_all.csv; detected ",
      length(blood_column_index)
    )
  )
}


blood_sample_names <- clean_column_names[
  blood_column_index
]


if (anyDuplicated(blood_sample_names) > 0) {
  stop(
    "Duplicated blood sample IDs detected in count matrix."
  )
}


cat("Blood columns detected:\n")

print(
  table(
    ifelse(
      grepl(
        "^BP",
        blood_sample_names
      ),
      "BP",
      "BC"
    )
  )
)

cat("\n")


# ==============================================================================
# 13. DETECT GENE IDENTIFIER COLUMN
# ==============================================================================

non_sample_columns <- setdiff(
  seq_along(raw_column_names),
  blood_column_index
)


candidate_gene_names <- c(
  "gene",
  "gene_symbol",
  "genesymbol",
  "symbol",
  "target",
  "target_id",
  "gene_id",
  "geneid"
)


normalized_raw_names <- tolower(
  gsub(
    "[^a-z0-9]",
    "",
    raw_column_names
  )
)


normalized_gene_candidates <- tolower(
  gsub(
    "[^a-z0-9]",
    "",
    candidate_gene_names
  )
)


candidate_hit <- which(
  normalized_raw_names %in%
    normalized_gene_candidates
)


candidate_hit <- intersect(
  candidate_hit,
  non_sample_columns
)


if (length(candidate_hit) >= 1) {
  
  gene_column_index <- candidate_hit[1]
  
} else {
  
  gene_column_index <- non_sample_columns[1]
}


gene_column_name <- raw_column_names[
  gene_column_index
]


cat("Detected gene identifier column:\n")
cat(
  gene_column_name,
  "\n\n"
)


# ==============================================================================
# 14. BUILD BLOOD COUNT MATRIX
# ==============================================================================

gene_symbols <- clean_gene_symbol(
  counts_raw[[gene_column_index]]
)


counts_df <- counts_raw[
  ,
  blood_column_index,
  drop = FALSE
]


counts_matrix <- do.call(
  cbind,
  lapply(
    counts_df,
    function(x) {
      suppressWarnings(
        as.numeric(x)
      )
    }
  )
)


counts_matrix <- as.matrix(
  counts_matrix
)


colnames(
  counts_matrix
) <- blood_sample_names


if (any(is.na(counts_matrix))) {
  
  stop(
    paste0(
      "NA values detected after conversion of blood counts to numeric. ",
      "Review counts_all.csv."
    )
  )
}


valid_gene <- !is.na(
  gene_symbols
) &
  gene_symbols != ""


gene_symbols <- gene_symbols[
  valid_gene
]


counts_matrix <- counts_matrix[
  valid_gene,
  ,
  drop = FALSE
]


# ==============================================================================
# 15. REMOVE ERCC SPIKE-INS
# ==============================================================================

ercc_flag <- grepl(
  "^ERCC[-_]",
  gene_symbols
)


cat(
  "ERCC spike-ins removed: ",
  sum(ercc_flag),
  "\n",
  sep = ""
)


gene_symbols <- gene_symbols[
  !ercc_flag
]


counts_matrix <- counts_matrix[
  !ercc_flag,
  ,
  drop = FALSE
]


# ==============================================================================
# 16. COLLAPSE DUPLICATE SYMBOLS
# ==============================================================================

counts_matrix <- rowsum(
  counts_matrix,
  group = gene_symbols,
  reorder = FALSE
)


if (
  any(
    counts_matrix < 0,
    na.rm = TRUE
  )
) {
  
  stop(
    "Negative raw counts detected."
  )
}


cat(
  "Genes after ERCC removal and duplicate collapsing: ",
  nrow(counts_matrix),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 17. TMM NORMALIZATION AND logCPM
# ==============================================================================

dge <- edgeR::DGEList(
  counts = counts_matrix
)


if (
  "normLibSizes" %in%
  getNamespaceExports(
    "edgeR"
  )
) {
  
  dge <- edgeR::normLibSizes(
    dge,
    method = "TMM"
  )
  
} else {
  
  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )
}


logcpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 2
)


rownames(
  logcpm
) <- clean_gene_symbol(
  rownames(logcpm)
)


# ==============================================================================
# 18. ALIGN SAMPLE ORDER
# ==============================================================================

sample_order <- blood_scores$sample_id


missing_samples <- setdiff(
  sample_order,
  colnames(logcpm)
)


if (length(missing_samples) > 0) {
  
  stop(
    paste0(
      "Samples missing from normalized blood count matrix:\n",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}


logcpm <- logcpm[
  ,
  sample_order,
  drop = FALSE
]


if (!identical(
  colnames(logcpm),
  blood_scores$sample_id
)) {
  
  stop(
    "Sample alignment between scores and logCPM failed."
  )
}


# ==============================================================================
# 19. GENE RESOLUTION HELPERS
# ==============================================================================

available_gene_symbols <- rownames(
  logcpm
)


resolve_gene <- function(canonical_gene) {
  
  if (!canonical_gene %in% names(gene_aliases)) {
    return(NA_character_)
  }
  
  candidates <- clean_gene_symbol(
    gene_aliases[[canonical_gene]]
  )
  
  hit <- candidates[
    candidates %in%
      available_gene_symbols
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  return(
    hit[1]
  )
}


get_gene_expression <- function(canonical_gene) {
  
  resolved <- resolve_gene(
    canonical_gene
  )
  
  if (is.na(resolved)) {
    
    return(
      rep(
        NA_real_,
        ncol(logcpm)
      )
    )
  }
  
  return(
    as.numeric(
      logcpm[
        resolved,
        ,
        drop = TRUE
      ]
    )
  )
}


signature_complete <- function(genes) {
  
  resolved <- vapply(
    genes,
    resolve_gene,
    character(1)
  )
  
  return(
    all(
      !is.na(resolved)
    )
  )
}


# ==============================================================================
# 20. SIGNATURE GENE COMPLETENESS
# ==============================================================================

signature_status_list <- list()

status_counter <- 1L


for (signature_name in names(signature_gene_sets)) {
  
  genes <- signature_gene_sets[[signature_name]]
  
  resolved <- vapply(
    genes,
    resolve_gene,
    character(1)
  )
  
  signature_status_list[[status_counter]] <- tibble(
    
    signature =
      signature_name,
    
    n_required_genes =
      length(genes),
    
    n_detected_genes =
      sum(
        !is.na(resolved)
      ),
    
    complete =
      all(
        !is.na(resolved)
      ),
    
    required_genes =
      paste(
        genes,
        collapse = "; "
      ),
    
    resolved_symbols =
      paste(
        ifelse(
          is.na(resolved),
          "MISSING",
          resolved
        ),
        collapse = "; "
      )
  )
  
  status_counter <- status_counter + 1L
}


signature_status <- bind_rows(
  signature_status_list
)


cat("Signature gene-set completeness:\n")

print(
  signature_status,
  n = Inf
)

cat("\n")


# ==============================================================================
# 21. GENE DIRECTION AUDIT
# ==============================================================================

all_benchmark_genes <- unique(
  unlist(
    signature_gene_sets,
    use.names = FALSE
  )
)


gene_audit_list <- list()

audit_counter <- 1L


for (gene in all_benchmark_genes) {
  
  resolved <- resolve_gene(
    gene
  )
  
  memberships <- names(
    signature_gene_sets
  )[
    vapply(
      signature_gene_sets,
      function(x) {
        gene %in% x
      },
      logical(1)
    )
  ]
  
  
  expected <- unname(
    expected_direction[
      gene
    ]
  )
  
  
  if (is.na(resolved)) {
    
    gene_audit_list[[audit_counter]] <- tibble(
      
      canonical_gene =
        gene,
      
      resolved_symbol =
        NA_character_,
      
      present =
        FALSE,
      
      signatures =
        paste(
          memberships,
          collapse = "; "
        ),
      
      expected_direction =
        expected,
      
      median_BC =
        NA_real_,
      
      median_BP =
        NA_real_,
      
      BP_minus_BC =
        NA_real_,
      
      observed_direction =
        NA_character_,
      
      direction_concordant =
        NA,
      
      wilcoxon_p =
        NA_real_
    )
    
    audit_counter <- audit_counter + 1L
    next
  }
  
  
  expr <- get_gene_expression(
    gene
  )
  
  
  bc_values <- expr[
    blood_scores$condition == "BC"
  ]
  
  
  bp_values <- expr[
    blood_scores$condition == "BP"
  ]
  
  
  median_bc <- median(
    bc_values,
    na.rm = TRUE
  )
  
  
  median_bp <- median(
    bp_values,
    na.rm = TRUE
  )
  
  
  difference <- median_bp -
    median_bc
  
  
  observed <- ifelse(
    difference >= 0,
    "UP",
    "DOWN"
  )
  
  
  ht <- stats::wilcox.test(
    x = bp_values,
    y = bc_values,
    exact = FALSE,
    paired = FALSE
  )
  
  
  gene_audit_list[[audit_counter]] <- tibble(
    
    canonical_gene =
      gene,
    
    resolved_symbol =
      resolved,
    
    present =
      TRUE,
    
    signatures =
      paste(
        memberships,
        collapse = "; "
      ),
    
    expected_direction =
      expected,
    
    median_BC =
      median_bc,
    
    median_BP =
      median_bp,
    
    BP_minus_BC =
      difference,
    
    observed_direction =
      observed,
    
    direction_concordant =
      observed == expected,
    
    wilcoxon_p =
      safe_htest_pvalue(
        ht
      )
  )
  
  
  audit_counter <- audit_counter + 1L
}


gene_audit <- bind_rows(
  gene_audit_list
)


gene_audit$BH_gene_audit <- NA_real_


valid_gene_p <- which(
  !is.na(
    gene_audit$wilcoxon_p
  ) &
    is.finite(
      gene_audit$wilcoxon_p
    )
)


if (length(valid_gene_p) > 0) {
  
  gene_audit$BH_gene_audit[
    valid_gene_p
  ] <- stats::p.adjust(
    gene_audit$wilcoxon_p[
      valid_gene_p
    ],
    method = "BH"
  )
}


cat("Benchmark gene availability and direction audit:\n")

print(
  gene_audit %>%
    select(
      canonical_gene,
      resolved_symbol,
      present,
      expected_direction,
      observed_direction,
      direction_concordant,
      wilcoxon_p,
      BH_gene_audit
    ) %>%
    tibble::as_tibble(),
  n = Inf
)

cat("\n")


# ==============================================================================
# 22. BUILD SAMPLE SCORE TABLE
# ==============================================================================

score_table <- blood_scores %>%
  select(
    sample_id,
    condition,
    primary_5gene_score,
    SRS,
    SRSq,
    CTS
  ) %>%
  left_join(
    demographics %>%
      select(
        sample_id,
        age_years,
        sex
      ),
    by = "sample_id"
  )


if (nrow(score_table) != 45) {
  
  stop(
    paste0(
      "Expected 45 rows in score table; observed ",
      nrow(score_table)
    )
  )
}


# ==============================================================================
# 23. DCAF17 ALTERNATIVE SCORE
# ==============================================================================

if (
  signature_complete(
    signature_gene_sets[["DCAF17_alternative"]]
  )
) {
  
  z_CD177 <- z_safe(
    get_gene_expression("CD177")
  )
  
  z_HK3 <- z_safe(
    get_gene_expression("HK3")
  )
  
  z_IRAK3 <- z_safe(
    get_gene_expression("IRAK3")
  )
  
  z_CARD11 <- z_safe(
    get_gene_expression("CARD11")
  )
  
  z_DCAF17 <- z_safe(
    get_gene_expression("DCAF17")
  )
  
  
  score_table$DCAF17_alt_score <-
    rowMeans(
      cbind(
        z_CD177,
        z_HK3,
        z_IRAK3
      ),
      na.rm = FALSE
    ) -
    rowMeans(
      cbind(
        z_CARD11,
        z_DCAF17
      ),
      na.rm = FALSE
    )
  
} else {
  
  score_table$DCAF17_alt_score <- NA_real_
}


# ==============================================================================
# 24. SEPTICYTE LAB-LIKE SCORE
# ==============================================================================

if (
  signature_complete(
    signature_gene_sets[["SeptiCyte_LAB_like"]]
  )
) {
  
  score_table$SeptiCyte_LAB_like <-
    get_gene_expression("PLAC8") +
    get_gene_expression("LAMP1") -
    get_gene_expression("PLA2G7") -
    get_gene_expression("CEACAM4")
  
} else {
  
  score_table$SeptiCyte_LAB_like <- NA_real_
}


# ==============================================================================
# 25. SEPTICYTE RAPID-LIKE TWO-GENE CONTRAST
# ==============================================================================

if (
  signature_complete(
    signature_gene_sets[["SeptiCyte_RAPID_like"]]
  )
) {
  
  score_table$SeptiCyte_RAPID_like <-
    get_gene_expression("PLAC8") -
    get_gene_expression("PLA2G7")
  
} else {
  
  score_table$SeptiCyte_RAPID_like <- NA_real_
}


# ==============================================================================
# 26. FAIM3:PLAC8-LIKE SCORE
# ==============================================================================

if (
  signature_complete(
    signature_gene_sets[["FAIM3_PLAC8_like"]]
  )
) {
  
  score_table$FAIM3_PLAC8_like <-
    get_gene_expression("PLAC8") -
    get_gene_expression("FAIM3")
  
} else {
  
  score_table$FAIM3_PLAC8_like <- NA_real_
}


# ==============================================================================
# 27. SEPSIS METASCORE-LIKE
# ==============================================================================

sms_up <- c(
  "CEACAM1",
  "ZDHHC19",
  "NMRK1",
  "GNA15",
  "BATF",
  "C3AR1"
)


sms_down <- c(
  "KIAA1370",
  "TGFBI",
  "MTCH1",
  "RPGRIP1",
  "HLA-DPB1"
)


if (
  signature_complete(
    c(
      sms_up,
      sms_down
    )
  )
) {
  
  sms_up_matrix <- sapply(
    sms_up,
    get_gene_expression
  )
  
  
  sms_down_matrix <- sapply(
    sms_down,
    get_gene_expression
  )
  
  
  score_table$Sepsis_MetaScore_like <-
    rowMeans(
      sms_up_matrix,
      na.rm = FALSE
    ) -
    rowMeans(
      sms_down_matrix,
      na.rm = FALSE
    )
  
} else {
  
  score_table$Sepsis_MetaScore_like <- NA_real_
}


# ==============================================================================
# 28. LIFTS-LIKE SCORE
# ==============================================================================

if (
  signature_complete(
    signature_gene_sets[["LIFTS_like"]]
  )
) {
  
  score_table$LIFTS_like <-
    -0.9305 * get_gene_expression("LRRN3") +
    -0.9692 * get_gene_expression("IL2RB") +
    -0.7378 * get_gene_expression("FCER1A") +
    0.8460 * get_gene_expression("TLR5") +
    0.8905 * get_gene_expression("S100A12")
  
} else {
  
  score_table$LIFTS_like <- NA_real_
}


# ==============================================================================
# 29. SCORE REGISTRY
# ==============================================================================

score_registry <- tibble(
  
  score_column = c(
    "primary_5gene_score",
    "DCAF17_alt_score",
    "SeptiCyte_LAB_like",
    "SeptiCyte_RAPID_like",
    "FAIM3_PLAC8_like",
    "Sepsis_MetaScore_like",
    "LIFTS_like"
  ),
  
  display_name = c(
    "Primary 5-gene",
    "DCAF17 alternative",
    "SeptiCyte LAB-like",
    "PLAC8-PLA2G7 contrast",
    "FAIM3:PLAC8-like",
    "Sepsis MetaScore-like",
    "LIFTS-like"
  ),
  
  panel_origin = c(
    "Current study",
    "Current study",
    "Published",
    "Published",
    "Published",
    "Published",
    "Published"
  ),
  
  same_cohort_selected = c(
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  )
)


score_registry$available <- vapply(
  
  score_registry$score_column,
  
  function(score_name) {
    
    score_name %in% names(score_table) &&
      any(
        !is.na(
          score_table[[score_name]]
        )
      )
  },
  
  logical(1)
)


cat("Signature score availability:\n")

print(
  score_registry,
  n = Inf
)

cat("\n")


available_scores <- score_registry$score_column[
  score_registry$available
]


if (length(available_scores) < 2) {
  
  stop(
    "Fewer than two biomarker scores available for benchmarking."
  )
}


# ==============================================================================
# 30. BP vs BC APPARENT AUC BENCHMARK
# ==============================================================================

auc_results_list <- list()

roc_objects <- list()

auc_counter <- 1L


for (score_name in available_scores) {
  
  x <- score_table[[score_name]]
  
  
  keep <- !is.na(x) &
    !is.na(
      score_table$condition
    )
  
  
  x2 <- x[keep]
  
  condition2 <- droplevels(
    score_table$condition[
      keep
    ]
  )
  
  
  if (
    nlevels(condition2) != 2 ||
    length(unique(x2)) < 2
  ) {
    
    next
  }
  
  
  roc_object <- pROC::roc(
    
    response =
      condition2,
    
    predictor =
      x2,
    
    levels =
      c(
        "BC",
        "BP"
      ),
    
    direction =
      "<",
    
    quiet =
      TRUE
  )
  
  
  auc_value <- as.numeric(
    pROC::auc(
      roc_object
    )
  )
  
  
  auc_ci <- tryCatch(
    
    as.numeric(
      pROC::ci.auc(
        roc_object,
        method = "delong"
      )
    ),
    
    error = function(e) {
      
      c(
        NA_real_,
        NA_real_,
        NA_real_
      )
    }
  )
  
  
  bp_values <- x[
    score_table$condition == "BP"
  ]
  
  
  bc_values <- x[
    score_table$condition == "BC"
  ]
  
  
  ht <- stats::wilcox.test(
    x = bp_values,
    y = bc_values,
    exact = FALSE,
    paired = FALSE
  )
  
  
  auc_results_list[[auc_counter]] <- tibble(
    
    score_column =
      score_name,
    
    n_total =
      sum(keep),
    
    n_BC =
      sum(
        condition2 == "BC"
      ),
    
    n_BP =
      sum(
        condition2 == "BP"
      ),
    
    median_BC =
      median(
        bc_values,
        na.rm = TRUE
      ),
    
    median_BP =
      median(
        bp_values,
        na.rm = TRUE
      ),
    
    median_difference_BP_minus_BC =
      median(
        bp_values,
        na.rm = TRUE
      ) -
      median(
        bc_values,
        na.rm = TRUE
      ),
    
    Wilcoxon_p =
      safe_htest_pvalue(
        ht
      ),
    
    AUC_fixed_direction =
      auc_value,
    
    AUC_CI_low =
      auc_ci[1],
    
    AUC_CI_mid =
      auc_ci[2],
    
    AUC_CI_high =
      auc_ci[3],
    
    interpretation =
      paste0(
        "Contextual BP-vs-BC discrimination only; ",
        "not independent clinical validation"
      )
  )
  
  
  roc_objects[[score_name]] <- roc_object
  
  auc_counter <- auc_counter + 1L
}


auc_results <- bind_rows(
  auc_results_list
) %>%
  left_join(
    score_registry,
    by = "score_column"
  ) %>%
  arrange(
    desc(
      AUC_fixed_direction
    )
  )


if (nrow(auc_results) == 0) {
  
  stop(
    "No valid AUC results were generated."
  )
}


auc_results$BH_Wilcoxon <- stats::p.adjust(
  auc_results$Wilcoxon_p,
  method = "BH"
)


cat("BP vs BC benchmark performance:\n")

print(
  auc_results %>%
    select(
      display_name,
      n_total,
      median_BC,
      median_BP,
      AUC_fixed_direction,
      AUC_CI_low,
      AUC_CI_high,
      Wilcoxon_p,
      BH_Wilcoxon
    ) %>%
    tibble::as_tibble(),
  n = Inf
)

cat("\n")


# ==============================================================================
# 31. AGE/SEX-ADJUSTED BP vs BC EFFECTS
# ==============================================================================

adjusted_results_list <- list()

adjusted_counter <- 1L


for (score_name in available_scores) {
  
  tmp <- tibble(
    
    score =
      score_table[[score_name]],
    
    condition_binary =
      ifelse(
        score_table$condition == "BP",
        1,
        0
      ),
    
    age_years =
      suppressWarnings(
        as.numeric(
          score_table$age_years
        )
      ),
    
    sex =
      factor(
        score_table$sex
      )
  )
  
  
  tmp <- tmp %>%
    filter(
      complete.cases(
        score,
        condition_binary,
        age_years,
        sex
      )
    )
  
  
  if (
    nrow(tmp) < 20 ||
    !is.finite(
      stats::sd(
        tmp$score
      )
    ) ||
    stats::sd(
      tmp$score
    ) == 0
  ) {
    
    next
  }
  
  
  tmp$score_z <- z_safe(
    tmp$score
  )
  
  
  fit <- stats::lm(
    score_z ~
      condition_binary +
      age_years +
      sex,
    data = tmp
  )
  
  
  coefficients <- summary(
    fit
  )$coefficients
  
  
  if (
    !"condition_binary" %in%
    rownames(coefficients)
  ) {
    
    next
  }
  
  
  condition_coef <- coefficients[
    "condition_binary",
    ,
    drop = FALSE
  ]
  
  
  adjusted_results_list[[adjusted_counter]] <- tibble(
    
    score_column =
      score_name,
    
    n =
      nrow(tmp),
    
    adjusted_standardized_BP_vs_BC_effect =
      condition_coef[
        1,
        "Estimate"
      ],
    
    standard_error =
      condition_coef[
        1,
        "Std. Error"
      ],
    
    t_value =
      condition_coef[
        1,
        "t value"
      ],
    
    p_value =
      condition_coef[
        1,
        "Pr(>|t|)"
      ],
    
    adjusted_R2 =
      summary(
        fit
      )$adj.r.squared
  )
  
  
  adjusted_counter <- adjusted_counter + 1L
}


adjusted_results <- bind_rows(
  adjusted_results_list
) %>%
  left_join(
    score_registry,
    by = "score_column"
  )


if (nrow(adjusted_results) > 0) {
  
  adjusted_results$BH_adjusted <- stats::p.adjust(
    adjusted_results$p_value,
    method = "BH"
  )
  
  
  adjusted_results <- adjusted_results %>%
    arrange(
      p_value
    )
}


cat("Age/sex-adjusted BP vs BC score effects:\n")

if (nrow(adjusted_results) > 0) {
  
  print(
    adjusted_results %>%
      select(
        display_name,
        n,
        adjusted_standardized_BP_vs_BC_effect,
        p_value,
        BH_adjusted,
        adjusted_R2
      ) %>%
      tibble::as_tibble(),
    n = Inf
  )
  
} else {
  
  cat("No adjusted score results generated.\n")
}

cat("\n")


# ==============================================================================
# 32. BP-ONLY DATA FOR ENDOTYPE ANALYSES
# ==============================================================================

bp_scores <- score_table %>%
  filter(
    condition == "BP"
  )


if (nrow(bp_scores) != 35) {
  
  stop(
    paste0(
      "Expected 35 BP samples; observed ",
      nrow(bp_scores)
    )
  )
}


cat("BP SRS distribution:\n")

print(
  table(
    bp_scores$SRS,
    useNA = "ifany"
  )
)

cat("\n")


cat("BP CTS distribution:\n")

print(
  table(
    bp_scores$CTS,
    useNA = "ifany"
  )
)

cat("\n")


# ==============================================================================
# 33. CORRELATIONS WITH PRIMARY SCORE AND SRSq
# ==============================================================================

correlation_results_list <- list()

correlation_counter <- 1L


for (score_name in available_scores) {
  
  for (
    target_name in c(
      "primary_5gene_score",
      "SRSq"
    )
  ) {
    
    x <- bp_scores[[score_name]]
    y <- bp_scores[[target_name]]
    
    
    keep <- complete.cases(
      x,
      y
    )
    
    
    x2 <- x[keep]
    y2 <- y[keep]
    
    
    if (
      length(x2) < 8 ||
      length(unique(x2)) < 3 ||
      length(unique(y2)) < 3
    ) {
      
      next
    }
    
    
    ht <- suppressWarnings(
      stats::cor.test(
        x = x2,
        y = y2,
        method = "spearman",
        exact = FALSE
      )
    )
    
    
    correlation_results_list[[correlation_counter]] <- tibble(
      
      score_column =
        score_name,
      
      target =
        target_name,
      
      n =
        length(x2),
      
      Spearman_rho =
        safe_htest_estimate(
          ht
        ),
      
      p_value =
        safe_htest_pvalue(
          ht
        )
    )
    
    
    correlation_counter <- correlation_counter + 1L
  }
}


correlation_results <- bind_rows(
  correlation_results_list
) %>%
  left_join(
    score_registry,
    by = "score_column"
  )


if (nrow(correlation_results) > 0) {
  
  correlation_results$BH_correlation <-
    stats::p.adjust(
      correlation_results$p_value,
      method = "BH"
    )
  
  
  correlation_results <- correlation_results %>%
    arrange(
      target,
      p_value
    )
}


cat("Correlations with primary score and SRSq in BP:\n")

if (nrow(correlation_results) > 0) {
  
  print(
    correlation_results %>%
      select(
        display_name,
        target,
        n,
        Spearman_rho,
        p_value,
        BH_correlation
      ) %>%
      tibble::as_tibble(),
    n = Inf
  )
  
} else {
  
  cat("No correlation results generated.\n")
}

cat("\n")


# ==============================================================================
# 34. ENDOTYPE ASSOCIATIONS
# ==============================================================================

endotype_results_list <- list()

endotype_counter <- 1L


for (score_name in available_scores) {
  
  score_values <- bp_scores[[score_name]]
  
  
  # ---------------------------------------------------------------------------
  # SRS
  # ---------------------------------------------------------------------------
  
  keep_srs <- complete.cases(
    score_values,
    bp_scores$SRS
  )
  
  
  x_srs <- score_values[
    keep_srs
  ]
  
  
  srs_group <- droplevels(
    factor(
      bp_scores$SRS[
        keep_srs
      ]
    )
  )
  
  
  if (
    length(x_srs) >= 6 &&
    nlevels(srs_group) == 2
  ) {
    
    level_names <- levels(
      srs_group
    )
    
    
    values_1 <- x_srs[
      srs_group == level_names[1]
    ]
    
    
    values_2 <- x_srs[
      srs_group == level_names[2]
    ]
    
    
    ht <- stats::wilcox.test(
      x = values_1,
      y = values_2,
      exact = FALSE,
      paired = FALSE
    )
    
    
    endotype_results_list[[endotype_counter]] <- tibble(
      
      score_column =
        score_name,
      
      framework =
        "SRS",
      
      test =
        "Wilcoxon_rank_sum",
      
      n =
        length(x_srs),
      
      effect =
        median(
          values_1,
          na.rm = TRUE
        ) -
        median(
          values_2,
          na.rm = TRUE
        ),
      
      effect_name =
        paste0(
          "median_",
          level_names[1],
          "_minus_",
          level_names[2]
        ),
      
      p_value =
        safe_htest_pvalue(
          ht
        )
    )
    
    
    endotype_counter <- endotype_counter + 1L
  }
  
  
  # ---------------------------------------------------------------------------
  # CTS
  # ---------------------------------------------------------------------------
  
  keep_cts <- complete.cases(
    score_values,
    bp_scores$CTS
  )
  
  
  x_cts <- score_values[
    keep_cts
  ]
  
  
  cts_group <- droplevels(
    factor(
      bp_scores$CTS[
        keep_cts
      ]
    )
  )
  
  
  if (
    length(x_cts) >= 6 &&
    nlevels(cts_group) >= 2
  ) {
    
    ht <- stats::kruskal.test(
      x_cts ~ cts_group
    )
    
    
    epsilon2 <- epsilon_squared_kw(
      x = x_cts,
      group = cts_group
    )
    
    
    endotype_results_list[[endotype_counter]] <- tibble(
      
      score_column =
        score_name,
      
      framework =
        "CTS",
      
      test =
        "Kruskal_Wallis",
      
      n =
        length(x_cts),
      
      effect =
        epsilon2,
      
      effect_name =
        "epsilon_squared",
      
      p_value =
        safe_htest_pvalue(
          ht
        )
    )
    
    
    endotype_counter <- endotype_counter + 1L
  }
}


endotype_results <- bind_rows(
  endotype_results_list
) %>%
  left_join(
    score_registry,
    by = "score_column"
  )


if (nrow(endotype_results) > 0) {
  
  endotype_results$BH_endotype <-
    stats::p.adjust(
      endotype_results$p_value,
      method = "BH"
    )
  
  
  endotype_results <- endotype_results %>%
    arrange(
      framework,
      p_value
    )
}


cat("Signature associations with SRS and CTS:\n")

if (nrow(endotype_results) > 0) {
  
  print(
    endotype_results %>%
      select(
        display_name,
        framework,
        test,
        n,
        effect,
        effect_name,
        p_value,
        BH_endotype
      ) %>%
      tibble::as_tibble(),
    n = Inf
  )
  
} else {
  
  cat("No endotype association results generated.\n")
}

cat("\n")


# ==============================================================================
# 35. PAIRWISE SIGNATURE CORRELATION MATRIX IN BP
# ==============================================================================

bp_score_matrix <- bp_scores %>%
  select(
    all_of(
      available_scores
    )
  )


signature_correlation_matrix <- stats::cor(
  bp_score_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)


signature_correlation_df <- as.data.frame(
  signature_correlation_matrix
)


signature_correlation_df <- tibble::rownames_to_column(
  signature_correlation_df,
  var = "score_column"
)


# ==============================================================================
# 36. FIGURE A — AUC BENCHMARK
# ==============================================================================

auc_plot_data <- auc_results %>%
  filter(
    !is.na(
      AUC_fixed_direction
    )
  )


if (nrow(auc_plot_data) > 0) {
  
  p_auc <- ggplot(
    auc_plot_data,
    aes(
      x =
        reorder(
          display_name,
          AUC_fixed_direction
        ),
      
      y =
        AUC_fixed_direction
    )
  ) +
    
    geom_point(
      size = 3
    ) +
    
    geom_errorbar(
      aes(
        ymin =
          AUC_CI_low,
        
        ymax =
          AUC_CI_high
      ),
      width = 0.15
    ) +
    
    geom_hline(
      yintercept = 0.5,
      linetype = 2
    ) +
    
    coord_flip() +
    
    labs(
      title =
        "Contextual blood transcriptomic biomarker benchmarking",
      
      subtitle =
        paste0(
          "Sepsis versus healthy controls; ",
          "not independent diagnostic validation"
        ),
      
      x =
        NULL,
      
      y =
        "AUC"
    ) +
    
    theme_bw(
      base_size = 11
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "137_Figure_A_AUC_benchmark.png"
    ),
    plot = p_auc,
    width = 8,
    height = 5.5,
    dpi = 300
  )
}


# ==============================================================================
# 37. FIGURE B — AGE/SEX-ADJUSTED EFFECTS
# ==============================================================================

if (nrow(adjusted_results) > 0) {
  
  adjusted_plot_data <- adjusted_results %>%
    filter(
      !is.na(
        adjusted_standardized_BP_vs_BC_effect
      )
    )
  
  
  if (nrow(adjusted_plot_data) > 0) {
    
    p_adjusted <- ggplot(
      adjusted_plot_data,
      aes(
        x =
          adjusted_standardized_BP_vs_BC_effect,
        
        y =
          reorder(
            display_name,
            adjusted_standardized_BP_vs_BC_effect
          )
      )
    ) +
      
      geom_vline(
        xintercept = 0,
        linetype = 2
      ) +
      
      geom_point(
        size = 3
      ) +
      
      labs(
        title =
          "Age- and sex-adjusted biomarker score differences",
        
        subtitle =
          "Standardized score difference: sepsis minus healthy controls",
        
        x =
          "Adjusted standardized BP-BC effect",
        
        y =
          NULL
      ) +
      
      theme_bw(
        base_size = 11
      )
    
    
    ggsave(
      filename = file.path(
        figures_dir,
        "137_Figure_B_age_sex_adjusted_effects.png"
      ),
      plot = p_adjusted,
      width = 8,
      height = 5.5,
      dpi = 300
    )
  }
}


# ==============================================================================
# 38. FIGURE C — CORRELATIONS WITH PRIMARY SCORE / SRSq
# ==============================================================================

if (nrow(correlation_results) > 0) {
  
  correlation_plot_data <- correlation_results %>%
    filter(
      !is.na(
        Spearman_rho
      )
    )
  
  
  if (nrow(correlation_plot_data) > 0) {
    
    p_corr <- ggplot(
      correlation_plot_data,
      aes(
        x =
          Spearman_rho,
        
        y =
          reorder(
            display_name,
            Spearman_rho
          )
      )
    ) +
      
      geom_vline(
        xintercept = 0,
        linetype = 2
      ) +
      
      geom_point(
        size = 3
      ) +
      
      facet_wrap(
        ~ target
      ) +
      
      labs(
        title =
          "Alignment of biomarker signatures with the host-response axis",
        
        x =
          "Spearman rho",
        
        y =
          NULL
      ) +
      
      theme_bw(
        base_size = 10
      )
    
    
    ggsave(
      filename = file.path(
        figures_dir,
        "137_Figure_C_signature_correlations.png"
      ),
      plot = p_corr,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
}


# ==============================================================================
# 39. FIGURE D — STANDARDIZED SCORES BY CTS
# ==============================================================================

cts_plot_data <- bp_scores %>%
  select(
    sample_id,
    CTS,
    all_of(
      available_scores
    )
  ) %>%
  pivot_longer(
    cols =
      all_of(
        available_scores
      ),
    
    names_to =
      "score_column",
    
    values_to =
      "score"
  ) %>%
  filter(
    !is.na(
      CTS
    ),
    !is.na(
      score
    )
  ) %>%
  group_by(
    score_column
  ) %>%
  mutate(
    score_z =
      z_safe(
        score
      )
  ) %>%
  ungroup() %>%
  left_join(
    score_registry %>%
      select(
        score_column,
        display_name
      ),
    by = "score_column"
  )


if (nrow(cts_plot_data) > 0) {
  
  p_cts <- ggplot(
    cts_plot_data,
    aes(
      x =
        CTS,
      
      y =
        score_z
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      width = 0.55
    ) +
    
    geom_jitter(
      width = 0.12,
      size = 1.5
    ) +
    
    facet_wrap(
      ~ display_name,
      scales = "free_y"
    ) +
    
    labs(
      title =
        "Blood biomarker signatures across CTS classes",
      
      subtitle =
        "Scores standardized within the sepsis cohort",
      
      x =
        "Consensus Transcriptomic Subtype",
      
      y =
        "Standardized score"
    ) +
    
    theme_bw(
      base_size = 9
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "137_Figure_D_scores_by_CTS.png"
    ),
    plot = p_cts,
    width = 12,
    height = 8,
    dpi = 300
  )
}


# ==============================================================================
# 40. FIGURE E — GENE DIRECTION AUDIT
# ==============================================================================

gene_plot_data <- gene_audit %>%
  filter(
    present,
    !is.na(
      BP_minus_BC
    )
  )


if (nrow(gene_plot_data) > 0) {
  
  p_gene <- ggplot(
    gene_plot_data,
    aes(
      x =
        BP_minus_BC,
      
      y =
        reorder(
          canonical_gene,
          BP_minus_BC
        )
    )
  ) +
    
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    
    geom_point(
      size = 2.5
    ) +
    
    labs(
      title =
        "Direction audit of predefined biomarker genes",
      
      subtitle =
        "Median TMM logCPM difference: sepsis minus healthy controls",
      
      x =
        "Median BP - BC logCPM",
      
      y =
        NULL
    ) +
    
    theme_bw(
      base_size = 9
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "137_Figure_E_gene_direction_audit.png"
    ),
    plot = p_gene,
    width = 8,
    height = 8,
    dpi = 300
  )
}


# ==============================================================================
# 41. FIGURE F — ROC CURVES
# ==============================================================================

roc_coordinate_list <- list()

roc_counter <- 1L


for (score_name in names(roc_objects)) {
  
  roc_object <- roc_objects[[score_name]]
  
  
  roc_coordinate_list[[roc_counter]] <- tibble(
    
    score_column =
      score_name,
    
    false_positive_rate =
      1 -
      as.numeric(
        roc_object$specificities
      ),
    
    sensitivity =
      as.numeric(
        roc_object$sensitivities
      )
  )
  
  
  roc_counter <- roc_counter + 1L
}


roc_coordinates <- bind_rows(
  roc_coordinate_list
) %>%
  left_join(
    score_registry %>%
      select(
        score_column,
        display_name
      ),
    by = "score_column"
  )


if (nrow(roc_coordinates) > 0) {
  
  p_roc <- ggplot(
    roc_coordinates,
    aes(
      x =
        false_positive_rate,
      
      y =
        sensitivity,
      
      group =
        display_name,
      
      linetype =
        display_name
    )
  ) +
    
    geom_line(
      linewidth = 0.8
    ) +
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2
    ) +
    
    coord_equal() +
    
    labs(
      title =
        "Contextual BP-versus-BC ROC curves",
      
      subtitle =
        paste0(
          "Healthy-control comparison; ",
          "not independent diagnostic validation"
        ),
      
      x =
        "1 - specificity",
      
      y =
        "Sensitivity",
      
      linetype =
        "Signature"
    ) +
    
    theme_bw(
      base_size = 10
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "137_Figure_F_ROC_curves.png"
    ),
    plot = p_roc,
    width = 8,
    height = 7,
    dpi = 300
  )
}


# ==============================================================================
# 42. SAVE CSV TABLES
# ==============================================================================

write.csv(
  reference_table,
  file.path(
    tables_dir,
    "137_signature_reference_registry.csv"
  ),
  row.names = FALSE
)


write.csv(
  signature_status,
  file.path(
    tables_dir,
    "137_signature_gene_completeness.csv"
  ),
  row.names = FALSE
)


write.csv(
  gene_audit,
  file.path(
    tables_dir,
    "137_benchmark_gene_direction_audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  score_table,
  file.path(
    tables_dir,
    "137_all_blood_biomarker_scores.csv"
  ),
  row.names = FALSE
)


write.csv(
  auc_results,
  file.path(
    tables_dir,
    "137_BP_BC_AUC_benchmark.csv"
  ),
  row.names = FALSE
)


write.csv(
  adjusted_results,
  file.path(
    tables_dir,
    "137_age_sex_adjusted_benchmark.csv"
  ),
  row.names = FALSE
)


write.csv(
  correlation_results,
  file.path(
    tables_dir,
    "137_signature_correlations_primary_SRSq.csv"
  ),
  row.names = FALSE
)


write.csv(
  endotype_results,
  file.path(
    tables_dir,
    "137_signature_endotype_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  signature_correlation_df,
  file.path(
    tables_dir,
    "137_BP_signature_correlation_matrix.csv"
  ),
  row.names = FALSE
)


write.csv(
  roc_coordinates,
  file.path(
    tables_dir,
    "137_ROC_coordinates.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 43. PUBLICATION WORKBOOK
# ==============================================================================

wb <- openxlsx::createWorkbook()


run_info <- tibble(
  
  parameter = c(
    "script",
    "run_date",
    "cohort",
    "n_blood",
    "n_BP",
    "n_BC",
    "normalization",
    "AUC_direction",
    "primary_panel_status",
    "external_signature_status",
    "interpretation"
  ),
  
  value = c(
    script_name,
    as.character(run_date),
    "Blood RNA-seq only",
    "45",
    "35",
    "10",
    "edgeR TMM + logCPM, prior.count=2",
    "Higher predefined score interpreted as more sepsis-like",
    "Frozen from Script 135; selected in current cohort",
    "Predefined published signatures; platform-adapted implementations",
    paste0(
      "Contextual benchmarking only. ",
      "Healthy controls are not equivalent to noninfectious ICU/SIRS controls."
    )
  )
)


sheet_list <- list(
  
  "00_run_info" =
    run_info,
  
  "01_references" =
    reference_table,
  
  "02_signature_status" =
    signature_status,
  
  "03_gene_audit" =
    gene_audit,
  
  "04_sample_scores" =
    score_table,
  
  "05_AUC_benchmark" =
    auc_results,
  
  "06_age_sex_adjusted" =
    adjusted_results,
  
  "07_corr_primary_SRSq" =
    correlation_results,
  
  "08_endotype_tests" =
    endotype_results,
  
  "09_BP_score_correlations" =
    signature_correlation_df
)


for (sheet_name in names(sheet_list)) {
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  openxlsx::writeData(
    wb,
    sheet_name,
    sheet_list[[sheet_name]]
  )
}


workbook_file <- file.path(
  tables_dir,
  "137_blood_biomarker_benchmarking.xlsx"
)


openxlsx::saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)


# ==============================================================================
# 44. AUTOMATED TEXT SUMMARY
# ==============================================================================

top_auc <- auc_results %>%
  arrange(
    desc(
      AUC_fixed_direction
    )
  )


summary_ru <- c(
  
  "SCRIPT 137 — BLOOD BIOMARKER BENCHMARKING",
  
  "====================================================================",
  
  "",
  
  "АНАЛИЗ:",
  
  "Blood RNA-seq only: BP n=35, BC n=10.",
  
  paste0(
    "Доступно biomarker scores: ",
    length(
      available_scores
    ),
    "."
  ),
  
  "",
  
  "ВАЖНО:",
  
  paste0(
    "Все BP-vs-BC AUC отражают apparent/contextual discrimination ",
    "в текущей когорте."
  ),
  
  paste0(
    "Healthy controls не эквивалентны noninfectious ICU/SIRS controls ",
    "для которых разрабатывалась часть опубликованных signatures."
  ),
  
  paste0(
    "Primary five-gene panel был выбран в этой же когорте и ",
    "не является independently validated diagnostic test."
  ),
  
  "",
  
  "AUC BENCHMARK:"
)


for (i in seq_len(nrow(top_auc))) {
  
  summary_ru <- c(
    
    summary_ru,
    
    paste0(
      top_auc$display_name[i],
      ": AUC=",
      signif(
        top_auc$AUC_fixed_direction[i],
        4
      ),
      "; 95% CI ",
      signif(
        top_auc$AUC_CI_low[i],
        4
      ),
      "-",
      signif(
        top_auc$AUC_CI_high[i],
        4
      )
    )
  )
}


summary_ru <- c(
  
  summary_ru,
  
  "",
  
  "INTERPRETATION:",
  
  paste0(
    "Основная цель Script 137 — оценить biological concordance ",
    "нашей фиксированной пятигенной оси с опубликованными ",
    "host-response signatures, SRSq и CTS, а не доказать ",
    "превосходство диагностической точности."
  )
)


summary_ru_file <- file.path(
  text_dir,
  "137_summary_RU.txt"
)


writeLines(
  summary_ru,
  summary_ru_file
)


summary_en <- c(
  
  "SCRIPT 137 — BLOOD BIOMARKER BENCHMARKING",
  
  "====================================================================",
  
  "",
  
  "Blood RNA-seq only: BP n=35, BC n=10.",
  
  "",
  
  paste0(
    "The analysis provides contextual benchmarking against predefined ",
    "published host-response signatures."
  ),
  
  paste0(
    "BP-versus-BC AUC values represent apparent discrimination against ",
    "healthy controls and should not be interpreted as independent ",
    "diagnostic validation."
  ),
  
  paste0(
    "For signatures originally developed to distinguish sepsis or infection ",
    "from noninfectious systemic inflammation, the present healthy-control ",
    "comparison represents a different classification task."
  ),
  
  paste0(
    "The main purpose of the benchmarking analysis is therefore biological ",
    "concordance with the five-gene host-response axis, SRSq and CTS."
  )
)


summary_en_file <- file.path(
  text_dir,
  "137_summary_EN.txt"
)


writeLines(
  summary_en,
  summary_en_file
)


# ==============================================================================
# 45. INPUT MANIFEST
# ==============================================================================

input_manifest <- tibble(
  
  input = c(
    "raw_counts",
    "Script135_scores",
    "Script136b_demographics"
  ),
  
  path = c(
    counts_file,
    script135_scores_file,
    demographics_file
  )
)


file_information <- file.info(
  input_manifest$path
)


input_manifest$file_size_bytes <-
  file_information$size


input_manifest$modified_time <-
  as.character(
    file_information$mtime
  )


input_manifest$md5 <-
  unname(
    tools::md5sum(
      input_manifest$path
    )
  )


write.csv(
  input_manifest,
  file.path(
    logs_dir,
    "137_input_file_manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 46. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "137_sessionInfo.txt"
  )
)


# ==============================================================================
# 47. FINAL REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 137 completed successfully.\n")
cat("====================================================================\n\n")


cat("SIGNATURE COMPLETENESS:\n")

print(
  signature_status,
  n = Inf
)

cat("\n")


cat("GENE DIRECTION AUDIT:\n")

print(
  gene_audit %>%
    select(
      canonical_gene,
      resolved_symbol,
      expected_direction,
      observed_direction,
      direction_concordant,
      BP_minus_BC,
      wilcoxon_p,
      BH_gene_audit
    ) %>%
    tibble::as_tibble(),
  n = Inf
)

cat("\n")


cat("BP vs BC AUC benchmark:\n")

print(
  auc_results %>%
    select(
      display_name,
      n_total,
      AUC_fixed_direction,
      AUC_CI_low,
      AUC_CI_high,
      Wilcoxon_p,
      BH_Wilcoxon
    ) %>%
    tibble::as_tibble(),
  n = Inf
)

cat("\n")


cat("Age/sex-adjusted score effects:\n")

if (nrow(adjusted_results) > 0) {
  
  print(
    adjusted_results %>%
      select(
        display_name,
        n,
        adjusted_standardized_BP_vs_BC_effect,
        p_value,
        BH_adjusted,
        adjusted_R2
      ) %>%
      tibble::as_tibble(),
    n = Inf
  )
  
} else {
  
  cat("No adjusted score results generated.\n")
}

cat("\n")


cat("Correlations with primary score and SRSq in BP:\n")

if (nrow(correlation_results) > 0) {
  
  print(
    correlation_results %>%
      select(
        display_name,
        target,
        n,
        Spearman_rho,
        p_value,
        BH_correlation
      ) %>%
      tibble::as_tibble(),
    n = Inf
  )
  
} else {
  
  cat("No correlation results generated.\n")
}

cat("\n")


cat("Endotype associations:\n")

if (nrow(endotype_results) > 0) {
  
  print(
    endotype_results %>%
      select(
        display_name,
        framework,
        n,
        effect,
        effect_name,
        p_value,
        BH_endotype
      ) %>%
      tibble::as_tibble(),
    n = Inf
  )
  
} else {
  
  cat("No endotype association results generated.\n")
}

cat("\n")


cat("Main output workbook:\n")

cat(
  normalizePath(
    workbook_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat("Key output tables:\n")

cat(
  "1) ",
  normalizePath(
    file.path(
      tables_dir,
      "137_BP_BC_AUC_benchmark.csv"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "2) ",
  normalizePath(
    file.path(
      tables_dir,
      "137_age_sex_adjusted_benchmark.csv"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "3) ",
  normalizePath(
    file.path(
      tables_dir,
      "137_signature_correlations_primary_SRSq.csv"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "4) ",
  normalizePath(
    file.path(
      tables_dir,
      "137_signature_endotype_associations.csv"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "5) ",
  normalizePath(
    file.path(
      tables_dir,
      "137_benchmark_gene_direction_audit.csv"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat("Figures:\n")
cat("- 137_Figure_A_AUC_benchmark.png\n")
cat("- 137_Figure_B_age_sex_adjusted_effects.png\n")
cat("- 137_Figure_C_signature_correlations.png\n")
cat("- 137_Figure_D_scores_by_CTS.png\n")
cat("- 137_Figure_E_gene_direction_audit.png\n")
cat("- 137_Figure_F_ROC_curves.png\n\n")


cat("IMPORTANT:\n")
cat("- Primary five-gene panel remains frozen from Script 135.\n")
cat("- No additional feature selection was performed.\n")
cat("- External signatures are predefined literature comparators.\n")
cat("- SeptiCyte-derived RNA-seq contrasts are NOT official SeptiScores.\n")
cat("- BP vs BC means sepsis versus healthy controls.\n")
cat("- AUCs are apparent/contextual, not independent validation.\n")
cat("- Do not interpret these AUCs as sepsis-versus-SIRS validation.\n")
cat("- Age/sex-adjusted sensitivity analyses are included.\n")
cat("- Main biological question is concordance with SRSq and CTS.\n")
cat("- Blood only; no urine; no lncRNA.\n\n")


cat("Done.\n")