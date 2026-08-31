################################################################################
# Script 153
# FINAL v2
#
# Supplementary Figure S1
#
# Discovery cohort quality control and sequencing structure
#
# Project:
#   Sepsis_DESeq2
#
# Manuscript:
#   Blood-only sepsis transcriptomic endotypes /
#   five-gene host-response signature
#
#
# IMPORTANT PROVENANCE CORRECTION
# -------------------------------
#
# Original blood DESeq2 Scripts 08 and 09:
#
#   1. selected BP + BC samples;
#   2. applied expression prefilter:
#
#        count >= 10 in >= 3 blood samples
#
#   3. did NOT remove ERCC-prefixed technical features before this filter.
#
# Therefore:
#
#   original frozen DESeq2 universe = 12,400 retained FEATURES
#
# comprising:
#
#   12,393 biological/human targets
#        7 ERCC-prefixed technical features
#
# Ten ERCC-prefixed rows exist in the complete input matrix; seven passed
# the original expression prefilter.
#
# None of these seven ERCC-prefixed retained features is a primary DEG,
# batch-adjusted DEG, or robust-core DEG.
#
#
# FIGURE S1 BIOLOGICAL QC MATRIX
# ------------------------------
#
# For unsupervised biological QC in this supplementary figure:
#
#   - all ERCC-prefixed features are removed;
#   - the same expression criterion is applied;
#   - 12,393 biological targets are retained.
#
#
# FIGURE S1 PANELS
# ----------------
#
# A. Human-target library-size distribution by study group
#
# B. PCA of VST-transformed biological targets,
#    colored by sequencing batch and shaped by study group
#
# C. VST sample-distance heatmap with condition and batch annotations
#
# D. Number of blood samples per sequencing batch stratified by group
#
#
# THIS SCRIPT DOES NOT:
# ---------------------
#
#   - re-run blood differential-expression testing;
#   - modify frozen DESeq2 results;
#   - claim absence of sequencing-batch effects;
#   - remove biological batch-associated variation;
#   - perform hypothesis testing on the QC panels.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 153 FINAL v2\n")
cat("Supplementary Figure S1\n")
cat("Discovery blood cohort QC and sequencing structure\n")
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
  "ggplot2",
  "pheatmap",
  "gridExtra",
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
  
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(gridExtra)
  library(openxlsx)
  
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "153_FigureS1_discovery_QC"
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

find_counts_file <- function() {
  
  preferred_candidates <- c(
    
    file.path(
      project_dir,
      "data",
      "counts_all.csv"
    ),
    
    file.path(
      project_dir,
      "counts_all.csv"
    ),
    
    file.path(
      project_dir,
      "data",
      "counts",
      "counts_all.csv"
    ),
    
    file.path(
      project_dir,
      "data",
      "expression",
      "counts_all.csv"
    ),
    
    file.path(
      project_dir,
      "data",
      "processed",
      "counts_all.csv"
    )
  )
  
  
  hits <- preferred_candidates[
    file.exists(
      preferred_candidates
    )
  ]
  
  
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
    ) >
    1
  ) {
    
    cat(
      "\nMultiple preferred counts_all.csv files detected:\n"
    )
    
    
    print(
      hits
    )
    
    
    stop(
      "Resolve counts_all.csv provenance before continuing."
    )
  }
  
  
  all_hits <- list.files(
    project_dir,
    pattern = "^counts_all\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  
  all_hits <- all_hits[
    !grepl(
      "[/\\\\]results[/\\\\]",
      all_hits,
      ignore.case = TRUE
    )
  ]
  
  
  all_hits <- sort(
    unique(
      all_hits
    )
  )
  
  
  if (
    length(
      all_hits
    ) ==
    1
  ) {
    
    return(
      all_hits[1]
    )
  }
  
  
  if (
    length(
      all_hits
    ) ==
    0
  ) {
    
    stop(
      "counts_all.csv was not found."
    )
  }
  
  
  cat(
    "\nMultiple counts_all.csv files found:\n"
  )
  
  
  print(
    all_hits
  )
  
  
  stop(
    "Resolve counts_all.csv provenance before continuing."
  )
}


find_gene_column <- function(data) {
  
  candidates <- c(
    "Gene",
    "gene",
    "GeneSymbol",
    "gene_symbol",
    "symbol"
  )
  
  
  hit <- candidates[
    candidates %in%
      names(data)
  ]
  
  
  if (
    length(
      hit
    ) ==
    0
  ) {
    
    stop(
      paste0(
        "Could not identify gene column. First columns: ",
        paste(
          head(
            names(data),
            10
          ),
          collapse = ", "
        )
      )
    )
  }
  
  
  hit[1]
}


normalize_gene <- function(x) {
  
  toupper(
    trimws(
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
    trimws(
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


# =============================================================================
# 5. FROZEN BLOOD SAMPLE METADATA FROM SCRIPT 151
# =============================================================================

metadata_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "151_Table1_TableS1_discovery_cohort",
  "tables",
  "151_TableS1_blood_sample_metadata.csv"
)


if (
  !file.exists(
    metadata_file
  )
) {
  
  stop(
    paste0(
      "Frozen blood metadata from Script 151 not found:\n",
      metadata_file
    )
  )
}


sample_meta <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_metadata_columns <- c(
  "sample_id",
  "participant_id",
  "cohort",
  "biofluid",
  "sample_group",
  "sequencing_batch"
)


missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(
    sample_meta
  )
)


if (
  length(
    missing_metadata_columns
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing metadata column(s): ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    )
  )
}


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
      "Expected 45 blood samples; observed ",
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
    "Duplicate blood sample IDs detected."
  )
}


if (
  sum(
    sample_meta$sample_group ==
    "BP"
  ) !=
  35
) {
  
  stop(
    "Expected BP n=35."
  )
}


if (
  sum(
    sample_meta$sample_group ==
    "BC"
  ) !=
  10
) {
  
  stop(
    "Expected BC n=10."
  )
}


if (
  any(
    is.na(
      sample_meta$sequencing_batch
    ) |
    trimws(
      sample_meta$sequencing_batch
    ) ==
    ""
  )
) {
  
  stop(
    "Missing sequencing-batch values detected."
  )
}


sample_meta <- sample_meta %>%
  
  dplyr::mutate(
    
    condition_label =
      dplyr::case_when(
        
        sample_group ==
          "BP" ~
          "Sepsis (BP)",
        
        sample_group ==
          "BC" ~
          "Healthy control (BC)",
        
        TRUE ~
          NA_character_
      ),
    
    condition_label =
      factor(
        condition_label,
        levels = c(
          "Sepsis (BP)",
          "Healthy control (BC)"
        )
      ),
    
    sequencing_batch =
      factor(
        as.character(
          sequencing_batch
        ),
        levels = sort(
          unique(
            as.character(
              sequencing_batch
            )
          )
        )
      )
  )


sample_ids <- as.character(
  sample_meta$sample_id
)


cat("\nFrozen metadata source:\n")

cat(
  normalizePath(
    metadata_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


cat(
  "Blood samples = ",
  nrow(
    sample_meta
  ),
  " (BP ",
  sum(
    sample_meta$sample_group ==
      "BP"
  ),
  "; BC ",
  sum(
    sample_meta$sample_group ==
      "BC"
  ),
  ")\n",
  sep = ""
)


# =============================================================================
# 6. RAW COUNT MATRIX
# =============================================================================

counts_file <- find_counts_file()


cat("\nRaw counts source:\n")

cat(
  normalizePath(
    counts_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


counts_all <- read.csv(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


gene_col <- find_gene_column(
  counts_all
)


genes_all <- as.character(
  counts_all[[gene_col]]
)


missing_samples <- setdiff(
  sample_ids,
  names(
    counts_all
  )
)


if (
  length(
    missing_samples
  ) >
  0
) {
  
  stop(
    paste0(
      "Blood sample(s) missing from counts_all.csv: ",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}


counts_df <- counts_all[
  sample_ids
]


counts_df[] <- lapply(
  counts_df,
  function(x) {
    
    suppressWarnings(
      as.numeric(x)
    )
  }
)


counts_blood <- as.matrix(
  counts_df
)


if (
  any(
    !is.finite(
      counts_blood
    )
  )
) {
  
  stop(
    "Non-finite raw-count values detected."
  )
}


if (
  any(
    counts_blood <
    0
  )
) {
  
  stop(
    "Negative raw-count values detected."
  )
}


if (
  max(
    abs(
      counts_blood -
      round(
        counts_blood
      )
    )
  ) >
  1e-6
) {
  
  stop(
    "Non-integer values detected in raw count matrix."
  )
}


counts_blood <- round(
  counts_blood
)


storage.mode(
  counts_blood
) <- "integer"


rownames(
  counts_blood
) <- make.unique(
  genes_all
)


n_raw_features <- nrow(
  counts_blood
)


cat(
  "Raw matrix features = ",
  n_raw_features,
  "\n",
  sep = ""
)


if (
  n_raw_features !=
  20812
) {
  
  stop(
    paste0(
      "Expected 20,812 raw rows; observed ",
      n_raw_features
    )
  )
}


# =============================================================================
# 7. ERCC AUDIT
# =============================================================================

ercc_flag <- grepl(
  "^ERCC[-_]",
  genes_all,
  ignore.case = TRUE
)


n_ercc_total <- sum(
  ercc_flag
)


cat(
  "ERCC-prefixed technical features in complete matrix = ",
  n_ercc_total,
  "\n",
  sep = ""
)


if (
  n_ercc_total !=
  10
) {
  
  stop(
    paste0(
      "Expected 10 ERCC-prefixed rows; observed ",
      n_ercc_total
    )
  )
}


# =============================================================================
# 8. REPRODUCE ORIGINAL BLOOD DE EXPRESSION PREFILTER
# =============================================================================
#
# IMPORTANT:
#
# This reproduces Scripts 08/09 exactly with respect to filtering:
#
#   count >=10 in >=3 blood samples
#
# applied BEFORE exclusion of ERCC-prefixed technical rows.
#
# =============================================================================

original_prefilter_keep <- rowSums(
  counts_blood >=
    10
) >=
  3


n_original_prefilter <- sum(
  original_prefilter_keep
)


n_ercc_passing_prefilter <- sum(
  original_prefilter_keep &
    ercc_flag
)


n_biological_passing_prefilter <- sum(
  original_prefilter_keep &
    !ercc_flag
)


cat("\nORIGINAL BLOOD DE PREFILTER AUDIT\n")
cat("---------------------------------\n")


cat(
  "All retained features = ",
  n_original_prefilter,
  "\n",
  sep = ""
)


cat(
  "Retained ERCC-prefixed technical features = ",
  n_ercc_passing_prefilter,
  "\n",
  sep = ""
)


cat(
  "Retained biological/human targets = ",
  n_biological_passing_prefilter,
  "\n",
  sep = ""
)


if (
  n_original_prefilter !=
  12400
) {
  
  stop(
    paste0(
      "Original prefilter should retain 12,400 features; observed ",
      n_original_prefilter
    )
  )
}


if (
  n_ercc_passing_prefilter !=
  7
) {
  
  stop(
    paste0(
      "Expected 7 ERCC-prefixed retained technical features; observed ",
      n_ercc_passing_prefilter
    )
  )
}


if (
  n_biological_passing_prefilter !=
  12393
) {
  
  stop(
    paste0(
      "Expected 12,393 retained biological targets; observed ",
      n_biological_passing_prefilter
    )
  )
}


cat(
  "Original DE prefilter provenance reproduced exactly.\n"
)


# =============================================================================
# 9. AUDIT AGAINST FROZEN MODEL-COMPARISON TABLE
# =============================================================================

model_comparison_file <- file.path(
  project_dir,
  "results",
  "blood_BP_vs_BC_model_comparison",
  "blood_model_comparison_all_genes.csv"
)


if (
  !file.exists(
    model_comparison_file
  )
) {
  
  stop(
    paste0(
      "Frozen model-comparison file not found:\n",
      model_comparison_file
    )
  )
}


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
      "Expected 12,400 rows in frozen model comparison; observed ",
      nrow(
        model_df
      )
    )
  )
}


model_gene_col <- find_gene_column(
  model_df
)


model_genes <- as.character(
  model_df[[model_gene_col]]
)


model_ercc_flag <- grepl(
  "^ERCC[-_]",
  model_genes,
  ignore.case = TRUE
)


n_model_ercc <- sum(
  model_ercc_flag
)


cat("\nFROZEN MODEL-COMPARISON AUDIT\n")
cat("-----------------------------\n")


cat(
  "Frozen model rows = ",
  nrow(
    model_df
  ),
  "\n",
  sep = ""
)


cat(
  "ERCC-prefixed rows in frozen model = ",
  n_model_ercc,
  "\n",
  sep = ""
)


cat(
  "Biological rows in frozen model = ",
  sum(
    !model_ercc_flag
  ),
  "\n",
  sep = ""
)


if (
  n_model_ercc !=
  7
) {
  
  stop(
    paste0(
      "Expected 7 retained ERCC rows in frozen model table; observed ",
      n_model_ercc
    )
  )
}


raw_prefilter_genes <- normalize_gene(
  genes_all[
    original_prefilter_keep
  ]
)


frozen_model_genes <- normalize_gene(
  model_genes
)


if (
  !setequal(
    raw_prefilter_genes,
    frozen_model_genes
  )
) {
  
  missing_from_model <- setdiff(
    raw_prefilter_genes,
    frozen_model_genes
  )
  
  
  extra_in_model <- setdiff(
    frozen_model_genes,
    raw_prefilter_genes
  )
  
  
  cat(
    "\nFeatures missing from frozen model table:\n"
  )
  
  print(
    missing_from_model
  )
  
  
  cat(
    "\nUnexpected features in frozen model table:\n"
  )
  
  print(
    extra_in_model
  )
  
  
  stop(
    "Frozen model feature universe does not match the reproduced raw-count prefilter."
  )
}


cat(
  "Frozen model feature set exactly matches reproduced original prefilter.\n"
)


# =============================================================================
# 10. AUDIT ERCC DIFFERENTIAL-EXPRESSION FLAGS
# =============================================================================

primary_deg_col <- if (
  "DEG_simple" %in%
  names(
    model_df
  )
) {
  "DEG_simple"
} else {
  NA_character_
}


batch_deg_col <- if (
  "DEG_batch" %in%
  names(
    model_df
  )
) {
  "DEG_batch"
} else {
  NA_character_
}


core_deg_col <- if (
  "core_DEG" %in%
  names(
    model_df
  )
) {
  "core_DEG"
} else {
  NA_character_
}


if (
  any(
    is.na(
      c(
        primary_deg_col,
        batch_deg_col,
        core_deg_col
      )
    )
  )
) {
  
  stop(
    "Could not identify frozen DEG flag columns."
  )
}


model_primary_deg <- as_logical_flag(
  model_df[[primary_deg_col]]
)


model_batch_deg <- as_logical_flag(
  model_df[[batch_deg_col]]
)


model_core_deg <- as_logical_flag(
  model_df[[core_deg_col]]
)


ercc_primary_deg_n <- sum(
  model_primary_deg[
    model_ercc_flag
  ]
)


ercc_batch_deg_n <- sum(
  model_batch_deg[
    model_ercc_flag
  ]
)


ercc_core_deg_n <- sum(
  model_core_deg[
    model_ercc_flag
  ]
)


cat(
  "ERCC primary DEG = ",
  ercc_primary_deg_n,
  "\n",
  sep = ""
)


cat(
  "ERCC batch-adjusted DEG = ",
  ercc_batch_deg_n,
  "\n",
  sep = ""
)


cat(
  "ERCC robust-core DEG = ",
  ercc_core_deg_n,
  "\n",
  sep = ""
)


if (
  ercc_primary_deg_n !=
  0 ||
  ercc_batch_deg_n !=
  0 ||
  ercc_core_deg_n !=
  0
) {
  
  stop(
    "At least one retained ERCC technical feature is classified as a DEG."
  )
}


cat(
  "All retained ERCC technical features are non-DE in both models.\n"
)


# =============================================================================
# 11. ERCC AUDIT TABLE
# =============================================================================

ercc_model_audit <- model_df[
  model_ercc_flag,
  ,
  drop = FALSE
]


# =============================================================================
# 12. BIOLOGICAL MATRIX FOR FIGURE S1
# =============================================================================
#
# Use only biological/human targets that passed the original blood
# expression threshold.
#
# =============================================================================

biological_qc_keep <-
  original_prefilter_keep &
  !ercc_flag


counts_qc <- counts_blood[
  biological_qc_keep,
  ,
  drop = FALSE
]


genes_qc <- genes_all[
  biological_qc_keep
]


rownames(
  counts_qc
) <- make.unique(
  genes_qc
)


if (
  nrow(
    counts_qc
  ) !=
  12393
) {
  
  stop(
    "Biological QC matrix should contain 12,393 targets."
  )
}


# =============================================================================
# 13. HUMAN-TARGET LIBRARY SIZE
# =============================================================================
#
# Library-size QC excludes all ERCC-prefixed technical rows,
# but is calculated before the expression threshold.
#
# =============================================================================

counts_human_all <- counts_blood[
  !ercc_flag,
  ,
  drop = FALSE
]


human_library_size <- colSums(
  counts_human_all
)


if (
  any(
    human_library_size <=
    0
  )
) {
  
  stop(
    "Zero human-target library size detected."
  )
}


sample_qc <- sample_meta %>%
  
  dplyr::mutate(
    
    human_target_counts =
      as.numeric(
        human_library_size[
          sample_id
        ]
      ),
    
    human_target_counts_million =
      human_target_counts /
      1e6
  )


# =============================================================================
# 14. DESEQ2 VST FOR UNSUPERVISED BIOLOGICAL QC
# =============================================================================

col_data <- data.frame(
  
  condition =
    factor(
      sample_meta$sample_group,
      levels = c(
        "BC",
        "BP"
      )
    ),
  
  sequencing_batch =
    sample_meta$sequencing_batch,
  
  row.names =
    sample_meta$sample_id
)


if (
  !identical(
    colnames(
      counts_qc
    ),
    rownames(
      col_data
    )
  )
) {
  
  stop(
    "Count matrix and metadata sample order do not match."
  )
}


dds_qc <- DESeq2::DESeqDataSetFromMatrix(
  
  countData =
    counts_qc,
  
  colData =
    col_data,
  
  design =
    ~ condition
)


vsd_qc <- DESeq2::vst(
  dds_qc,
  blind = TRUE
)


vst_matrix <- SummarizedExperiment::assay(
  vsd_qc
)


if (
  nrow(
    vst_matrix
  ) !=
  12393
) {
  
  stop(
    "Unexpected VST matrix row count."
  )
}


# =============================================================================
# 15. PCA
# =============================================================================

pca_fit <- stats::prcomp(
  t(
    vst_matrix
  ),
  center = TRUE,
  scale. = FALSE
)


percent_variance <- 100 *
  (
    pca_fit$sdev^2 /
      sum(
        pca_fit$sdev^2
      )
  )


pc1_percent <- percent_variance[1]

pc2_percent <- percent_variance[2]


pca_df <- data.frame(
  
  sample_id =
    rownames(
      pca_fit$x
    ),
  
  PC1 =
    pca_fit$x[
      ,
      1
    ],
  
  PC2 =
    pca_fit$x[
      ,
      2
    ],
  
  stringsAsFactors = FALSE
) %>%
  
  dplyr::left_join(
    
    sample_meta %>%
      
      dplyr::select(
        sample_id,
        participant_id,
        cohort,
        sample_group,
        condition_label,
        sequencing_batch
      ),
    
    by =
      "sample_id"
  )


cat("\nPCA variance:\n")


cat(
  "PC1 = ",
  sprintf(
    "%.3f",
    pc1_percent
  ),
  "%\n",
  sep = ""
)


cat(
  "PC2 = ",
  sprintf(
    "%.3f",
    pc2_percent
  ),
  "%\n",
  sep = ""
)


# =============================================================================
# 16. SAMPLE-DISTANCE MATRIX
# =============================================================================

sample_distance <- stats::dist(
  t(
    vst_matrix
  ),
  method = "euclidean"
)


sample_distance_matrix <- as.matrix(
  sample_distance
)


sample_distance_export <- data.frame(
  
  sample_id =
    rownames(
      sample_distance_matrix
    ),
  
  sample_distance_matrix,
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 17. SEQUENCING BATCH × CONDITION
# =============================================================================

batch_counts <- sample_meta %>%
  
  dplyr::count(
    sequencing_batch,
    condition_label,
    name = "n"
  ) %>%
  
  dplyr::arrange(
    sequencing_batch,
    condition_label
  )


batch_condition_table <- table(
  sample_meta$sequencing_batch,
  sample_meta$condition_label
)


batch_condition_export <- data.frame(
  
  sequencing_batch =
    rownames(
      batch_condition_table
    ),
  
  as.data.frame.matrix(
    batch_condition_table
  ),
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
)


cat("\nSequencing batch x study group:\n")

print(
  batch_condition_export,
  row.names = FALSE
)


# =============================================================================
# 18. COLORS AND SHAPES
# =============================================================================

condition_colors <- c(
  
  "Sepsis (BP)" =
    "#B2182B",
  
  "Healthy control (BC)" =
    "#2166AC"
)


condition_shapes <- c(
  
  "Sepsis (BP)" =
    16,
  
  "Healthy control (BC)" =
    17
)


batch_levels <- levels(
  sample_meta$sequencing_batch
)


n_batches <- length(
  batch_levels
)


batch_colors <- grDevices::hcl.colors(
  n_batches,
  palette = "Dark 3"
)


names(
  batch_colors
) <- batch_levels


# =============================================================================
# 19. PANEL A — HUMAN-TARGET LIBRARY SIZE
# =============================================================================

panel_A <- ggplot2::ggplot(
  
  sample_qc,
  
  ggplot2::aes(
    x = condition_label,
    y = human_target_counts_million,
    fill = condition_label
  )
  
) +
  
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.70
  ) +
  
  ggplot2::geom_jitter(
    
    ggplot2::aes(
      color = condition_label
    ),
    
    width = 0.12,
    size = 2.2,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values = condition_colors
  ) +
  
  ggplot2::scale_color_manual(
    values = condition_colors
  ) +
  
  ggplot2::labs(
    title = "A  Human-target library size",
    x = NULL,
    y = "Total human-target counts (millions)"
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
    
    axis.text.x =
      ggplot2::element_text(
        angle = 15,
        hjust = 1
      )
  )


# =============================================================================
# 20. PANEL B — PCA COLORED BY SEQUENCING BATCH
# =============================================================================

panel_B <- ggplot2::ggplot(
  
  pca_df,
  
  ggplot2::aes(
    x = PC1,
    y = PC2,
    color = sequencing_batch,
    shape = condition_label
  )
  
) +
  
  ggplot2::geom_point(
    size = 3.2,
    alpha = 0.90
  ) +
  
  ggplot2::scale_color_manual(
    values = batch_colors,
    name = "Sequencing batch"
  ) +
  
  ggplot2::scale_shape_manual(
    values = condition_shapes,
    name = "Study group"
  ) +
  
  ggplot2::labs(
    
    title =
      "B  PCA colored by sequencing batch",
    
    x =
      paste0(
        "PC1 (",
        sprintf(
          "%.1f",
          pc1_percent
        ),
        "%)"
      ),
    
    y =
      paste0(
        "PC2 (",
        sprintf(
          "%.1f",
          pc2_percent
        ),
        "%)"
      )
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
      "right"
  )


# =============================================================================
# 21. PANEL C — SAMPLE-DISTANCE HEATMAP
# =============================================================================

heatmap_annotation <- data.frame(
  
  Study_group =
    sample_meta$condition_label,
  
  Sequencing_batch =
    sample_meta$sequencing_batch,
  
  row.names =
    sample_meta$sample_id
)


annotation_colors <- list(
  
  Study_group =
    condition_colors,
  
  Sequencing_batch =
    batch_colors
)


heatmap_colors <- grDevices::colorRampPalette(
  c(
    "#FFFFFF",
    "#D9EAF7",
    "#74A9CF",
    "#2B8CBE",
    "#045A8D"
  )
)(
  100
)


heatmap_result <- pheatmap::pheatmap(
  
  sample_distance_matrix,
  
  color =
    heatmap_colors,
  
  cluster_rows =
    TRUE,
  
  cluster_cols =
    TRUE,
  
  annotation_row =
    heatmap_annotation,
  
  annotation_col =
    heatmap_annotation,
  
  annotation_colors =
    annotation_colors,
  
  show_rownames =
    FALSE,
  
  show_colnames =
    FALSE,
  
  border_color =
    NA,
  
  main =
    "C  VST sample-distance matrix",
  
  fontsize =
    9,
  
  silent =
    TRUE
)


# =============================================================================
# 22. PANEL D — SEQUENCING-BATCH COMPOSITION
# =============================================================================

panel_D <- ggplot2::ggplot(
  
  batch_counts,
  
  ggplot2::aes(
    x = sequencing_batch,
    y = n,
    fill = condition_label
  )
  
) +
  
  ggplot2::geom_col(
    width = 0.72
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        ifelse(
          n >
            0,
          n,
          ""
        )
    ),
    
    position =
      ggplot2::position_stack(
        vjust = 0.5
      ),
    
    size =
      3.2
  ) +
  
  ggplot2::scale_fill_manual(
    values = condition_colors,
    name = "Study group"
  ) +
  
  ggplot2::labs(
    title = "D  Samples per sequencing batch",
    x = "Sequencing batch",
    y = "Number of blood samples"
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
      "right",
    
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
  )


# =============================================================================
# 23. COMBINE FIGURE S1
# =============================================================================

figure_S1_grob <- gridExtra::arrangeGrob(
  
  panel_A,
  panel_B,
  heatmap_result$gtable,
  panel_D,
  
  ncol =
    2,
  
  widths =
    c(
      1,
      1
    ),
  
  heights =
    c(
      1,
      1.15
    )
)


# =============================================================================
# 24. SAVE COMBINED FIGURE
# =============================================================================

figure_pdf <- file.path(
  figures_dir,
  "153_FigureS1_discovery_blood_QC.pdf"
)


figure_png <- file.path(
  figures_dir,
  "153_FigureS1_discovery_blood_QC.png"
)


figure_tiff <- file.path(
  figures_dir,
  "153_FigureS1_discovery_blood_QC.tiff"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure_S1_grob,
  width = 14,
  height = 11,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure_S1_grob,
  width = 14,
  height = 11,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure_S1_grob,
  width = 14,
  height = 11,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 25. SAVE INDIVIDUAL PANELS
# =============================================================================

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "153_FigureS1A_library_size.pdf"
  ),
  plot = panel_A,
  width = 6.5,
  height = 5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "153_FigureS1B_PCA_by_batch.pdf"
  ),
  plot = panel_B,
  width = 7,
  height = 5.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "153_FigureS1D_batch_composition.pdf"
  ),
  plot = panel_D,
  width = 7,
  height = 5.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


heatmap_pdf <- file.path(
  figures_dir,
  "153_FigureS1C_sample_distance_heatmap.pdf"
)


grDevices::cairo_pdf(
  heatmap_pdf,
  width = 8,
  height = 7
)


grid::grid.newpage()


grid::grid.draw(
  heatmap_result$gtable
)


grDevices::dev.off()


# =============================================================================
# 26. LIBRARY-SIZE SUMMARY
# =============================================================================

library_summary <- sample_qc %>%
  
  dplyr::group_by(
    condition_label
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    median_million =
      stats::median(
        human_target_counts_million
      ),
    
    q1_million =
      as.numeric(
        stats::quantile(
          human_target_counts_million,
          0.25,
          type = 7
        )
      ),
    
    q3_million =
      as.numeric(
        stats::quantile(
          human_target_counts_million,
          0.75,
          type = 7
        )
      ),
    
    .groups =
      "drop"
  )


# =============================================================================
# 27. FEATURE-UNIVERSE AUDIT TABLE
# =============================================================================

feature_universe_audit <- data.frame(
  
  metric = c(
    
    "Complete raw matrix rows",
    
    "ERCC-prefixed technical rows in complete matrix",
    
    "Biological rows before expression prefilter",
    
    "Original DE prefilter criterion",
    
    "Features passing original DE prefilter",
    
    "ERCC-prefixed features passing original DE prefilter",
    
    "Biological targets passing original DE prefilter",
    
    "Frozen model-comparison rows",
    
    "ERCC-prefixed rows in frozen model comparison",
    
    "Biological rows in frozen model comparison",
    
    "ERCC primary DEG",
    
    "ERCC batch-adjusted DEG",
    
    "ERCC robust-core DEG",
    
    "Biological targets used for Figure S1 VST/PCA"
  ),
  
  value = c(
    
    n_raw_features,
    
    n_ercc_total,
    
    sum(
      !ercc_flag
    ),
    
    "count >=10 in >=3 blood samples",
    
    n_original_prefilter,
    
    n_ercc_passing_prefilter,
    
    n_biological_passing_prefilter,
    
    nrow(
      model_df
    ),
    
    n_model_ercc,
    
    sum(
      !model_ercc_flag
    ),
    
    ercc_primary_deg_n,
    
    ercc_batch_deg_n,
    
    ercc_core_deg_n,
    
    nrow(
      counts_qc
    )
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 28. FIGURE SOURCE-DATA WORKBOOK
# =============================================================================

source_data_file <- file.path(
  tables_dir,
  "153_FigureS1_source_data.xlsx"
)


wb <- openxlsx::createWorkbook()


source_objects <- list(
  
  Feature_universe_audit =
    feature_universe_audit,
  
  ERCC_frozen_model_audit =
    ercc_model_audit,
  
  Sample_QC =
    sample_qc,
  
  Library_summary =
    library_summary,
  
  PCA_coordinates =
    pca_df,
  
  Batch_counts =
    batch_counts,
  
  Batch_condition_table =
    batch_condition_export,
  
  Sample_distance =
    sample_distance_export
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
  
  data_object <- source_objects[[sheet_name]]
  
  
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
# 29. FIGURE LEGEND
# =============================================================================

figure_legend <- paste0(
  
  "Supplementary Figure S1. Discovery-cohort quality control and sequencing ",
  "structure. ",
  
  "(A) Distribution of total human-target counts across 35 sepsis blood ",
  "samples (BP) and 10 healthy-control blood samples (BC). ERCC-prefixed ",
  "technical features were excluded from this descriptive library-size metric. ",
  
  "(B) Principal-component analysis of variance-stabilized expression values ",
  "for 12,393 biological targets satisfying the blood expression prefilter ",
  "(at least 10 raw counts in at least three samples), after exclusion of ",
  "ERCC-prefixed technical features. Points are colored by sequencing batch ",
  "and shaped by study group. ",
  
  "(C) Unsupervised Euclidean sample-distance matrix calculated from the same ",
  "variance-stabilized biological expression matrix, with study-group and ",
  "sequencing-batch annotations. ",
  
  "(D) Distribution of blood samples across sequencing batches according to ",
  "study group. The original blood DESeq2 workflows retained 12,400 features ",
  "after expression prefiltering, comprising 12,393 biological targets and ",
  "seven ERCC-prefixed technical features; none of the retained ERCC features ",
  "was differentially expressed in the primary or batch-adjusted analysis. ",
  "The figure is descriptive and does not imply absence of sequencing-batch ",
  "effects."
)


legend_file <- file.path(
  text_dir,
  "153_FigureS1_legend_EN.txt"
)


writeLines(
  figure_legend,
  legend_file
)


# =============================================================================
# 30. PROPOSED RESULTS 3.1 INSERT
# =============================================================================

results_insert <- paste0(
  
  "Targeted whole-blood transcriptomic profiles were available for all 35 ",
  "patients with sepsis and 10 healthy controls. Quality-control assessment ",
  "included human-target library size, global variance-stabilized expression ",
  "structure, pairwise sample distances, and sequencing-batch composition. ",
  "The original blood DESeq2 expression prefilter retained 12,400 features, ",
  "including 12,393 biological targets and seven ERCC-prefixed technical ",
  "features. None of the retained ERCC-prefixed features was classified as ",
  "differentially expressed in either the primary or batch-adjusted model. ",
  "Accordingly, unsupervised biological QC in Supplementary Fig. S1 was ",
  "performed using the 12,393 retained biological targets after exclusion of ",
  "ERCC-prefixed technical features."
)


results_insert_file <- file.path(
  text_dir,
  "153_proposed_Results_3.1_QC_insert_EN.txt"
)


writeLines(
  results_insert,
  results_insert_file
)


# =============================================================================
# 31. PROPOSED METHODS CORRECTION
# =============================================================================

methods_insert <- paste0(
  
  "For the original blood sepsis-versus-control DESeq2 analyses, the ",
  "prespecified expression filter of at least 10 raw counts in at least ",
  "three blood samples retained 12,400 features. Seven ERCC-prefixed ",
  "technical-control features passed this filter; none was differentially ",
  "expressed in either the primary or batch-adjusted model, and none entered ",
  "the robust-core DEG set. ERCC-prefixed features were excluded from ",
  "downstream biological interpretation and from unsupervised biological ",
  "quality-control visualizations."
)


methods_insert_file <- file.path(
  text_dir,
  "153_proposed_Methods_ERCC_clarification_EN.txt"
)


writeLines(
  methods_insert,
  methods_insert_file
)


# =============================================================================
# 32. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "153_sessionInfo.txt"
  )
)


# =============================================================================
# 33. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 153 FINAL v2 completed successfully.\n")
cat("====================================================================\n\n")


cat("DISCOVERY BLOOD DATASET\n")
cat("-----------------------\n")


cat(
  "Blood samples = ",
  ncol(
    counts_blood
  ),
  "\n",
  sep = ""
)


cat(
  "BP = ",
  sum(
    sample_meta$sample_group ==
      "BP"
  ),
  "\n",
  sep = ""
)


cat(
  "BC = ",
  sum(
    sample_meta$sample_group ==
      "BC"
  ),
  "\n",
  sep = ""
)


cat("\nFEATURE-UNIVERSE PROVENANCE\n")
cat("---------------------------\n")


print(
  feature_universe_audit,
  row.names = FALSE
)


cat("\nERCC FEATURES RETAINED IN ORIGINAL DE UNIVERSE\n")
cat("-----------------------------------------------\n")


print(
  ercc_model_audit[
    ,
    intersect(
      c(
        "Gene",
        "log2FC_simple",
        "padj_simple",
        "DEG_simple",
        "log2FC_batch",
        "padj_batch",
        "DEG_batch",
        "core_DEG"
      ),
      names(
        ercc_model_audit
      )
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


cat("\nPCA\n")
cat("---\n")


cat(
  "Biological targets used for PCA = ",
  nrow(
    vst_matrix
  ),
  "\n",
  sep = ""
)


cat(
  "PC1 variance = ",
  sprintf(
    "%.3f",
    pc1_percent
  ),
  "%\n",
  sep = ""
)


cat(
  "PC2 variance = ",
  sprintf(
    "%.3f",
    pc2_percent
  ),
  "%\n",
  sep = ""
)


cat("\nLIBRARY SIZE\n")
cat("------------\n")


print(
  library_summary,
  row.names = FALSE
)


cat("\nSEQUENCING BATCH COMPOSITION\n")
cat("----------------------------\n")


print(
  batch_condition_export,
  row.names = FALSE
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Combined Figure S1 PDF:\n  ",
  normalizePath(
    figure_pdf,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Combined Figure S1 PNG:\n  ",
  normalizePath(
    figure_png,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Combined Figure S1 TIFF:\n  ",
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
  "Methods ERCC clarification:\n  ",
  normalizePath(
    methods_insert_file,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nINTERPRETATION GUARDRAILS\n")
cat("-------------------------\n")


cat(
  "- Original blood DESeq2 universe = 12,400 retained features.\n"
)


cat(
  "- This comprises 12,393 biological targets + 7 ERCC-prefixed technical features.\n"
)


cat(
  "- None of the 7 retained ERCC features is a primary, batch-adjusted, or robust-core DEG.\n"
)


cat(
  "- Figure S1 biological PCA and sample-distance analyses use 12,393 targets after ERCC exclusion.\n"
)


cat(
  "- Script 153 does not re-run or replace the frozen differential-expression results.\n"
)


cat(
  "- Sequencing batch is displayed explicitly; absence of batch effects is not claimed.\n"
)


cat("\nDone.\n")