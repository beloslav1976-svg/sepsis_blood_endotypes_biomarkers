# ==============================================================================
# Script 142b
# External validation of the frozen five-gene host-response signature
# Dataset: GSE185263
# Project: Sepsis_DESeq2
#
# PRIMARY QUESTION
# Does the frozen five-gene myeloid-adaptive host-response score track
# organ-dysfunction severity, measured by continuous 24-h SOFA,
# in an independent whole-blood RNA-seq cohort?
#
# FROZEN PANEL
#   UP:   CD177, HK3, IRAK3
#   DOWN: CARD11, IKZF2
#
# FROZEN SCORE
#   mean[z(CD177), z(HK3), z(IRAK3)] -
#   mean[z(CARD11), z(IKZF2)]
#
# PRIMARY EXTERNAL ENDPOINT
#   Spearman correlation between five-gene score and continuous 24-h SOFA
#   among patients with sepsis.
#
# SECONDARY PREDECLARED ANALYSES
#   1. SOFA >=2 versus SOFA 0-1
#   2. In-hospital mortality: Died versus Survived
#   3. ICU versus Emergency Room
#   4. Sepsis versus healthy controls - contextual only
#   5. Individual frozen genes versus SOFA
#   6. Location-specific score-SOFA correlations
#   7. Age-, sex-, and location-adjusted score-SOFA association
#
# IMPORTANT
#   - design frozen after Script 142a metadata audit
#   - design frozen BEFORE inspecting five-gene expression in GSE185263
#   - no feature selection
#   - no gene substitution
#   - no coefficient refitting
#   - no cutoff optimization
#   - no post hoc score-direction flipping
#   - no endotype reconstruction
#
# NOTE ON SCRIPT 142a
# The exploratory "timepoint suffix" parser in Script 142a incorrectly
# interpreted endings of ordinary sample IDs such as sepnet0186 as timepoints.
# Therefore NO samples are excluded based on that parser.
#
# Interpretation:
# External transcriptomic replication of a frozen host-response score,
# not validation of a calibrated clinical assay.
# ==============================================================================


# ==============================================================================
# 0. GLOBAL SETTINGS
# ==============================================================================

options(stringsAsFactors = FALSE)

set.seed(20260817)

project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))

script_name <- "142b_external_validation_GSE185263.R"

gse_id <- "GSE185263"

run_date <- Sys.time()


if (!dir.exists(project_dir)) {
  stop(
    "Project directory not found: ",
    project_dir
  )
}

setwd(project_dir)


cat("\n")
cat("====================================================================\n")
cat("Running Script 142b\n")
cat("External validation of frozen five-gene signature: GSE185263\n")
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
# 1. INPUT FILES
# ==============================================================================

audit_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "142a_GSE185263_metadata_audit"
)


metadata_file <- file.path(
  audit_dir,
  "tables",
  "142a_GSE185263_metadata_parsed.csv"
)


raw_counts_file <- file.path(
  audit_dir,
  "raw_download",
  "GSE185263_raw_counts.csv.gz"
)


required_files <- c(
  metadata_file,
  raw_counts_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "Missing required Script 142a output(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


cat(
  "Required Script 142a inputs found.\n\n"
)


# ==============================================================================
# 2. OUTPUT DIRECTORIES
# ==============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "142b_external_validation_GSE185263"
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


for (
  dir_path in c(
    output_dir,
    tables_dir,
    figures_dir,
    text_dir,
    logs_dir
  )
) {
  
  dir.create(
    dir_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


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
# 3. PACKAGES
# ==============================================================================

cran_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "ggplot2",
  "pROC",
  "patchwork",
  "openxlsx"
)


missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_cran) > 0) {
  
  install.packages(
    missing_cran
  )
}


if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages(
    "BiocManager"
  )
}


bioc_packages <- c(
  "edgeR",
  "AnnotationDbi",
  "org.Hs.eg.db"
)


missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_bioc) > 0) {
  
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}


suppressPackageStartupMessages({
  
  library(data.table)
  
  library(dplyr)
  
  library(tidyr)
  
  library(stringr)
  
  library(tibble)
  
  library(ggplot2)
  
  library(pROC)
  
  library(patchwork)
  
  library(openxlsx)
  
  library(edgeR)
  
  library(AnnotationDbi)
  
  library(org.Hs.eg.db)
})


cat(
  "Required packages loaded successfully.\n\n"
)


# ==============================================================================
# 4. FROZEN PANEL
# ==============================================================================

up_genes <- c(
  "CD177",
  "HK3",
  "IRAK3"
)


down_genes <- c(
  "CARD11",
  "IKZF2"
)


five_genes <- c(
  up_genes,
  down_genes
)


expected_sofa_direction <- c(
  CD177 = "POSITIVE",
  HK3 = "POSITIVE",
  IRAK3 = "POSITIVE",
  CARD11 = "NEGATIVE",
  IKZF2 = "NEGATIVE"
)


# ==============================================================================
# 5. WRITE PREDECLARED ANALYSIS PLAN
#
# Written before reading five-gene expression from the raw-count matrix.
# ==============================================================================

predeclared_plan <- c(
  
  "GSE185263 EXTERNAL VALIDATION - PREDECLARED ANALYSIS PLAN",
  
  "====================================================================",
  
  "",
  
  paste0(
    "Plan frozen at: ",
    run_date
  ),
  
  "",
  
  "FROZEN PANEL:",
  
  "UP: CD177, HK3, IRAK3",
  
  "DOWN: CARD11, IKZF2",
  
  "",
  
  "FROZEN SCORE:",
  
  paste0(
    "mean[z(CD177), z(HK3), z(IRAK3)] - ",
    "mean[z(CARD11), z(IKZF2)]"
  ),
  
  "",
  
  "NORMALIZATION:",
  
  paste0(
    "Raw counts -> Ensembl-to-symbol mapping -> ",
    "duplicate-symbol collapse -> edgeR TMM -> logCPM."
  ),
  
  "",
  
  "PRIMARY STANDARDIZATION:",
  
  paste0(
    "Each frozen gene is z-standardized using all sepsis samples. ",
    "No clinical outcome is used in score construction."
  ),
  
  "",
  
  "PRIMARY ENDPOINT:",
  
  paste0(
    "Spearman correlation between frozen five-gene score and ",
    "continuous 24-h SOFA among sepsis samples with available SOFA."
  ),
  
  "",
  
  "PRIMARY DIRECTION:",
  
  paste0(
    "Higher five-gene score is expected to associate with higher SOFA."
  ),
  
  "",
  
  "SECONDARY ANALYSES:",
  
  "1. SOFA >=2 versus SOFA 0-1.",
  
  "2. In-hospital mortality: Died versus Survived.",
  
  "3. ICU versus Emergency Room.",
  
  "4. Sepsis versus healthy controls - contextual only.",
  
  "5. Individual frozen genes versus continuous SOFA.",
  
  "6. Collection-location-specific score-SOFA associations.",
  
  paste0(
    "7. Age-, sex-, and collection-location-adjusted ",
    "score-SOFA association."
  ),
  
  "",
  
  "MULTIPLICITY:",
  
  paste0(
    "The score-SOFA correlation is the single primary endpoint."
  ),
  
  paste0(
    "Secondary score-level tests are BH-adjusted as one family."
  ),
  
  paste0(
    "Five gene-level SOFA tests are BH-adjusted as a separate family."
  ),
  
  "",
  
  "FROZEN RULES:",
  
  "- no feature selection",
  
  "- no gene addition or deletion",
  
  "- no gene substitution",
  
  "- no coefficient refitting",
  
  "- no cutoff optimization",
  
  "- no post hoc score-direction reversal",
  
  "- no endotype reconstruction",
  
  "",
  
  "SAMPLE-HANDLING NOTE:",
  
  paste0(
    "The exploratory timepoint parser from Script 142a is not used ",
    "because ordinary sample IDs such as sepnet0186 were incorrectly ",
    "interpreted as timepoint suffixes."
  ),
  
  "",
  
  "INTERPRETATION:",
  
  paste0(
    "External transcriptomic replication of a frozen molecular ",
    "host-response score, not validation of a calibrated clinical assay."
  )
)


plan_file <- file.path(
  text_dir,
  "142b_PREDECLARED_external_validation_plan.txt"
)


writeLines(
  predeclared_plan,
  con = plan_file,
  useBytes = TRUE
)


cat(
  "Predeclared validation plan written before expression analysis.\n\n"
)


# ==============================================================================
# 6. HELPER FUNCTIONS
# ==============================================================================

normalize_key <- function(x) {
  
  x <- as.character(x)
  
  x <- tolower(x)
  
  x <- gsub(
    "[^a-z0-9]",
    "",
    x
  )
  
  return(x)
}


safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


standardize_by_reference <- function(
    mat,
    reference_samples
) {
  
  mat <- as.matrix(mat)
  
  storage.mode(mat) <- "numeric"
  
  
  if (!all(
    reference_samples %in%
    colnames(mat)
  )) {
    
    stop(
      "Some reference samples are absent from expression matrix."
    )
  }
  
  
  reference_matrix <- mat[
    ,
    reference_samples,
    drop = FALSE
  ]
  
  
  reference_mean <- rowMeans(
    reference_matrix,
    na.rm = TRUE
  )
  
  
  reference_sd <- apply(
    reference_matrix,
    1,
    stats::sd,
    na.rm = TRUE
  )
  
  
  bad <- !is.finite(reference_sd) |
    reference_sd == 0
  
  
  if (any(bad)) {
    
    stop(
      paste0(
        "Cannot standardize gene(s) with zero/non-finite SD: ",
        paste(
          rownames(mat)[bad],
          collapse = ", "
        )
      )
    )
  }
  
  
  out <- sweep(
    mat,
    1,
    reference_mean,
    FUN = "-"
  )
  
  
  out <- sweep(
    out,
    1,
    reference_sd,
    FUN = "/"
  )
  
  
  return(out)
}


calculate_five_gene_score <- function(
    z_matrix
) {
  
  if (!all(
    five_genes %in%
    rownames(z_matrix)
  )) {
    
    stop(
      "Frozen five-gene panel incomplete in standardized matrix."
    )
  }
  
  
  up_component <- colMeans(
    z_matrix[
      up_genes,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  
  down_component <- colMeans(
    z_matrix[
      down_genes,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  
  return(
    up_component -
      down_component
  )
}


safe_spearman <- function(
    x,
    y
) {
  
  keep <- is.finite(x) &
    is.finite(y)
  
  
  x <- x[keep]
  
  y <- y[keep]
  
  
  if (length(x) < 5) {
    
    return(
      tibble::tibble(
        n = length(x),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  test <- suppressWarnings(
    stats::cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  
  return(
    tibble::tibble(
      n = length(x),
      rho = unname(
        as.numeric(
          test$estimate
        )
      ),
      p_value = test$p.value
    )
  )
}


safe_wilcox <- function(
    case_values,
    control_values
) {
  
  case_values <- case_values[
    is.finite(case_values)
  ]
  
  
  control_values <- control_values[
    is.finite(control_values)
  ]
  
  
  if (
    length(case_values) < 2 ||
    length(control_values) < 2
  ) {
    
    return(
      tibble::tibble(
        W = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  test <- stats::wilcox.test(
    case_values,
    control_values,
    exact = FALSE,
    paired = FALSE
  )
  
  
  return(
    tibble::tibble(
      W = unname(
        as.numeric(
          test$statistic
        )
      ),
      p_value = test$p.value
    )
  )
}


fixed_direction_auc <- function(
    case_values,
    control_values
) {
  
  case_values <- case_values[
    is.finite(case_values)
  ]
  
  
  control_values <- control_values[
    is.finite(control_values)
  ]
  
  
  if (
    length(case_values) < 3 ||
    length(control_values) < 3
  ) {
    
    return(
      tibble::tibble(
        AUC = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
    )
  }
  
  
  response <- factor(
    c(
      rep(
        "control",
        length(control_values)
      ),
      rep(
        "case",
        length(case_values)
      )
    ),
    levels = c(
      "control",
      "case"
    )
  )
  
  
  predictor <- c(
    control_values,
    case_values
  )
  
  
  roc_object <- pROC::roc(
    response = response,
    predictor = predictor,
    levels = c(
      "control",
      "case"
    ),
    direction = "<",
    quiet = TRUE
  )
  
  
  auc_ci <- suppressWarnings(
    pROC::ci.auc(
      roc_object,
      method = "delong"
    )
  )
  
  
  return(
    tibble::tibble(
      AUC = as.numeric(
        pROC::auc(
          roc_object
        )
      ),
      CI_low = as.numeric(
        auc_ci[1]
      ),
      CI_high = as.numeric(
        auc_ci[3]
      )
    )
  )
}


two_group_score_test <- function(
    df,
    group_column,
    case_value,
    control_value,
    analysis_name,
    score_column = "five_gene_score"
) {
  
  case_values <- df %>%
    
    dplyr::filter(
      .data[[group_column]] ==
        case_value
    ) %>%
    
    dplyr::pull(
      dplyr::all_of(
        score_column
      )
    )
  
  
  control_values <- df %>%
    
    dplyr::filter(
      .data[[group_column]] ==
        control_value
    ) %>%
    
    dplyr::pull(
      dplyr::all_of(
        score_column
      )
    )
  
  
  wt <- safe_wilcox(
    case_values,
    control_values
  )
  
  
  auc_result <- fixed_direction_auc(
    case_values,
    control_values
  )
  
  
  return(
    tibble::tibble(
      
      analysis =
        analysis_name,
      
      score =
        score_column,
      
      case =
        case_value,
      
      control =
        control_value,
      
      n_case =
        sum(
          is.finite(
            case_values
          )
        ),
      
      n_control =
        sum(
          is.finite(
            control_values
          )
        ),
      
      median_case =
        stats::median(
          case_values,
          na.rm = TRUE
        ),
      
      median_control =
        stats::median(
          control_values,
          na.rm = TRUE
        ),
      
      median_difference =
        stats::median(
          case_values,
          na.rm = TRUE
        ) -
        stats::median(
          control_values,
          na.rm = TRUE
        ),
      
      W =
        wt$W,
      
      p_value =
        wt$p_value,
      
      AUC =
        auc_result$AUC,
      
      CI_low =
        auc_result$CI_low,
      
      CI_high =
        auc_result$CI_high
    )
  )
}


save_plot_all_formats <- function(
    plot_object,
    filename_base,
    width,
    height
) {
  
  png_file <- file.path(
    figures_dir,
    paste0(
      filename_base,
      ".png"
    )
  )
  
  
  tiff_file <- file.path(
    figures_dir,
    paste0(
      filename_base,
      ".tiff"
    )
  )
  
  
  pdf_file <- file.path(
    figures_dir,
    paste0(
      filename_base,
      ".pdf"
    )
  )
  
  
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  
  ggplot2::ggsave(
    filename = tiff_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  
  
  if (isTRUE(
    capabilities("cairo")
  )) {
    
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot_object,
      width = width,
      height = height,
      device = grDevices::cairo_pdf,
      bg = "white"
    )
    
  } else {
    
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot_object,
      width = width,
      height = height,
      bg = "white"
    )
  }
  
  
  invisible(
    c(
      png_file,
      tiff_file,
      pdf_file
    )
  )
}


# ==============================================================================
# 7. READ AND STANDARDIZE METADATA
# ==============================================================================

metadata <- utils::read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_metadata_columns <- c(
  "geo_accession",
  "sample_title",
  "age",
  "sex",
  "collection_location",
  "collection_site",
  "disease_state",
  "in_hospital_mortality",
  "sofa_24h_post_admisssion"
)


missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(metadata)
)


if (length(
  missing_metadata_columns
) > 0) {
  
  stop(
    paste0(
      "Missing required metadata column(s): ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    )
  )
}


metadata <- metadata %>%
  
  dplyr::mutate(
    
    geo_accession =
      as.character(
        geo_accession
      ),
    
    sample_title =
      as.character(
        sample_title
      ),
    
    sample_key =
      normalize_key(
        sample_title
      ),
    
    disease_state =
      tolower(
        trimws(
          as.character(
            disease_state
          )
        )
      ),
    
    age_numeric =
      safe_numeric(
        age
      ),
    
    sex_standardized =
      toupper(
        trimws(
          as.character(
            sex
          )
        )
      ),
    
    collection_location =
      tolower(
        trimws(
          as.character(
            collection_location
          )
        )
      ),
    
    collection_site =
      trimws(
        as.character(
          collection_site
        )
      ),
    
    mortality =
      trimws(
        as.character(
          in_hospital_mortality
        )
      ),
    
    sofa =
      safe_numeric(
        sofa_24h_post_admisssion
      )
  )


cat(
  "Metadata disease-state distribution:\n"
)


print(
  table(
    metadata$disease_state,
    useNA = "ifany"
  )
)


cat("\n")


if (nrow(metadata) != 392) {
  
  stop(
    paste0(
      "Expected 392 metadata rows; observed ",
      nrow(metadata)
    )
  )
}


if (
  sum(
    metadata$disease_state ==
    "sepsis",
    na.rm = TRUE
  ) != 348
) {
  
  warning(
    "Sepsis sample count differs from expected n=348."
  )
}


# ==============================================================================
# 8. READ RAW COUNT MATRIX
# ==============================================================================

cat(
  "Reading raw-count matrix...\n"
)


counts_raw <- data.table::fread(
  raw_counts_file,
  data.table = FALSE,
  check.names = FALSE
)


if (
  !"ensembl_gene_id" %in%
  names(counts_raw)
) {
  
  stop(
    "Column 'ensembl_gene_id' not found in raw-count matrix."
  )
}


count_sample_columns <- setdiff(
  names(counts_raw),
  "ensembl_gene_id"
)


cat(
  "Raw count matrix: ",
  nrow(counts_raw),
  " genes x ",
  length(count_sample_columns),
  " samples\n\n",
  sep = ""
)


if (
  length(
    count_sample_columns
  ) != 392
) {
  
  stop(
    paste0(
      "Expected 392 sample columns; observed ",
      length(
        count_sample_columns
      )
    )
  )
}


# ==============================================================================
# 9. MAP COUNT-MATRIX COLUMNS TO GEO METADATA
# ==============================================================================

count_mapping <- tibble::tibble(
  
  count_column =
    count_sample_columns,
  
  sample_key =
    normalize_key(
      count_sample_columns
    )
) %>%
  
  dplyr::left_join(
    
    metadata %>%
      
      dplyr::select(
        geo_accession,
        sample_title,
        sample_key
      ),
    
    by =
      "sample_key"
  )


if (
  any(
    is.na(
      count_mapping$geo_accession
    )
  )
) {
  
  failed <- count_mapping %>%
    
    dplyr::filter(
      is.na(
        geo_accession
      )
    )
  
  
  print(
    failed,
    n = Inf
  )
  
  
  stop(
    "Some count-matrix sample columns could not be matched to metadata."
  )
}


if (
  anyDuplicated(
    count_mapping$geo_accession
  ) > 0
) {
  
  stop(
    "Duplicate GEO accessions after count-column mapping."
  )
}


metadata_counts_order <- metadata[
  match(
    count_mapping$geo_accession,
    metadata$geo_accession
  ),
  ,
  drop = FALSE
]


if (
  any(
    is.na(
      metadata_counts_order$geo_accession
    )
  )
) {
  
  stop(
    "Failed to reorder metadata according to count matrix."
  )
}


cat(
  "Count-matrix sample mapping: 392/392 successful.\n\n"
)


# ==============================================================================
# 10. ENSEMBL -> GENE SYMBOL MAPPING
# ==============================================================================

ensembl_raw <- as.character(
  counts_raw$ensembl_gene_id
)


ensembl_clean <- sub(
  "\\.[0-9]+$",
  "",
  ensembl_raw
)


cat(
  "Mapping Ensembl identifiers to HGNC gene symbols...\n"
)


symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(
    ensembl_clean
  ),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)


gene_symbol <- unname(
  symbol_map[
    ensembl_clean
  ]
)


gene_symbol <- toupper(
  trimws(
    as.character(
      gene_symbol
    )
  )
)


mapping_keep <- !is.na(
  gene_symbol
) &
  gene_symbol != ""


cat(
  "Mapped Ensembl rows: ",
  sum(mapping_keep),
  " / ",
  length(mapping_keep),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 11. BUILD GENE-SYMBOL RAW COUNT MATRIX
# ==============================================================================

counts_matrix <- as.matrix(
  counts_raw[
    ,
    count_sample_columns,
    drop = FALSE
  ]
)


storage.mode(
  counts_matrix
) <- "numeric"


counts_matrix <- counts_matrix[
  mapping_keep,
  ,
  drop = FALSE
]


mapped_symbols <- gene_symbol[
  mapping_keep
]


# Duplicate HGNC symbols are collapsed by SUM before normalization.

counts_symbol <- rowsum(
  counts_matrix,
  group = mapped_symbols,
  reorder = FALSE
)


nonzero_rows <- rowSums(
  counts_symbol,
  na.rm = TRUE
) > 0


counts_symbol <- counts_symbol[
  nonzero_rows,
  ,
  drop = FALSE
]


cat(
  "Collapsed count matrix: ",
  nrow(counts_symbol),
  " gene symbols x ",
  ncol(counts_symbol),
  " samples\n\n",
  sep = ""
)


# ==============================================================================
# 12. FROZEN FIVE-GENE COVERAGE
# ==============================================================================

gene_coverage <- tibble::tibble(
  
  gene =
    five_genes,
  
  expected_sofa_direction =
    unname(
      expected_sofa_direction[
        five_genes
      ]
    ),
  
  present =
    five_genes %in%
    rownames(counts_symbol)
)


cat(
  "Frozen five-gene coverage:\n"
)


print(
  gene_coverage,
  n = Inf
)


cat("\n")


if (!all(
  gene_coverage$present
)) {
  
  stop(
    paste0(
      "Frozen panel incomplete. Missing: ",
      paste(
        gene_coverage$gene[
          !gene_coverage$present
        ],
        collapse = ", "
      )
    )
  )
}


cat(
  "All five frozen genes detected.\n\n"
)


# ==============================================================================
# 13. edgeR TMM NORMALIZATION
# ==============================================================================

cat(
  "Performing edgeR TMM normalization...\n"
)


dge <- edgeR::DGEList(
  counts = counts_symbol
)


dge <- edgeR::calcNormFactors(
  dge,
  method = "TMM"
)


logcpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 2
)


five_gene_logcpm <- logcpm[
  five_genes,
  ,
  drop = FALSE
]


colnames(
  five_gene_logcpm
) <- metadata_counts_order$geo_accession


cat(
  "TMM logCPM calculated successfully.\n\n"
)


# ==============================================================================
# 14. DEFINE SEPSIS AND HEALTHY COHORTS
# ==============================================================================

sepsis_meta <- metadata_counts_order %>%
  
  dplyr::filter(
    disease_state ==
      "sepsis"
  )


healthy_meta <- metadata_counts_order %>%
  
  dplyr::filter(
    disease_state ==
      "healthy"
  )


sepsis_samples <- sepsis_meta$geo_accession

healthy_samples <- healthy_meta$geo_accession


cat(
  "GSE185263 cohort:\n"
)


cat(
  "Sepsis = ",
  length(
    sepsis_samples
  ),
  "\n",
  sep = ""
)


cat(
  "Healthy = ",
  length(
    healthy_samples
  ),
  "\n\n",
  sep = ""
)


if (
  length(
    sepsis_samples
  ) != 348
) {
  
  warning(
    "Expected 348 sepsis samples."
  )
}


if (
  length(
    healthy_samples
  ) != 44
) {
  
  warning(
    "Expected 44 healthy samples."
  )
}


# ==============================================================================
# 15. FROZEN SCORE STANDARDIZATION
#
# PRIMARY:
# z-standardization using all sepsis samples.
#
# SENSITIVITY:
# all-sample reference and healthy-control reference.
# ==============================================================================

z_sepsis_reference <- standardize_by_reference(
  five_gene_logcpm,
  reference_samples = sepsis_samples
)


score_sepsis_reference <- calculate_five_gene_score(
  z_sepsis_reference
)


z_all_reference <- standardize_by_reference(
  five_gene_logcpm,
  reference_samples = colnames(
    five_gene_logcpm
  )
)


score_all_reference <- calculate_five_gene_score(
  z_all_reference
)


z_healthy_reference <- standardize_by_reference(
  five_gene_logcpm,
  reference_samples = healthy_samples
)


score_healthy_reference <- calculate_five_gene_score(
  z_healthy_reference
)


cat(
  "Frozen five-gene scores calculated.\n\n"
)


# ==============================================================================
# 16. SAMPLE-LEVEL SCORE TABLE
# ==============================================================================

scores <- metadata_counts_order %>%
  
  dplyr::select(
    geo_accession,
    sample_title,
    disease_state,
    age_numeric,
    sex_standardized,
    collection_location,
    collection_site,
    mortality,
    sofa
  ) %>%
  
  dplyr::mutate(
    
    five_gene_score =
      unname(
        score_sepsis_reference[
          geo_accession
        ]
      ),
    
    five_gene_score_all_reference =
      unname(
        score_all_reference[
          geo_accession
        ]
      ),
    
    five_gene_score_healthy_reference =
      unname(
        score_healthy_reference[
          geo_accession
        ]
      ),
    
    sofa_group =
      dplyr::case_when(
        
        !is.finite(
          sofa
        ) ~
          NA_character_,
        
        sofa >= 2 ~
          "SOFA_ge2",
        
        sofa < 2 ~
          "SOFA_0_1",
        
        TRUE ~
          NA_character_
      )
  )


for (
  gene in five_genes
) {
  
  scores[[paste0(gene, "_logCPM")]] <-
    five_gene_logcpm[
      gene,
      scores$geo_accession
    ]
  
  
  scores[[paste0(gene, "_z_sepsis")]] <-
    z_sepsis_reference[
      gene,
      scores$geo_accession
    ]
}


sepsis_scores <- scores %>%
  
  dplyr::filter(
    disease_state ==
      "sepsis"
  )


sofa_available_n <- sum(
  is.finite(
    sepsis_scores$sofa
  )
)


mortality_available_n <- sum(
  sepsis_scores$mortality %in%
    c(
      "Died",
      "Survived"
    )
)


cat(
  "SEPSIS VALIDATION COHORT:\n"
)


cat(
  "Total sepsis samples = ",
  nrow(
    sepsis_scores
  ),
  "\n",
  sep = ""
)


cat(
  "SOFA available = ",
  sofa_available_n,
  "\n",
  sep = ""
)


cat(
  "Mortality available = ",
  mortality_available_n,
  "\n\n",
  sep = ""
)


# ==============================================================================
# 17. PRIMARY ENDPOINT
# FIVE-GENE SCORE vs CONTINUOUS SOFA
# ==============================================================================

primary_sofa <- safe_spearman(
  sepsis_scores$five_gene_score,
  sepsis_scores$sofa
) %>%
  
  dplyr::mutate(
    
    analysis =
      "PRIMARY_five_gene_score_vs_SOFA",
    
    expected_direction =
      "POSITIVE"
  ) %>%
  
  dplyr::select(
    analysis,
    expected_direction,
    dplyr::everything()
  )


cat(
  "PRIMARY EXTERNAL ENDPOINT:\n"
)


print(
  primary_sofa,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 18. SECONDARY SCORE-LEVEL ANALYSES
# ==============================================================================

secondary_sofa <- two_group_score_test(
  
  df =
    sepsis_scores %>%
    
    dplyr::filter(
      !is.na(
        sofa_group
      )
    ),
  
  group_column =
    "sofa_group",
  
  case_value =
    "SOFA_ge2",
  
  control_value =
    "SOFA_0_1",
  
  analysis_name =
    "SOFA_ge2_vs_SOFA_0_1",
  
  score_column =
    "five_gene_score"
)


secondary_mortality <- two_group_score_test(
  
  df =
    sepsis_scores %>%
    
    dplyr::filter(
      mortality %in%
        c(
          "Died",
          "Survived"
        )
    ),
  
  group_column =
    "mortality",
  
  case_value =
    "Died",
  
  control_value =
    "Survived",
  
  analysis_name =
    "Died_vs_Survived",
  
  score_column =
    "five_gene_score"
)


secondary_site <- two_group_score_test(
  
  df =
    sepsis_scores %>%
    
    dplyr::filter(
      collection_site %in%
        c(
          "ICU",
          "Emergency Room"
        )
    ),
  
  group_column =
    "collection_site",
  
  case_value =
    "ICU",
  
  control_value =
    "Emergency Room",
  
  analysis_name =
    "ICU_vs_Emergency_Room",
  
  score_column =
    "five_gene_score"
)


disease_context <- scores %>%
  
  dplyr::filter(
    disease_state %in%
      c(
        "sepsis",
        "healthy"
      )
  )


secondary_disease <- two_group_score_test(
  
  df =
    disease_context,
  
  group_column =
    "disease_state",
  
  case_value =
    "sepsis",
  
  control_value =
    "healthy",
  
  analysis_name =
    "Sepsis_vs_healthy_contextual",
  
  score_column =
    "five_gene_score_all_reference"
)


secondary_results <- dplyr::bind_rows(
  secondary_sofa,
  secondary_mortality,
  secondary_site,
  secondary_disease
) %>%
  
  dplyr::mutate(
    
    BH_secondary =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )


cat(
  "SECONDARY SCORE-LEVEL ANALYSES:\n"
)


print(
  secondary_results,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 19. INDIVIDUAL FROZEN GENES vs SOFA
# ==============================================================================

gene_sofa_results <- dplyr::bind_rows(
  
  lapply(
    five_genes,
    function(gene) {
      
      gene_column <- paste0(
        gene,
        "_z_sepsis"
      )
      
      
      gene_values <- sepsis_scores[[gene_column]]
      
      
      result <- safe_spearman(
        gene_values,
        sepsis_scores$sofa
      )
      
      
      expected <- expected_sofa_direction[[gene]]
      
      
      observed <- dplyr::case_when(
        
        is.na(
          result$rho
        ) ~
          NA_character_,
        
        result$rho > 0 ~
          "POSITIVE",
        
        result$rho < 0 ~
          "NEGATIVE",
        
        TRUE ~
          "ZERO"
      )
      
      
      result %>%
        
        dplyr::mutate(
          
          gene =
            gene,
          
          expected_direction =
            expected,
          
          observed_direction =
            observed,
          
          direction_concordant =
            expected ==
            observed
        ) %>%
        
        dplyr::select(
          gene,
          expected_direction,
          observed_direction,
          direction_concordant,
          dplyr::everything()
        )
    }
  )
) %>%
  
  dplyr::mutate(
    
    BH_five_genes =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )


cat(
  "INDIVIDUAL FIVE-GENE SOFA ASSOCIATIONS:\n"
)


print(
  gene_sofa_results,
  n = Inf,
  width = Inf
)


cat("\n")


gene_concordant_n <- sum(
  gene_sofa_results$direction_concordant,
  na.rm = TRUE
)


# ==============================================================================
# 20. AGE / SEX / LOCATION-ADJUSTED SOFA MODEL
#
# Outcome:
#   five_gene_score
#
# Predictor of interest:
#   SOFA
#
# Covariates:
#   age
#   sex
#   collection location
# ==============================================================================

adjusted_data <- sepsis_scores %>%
  
  dplyr::filter(
    
    is.finite(
      five_gene_score
    ),
    
    is.finite(
      sofa
    ),
    
    is.finite(
      age_numeric
    ),
    
    sex_standardized %in%
      c(
        "M",
        "F"
      ),
    
    !is.na(
      collection_location
    ),
    
    collection_location != ""
  ) %>%
  
  dplyr::mutate(
    
    sex_standardized =
      factor(
        sex_standardized
      ),
    
    collection_location =
      factor(
        collection_location
      )
  )


adjusted_model <- stats::lm(
  
  five_gene_score ~
    sofa +
    age_numeric +
    sex_standardized +
    collection_location,
  
  data =
    adjusted_data
)


adjusted_coefficients_matrix <- summary(
  adjusted_model
)$coefficients


adjusted_coefficients <- tibble::as_tibble(
  adjusted_coefficients_matrix,
  rownames = "term"
)


names(
  adjusted_coefficients
) <- c(
  "term",
  "estimate",
  "SE",
  "t_value",
  "p_value"
)


adjusted_sofa_row <- adjusted_coefficients %>%
  
  dplyr::filter(
    term ==
      "sofa"
  )


if (
  nrow(
    adjusted_sofa_row
  ) != 1
) {
  
  stop(
    "SOFA coefficient was not uniquely identified in adjusted model."
  )
}


adjusted_model_summary <- tibble::tibble(
  
  metric = c(
    "n",
    "R_squared",
    "adjusted_R_squared",
    "residual_SD"
  ),
  
  value = c(
    nrow(
      adjusted_data
    ),
    
    summary(
      adjusted_model
    )$r.squared,
    
    summary(
      adjusted_model
    )$adj.r.squared,
    
    summary(
      adjusted_model
    )$sigma
  )
)


cat(
  "AGE/SEX/LOCATION-ADJUSTED SOFA ASSOCIATION:\n"
)


print(
  adjusted_sofa_row,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "Adjusted model summary:\n"
)


print(
  adjusted_model_summary,
  n = Inf
)


cat("\n")


# ==============================================================================
# 21. LOCATION-SPECIFIC SCORE-SOFA ASSOCIATIONS
# ==============================================================================

location_counts <- sepsis_scores %>%
  
  dplyr::filter(
    is.finite(
      sofa
    )
  ) %>%
  
  dplyr::count(
    collection_location,
    name = "n_with_SOFA"
  )


eligible_locations <- location_counts %>%
  
  dplyr::filter(
    n_with_SOFA >= 10
  ) %>%
  
  dplyr::pull(
    collection_location
  )


location_sofa_results <- dplyr::bind_rows(
  
  lapply(
    eligible_locations,
    function(location_value) {
      
      location_data <- sepsis_scores %>%
        
        dplyr::filter(
          collection_location ==
            location_value
        )
      
      
      result <- safe_spearman(
        location_data$five_gene_score,
        location_data$sofa
      )
      
      
      result %>%
        
        dplyr::mutate(
          
          collection_location =
            location_value,
          
          direction_concordant =
            rho > 0
        ) %>%
        
        dplyr::select(
          collection_location,
          direction_concordant,
          dplyr::everything()
        )
    }
  )
)


if (
  nrow(
    location_sofa_results
  ) > 0
) {
  
  location_sofa_results <- location_sofa_results %>%
    
    dplyr::mutate(
      
      BH_location =
        stats::p.adjust(
          p_value,
          method = "BH"
        )
    )
}


cat(
  "LOCATION-STRATIFIED SCORE-SOFA ASSOCIATIONS:\n"
)


print(
  location_sofa_results,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 22. CROSS-LOCATION FISHER-Z SYNTHESIS
#
# Descriptive fixed-effect synthesis of correlations.
# ==============================================================================

meta_locations <- location_sofa_results %>%
  
  dplyr::filter(
    
    is.finite(
      rho
    ),
    
    n > 3,
    
    abs(
      rho
    ) < 1
  ) %>%
  
  dplyr::mutate(
    
    fisher_z =
      atanh(
        rho
      ),
    
    weight =
      n -
      3
  )


if (
  nrow(
    meta_locations
  ) >= 2
) {
  
  pooled_z <- sum(
    meta_locations$weight *
      meta_locations$fisher_z
  ) /
    sum(
      meta_locations$weight
    )
  
  
  pooled_se <- sqrt(
    1 /
      sum(
        meta_locations$weight
      )
  )
  
  
  pooled_z_stat <- pooled_z /
    pooled_se
  
  
  pooled_p <- 2 *
    stats::pnorm(
      -abs(
        pooled_z_stat
      )
    )
  
  
  pooled_ci_low_z <- pooled_z -
    1.96 *
    pooled_se
  
  
  pooled_ci_high_z <- pooled_z +
    1.96 *
    pooled_se
  
  
  pooled_location_result <- tibble::tibble(
    
    locations =
      nrow(
        meta_locations
      ),
    
    total_n =
      sum(
        meta_locations$n
      ),
    
    pooled_rho =
      tanh(
        pooled_z
      ),
    
    CI_low =
      tanh(
        pooled_ci_low_z
      ),
    
    CI_high =
      tanh(
        pooled_ci_high_z
      ),
    
    p_value =
      pooled_p
  )
  
} else {
  
  pooled_location_result <- tibble::tibble(
    
    locations =
      NA_integer_,
    
    total_n =
      NA_integer_,
    
    pooled_rho =
      NA_real_,
    
    CI_low =
      NA_real_,
    
    CI_high =
      NA_real_,
    
    p_value =
      NA_real_
  )
}


cat(
  "CROSS-LOCATION FIXED-EFFECT SUMMARY:\n"
)


print(
  pooled_location_result,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 23. SCORE SCALING SENSITIVITY
# ==============================================================================

sepsis_scaling <- scores %>%
  
  dplyr::filter(
    disease_state ==
      "sepsis"
  )


scaling_all <- safe_spearman(
  sepsis_scaling$five_gene_score,
  sepsis_scaling$five_gene_score_all_reference
) %>%
  
  dplyr::mutate(
    
    comparison =
      "sepsis_reference_vs_all_sample_reference"
  )


scaling_healthy <- safe_spearman(
  sepsis_scaling$five_gene_score,
  sepsis_scaling$five_gene_score_healthy_reference
) %>%
  
  dplyr::mutate(
    
    comparison =
      "sepsis_reference_vs_healthy_reference"
  )


scaling_results <- dplyr::bind_rows(
  scaling_all,
  scaling_healthy
) %>%
  
  dplyr::select(
    comparison,
    dplyr::everything()
  )


cat(
  "SCALING SENSITIVITY:\n"
)


print(
  scaling_results,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 24. BASIC SUMMARIES
# ==============================================================================

sofa_summary <- sepsis_scores %>%
  
  dplyr::filter(
    is.finite(
      sofa
    )
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    min =
      min(
        sofa
      ),
    
    q1 =
      stats::quantile(
        sofa,
        0.25
      ),
    
    median =
      stats::median(
        sofa
      ),
    
    mean =
      mean(
        sofa
      ),
    
    q3 =
      stats::quantile(
        sofa,
        0.75
      ),
    
    max =
      max(
        sofa
      )
  )


sofa_group_summary <- sepsis_scores %>%
  
  dplyr::filter(
    !is.na(
      sofa_group
    )
  ) %>%
  
  dplyr::group_by(
    sofa_group
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    score_median =
      stats::median(
        five_gene_score
      ),
    
    score_q1 =
      stats::quantile(
        five_gene_score,
        0.25
      ),
    
    score_q3 =
      stats::quantile(
        five_gene_score,
        0.75
      ),
    
    .groups =
      "drop"
  )


mortality_summary <- sepsis_scores %>%
  
  dplyr::filter(
    mortality %in%
      c(
        "Died",
        "Survived"
      )
  ) %>%
  
  dplyr::group_by(
    mortality
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    score_median =
      stats::median(
        five_gene_score
      ),
    
    score_q1 =
      stats::quantile(
        five_gene_score,
        0.25
      ),
    
    score_q3 =
      stats::quantile(
        five_gene_score,
        0.75
      ),
    
    sofa_median =
      stats::median(
        sofa,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


site_summary <- sepsis_scores %>%
  
  dplyr::count(
    collection_site,
    name = "n"
  )


location_summary <- sepsis_scores %>%
  
  dplyr::count(
    collection_location,
    name = "n"
  )


# ==============================================================================
# 25. PUBLICATION COLORS
# ==============================================================================

sofa_colors <- c(
  "SOFA_0_1" = "#56B4E9",
  "SOFA_ge2" = "#D55E00"
)


mortality_colors <- c(
  "Survived" = "#0072B2",
  "Died" = "#D55E00"
)


disease_colors <- c(
  "healthy" = "#56B4E9",
  "sepsis" = "#D55E00"
)


direction_colors <- c(
  "NEGATIVE" = "#0072B2",
  "POSITIVE" = "#D55E00"
)


# ==============================================================================
# 26. FIGURE A - PRIMARY SCORE vs SOFA
# ==============================================================================

pA_data <- sepsis_scores %>%
  
  dplyr::filter(
    is.finite(
      sofa
    )
  )


pA <- ggplot2::ggplot(
  
  pA_data,
  
  ggplot2::aes(
    x =
      sofa,
    y =
      five_gene_score
  )
) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.50,
    color = "#D55E00"
  ) +
  
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey30",
    fill = "grey82",
    linewidth = 0.9
  ) +
  
  ggplot2::labs(
    
    title =
      "Frozen five-gene score and organ-dysfunction severity",
    
    subtitle =
      paste0(
        "GSE185263 sepsis cohort; Spearman rho = ",
        sprintf(
          "%.3f",
          primary_sofa$rho
        ),
        "; p = ",
        format.pval(
          primary_sofa$p_value,
          digits = 3
        ),
        "; n = ",
        primary_sofa$n
      ),
    
    x =
      "SOFA score, 24 h after admission",
    
    y =
      "Five-gene host-response score"
  ) +
  
  ggplot2::theme_bw(
    base_size = 14
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 27. FIGURE B - SCORE BY SOFA >=2
# ==============================================================================

pB_data <- sepsis_scores %>%
  
  dplyr::filter(
    !is.na(
      sofa_group
    )
  )


pB <- ggplot2::ggplot(
  
  pB_data,
  
  ggplot2::aes(
    x =
      sofa_group,
    y =
      five_gene_score,
    fill =
      sofa_group,
    color =
      sofa_group
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.62,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.45
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      sofa_colors
  ) +
  
  ggplot2::scale_color_manual(
    values =
      sofa_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Five-gene score by SOFA threshold",
    
    subtitle =
      paste0(
        "SOFA >=2 vs 0-1; p = ",
        format.pval(
          secondary_sofa$p_value,
          digits = 3
        ),
        "; AUC = ",
        sprintf(
          "%.3f",
          secondary_sofa$AUC
        )
      ),
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  ) +
  
  ggplot2::theme_bw(
    base_size = 14
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 28. FIGURE C - MORTALITY
# ==============================================================================

pC_data <- sepsis_scores %>%
  
  dplyr::filter(
    mortality %in%
      c(
        "Survived",
        "Died"
      )
  )


pC <- ggplot2::ggplot(
  
  pC_data,
  
  ggplot2::aes(
    x =
      mortality,
    y =
      five_gene_score,
    fill =
      mortality,
    color =
      mortality
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.62,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.45
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      mortality_colors
  ) +
  
  ggplot2::scale_color_manual(
    values =
      mortality_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Five-gene score and in-hospital mortality",
    
    subtitle =
      paste0(
        "Died vs Survived; p = ",
        format.pval(
          secondary_mortality$p_value,
          digits = 3
        ),
        "; AUC = ",
        sprintf(
          "%.3f",
          secondary_mortality$AUC
        )
      ),
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  ) +
  
  ggplot2::theme_bw(
    base_size = 14
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 29. FIGURE D - INDIVIDUAL GENE-SOFA CORRELATIONS
# ==============================================================================

gene_plot_data <- gene_sofa_results %>%
  
  dplyr::mutate(
    
    gene =
      factor(
        gene,
        levels =
          rev(
            five_genes
          )
      ),
    
    expected_direction =
      factor(
        expected_direction,
        levels = c(
          "NEGATIVE",
          "POSITIVE"
        )
      )
  )


pD <- ggplot2::ggplot(
  
  gene_plot_data,
  
  ggplot2::aes(
    x =
      rho,
    y =
      gene,
    color =
      expected_direction
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend =
        rho,
      y =
        gene,
      yend =
        gene
    ),
    
    linewidth = 1
  ) +
  
  ggplot2::geom_point(
    size = 4
  ) +
  
  ggplot2::scale_color_manual(
    values =
      direction_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Frozen component genes and SOFA",
    
    subtitle =
      "Spearman correlations within the independent sepsis cohort",
    
    x =
      "Spearman rho with SOFA",
    
    y =
      NULL,
    
    color =
      "Expected direction"
  ) +
  
  ggplot2::theme_bw(
    base_size = 14
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 30. FIGURE E - LOCATION-SPECIFIC CORRELATIONS
# ==============================================================================

pE_data <- location_sofa_results %>%
  
  dplyr::mutate(
    
    collection_location =
      reorder(
        collection_location,
        rho
      )
  )


pE <- ggplot2::ggplot(
  
  pE_data,
  
  ggplot2::aes(
    x =
      rho,
    y =
      collection_location
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend =
        rho,
      y =
        collection_location,
      yend =
        collection_location
    ),
    
    linewidth = 1,
    color = "grey55"
  ) +
  
  ggplot2::geom_point(
    size = 4,
    color = "#D55E00"
  ) +
  
  ggplot2::labs(
    
    title =
      "Cross-location replication",
    
    subtitle =
      "Five-gene score versus SOFA within geographic cohorts",
    
    x =
      "Spearman rho",
    
    y =
      NULL
  ) +
  
  ggplot2::theme_bw(
    base_size = 14
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 31. FIGURE F - CONTEXTUAL SEPSIS vs HEALTHY
# ==============================================================================

pF <- ggplot2::ggplot(
  
  disease_context,
  
  ggplot2::aes(
    x =
      disease_state,
    y =
      five_gene_score_all_reference,
    fill =
      disease_state,
    color =
      disease_state
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.62,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    size = 1.4,
    alpha = 0.35
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      disease_colors
  ) +
  
  ggplot2::scale_color_manual(
    values =
      disease_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Contextual sepsis-versus-healthy comparison",
    
    subtitle =
      paste0(
        "Descriptive only; AUC = ",
        sprintf(
          "%.3f",
          secondary_disease$AUC
        )
      ),
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  ) +
  
  ggplot2::theme_bw(
    base_size = 14
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 32. SAVE FIGURES
# ==============================================================================

save_plot_all_formats(
  pA,
  "142b_Figure_A_score_vs_SOFA",
  width = 9,
  height = 7
)


save_plot_all_formats(
  pB,
  "142b_Figure_B_score_by_SOFA_ge2",
  width = 8,
  height = 7
)


save_plot_all_formats(
  pC,
  "142b_Figure_C_score_by_mortality",
  width = 8,
  height = 7
)


save_plot_all_formats(
  pD,
  "142b_Figure_D_component_gene_SOFA_correlations",
  width = 9,
  height = 6
)


save_plot_all_formats(
  pE,
  "142b_Figure_E_location_specific_SOFA_correlations",
  width = 9,
  height = 6
)


save_plot_all_formats(
  pF,
  "142b_Supplementary_Figure_sepsis_vs_healthy",
  width = 8,
  height = 7
)


integrated_figure <- (
  
  pA +
    pB
  
) / (
  
  pD +
    pE
  
) +
  
  patchwork::plot_annotation(
    
    title =
      paste0(
        "External severity validation of the frozen ",
        "five-gene host-response score"
      ),
    
    subtitle =
      "GSE185263 whole-blood RNA-seq"
  )


save_plot_all_formats(
  integrated_figure,
  "142b_Figure_external_severity_validation_integrated",
  width = 16,
  height = 12
)


cat(
  "Publication-quality figures saved.\n\n"
)


# ==============================================================================
# 33. EXPORT CSV TABLES
# ==============================================================================

write.csv(
  scores,
  file.path(
    tables_dir,
    "142b_GSE185263_sample_scores.csv"
  ),
  row.names = FALSE
)


write.csv(
  gene_coverage,
  file.path(
    tables_dir,
    "142b_frozen_gene_coverage.csv"
  ),
  row.names = FALSE
)


write.csv(
  primary_sofa,
  file.path(
    tables_dir,
    "142b_PRIMARY_score_vs_SOFA.csv"
  ),
  row.names = FALSE
)


write.csv(
  secondary_results,
  file.path(
    tables_dir,
    "142b_secondary_score_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  gene_sofa_results,
  file.path(
    tables_dir,
    "142b_component_gene_SOFA_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  adjusted_coefficients,
  file.path(
    tables_dir,
    "142b_age_sex_location_adjusted_model.csv"
  ),
  row.names = FALSE
)


write.csv(
  adjusted_model_summary,
  file.path(
    tables_dir,
    "142b_adjusted_model_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  location_sofa_results,
  file.path(
    tables_dir,
    "142b_location_specific_SOFA_correlations.csv"
  ),
  row.names = FALSE
)


write.csv(
  pooled_location_result,
  file.path(
    tables_dir,
    "142b_cross_location_pooled_correlation.csv"
  ),
  row.names = FALSE
)


write.csv(
  scaling_results,
  file.path(
    tables_dir,
    "142b_scaling_sensitivity.csv"
  ),
  row.names = FALSE
)


write.csv(
  sofa_summary,
  file.path(
    tables_dir,
    "142b_SOFA_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  sofa_group_summary,
  file.path(
    tables_dir,
    "142b_SOFA_group_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  mortality_summary,
  file.path(
    tables_dir,
    "142b_mortality_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  site_summary,
  file.path(
    tables_dir,
    "142b_collection_site_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  location_summary,
  file.path(
    tables_dir,
    "142b_collection_location_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 34. EXCEL WORKBOOK
# ==============================================================================

run_info <- tibble::tibble(
  
  parameter = c(
    "script",
    "run_date",
    "GEO_accession",
    "total_samples",
    "sepsis_samples",
    "healthy_samples",
    "SOFA_available_sepsis",
    "mortality_available_sepsis",
    "primary_panel",
    "primary_score",
    "primary_endpoint",
    "feature_selection_external",
    "coefficient_refitting_external",
    "cutoff_optimization_external",
    "score_direction_flipping_external",
    "endotype_reconstruction"
  ),
  
  value = c(
    
    script_name,
    
    as.character(
      run_date
    ),
    
    gse_id,
    
    as.character(
      nrow(
        scores
      )
    ),
    
    as.character(
      sum(
        scores$disease_state ==
          "sepsis",
        na.rm = TRUE
      )
    ),
    
    as.character(
      sum(
        scores$disease_state ==
          "healthy",
        na.rm = TRUE
      )
    ),
    
    as.character(
      sofa_available_n
    ),
    
    as.character(
      mortality_available_n
    ),
    
    paste(
      five_genes,
      collapse = "; "
    ),
    
    paste0(
      "mean z(CD177,HK3,IRAK3) - ",
      "mean z(CARD11,IKZF2)"
    ),
    
    paste0(
      "Spearman five-gene score versus ",
      "continuous 24-h SOFA in sepsis"
    ),
    
    "NO",
    
    "NO",
    
    "NO",
    
    "NO",
    
    "NO"
  )
)


workbook_file <- file.path(
  tables_dir,
  "142b_GSE185263_external_validation.xlsx"
)


wb <- openxlsx::createWorkbook()


sheet_list <- list(
  
  "00_run_info" =
    run_info,
  
  "01_PRIMARY_SOFA" =
    primary_sofa,
  
  "02_secondary" =
    secondary_results,
  
  "03_gene_SOFA" =
    gene_sofa_results,
  
  "04_adjusted_model" =
    adjusted_coefficients,
  
  "05_adjusted_summary" =
    adjusted_model_summary,
  
  "06_location_SOFA" =
    location_sofa_results,
  
  "07_location_pooled" =
    pooled_location_result,
  
  "08_scaling" =
    scaling_results,
  
  "09_SOFA_summary" =
    sofa_summary,
  
  "10_SOFA_groups" =
    sofa_group_summary,
  
  "11_mortality" =
    mortality_summary,
  
  "12_site_summary" =
    site_summary,
  
  "13_location_summary" =
    location_summary,
  
  "14_gene_coverage" =
    gene_coverage,
  
  "15_sample_scores" =
    scores
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)


for (
  sheet_name in names(
    sheet_list
  )
) {
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  
  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = sheet_list[[sheet_name]],
    headerStyle = header_style
  )
  
  
  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  
  openxlsx::setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_len(
      max(
        1,
        ncol(
          sheet_list[[sheet_name]]
        )
      )
    ),
    widths = "auto"
  )
}


openxlsx::saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)


# ==============================================================================
# 35. AUTOMATED MANUSCRIPT SUMMARY
# ==============================================================================

sofa_secondary_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "SOFA_ge2_vs_SOFA_0_1"
  )


mortality_secondary_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "Died_vs_Survived"
  )


site_secondary_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "ICU_vs_Emergency_Room"
  )


context_secondary_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "Sepsis_vs_healthy_contextual"
  )


summary_en <- c(
  
  "GSE185263 EXTERNAL SEVERITY VALIDATION",
  
  "====================================================================",
  
  "",
  
  paste0(
    "The frozen five-gene host-response score was evaluated in ",
    nrow(
      sepsis_scores
    ),
    " sepsis samples from GSE185263 without feature selection, ",
    "coefficient refitting, cutoff optimization, or score-direction reversal."
  ),
  
  "",
  
  "PRIMARY EXTERNAL ENDPOINT:",
  
  paste0(
    "Among ",
    primary_sofa$n,
    " sepsis samples with available 24-h SOFA, the five-gene score ",
    "showed Spearman rho = ",
    sprintf(
      "%.3f",
      primary_sofa$rho
    ),
    " with SOFA (p = ",
    format.pval(
      primary_sofa$p_value,
      digits = 4
    ),
    ")."
  ),
  
  "",
  
  "SECONDARY SOFA ANALYSIS:",
  
  paste0(
    "SOFA >=2 versus SOFA 0-1: median score ",
    sprintf(
      "%.3f",
      sofa_secondary_row$median_case
    ),
    " versus ",
    sprintf(
      "%.3f",
      sofa_secondary_row$median_control
    ),
    "; p = ",
    format.pval(
      sofa_secondary_row$p_value,
      digits = 4
    ),
    "; AUC = ",
    sprintf(
      "%.3f",
      sofa_secondary_row$AUC
    ),
    "."
  ),
  
  "",
  
  "MORTALITY CONTEXT:",
  
  paste0(
    "Died versus Survived: median score ",
    sprintf(
      "%.3f",
      mortality_secondary_row$median_case
    ),
    " versus ",
    sprintf(
      "%.3f",
      mortality_secondary_row$median_control
    ),
    "; p = ",
    format.pval(
      mortality_secondary_row$p_value,
      digits = 4
    ),
    "; AUC = ",
    sprintf(
      "%.3f",
      mortality_secondary_row$AUC
    ),
    "."
  ),
  
  "",
  
  "COLLECTION-SITE CONTEXT:",
  
  paste0(
    "ICU versus Emergency Room: p = ",
    format.pval(
      site_secondary_row$p_value,
      digits = 4
    ),
    "; AUC = ",
    sprintf(
      "%.3f",
      site_secondary_row$AUC
    ),
    "."
  ),
  
  "",
  
  "GENE-LEVEL DIRECTION:",
  
  paste0(
    gene_concordant_n,
    " of 5 frozen component genes showed the prespecified direction ",
    "of association with SOFA."
  ),
  
  "",
  
  "COVARIATE SENSITIVITY:",
  
  paste0(
    "After adjustment for age, sex, and collection location, ",
    "the SOFA coefficient was ",
    sprintf(
      "%.4f",
      adjusted_sofa_row$estimate
    ),
    " five-gene score units per SOFA point (p = ",
    format.pval(
      adjusted_sofa_row$p_value,
      digits = 4
    ),
    ")."
  ),
  
  "",
  
  "CROSS-LOCATION SYNTHESIS:",
  
  paste0(
    "The descriptive fixed-effect Fisher-z synthesis across eligible ",
    "locations gave pooled rho = ",
    sprintf(
      "%.3f",
      pooled_location_result$pooled_rho
    ),
    " (95% CI ",
    sprintf(
      "%.3f",
      pooled_location_result$CI_low
    ),
    " to ",
    sprintf(
      "%.3f",
      pooled_location_result$CI_high
    ),
    ")."
  ),
  
  "",
  
  "CONTEXTUAL SEPSIS-vs-HEALTHY COMPARISON:",
  
  paste0(
    "AUC = ",
    sprintf(
      "%.3f",
      context_secondary_row$AUC
    ),
    ". This comparison is descriptive only."
  ),
  
  "",
  
  "INTERPRETATION:",
  
  paste0(
    "The primary question is whether the frozen five-gene score tracks ",
    "organ-dysfunction severity within sepsis. Diagnostic, collection-site, ",
    "and mortality analyses are secondary and are not clinical assay validation."
  )
)


summary_ru <- c(
  
  "GSE185263 - ВНЕШНЯЯ ВАЛИДАЦИЯ СВЯЗИ С ТЯЖЕСТЬЮ",
  
  "====================================================================",
  
  "",
  
  paste0(
    "Замороженный пятигенный host-response score был исследован у ",
    nrow(
      sepsis_scores
    ),
    " образцов сепсиса GSE185263 без нового отбора генов, ",
    "изменения коэффициентов, оптимизации порога или изменения ",
    "направления score."
  ),
  
  "",
  
  "ОСНОВНАЯ ВНЕШНЯЯ КОНЕЧНАЯ ТОЧКА:",
  
  paste0(
    "У ",
    primary_sofa$n,
    " образцов сепсиса с доступным 24-часовым SOFA корреляция ",
    "пятигенного score с SOFA составила rho = ",
    sprintf(
      "%.3f",
      primary_sofa$rho
    ),
    ", p = ",
    format.pval(
      primary_sofa$p_value,
      digits = 4
    ),
    "."
  ),
  
  "",
  
  "ВТОРИЧНЫЙ АНАЛИЗ SOFA:",
  
  paste0(
    "SOFA >=2 против SOFA 0-1: медианы score ",
    sprintf(
      "%.3f",
      sofa_secondary_row$median_case
    ),
    " и ",
    sprintf(
      "%.3f",
      sofa_secondary_row$median_control
    ),
    "; p = ",
    format.pval(
      sofa_secondary_row$p_value,
      digits = 4
    ),
    "; AUC = ",
    sprintf(
      "%.3f",
      sofa_secondary_row$AUC
    ),
    "."
  ),
  
  "",
  
  "НАПРАВЛЕНИЕ ОТДЕЛЬНЫХ ГЕНОВ:",
  
  paste0(
    gene_concordant_n,
    " из 5 компонентов панели показали заранее ожидаемое направление ",
    "связи с SOFA."
  ),
  
  "",
  
  "КОВАРИАТНАЯ ЧУВСТВИТЕЛЬНОСТЬ:",
  
  paste0(
    "После коррекции по возрасту, полу и месту набора коэффициент SOFA ",
    "составил ",
    sprintf(
      "%.4f",
      adjusted_sofa_row$estimate
    ),
    " единиц score на один балл SOFA; p = ",
    format.pval(
      adjusted_sofa_row$p_value,
      digits = 4
    ),
    "."
  ),
  
  "",
  
  "ИНТЕРПРЕТАЦИЯ:",
  
  paste0(
    "Главный вопрос анализа - отражает ли замороженный пятигенный score ",
    "тяжесть органной дисфункции внутри сепсиса. Анализы смертности, ",
    "места госпитализации и sepsis-vs-healthy являются вторичными и ",
    "не рассматриваются как валидация клинического теста."
  )
)


summary_en_file <- file.path(
  text_dir,
  "142b_external_validation_summary_EN.txt"
)


summary_ru_file <- file.path(
  text_dir,
  "142b_external_validation_summary_RU.txt"
)


writeLines(
  summary_en,
  con = summary_en_file,
  useBytes = TRUE
)


writeLines(
  summary_ru,
  con = summary_ru_file,
  useBytes = TRUE
)


# ==============================================================================
# 36. INPUT MANIFEST
# ==============================================================================

manifest <- tibble::tibble(
  
  item = c(
    "metadata",
    "raw_counts",
    "predeclared_plan"
  ),
  
  path = c(
    metadata_file,
    raw_counts_file,
    plan_file
  )
)


manifest_info <- file.info(
  manifest$path
)


manifest$file_size_bytes <-
  manifest_info$size


manifest$modified_time <-
  as.character(
    manifest_info$mtime
  )


manifest$md5 <-
  unname(
    tools::md5sum(
      manifest$path
    )
  )


write.csv(
  manifest,
  file.path(
    logs_dir,
    "142b_input_manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 37. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "142b_sessionInfo.txt"
  )
)


# ==============================================================================
# 38. FINAL REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 142b completed successfully.\n")
cat("====================================================================\n\n")


cat(
  "GSE185263 COHORT:\n"
)


cat(
  "Sepsis = ",
  sum(
    scores$disease_state ==
      "sepsis",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Healthy = ",
  sum(
    scores$disease_state ==
      "healthy",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Sepsis with SOFA = ",
  primary_sofa$n,
  "\n\n",
  sep = ""
)


cat(
  "PRIMARY EXTERNAL ENDPOINT:\n"
)


cat(
  "Five-gene score vs continuous 24-h SOFA\n"
)


cat(
  "Spearman rho = ",
  sprintf(
    "%.3f",
    primary_sofa$rho
  ),
  "\n",
  sep = ""
)


cat(
  "p = ",
  format.pval(
    primary_sofa$p_value,
    digits = 4
  ),
  "\n\n",
  sep = ""
)


cat(
  "SECONDARY SCORE ASSOCIATIONS:\n"
)


print(
  secondary_results,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "INDIVIDUAL GENE-SOFA ASSOCIATIONS:\n"
)


print(
  gene_sofa_results,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "Expected gene-level direction = ",
  gene_concordant_n,
  "/5 concordant\n\n",
  sep = ""
)


cat(
  "AGE/SEX/LOCATION-ADJUSTED SOFA ASSOCIATION:\n"
)


print(
  adjusted_sofa_row,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "LOCATION-SPECIFIC SOFA ASSOCIATIONS:\n"
)


print(
  location_sofa_results,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "CROSS-LOCATION POOLED RESULT:\n"
)


print(
  pooled_location_result,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "SCALING SENSITIVITY:\n"
)


print(
  scaling_results,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "Main workbook:\n"
)


cat(
  normalizePath(
    workbook_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "Main integrated figure:\n"
)


cat(
  normalizePath(
    file.path(
      figures_dir,
      "142b_Figure_external_severity_validation_integrated.png"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "IMPORTANT:\n"
)


cat(
  "- Primary endpoint frozen before five-gene expression analysis.\n"
)


cat(
  "- No feature selection.\n"
)


cat(
  "- No gene substitution.\n"
)


cat(
  "- No coefficient refitting.\n"
)


cat(
  "- No cutoff optimization.\n"
)


cat(
  "- No score-direction flipping.\n"
)


cat(
  "- No endotype reconstruction.\n"
)


cat(
  "- Sepsis-vs-healthy comparison is contextual only.\n"
)


cat(
  "- Mortality analysis is secondary and not a prognostic validation claim.\n\n"
)


cat(
  "Done.\n"
)