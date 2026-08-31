# ==============================================================================
# Script 141
# External validation of the frozen five-gene blood host-response signature
# Dataset: GSE154918
# Project: Sepsis_DESeq2
#
# FROZEN PANEL
#   UP:   CD177, HK3, IRAK3
#   DOWN: CARD11, IKZF2
#
# SCORE
#   mean[z(CD177), z(HK3), z(IRAK3)] -
#   mean[z(CARD11), z(IKZF2)]
#
# PRIMARY EXTERNAL COMPARISON
#   Baseline Seps_P + Shock_P versus Inf1_P
#
# SECONDARY ANALYSES
#   1. Seps_P + Shock_P versus healthy controls
#   2. Seps_P versus Inf1_P
#   3. Shock_P versus Inf1_P
#   4. Shock_P versus Seps_P
#   5. Ordered baseline states:
#      Hlty -> Inf1_P -> Seps_P -> Shock_P
#
# IMPORTANT
#   - frozen gene composition
#   - no gene substitution
#   - no new feature selection
#   - no coefficient fitting
#   - no cutoff optimization
#   - no use of follow-up samples in primary validation
#
# Interpretation:
# External transcriptomic replication of a frozen molecular signature,
# not validation of a pre-calibrated clinical diagnostic assay.
# ==============================================================================


# ==============================================================================
# 0. GLOBAL SETTINGS
# ==============================================================================

options(
  stringsAsFactors = FALSE
)

set.seed(
  20260817
)


project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))

gse_id <- "GSE154918"

script_name <- "141_external_validation_GSE154918.R"

run_date <- Sys.time()


if (!dir.exists(project_dir)) {
  
  stop(
    paste0(
      "Project directory not found: ",
      project_dir
    )
  )
}


setwd(
  project_dir
)


cat("\n")
cat("====================================================================\n")
cat("Running Script 141\n")
cat("External validation of frozen five-gene signature: GSE154918\n")
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
  as.character(
    run_date
  ),
  "\n\n"
)


# ==============================================================================
# 1. OUTPUT DIRECTORIES
# ==============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "141_external_validation_GSE154918"
)


raw_dir <- file.path(
  output_dir,
  "raw_download"
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
    raw_dir,
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
# 2. PACKAGES
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
  
  cat(
    "Installing missing CRAN packages:\n",
    paste(
      missing_cran,
      collapse = ", "
    ),
    "\n\n"
  )
  
  
  install.packages(
    missing_cran
  )
}


if (
  !requireNamespace(
    "BiocManager",
    quietly = TRUE
  )
) {
  
  install.packages(
    "BiocManager"
  )
}


bioc_packages <- c(
  "GEOquery",
  "Biobase",
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
  
  cat(
    "Installing missing Bioconductor packages:\n",
    paste(
      missing_bioc,
      collapse = ", "
    ),
    "\n\n"
  )
  
  
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
  
  library(GEOquery)
  
  library(Biobase)
  
  library(AnnotationDbi)
  
  library(org.Hs.eg.db)
})


cat(
  "Required packages loaded successfully.\n\n"
)


# ==============================================================================
# 3. FROZEN PRE-ANALYSIS SPECIFICATION
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


expected_direction <- c(
  
  "CD177" =
    "UP",
  
  "HK3" =
    "UP",
  
  "IRAK3" =
    "UP",
  
  "CARD11" =
    "DOWN",
  
  "IKZF2" =
    "DOWN"
)


baseline_status_levels <- c(
  "Hlty",
  "Inf1_P",
  "Seps_P",
  "Shock_P"
)


all_status_levels <- c(
  baseline_status_levels,
  "Seps_FU",
  "Shock_FU"
)


# Expected published group sizes.
# These are used ONLY as an audit.
# They do not determine inclusion.

expected_status_counts <- tibble::tibble(
  
  status = c(
    "Hlty",
    "Inf1_P",
    "Seps_P",
    "Shock_P",
    "Seps_FU",
    "Shock_FU"
  ),
  
  expected_n = c(
    40,
    12,
    20,
    19,
    4,
    10
  )
)


# ==============================================================================
# 4. WRITE PREDECLARED ANALYSIS PLAN
# ==============================================================================

predeclared_plan <- c(
  
  "GSE154918 EXTERNAL VALIDATION - PREDECLARED ANALYSIS PLAN",
  
  "====================================================================",
  
  "",
  
  paste0(
    "Created before expression-level validation: ",
    run_date
  ),
  
  "",
  
  "FROZEN FIVE-GENE PANEL:",
  
  "UP: CD177, HK3, IRAK3",
  
  "DOWN: CARD11, IKZF2",
  
  "",
  
  "FROZEN SCORE:",
  
  paste0(
    "mean[z(CD177), z(HK3), z(IRAK3)] - ",
    "mean[z(CARD11), z(IKZF2)]"
  ),
  
  "",
  
  "PRIMARY EXTERNAL COMPARISON:",
  
  paste0(
    "Baseline sepsis/septic shock (Seps_P + Shock_P) ",
    "versus uncomplicated infection (Inf1_P)."
  ),
  
  "",
  
  "SECONDARY COMPARISONS:",
  
  "1. Seps_P + Shock_P versus Hlty",
  
  "2. Seps_P versus Inf1_P",
  
  "3. Shock_P versus Inf1_P",
  
  "4. Shock_P versus Seps_P",
  
  paste0(
    "5. Ordered baseline states: ",
    "Hlty -> Inf1_P -> Seps_P -> Shock_P"
  ),
  
  "",
  
  "PRIMARY SCALING:",
  
  paste0(
    "Each frozen gene is standardized across all baseline samples ",
    "in the external cohort."
  ),
  
  paste0(
    "No phenotype labels are used to estimate gene weights."
  ),
  
  "",
  
  "SCALING SENSITIVITY:",
  
  paste0(
    "A second implementation standardizes expression using ",
    "healthy-control mean and standard deviation."
  ),
  
  "",
  
  "FROZEN RULES:",
  
  "- no gene substitution",
  
  "- no gene addition",
  
  "- no gene deletion",
  
  "- no feature selection",
  
  "- no coefficient refitting",
  
  "- no cutoff optimization",
  
  "- no post hoc score-direction flipping",
  
  "- follow-up samples excluded from primary validation",
  
  "",
  
  "INTERPRETATION:",
  
  paste0(
    "This is external transcriptomic replication of a frozen gene ",
    "composition and score direction."
  ),
  
  paste0(
    "It is not validation of a pre-calibrated clinical diagnostic ",
    "assay or clinical decision threshold."
  )
)


plan_file <- file.path(
  text_dir,
  "141_PREDECLARED_external_validation_plan.txt"
)


writeLines(
  predeclared_plan,
  con = plan_file,
  useBytes = TRUE
)


cat(
  "Frozen pre-analysis plan written before validation.\n\n"
)


# ==============================================================================
# 5. HELPER FUNCTIONS
# ==============================================================================

normalize_key <- function(x) {
  
  x <- as.character(
    x
  )
  
  x <- tolower(
    x
  )
  
  x <- gsub(
    "[^a-z0-9]",
    "",
    x
  )
  
  return(
    x
  )
}


safe_numeric <- function(x) {
  
  suppressWarnings(
    
    as.numeric(
      as.character(
        x
      )
    )
  )
}


safe_z_rows <- function(mat) {
  
  mat <- as.matrix(
    mat
  )
  
  
  storage.mode(
    mat
  ) <- "numeric"
  
  
  row_means <- rowMeans(
    mat,
    na.rm = TRUE
  )
  
  
  row_sds <- apply(
    mat,
    1,
    stats::sd,
    na.rm = TRUE
  )
  
  
  bad <- !is.finite(
    row_sds
  ) |
    row_sds == 0
  
  
  if (any(bad)) {
    
    stop(
      paste0(
        "Cannot z-standardize gene(s) with zero/non-finite SD: ",
        paste(
          rownames(
            mat
          )[bad],
          collapse = ", "
        )
      )
    )
  }
  
  
  out <- sweep(
    mat,
    1,
    row_means,
    FUN = "-"
  )
  
  
  out <- sweep(
    out,
    1,
    row_sds,
    FUN = "/"
  )
  
  
  return(
    out
  )
}


healthy_reference_z <- function(
    mat,
    healthy_samples
) {
  
  mat <- as.matrix(
    mat
  )
  
  
  storage.mode(
    mat
  ) <- "numeric"
  
  
  if (
    !all(
      healthy_samples %in%
      colnames(
        mat
      )
    )
  ) {
    
    stop(
      paste0(
        "Healthy-reference samples are not all present ",
        "in the expression matrix."
      )
    )
  }
  
  
  ref_mat <- mat[
    ,
    healthy_samples,
    drop = FALSE
  ]
  
  
  ref_mean <- rowMeans(
    ref_mat,
    na.rm = TRUE
  )
  
  
  ref_sd <- apply(
    ref_mat,
    1,
    stats::sd,
    na.rm = TRUE
  )
  
  
  bad <- !is.finite(
    ref_sd
  ) |
    ref_sd == 0
  
  
  if (any(bad)) {
    
    stop(
      paste0(
        "Healthy-reference SD is zero/non-finite for: ",
        paste(
          rownames(
            mat
          )[bad],
          collapse = ", "
        )
      )
    )
  }
  
  
  out <- sweep(
    mat,
    1,
    ref_mean,
    FUN = "-"
  )
  
  
  out <- sweep(
    out,
    1,
    ref_sd,
    FUN = "/"
  )
  
  
  return(
    out
  )
}


calculate_five_gene_score <- function(z_mat) {
  
  if (
    !all(
      five_genes %in%
      rownames(
        z_mat
      )
    )
  ) {
    
    stop(
      "Not all five frozen genes are present in standardized matrix."
    )
  }
  
  
  up_component <- colMeans(
    
    z_mat[
      up_genes,
      ,
      drop = FALSE
    ],
    
    na.rm = TRUE
  )
  
  
  down_component <- colMeans(
    
    z_mat[
      down_genes,
      ,
      drop = FALSE
    ],
    
    na.rm = TRUE
  )
  
  
  score <- up_component -
    down_component
  
  
  return(
    score
  )
}


extract_status_from_row <- function(
    row_values
) {
  
  txt <- paste(
    as.character(
      row_values
    ),
    collapse = " | "
  )
  
  
  hits <- all_status_levels[
    vapply(
      all_status_levels,
      function(x) {
        
        stringr::str_detect(
          txt,
          stringr::fixed(
            x
          )
        )
      },
      logical(1)
    )
  ]
  
  
  if (length(hits) == 0) {
    
    return(
      NA_character_
    )
  }
  
  
  if (length(hits) == 1) {
    
    return(
      hits[1]
    )
  }
  
  
  # In case multiple substrings are found,
  # retain the longest matching label.
  
  hits <- hits[
    order(
      nchar(
        hits
      ),
      decreasing = TRUE
    )
  ]
  
  
  return(
    hits[1]
  )
}


extract_sex_from_row <- function(
    row_values
) {
  
  txt <- paste(
    as.character(
      row_values
    ),
    collapse = " | "
  )
  
  
  txt_lower <- tolower(
    txt
  )
  
  
  if (
    grepl(
      "sex:[[:space:]]*female",
      txt_lower,
      perl = TRUE
    )
  ) {
    
    return(
      "Female"
    )
  }
  
  
  if (
    grepl(
      "sex:[[:space:]]*male",
      txt_lower,
      perl = TRUE
    )
  ) {
    
    return(
      "Male"
    )
  }
  
  
  if (
    grepl(
      "sex:[[:space:]]*f([^a-z]|$)",
      txt_lower,
      perl = TRUE
    )
  ) {
    
    return(
      "Female"
    )
  }
  
  
  if (
    grepl(
      "sex:[[:space:]]*m([^a-z]|$)",
      txt_lower,
      perl = TRUE
    )
  ) {
    
    return(
      "Male"
    )
  }
  
  
  return(
    NA_character_
  )
}


safe_wilcox <- function(
    x_case,
    x_control
) {
  
  x_case <- x_case[
    is.finite(
      x_case
    )
  ]
  
  
  x_control <- x_control[
    is.finite(
      x_control
    )
  ]
  
  
  if (
    length(
      x_case
    ) < 2 ||
    length(
      x_control
    ) < 2
  ) {
    
    return(
      list(
        p = NA_real_,
        statistic = NA_real_
      )
    )
  }
  
  
  ht <- stats::wilcox.test(
    
    x =
      x_case,
    
    y =
      x_control,
    
    exact =
      FALSE,
    
    paired =
      FALSE
  )
  
  
  return(
    list(
      
      p =
        as.numeric(
          ht$p.value
        ),
      
      statistic =
        unname(
          as.numeric(
            ht$statistic
          )
        )
    )
  )
}


fixed_direction_auc <- function(
    x_case,
    x_control
) {
  
  x_case <- x_case[
    is.finite(
      x_case
    )
  ]
  
  
  x_control <- x_control[
    is.finite(
      x_control
    )
  ]
  
  
  if (
    length(
      x_case
    ) < 3 ||
    length(
      x_control
    ) < 3
  ) {
    
    return(
      tibble::tibble(
        
        auc =
          NA_real_,
        
        ci_low =
          NA_real_,
        
        ci_high =
          NA_real_
      )
    )
  }
  
  
  response <- factor(
    
    c(
      
      rep(
        "control",
        length(
          x_control
        )
      ),
      
      rep(
        "case",
        length(
          x_case
        )
      )
    ),
    
    levels = c(
      "control",
      "case"
    )
  )
  
  
  predictor <- c(
    x_control,
    x_case
  )
  
  
  roc_object <- pROC::roc(
    
    response =
      response,
    
    predictor =
      predictor,
    
    levels = c(
      "control",
      "case"
    ),
    
    direction =
      "<",
    
    quiet =
      TRUE
  )
  
  
  roc_ci <- suppressWarnings(
    
    pROC::ci.auc(
      roc_object,
      method = "delong"
    )
  )
  
  
  return(
    tibble::tibble(
      
      auc =
        as.numeric(
          pROC::auc(
            roc_object
          )
        ),
      
      ci_low =
        as.numeric(
          roc_ci[1]
        ),
      
      ci_high =
        as.numeric(
          roc_ci[3]
        )
    )
  )
}


compare_groups <- function(
    df,
    score_column,
    case_status,
    control_status,
    comparison_name
) {
  
  case_values <- df %>%
    
    dplyr::filter(
      status %in%
        case_status
    ) %>%
    
    dplyr::pull(
      dplyr::all_of(
        score_column
      )
    )
  
  
  control_values <- df %>%
    
    dplyr::filter(
      status %in%
        control_status
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
  
  
  auc_table <- fixed_direction_auc(
    case_values,
    control_values
  )
  
  
  result <- tibble::tibble(
    
    comparison =
      comparison_name,
    
    score =
      score_column,
    
    case_status =
      paste(
        case_status,
        collapse = "+"
      ),
    
    control_status =
      paste(
        control_status,
        collapse = "+"
      ),
    
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
    
    case_median =
      stats::median(
        case_values,
        na.rm = TRUE
      ),
    
    control_median =
      stats::median(
        control_values,
        na.rm = TRUE
      ),
    
    median_difference_case_minus_control =
      stats::median(
        case_values,
        na.rm = TRUE
      ) -
      stats::median(
        control_values,
        na.rm = TRUE
      ),
    
    wilcoxon_W =
      wt$statistic,
    
    p_value =
      wt$p,
    
    auc_fixed_direction =
      auc_table$auc,
    
    auc_ci_low =
      auc_table$ci_low,
    
    auc_ci_high =
      auc_table$ci_high
  )
  
  
  return(
    result
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
    
    filename =
      png_file,
    
    plot =
      plot_object,
    
    width =
      width,
    
    height =
      height,
    
    dpi =
      600,
    
    bg =
      "white"
  )
  
  
  ggplot2::ggsave(
    
    filename =
      tiff_file,
    
    plot =
      plot_object,
    
    width =
      width,
    
    height =
      height,
    
    dpi =
      600,
    
    compression =
      "lzw",
    
    bg =
      "white"
  )
  
  
  if (
    isTRUE(
      capabilities(
        "cairo"
      )
    )
  ) {
    
    ggplot2::ggsave(
      
      filename =
        pdf_file,
      
      plot =
        plot_object,
      
      width =
        width,
      
      height =
        height,
      
      device =
        grDevices::cairo_pdf,
      
      bg =
        "white"
    )
    
  } else {
    
    ggplot2::ggsave(
      
      filename =
        pdf_file,
      
      plot =
        plot_object,
      
      width =
        width,
      
      height =
        height,
      
      bg =
        "white"
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
# 6. DOWNLOAD GEO METADATA
# ==============================================================================

cat(
  "Downloading GEO metadata for ",
  gse_id,
  " ...\n",
  sep = ""
)


gse_raw <- GEOquery::getGEO(
  
  gse_id,
  
  GSEMatrix =
    TRUE,
  
  getGPL =
    FALSE
)


if (
  is.list(
    gse_raw
  ) &&
  !inherits(
    gse_raw,
    "ExpressionSet"
  )
) {
  
  geo_object <- gse_raw[[1]]
  
} else {
  
  geo_object <- gse_raw
}


if (
  !inherits(
    geo_object,
    "ExpressionSet"
  )
) {
  
  stop(
    paste0(
      "Expected ExpressionSet from GEOquery. Observed class: ",
      paste(
        class(
          geo_object
        ),
        collapse = ", "
      )
    )
  )
}


pheno_raw <- as.data.frame(
  
  Biobase::pData(
    geo_object
  ),
  
  stringsAsFactors =
    FALSE
)


pheno_raw$.geo_rowname <- rownames(
  pheno_raw
)


cat(
  "GEO phenotype table: ",
  nrow(
    pheno_raw
  ),
  " samples x ",
  ncol(
    pheno_raw
  ),
  " columns\n\n",
  sep = ""
)


# ==============================================================================
# 7. STANDARDIZE GEO METADATA
# ==============================================================================

if (
  "geo_accession" %in%
  names(
    pheno_raw
  )
) {
  
  geo_accession_vector <- as.character(
    pheno_raw$geo_accession
  )
  
} else {
  
  geo_accession_vector <- as.character(
    pheno_raw$.geo_rowname
  )
}


if (
  "title" %in%
  names(
    pheno_raw
  )
) {
  
  sample_title_vector <- as.character(
    pheno_raw$title
  )
  
} else {
  
  sample_title_vector <-
    geo_accession_vector
}


status_vector <- apply(
  pheno_raw,
  1,
  extract_status_from_row
)


sex_vector <- apply(
  pheno_raw,
  1,
  extract_sex_from_row
)


metadata <- pheno_raw %>%
  
  dplyr::mutate(
    
    geo_accession =
      geo_accession_vector,
    
    sample_title =
      sample_title_vector,
    
    status =
      status_vector,
    
    sex_standardized =
      sex_vector
  ) %>%
  
  dplyr::select(
    
    geo_accession,
    
    sample_title,
    
    status,
    
    sex_standardized,
    
    dplyr::everything()
  )


cat(
  "Detected GEO status distribution:\n"
)


print(
  table(
    metadata$status,
    useNA = "ifany"
  )
)


cat("\n")


unresolved_status <- metadata %>%
  
  dplyr::filter(
    is.na(
      status
    )
  ) %>%
  
  dplyr::select(
    geo_accession,
    sample_title
  )


if (
  nrow(
    unresolved_status
  ) > 0
) {
  
  cat(
    "WARNING: unresolved status for these samples:\n"
  )
  
  
  print(
    unresolved_status
  )
  
  
  cat("\n")
}


observed_status_counts <- metadata %>%
  
  dplyr::filter(
    !is.na(
      status
    )
  ) %>%
  
  dplyr::count(
    status,
    name = "observed_n"
  )


status_count_check <- expected_status_counts %>%
  
  dplyr::left_join(
    
    observed_status_counts,
    
    by =
      "status"
  ) %>%
  
  dplyr::mutate(
    
    observed_n =
      tidyr::replace_na(
        observed_n,
        0L
      ),
    
    count_matches_expected =
      expected_n ==
      observed_n
  )


cat(
  "Status-count audit:\n"
)


print(
  status_count_check,
  n = Inf
)


cat("\n")


if (
  !all(
    status_count_check$count_matches_expected
  )
) {
  
  warning(
    paste0(
      "Observed status counts differ from expected published counts. ",
      "Inspect metadata before final interpretation."
    )
  )
}


write.csv(
  
  metadata,
  
  file.path(
    tables_dir,
    "141_GSE154918_metadata_all_samples.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  status_count_check,
  
  file.path(
    tables_dir,
    "141_GSE154918_status_count_check.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 8. DOWNLOAD SUPPLEMENTARY PROCESSED EXPRESSION FILE
# ==============================================================================

supp_filename <-
  "GSE154918_Schughart_Sepsis_200320.txt.gz"


supp_file <- file.path(
  raw_dir,
  supp_filename
)


if (
  !file.exists(
    supp_file
  )
) {
  
  cat(
    "Downloading GEO supplementary files...\n"
  )
  
  
  GEOquery::getGEOSuppFiles(
    
    gse_id,
    
    makeDirectory =
      TRUE,
    
    baseDir =
      raw_dir
  )
  
  
  candidate_files <- list.files(
    
    raw_dir,
    
    pattern =
      "GSE154918_Schughart_Sepsis_200320\\.txt\\.gz$",
    
    recursive =
      TRUE,
    
    full.names =
      TRUE
  )
  
  
  if (
    length(
      candidate_files
    ) == 0
  ) {
    
    stop(
      paste0(
        "Could not locate supplementary expression file: ",
        supp_filename
      )
    )
  }
  
  
  candidate_file <- candidate_files[1]
  
  
  candidate_norm <- normalizePath(
    candidate_file,
    winslash = "/",
    mustWork = FALSE
  )
  
  
  target_norm <- normalizePath(
    supp_file,
    winslash = "/",
    mustWork = FALSE
  )
  
  
  if (
    !identical(
      candidate_norm,
      target_norm
    )
  ) {
    
    copy_ok <- file.copy(
      
      from =
        candidate_file,
      
      to =
        supp_file,
      
      overwrite =
        TRUE
    )
    
    
    if (!isTRUE(copy_ok)) {
      
      stop(
        "Failed to copy supplementary expression file into raw_download."
      )
    }
  }
}


if (
  !file.exists(
    supp_file
  )
) {
  
  stop(
    paste0(
      "Processed GEO expression file not found: ",
      supp_file
    )
  )
}


cat(
  "Processed expression file:\n"
)


cat(
  normalizePath(
    supp_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)


cat(
  "File size: ",
  round(
    file.info(
      supp_file
    )$size /
      1024^2,
    2
  ),
  " MB\n\n",
  sep = ""
)


# ==============================================================================
# 9. READ PROCESSED EXPRESSION DATA
# ==============================================================================

cat(
  "Reading processed expression matrix...\n"
)


expr_raw <- tryCatch(
  
  {
    
    data.table::fread(
      
      supp_file,
      
      data.table =
        FALSE,
      
      check.names =
        FALSE
    )
  },
  
  error = function(e) {
    
    message(
      paste0(
        "fread() failed: ",
        e$message,
        "\nUsing read.delim(gzfile()) fallback."
      )
    )
    
    
    utils::read.delim(
      
      gzfile(
        supp_file
      ),
      
      header =
        TRUE,
      
      check.names =
        FALSE,
      
      stringsAsFactors =
        FALSE
    )
  }
)


cat(
  "Raw processed table: ",
  nrow(
    expr_raw
  ),
  " rows x ",
  ncol(
    expr_raw
  ),
  " columns\n\n",
  sep = ""
)


cat(
  "First column names:\n"
)


print(
  head(
    names(
      expr_raw
    ),
    20
  )
)


cat("\n")


# ==============================================================================
# 10. MATCH EXPRESSION COLUMNS TO GEO SAMPLES
# ==============================================================================

meta_lookup <- metadata %>%
  
  dplyr::select(
    geo_accession,
    sample_title,
    status
  ) %>%
  
  dplyr::mutate(
    
    title_key =
      normalize_key(
        sample_title
      ),
    
    geo_key =
      normalize_key(
        geo_accession
      )
  )


title_to_geo <- stats::setNames(
  meta_lookup$geo_accession,
  meta_lookup$title_key
)


geo_to_geo <- stats::setNames(
  meta_lookup$geo_accession,
  meta_lookup$geo_key
)


expression_column_mapping <- tibble::tibble(
  
  expression_column =
    names(
      expr_raw
    ),
  
  expression_key =
    normalize_key(
      names(
        expr_raw
      )
    )
) %>%
  
  dplyr::mutate(
    
    geo_from_title =
      unname(
        title_to_geo[
          expression_key
        ]
      ),
    
    geo_from_accession =
      unname(
        geo_to_geo[
          expression_key
        ]
      ),
    
    geo_accession =
      dplyr::coalesce(
        geo_from_title,
        geo_from_accession
      )
  )


sample_column_mapping <- expression_column_mapping %>%
  
  dplyr::filter(
    !is.na(
      geo_accession
    )
  ) %>%
  
  dplyr::left_join(
    
    meta_lookup %>%
      
      dplyr::select(
        geo_accession,
        sample_title,
        status
      ),
    
    by =
      "geo_accession"
  )


cat(
  "Expression columns matched to GEO samples: ",
  nrow(
    sample_column_mapping
  ),
  "\n\n",
  sep = ""
)


if (
  nrow(
    sample_column_mapping
  ) < 90
) {
  
  cat(
    "Matched expression columns:\n"
  )
  
  
  print(
    sample_column_mapping,
    n = Inf
  )
  
  
  cat("\n")
  
  
  stop(
    paste0(
      "Too few sample columns were matched to GEO metadata. ",
      "Inspect supplementary matrix structure."
    )
  )
}


if (
  anyDuplicated(
    sample_column_mapping$geo_accession
  ) > 0
) {
  
  stop(
    "Duplicate GEO accessions detected after expression-column matching."
  )
}


sample_cols <-
  sample_column_mapping$expression_column


non_sample_cols <- setdiff(
  names(
    expr_raw
  ),
  sample_cols
)


cat(
  "Non-sample columns in processed file:\n"
)


print(
  non_sample_cols
)


cat("\n")


write.csv(
  
  sample_column_mapping,
  
  file.path(
    tables_dir,
    "141_expression_sample_column_mapping.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 11. DETECT AND MAP GENE IDENTIFIER
# ==============================================================================

symbol_candidates <- non_sample_cols[
  grepl(
    "symbol|external.*gene.*name|gene.*name",
    non_sample_cols,
    ignore.case = TRUE
  )
]


ensembl_candidates <- non_sample_cols[
  grepl(
    "ensembl|gene.*id",
    non_sample_cols,
    ignore.case = TRUE
  )
]


generic_gene_candidates <- non_sample_cols[
  grepl(
    "gene|id",
    non_sample_cols,
    ignore.case = TRUE
  )
]


gene_id_column_used <-
  NA_character_


gene_id_type <-
  NA_character_


gene_symbol <-
  NULL


if (
  length(
    symbol_candidates
  ) > 0
) {
  
  gene_id_column_used <-
    symbol_candidates[1]
  
  
  gene_symbol <- as.character(
    expr_raw[[gene_id_column_used]]
  )
  
  
  gene_id_type <-
    "gene_symbol_column"
  
} else {
  
  if (
    length(
      ensembl_candidates
    ) > 0
  ) {
    
    gene_id_column_used <-
      ensembl_candidates[1]
    
  } else if (
    length(
      generic_gene_candidates
    ) > 0
  ) {
    
    gene_id_column_used <-
      generic_gene_candidates[1]
    
  } else if (
    length(
      non_sample_cols
    ) > 0
  ) {
    
    gene_id_column_used <-
      non_sample_cols[1]
    
  } else {
    
    stop(
      "No non-sample column is available as a gene identifier."
    )
  }
  
  
  raw_gene_id <- as.character(
    expr_raw[[gene_id_column_used]]
  )
  
  
  ensembl_fraction <- mean(
    
    grepl(
      "^ENSG[0-9]+",
      raw_gene_id
    ),
    
    na.rm =
      TRUE
  )
  
  
  if (
    is.finite(
      ensembl_fraction
    ) &&
    ensembl_fraction > 0.5
  ) {
    
    gene_id_type <-
      "ENSEMBL"
    
    
    ensembl_clean <- sub(
      "\\.[0-9]+$",
      "",
      raw_gene_id
    )
    
    
    symbol_map <- AnnotationDbi::mapIds(
      
      org.Hs.eg.db,
      
      keys =
        unique(
          ensembl_clean
        ),
      
      column =
        "SYMBOL",
      
      keytype =
        "ENSEMBL",
      
      multiVals =
        "first"
    )
    
    
    gene_symbol <- unname(
      symbol_map[
        ensembl_clean
      ]
    )
    
  } else {
    
    gene_id_type <-
      "assumed_gene_symbol"
    
    
    gene_symbol <-
      raw_gene_id
  }
}


gene_symbol <- toupper(
  trimws(
    as.character(
      gene_symbol
    )
  )
)


cat(
  "Gene identifier column used:\n"
)


cat(
  gene_id_column_used,
  "\n"
)


cat(
  "Gene identifier interpretation:\n"
)


cat(
  gene_id_type,
  "\n\n"
)


# ==============================================================================
# 12. BUILD GENE x SAMPLE EXPRESSION MATRIX
# ==============================================================================

expr_sample_df <- expr_raw[
  ,
  sample_cols,
  drop = FALSE
]


expr_sample_df[] <- lapply(
  expr_sample_df,
  safe_numeric
)


expr_dt <- data.table::as.data.table(
  expr_sample_df
)


expr_dt[["gene_symbol"]] <-
  gene_symbol


expr_dt <- expr_dt[
  !is.na(
    gene_symbol
  ) &
    gene_symbol != ""
]


expr_collapsed <- expr_dt[
  ,
  lapply(
    .SD,
    function(x) {
      
      value <- mean(
        as.numeric(
          x
        ),
        na.rm = TRUE
      )
      
      
      if (
        is.nan(
          value
        )
      ) {
        
        return(
          NA_real_
        )
        
      } else {
        
        return(
          value
        )
      }
    }
  ),
  by =
    gene_symbol,
  .SDcols =
    sample_cols
]


gene_mapping_summary <- tibble::tibble(
  
  metric = c(
    "raw_expression_rows",
    "mapped_nonempty_gene_symbols",
    "unique_collapsed_gene_symbols"
  ),
  
  value = c(
    nrow(
      expr_raw
    ),
    nrow(
      expr_dt
    ),
    nrow(
      expr_collapsed
    )
  )
)


cat(
  "Gene mapping summary:\n"
)


print(
  gene_mapping_summary,
  n = Inf
)


cat("\n")


rename_vector <- stats::setNames(
  sample_column_mapping$geo_accession,
  sample_column_mapping$expression_column
)


for (
  old_name in names(
    rename_vector
  )
) {
  
  data.table::setnames(
    
    expr_collapsed,
    
    old =
      old_name,
    
    new =
      rename_vector[[old_name]]
  )
}


expression_columns_after_rename <- setdiff(
  names(
    expr_collapsed
  ),
  "gene_symbol"
)


expr_matrix <- as.matrix(
  
  expr_collapsed[
    ,
    expression_columns_after_rename,
    with = FALSE
  ]
)


rownames(
  expr_matrix
) <- expr_collapsed$gene_symbol


storage.mode(
  expr_matrix
) <- "numeric"


cat(
  "Collapsed expression matrix: ",
  nrow(
    expr_matrix
  ),
  " genes x ",
  ncol(
    expr_matrix
  ),
  " samples\n\n",
  sep = ""
)


# ==============================================================================
# 13. FROZEN FIVE-GENE COVERAGE
# ==============================================================================

gene_coverage <- tibble::tibble(
  
  gene =
    five_genes,
  
  expected_direction =
    unname(
      expected_direction[
        five_genes
      ]
    ),
  
  present =
    five_genes %in%
    rownames(
      expr_matrix
    )
)


cat(
  "Frozen five-gene coverage:\n"
)


print(
  gene_coverage,
  n = Inf
)


cat("\n")


if (
  !all(
    gene_coverage$present
  )
) {
  
  missing_frozen_genes <- gene_coverage$gene[
    !gene_coverage$present
  ]
  
  
  stop(
    paste0(
      "Frozen five-gene panel is incomplete in GSE154918. Missing: ",
      paste(
        missing_frozen_genes,
        collapse = ", "
      )
    )
  )
}


cat(
  "All five frozen genes detected.\n\n"
)


# ==============================================================================
# 14. DEFINE BASELINE EXTERNAL VALIDATION COHORT
# ==============================================================================

baseline_meta <- metadata %>%
  
  dplyr::filter(
    status %in%
      baseline_status_levels
  ) %>%
  
  dplyr::mutate(
    
    status =
      factor(
        status,
        levels =
          baseline_status_levels,
        ordered =
          TRUE
      )
  )


missing_baseline_expression <- setdiff(
  baseline_meta$geo_accession,
  colnames(
    expr_matrix
  )
)


if (
  length(
    missing_baseline_expression
  ) > 0
) {
  
  stop(
    paste0(
      "Baseline samples missing from expression matrix: ",
      paste(
        missing_baseline_expression,
        collapse = ", "
      )
    )
  )
}


baseline_expr <- expr_matrix[
  five_genes,
  baseline_meta$geo_accession,
  drop = FALSE
]


cat(
  "Baseline validation cohort:\n"
)


print(
  table(
    baseline_meta$status
  )
)


cat("\n")


cat(
  "Baseline n = ",
  nrow(
    baseline_meta
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 15. PRIMARY FROZEN SCORE
# COHORT-WIDE UNSUPERVISED Z-STANDARDIZATION
# ==============================================================================

z_main <- safe_z_rows(
  baseline_expr
)


score_main <- calculate_five_gene_score(
  z_main
)


# ==============================================================================
# 16. HEALTHY-REFERENCE SCALING SENSITIVITY
# ==============================================================================

healthy_samples <- baseline_meta %>%
  
  dplyr::filter(
    status ==
      "Hlty"
  ) %>%
  
  dplyr::pull(
    geo_accession
  )


if (
  length(
    healthy_samples
  ) < 5
) {
  
  stop(
    "Too few healthy-control samples for healthy-reference scaling."
  )
}


z_healthy_reference <- healthy_reference_z(
  
  baseline_expr,
  
  healthy_samples =
    healthy_samples
)


score_healthy_reference <- calculate_five_gene_score(
  z_healthy_reference
)


# ==============================================================================
# 17. SAMPLE-LEVEL VALIDATION TABLE
# ==============================================================================

scores_df <- baseline_meta %>%
  
  dplyr::select(
    geo_accession,
    sample_title,
    status,
    sex_standardized
  ) %>%
  
  dplyr::mutate(
    
    five_gene_score =
      unname(
        score_main[
          geo_accession
        ]
      ),
    
    five_gene_score_healthy_reference =
      unname(
        score_healthy_reference[
          geo_accession
        ]
      ),
    
    severity_order =
      dplyr::case_when(
        
        status ==
          "Hlty" ~
          0,
        
        status ==
          "Inf1_P" ~
          1,
        
        status ==
          "Seps_P" ~
          2,
        
        status ==
          "Shock_P" ~
          3,
        
        TRUE ~
          NA_real_
      )
  )


for (
  gene in five_genes
) {
  
  scores_df[[paste0(gene, "_expression")]] <-
    baseline_expr[
      gene,
      scores_df$geo_accession
    ]
  
  
  scores_df[[paste0(gene, "_z")]] <-
    z_main[
      gene,
      scores_df$geo_accession
    ]
}


cat(
  "Frozen five-gene scores calculated.\n\n"
)


# ==============================================================================
# 18. PRIMARY AND SECONDARY EXTERNAL COMPARISONS
# ==============================================================================

comparison_specs <- list(
  
  list(
    
    name =
      "Sepsis_or_shock_vs_uncomplicated",
    
    case =
      c(
        "Seps_P",
        "Shock_P"
      ),
    
    control =
      "Inf1_P"
  ),
  
  list(
    
    name =
      "Sepsis_or_shock_vs_healthy",
    
    case =
      c(
        "Seps_P",
        "Shock_P"
      ),
    
    control =
      "Hlty"
  ),
  
  list(
    
    name =
      "Sepsis_vs_uncomplicated",
    
    case =
      "Seps_P",
    
    control =
      "Inf1_P"
  ),
  
  list(
    
    name =
      "Shock_vs_uncomplicated",
    
    case =
      "Shock_P",
    
    control =
      "Inf1_P"
  ),
  
  list(
    
    name =
      "Shock_vs_sepsis",
    
    case =
      "Shock_P",
    
    control =
      "Seps_P"
  )
)


comparison_results <- dplyr::bind_rows(
  
  lapply(
    comparison_specs,
    function(spec) {
      
      compare_groups(
        
        df =
          scores_df,
        
        score_column =
          "five_gene_score",
        
        case_status =
          spec$case,
        
        control_status =
          spec$control,
        
        comparison_name =
          spec$name
      )
    }
  )
) %>%
  
  dplyr::mutate(
    
    p_BH_across_score_comparisons =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )


comparison_results_healthy_reference <- dplyr::bind_rows(
  
  lapply(
    comparison_specs,
    function(spec) {
      
      compare_groups(
        
        df =
          scores_df,
        
        score_column =
          "five_gene_score_healthy_reference",
        
        case_status =
          spec$case,
        
        control_status =
          spec$control,
        
        comparison_name =
          spec$name
      )
    }
  )
) %>%
  
  dplyr::mutate(
    
    p_BH_across_score_comparisons =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )


cat(
  "Frozen-score external comparisons:\n"
)


print(
  comparison_results,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 19. ORDERED BASELINE HOST-RESPONSE ANALYSIS
# ==============================================================================

ordered_spearman <- suppressWarnings(
  
  stats::cor.test(
    
    x =
      scores_df$five_gene_score,
    
    y =
      scores_df$severity_order,
    
    method =
      "spearman",
    
    exact =
      FALSE
  )
)


ordered_kruskal <- stats::kruskal.test(
  
  five_gene_score ~ status,
  
  data =
    scores_df
)


ordered_results <- tibble::tibble(
  
  analysis = c(
    "Spearman_score_vs_ordered_status",
    "Kruskal_Wallis_across_four_baseline_groups"
  ),
  
  statistic = c(
    
    unname(
      as.numeric(
        ordered_spearman$estimate
      )
    ),
    
    unname(
      as.numeric(
        ordered_kruskal$statistic
      )
    )
  ),
  
  statistic_name = c(
    "rho",
    "chi_squared"
  ),
  
  p_value = c(
    ordered_spearman$p.value,
    ordered_kruskal$p.value
  )
)


pairwise_wilcox <- stats::pairwise.wilcox.test(
  
  x =
    scores_df$five_gene_score,
  
  g =
    scores_df$status,
  
  p.adjust.method =
    "BH",
  
  exact =
    FALSE
)


pairwise_matrix <- as.data.frame(
  
  as.table(
    pairwise_wilcox$p.value
  )
)


names(
  pairwise_matrix
) <- c(
  "group_1",
  "group_2",
  "BH_adjusted_p"
)


pairwise_matrix <- pairwise_matrix %>%
  
  dplyr::filter(
    !is.na(
      BH_adjusted_p
    )
  )


cat(
  "Ordered baseline analysis:\n"
)


print(
  ordered_results,
  n = Inf
)


cat("\n")


# ==============================================================================
# 20. INDIVIDUAL FIVE-GENE DIRECTION AUDIT
#
# PRIMARY EXTERNAL CONTRAST:
# Seps_P + Shock_P versus Inf1_P
# ==============================================================================

primary_case_samples <- scores_df %>%
  
  dplyr::filter(
    status %in%
      c(
        "Seps_P",
        "Shock_P"
      )
  ) %>%
  
  dplyr::pull(
    geo_accession
  )


primary_control_samples <- scores_df %>%
  
  dplyr::filter(
    status ==
      "Inf1_P"
  ) %>%
  
  dplyr::pull(
    geo_accession
  )


gene_direction_results <- dplyr::bind_rows(
  
  lapply(
    five_genes,
    function(gene) {
      
      x_case <- baseline_expr[
        gene,
        primary_case_samples
      ]
      
      
      x_control <- baseline_expr[
        gene,
        primary_control_samples
      ]
      
      
      wt <- safe_wilcox(
        x_case,
        x_control
      )
      
      
      difference <- stats::median(
        x_case,
        na.rm = TRUE
      ) -
        stats::median(
          x_control,
          na.rm = TRUE
        )
      
      
      expected <-
        expected_direction[[gene]]
      
      
      observed <- dplyr::case_when(
        
        difference > 0 ~
          "UP",
        
        difference < 0 ~
          "DOWN",
        
        TRUE ~
          "NO_CHANGE"
      )
      
      
      tibble::tibble(
        
        gene =
          gene,
        
        expected_direction =
          expected,
        
        observed_direction =
          observed,
        
        direction_concordant =
          expected ==
          observed,
        
        n_case =
          sum(
            is.finite(
              x_case
            )
          ),
        
        n_control =
          sum(
            is.finite(
              x_control
            )
          ),
        
        case_median =
          stats::median(
            x_case,
            na.rm = TRUE
          ),
        
        control_median =
          stats::median(
            x_control,
            na.rm = TRUE
          ),
        
        median_difference_case_minus_control =
          difference,
        
        wilcoxon_W =
          wt$statistic,
        
        p_value =
          wt$p
      )
    }
  )
) %>%
  
  dplyr::mutate(
    
    p_BH_five_genes =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )


cat(
  "Individual-gene direction audit:\n"
)


print(
  gene_direction_results,
  n = Inf,
  width = Inf
)


cat("\n")


n_direction_concordant <- sum(
  gene_direction_results$direction_concordant,
  na.rm = TRUE
)


cat(
  "Directionally concordant frozen genes: ",
  n_direction_concordant,
  "/5\n\n",
  sep = ""
)


# ==============================================================================
# 21. EXPRESSION AND SCORE GROUP SUMMARIES
# ==============================================================================

expression_long <- scores_df %>%
  
  dplyr::select(
    
    geo_accession,
    
    sample_title,
    
    status,
    
    dplyr::ends_with(
      "_expression"
    )
  ) %>%
  
  tidyr::pivot_longer(
    
    cols =
      dplyr::ends_with(
        "_expression"
      ),
    
    names_to =
      "gene",
    
    values_to =
      "expression"
  ) %>%
  
  dplyr::mutate(
    
    gene =
      stringr::str_remove(
        gene,
        "_expression$"
      ),
    
    gene =
      factor(
        gene,
        levels =
          five_genes
      )
  )


expression_group_summary <- expression_long %>%
  
  dplyr::group_by(
    gene,
    status
  ) %>%
  
  dplyr::summarise(
    
    n =
      sum(
        is.finite(
          expression
        )
      ),
    
    median =
      stats::median(
        expression,
        na.rm = TRUE
      ),
    
    mean =
      mean(
        expression,
        na.rm = TRUE
      ),
    
    sd =
      stats::sd(
        expression,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


score_group_summary <- scores_df %>%
  
  dplyr::group_by(
    status
  ) %>%
  
  dplyr::summarise(
    
    n =
      sum(
        is.finite(
          five_gene_score
        )
      ),
    
    median =
      stats::median(
        five_gene_score,
        na.rm = TRUE
      ),
    
    q1 =
      stats::quantile(
        five_gene_score,
        0.25,
        na.rm = TRUE
      ),
    
    q3 =
      stats::quantile(
        five_gene_score,
        0.75,
        na.rm = TRUE
      ),
    
    mean =
      mean(
        five_gene_score,
        na.rm = TRUE
      ),
    
    sd =
      stats::sd(
        five_gene_score,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


cat(
  "Five-gene score by baseline group:\n"
)


print(
  score_group_summary,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 22. SCORE-SCALING SENSITIVITY
# ==============================================================================

scaling_cor <- suppressWarnings(
  
  stats::cor.test(
    
    x =
      scores_df$five_gene_score,
    
    y =
      scores_df$five_gene_score_healthy_reference,
    
    method =
      "spearman",
    
    exact =
      FALSE
  )
)


scaling_sensitivity <- tibble::tibble(
  
  comparison =
    "cohort_standardized_vs_healthy_reference",
  
  n =
    sum(
      complete.cases(
        scores_df$five_gene_score,
        scores_df$five_gene_score_healthy_reference
      )
    ),
  
  spearman_rho =
    unname(
      as.numeric(
        scaling_cor$estimate
      )
    ),
  
  p_value =
    scaling_cor$p.value
)


cat(
  "Score-scaling sensitivity:\n"
)


print(
  scaling_sensitivity,
  n = Inf
)


cat("\n")


# ==============================================================================
# 23. PRIMARY EXTERNAL ROC
# ==============================================================================

primary_df <- scores_df %>%
  
  dplyr::filter(
    status %in%
      c(
        "Inf1_P",
        "Seps_P",
        "Shock_P"
      )
  ) %>%
  
  dplyr::mutate(
    
    primary_group =
      ifelse(
        status ==
          "Inf1_P",
        "Uncomplicated infection",
        "Sepsis / septic shock"
      ),
    
    primary_group =
      factor(
        primary_group,
        levels = c(
          "Uncomplicated infection",
          "Sepsis / septic shock"
        )
      )
  )


primary_roc <- pROC::roc(
  
  response =
    primary_df$primary_group,
  
  predictor =
    primary_df$five_gene_score,
  
  levels = c(
    "Uncomplicated infection",
    "Sepsis / septic shock"
  ),
  
  direction =
    "<",
  
  quiet =
    TRUE
)


primary_roc_ci <- suppressWarnings(
  
  pROC::ci.auc(
    primary_roc,
    method = "delong"
  )
)


primary_auc <- as.numeric(
  pROC::auc(
    primary_roc
  )
)


cat(
  "PRIMARY EXTERNAL VALIDATION:\n"
)


cat(
  "Sepsis/septic shock vs uncomplicated infection\n"
)


cat(
  "AUC = ",
  signif(
    primary_auc,
    4
  ),
  "\n",
  sep = ""
)


cat(
  "95% CI = ",
  signif(
    primary_roc_ci[1],
    4
  ),
  " - ",
  signif(
    primary_roc_ci[3],
    4
  ),
  "\n\n",
  sep = ""
)


roc_coords <- as.data.frame(
  
  pROC::coords(
    
    primary_roc,
    
    x =
      "all",
    
    ret = c(
      "specificity",
      "sensitivity"
    ),
    
    transpose =
      FALSE
  )
) %>%
  
  dplyr::mutate(
    
    false_positive_rate =
      1 -
      specificity
  )


# ==============================================================================
# 24. PUBLICATION COLOR PALETTE
# ==============================================================================

status_colors <- c(
  
  "Hlty" =
    "#56B4E9",
  
  "Inf1_P" =
    "#7F8C8D",
  
  "Seps_P" =
    "#E69F00",
  
  "Shock_P" =
    "#D55E00"
)


direction_colors <- c(
  
  "UP" =
    "#D55E00",
  
  "DOWN" =
    "#0072B2"
)


# ==============================================================================
# 25. FIGURE A
# FROZEN SCORE ACROSS BASELINE GROUPS
# ==============================================================================

pA <- ggplot2::ggplot(
  
  scores_df,
  
  ggplot2::aes(
    x =
      status,
    y =
      five_gene_score,
    fill =
      status,
    color =
      status
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.62,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.4,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      status_colors
  ) +
  
  ggplot2::scale_color_manual(
    values =
      status_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Frozen five-gene host-response score in GSE154918",
    
    subtitle =
      "Independent whole-blood RNA-seq cohort",
    
    x =
      NULL,
    
    y =
      "Five-gene score"
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
# 26. FIGURE B
# COMPONENT GENES ACROSS BASELINE GROUPS
# ==============================================================================

pB <- ggplot2::ggplot(
  
  expression_long,
  
  ggplot2::aes(
    x =
      status,
    y =
      expression,
    fill =
      status,
    color =
      status
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.62,
    alpha = 0.50,
    outlier.shape = NA
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    size = 1.4,
    alpha = 0.72
  ) +
  
  ggplot2::facet_wrap(
    ~ gene,
    nrow = 1,
    scales = "free_y"
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      status_colors
  ) +
  
  ggplot2::scale_color_manual(
    values =
      status_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Expression of the five frozen component genes",
    
    x =
      NULL,
    
    y =
      "GEO processed expression"
  ) +
  
  ggplot2::theme_bw(
    base_size = 12
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    strip.text =
      ggplot2::element_text(
        face = "bold"
      ),
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 27. FIGURE C
# PRIMARY EXTERNAL ROC
# ==============================================================================

auc_label <- paste0(
  
  "AUC = ",
  sprintf(
    "%.3f",
    primary_auc
  ),
  
  "\n95% CI ",
  
  sprintf(
    "%.3f",
    as.numeric(
      primary_roc_ci[1]
    )
  ),
  
  "-",
  
  sprintf(
    "%.3f",
    as.numeric(
      primary_roc_ci[3]
    )
  )
)


pC <- ggplot2::ggplot(
  
  roc_coords,
  
  ggplot2::aes(
    x =
      false_positive_rate,
    y =
      sensitivity
  )
) +
  
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  ggplot2::geom_line(
    linewidth = 1.2,
    color = "#D55E00"
  ) +
  
  ggplot2::coord_equal() +
  
  ggplot2::annotate(
    "text",
    x = 0.58,
    y = 0.18,
    label = auc_label,
    hjust = 0,
    size = 4.5
  ) +
  
  ggplot2::labs(
    
    title =
      "Primary external comparison",
    
    subtitle =
      "Sepsis / septic shock versus uncomplicated infection",
    
    x =
      "1 - Specificity",
    
    y =
      "Sensitivity"
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
# 28. FIGURE D
# SCALING SENSITIVITY
# ==============================================================================

pD <- ggplot2::ggplot(
  
  scores_df,
  
  ggplot2::aes(
    x =
      five_gene_score,
    y =
      five_gene_score_healthy_reference,
    color =
      status
  )
) +
  
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey35",
    fill = "grey85",
    linewidth = 0.8
  ) +
  
  ggplot2::scale_color_manual(
    values =
      status_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Score robustness to scaling strategy",
    
    subtitle =
      paste0(
        "Spearman rho = ",
        sprintf(
          "%.3f",
          unname(
            as.numeric(
              scaling_cor$estimate
            )
          )
        ),
        "; p = ",
        format.pval(
          scaling_cor$p.value,
          digits = 3
        )
      ),
    
    x =
      "Baseline cohort-standardized score",
    
    y =
      "Healthy-reference standardized score",
    
    color =
      "Clinical group"
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
# 29. FIGURE E
# EXTERNAL DIRECTION AUDIT
# ==============================================================================

pE_data <- gene_direction_results %>%
  
  dplyr::mutate(
    
    gene =
      factor(
        gene,
        levels =
          rev(
            five_genes
          )
      )
  )


pE <- ggplot2::ggplot(
  
  pE_data,
  
  ggplot2::aes(
    x =
      median_difference_case_minus_control,
    y =
      gene,
    color =
      expected_direction
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend =
        median_difference_case_minus_control,
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
      "External direction audit of the frozen five-gene panel",
    
    subtitle =
      paste0(
        "Median expression difference: sepsis/septic shock ",
        "minus uncomplicated infection"
      ),
    
    x =
      "Median expression difference",
    
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
# 30. SAVE FIGURES
# ==============================================================================

save_plot_all_formats(
  
  plot_object =
    pA,
  
  filename_base =
    "141_Figure_A_five_gene_score_baseline_groups",
  
  width =
    9,
  
  height =
    7
)


save_plot_all_formats(
  
  plot_object =
    pB,
  
  filename_base =
    "141_Figure_B_five_component_genes_baseline_groups",
  
  width =
    16,
  
  height =
    6.5
)


save_plot_all_formats(
  
  plot_object =
    pC,
  
  filename_base =
    "141_Figure_C_primary_external_ROC",
  
  width =
    8,
  
  height =
    7
)


save_plot_all_formats(
  
  plot_object =
    pD,
  
  filename_base =
    "141_Figure_D_score_scaling_sensitivity",
  
  width =
    8.5,
  
  height =
    7
)


save_plot_all_formats(
  
  plot_object =
    pE,
  
  filename_base =
    "141_Supplementary_Figure_gene_direction_audit",
  
  width =
    9,
  
  height =
    6
)


combined_figure <- (
  
  pA +
    pC
  
) / (
  
  pD +
    pE
  
) +
  
  patchwork::plot_annotation(
    
    title =
      paste0(
        "External validation of the frozen ",
        "five-gene host-response signature"
      ),
    
    subtitle =
      "GSE154918 whole-blood RNA-seq"
  )


save_plot_all_formats(
  
  plot_object =
    combined_figure,
  
  filename_base =
    "141_Figure_external_validation_integrated",
  
  width =
    16,
  
  height =
    12
)


cat(
  "Publication-quality color figures saved.\n\n"
)


# ==============================================================================
# 31. EXPORT TABLES
# ==============================================================================

write.csv(
  
  scores_df,
  
  file.path(
    tables_dir,
    "141_GSE154918_five_gene_scores.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  comparison_results,
  
  file.path(
    tables_dir,
    "141_external_score_comparisons.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  comparison_results_healthy_reference,
  
  file.path(
    tables_dir,
    "141_external_score_comparisons_healthy_reference.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  gene_direction_results,
  
  file.path(
    tables_dir,
    "141_external_five_gene_direction_audit.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  ordered_results,
  
  file.path(
    tables_dir,
    "141_ordered_baseline_status_tests.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  pairwise_matrix,
  
  file.path(
    tables_dir,
    "141_pairwise_baseline_group_tests.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  score_group_summary,
  
  file.path(
    tables_dir,
    "141_score_group_summary.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  expression_group_summary,
  
  file.path(
    tables_dir,
    "141_gene_expression_group_summary.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  gene_coverage,
  
  file.path(
    tables_dir,
    "141_frozen_gene_coverage.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  gene_mapping_summary,
  
  file.path(
    tables_dir,
    "141_gene_mapping_summary.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  scaling_sensitivity,
  
  file.path(
    tables_dir,
    "141_score_scaling_sensitivity.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 32. EXCEL WORKBOOK
# ==============================================================================

run_info <- tibble::tibble(
  
  parameter = c(
    "script",
    "run_date",
    "GEO_accession",
    "data_type",
    "primary_panel",
    "score_definition",
    "primary_external_comparison",
    "feature_selection_external",
    "coefficient_refitting_external",
    "cutoff_optimization_external",
    "direction_flipping_external",
    "followup_in_primary_analysis"
  ),
  
  value = c(
    
    script_name,
    
    as.character(
      run_date
    ),
    
    gse_id,
    
    "whole-blood RNA-seq; GEO processed expression",
    
    paste(
      five_genes,
      collapse = "; "
    ),
    
    paste0(
      "mean z(CD177,HK3,IRAK3) - ",
      "mean z(CARD11,IKZF2)"
    ),
    
    "Seps_P + Shock_P versus Inf1_P",
    
    "NO",
    
    "NO",
    
    "NO",
    
    "NO",
    
    "NO"
  )
)


workbook_file <- file.path(
  tables_dir,
  "141_GSE154918_external_validation.xlsx"
)


wb <- openxlsx::createWorkbook()


sheet_list <- list(
  
  "00_run_info" =
    run_info,
  
  "01_status_count_check" =
    status_count_check,
  
  "02_metadata" =
    metadata,
  
  "03_gene_coverage" =
    gene_coverage,
  
  "04_sample_scores" =
    scores_df,
  
  "05_score_comparisons" =
    comparison_results,
  
  "06_healthyref_comparisons" =
    comparison_results_healthy_reference,
  
  "07_gene_direction_audit" =
    gene_direction_results,
  
  "08_ordered_tests" =
    ordered_results,
  
  "09_pairwise_tests" =
    pairwise_matrix,
  
  "10_score_group_summary" =
    score_group_summary,
  
  "11_gene_group_summary" =
    expression_group_summary,
  
  "12_scaling_sensitivity" =
    scaling_sensitivity,
  
  "13_gene_mapping_summary" =
    gene_mapping_summary,
  
  "14_sample_column_mapping" =
    sample_column_mapping
)


header_style <- openxlsx::createStyle(
  
  textDecoration =
    "bold",
  
  halign =
    "center",
  
  valign =
    "center",
  
  border =
    "Bottom"
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
    
    sheet =
      sheet_name,
    
    x =
      sheet_list[[sheet_name]],
    
    headerStyle =
      header_style
  )
  
  
  openxlsx::freezePane(
    
    wb,
    
    sheet =
      sheet_name,
    
    firstRow =
      TRUE
  )
  
  
  openxlsx::setColWidths(
    
    wb,
    
    sheet =
      sheet_name,
    
    cols =
      seq_len(
        max(
          1,
          ncol(
            sheet_list[[sheet_name]]
          )
        )
      ),
    
    widths =
      "auto"
  )
}


openxlsx::saveWorkbook(
  
  wb,
  
  workbook_file,
  
  overwrite =
    TRUE
)


# ==============================================================================
# 33. AUTOMATED MANUSCRIPT SUMMARY
# ==============================================================================

primary_row <- comparison_results %>%
  
  dplyr::filter(
    comparison ==
      "Sepsis_or_shock_vs_uncomplicated"
  )


healthy_row <- comparison_results %>%
  
  dplyr::filter(
    comparison ==
      "Sepsis_or_shock_vs_healthy"
  )


if (
  nrow(
    primary_row
  ) != 1
) {
  
  stop(
    "Primary external comparison was not uniquely identified."
  )
}


if (
  nrow(
    healthy_row
  ) != 1
) {
  
  stop(
    "Healthy-control comparison was not uniquely identified."
  )
}


ordered_rho <- ordered_results %>%
  
  dplyr::filter(
    analysis ==
      "Spearman_score_vs_ordered_status"
  ) %>%
  
  dplyr::pull(
    statistic
  )


ordered_p <- ordered_results %>%
  
  dplyr::filter(
    analysis ==
      "Spearman_score_vs_ordered_status"
  ) %>%
  
  dplyr::pull(
    p_value
  )


summary_en <- c(
  
  "GSE154918 EXTERNAL VALIDATION - MANUSCRIPT SUMMARY",
  
  "====================================================================",
  
  "",
  
  paste0(
    "The frozen five-gene host-response signature ",
    "(CD177, HK3, IRAK3, CARD11, and IKZF2) was evaluated ",
    "in GSE154918 without feature selection, coefficient refitting, ",
    "score-direction flipping, or cutoff optimization."
  ),
  
  "",
  
  paste0(
    "Primary external comparison: baseline sepsis/septic shock ",
    "versus uncomplicated infection; n_case = ",
    primary_row$n_case,
    ", n_control = ",
    primary_row$n_control,
    "."
  ),
  
  paste0(
    "Median five-gene score: case = ",
    signif(
      primary_row$case_median,
      4
    ),
    "; control = ",
    signif(
      primary_row$control_median,
      4
    ),
    "."
  ),
  
  paste0(
    "Wilcoxon p = ",
    format.pval(
      primary_row$p_value,
      digits = 4
    ),
    "."
  ),
  
  paste0(
    "Fixed-direction AUC = ",
    sprintf(
      "%.3f",
      primary_row$auc_fixed_direction
    ),
    " (95% CI ",
    sprintf(
      "%.3f",
      primary_row$auc_ci_low
    ),
    "-",
    sprintf(
      "%.3f",
      primary_row$auc_ci_high
    ),
    ")."
  ),
  
  "",
  
  paste0(
    "Secondary sepsis/septic shock versus healthy-control AUC = ",
    sprintf(
      "%.3f",
      healthy_row$auc_fixed_direction
    ),
    "."
  ),
  
  "",
  
  paste0(
    n_direction_concordant,
    " of 5 frozen component genes showed the expected direction ",
    "in the primary external comparison."
  ),
  
  "",
  
  paste0(
    "Across the ordered baseline states Hlty -> Inf1_P -> ",
    "Seps_P -> Shock_P, the five-gene score showed Spearman rho = ",
    sprintf(
      "%.3f",
      ordered_rho
    ),
    ", p = ",
    format.pval(
      ordered_p,
      digits = 4
    ),
    "."
  ),
  
  "",
  
  paste0(
    "Cohort-standardized and healthy-reference implementations ",
    "were concordant (Spearman rho = ",
    sprintf(
      "%.3f",
      scaling_sensitivity$spearman_rho
    ),
    ")."
  ),
  
  "",
  
  "INTERPRETATION:",
  
  paste0(
    "This analysis evaluates external transcriptomic replication ",
    "of a frozen five-gene composition and score direction. ",
    "It does not constitute validation of a pre-calibrated clinical ",
    "assay or clinical decision threshold."
  )
)


summary_ru <- c(
  
  "GSE154918 - ВНЕШНЯЯ ВАЛИДАЦИЯ: РЕЗЮМЕ",
  
  "====================================================================",
  
  "",
  
  paste0(
    "Замороженная пятигенная host-response сигнатура ",
    "(CD177, HK3, IRAK3, CARD11 и IKZF2) была проверена ",
    "в GSE154918 без нового отбора признаков, изменения коэффициентов, ",
    "изменения направления score или оптимизации диагностического порога."
  ),
  
  "",
  
  paste0(
    "Основное внешнее сравнение: исходные sepsis/septic shock ",
    "против uncomplicated infection; n_case = ",
    primary_row$n_case,
    ", n_control = ",
    primary_row$n_control,
    "."
  ),
  
  paste0(
    "Медиана пятигенного score: sepsis/shock = ",
    signif(
      primary_row$case_median,
      4
    ),
    "; uncomplicated infection = ",
    signif(
      primary_row$control_median,
      4
    ),
    "."
  ),
  
  paste0(
    "Wilcoxon p = ",
    format.pval(
      primary_row$p_value,
      digits = 4
    ),
    "."
  ),
  
  paste0(
    "AUC с заранее фиксированным направлением = ",
    sprintf(
      "%.3f",
      primary_row$auc_fixed_direction
    ),
    " (95% ДИ ",
    sprintf(
      "%.3f",
      primary_row$auc_ci_low
    ),
    "-",
    sprintf(
      "%.3f",
      primary_row$auc_ci_high
    ),
    ")."
  ),
  
  "",
  
  paste0(
    n_direction_concordant,
    " из 5 компонентов панели сохранили ожидаемое направление ",
    "в основном внешнем сравнении."
  ),
  
  "",
  
  paste0(
    "В последовательности Hlty -> Inf1_P -> Seps_P -> Shock_P ",
    "корреляция пятигенного score с порядком групп составила rho = ",
    sprintf(
      "%.3f",
      ordered_rho
    ),
    ", p = ",
    format.pval(
      ordered_p,
      digits = 4
    ),
    "."
  ),
  
  "",
  
  "ИНТЕРПРЕТАЦИЯ:",
  
  paste0(
    "Этот анализ представляет внешнюю транскриптомную проверку ",
    "заранее замороженного состава и направления пятигенной сигнатуры. ",
    "Он не является валидацией готового клинического теста или ",
    "заранее калиброванного диагностического порога."
  )
)


summary_en_file <- file.path(
  text_dir,
  "141_external_validation_summary_EN.txt"
)


summary_ru_file <- file.path(
  text_dir,
  "141_external_validation_summary_RU.txt"
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
# 34. INPUT / DOWNLOAD MANIFEST
# ==============================================================================

manifest <- tibble::tibble(
  
  item = c(
    "predeclared_plan",
    "GEO_supplementary_expression_file"
  ),
  
  path = c(
    plan_file,
    supp_file
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
    "141_input_download_manifest.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 35. SESSION INFO
# ==============================================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      logs_dir,
      "141_sessionInfo.txt"
    )
)


# ==============================================================================
# 36. FINAL REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 141 completed successfully.\n")
cat("====================================================================\n\n")


cat(
  "GSE154918 BASELINE COHORT:\n"
)


print(
  table(
    scores_df$status
  )
)


cat("\n")


cat(
  "FROZEN PANEL:\n"
)


cat(
  "UP:   ",
  paste(
    up_genes,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "DOWN: ",
  paste(
    down_genes,
    collapse = ", "
  ),
  "\n\n",
  sep = ""
)


cat(
  "PRIMARY EXTERNAL COMPARISON:\n"
)


cat(
  "Seps_P + Shock_P vs Inf1_P\n"
)


cat(
  "n = ",
  primary_row$n_case,
  " vs ",
  primary_row$n_control,
  "\n",
  sep = ""
)


cat(
  "Median score = ",
  signif(
    primary_row$case_median,
    4
  ),
  " vs ",
  signif(
    primary_row$control_median,
    4
  ),
  "\n",
  sep = ""
)


cat(
  "Wilcoxon p = ",
  format.pval(
    primary_row$p_value,
    digits = 4
  ),
  "\n",
  sep = ""
)


cat(
  "AUC = ",
  sprintf(
    "%.3f",
    primary_row$auc_fixed_direction
  ),
  " [",
  sprintf(
    "%.3f",
    primary_row$auc_ci_low
  ),
  ", ",
  sprintf(
    "%.3f",
    primary_row$auc_ci_high
  ),
  "]\n\n",
  sep = ""
)


cat(
  "Expected individual-gene direction: ",
  n_direction_concordant,
  "/5\n\n",
  sep = ""
)


cat(
  "ORDERED BASELINE GRADIENT:\n"
)


cat(
  "Hlty -> Inf1_P -> Seps_P -> Shock_P\n"
)


cat(
  "Spearman rho = ",
  sprintf(
    "%.3f",
    ordered_rho
  ),
  "; p = ",
  format.pval(
    ordered_p,
    digits = 4
  ),
  "\n\n",
  sep = ""
)


cat(
  "SCALING SENSITIVITY:\n"
)


cat(
  "Cohort-z vs healthy-reference rho = ",
  sprintf(
    "%.3f",
    scaling_sensitivity$spearman_rho
  ),
  "\n\n",
  sep = ""
)


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
  "Main integrated color figure:\n"
)


cat(
  normalizePath(
    file.path(
      figures_dir,
      "141_Figure_external_validation_integrated.png"
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
  "- Frozen five-gene composition unchanged.\n"
)


cat(
  "- No feature selection in GSE154918.\n"
)


cat(
  "- No coefficient refitting.\n"
)


cat(
  "- No cutoff optimization.\n"
)


cat(
  "- No post hoc score-direction flipping.\n"
)


cat(
  "- Primary comparison was declared before expression analysis.\n"
)


cat(
  "- Follow-up samples excluded from primary analysis.\n"
)


cat(
  paste0(
    "- This is external transcriptomic replication, ",
    "not clinical assay validation.\n\n"
  )
)


cat(
  "Done.\n"
)