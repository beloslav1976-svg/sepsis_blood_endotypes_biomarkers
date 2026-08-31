# ==============================================================================
# Script 143
# Multicohort integration of the frozen five-gene host-response signature
#
# Project: Sepsis_DESeq2
#
# PURPOSE
# Integrate:
#
#   DISCOVERY / INTERNAL EVIDENCE
#     Script 135:
#       - frozen five-gene score
#       - SRS
#       - SRSq
#       - CTS
#
#     Script 136b:
#       - demographic sensitivity result
#
#     Script 137:
#       - convergence with previously published transcriptomic signatures
#
#   EXTERNAL EVIDENCE
#     Script 141:
#       - GSE154918
#       - independent whole-blood RNA-seq
#       - directional replication
#       - clinical-state gradient
#
#     Script 142b:
#       - GSE185263
#       - independent whole-blood RNA-seq
#       - primary score-SOFA association
#       - individual gene-SOFA associations
#       - covariate-adjusted association
#       - geographic replication
#
#
# FROZEN PANEL
#   UP:
#     CD177
#     HK3
#     IRAK3
#
#   DOWN:
#     CARD11
#     IKZF2
#
#
# CENTRAL BIOLOGICAL INTERPRETATION
#
# The five-gene signature is treated as a compact molecular representation
# of a myeloid-adaptive host-response axis.
#
# It is NOT treated as a fully calibrated clinical diagnostic or prognostic
# assay.
#
#
# IMPORTANT
#
#   - NO new feature selection
#   - NO gene substitution
#   - NO gene addition/deletion
#   - NO coefficient refitting
#   - NO cutoff optimization
#   - NO post hoc score-direction flipping
#
#   - GSE154918 and GSE185263 have different prespecified endpoints.
#     They are NOT pooled into one cross-dataset effect estimate.
#
#
# NEW ANALYSIS IN SCRIPT 143
#
# For the five geographic subcohorts within GSE185263:
#
#   - Fisher-z transformation of Spearman correlations
#   - fixed-effect pooled correlation
#   - Cochran Q
#   - heterogeneity p value
#   - I^2
#   - DerSimonian-Laird tau^2
#   - random-effects pooled correlation
#
#
# MAIN OUTPUTS
#
#   tables/
#     143_multicohort_five_gene_evidence_package.xlsx
#     143_manuscript_key_numbers.csv
#     143_multicohort_evidence_summary.csv
#
#   figures/
#     143_MAIN_Figure_multicohort_five_gene_evidence.png
#     + individual publication panels
#
#   text/
#     143_multicohort_evidence_summary_EN.txt
#     143_multicohort_evidence_summary_RU.txt
#     143_proposed_external_validation_results_paragraph_EN.txt
#
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

script_name <-
  "143_multicohort_five_gene_evidence_integration.R"

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
cat("Running Script 143\n")
cat("Multicohort five-gene evidence integration\n")
cat("====================================================================\n\n")


cat(
  "Project directory:\n"
)


cat(
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "Run date:\n"
)


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
  "143_multicohort_integration"
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


cat(
  "Output folder:\n"
)


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
  "ggrepel",
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


if (
  length(
    missing_cran
  ) > 0
) {
  
  install.packages(
    missing_cran
  )
}


suppressPackageStartupMessages({
  
  library(data.table)
  
  library(dplyr)
  
  library(tidyr)
  
  library(stringr)
  
  library(tibble)
  
  library(ggplot2)
  
  library(ggrepel)
  
  library(patchwork)
  
  library(openxlsx)
})


cat(
  "Required packages loaded successfully.\n\n"
)


# ==============================================================================
# 3. FROZEN PANEL
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


discovery_direction <- c(
  
  CD177 =
    "POSITIVE",
  
  HK3 =
    "POSITIVE",
  
  IRAK3 =
    "POSITIVE",
  
  CARD11 =
    "NEGATIVE",
  
  IKZF2 =
    "NEGATIVE"
)


# ==============================================================================
# 4. HELPER FUNCTIONS
# ==============================================================================

read_csv_fast <- function(
    path
) {
  
  data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )
}


find_column <- function(
    df,
    patterns,
    exclude_pattern = NULL,
    label = "column"
) {
  
  nms <- names(
    df
  )
  
  
  for (
    pattern in patterns
  ) {
    
    hits <- nms[
      grepl(
        pattern,
        nms,
        ignore.case = TRUE,
        perl = TRUE
      )
    ]
    
    
    if (!is.null(
      exclude_pattern
    )) {
      
      hits <- hits[
        !grepl(
          exclude_pattern,
          hits,
          ignore.case = TRUE,
          perl = TRUE
        )
      ]
    }
    
    
    if (
      length(
        hits
      ) > 0
    ) {
      
      return(
        hits[1]
      )
    }
  }
  
  
  stop(
    paste0(
      "Unable to detect ",
      label,
      ".\nAvailable columns:\n",
      paste(
        nms,
        collapse = ", "
      )
    )
  )
}


find_csv_by_required_columns <- function(
    directory,
    required_columns,
    recursive = TRUE,
    label = "CSV file"
) {
  
  files <- list.files(
    directory,
    pattern = "\\.csv$",
    recursive = recursive,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  
  if (
    length(
      files
    ) == 0
  ) {
    
    stop(
      paste0(
        "No CSV files found in: ",
        directory
      )
    )
  }
  
  
  matched_files <- character()
  
  
  for (
    file_i in files
  ) {
    
    header_i <- tryCatch(
      
      data.table::fread(
        file_i,
        nrows = 0,
        data.table = FALSE,
        check.names = FALSE
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    if (is.null(
      header_i
    )) {
      
      next
    }
    
    
    if (
      all(
        required_columns %in%
        names(
          header_i
        )
      )
    ) {
      
      matched_files <- c(
        matched_files,
        file_i
      )
    }
  }
  
  
  if (
    length(
      matched_files
    ) == 0
  ) {
    
    stop(
      paste0(
        "Unable to locate ",
        label,
        " containing columns: ",
        paste(
          required_columns,
          collapse = ", "
        ),
        "\nSearch directory: ",
        directory
      )
    )
  }
  
  
  if (
    length(
      matched_files
    ) > 1
  ) {
    
    cat(
      "Multiple candidate files found for ",
      label,
      ":\n",
      sep = ""
    )
    
    print(
      matched_files
    )
    
    cat(
      "Using first candidate:\n",
      matched_files[1],
      "\n\n",
      sep = ""
    )
  }
  
  
  return(
    matched_files[1]
  )
}


safe_spearman <- function(
    x,
    y
) {
  
  x <- suppressWarnings(
    as.numeric(
      as.character(
        x
      )
    )
  )
  
  
  y <- suppressWarnings(
    as.numeric(
      as.character(
        y
      )
    )
  )
  
  
  keep <- is.finite(
    x
  ) &
    is.finite(
      y
    )
  
  
  x <- x[
    keep
  ]
  
  
  y <- y[
    keep
  ]
  
  
  if (
    length(
      x
    ) < 5
  ) {
    
    return(
      tibble::tibble(
        
        n =
          length(
            x
          ),
        
        rho =
          NA_real_,
        
        p_value =
          NA_real_
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
  
  
  tibble::tibble(
    
    n =
      length(
        x
      ),
    
    rho =
      unname(
        as.numeric(
          test$estimate
        )
      ),
    
    p_value =
      test$p.value
  )
}


kruskal_epsilon2 <- function(
    value,
    group
) {
  
  value <- suppressWarnings(
    as.numeric(
      as.character(
        value
      )
    )
  )
  
  
  group <- as.character(
    group
  )
  
  
  keep <- is.finite(
    value
  ) &
    !is.na(
      group
    ) &
    group != ""
  
  
  value <- value[
    keep
  ]
  
  
  group <- factor(
    group[
      keep
    ]
  )
  
  
  if (
    nlevels(
      group
    ) < 2
  ) {
    
    return(
      tibble::tibble(
        
        n =
          length(
            value
          ),
        
        groups =
          nlevels(
            group
          ),
        
        H =
          NA_real_,
        
        p_value =
          NA_real_,
        
        epsilon2 =
          NA_real_
      )
    )
  }
  
  
  test <- stats::kruskal.test(
    value ~ group
  )
  
  
  n <- length(
    value
  )
  
  
  k <- nlevels(
    group
  )
  
  
  H <- unname(
    as.numeric(
      test$statistic
    )
  )
  
  
  epsilon2 <- (
    H - k + 1
  ) / (
    n - k
  )
  
  
  epsilon2 <- max(
    0,
    epsilon2
  )
  
  
  tibble::tibble(
    
    n =
      n,
    
    groups =
      k,
    
    H =
      H,
    
    p_value =
      test$p.value,
    
    epsilon2 =
      epsilon2
  )
}


correlation_ci_fisher <- function(
    rho,
    n,
    conf_level = 0.95
) {
  
  if (
    !is.finite(
      rho
    ) ||
    !is.finite(
      n
    ) ||
    n <= 3 ||
    abs(
      rho
    ) >= 1
  ) {
    
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }
  
  
  alpha <- 1 -
    conf_level
  
  
  zcrit <- stats::qnorm(
    1 -
      alpha /
      2
  )
  
  
  z <- atanh(
    rho
  )
  
  
  se <- 1 /
    sqrt(
      n - 3
    )
  
  
  c(
    
    lower =
      tanh(
        z -
          zcrit *
          se
      ),
    
    upper =
      tanh(
        z +
          zcrit *
          se
      )
  )
}


meta_correlations_fisher <- function(
    rho,
    n
) {
  
  keep <- is.finite(
    rho
  ) &
    is.finite(
      n
    ) &
    n > 3 &
    abs(
      rho
    ) < 1
  
  
  rho <- rho[
    keep
  ]
  
  
  n <- n[
    keep
  ]
  
  
  k <- length(
    rho
  )
  
  
  if (
    k < 2
  ) {
    
    stop(
      "At least two valid correlations are required for meta-analysis."
    )
  }
  
  
  fisher_z <- atanh(
    rho
  )
  
  
  variance <- 1 /
    (
      n - 3
    )
  
  
  weight_fixed <- 1 /
    variance
  
  
  pooled_z_fixed <- sum(
    weight_fixed *
      fisher_z
  ) /
    sum(
      weight_fixed
    )
  
  
  se_fixed <- sqrt(
    1 /
      sum(
        weight_fixed
      )
  )
  
  
  Q <- sum(
    weight_fixed *
      (
        fisher_z -
          pooled_z_fixed
      )^2
  )
  
  
  Q_df <- k -
    1
  
  
  Q_p <- stats::pchisq(
    Q,
    df = Q_df,
    lower.tail = FALSE
  )
  
  
  C_value <- sum(
    weight_fixed
  ) -
    (
      sum(
        weight_fixed^2
      ) /
        sum(
          weight_fixed
        )
    )
  
  
  tau2 <- max(
    0,
    (
      Q -
        Q_df
    ) /
      C_value
  )
  
  
  I2 <- if (
    Q > 0
  ) {
    
    max(
      0,
      (
        Q -
          Q_df
      ) /
        Q
    ) *
      100
    
  } else {
    
    0
  }
  
  
  weight_random <- 1 /
    (
      variance +
        tau2
    )
  
  
  pooled_z_random <- sum(
    weight_random *
      fisher_z
  ) /
    sum(
      weight_random
    )
  
  
  se_random <- sqrt(
    1 /
      sum(
        weight_random
      )
  )
  
  
  zcrit <- stats::qnorm(
    0.975
  )
  
  
  fixed_low <- pooled_z_fixed -
    zcrit *
    se_fixed
  
  
  fixed_high <- pooled_z_fixed +
    zcrit *
    se_fixed
  
  
  random_low <- pooled_z_random -
    zcrit *
    se_random
  
  
  random_high <- pooled_z_random +
    zcrit *
    se_random
  
  
  fixed_p <- 2 *
    stats::pnorm(
      -abs(
        pooled_z_fixed /
          se_fixed
      )
    )
  
  
  random_p <- 2 *
    stats::pnorm(
      -abs(
        pooled_z_random /
          se_random
      )
    )
  
  
  tibble::tibble(
    
    k =
      k,
    
    total_n =
      sum(
        n
      ),
    
    Q =
      Q,
    
    Q_df =
      Q_df,
    
    Q_p =
      Q_p,
    
    I2_percent =
      I2,
    
    tau2_fisher_z =
      tau2,
    
    fixed_rho =
      tanh(
        pooled_z_fixed
      ),
    
    fixed_CI_low =
      tanh(
        fixed_low
      ),
    
    fixed_CI_high =
      tanh(
        fixed_high
      ),
    
    fixed_p =
      fixed_p,
    
    random_rho =
      tanh(
        pooled_z_random
      ),
    
    random_CI_low =
      tanh(
        random_low
      ),
    
    random_CI_high =
      tanh(
        random_high
      ),
    
    random_p =
      random_p
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
  
  
  if (
    isTRUE(
      capabilities(
        "cairo"
      )
    )
  ) {
    
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
# 5. INPUT DIRECTORIES
# ==============================================================================

dir_135 <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "135_validation",
  "tables"
)


dir_141 <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "141_external_validation_GSE154918",
  "tables"
)


dir_142 <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "142b_external_validation_GSE185263",
  "tables"
)


for (
  dir_required in c(
    dir_135,
    dir_141,
    dir_142
  )
) {
  
  if (!dir.exists(
    dir_required
  )) {
    
    stop(
      paste0(
        "Required upstream directory not found:\n",
        dir_required
      )
    )
  }
}


# ==============================================================================
# 6. LOCATE SCRIPT 135 FILE
# ==============================================================================

file_135_preferred <- file.path(
  dir_135,
  "135_sepsis_blood_scores_with_SRS_CTS.csv"
)


if (
  file.exists(
    file_135_preferred
  )
) {
  
  file_135 <- file_135_preferred
  
} else {
  
  candidate_135_files <- list.files(
    dir_135,
    pattern = "\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  
  if (
    length(
      candidate_135_files
    ) == 0
  ) {
    
    stop(
      "No Script 135 CSV files found."
    )
  }
  
  
  file_135 <- candidate_135_files[1]
  
  
  cat(
    "Preferred Script 135 score file not found.\n"
  )
  
  cat(
    "Using candidate:\n",
    file_135,
    "\n\n",
    sep = ""
  )
}


# ==============================================================================
# 7. LOCATE SCRIPT 141 FILES ROBUSTLY
# ==============================================================================

file_141_scores <- find_csv_by_required_columns(
  
  directory =
    dir_141,
  
  required_columns = c(
    "status",
    "five_gene_score"
  ),
  
  label =
    "Script 141 sample-score table"
)


file_141_comparisons <- find_csv_by_required_columns(
  
  directory =
    dir_141,
  
  required_columns = c(
    "comparison",
    "auc_fixed_direction"
  ),
  
  label =
    "Script 141 comparison table"
)


file_141_direction <- find_csv_by_required_columns(
  
  directory =
    dir_141,
  
  required_columns = c(
    "gene",
    "observed_direction",
    "direction_concordant"
  ),
  
  label =
    "Script 141 gene-direction table"
)


# ==============================================================================
# 8. SCRIPT 142b FILES
# ==============================================================================

file_142_scores <- file.path(
  dir_142,
  "142b_GSE185263_sample_scores.csv"
)


file_142_primary <- file.path(
  dir_142,
  "142b_PRIMARY_score_vs_SOFA.csv"
)


file_142_secondary <- file.path(
  dir_142,
  "142b_secondary_score_associations.csv"
)


file_142_gene_sofa <- file.path(
  dir_142,
  "142b_component_gene_SOFA_associations.csv"
)


file_142_adjusted <- file.path(
  dir_142,
  "142b_age_sex_location_adjusted_model.csv"
)


file_142_location <- file.path(
  dir_142,
  "142b_location_specific_SOFA_correlations.csv"
)


required_files <- c(
  file_135,
  file_141_scores,
  file_141_comparisons,
  file_141_direction,
  file_142_scores,
  file_142_primary,
  file_142_secondary,
  file_142_gene_sofa,
  file_142_adjusted,
  file_142_location
)


missing_files <- required_files[
  !file.exists(
    required_files
  )
]


if (
  length(
    missing_files
  ) > 0
) {
  
  stop(
    paste0(
      "Required upstream result file(s) missing:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


cat(
  "All required upstream result files found.\n\n"
)


cat(
  "Script 135 input:\n",
  file_135,
  "\n\n",
  sep = ""
)


cat(
  "Script 141 score input:\n",
  file_141_scores,
  "\n\n",
  sep = ""
)


cat(
  "Script 141 comparison input:\n",
  file_141_comparisons,
  "\n\n",
  sep = ""
)


cat(
  "Script 141 direction input:\n",
  file_141_direction,
  "\n\n",
  sep = ""
)


# ==============================================================================
# 9. READ UPSTREAM TABLES
# ==============================================================================

data_135 <- read_csv_fast(
  file_135
)


data_141_scores <- read_csv_fast(
  file_141_scores
)


data_141_comparisons <- read_csv_fast(
  file_141_comparisons
)


data_141_direction <- read_csv_fast(
  file_141_direction
)


data_142_scores <- read_csv_fast(
  file_142_scores
)


data_142_primary <- read_csv_fast(
  file_142_primary
)


data_142_secondary <- read_csv_fast(
  file_142_secondary
)


data_142_gene_sofa <- read_csv_fast(
  file_142_gene_sofa
)


data_142_adjusted <- read_csv_fast(
  file_142_adjusted
)


data_142_location <- read_csv_fast(
  file_142_location
)


cat(
  "Upstream result tables loaded successfully.\n\n"
)


# ==============================================================================
# 10. DISCOVERY DATA: DETECT IMPORTANT COLUMNS
# ==============================================================================

score_135_col <- find_column(
  
  data_135,
  
  patterns = c(
    "^primary_score$",
    "^Primary_5_gene_score$",
    "primary.*five.*gene.*score",
    "primary.*5.*gene.*score",
    "five.*gene.*score",
    "primary.*score",
    "host.*response.*score"
  ),
  
  exclude_pattern =
    "DCAF|Septi|SRS|CTS",
  
  label =
    "Script 135 primary five-gene score"
)


srs_135_col <- find_column(
  
  data_135,
  
  patterns = c(
    "^SRS$",
    "^SRS_class$",
    "^SRS_final$",
    "FINAL.*SRS",
    "SRS.*class",
    "SRS.*label"
  ),
  
  exclude_pattern =
    "SRSq|score|prob|distance|corr",
  
  label =
    "Script 135 SRS class"
)


cts_135_col <- find_column(
  
  data_135,
  
  patterns = c(
    "^CTS$",
    "^CTS_class$",
    "^CTS_final$",
    "FINAL.*CTS",
    "CTS.*class",
    "CTS.*label"
  ),
  
  exclude_pattern =
    "score|prob|distance|corr",
  
  label =
    "Script 135 CTS class"
)


srsq_candidates <- names(
  data_135
)[
  grepl(
    "SRSq",
    names(
      data_135
    ),
    ignore.case = TRUE
  )
]


if (
  length(
    srsq_candidates
  ) > 0
) {
  
  srsq_135_col <- srsq_candidates[1]
  
} else {
  
  srsq_135_col <- NA_character_
}


cat(
  "Detected Script 135 columns:\n"
)


cat(
  "Primary score: ",
  score_135_col,
  "\n",
  sep = ""
)


cat(
  "SRS class: ",
  srs_135_col,
  "\n",
  sep = ""
)


cat(
  "CTS class: ",
  cts_135_col,
  "\n",
  sep = ""
)


cat(
  "SRSq: ",
  ifelse(
    is.na(
      srsq_135_col
    ),
    "<not detected - locked upstream result will be used>",
    srsq_135_col
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 11. DISCOVERY ENDOTYPE HIERARCHY
# ==============================================================================

discovery_endotype <- data_135 %>%
  
  dplyr::transmute(
    
    five_gene_score =
      suppressWarnings(
        as.numeric(
          as.character(
            .data[[score_135_col]]
          )
        )
      ),
    
    SRS =
      as.character(
        .data[[srs_135_col]]
      ),
    
    CTS =
      as.character(
        .data[[cts_135_col]]
      )
  ) %>%
  
  dplyr::filter(
    
    is.finite(
      five_gene_score
    ),
    
    !is.na(
      SRS
    ),
    
    SRS != "",
    
    !is.na(
      CTS
    ),
    
    CTS != ""
  ) %>%
  
  dplyr::mutate(
    
    integrated_group =
      paste0(
        CTS,
        "/",
        SRS
      )
  )


preferred_integrated_order <- c(
  "CTS1/SRS1",
  "CTS2/SRS1",
  "CTS3/SRS1",
  "CTS3/SRS2"
)


observed_integrated <- unique(
  discovery_endotype$integrated_group
)


final_integrated_order <- c(
  
  preferred_integrated_order[
    preferred_integrated_order %in%
      observed_integrated
  ],
  
  setdiff(
    observed_integrated,
    preferred_integrated_order
  )
)


discovery_endotype$integrated_group <- factor(
  
  discovery_endotype$integrated_group,
  
  levels =
    final_integrated_order
)


integrated_group_summary <- discovery_endotype %>%
  
  dplyr::group_by(
    integrated_group
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    median =
      stats::median(
        five_gene_score
      ),
    
    q1 =
      stats::quantile(
        five_gene_score,
        0.25
      ),
    
    q3 =
      stats::quantile(
        five_gene_score,
        0.75
      ),
    
    mean =
      mean(
        five_gene_score
      ),
    
    sd =
      stats::sd(
        five_gene_score
      ),
    
    .groups =
      "drop"
  )


cts_result_135 <- kruskal_epsilon2(
  
  value =
    data_135[[score_135_col]],
  
  group =
    data_135[[cts_135_col]]
)


if (
  !is.na(
    srsq_135_col
  )
) {
  
  srsq_result_135 <- safe_spearman(
    
    data_135[[score_135_col]],
    
    data_135[[srsq_135_col]]
  )
  
} else {
  
  # Locked upstream Script 135 result.
  # Used only if SRSq is not exported in the selected sample-level table.
  
  srsq_result_135 <- tibble::tibble(
    
    n =
      35L,
    
    rho =
      0.765,
    
    p_value =
      8.74e-08
  )
}


cat(
  "Discovery endotype integration:\n"
)


print(
  integrated_group_summary,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "Discovery SRSq association:\n"
)


print(
  srsq_result_135,
  n = Inf
)


cat("\n")


cat(
  "Discovery CTS association:\n"
)


print(
  cts_result_135,
  n = Inf
)


cat("\n")


# ==============================================================================
# 12. GSE154918 CLINICAL-STATE GRADIENT
# ==============================================================================

required_141_score_columns <- c(
  "status",
  "five_gene_score"
)


if (
  !all(
    required_141_score_columns %in%
    names(
      data_141_scores
    )
  )
) {
  
  stop(
    "Required columns missing from Script 141 score table."
  )
}


gse154918_scores <- data_141_scores %>%
  
  dplyr::filter(
    status %in%
      c(
        "Hlty",
        "Inf1_P",
        "Seps_P",
        "Shock_P"
      )
  ) %>%
  
  dplyr::mutate(
    
    status =
      factor(
        status,
        levels = c(
          "Hlty",
          "Inf1_P",
          "Seps_P",
          "Shock_P"
        ),
        ordered = TRUE
      )
  )


gse154918_group_summary <- gse154918_scores %>%
  
  dplyr::group_by(
    status
  ) %>%
  
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
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


row_141_primary <- data_141_comparisons %>%
  
  dplyr::filter(
    comparison ==
      "Sepsis_or_shock_vs_uncomplicated"
  )


row_141_shock <- data_141_comparisons %>%
  
  dplyr::filter(
    comparison ==
      "Shock_vs_uncomplicated"
  )


if (
  nrow(
    row_141_primary
  ) != 1
) {
  
  stop(
    "Could not uniquely identify GSE154918 primary comparison."
  )
}


if (
  nrow(
    row_141_shock
  ) != 1
) {
  
  stop(
    "Could not uniquely identify GSE154918 Shock vs uncomplicated comparison."
  )
}


# ==============================================================================
# 13. GSE185263 KEY RESULTS
# ==============================================================================

if (
  nrow(
    data_142_primary
  ) != 1
) {
  
  stop(
    "Script 142b primary endpoint table should contain exactly one row."
  )
}


primary_142 <- tibble::as_tibble(
  data_142_primary
)


adjusted_142 <- tibble::as_tibble(
  data_142_adjusted
) %>%
  
  dplyr::filter(
    term ==
      "sofa"
  )


if (
  nrow(
    adjusted_142
  ) != 1
) {
  
  stop(
    "Could not uniquely identify adjusted SOFA coefficient from Script 142b."
  )
}


# ==============================================================================
# 14. CROSS-COHORT FIVE-GENE DIRECTION MATRIX
# ==============================================================================

direction_141 <- data_141_direction %>%
  
  dplyr::filter(
    gene %in%
      five_genes
  ) %>%
  
  dplyr::transmute(
    
    gene =
      gene,
    
    GSE154918 =
      dplyr::case_when(
        
        toupper(
          observed_direction
        ) ==
          "UP" ~
          "POSITIVE",
        
        toupper(
          observed_direction
        ) ==
          "DOWN" ~
          "NEGATIVE",
        
        toupper(
          observed_direction
        ) ==
          "POSITIVE" ~
          "POSITIVE",
        
        toupper(
          observed_direction
        ) ==
          "NEGATIVE" ~
          "NEGATIVE",
        
        TRUE ~
          NA_character_
      )
  )


direction_142 <- data_142_gene_sofa %>%
  
  dplyr::filter(
    gene %in%
      five_genes
  ) %>%
  
  dplyr::transmute(
    
    gene =
      gene,
    
    GSE185263 =
      toupper(
        observed_direction
      )
  )


direction_matrix <- tibble::tibble(
  
  gene =
    five_genes,
  
  Discovery =
    unname(
      discovery_direction[
        five_genes
      ]
    )
) %>%
  
  dplyr::left_join(
    direction_141,
    by = "gene"
  ) %>%
  
  dplyr::left_join(
    direction_142,
    by = "gene"
  )


if (
  any(
    is.na(
      direction_matrix$GSE154918
    )
  )
) {
  
  warning(
    "One or more GSE154918 gene directions were not detected."
  )
}


if (
  any(
    is.na(
      direction_matrix$GSE185263
    )
  )
) {
  
  warning(
    "One or more GSE185263 gene directions were not detected."
  )
}


direction_concordance_summary <- direction_matrix %>%
  
  dplyr::mutate(
    
    GSE154918_concordant =
      Discovery ==
      GSE154918,
    
    GSE185263_concordant =
      Discovery ==
      GSE185263
  )


direction_long <- direction_matrix %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      Discovery,
      GSE154918,
      GSE185263
    ),
    
    names_to =
      "evidence_source",
    
    values_to =
      "direction"
  ) %>%
  
  dplyr::mutate(
    
    evidence_source =
      factor(
        evidence_source,
        levels = c(
          "Discovery",
          "GSE154918",
          "GSE185263"
        )
      ),
    
    gene =
      factor(
        gene,
        levels =
          rev(
            five_genes
          )
      ),
    
    direction =
      factor(
        direction,
        levels = c(
          "NEGATIVE",
          "POSITIVE"
        )
      ),
    
    symbol =
      ifelse(
        !is.na(
          direction
        ),
        "\u2713",
        ""
      )
  )


cat(
  "Cross-cohort gene direction matrix:\n"
)


print(
  direction_matrix,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 15. GSE185263 LOCATION-SPECIFIC RANDOM-EFFECTS META-ANALYSIS
# ==============================================================================

location_meta_input <- data_142_location %>%
  
  dplyr::filter(
    
    is.finite(
      rho
    ),
    
    is.finite(
      n
    ),
    
    n > 3,
    
    abs(
      rho
    ) < 1
  )


if (
  nrow(
    location_meta_input
  ) < 2
) {
  
  stop(
    "Too few valid GSE185263 location-specific correlations."
  )
}


location_ci_list <- lapply(
  
  seq_len(
    nrow(
      location_meta_input
    )
  ),
  
  function(i) {
    
    correlation_ci_fisher(
      
      rho =
        location_meta_input$rho[i],
      
      n =
        location_meta_input$n[i]
    )
  }
)


location_meta_input$CI_low <- vapply(
  
  location_ci_list,
  
  function(x) {
    
    as.numeric(
      x["lower"]
    )
  },
  
  numeric(1)
)


location_meta_input$CI_high <- vapply(
  
  location_ci_list,
  
  function(x) {
    
    as.numeric(
      x["upper"]
    )
  },
  
  numeric(1)
)


location_meta_result <- meta_correlations_fisher(
  
  rho =
    location_meta_input$rho,
  
  n =
    location_meta_input$n
)


cat(
  "GSE185263 geographic random-effects synthesis:\n"
)


print(
  location_meta_result,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 16. FOREST-PLOT DATA
# ==============================================================================

forest_data <- location_meta_input %>%
  
  dplyr::transmute(
    
    label =
      stringr::str_to_title(
        collection_location
      ),
    
    n =
      n,
    
    rho =
      rho,
    
    CI_low =
      CI_low,
    
    CI_high =
      CI_high,
    
    type =
      "Location"
  )


forest_data <- dplyr::bind_rows(
  
  forest_data,
  
  tibble::tibble(
    
    label =
      "Random-effects pooled",
    
    n =
      location_meta_result$total_n,
    
    rho =
      location_meta_result$random_rho,
    
    CI_low =
      location_meta_result$random_CI_low,
    
    CI_high =
      location_meta_result$random_CI_high,
    
    type =
      "Pooled"
  )
)


forest_order <- c(
  as.character(
    forest_data$label[
      forest_data$type ==
        "Location"
    ]
  ),
  "Random-effects pooled"
)


forest_data$label <- factor(
  
  forest_data$label,
  
  levels =
    rev(
      forest_order
    )
)


# ==============================================================================
# 17. LOCKED SCRIPT 137 BENCHMARKING RESULTS
#
# These values come from the completed Script 137.
# Script 143 performs no new benchmarking or model fitting.
# ==============================================================================

benchmark_results <- tibble::tibble(
  
  signature = c(
    "Primary five-gene",
    "DCAF17 alternative",
    "LIFTS",
    "FAIM3:PLAC8-like",
    "Sepsis MetaScore-like",
    "RAPID-like",
    "SeptiCyte LAB-like"
  ),
  
  category = c(
    "Current study",
    "Current study",
    "Published comparator",
    "Published comparator",
    "Published comparator",
    "Published comparator",
    "Published comparator"
  ),
  
  SRSq_rho = c(
    0.765,
    0.794,
    0.852,
    0.716,
    0.685,
    0.495,
    0.462
  ),
  
  CTS_epsilon2 = c(
    0.661,
    0.610,
    0.698,
    0.649,
    0.529,
    0.436,
    0.328
  ),
  
  BP_vs_BC_AUC = c(
    1.000,
    1.000,
    1.000,
    0.946,
    0.974,
    0.900,
    0.720
  )
)


# ==============================================================================
# 18. LOCKED DISCOVERY / DEMOGRAPHIC RESULTS
# ==============================================================================

discovery_locked <- tibble::tibble(
  
  result = c(
    "Robust core blood DEGs",
    "Primary five-gene apparent BP-vs-BC AUC",
    "Primary score vs SRSq",
    "Primary score across CTS",
    "Age/sex-adjusted BP-vs-BC score effect"
  ),
  
  estimate = c(
    1796,
    1.000,
    srsq_result_135$rho,
    cts_result_135$epsilon2,
    3.88
  ),
  
  estimate_type = c(
    "n genes",
    "AUC",
    "Spearman rho",
    "Kruskal-Wallis epsilon2",
    "Adjusted score difference"
  ),
  
  p_value = c(
    NA_real_,
    NA_real_,
    srsq_result_135$p_value,
    cts_result_135$p_value,
    3.73e-11
  ),
  
  source = c(
    "Upstream robust-core analysis",
    "Script 135",
    "Script 135",
    "Script 135",
    "Script 136b"
  )
)


# ==============================================================================
# 19. GSE185263 SECONDARY RESULTS
# ==============================================================================

secondary_142_sofa <- data_142_secondary %>%
  
  dplyr::filter(
    analysis ==
      "SOFA_ge2_vs_SOFA_0_1"
  )


secondary_142_mortality <- data_142_secondary %>%
  
  dplyr::filter(
    analysis ==
      "Died_vs_Survived"
  )


secondary_142_site <- data_142_secondary %>%
  
  dplyr::filter(
    analysis ==
      "ICU_vs_Emergency_Room"
  )


secondary_142_disease <- data_142_secondary %>%
  
  dplyr::filter(
    analysis ==
      "Sepsis_vs_healthy_contextual"
  )


for (
  object_name in c(
    "secondary_142_sofa",
    "secondary_142_mortality",
    "secondary_142_site",
    "secondary_142_disease"
  )
) {
  
  object_value <- get(
    object_name
  )
  
  
  if (
    nrow(
      object_value
    ) != 1
  ) {
    
    stop(
      paste0(
        "Could not uniquely identify Script 142b result: ",
        object_name
      )
    )
  }
}


# ==============================================================================
# 20. MULTICOHORT EVIDENCE SUMMARY
# ==============================================================================

evidence_summary <- tibble::tibble(
  
  evidence_level = c(
    
    "Discovery",
    
    "Discovery",
    
    "Discovery",
    
    "External GSE154918",
    
    "External GSE154918",
    
    "External GSE154918",
    
    "External GSE185263",
    
    "External GSE185263",
    
    "External GSE185263",
    
    "External GSE185263",
    
    "External GSE185263",
    
    "External GSE185263",
    
    "External GSE185263 geographic"
  ),
  
  analysis = c(
    
    "Primary score vs SRSq",
    
    "Primary score across CTS",
    
    "Age/sex-adjusted sepsis-vs-healthy score",
    
    "Frozen genes with expected direction",
    
    "Sepsis/shock vs uncomplicated infection",
    
    "Shock vs uncomplicated infection",
    
    "Primary score vs continuous SOFA",
    
    "Frozen genes with expected SOFA direction",
    
    "Age/sex/location-adjusted SOFA coefficient",
    
    "SOFA >=2 vs SOFA 0-1",
    
    "Died vs Survived",
    
    "ICU vs Emergency Room",
    
    "Random-effects geographic score-SOFA correlation"
  ),
  
  effect_type = c(
    
    "Spearman rho",
    
    "epsilon2",
    
    "adjusted score difference",
    
    "directional concordance",
    
    "AUC",
    
    "AUC",
    
    "Spearman rho",
    
    "directional concordance",
    
    "beta per SOFA point",
    
    "AUC",
    
    "AUC",
    
    "AUC",
    
    "random-effects pooled rho"
  ),
  
  estimate = c(
    
    srsq_result_135$rho,
    
    cts_result_135$epsilon2,
    
    3.88,
    
    sum(
      direction_concordance_summary$GSE154918_concordant,
      na.rm = TRUE
    ),
    
    row_141_primary$auc_fixed_direction,
    
    row_141_shock$auc_fixed_direction,
    
    primary_142$rho,
    
    sum(
      direction_concordance_summary$GSE185263_concordant,
      na.rm = TRUE
    ),
    
    adjusted_142$estimate,
    
    secondary_142_sofa$AUC,
    
    secondary_142_mortality$AUC,
    
    secondary_142_site$AUC,
    
    location_meta_result$random_rho
  ),
  
  p_value = c(
    
    srsq_result_135$p_value,
    
    cts_result_135$p_value,
    
    3.73e-11,
    
    NA_real_,
    
    row_141_primary$p_value,
    
    row_141_shock$p_value,
    
    primary_142$p_value,
    
    NA_real_,
    
    adjusted_142$p_value,
    
    secondary_142_sofa$p_value,
    
    secondary_142_mortality$p_value,
    
    secondary_142_site$p_value,
    
    location_meta_result$random_p
  ),
  
  interpretation = c(
    
    "Strong alignment with the SRS immune-dysfunction axis.",
    
    "Strong separation across consensus transcriptomic subtypes.",
    
    "Discovery sepsis association persists after demographic adjustment.",
    
    "All five frozen genes reproduce the prespecified direction.",
    
    "Prespecified infection-control comparison shows modest discrimination and is not statistically significant.",
    
    "Score is higher in septic shock than uncomplicated infection.",
    
    "Prespecified independent severity endpoint is positive.",
    
    "All five components reproduce prespecified SOFA directions and are FDR-significant.",
    
    "Association with SOFA persists after age, sex and geographic adjustment.",
    
    "Higher score is associated with greater organ dysfunction.",
    
    "Higher score is observed among nonsurvivors; secondary association only.",
    
    "Higher score is observed in ICU samples; secondary contextual association.",
    
    "Positive score-SOFA association is reproducible across geographic subcohorts."
  )
)


# ==============================================================================
# 21. FIGURE A
# CROSS-COHORT DIRECTIONAL REPLICATION
# ==============================================================================

direction_colors <- c(
  "NEGATIVE" = "#0072B2",
  "POSITIVE" = "#D55E00"
)


pA <- ggplot2::ggplot(
  
  direction_long,
  
  ggplot2::aes(
    x =
      evidence_source,
    y =
      gene,
    fill =
      direction
  )
) +
  
  ggplot2::geom_tile(
    color = "white",
    linewidth = 1.2
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        symbol
    ),
    
    size = 6,
    color = "white",
    fontface = "bold"
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      direction_colors,
    drop = FALSE,
    na.value = "grey80"
  ) +
  
  ggplot2::labs(
    
    title =
      "Directional replication of the frozen five-gene panel",
    
    subtitle =
      paste0(
        "Discovery direction, GSE154918 disease-state contrast, ",
        "and GSE185263 SOFA association"
      ),
    
    x =
      NULL,
    
    y =
      NULL,
    
    fill =
      "Direction"
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid =
      ggplot2::element_blank(),
    
    axis.text.x =
      ggplot2::element_text(
        face = "bold"
      ),
    
    axis.text.y =
      ggplot2::element_text(
        face = "italic"
      )
  )


# ==============================================================================
# 22. FIGURE B
# DISCOVERY CTS / SRS HIERARCHY
# ==============================================================================

integrated_colors <- c(
  
  "CTS1/SRS1" =
    "#D55E00",
  
  "CTS2/SRS1" =
    "#E69F00",
  
  "CTS3/SRS1" =
    "#56B4E9",
  
  "CTS3/SRS2" =
    "#0072B2"
)


missing_integrated_colors <- setdiff(
  levels(
    discovery_endotype$integrated_group
  ),
  names(
    integrated_colors
  )
)


if (
  length(
    missing_integrated_colors
  ) > 0
) {
  
  extra_colors <- rep(
    "#7F7F7F",
    length(
      missing_integrated_colors
    )
  )
  
  
  names(
    extra_colors
  ) <- missing_integrated_colors
  
  
  integrated_colors <- c(
    integrated_colors,
    extra_colors
  )
}


pB <- ggplot2::ggplot(
  
  discovery_endotype,
  
  ggplot2::aes(
    x =
      integrated_group,
    y =
      five_gene_score,
    fill =
      integrated_group,
    color =
      integrated_group
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.62,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  
  ggplot2::geom_jitter(
    width = 0.10,
    height = 0,
    size = 2,
    alpha = 0.75
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      integrated_colors
  ) +
  
  ggplot2::scale_color_manual(
    values =
      integrated_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Discovery endotype hierarchy",
    
    subtitle =
      paste0(
        "CTS epsilon² = ",
        sprintf(
          "%.3f",
          cts_result_135$epsilon2
        ),
        "; SRSq rho = ",
        sprintf(
          "%.3f",
          srsq_result_135$rho
        )
      ),
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 25,
        hjust = 1
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 23. FIGURE C
# GSE154918 CLINICAL-STATE GRADIENT
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


pC <- ggplot2::ggplot(
  
  gse154918_scores,
  
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
    width = 0.11,
    size = 1.8,
    alpha = 0.70
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
      "Independent clinical-state replication",
    
    subtitle =
      paste0(
        "GSE154918; primary AUC = ",
        sprintf(
          "%.3f",
          row_141_primary$auc_fixed_direction
        ),
        "; shock vs uncomplicated AUC = ",
        sprintf(
          "%.3f",
          row_141_shock$auc_fixed_direction
        )
      ),
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
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
# 24. FIGURE D
# GSE185263 SCORE vs SOFA
# ==============================================================================

gse185263_sofa <- data_142_scores %>%
  
  dplyr::filter(
    
    disease_state ==
      "sepsis",
    
    is.finite(
      sofa
    ),
    
    is.finite(
      five_gene_score
    )
  )


pD <- ggplot2::ggplot(
  
  gse185263_sofa,
  
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
    size = 1.8,
    alpha = 0.45,
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
      "Independent association with organ dysfunction",
    
    subtitle =
      paste0(
        "GSE185263; rho = ",
        sprintf(
          "%.3f",
          primary_142$rho
        ),
        "; p = ",
        format.pval(
          primary_142$p_value,
          digits = 3
        ),
        "; n = ",
        primary_142$n
      ),
    
    x =
      "SOFA score, 24 h after admission",
    
    y =
      "Five-gene host-response score"
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
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
# 25. FIGURE E
# GSE185263 GEOGRAPHIC FOREST PLOT
# ==============================================================================

pE <- ggplot2::ggplot(
  
  forest_data,
  
  ggplot2::aes(
    x =
      rho,
    y =
      label
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x =
        CI_low,
      xend =
        CI_high,
      y =
        label,
      yend =
        label,
      color =
        type
    ),
    
    linewidth = 1
  ) +
  
  ggplot2::geom_point(
    
    ggplot2::aes(
      color =
        type,
      shape =
        type
    ),
    
    size = 3.7
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "Location" = "#D55E00",
      "Pooled" = "#000000"
    )
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "Location" = 16,
      "Pooled" = 18
    )
  ) +
  
  ggplot2::labs(
    
    title =
      "Geographic replication of the SOFA association",
    
    subtitle =
      paste0(
        "Random-effects rho = ",
        sprintf(
          "%.3f",
          location_meta_result$random_rho
        ),
        " [",
        sprintf(
          "%.3f",
          location_meta_result$random_CI_low
        ),
        ", ",
        sprintf(
          "%.3f",
          location_meta_result$random_CI_high
        ),
        "]; I² = ",
        sprintf(
          "%.1f",
          location_meta_result$I2_percent
        ),
        "%"
      ),
    
    x =
      "Spearman rho with SOFA",
    
    y =
      NULL,
    
    color =
      NULL,
    
    shape =
      NULL
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
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
# 26. FIGURE F
# CONVERGENCE WITH PUBLISHED TRANSCRIPTOMIC SIGNATURES
# ==============================================================================

benchmark_colors <- c(
  
  "Current study" =
    "#D55E00",
  
  "Published comparator" =
    "#0072B2"
)


pF <- ggplot2::ggplot(
  
  benchmark_results,
  
  ggplot2::aes(
    x =
      SRSq_rho,
    y =
      CTS_epsilon2,
    color =
      category
  )
) +
  
  ggplot2::geom_point(
    size = 4,
    alpha = 0.90
  ) +
  
  ggrepel::geom_text_repel(
    
    ggplot2::aes(
      label =
        signature
    ),
    
    size = 3.3,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  
  ggplot2::scale_color_manual(
    values =
      benchmark_colors
  ) +
  
  ggplot2::labs(
    
    title =
      "Published signatures converge on the same host-response axis",
    
    subtitle =
      paste0(
        "Relationships with SRSq and consensus transcriptomic ",
        "subtypes in the discovery cohort"
      ),
    
    x =
      "Spearman rho with SRSq",
    
    y =
      "CTS effect size (epsilon²)",
    
    color =
      NULL
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "bottom",
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


# ==============================================================================
# 27. SAVE INDIVIDUAL FIGURES
# ==============================================================================

save_plot_all_formats(
  pA,
  "143_Figure_A_cross_cohort_gene_direction",
  width = 8.5,
  height = 6
)


save_plot_all_formats(
  pB,
  "143_Figure_B_discovery_endotype_hierarchy",
  width = 8.5,
  height = 6
)


save_plot_all_formats(
  pC,
  "143_Figure_C_GSE154918_state_gradient",
  width = 8.5,
  height = 6
)


save_plot_all_formats(
  pD,
  "143_Figure_D_GSE185263_score_SOFA",
  width = 8.5,
  height = 6
)


save_plot_all_formats(
  pE,
  "143_Figure_E_GSE185263_geographic_forest",
  width = 9,
  height = 6
)


save_plot_all_formats(
  pF,
  "143_Figure_F_published_signature_convergence",
  width = 9,
  height = 7
)


# ==============================================================================
# 28. MAIN MULTICOHORT FIGURE
# ==============================================================================

main_figure <- (
  
  pA +
    pB
  
) / (
  
  pC +
    pD
  
) / (
  
  pE +
    pF
  
) +
  
  patchwork::plot_annotation(
    
    title =
      paste0(
        "A frozen five-gene score captures a reproducible ",
        "sepsis host-response axis across cohorts"
      ),
    
    subtitle =
      paste0(
        "Discovery endotypes, independent clinical-state replication, ",
        "organ-dysfunction severity, and geographic reproducibility"
      ),
    
    tag_levels =
      "A"
  )


save_plot_all_formats(
  main_figure,
  "143_MAIN_Figure_multicohort_five_gene_evidence",
  width = 18,
  height = 21
)


cat(
  "Integrated publication figures saved.\n\n"
)


# ==============================================================================
# 29. SUPPLEMENTARY EXTERNAL ENDPOINT FIGURE
# ==============================================================================

external_endpoint_plot_data <- tibble::tibble(
  
  analysis = c(
    
    "GSE154918:\nSepsis/shock vs uncomplicated",
    
    "GSE154918:\nShock vs uncomplicated",
    
    "GSE185263:\nSOFA >=2 vs 0-1",
    
    "GSE185263:\nDied vs Survived",
    
    "GSE185263:\nICU vs Emergency Room",
    
    "GSE185263:\nSepsis vs healthy"
  ),
  
  AUC = c(
    
    row_141_primary$auc_fixed_direction,
    
    row_141_shock$auc_fixed_direction,
    
    secondary_142_sofa$AUC,
    
    secondary_142_mortality$AUC,
    
    secondary_142_site$AUC,
    
    secondary_142_disease$AUC
  ),
  
  CI_low = c(
    
    row_141_primary$auc_ci_low,
    
    row_141_shock$auc_ci_low,
    
    secondary_142_sofa$CI_low,
    
    secondary_142_mortality$CI_low,
    
    secondary_142_site$CI_low,
    
    secondary_142_disease$CI_low
  ),
  
  CI_high = c(
    
    row_141_primary$auc_ci_high,
    
    row_141_shock$auc_ci_high,
    
    secondary_142_sofa$CI_high,
    
    secondary_142_mortality$CI_high,
    
    secondary_142_site$CI_high,
    
    secondary_142_disease$CI_high
  ),
  
  role = c(
    
    "Prespecified primary",
    
    "Secondary",
    
    "Secondary",
    
    "Secondary",
    
    "Secondary",
    
    "Contextual"
  )
)


external_endpoint_plot_data$analysis <- factor(
  
  external_endpoint_plot_data$analysis,
  
  levels =
    rev(
      external_endpoint_plot_data$analysis
    )
)


endpoint_colors <- c(
  
  "Prespecified primary" =
    "#D55E00",
  
  "Secondary" =
    "#0072B2",
  
  "Contextual" =
    "#7F8C8D"
)


p_external_endpoints <- ggplot2::ggplot(
  
  external_endpoint_plot_data,
  
  ggplot2::aes(
    x =
      AUC,
    y =
      analysis,
    color =
      role
  )
) +
  
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dashed"
  ) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x =
        CI_low,
      xend =
        CI_high,
      y =
        analysis,
      yend =
        analysis
    ),
    
    linewidth = 1
  ) +
  
  ggplot2::geom_point(
    size = 4
  ) +
  
  ggplot2::scale_color_manual(
    values =
      endpoint_colors
  ) +
  
  ggplot2::labs(
    
    title =
      paste0(
        "External binary comparisons are secondary to the ",
        "host-response-state hypothesis"
      ),
    
    x =
      "AUC",
    
    y =
      NULL,
    
    color =
      NULL
  ) +
  
  ggplot2::theme_bw(
    base_size = 13
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    legend.position =
      "bottom",
    
    panel.grid.minor =
      ggplot2::element_blank()
  )


save_plot_all_formats(
  p_external_endpoints,
  "143_Supplementary_external_binary_endpoint_summary",
  width = 10,
  height = 7
)


# ==============================================================================
# 30. MANUSCRIPT KEY NUMBERS
# ==============================================================================

key_numbers <- tibble::tibble(
  
  domain = c(
    
    "Discovery",
    
    "Discovery",
    
    "Discovery",
    
    "Discovery",
    
    "GSE154918",
    
    "GSE154918",
    
    "GSE154918",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263 geographic",
    
    "GSE185263 geographic"
  ),
  
  metric = c(
    
    "Robust core DEGs",
    
    "Primary score vs SRSq",
    
    "CTS effect size",
    
    "Age/sex-adjusted BP-vs-BC score difference",
    
    "Expected gene directions",
    
    "Primary sepsis/shock vs uncomplicated AUC",
    
    "Shock vs uncomplicated AUC",
    
    "Primary score vs SOFA",
    
    "Expected gene-SOFA directions",
    
    "Adjusted SOFA coefficient",
    
    "SOFA >=2 vs 0-1 AUC",
    
    "Mortality AUC",
    
    "Random-effects pooled rho",
    
    "I2"
  ),
  
  value = c(
    
    "1796",
    
    sprintf(
      "%.3f",
      srsq_result_135$rho
    ),
    
    sprintf(
      "%.3f",
      cts_result_135$epsilon2
    ),
    
    "3.88",
    
    paste0(
      sum(
        direction_concordance_summary$GSE154918_concordant,
        na.rm = TRUE
      ),
      "/5"
    ),
    
    sprintf(
      "%.3f",
      row_141_primary$auc_fixed_direction
    ),
    
    sprintf(
      "%.3f",
      row_141_shock$auc_fixed_direction
    ),
    
    sprintf(
      "%.3f",
      primary_142$rho
    ),
    
    paste0(
      sum(
        direction_concordance_summary$GSE185263_concordant,
        na.rm = TRUE
      ),
      "/5"
    ),
    
    sprintf(
      "%.3f",
      adjusted_142$estimate
    ),
    
    sprintf(
      "%.3f",
      secondary_142_sofa$AUC
    ),
    
    sprintf(
      "%.3f",
      secondary_142_mortality$AUC
    ),
    
    sprintf(
      "%.3f",
      location_meta_result$random_rho
    ),
    
    paste0(
      sprintf(
        "%.1f",
        location_meta_result$I2_percent
      ),
      "%"
    )
  ),
  
  p_value = c(
    
    NA_real_,
    
    srsq_result_135$p_value,
    
    cts_result_135$p_value,
    
    3.73e-11,
    
    NA_real_,
    
    row_141_primary$p_value,
    
    row_141_shock$p_value,
    
    primary_142$p_value,
    
    NA_real_,
    
    adjusted_142$p_value,
    
    secondary_142_sofa$p_value,
    
    secondary_142_mortality$p_value,
    
    location_meta_result$random_p,
    
    location_meta_result$Q_p
  )
)


# ==============================================================================
# 31. EXPORT CSV TABLES
# ==============================================================================

write.csv(
  
  evidence_summary,
  
  file.path(
    tables_dir,
    "143_multicohort_evidence_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  direction_matrix,
  
  file.path(
    tables_dir,
    "143_cross_cohort_gene_direction_matrix.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  direction_concordance_summary,
  
  file.path(
    tables_dir,
    "143_gene_direction_concordance_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  integrated_group_summary,
  
  file.path(
    tables_dir,
    "143_discovery_integrated_endotype_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  benchmark_results,
  
  file.path(
    tables_dir,
    "143_published_signature_benchmark_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  gse154918_group_summary,
  
  file.path(
    tables_dir,
    "143_GSE154918_group_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  location_meta_input,
  
  file.path(
    tables_dir,
    "143_GSE185263_location_correlations_with_CI.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  location_meta_result,
  
  file.path(
    tables_dir,
    "143_GSE185263_location_random_effects_meta.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  external_endpoint_plot_data,
  
  file.path(
    tables_dir,
    "143_external_binary_endpoint_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  discovery_locked,
  
  file.path(
    tables_dir,
    "143_discovery_locked_key_results.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  key_numbers,
  
  file.path(
    tables_dir,
    "143_manuscript_key_numbers.csv"
  ),
  
  row.names = FALSE
)


# ==============================================================================
# 32. MASTER EXCEL WORKBOOK
# ==============================================================================

run_info <- tibble::tibble(
  
  parameter = c(
    
    "script",
    
    "run_date",
    
    "purpose",
    
    "frozen_panel",
    
    "new_feature_selection",
    
    "new_gene_weighting",
    
    "cross_dataset_endpoint_meta_analysis",
    
    "new_statistical_analysis"
  ),
  
  value = c(
    
    script_name,
    
    as.character(
      run_date
    ),
    
    "Multicohort publication integration",
    
    paste(
      five_genes,
      collapse = "; "
    ),
    
    "NO",
    
    "NO",
    
    paste0(
      "NO - GSE154918 and GSE185263 have different ",
      "prespecified endpoints"
    ),
    
    paste0(
      "Random-effects Fisher-z meta-analysis of geographic ",
      "score-SOFA correlations within GSE185263"
    )
  )
)


workbook_file <- file.path(
  tables_dir,
  "143_multicohort_five_gene_evidence_package.xlsx"
)


wb <- openxlsx::createWorkbook()


sheet_list <- list(
  
  "00_run_info" =
    run_info,
  
  "01_key_numbers" =
    key_numbers,
  
  "02_evidence_summary" =
    evidence_summary,
  
  "03_gene_direction" =
    direction_matrix,
  
  "04_direction_concordance" =
    direction_concordance_summary,
  
  "05_discovery_endotypes" =
    integrated_group_summary,
  
  "06_benchmark_signatures" =
    benchmark_results,
  
  "07_GSE154918_groups" =
    gse154918_group_summary,
  
  "08_GSE154918_comparisons" =
    data_141_comparisons,
  
  "09_GSE185263_primary" =
    data_142_primary,
  
  "10_GSE185263_secondary" =
    data_142_secondary,
  
  "11_GSE185263_gene_SOFA" =
    data_142_gene_sofa,
  
  "12_GSE185263_adjusted" =
    data_142_adjusted,
  
  "13_location_correlations" =
    location_meta_input,
  
  "14_location_meta" =
    location_meta_result,
  
  "15_external_AUCs" =
    external_endpoint_plot_data,
  
  "16_discovery_locked" =
    discovery_locked
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
# 33. MANUSCRIPT-READY ENGLISH SUMMARY
# ==============================================================================

summary_en <- c(
  
  "MULTICOHORT FIVE-GENE HOST-RESPONSE EVIDENCE SUMMARY",
  
  "====================================================================",
  
  "",
  
  "DISCOVERY COHORT",
  
  paste0(
    "The five-gene score tracked the SRS immune-dysfunction axis ",
    "(Spearman rho = ",
    sprintf(
      "%.3f",
      srsq_result_135$rho
    ),
    ", p = ",
    format.pval(
      srsq_result_135$p_value,
      digits = 4
    ),
    ") and differed strongly across consensus transcriptomic ",
    "subtypes (epsilon2 = ",
    sprintf(
      "%.3f",
      cts_result_135$epsilon2
    ),
    ", p = ",
    format.pval(
      cts_result_135$p_value,
      digits = 4
    ),
    ")."
  ),
  
  "",
  
  "GSE154918 EXTERNAL REPLICATION",
  
  paste0(
    "All five frozen component genes reproduced the prespecified ",
    "direction in the independent GSE154918 whole-blood RNA-seq cohort."
  ),
  
  paste0(
    "The prespecified sepsis/septic-shock versus uncomplicated-infection ",
    "comparison showed modest discrimination (AUC = ",
    sprintf(
      "%.3f",
      row_141_primary$auc_fixed_direction
    ),
    ", 95% CI ",
    sprintf(
      "%.3f",
      row_141_primary$auc_ci_low
    ),
    "-",
    sprintf(
      "%.3f",
      row_141_primary$auc_ci_high
    ),
    "; p = ",
    format.pval(
      row_141_primary$p_value,
      digits = 4
    ),
    ")."
  ),
  
  paste0(
    "Patients with septic shock showed higher scores than patients ",
    "with uncomplicated infection (AUC = ",
    sprintf(
      "%.3f",
      row_141_shock$auc_fixed_direction
    ),
    "; p = ",
    format.pval(
      row_141_shock$p_value,
      digits = 4
    ),
    ")."
  ),
  
  "",
  
  "GSE185263 EXTERNAL SEVERITY VALIDATION",
  
  paste0(
    "In GSE185263, the prespecified primary endpoint was positive. ",
    "Among ",
    primary_142$n,
    " sepsis samples with available 24-h SOFA data, the frozen ",
    "five-gene score correlated with SOFA (rho = ",
    sprintf(
      "%.3f",
      primary_142$rho
    ),
    ", p = ",
    format.pval(
      primary_142$p_value,
      digits = 4
    ),
    ")."
  ),
  
  paste0(
    "All five component genes reproduced their prespecified direction ",
    "of association with SOFA and remained significant after ",
    "gene-level FDR correction."
  ),
  
  paste0(
    "The score-SOFA association persisted after adjustment for age, ",
    "sex and collection location (beta = ",
    sprintf(
      "%.3f",
      adjusted_142$estimate
    ),
    " score units per SOFA point; p = ",
    format.pval(
      adjusted_142$p_value,
      digits = 4
    ),
    ")."
  ),
  
  "",
  
  "GEOGRAPHIC REPRODUCIBILITY",
  
  paste0(
    "The score-SOFA association was directionally positive in all five ",
    "geographic subcohorts. Random-effects Fisher-z synthesis yielded ",
    "rho = ",
    sprintf(
      "%.3f",
      location_meta_result$random_rho
    ),
    " (95% CI ",
    sprintf(
      "%.3f",
      location_meta_result$random_CI_low
    ),
    "-",
    sprintf(
      "%.3f",
      location_meta_result$random_CI_high
    ),
    "; p = ",
    format.pval(
      location_meta_result$random_p,
      digits = 4
    ),
    "), with I2 = ",
    sprintf(
      "%.1f",
      location_meta_result$I2_percent
    ),
    "%."
  ),
  
  "",
  
  "INTEGRATED INTERPRETATION",
  
  paste0(
    "Together, the discovery and external-cohort analyses support ",
    "the five-gene signature as a compact and reproducible molecular ",
    "readout of a myeloid-adaptive host-response axis associated with ",
    "established transcriptomic endotypes and organ-dysfunction severity."
  ),
  
  paste0(
    "These results do not establish a calibrated clinical diagnostic ",
    "or prognostic assay."
  )
)


summary_en_file <- file.path(
  text_dir,
  "143_multicohort_evidence_summary_EN.txt"
)


writeLines(
  summary_en,
  con = summary_en_file,
  useBytes = TRUE
)


# ==============================================================================
# 34. MANUSCRIPT-READY RUSSIAN SUMMARY
# ==============================================================================

summary_ru <- c(
  
  "МУЛЬТИКОГОРТНОЕ РЕЗЮМЕ ДОКАЗАТЕЛЬСТВ ДЛЯ ПЯТИГЕННОЙ СИГНАТУРЫ",
  
  "====================================================================",
  
  "",
  
  "ИСХОДНАЯ КОГОРТА",
  
  paste0(
    "Пятигенный score был связан с осью SRS ",
    "(Spearman rho = ",
    sprintf(
      "%.3f",
      srsq_result_135$rho
    ),
    ", p = ",
    format.pval(
      srsq_result_135$p_value,
      digits = 4
    ),
    ") и выраженно различался между CTS-подтипами ",
    "(epsilon2 = ",
    sprintf(
      "%.3f",
      cts_result_135$epsilon2
    ),
    ")."
  ),
  
  "",
  
  "GSE154918",
  
  paste0(
    "Все 5 из 5 компонентов замороженной панели воспроизвели ",
    "заранее ожидаемое направление."
  ),
  
  paste0(
    "В заранее заданном сравнении sepsis/septic shock против ",
    "uncomplicated infection разделение было умеренным ",
    "(AUC = ",
    sprintf(
      "%.3f",
      row_141_primary$auc_fixed_direction
    ),
    "; p = ",
    format.pval(
      row_141_primary$p_value,
      digits = 4
    ),
    "), однако score был выше при septic shock."
  ),
  
  "",
  
  "GSE185263",
  
  paste0(
    "Основная заранее заданная внешняя конечная точка была подтверждена: ",
    "пятигенный score коррелировал с 24-часовым SOFA ",
    "(rho = ",
    sprintf(
      "%.3f",
      primary_142$rho
    ),
    ", p = ",
    format.pval(
      primary_142$p_value,
      digits = 4
    ),
    "; n = ",
    primary_142$n,
    ")."
  ),
  
  paste0(
    "Все пять отдельных генов воспроизвели ожидаемое направление ",
    "связи с SOFA и сохранили значимость после FDR-коррекции."
  ),
  
  paste0(
    "После поправки на возраст, пол и географическое место набора ",
    "связь с SOFA сохранялась (beta = ",
    sprintf(
      "%.3f",
      adjusted_142$estimate
    ),
    "; p = ",
    format.pval(
      adjusted_142$p_value,
      digits = 4
    ),
    ")."
  ),
  
  "",
  
  "ГЕОГРАФИЧЕСКАЯ ВОСПРОИЗВОДИМОСТЬ",
  
  paste0(
    "Во всех пяти географических подкогортах направление связи ",
    "score-SOFA было положительным. Random-effects pooled rho = ",
    sprintf(
      "%.3f",
      location_meta_result$random_rho
    ),
    " (95% ДИ ",
    sprintf(
      "%.3f",
      location_meta_result$random_CI_low
    ),
    "-",
    sprintf(
      "%.3f",
      location_meta_result$random_CI_high
    ),
    "), I2 = ",
    sprintf(
      "%.1f",
      location_meta_result$I2_percent
    ),
    "%."
  ),
  
  "",
  
  "ИНТЕГРИРОВАННАЯ ИНТЕРПРЕТАЦИЯ",
  
  paste0(
    "Совокупность результатов исходной и независимых внешних когорт ",
    "поддерживает интерпретацию пятигенной сигнатуры как компактного ",
    "и воспроизводимого молекулярного отражения myeloid-adaptive ",
    "host-response axis, связанного с транскриптомными эндотипами ",
    "и тяжестью органной дисфункции."
  ),
  
  paste0(
    "Эти данные не являются валидацией готового диагностического ",
    "или прогностического клинического теста."
  )
)


summary_ru_file <- file.path(
  text_dir,
  "143_multicohort_evidence_summary_RU.txt"
)


writeLines(
  summary_ru,
  con = summary_ru_file,
  useBytes = TRUE
)


# ==============================================================================
# 35. PROPOSED MANUSCRIPT RESULTS PARAGRAPH
# ==============================================================================

results_paragraph <- c(
  
  "PROPOSED MANUSCRIPT RESULTS PARAGRAPH",
  
  "====================================================================",
  
  "",
  
  paste0(
    "The frozen five-gene signature was subsequently evaluated across ",
    "two independent whole-blood RNA-seq cohorts. In GSE154918, all five ",
    "component genes reproduced their prespecified direction of change. ",
    "The prespecified comparison of sepsis or septic shock with ",
    "uncomplicated infection showed modest discrimination (AUC ",
    sprintf(
      "%.3f",
      row_141_primary$auc_fixed_direction
    ),
    ", 95% CI ",
    sprintf(
      "%.3f",
      row_141_primary$auc_ci_low
    ),
    "-",
    sprintf(
      "%.3f",
      row_141_primary$auc_ci_high
    ),
    "; p = ",
    format.pval(
      row_141_primary$p_value,
      digits = 3
    ),
    "), whereas patients with septic shock showed higher scores than ",
    "those with uncomplicated infection (AUC ",
    sprintf(
      "%.3f",
      row_141_shock$auc_fixed_direction
    ),
    "; p = ",
    format.pval(
      row_141_shock$p_value,
      digits = 3
    ),
    ")."
  ),
  
  "",
  
  paste0(
    "In the larger GSE185263 cohort, the prespecified primary external ",
    "severity endpoint was met. Among ",
    primary_142$n,
    " sepsis samples with available 24-h SOFA data, the five-gene score ",
    "was positively associated with organ-dysfunction severity ",
    "(Spearman rho = ",
    sprintf(
      "%.3f",
      primary_142$rho
    ),
    ", p = ",
    format.pval(
      primary_142$p_value,
      digits = 3
    ),
    "). Each of the five component genes independently reproduced its ",
    "expected direction of association with SOFA and remained significant ",
    "after gene-level false-discovery-rate correction. The score-SOFA ",
    "association persisted after adjustment for age, sex and geographic ",
    "collection location (beta = ",
    sprintf(
      "%.3f",
      adjusted_142$estimate
    ),
    " score units per SOFA point, p = ",
    format.pval(
      adjusted_142$p_value,
      digits = 3
    ),
    ")."
  ),
  
  "",
  
  paste0(
    "The direction of the score-SOFA relationship was positive in all ",
    "five geographic subcohorts. Random-effects Fisher-z synthesis of ",
    "the location-specific correlations yielded a pooled rho of ",
    sprintf(
      "%.3f",
      location_meta_result$random_rho
    ),
    " (95% CI ",
    sprintf(
      "%.3f",
      location_meta_result$random_CI_low
    ),
    "-",
    sprintf(
      "%.3f",
      location_meta_result$random_CI_high
    ),
    "), with I2 = ",
    sprintf(
      "%.1f",
      location_meta_result$I2_percent
    ),
    "%. Together with the discovery-cohort associations with SRS and CTS, ",
    "these findings support the five-gene score as a compact molecular ",
    "representation of a reproducible host-response axis rather than as ",
    "a simple binary diagnostic classifier."
  )
)


results_paragraph_file <- file.path(
  text_dir,
  "143_proposed_external_validation_results_paragraph_EN.txt"
)


writeLines(
  results_paragraph,
  con = results_paragraph_file,
  useBytes = TRUE
)


# ==============================================================================
# 36. INPUT MANIFEST
# ==============================================================================

manifest <- tibble::tibble(
  
  item = c(
    "Script135 scores/SRS/CTS",
    "Script141 scores",
    "Script141 comparisons",
    "Script141 gene direction",
    "Script142b sample scores",
    "Script142b primary SOFA",
    "Script142b secondary",
    "Script142b gene-SOFA",
    "Script142b adjusted model",
    "Script142b locations"
  ),
  
  path = c(
    file_135,
    file_141_scores,
    file_141_comparisons,
    file_141_direction,
    file_142_scores,
    file_142_primary,
    file_142_secondary,
    file_142_gene_sofa,
    file_142_adjusted,
    file_142_location
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
    "143_input_manifest.csv"
  ),
  
  row.names = FALSE
)


# ==============================================================================
# 37. SESSION INFO
# ==============================================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      logs_dir,
      "143_sessionInfo.txt"
    )
)


# ==============================================================================
# 38. FINAL REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 143 completed successfully.\n")
cat("====================================================================\n\n")


cat(
  "DISCOVERY ENDOTYPE EVIDENCE:\n"
)


cat(
  "SRSq rho = ",
  sprintf(
    "%.3f",
    srsq_result_135$rho
  ),
  "; p = ",
  format.pval(
    srsq_result_135$p_value,
    digits = 4
  ),
  "\n",
  sep = ""
)


cat(
  "CTS epsilon2 = ",
  sprintf(
    "%.3f",
    cts_result_135$epsilon2
  ),
  "; p = ",
  format.pval(
    cts_result_135$p_value,
    digits = 4
  ),
  "\n\n",
  sep = ""
)


cat(
  "GSE154918:\n"
)


cat(
  "Expected gene directions = ",
  sum(
    direction_concordance_summary$GSE154918_concordant,
    na.rm = TRUE
  ),
  "/5\n",
  sep = ""
)


cat(
  "Primary sepsis/shock vs uncomplicated AUC = ",
  sprintf(
    "%.3f",
    row_141_primary$auc_fixed_direction
  ),
  "; p = ",
  format.pval(
    row_141_primary$p_value,
    digits = 4
  ),
  "\n",
  sep = ""
)


cat(
  "Shock vs uncomplicated AUC = ",
  sprintf(
    "%.3f",
    row_141_shock$auc_fixed_direction
  ),
  "; p = ",
  format.pval(
    row_141_shock$p_value,
    digits = 4
  ),
  "\n\n",
  sep = ""
)


cat(
  "GSE185263:\n"
)


cat(
  "Primary score vs SOFA rho = ",
  sprintf(
    "%.3f",
    primary_142$rho
  ),
  "; p = ",
  format.pval(
    primary_142$p_value,
    digits = 4
  ),
  "\n",
  sep = ""
)


cat(
  "Expected gene-SOFA directions = ",
  sum(
    direction_concordance_summary$GSE185263_concordant,
    na.rm = TRUE
  ),
  "/5\n",
  sep = ""
)


cat(
  "Adjusted SOFA beta = ",
  sprintf(
    "%.3f",
    adjusted_142$estimate
  ),
  "; p = ",
  format.pval(
    adjusted_142$p_value,
    digits = 4
  ),
  "\n\n",
  sep = ""
)


cat(
  "GSE185263 GEOGRAPHIC RANDOM-EFFECTS META-ANALYSIS:\n"
)


cat(
  "k = ",
  location_meta_result$k,
  "\n",
  sep = ""
)


cat(
  "Random-effects pooled rho = ",
  sprintf(
    "%.3f",
    location_meta_result$random_rho
  ),
  " [",
  sprintf(
    "%.3f",
    location_meta_result$random_CI_low
  ),
  ", ",
  sprintf(
    "%.3f",
    location_meta_result$random_CI_high
  ),
  "]\n",
  sep = ""
)


cat(
  "p = ",
  format.pval(
    location_meta_result$random_p,
    digits = 4
  ),
  "\n",
  sep = ""
)


cat(
  "Q = ",
  sprintf(
    "%.3f",
    location_meta_result$Q
  ),
  "; df = ",
  location_meta_result$Q_df,
  "; heterogeneity p = ",
  format.pval(
    location_meta_result$Q_p,
    digits = 4
  ),
  "\n",
  sep = ""
)


cat(
  "I2 = ",
  sprintf(
    "%.1f",
    location_meta_result$I2_percent
  ),
  "%\n",
  sep = ""
)


cat(
  "tau2 = ",
  signif(
    location_meta_result$tau2_fisher_z,
    4
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
  "Main multicohort figure:\n"
)


cat(
  normalizePath(
    file.path(
      figures_dir,
      "143_MAIN_Figure_multicohort_five_gene_evidence.png"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "English evidence summary:\n"
)


cat(
  normalizePath(
    summary_en_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "Russian evidence summary:\n"
)


cat(
  normalizePath(
    summary_ru_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "Manuscript-ready Results paragraph:\n"
)


cat(
  normalizePath(
    results_paragraph_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "IMPORTANT INTERPRETATION:\n"
)


cat(
  "- GSE154918 and GSE185263 have different prespecified endpoints.\n"
)


cat(
  "- They are integrated narratively, not pooled into one effect estimate.\n"
)


cat(
  paste0(
    "- Random-effects pooling is performed only across geographic ",
    "score-SOFA correlations within GSE185263.\n"
  )
)


cat(
  "- No new feature selection was performed.\n"
)


cat(
  "- No gene substitution was performed.\n"
)


cat(
  "- No score refitting was performed.\n"
)


cat(
  paste0(
    "- The evidence supports a reproducible host-response axis, ",
    "not a calibrated clinical diagnostic/prognostic assay.\n\n"
  )
)


cat(
  "Done.\n"
)