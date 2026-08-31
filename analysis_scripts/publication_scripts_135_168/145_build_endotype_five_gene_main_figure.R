################################################################################
# Script 145 v3
# Main publication figure:
# Blood endotype hierarchy and five-gene host-response score
#
# Project: Sepsis_DESeq2
#
# Purpose:
#   Build the central publication figure demonstrating that:
#
#   1) SRS and CTS define a hierarchical structure of blood host-response states;
#   2) the biology-guided five-gene score differs across SRS and CTS;
#   3) integrated CTS/SRS groups form an ordered host-response continuum;
#   4) the five component genes reproduce the expected
#      myeloid-up / adaptive-down expression pattern;
#   5) the five-gene score tracks quantitative SRSq.
#
# Current manuscript primary signature:
#
#   UP:
#     CD177
#     HK3
#     IRAK3
#
#   DOWN:
#     CARD11
#     IKZF2
#
#   Score:
#     mean[z(CD177), z(HK3), z(IRAK3)]
#       -
#     mean[z(CARD11), z(IKZF2)]
#
# IMPORTANT:
#   - This script DOES NOT reassign SRS or CTS.
#   - This script DOES NOT perform feature selection.
#   - The five-gene score is read from the finalized Script 135 table.
#   - Expression recalculation from counts is used for HEATMAP VISUALIZATION ONLY.
#   - Historical erroneous SRS branches (127/127b) are explicitly rejected
#     by distribution QC.
#
# Expected final BP distributions:
#   SRS1 = 28
#   SRS2 = 7
#
#   CTS1 = 14
#   CTS2 = 6
#   CTS3 = 15
#
# Output:
#   results/blood_endotypes_biomarkers/145_main_endotype_figure/
#
################################################################################


cat("====================================================================\n")
cat("Running Script 145 v3\n")
cat("Main endotype + five-gene host-response publication figure\n")
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
# 2. OUTPUT DIRECTORIES
# =============================================================================

out_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "145_main_endotype_figure"
)

table_dir <- file.path(
  out_dir,
  "tables"
)

figure_dir <- file.path(
  out_dir,
  "figures"
)

text_dir <- file.path(
  out_dir,
  "text"
)

for (d in c(
  out_dir,
  table_dir,
  figure_dir,
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
# 3. PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "scales",
  "openxlsx",
  "readr"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    stop(
      "Required package is not installed: ",
      pkg
    )
  }
}

suppressPackageStartupMessages({
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(openxlsx)
  library(readr)
})


package_versions <- data.frame(
  
  package = c(
    "R",
    required_packages
  ),
  
  version = c(
    
    paste0(
      R.version$major,
      ".",
      strsplit(
        R.version$minor,
        "\\."
      )[[1]][1]
    ),
    
    vapply(
      required_packages,
      function(x) {
        
        as.character(
          packageVersion(x)
        )
        
      },
      character(1)
    )
  ),
  
  stringsAsFactors = FALSE
)


cat("\nPackage versions:\n")
print(package_versions)


write.csv(
  package_versions,
  file.path(
    table_dir,
    "145_package_versions.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file = file.path(
    text_dir,
    "145_sessionInfo.txt"
  )
)


# =============================================================================
# 4. SETTINGS
# =============================================================================

primary_genes_up <- c(
  "CD177",
  "HK3",
  "IRAK3"
)

primary_genes_down <- c(
  "CARD11",
  "IKZF2"
)

primary_genes <- c(
  primary_genes_up,
  primary_genes_down
)


expected_srs <- c(
  SRS1 = 28,
  SRS2 = 7
)

expected_cts <- c(
  CTS1 = 14,
  CTS2 = 6,
  CTS3 = 15
)


# =============================================================================
# 5. PUBLICATION COLORS
# =============================================================================

col_srs <- c(
  SRS1 = "#B2182B",
  SRS2 = "#2166AC"
)

col_cts <- c(
  CTS1 = "#D6604D",
  CTS2 = "#FDB863",
  CTS3 = "#4393C3"
)

col_integrated <- c(
  "CTS1/SRS1" = "#B2182B",
  "CTS2/SRS1" = "#EF8A62",
  "CTS3/SRS1" = "#67A9CF",
  "CTS3/SRS2" = "#2166AC"
)

col_expression <- c(
  "#2166AC",
  "#F7F7F7",
  "#B2182B"
)


# =============================================================================
# 6. HELPER FUNCTIONS
# =============================================================================

clean_names_simple <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    x
  )
  
  x
}


find_exact_file_recursive <- function(
    root,
    basename_pattern
) {
  
  hits <- list.files(
    root,
    pattern = basename_pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  hits <- hits[
    file.exists(hits)
  ]
  
  if (length(hits) == 0) {
    return(NA_character_)
  }
  
  if (length(hits) > 1) {
    
    info <- file.info(hits)
    
    hits <- hits[
      order(
        info$mtime,
        decreasing = TRUE
      )
    ]
    
    warning(
      "Multiple files matched pattern:\n",
      basename_pattern,
      "\nUsing most recently modified:\n",
      hits[1]
    )
  }
  
  hits[1]
}


read_table_auto <- function(file) {
  
  ext <- tolower(
    tools::file_ext(file)
  )
  
  if (ext == "csv") {
    
    return(
      readr::read_csv(
        file,
        show_col_types = FALSE
      ) %>%
        as.data.frame()
    )
  }
  
  if (ext %in% c(
    "xlsx",
    "xls"
  )) {
    
    return(
      openxlsx::read.xlsx(
        file
      )
    )
  }
  
  stop(
    "Unsupported file format: ",
    file
  )
}


find_column <- function(
    df,
    exact = character(0),
    regex = character(0),
    required = TRUE,
    label = "column"
) {
  
  nms <- names(df)
  
  for (candidate in exact) {
    
    hit <- nms[
      tolower(nms) ==
        tolower(candidate)
    ]
    
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  
  
  for (pattern in regex) {
    
    hit <- nms[
      grepl(
        pattern,
        nms,
        ignore.case = TRUE
      )
    ]
    
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  
  
  if (required) {
    
    stop(
      "Could not identify ",
      label,
      ".\nAvailable columns:\n",
      paste(
        nms,
        collapse = ", "
      )
    )
  }
  
  NA_character_
}


normalize_srs <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  out <- rep(
    NA_character_,
    length(x)
  )
  
  out[
    grepl(
      "SRS.?1|^1$",
      x
    )
  ] <- "SRS1"
  
  out[
    grepl(
      "SRS.?2|^2$",
      x
    )
  ] <- "SRS2"
  
  out[
    grepl(
      "SRS.?3|^3$",
      x
    )
  ] <- "SRS3"
  
  out
}


normalize_cts <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  out <- rep(
    NA_character_,
    length(x)
  )
  
  out[
    grepl(
      "CTS.?1|^1$",
      x
    )
  ] <- "CTS1"
  
  out[
    grepl(
      "CTS.?2|^2$",
      x
    )
  ] <- "CTS2"
  
  out[
    grepl(
      "CTS.?3|^3$",
      x
    )
  ] <- "CTS3"
  
  out
}


format_p <- function(p) {
  
  if (is.na(p)) {
    return("NA")
  }
  
  if (p < 2.2e-16) {
    return("<2.2×10⁻¹⁶")
  }
  
  if (p < 0.001) {
    
    return(
      format(
        p,
        scientific = TRUE,
        digits = 2
      )
    )
  }
  
  format(
    round(
      p,
      3
    ),
    nsmall = 3
  )
}


epsilon_squared_kw <- function(
    x,
    group
) {
  
  dat <- data.frame(
    x = x,
    group = group
  )
  
  dat <- dat[
    complete.cases(dat),
    ,
    drop = FALSE
  ]
  
  kw <- kruskal.test(
    x ~ group,
    data = dat
  )
  
  H <- as.numeric(
    kw$statistic
  )
  
  k <- length(
    unique(dat$group)
  )
  
  n <- nrow(dat)
  
  eps <- (
    H - k + 1
  ) / (
    n - k
  )
  
  max(
    0,
    eps
  )
}


theme_publication <- function(
    base_size = 10
) {
  
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    
    ggplot2::theme(
      
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 1
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.5
      ),
      
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      
      legend.title = ggplot2::element_blank(),
      
      legend.position = "top",
      
      strip.background = ggplot2::element_blank(),
      
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      
      plot.margin = ggplot2::margin(
        7,
        7,
        7,
        7
      )
    )
}


# =============================================================================
# 7. LOCATE FINAL SCRIPT 135 SCORE TABLE
# =============================================================================

score_file <- find_exact_file_recursive(
  "results",
  "^135_sepsis_blood_scores_with_SRS_CTS\\.csv$"
)


if (is.na(score_file)) {
  
  stop(
    paste0(
      "\nCould not find:\n",
      "135_sepsis_blood_scores_with_SRS_CTS.csv\n\n",
      "Do NOT substitute an older Script 127/130 table automatically.\n",
      "The main figure must use the finalized Script 135 score table."
    )
  )
}


cat("\nFinal Script 135 score table:\n")

print(
  normalizePath(
    score_file,
    winslash = "\\",
    mustWork = TRUE
  )
)


# =============================================================================
# 8. READ FINAL SCORE TABLE
# =============================================================================

score_raw <- read_table_auto(
  score_file
)

names(score_raw) <- clean_names_simple(
  names(score_raw)
)


cat("\nScript 135 table dimensions:\n")
print(dim(score_raw))

cat("\nScript 135 columns:\n")
print(names(score_raw))


# =============================================================================
# 9. IDENTIFY REQUIRED COLUMNS
# =============================================================================

sample_col <- find_column(
  
  score_raw,
  
  exact = c(
    "sample_id",
    "Sample_ID",
    "sample"
  ),
  
  regex = c(
    "^sample.*id$"
  ),
  
  label = "sample ID column"
)


score_col <- find_column(
  
  score_raw,
  
  exact = c(
    "primary_score",
    "Primary_5_gene",
    "Primary_5_gene_score",
    "primary_5_gene_score",
    "primary_5gene_score",
    "five_gene_score"
  ),
  
  regex = c(
    "primary.*5.*gene.*score",
    "primary.*5gene.*score",
    "primary.*score",
    "five.*gene.*score"
  ),
  
  label = "primary five-gene score column"
)


srs_col <- find_column(
  
  score_raw,
  
  exact = c(
    "SRS",
    "SRS_class",
    "srs_label",
    "corrected_SRS_endotype"
  ),
  
  regex = c(
    "^SRS$",
    "SRS.*class",
    "SRS.*label",
    "corrected.*SRS"
  ),
  
  label = "SRS class column"
)


srsq_col <- find_column(
  
  score_raw,
  
  exact = c(
    "SRSq",
    "srsq"
  ),
  
  regex = c(
    "^SRSq$",
    "SRS.*q"
  ),
  
  label = "SRSq column"
)


cts_col <- find_column(
  
  score_raw,
  
  exact = c(
    "CTS",
    "CTS_class",
    "CTS_endotype",
    "cts_label"
  ),
  
  regex = c(
    "^CTS$",
    "CTS.*class",
    "CTS.*endotype",
    "CTS.*label"
  ),
  
  label = "CTS class column"
)


cat("\nDetected columns:\n")
cat("sample_id:", sample_col, "\n")
cat("five-gene score:", score_col, "\n")
cat("SRS:", srs_col, "\n")
cat("SRSq:", srsq_col, "\n")
cat("CTS:", cts_col, "\n")


# =============================================================================
# 10. BUILD FINAL BP DATASET
# =============================================================================

blood <- score_raw %>%
  
  dplyr::transmute(
    
    sample_id =
      as.character(
        .data[[sample_col]]
      ),
    
    five_gene_score =
      as.numeric(
        .data[[score_col]]
      ),
    
    SRS =
      normalize_srs(
        .data[[srs_col]]
      ),
    
    SRSq =
      as.numeric(
        .data[[srsq_col]]
      ),
    
    CTS =
      normalize_cts(
        .data[[cts_col]]
      )
  ) %>%
  
  dplyr::filter(
    grepl(
      "^BP",
      sample_id,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::distinct(
    sample_id,
    .keep_all = TRUE
  )


if (nrow(blood) != 35) {
  
  stop(
    "Expected 35 BP samples in final Script 135 table; found ",
    nrow(blood),
    "."
  )
}


if (anyNA(
  blood[
    ,
    c(
      "five_gene_score",
      "SRS",
      "SRSq",
      "CTS"
    )
  ]
)) {
  
  stop(
    "Missing five-gene score / SRS / SRSq / CTS values detected."
  )
}


blood <- blood %>%
  
  dplyr::mutate(
    
    SRS = factor(
      SRS,
      levels = c(
        "SRS1",
        "SRS2"
      )
    ),
    
    CTS = factor(
      CTS,
      levels = c(
        "CTS1",
        "CTS2",
        "CTS3"
      )
    )
  )


# =============================================================================
# 11. CRITICAL ENDOTYPE DISTRIBUTION QC
# =============================================================================

observed_srs <- table(
  blood$SRS
)

observed_cts <- table(
  blood$CTS
)


cat("\nFinal SRS distribution:\n")
print(observed_srs)

cat("\nFinal CTS distribution:\n")
print(observed_cts)


if (!all(
  as.integer(
    observed_srs[
      names(expected_srs)
    ]
  ) ==
  as.integer(expected_srs)
)) {
  
  stop(
    paste0(
      "\nSRS DISTRIBUTION QC FAILED.\n",
      "Expected SRS1=28, SRS2=7.\n",
      "This may indicate that an obsolete/incorrect SRS source was loaded.\n",
      "Figure generation stopped deliberately."
    )
  )
}


if (!all(
  as.integer(
    observed_cts[
      names(expected_cts)
    ]
  ) ==
  as.integer(expected_cts)
)) {
  
  stop(
    paste0(
      "\nCTS DISTRIBUTION QC FAILED.\n",
      "Expected CTS1=14, CTS2=6, CTS3=15.\n",
      "Figure generation stopped deliberately."
    )
  )
}


cat("\nEndotype distribution QC: PASS\n")


# =============================================================================
# 12. INTEGRATED CTS/SRS GROUP
# =============================================================================

blood <- blood %>%
  
  dplyr::mutate(
    
    integrated_group = paste0(
      as.character(CTS),
      "/",
      as.character(SRS)
    )
  )


expected_integrated <- c(
  "CTS1/SRS1",
  "CTS2/SRS1",
  "CTS3/SRS1",
  "CTS3/SRS2"
)


unexpected_integrated <- setdiff(
  unique(
    blood$integrated_group
  ),
  expected_integrated
)


if (length(unexpected_integrated) > 0) {
  
  warning(
    "Unexpected integrated groups detected: ",
    paste(
      unexpected_integrated,
      collapse = ", "
    )
  )
}


blood <- blood %>%
  
  dplyr::mutate(
    
    integrated_group = factor(
      integrated_group,
      levels = expected_integrated
    )
  )


cat("\nIntegrated CTS/SRS distribution:\n")

print(
  table(
    blood$integrated_group,
    useNA = "ifany"
  )
)


# =============================================================================
# 13. STATISTICS FOR FIGURE
# =============================================================================

srs_test <- wilcox.test(
  five_gene_score ~ SRS,
  data = blood,
  exact = FALSE
)


cts_test <- kruskal.test(
  five_gene_score ~ CTS,
  data = blood
)


cts_epsilon <- epsilon_squared_kw(
  blood$five_gene_score,
  blood$CTS
)


integrated_test <- kruskal.test(
  five_gene_score ~ integrated_group,
  data = blood
)


integrated_epsilon <- epsilon_squared_kw(
  blood$five_gene_score,
  blood$integrated_group
)


srsq_cor <- cor.test(
  blood$five_gene_score,
  blood$SRSq,
  method = "spearman",
  exact = FALSE
)


stats_summary <- data.frame(
  
  analysis = c(
    "Five-gene score by SRS",
    "Five-gene score by CTS",
    "Five-gene score by integrated CTS/SRS",
    "Five-gene score vs SRSq"
  ),
  
  test = c(
    "Wilcoxon rank-sum",
    "Kruskal-Wallis",
    "Kruskal-Wallis",
    "Spearman correlation"
  ),
  
  statistic = c(
    
    unname(
      srs_test$statistic
    ),
    
    unname(
      cts_test$statistic
    ),
    
    unname(
      integrated_test$statistic
    ),
    
    unname(
      srsq_cor$estimate
    )
  ),
  
  p_value = c(
    srs_test$p.value,
    cts_test$p.value,
    integrated_test$p.value,
    srsq_cor$p.value
  ),
  
  effect_size = c(
    NA_real_,
    cts_epsilon,
    integrated_epsilon,
    unname(
      srsq_cor$estimate
    )
  ),
  
  effect_size_type = c(
    NA_character_,
    "epsilon_squared",
    "epsilon_squared",
    "Spearman_rho"
  ),
  
  stringsAsFactors = FALSE
)


cat("\nFigure statistics:\n")
print(stats_summary)


write.csv(
  stats_summary,
  file.path(
    table_dir,
    "145_main_figure_statistics.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 14. SUMMARY TABLES
# =============================================================================

group_summary <- dplyr::bind_rows(
  
  
  blood %>%
    
    dplyr::group_by(SRS) %>%
    
    dplyr::summarise(
      
      n = dplyr::n(),
      
      median_score = median(
        five_gene_score,
        na.rm = TRUE
      ),
      
      IQR_score = IQR(
        five_gene_score,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    dplyr::transmute(
      
      framework = "SRS",
      
      group = as.character(SRS),
      
      n = n,
      
      median_score = median_score,
      
      IQR_score = IQR_score
    ),
  
  
  blood %>%
    
    dplyr::group_by(CTS) %>%
    
    dplyr::summarise(
      
      n = dplyr::n(),
      
      median_score = median(
        five_gene_score,
        na.rm = TRUE
      ),
      
      IQR_score = IQR(
        five_gene_score,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    dplyr::transmute(
      
      framework = "CTS",
      
      group = as.character(CTS),
      
      n = n,
      
      median_score = median_score,
      
      IQR_score = IQR_score
    ),
  
  
  blood %>%
    
    dplyr::group_by(
      integrated_group
    ) %>%
    
    dplyr::summarise(
      
      n = dplyr::n(),
      
      median_score = median(
        five_gene_score,
        na.rm = TRUE
      ),
      
      IQR_score = IQR(
        five_gene_score,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    dplyr::transmute(
      
      framework = "CTS/SRS",
      
      group = as.character(
        integrated_group
      ),
      
      n = n,
      
      median_score = median_score,
      
      IQR_score = IQR_score
    )
)


write.csv(
  group_summary,
  file.path(
    table_dir,
    "145_endotype_score_group_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 15. PANEL A
# CTS x SRS HIERARCHICAL OVERLAP
# =============================================================================

cross_tab <- blood %>%
  
  dplyr::count(
    CTS,
    SRS,
    name = "n"
  ) %>%
  
  tidyr::complete(
    CTS,
    SRS,
    fill = list(
      n = 0
    )
  ) %>%
  
  dplyr::group_by(SRS) %>%
  
  dplyr::mutate(
    
    SRS_total = sum(n),
    
    row_percent =
      100 *
      n /
      SRS_total
  ) %>%
  
  dplyr::ungroup()


write.csv(
  cross_tab,
  file.path(
    table_dir,
    "145_CTS_by_SRS_cross_table.csv"
  ),
  row.names = FALSE
)


p_A <- ggplot2::ggplot(
  
  cross_tab,
  
  ggplot2::aes(
    x = CTS,
    y = SRS,
    fill = n
  )
  
) +
  
  ggplot2::geom_tile(
    linewidth = 1.2,
    color = "white"
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = paste0(
        "n=",
        n,
        "\n",
        round(
          row_percent,
          0
        ),
        "%"
      )
    ),
    
    size = 3.4,
    fontface = "bold"
  ) +
  
  ggplot2::scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2166AC"
  ) +
  
  theme_publication(
    base_size = 10
  ) +
  
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  ) +
  
  ggplot2::labs(
    
    title =
      "Hierarchical overlap of SRS and CTS",
    
    subtitle =
      "CTS1 and CTS2 occur within the SRS1-dominant state",
    
    x =
      "Consensus Transcriptomic Subtype",
    
    y =
      "Sepsis Response Signature",
    
    fill =
      "Patients"
  )


# =============================================================================
# 16. PANEL B
# FIVE-GENE SCORE BY SRS
# =============================================================================

p_B <- ggplot2::ggplot(
  
  blood,
  
  ggplot2::aes(
    x = SRS,
    y = five_gene_score,
    fill = SRS
  )
  
) +
  
  ggplot2::geom_violin(
    width = 0.85,
    alpha = 0.25,
    trim = FALSE,
    color = NA
  ) +
  
  ggplot2::geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    alpha = 0.70,
    linewidth = 0.45
  ) +
  
  ggplot2::geom_jitter(
    width = 0.10,
    size = 2.0,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_srs
  ) +
  
  theme_publication(
    base_size = 10
  ) +
  
  ggplot2::theme(
    legend.position = "none"
  ) +
  
  ggplot2::labs(
    
    title =
      "Five-gene score separates SRS states",
    
    subtitle =
      paste0(
        "Wilcoxon p=",
        format_p(
          srs_test$p.value
        )
      ),
    
    x = NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 17. PANEL C
# FIVE-GENE SCORE BY CTS
# =============================================================================

p_C <- ggplot2::ggplot(
  
  blood,
  
  ggplot2::aes(
    x = CTS,
    y = five_gene_score,
    fill = CTS
  )
  
) +
  
  ggplot2::geom_violin(
    width = 0.85,
    alpha = 0.25,
    trim = FALSE,
    color = NA
  ) +
  
  ggplot2::geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    alpha = 0.70,
    linewidth = 0.45
  ) +
  
  ggplot2::geom_jitter(
    width = 0.10,
    size = 2.0,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_cts
  ) +
  
  theme_publication(
    base_size = 10
  ) +
  
  ggplot2::theme(
    legend.position = "none"
  ) +
  
  ggplot2::labs(
    
    title =
      "The same molecular axis spans CTS classes",
    
    subtitle =
      paste0(
        "Kruskal-Wallis p=",
        format_p(
          cts_test$p.value
        ),
        "; ε²=",
        format(
          round(
            cts_epsilon,
            3
          ),
          nsmall = 3
        )
      ),
    
    x = NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 18. PANEL D
# INTEGRATED CTS/SRS CONTINUUM
# =============================================================================

p_D <- ggplot2::ggplot(
  
  blood,
  
  ggplot2::aes(
    x = integrated_group,
    y = five_gene_score,
    fill = integrated_group
  )
  
) +
  
  ggplot2::geom_violin(
    width = 0.85,
    alpha = 0.22,
    trim = FALSE,
    color = NA
  ) +
  
  ggplot2::geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    alpha = 0.72,
    linewidth = 0.45
  ) +
  
  ggplot2::geom_jitter(
    width = 0.10,
    size = 2.0,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_integrated,
    drop = FALSE
  ) +
  
  theme_publication(
    base_size = 9
  ) +
  
  ggplot2::theme(
    
    legend.position = "none",
    
    axis.text.x = ggplot2::element_text(
      angle = 25,
      hjust = 1
    )
  ) +
  
  ggplot2::labs(
    
    title =
      "Integrated endotypes reveal an ordered continuum",
    
    subtitle =
      paste0(
        "Kruskal-Wallis p=",
        format_p(
          integrated_test$p.value
        ),
        "; ε²=",
        format(
          round(
            integrated_epsilon,
            3
          ),
          nsmall = 3
        )
      ),
    
    x = NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 19. LOAD RAW COUNTS FOR FIVE-GENE HEATMAP
#
# VISUALIZATION ONLY.
# Primary score remains the finalized Script 135 score.
# =============================================================================

counts_csv <- file.path(
  "data",
  "counts_all.csv"
)

counts_xlsx <- file.path(
  "data",
  "counts_all.xlsx"
)


if (file.exists(counts_csv)) {
  
  counts_raw <- readr::read_csv(
    counts_csv,
    show_col_types = FALSE
  ) %>%
    as.data.frame()
  
  counts_source <- counts_csv
  
} else if (file.exists(counts_xlsx)) {
  
  counts_raw <- openxlsx::read.xlsx(
    counts_xlsx
  )
  
  counts_source <- counts_xlsx
  
} else {
  
  stop(
    "Neither data/counts_all.csv nor data/counts_all.xlsx was found."
  )
}


cat("\nExpression source for heatmap:\n")

print(
  normalizePath(
    counts_source,
    winslash = "\\",
    mustWork = TRUE
  )
)


gene_col <- names(
  counts_raw
)[1]


counts_raw[[gene_col]] <- as.character(
  counts_raw[[gene_col]]
)


# Historical alias harmonization
counts_raw[[gene_col]][
  counts_raw[[gene_col]] == "EMR3"
] <- "ADGRE3"


sample_numeric_cols <- names(counts_raw)[
  vapply(
    counts_raw,
    is.numeric,
    logical(1)
  )
]


counts_clean <- counts_raw %>%
  
  dplyr::rename(
    gene_symbol =
      dplyr::all_of(gene_col)
  ) %>%
  
  dplyr::group_by(
    gene_symbol
  ) %>%
  
  dplyr::summarise(
    
    dplyr::across(
      dplyr::all_of(
        sample_numeric_cols
      ),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),
    
    .groups = "drop"
  )


missing_panel_genes <- setdiff(
  primary_genes,
  counts_clean$gene_symbol
)


if (length(
  missing_panel_genes
) > 0) {
  
  stop(
    "Five-gene panel genes missing from counts matrix: ",
    paste(
      missing_panel_genes,
      collapse = ", "
    )
  )
}


# =============================================================================
# 20. RPM -> log2(RPM + 1)
# =============================================================================

count_matrix <- counts_clean %>%
  
  tibble::column_to_rownames(
    "gene_symbol"
  ) %>%
  
  as.matrix()


storage.mode(
  count_matrix
) <- "numeric"


bp_samples <- intersect(
  blood$sample_id,
  colnames(
    count_matrix
  )
)


if (length(bp_samples) != 35) {
  
  stop(
    "Expected all 35 BP samples in count matrix; found ",
    length(bp_samples)
  )
}


library_sizes <- colSums(
  count_matrix[
    ,
    bp_samples,
    drop = FALSE
  ],
  na.rm = TRUE
)


if (any(
  library_sizes <= 0
)) {
  
  stop(
    "Zero/negative library size detected."
  )
}


rpm_matrix <- sweep(
  
  count_matrix[
    ,
    bp_samples,
    drop = FALSE
  ],
  
  2,
  
  library_sizes,
  
  "/"
  
) * 1e6


log2rpm_matrix <- log2(
  rpm_matrix + 1
)


panel_expression <- log2rpm_matrix[
  primary_genes,
  bp_samples,
  drop = FALSE
]


panel_z <- t(
  scale(
    t(
      panel_expression
    )
  )
)


panel_z[
  !is.finite(panel_z)
] <- 0


# =============================================================================
# 21. ORDER PATIENTS BY FIVE-GENE SCORE
# =============================================================================

sample_order <- blood %>%
  
  dplyr::arrange(
    dplyr::desc(
      five_gene_score
    )
  ) %>%
  
  dplyr::pull(
    sample_id
  )


panel_z <- panel_z[
  ,
  sample_order,
  drop = FALSE
]


heat_long <- as.data.frame(
  panel_z
) %>%
  
  tibble::rownames_to_column(
    "gene"
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = -gene,
    
    names_to = "sample_id",
    
    values_to = "z_expression"
  ) %>%
  
  dplyr::mutate(
    
    sample_id = factor(
      sample_id,
      levels = sample_order
    ),
    
    gene = factor(
      gene,
      levels = rev(
        primary_genes
      )
    )
  )


# =============================================================================
# 22. HEATMAP ANNOTATION TRACKS
# =============================================================================

annotation_long <- blood %>%
  
  dplyr::select(
    sample_id,
    SRS,
    CTS
  ) %>%
  
  dplyr::filter(
    sample_id %in%
      sample_order
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      SRS,
      CTS
    ),
    
    names_to = "track",
    
    values_to = "class"
  ) %>%
  
  dplyr::mutate(
    
    sample_id = factor(
      sample_id,
      levels = sample_order
    ),
    
    track = factor(
      track,
      levels = c(
        "SRS",
        "CTS"
      )
    )
  )


annotation_colors <- c(
  col_srs,
  col_cts
)


p_E_annotation <- ggplot2::ggplot(
  
  annotation_long,
  
  ggplot2::aes(
    x = sample_id,
    y = track,
    fill = class
  )
  
) +
  
  ggplot2::geom_tile() +
  
  ggplot2::scale_fill_manual(
    values = annotation_colors
  ) +
  
  ggplot2::theme_void(
    base_size = 8
  ) +
  
  ggplot2::theme(
    
    legend.position = "top",
    
    legend.direction = "horizontal",
    
    legend.key.size =
      grid::unit(
        0.30,
        "cm"
      ),
    
    plot.margin = ggplot2::margin(
      0,
      3,
      0,
      3
    )
  ) +
  
  ggplot2::labs(
    fill = NULL
  )


p_E_heatmap <- ggplot2::ggplot(
  
  heat_long,
  
  ggplot2::aes(
    x = sample_id,
    y = gene,
    fill = z_expression
  )
  
) +
  
  ggplot2::geom_tile() +
  
  ggplot2::scale_fill_gradient2(
    
    low = col_expression[1],
    
    mid = col_expression[2],
    
    high = col_expression[3],
    
    midpoint = 0,
    
    limits = c(
      -2.5,
      2.5
    ),
    
    oob = scales::squish
  ) +
  
  ggplot2::theme_classic(
    base_size = 8
  ) +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_blank(),
    
    axis.ticks.x =
      ggplot2::element_blank(),
    
    axis.title =
      ggplot2::element_blank(),
    
    legend.position =
      "right",
    
    panel.grid =
      ggplot2::element_blank(),
    
    plot.margin = ggplot2::margin(
      0,
      3,
      3,
      3
    )
  ) +
  
  ggplot2::labs(
    fill =
      "Gene-wise\nz-score"
  )


p_E <- (
  
  p_E_annotation /
    
    p_E_heatmap
  
) +
  
  patchwork::plot_layout(
    
    heights = c(
      0.25,
      1
    ),
    
    guides = "collect"
  ) +
  
  patchwork::plot_annotation(
    
    title =
      "Five-gene expression across the host-response continuum",
    
    subtitle =
      "Patients ordered from highest to lowest five-gene score"
  ) &
  
  ggplot2::theme(
    
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    
    plot.subtitle = ggplot2::element_text(
      size = 9
    )
  )


# =============================================================================
# 23. PANEL F
# FIVE-GENE SCORE vs SRSq
# =============================================================================

rho_srsq <- unname(
  srsq_cor$estimate
)


p_F <- ggplot2::ggplot(
  
  blood,
  
  ggplot2::aes(
    x = SRSq,
    y = five_gene_score,
    fill = SRS
  )
  
) +
  
  ggplot2::geom_smooth(
    
    method = "lm",
    
    formula = y ~ x,
    
    se = TRUE,
    
    linewidth = 0.7,
    
    color = "grey35",
    
    fill = "grey80"
  ) +
  
  ggplot2::geom_point(
    
    shape = 21,
    
    size = 2.8,
    
    stroke = 0.4,
    
    alpha = 0.90
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_srs
  ) +
  
  theme_publication(
    base_size = 10
  ) +
  
  ggplot2::theme(
    legend.position = "top"
  ) +
  
  ggplot2::labs(
    
    title =
      "Five-gene score tracks quantitative SRSq",
    
    subtitle =
      paste0(
        "Spearman ρ=",
        format(
          round(
            rho_srsq,
            3
          ),
          nsmall = 3
        ),
        "; p=",
        format_p(
          srsq_cor$p.value
        )
      ),
    
    x =
      "Quantitative SRS score (SRSq)",
    
    y =
      "Five-gene host-response score",
    
    fill =
      "SRS"
  )


# =============================================================================
# 24. ASSEMBLE MAIN FIGURE
# =============================================================================

main_figure <- (
  
  (
    p_A |
      p_B |
      p_C
  ) /
    
    (
      p_D |
        p_E |
        p_F
    )
  
) +
  
  patchwork::plot_layout(
    
    heights = c(
      1,
      1.20
    ),
    
    widths = c(
      1,
      1.35,
      1
    ),
    
    guides = "collect"
  ) +
  
  patchwork::plot_annotation(
    
    title =
      "Blood transcriptomic endotypes converge on a five-gene host-response continuum",
    
    subtitle =
      paste0(
        "Primary signature: CD177 + HK3 + IRAK3 − CARD11 − IKZF2; ",
        "n=35 patients with sepsis"
      ),
    
    tag_levels = "A",
    
    theme = ggplot2::theme(
      
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 16
      ),
      
      plot.subtitle = ggplot2::element_text(
        
        size = 11,
        
        margin = ggplot2::margin(
          b = 8
        )
      ),
      
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 14
      )
    )
  )


# =============================================================================
# 25. EXPORT MAIN FIGURE
# =============================================================================

main_png <- file.path(
  figure_dir,
  "145_Figure2_endotype_five_gene_host_response.png"
)

main_pdf <- file.path(
  figure_dir,
  "145_Figure2_endotype_five_gene_host_response.pdf"
)

main_tiff <- file.path(
  figure_dir,
  "145_Figure2_endotype_five_gene_host_response.tiff"
)


ggplot2::ggsave(
  
  filename = main_png,
  
  plot = main_figure,
  
  width = 15,
  
  height = 9.5,
  
  dpi = 600,
  
  bg = "white"
)


ggplot2::ggsave(
  
  filename = main_pdf,
  
  plot = main_figure,
  
  width = 15,
  
  height = 9.5,
  
  device = grDevices::cairo_pdf,
  
  bg = "white"
)


ggplot2::ggsave(
  
  filename = main_tiff,
  
  plot = main_figure,
  
  width = 15,
  
  height = 9.5,
  
  dpi = 600,
  
  compression = "lzw",
  
  bg = "white"
)


# =============================================================================
# 26. EXPORT INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  
  A_CTS_by_SRS =
    p_A,
  
  B_score_by_SRS =
    p_B,
  
  C_score_by_CTS =
    p_C,
  
  D_integrated_continuum =
    p_D,
  
  E_five_gene_heatmap =
    p_E,
  
  F_score_vs_SRSq =
    p_F
)


for (nm in names(
  individual_panels
)) {
  
  p <- individual_panels[[nm]]
  
  
  width_i <- ifelse(
    nm ==
      "E_five_gene_heatmap",
    9,
    6
  )
  
  
  height_i <- ifelse(
    nm ==
      "E_five_gene_heatmap",
    5,
    5
  )
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      figure_dir,
      paste0(
        "145_panel_",
        nm,
        ".png"
      )
    ),
    
    plot = p,
    
    width = width_i,
    
    height = height_i,
    
    dpi = 600,
    
    bg = "white"
  )
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      figure_dir,
      paste0(
        "145_panel_",
        nm,
        ".pdf"
      )
    ),
    
    plot = p,
    
    width = width_i,
    
    height = height_i,
    
    device = grDevices::cairo_pdf,
    
    bg = "white"
  )
}


# =============================================================================
# 27. FIGURE CAPTION
# =============================================================================

caption_en <- paste0(
  
  "Figure 2. Blood transcriptomic endotypes converge on a compact ",
  "five-gene host-response continuum. ",
  
  "(A) Cross-classification of the 35 sepsis blood transcriptomes by ",
  "Sepsis Response Signature (SRS) and Consensus Transcriptomic Subtype ",
  "(CTS). CTS1 and CTS2 occurred exclusively within SRS1, whereas CTS3 ",
  "contained both SRS1 and SRS2 samples, indicating additional biological ",
  "resolution within the SRS1-dominant host-response state. ",
  
  "(B) Distribution of the five-gene host-response score across SRS classes. ",
  
  "(C) Distribution of the score across CTS1, CTS2, and CTS3. ",
  
  "(D) Integrated CTS/SRS categories reveal an ordered molecular continuum ",
  "from CTS1/SRS1 through CTS2/SRS1 and CTS3/SRS1 to CTS3/SRS2. ",
  
  "(E) Gene-wise standardized expression of the five component genes ",
  "(CD177, HK3, IRAK3, CARD11, and IKZF2), with patients ordered by the ",
  "composite score. The signature combines increased myeloid-associated ",
  "expression (CD177, HK3, IRAK3) with reduced adaptive-associated expression ",
  "(CARD11, IKZF2). ",
  
  "(F) Association between the five-gene host-response score and quantitative ",
  "SRSq. Group comparisons were evaluated using Wilcoxon rank-sum or ",
  "Kruskal-Wallis tests as appropriate; continuous association with SRSq ",
  "was evaluated using Spearman correlation. ",
  
  "The five-gene signature was developed independently of SRS and CTS ",
  "assignments."
)


caption_ru <- paste0(
  
  "Рисунок 2. Транскриптомные эндотипы крови сходятся на компактном ",
  "пятигенном континууме ответа хозяина. ",
  
  "(A) Совместная классификация 35 транскриптомов крови пациентов с сепсисом ",
  "по Sepsis Response Signature (SRS) и Consensus Transcriptomic Subtype ",
  "(CTS). CTS1 и CTS2 встречались только внутри SRS1, тогда как CTS3 ",
  "включал как SRS1, так и SRS2, что указывает на дополнительное ",
  "биологическое разрешение внутри доминирующего SRS1-состояния. ",
  
  "(B) Распределение пятигенного host-response score между классами SRS. ",
  
  "(C) Распределение score между CTS1, CTS2 и CTS3. ",
  
  "(D) Интегрированные CTS/SRS-группы формируют упорядоченный молекулярный ",
  "континуум CTS1/SRS1 → CTS2/SRS1 → CTS3/SRS1 → CTS3/SRS2. ",
  
  "(E) Стандартизированная экспрессия пяти генов CD177, HK3, IRAK3, CARD11 ",
  "и IKZF2; пациенты упорядочены по величине composite score. ",
  
  "Signature объединяет повышение myeloid-associated genes ",
  "(CD177, HK3, IRAK3) со снижением adaptive-associated genes ",
  "(CARD11, IKZF2). ",
  
  "(F) Связь пятигенного score с количественным SRSq. ",
  
  "Пятигенная signature была сформирована независимо от SRS и CTS; ",
  "эти классификации не использовались для отбора генов."
)


writeLines(
  caption_en,
  file.path(
    text_dir,
    "145_Figure2_caption_EN.txt"
  )
)


writeLines(
  caption_ru,
  file.path(
    text_dir,
    "145_Figure2_caption_RU.txt"
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
    "145_Figure2_caption_EN_RU.txt"
  )
)


# =============================================================================
# 28. RESULTS-READY NUMERICAL SUMMARY
# =============================================================================

srs_medians <- blood %>%
  
  dplyr::group_by(SRS) %>%
  
  dplyr::summarise(
    
    n = dplyr::n(),
    
    median = median(
      five_gene_score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


cts_medians <- blood %>%
  
  dplyr::group_by(CTS) %>%
  
  dplyr::summarise(
    
    n = dplyr::n(),
    
    median = median(
      five_gene_score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


integrated_medians <- blood %>%
  
  dplyr::group_by(
    integrated_group
  ) %>%
  
  dplyr::summarise(
    
    n = dplyr::n(),
    
    median = median(
      five_gene_score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


results_summary_lines <- c(
  
  "SCRIPT 145 RESULTS SUMMARY",
  "====================================================================",
  "",
  
  paste0(
    "BP samples: ",
    nrow(blood)
  ),
  
  "",
  
  "SRS distribution:",
  
  paste0(
    names(
      table(
        blood$SRS
      )
    ),
    "=",
    as.integer(
      table(
        blood$SRS
      )
    ),
    collapse = "; "
  ),
  
  "",
  
  "CTS distribution:",
  
  paste0(
    names(
      table(
        blood$CTS
      )
    ),
    "=",
    as.integer(
      table(
        blood$CTS
      )
    ),
    collapse = "; "
  ),
  
  "",
  
  "Integrated distribution:",
  
  paste0(
    names(
      table(
        blood$integrated_group
      )
    ),
    "=",
    as.integer(
      table(
        blood$integrated_group
      )
    ),
    collapse = "; "
  ),
  
  "",
  
  "SRS score comparison:",
  
  paste0(
    paste(
      srs_medians$SRS,
      "n=",
      srs_medians$n,
      "median=",
      round(
        srs_medians$median,
        3
      )
    ),
    collapse = "; "
  ),
  
  paste0(
    "Wilcoxon p=",
    signif(
      srs_test$p.value,
      5
    )
  ),
  
  "",
  
  "CTS score comparison:",
  
  paste0(
    paste(
      cts_medians$CTS,
      "n=",
      cts_medians$n,
      "median=",
      round(
        cts_medians$median,
        3
      )
    ),
    collapse = "; "
  ),
  
  paste0(
    "Kruskal-Wallis p=",
    signif(
      cts_test$p.value,
      5
    ),
    "; epsilon^2=",
    round(
      cts_epsilon,
      3
    )
  ),
  
  "",
  
  "Integrated CTS/SRS score comparison:",
  
  paste0(
    paste(
      integrated_medians$integrated_group,
      "n=",
      integrated_medians$n,
      "median=",
      round(
        integrated_medians$median,
        3
      )
    ),
    collapse = "; "
  ),
  
  paste0(
    "Kruskal-Wallis p=",
    signif(
      integrated_test$p.value,
      5
    ),
    "; epsilon^2=",
    round(
      integrated_epsilon,
      3
    )
  ),
  
  "",
  
  "Five-gene score vs SRSq:",
  
  paste0(
    "Spearman rho=",
    round(
      rho_srsq,
      3
    ),
    "; p=",
    signif(
      srsq_cor$p.value,
      5
    )
  ),
  
  "",
  
  "INTERPRETATION:",
  
  paste0(
    "The five-gene host-response score captures a shared quantitative ",
    "axis across two independently assigned transcriptomic endotype ",
    "frameworks. SRS and CTS were not used for five-gene feature selection."
  )
)


results_summary_file <- file.path(
  text_dir,
  "145_results_summary.txt"
)


writeLines(
  results_summary_lines,
  results_summary_file
)


# =============================================================================
# 29. MASTER EXCEL WORKBOOK
# =============================================================================

master_xlsx <- file.path(
  table_dir,
  "145_endotype_five_gene_main_figure_data.xlsx"
)


openxlsx::write.xlsx(
  
  list(
    
    Run_info = data.frame(
      
      item = c(
        "script",
        "score_source",
        "expression_source",
        "primary_genes_up",
        "primary_genes_down",
        "n_BP"
      ),
      
      value = c(
        
        "145_build_endotype_five_gene_main_figure.R",
        
        score_file,
        
        counts_source,
        
        paste(
          primary_genes_up,
          collapse = ";"
        ),
        
        paste(
          primary_genes_down,
          collapse = ";"
        ),
        
        nrow(blood)
      )
    ),
    
    
    Final_BP_data =
      blood,
    
    
    Endotype_group_summary =
      group_summary,
    
    
    Figure_statistics =
      stats_summary,
    
    
    CTS_by_SRS =
      cross_tab,
    
    
    Five_gene_heatmap_z =
      heat_long,
    
    
    Five_gene_log2RPM =
      as.data.frame(
        panel_expression
      ) %>%
      
      tibble::rownames_to_column(
        "gene"
      )
  ),
  
  master_xlsx,
  
  overwrite = TRUE
)


# =============================================================================
# 30. FINAL CONSOLE OUTPUT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 145 v3 completed successfully.\n")
cat("====================================================================\n\n")


cat("INPUT SCORE FILE:\n")

print(
  normalizePath(
    score_file,
    winslash = "\\",
    mustWork = TRUE
  )
)


cat("\nENDOTYPE DISTRIBUTIONS:\n")


cat("SRS:\n")

print(
  table(
    blood$SRS
  )
)


cat("\nCTS:\n")

print(
  table(
    blood$CTS
  )
)


cat("\nCTS x SRS:\n")

print(
  table(
    blood$CTS,
    blood$SRS
  )
)


cat("\nFIVE-GENE SCORE BY SRS:\n")
print(srs_medians)


cat(
  "Wilcoxon p = ",
  signif(
    srs_test$p.value,
    6
  ),
  "\n",
  sep = ""
)


cat("\nFIVE-GENE SCORE BY CTS:\n")
print(cts_medians)


cat(
  "Kruskal-Wallis p = ",
  signif(
    cts_test$p.value,
    6
  ),
  "\n",
  sep = ""
)


cat(
  "epsilon^2 = ",
  round(
    cts_epsilon,
    4
  ),
  "\n",
  sep = ""
)


cat("\nINTEGRATED CTS/SRS CONTINUUM:\n")
print(integrated_medians)


cat(
  "Kruskal-Wallis p = ",
  signif(
    integrated_test$p.value,
    6
  ),
  "\n",
  sep = ""
)


cat(
  "epsilon^2 = ",
  round(
    integrated_epsilon,
    4
  ),
  "\n",
  sep = ""
)


cat("\nFIVE-GENE SCORE vs SRSq:\n")


cat(
  "Spearman rho = ",
  round(
    rho_srsq,
    4
  ),
  "\n",
  sep = ""
)


cat(
  "p = ",
  signif(
    srsq_cor$p.value,
    6
  ),
  "\n",
  sep = ""
)


cat("\nMAIN FIGURES:\n")


print(
  normalizePath(
    main_png,
    winslash = "\\",
    mustWork = FALSE
  )
)


print(
  normalizePath(
    main_pdf,
    winslash = "\\",
    mustWork = FALSE
  )
)


print(
  normalizePath(
    main_tiff,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nMASTER TABLE:\n")


print(
  normalizePath(
    master_xlsx,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nRESULTS SUMMARY:\n")


print(
  normalizePath(
    results_summary_file,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nDone.\n")