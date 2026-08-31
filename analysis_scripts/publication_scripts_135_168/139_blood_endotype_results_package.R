# ==============================================================================
# Script 139
# Blood endotype and biomarker publication results package
#
# Project: Sepsis_DESeq2
#
# PURPOSE
# Consolidate validated results from Scripts 135, 136b, 137 and 138 into
# a single publication-oriented evidence package.
#
# IMPORTANT
# - NO feature selection
# - NO new biomarker discovery
# - NO new predictive model
# - NO urine
# - NO lncRNA
#
# This script:
# 1. validates the final blood cohort and endotype assignments
# 2. summarizes the frozen five-gene score across SRS and CTS
# 3. consolidates clinical associations
# 4. consolidates age/sex sensitivity results
# 5. consolidates external signature benchmarking
# 6. generates manuscript-ready key-number tables
# 7. generates draft figure legends
# 8. generates EN/RU evidence summaries
#
# It does not replace Script 140, which will generate the manuscript Results.
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

script_name <- "139_blood_endotype_results_package.R"
run_date <- Sys.time()


cat("\n")
cat("====================================================================\n")
cat("Running Script 139\n")
cat("Blood endotype publication results package\n")
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
  "dplyr",
  "tidyr",
  "tibble",
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
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
})


cat("Required packages loaded successfully.\n\n")


# ==============================================================================
# 2. HELPERS
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
  
  
  H <- safe_htest_statistic(
    ht
  )
  
  
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


format_p <- function(x) {
  
  if (
    length(x) == 0 ||
    is.na(x)
  ) {
    
    return("NA")
  }
  
  
  if (x < 0.001) {
    
    return(
      format(
        x,
        scientific = TRUE,
        digits = 3
      )
    )
  }
  
  
  return(
    format(
      x,
      scientific = FALSE,
      digits = 3,
      trim = TRUE
    )
  )
}


format_num <- function(
    x,
    digits = 3
) {
  
  if (
    length(x) == 0 ||
    is.na(x)
  ) {
    
    return("NA")
  }
  
  
  return(
    format(
      round(
        x,
        digits
      ),
      nsmall = 0,
      trim = TRUE
    )
  )
}


# ==============================================================================
# 3. INPUT FILES
# ==============================================================================

scores135_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "135_validation",
  "tables",
  "135_blood_scores_with_final_endotypes.csv"
)


clinical_tests136b_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity",
  "tables",
  "136b_all_clinical_tests_updated_FDR.csv"
)


demographics136b_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity",
  "tables",
  "136b_blood_demographics.csv"
)


demographic_balance136b_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity",
  "tables",
  "136b_BP_BC_demographic_balance.csv"
)


adjusted_score136b_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity",
  "tables",
  "136b_primary_score_age_sex_adjusted.csv"
)


adjusted_genes136b_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity",
  "tables",
  "136b_primary_genes_age_sex_adjusted.csv"
)


auc137_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "137_benchmarking",
  "tables",
  "137_BP_BC_AUC_benchmark.csv"
)


correlations137_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "137_benchmarking",
  "tables",
  "137_signature_correlations_primary_SRSq.csv"
)


endotypes137_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "137_benchmarking",
  "tables",
  "137_signature_endotype_associations.csv"
)


gene_audit137_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "137_benchmarking",
  "tables",
  "137_benchmark_gene_direction_audit.csv"
)


signature_status137_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "137_benchmarking",
  "tables",
  "137_signature_gene_completeness.csv"
)


figure_manifest138_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "138_integrated_figures",
  "tables",
  "138_figure_manifest.csv"
)


palette138_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "138_integrated_figures",
  "tables",
  "138_publication_color_palette.csv"
)


required_files <- c(
  scores135_file,
  clinical_tests136b_file,
  demographics136b_file,
  demographic_balance136b_file,
  adjusted_score136b_file,
  adjusted_genes136b_file,
  auc137_file,
  correlations137_file,
  endotypes137_file,
  gene_audit137_file,
  signature_status137_file,
  figure_manifest138_file,
  palette138_file
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


cat("All required Script 135-138 outputs found.\n\n")


# ==============================================================================
# 4. OUTPUT DIRECTORIES
# ==============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "139_results_package"
)


tables_dir <- file.path(
  output_dir,
  "tables"
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
# 5. READ INPUT TABLES
# ==============================================================================

scores135 <- read.csv(
  scores135_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


clinical_tests136b <- read.csv(
  clinical_tests136b_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


demographics136b <- read.csv(
  demographics136b_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


demographic_balance136b <- read.csv(
  demographic_balance136b_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


adjusted_score136b <- read.csv(
  adjusted_score136b_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


adjusted_genes136b <- read.csv(
  adjusted_genes136b_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


auc137 <- read.csv(
  auc137_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


correlations137 <- read.csv(
  correlations137_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


endotypes137 <- read.csv(
  endotypes137_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


gene_audit137 <- read.csv(
  gene_audit137_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


signature_status137 <- read.csv(
  signature_status137_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


figure_manifest138 <- read.csv(
  figure_manifest138_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


palette138 <- read.csv(
  palette138_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ==============================================================================
# 6. VALIDATE FINAL BLOOD COHORT
# ==============================================================================

scores135$sample_id <- clean_sample_id(
  scores135$sample_id
)


scores135$condition <- dplyr::case_when(
  
  grepl(
    "^BP[0-9]+$",
    scores135$sample_id
  ) ~ "BP",
  
  grepl(
    "^BC[0-9]+$",
    scores135$sample_id
  ) ~ "BC",
  
  TRUE ~ NA_character_
)


scores135$condition <- factor(
  scores135$condition,
  levels = c(
    "BC",
    "BP"
  )
)


if (
  sum(
    scores135$condition == "BP"
  ) != 35
) {
  
  stop(
    "Final cohort validation failed: BP must equal 35."
  )
}


if (
  sum(
    scores135$condition == "BC"
  ) != 10
) {
  
  stop(
    "Final cohort validation failed: BC must equal 10."
  )
}


bp <- scores135 %>%
  filter(
    condition == "BP"
  )


if (nrow(bp) != 35) {
  
  stop(
    "BP data frame does not contain 35 patients."
  )
}


bp$SRS <- factor(
  bp$SRS,
  levels = c(
    "SRS1",
    "SRS2"
  )
)


bp$CTS <- factor(
  bp$CTS,
  levels = c(
    "CTS1",
    "CTS2",
    "CTS3"
  )
)


cat("FINAL BLOOD COHORT:\n")

print(
  table(
    scores135$condition
  )
)

cat("\n")


cat("FINAL SRS:\n")

print(
  table(
    bp$SRS
  )
)

cat("\n")


cat("FINAL CTS:\n")

print(
  table(
    bp$CTS
  )
)

cat("\n")


cat("FINAL CTS x SRS:\n")

print(
  table(
    bp$CTS,
    bp$SRS
  )
)

cat("\n")


# ==============================================================================
# 7. COHORT SUMMARY
# ==============================================================================

demographics136b$sample_id <- clean_sample_id(
  demographics136b$sample_id
)


cohort_summary <- demographics136b %>%
  
  mutate(
    condition =
      factor(
        condition,
        levels = c(
          "BC",
          "BP"
        )
      )
  ) %>%
  
  group_by(
    condition
  ) %>%
  
  summarise(
    
    n =
      n(),
    
    age_available =
      sum(
        !is.na(
          age_years
        )
      ),
    
    age_median =
      median(
        age_years,
        na.rm = TRUE
      ),
    
    age_q1 =
      quantile(
        age_years,
        0.25,
        na.rm = TRUE
      ),
    
    age_q3 =
      quantile(
        age_years,
        0.75,
        na.rm = TRUE
      ),
    
    age_min =
      min(
        age_years,
        na.rm = TRUE
      ),
    
    age_max =
      max(
        age_years,
        na.rm = TRUE
      ),
    
    female_n =
      sum(
        sex == "Female",
        na.rm = TRUE
      ),
    
    male_n =
      sum(
        sex == "Male",
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


# ==============================================================================
# 8. FIVE-GENE SCORE BY SRS
# ==============================================================================

srs_score_summary <- bp %>%
  
  group_by(
    SRS
  ) %>%
  
  summarise(
    
    n =
      n(),
    
    median =
      median(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    q1 =
      quantile(
        primary_5gene_score,
        0.25,
        na.rm = TRUE
      ),
    
    q3 =
      quantile(
        primary_5gene_score,
        0.75,
        na.rm = TRUE
      ),
    
    mean =
      mean(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    sd =
      stats::sd(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


srs_levels <- levels(
  droplevels(
    bp$SRS
  )
)


srs_values_1 <- bp$primary_5gene_score[
  bp$SRS == srs_levels[1]
]


srs_values_2 <- bp$primary_5gene_score[
  bp$SRS == srs_levels[2]
]


srs_score_test <- stats::wilcox.test(
  x = srs_values_1,
  y = srs_values_2,
  exact = FALSE,
  paired = FALSE
)


srs_test_summary <- tibble(
  
  comparison =
    paste0(
      srs_levels[1],
      " vs ",
      srs_levels[2]
    ),
  
  test =
    "Wilcoxon rank-sum",
  
  n =
    length(
      c(
        srs_values_1,
        srs_values_2
      )
    ),
  
  median_difference =
    median(
      srs_values_1,
      na.rm = TRUE
    ) -
    median(
      srs_values_2,
      na.rm = TRUE
    ),
  
  p_value =
    safe_htest_pvalue(
      srs_score_test
    )
)


# ==============================================================================
# 9. FIVE-GENE SCORE vs SRSq
# ==============================================================================

srsq_keep <- complete.cases(
  bp$primary_5gene_score,
  bp$SRSq
)


srsq_test <- suppressWarnings(
  stats::cor.test(
    x =
      bp$primary_5gene_score[
        srsq_keep
      ],
    y =
      bp$SRSq[
        srsq_keep
      ],
    method =
      "spearman",
    exact =
      FALSE
  )
)


srsq_summary <- tibble(
  
  comparison =
    "Primary five-gene score vs SRSq",
  
  n =
    sum(
      srsq_keep
    ),
  
  Spearman_rho =
    safe_htest_estimate(
      srsq_test
    ),
  
  p_value =
    safe_htest_pvalue(
      srsq_test
    )
)


# ==============================================================================
# 10. FIVE-GENE SCORE BY CTS
# ==============================================================================

cts_score_summary <- bp %>%
  
  group_by(
    CTS
  ) %>%
  
  summarise(
    
    n =
      n(),
    
    median =
      median(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    q1 =
      quantile(
        primary_5gene_score,
        0.25,
        na.rm = TRUE
      ),
    
    q3 =
      quantile(
        primary_5gene_score,
        0.75,
        na.rm = TRUE
      ),
    
    mean =
      mean(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    sd =
      stats::sd(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


cts_kw <- stats::kruskal.test(
  primary_5gene_score ~ CTS,
  data = bp
)


cts_epsilon2 <- epsilon_squared_kw(
  x =
    bp$primary_5gene_score,
  
  group =
    bp$CTS
)


cts_test_summary <- tibble(
  
  comparison =
    "Primary five-gene score across CTS1/CTS2/CTS3",
  
  test =
    "Kruskal-Wallis",
  
  n =
    nrow(
      bp
    ),
  
  epsilon_squared =
    cts_epsilon2,
  
  p_value =
    safe_htest_pvalue(
      cts_kw
    )
)


cts_pairwise_raw <- stats::pairwise.wilcox.test(
  
  x =
    bp$primary_5gene_score,
  
  g =
    bp$CTS,
  
  p.adjust.method =
    "BH",
  
  exact =
    FALSE
)


pairwise_matrix <- cts_pairwise_raw$p.value


cts_pairwise <- tibble()


if (!is.null(pairwise_matrix)) {
  
  pairwise_rows <- list()
  
  pair_counter <- 1L
  
  
  for (i in seq_len(
    nrow(
      pairwise_matrix
    )
  )) {
    
    for (j in seq_len(
      ncol(
        pairwise_matrix
      )
    )) {
      
      value <- pairwise_matrix[
        i,
        j
      ]
      
      
      if (!is.na(value)) {
        
        pairwise_rows[[pair_counter]] <- tibble(
          
          group_1 =
            rownames(
              pairwise_matrix
            )[i],
          
          group_2 =
            colnames(
              pairwise_matrix
            )[j],
          
          BH_adjusted_p =
            as.numeric(
              value
            )
        )
        
        
        pair_counter <- pair_counter + 1L
      }
    }
  }
  
  
  if (length(pairwise_rows) > 0) {
    
    cts_pairwise <- bind_rows(
      pairwise_rows
    )
  }
}


# ==============================================================================
# 11. INTEGRATED CTS/SRS STATE
# ==============================================================================

bp$CTS_SRS_state <- paste(
  bp$CTS,
  bp$SRS,
  sep = "/"
)


integrated_levels <- c(
  "CTS1/SRS1",
  "CTS2/SRS1",
  "CTS3/SRS1",
  "CTS3/SRS2"
)


bp$CTS_SRS_state <- factor(
  bp$CTS_SRS_state,
  levels =
    integrated_levels
)


integrated_state_summary <- bp %>%
  
  filter(
    !is.na(
      CTS_SRS_state
    )
  ) %>%
  
  group_by(
    CTS_SRS_state
  ) %>%
  
  summarise(
    
    n =
      n(),
    
    median =
      median(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    q1 =
      quantile(
        primary_5gene_score,
        0.25,
        na.rm = TRUE
      ),
    
    q3 =
      quantile(
        primary_5gene_score,
        0.75,
        na.rm = TRUE
      ),
    
    mean =
      mean(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    sd =
      stats::sd(
        primary_5gene_score,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


integrated_kw <- stats::kruskal.test(
  primary_5gene_score ~ CTS_SRS_state,
  data =
    bp %>%
    filter(
      !is.na(
        CTS_SRS_state
      )
    )
)


integrated_epsilon2 <- epsilon_squared_kw(
  x =
    bp$primary_5gene_score[
      !is.na(
        bp$CTS_SRS_state
      )
    ],
  
  group =
    bp$CTS_SRS_state[
      !is.na(
        bp$CTS_SRS_state
      )
    ]
)


integrated_test_summary <- tibble(
  
  comparison =
    paste(
      integrated_levels,
      collapse = " -> "
    ),
  
  test =
    "Kruskal-Wallis",
  
  n =
    sum(
      !is.na(
        bp$CTS_SRS_state
      )
    ),
  
  epsilon_squared =
    integrated_epsilon2,
  
  p_value =
    safe_htest_pvalue(
      integrated_kw
    )
)


# ==============================================================================
# 12. CLINICAL ASSOCIATIONS
# ==============================================================================

clinical_ranked <- clinical_tests136b %>%
  
  arrange(
    p_value
  )


clinical_fdr <- clinical_tests136b %>%
  
  filter(
    !is.na(
      BH_global
    ),
    BH_global < 0.05
  ) %>%
  
  arrange(
    BH_global
  )


clinical_nominal <- clinical_tests136b %>%
  
  filter(
    !is.na(
      p_value
    ),
    p_value < 0.05
  ) %>%
  
  arrange(
    p_value
  )


clinical_summary_counts <- tibble(
  
  metric = c(
    "Total clinical association tests",
    "Nominal p < 0.05",
    "Global BH FDR < 0.05"
  ),
  
  n = c(
    nrow(
      clinical_tests136b
    ),
    
    nrow(
      clinical_nominal
    ),
    
    nrow(
      clinical_fdr
    )
  )
)


# ==============================================================================
# 13. DEMOGRAPHIC SENSITIVITY
# ==============================================================================

primary_adjusted_result <- adjusted_score136b %>%
  
  slice_head(
    n = 1
  )


primary_gene_adjusted <- adjusted_genes136b %>%
  
  arrange(
    p_value
  )


# ==============================================================================
# 14. EXTERNAL SIGNATURE BENCHMARKING
# ==============================================================================

external_auc <- auc137 %>%
  
  arrange(
    desc(
      AUC_fixed_direction
    )
  )


signature_primary_correlations <- correlations137 %>%
  
  filter(
    target ==
      "primary_5gene_score"
  ) %>%
  
  arrange(
    desc(
      Spearman_rho
    )
  )


signature_srsq_correlations <- correlations137 %>%
  
  filter(
    target ==
      "SRSq"
  ) %>%
  
  arrange(
    desc(
      Spearman_rho
    )
  )


signature_cts_associations <- endotypes137 %>%
  
  filter(
    framework ==
      "CTS"
  ) %>%
  
  arrange(
    desc(
      effect
    )
  )


signature_srs_associations <- endotypes137 %>%
  
  filter(
    framework ==
      "SRS"
  ) %>%
  
  arrange(
    p_value
  )


# ==============================================================================
# 15. GENE-DIRECTION CONCORDANCE
# ==============================================================================

gene_direction_summary <- gene_audit137 %>%
  
  summarise(
    
    genes_evaluated =
      sum(
        present,
        na.rm = TRUE
      ),
    
    direction_concordant =
      sum(
        direction_concordant == TRUE,
        na.rm = TRUE
      ),
    
    direction_discordant =
      sum(
        direction_concordant == FALSE,
        na.rm = TRUE
      )
  )


gene_direction_discordant <- gene_audit137 %>%
  
  filter(
    present,
    direction_concordant == FALSE
  )


# ==============================================================================
# 16. KEY RESULTS TABLE
# ==============================================================================

primary_crp <- clinical_tests136b %>%
  
  filter(
    framework ==
      "Primary_5gene_score",
    clinical_label ==
      "CRP",
    test ==
      "Spearman"
  ) %>%
  
  slice_head(
    n = 1
  )


srsq_crp <- clinical_tests136b %>%
  
  filter(
    framework ==
      "SRSq",
    clinical_label ==
      "CRP",
    test ==
      "Spearman"
  ) %>%
  
  slice_head(
    n = 1
  )


lifts_primary <- correlations137 %>%
  
  filter(
    display_name ==
      "LIFTS-like",
    target ==
      "primary_5gene_score"
  ) %>%
  
  slice_head(
    n = 1
  )


lifts_srsq <- correlations137 %>%
  
  filter(
    display_name ==
      "LIFTS-like",
    target ==
      "SRSq"
  ) %>%
  
  slice_head(
    n = 1
  )


lifts_cts <- endotypes137 %>%
  
  filter(
    display_name ==
      "LIFTS-like",
    framework ==
      "CTS"
  ) %>%
  
  slice_head(
    n = 1
  )


primary_auc <- auc137 %>%
  
  filter(
    display_name ==
      "Primary 5-gene"
  ) %>%
  
  slice_head(
    n = 1
  )


key_results <- bind_rows(
  
  tibble(
    section = "Cohort",
    result = "Blood RNA-seq cohort",
    value = "BP n=35; BC n=10"
  ),
  
  tibble(
    section = "Endotypes",
    result = "SRS distribution",
    value = "SRS1 n=28; SRS2 n=7"
  ),
  
  tibble(
    section = "Endotypes",
    result = "CTS distribution",
    value = "CTS1 n=14; CTS2 n=6; CTS3 n=15"
  ),
  
  tibble(
    section = "Five-gene score",
    result = "SRS1 vs SRS2",
    value = paste0(
      "p=",
      format_p(
        srs_test_summary$p_value[1]
      )
    )
  ),
  
  tibble(
    section = "Five-gene score",
    result = "Score vs SRSq",
    value = paste0(
      "rho=",
      format_num(
        srsq_summary$Spearman_rho[1]
      ),
      "; p=",
      format_p(
        srsq_summary$p_value[1]
      )
    )
  ),
  
  tibble(
    section = "Five-gene score",
    result = "Score across CTS",
    value = paste0(
      "epsilon^2=",
      format_num(
        cts_test_summary$epsilon_squared[1]
      ),
      "; p=",
      format_p(
        cts_test_summary$p_value[1]
      )
    )
  ),
  
  tibble(
    section = "Five-gene score",
    result = "Integrated CTS/SRS states",
    value = paste0(
      "epsilon^2=",
      format_num(
        integrated_test_summary$epsilon_squared[1]
      ),
      "; p=",
      format_p(
        integrated_test_summary$p_value[1]
      )
    )
  )
)


if (nrow(primary_crp) == 1) {
  
  key_results <- bind_rows(
    
    key_results,
    
    tibble(
      section = "Clinical",
      result = "Five-gene score vs CRP",
      value = paste0(
        "rho=",
        format_num(
          primary_crp$effect[1]
        ),
        "; p=",
        format_p(
          primary_crp$p_value[1]
        ),
        "; BH=",
        format_p(
          primary_crp$BH_global[1]
        )
      )
    )
  )
}


if (nrow(srsq_crp) == 1) {
  
  key_results <- bind_rows(
    
    key_results,
    
    tibble(
      section = "Clinical",
      result = "SRSq vs CRP",
      value = paste0(
        "rho=",
        format_num(
          srsq_crp$effect[1]
        ),
        "; p=",
        format_p(
          srsq_crp$p_value[1]
        ),
        "; BH=",
        format_p(
          srsq_crp$BH_global[1]
        )
      )
    )
  )
}


if (nrow(primary_adjusted_result) == 1) {
  
  key_results <- bind_rows(
    
    key_results,
    
    tibble(
      section = "Demographic sensitivity",
      result = "Age/sex-adjusted five-gene BP vs BC effect",
      value = paste0(
        "effect=",
        format_num(
          primary_adjusted_result$
            condition_effect_BP_vs_BC[1]
        ),
        "; p=",
        format_p(
          primary_adjusted_result$p_value[1]
        )
      )
    )
  )
}


if (nrow(primary_auc) == 1) {
  
  key_results <- bind_rows(
    
    key_results,
    
    tibble(
      section = "Contextual discrimination",
      result = "Primary five-gene apparent BP vs BC AUC",
      value = paste0(
        "AUC=",
        format_num(
          primary_auc$AUC_fixed_direction[1]
        ),
        "; CI ",
        format_num(
          primary_auc$AUC_CI_low[1]
        ),
        "-",
        format_num(
          primary_auc$AUC_CI_high[1]
        )
      )
    )
  )
}


if (nrow(lifts_primary) == 1) {
  
  key_results <- bind_rows(
    
    key_results,
    
    tibble(
      section = "External benchmarking",
      result = "LIFTS-like vs primary five-gene",
      value = paste0(
        "rho=",
        format_num(
          lifts_primary$Spearman_rho[1]
        ),
        "; p=",
        format_p(
          lifts_primary$p_value[1]
        )
      )
    )
  )
}


if (nrow(lifts_srsq) == 1) {
  
  key_results <- bind_rows(
    
    key_results,
    
    tibble(
      section = "External benchmarking",
      result = "LIFTS-like vs SRSq",
      value = paste0(
        "rho=",
        format_num(
          lifts_srsq$Spearman_rho[1]
        ),
        "; p=",
        format_p(
          lifts_srsq$p_value[1]
        )
      )
    )
  )
}


if (nrow(lifts_cts) == 1) {
  
  key_results <- bind_rows(
    
    key_results,
    
    tibble(
      section = "External benchmarking",
      result = "LIFTS-like across CTS",
      value = paste0(
        "epsilon^2=",
        format_num(
          lifts_cts$effect[1]
        ),
        "; p=",
        format_p(
          lifts_cts$p_value[1]
        ),
        "; BH=",
        format_p(
          lifts_cts$BH_endotype[1]
        )
      )
    )
  )
}


# ==============================================================================
# 17. FIGURE LEGENDS
# ==============================================================================

figure_legends_en <- c(
  
  "FIGURE 1",
  paste0(
    "Figure 1. Blood five-gene host-response signature in sepsis. ",
    "(A) Distribution of the frozen five-gene myeloid-adaptive balance ",
    "score in blood from patients with sepsis (BP) and healthy controls (BC). ",
    "(B) TMM-normalized logCPM expression of CD177, HK3, IRAK3, CARD11, ",
    "and IKZF2. (C) Sepsis-associated expression differences for the five ",
    "component genes after simultaneous adjustment for age and sex. ",
    "Positive effects indicate higher expression in sepsis and negative ",
    "effects indicate lower expression in sepsis."
  ),
  
  "",
  
  "FIGURE 2",
  paste0(
    "Figure 2. Hierarchical organization of blood transcriptomic endotypes. ",
    "(A) Five-gene score according to SRS class. ",
    "(B) Five-gene score across Consensus Transcriptomic Subtypes (CTS). ",
    "(C) Five-gene score across integrated CTS/SRS states, demonstrating ",
    "the gradient CTS1/SRS1 to CTS2/SRS1 to CTS3/SRS1 to CTS3/SRS2. ",
    "(D) Relationship between the five-gene score and continuous SRSq; ",
    "point color indicates CTS and point shape indicates SRS class."
  ),
  
  "",
  
  "FIGURE 3",
  paste0(
    "Figure 3. Five-gene expression architecture across sepsis endotypes. ",
    "Gene-wise standardized expression of CD177, HK3, IRAK3, CARD11, and ",
    "IKZF2 across 35 patients with sepsis ordered by CTS, SRS, and SRSq. ",
    "Annotation bars indicate final SRS and CTS assignments. Red indicates ",
    "higher and blue indicates lower expression relative to the gene-wise ",
    "mean within the sepsis cohort."
  ),
  
  "",
  
  "FIGURE 4",
  paste0(
    "Figure 4. Clinical context of the blood host-response axis. ",
    "(A) Association between the five-gene score and serum C-reactive ",
    "protein (CRP); colors denote CTS. ",
    "(B) Spearman correlations of the five-gene score and SRSq with ",
    "continuous clinical variables. Statistical support is indicated by ",
    "global FDR significance, nominal significance, or non-significance."
  ),
  
  "",
  
  "FIGURE 5",
  paste0(
    "Figure 5. Blood transcriptomic signatures converge on a common ",
    "host-response continuum. ",
    "(A) Spearman correlations between predefined published transcriptomic ",
    "signatures and the primary five-gene score. ",
    "(B) Corresponding correlations with continuous SRSq. ",
    "(C) CTS-associated effect sizes for the evaluated signatures, expressed ",
    "as Kruskal-Wallis epsilon-squared. ",
    "(D) Median standardized signature scores with interquartile ranges ",
    "across CTS1, CTS2, and CTS3."
  )
)


writeLines(
  figure_legends_en,
  file.path(
    text_dir,
    "139_figure_legends_EN.txt"
  )
)


# ==============================================================================
# 18. ENGLISH RESULTS EVIDENCE SUMMARY
# ==============================================================================

results_summary_en <- c(
  
  "SCRIPT 139 - BLOOD ENDOTYPE RESULTS EVIDENCE PACKAGE",
  
  "====================================================================",
  
  "",
  
  "COHORT",
  
  "The blood RNA-seq dataset comprised 35 patients with sepsis and 10 healthy controls.",
  
  "",
  
  "ENDOTYPE STRUCTURE",
  
  paste0(
    "Among patients with sepsis, final SRS classification identified ",
    "28 SRS1 and 7 SRS2 samples."
  ),
  
  paste0(
    "Final CTS classification identified 14 CTS1, 6 CTS2, and 15 CTS3 samples."
  ),
  
  paste0(
    "All CTS1 and CTS2 samples were SRS1, whereas CTS3 included ",
    "8 SRS1 and 7 SRS2 samples."
  ),
  
  "",
  
  "FIVE-GENE HOST-RESPONSE AXIS",
  
  paste0(
    "The frozen five-gene score differed between SRS classes ",
    "(Wilcoxon p=",
    format_p(
      srs_test_summary$p_value[1]
    ),
    ")."
  ),
  
  paste0(
    "The score correlated strongly with continuous SRSq ",
    "(Spearman rho=",
    format_num(
      srsq_summary$Spearman_rho[1]
    ),
    ", p=",
    format_p(
      srsq_summary$p_value[1]
    ),
    ")."
  ),
  
  paste0(
    "The score differed strongly across CTS classes ",
    "(Kruskal-Wallis epsilon^2=",
    format_num(
      cts_test_summary$epsilon_squared[1]
    ),
    ", p=",
    format_p(
      cts_test_summary$p_value[1]
    ),
    ")."
  ),
  
  paste0(
    "Integrated CTS/SRS states showed an even clearer hierarchical gradient ",
    "(epsilon^2=",
    format_num(
      integrated_test_summary$epsilon_squared[1]
    ),
    ", p=",
    format_p(
      integrated_test_summary$p_value[1]
    ),
    ")."
  ),
  
  "",
  
  "CLINICAL CONTEXT",
  
  paste0(
    "After global multiple-testing correction, ",
    nrow(
      clinical_fdr
    ),
    " of ",
    nrow(
      clinical_tests136b
    ),
    " clinical association tests remained significant."
  )
)


if (nrow(primary_crp) == 1) {
  
  results_summary_en <- c(
    
    results_summary_en,
    
    paste0(
      "The five-gene score was associated with CRP ",
      "(rho=",
      format_num(
        primary_crp$effect[1]
      ),
      ", p=",
      format_p(
        primary_crp$p_value[1]
      ),
      ", global BH=",
      format_p(
        primary_crp$BH_global[1]
      ),
      ")."
    )
  )
}


if (nrow(srsq_crp) == 1) {
  
  results_summary_en <- c(
    
    results_summary_en,
    
    paste0(
      "SRSq was also associated with CRP ",
      "(rho=",
      format_num(
        srsq_crp$effect[1]
      ),
      ", p=",
      format_p(
        srsq_crp$p_value[1]
      ),
      ", global BH=",
      format_p(
        srsq_crp$BH_global[1]
      ),
      ")."
    )
  )
}


results_summary_en <- c(
  
  results_summary_en,
  
  "",
  
  "DEMOGRAPHIC SENSITIVITY",
  
  paste0(
    "Although BP and BC differed in age and sex distribution, ",
    "neither age nor sex was significantly associated with the five-gene ",
    "score, SRSq, SRS class, or CTS class within the sepsis cohort."
  )
)


if (nrow(primary_adjusted_result) == 1) {
  
  results_summary_en <- c(
    
    results_summary_en,
    
    paste0(
      "The five-gene score remained strongly associated with sepsis after ",
      "simultaneous adjustment for age and sex ",
      "(adjusted BP-BC effect=",
      format_num(
        primary_adjusted_result$
          condition_effect_BP_vs_BC[1]
      ),
      ", p=",
      format_p(
        primary_adjusted_result$p_value[1]
      ),
      ")."
    )
  )
}


results_summary_en <- c(
  
  results_summary_en,
  
  "",
  
  "EXTERNAL TRANSCRIPTOMIC BENCHMARKING",
  
  paste0(
    "All genes required for the seven evaluated current-study and published ",
    "signature implementations were detected in the RNA-seq dataset."
  ),
  
  paste0(
    gene_direction_summary$direction_concordant[1],
    " of ",
    gene_direction_summary$genes_evaluated[1],
    " predefined biomarker genes showed the expected direction of change."
  )
)


if (nrow(lifts_primary) == 1) {
  
  results_summary_en <- c(
    
    results_summary_en,
    
    paste0(
      "LIFTS-like showed the strongest external concordance with the primary ",
      "five-gene score (rho=",
      format_num(
        lifts_primary$Spearman_rho[1]
      ),
      ", p=",
      format_p(
        lifts_primary$p_value[1]
      ),
      ")."
    )
  )
}


if (nrow(lifts_srsq) == 1) {
  
  results_summary_en <- c(
    
    results_summary_en,
    
    paste0(
      "LIFTS-like was also strongly associated with SRSq ",
      "(rho=",
      format_num(
        lifts_srsq$Spearman_rho[1]
      ),
      ", p=",
      format_p(
        lifts_srsq$p_value[1]
      ),
      ")."
    )
  )
}


if (nrow(lifts_cts) == 1) {
  
  results_summary_en <- c(
    
    results_summary_en,
    
    paste0(
      "Its CTS effect size was epsilon^2=",
      format_num(
        lifts_cts$effect[1]
      ),
      " (p=",
      format_p(
        lifts_cts$p_value[1]
      ),
      ")."
    )
  )
}


results_summary_en <- c(
  
  results_summary_en,
  
  "",
  
  "INTERPRETATION",
  
  paste0(
    "Together, the results support a hierarchical blood host-response ",
    "organization in which SRS describes a dominant immune-dysfunction axis, ",
    "CTS resolves biologically distinct positions within that axis, and the ",
    "five-gene score provides a compact molecular representation of the ",
    "resulting myeloid-adaptive continuum."
  ),
  
  paste0(
    "Concordance with independently developed published transcriptomic ",
    "signatures supports the biological reproducibility of this continuum."
  ),
  
  "",
  
  "IMPORTANT LIMITATIONS",
  
  paste0(
    "The primary five-gene panel was selected in the current cohort and ",
    "therefore does not constitute independent diagnostic validation."
  ),
  
  paste0(
    "The control group consists of healthy volunteers rather than ",
    "noninfectious critically ill or SIRS controls."
  ),
  
  paste0(
    "Contextual BP-versus-BC AUC values should therefore not be interpreted ",
    "as comparative clinical diagnostic performance."
  )
)


writeLines(
  results_summary_en,
  file.path(
    text_dir,
    "139_results_evidence_summary_EN.txt"
  )
)


# ==============================================================================
# 19. RUSSIAN RESULTS EVIDENCE SUMMARY
# ==============================================================================

results_summary_ru <- c(
  
  "SCRIPT 139 - ИТОГОВЫЙ ПАКЕТ РЕЗУЛЬТАТОВ ПО ЭНДОТИПАМ КРОВИ",
  
  "====================================================================",
  
  "",
  
  "КОГОРТА",
  
  "В анализ крови включены 35 пациентов с сепсисом и 10 здоровых контролей.",
  
  "",
  
  "ЭНДОТИПИЧЕСКАЯ СТРУКТУРА",
  
  "Финальное распределение SRS: SRS1 - 28, SRS2 - 7.",
  
  "Финальное распределение CTS: CTS1 - 14, CTS2 - 6, CTS3 - 15.",
  
  paste0(
    "Все CTS1 и CTS2 относились к SRS1; CTS3 включал 8 SRS1 и 7 SRS2."
  ),
  
  "",
  
  "ПЯТИГЕННАЯ HOST-RESPONSE ОСЬ",
  
  paste0(
    "Пятигенный score различался между SRS1 и SRS2: p=",
    format_p(
      srs_test_summary$p_value[1]
    ),
    "."
  ),
  
  paste0(
    "Связь score с SRSq: rho=",
    format_num(
      srsq_summary$Spearman_rho[1]
    ),
    "; p=",
    format_p(
      srsq_summary$p_value[1]
    ),
    "."
  ),
  
  paste0(
    "Различия между CTS: epsilon^2=",
    format_num(
      cts_test_summary$epsilon_squared[1]
    ),
    "; p=",
    format_p(
      cts_test_summary$p_value[1]
    ),
    "."
  ),
  
  paste0(
    "Для интегрированных CTS/SRS-состояний: epsilon^2=",
    format_num(
      integrated_test_summary$epsilon_squared[1]
    ),
    "; p=",
    format_p(
      integrated_test_summary$p_value[1]
    ),
    "."
  ),
  
  "",
  
  "КЛИНИЧЕСКИЕ АССОЦИАЦИИ",
  
  paste0(
    "Всего оценено ",
    nrow(
      clinical_tests136b
    ),
    " клинических ассоциаций; после глобальной BH-коррекции ",
    "значимыми остались ",
    nrow(
      clinical_fdr
    ),
    "."
  )
)


if (nrow(primary_crp) == 1) {
  
  results_summary_ru <- c(
    
    results_summary_ru,
    
    paste0(
      "Пятигенный score коррелировал с CRP: rho=",
      format_num(
        primary_crp$effect[1]
      ),
      "; p=",
      format_p(
        primary_crp$p_value[1]
      ),
      "; global BH=",
      format_p(
        primary_crp$BH_global[1]
      ),
      "."
    )
  )
}


if (nrow(srsq_crp) == 1) {
  
  results_summary_ru <- c(
    
    results_summary_ru,
    
    paste0(
      "SRSq также коррелировал с CRP: rho=",
      format_num(
        srsq_crp$effect[1]
      ),
      "; p=",
      format_p(
        srsq_crp$p_value[1]
      ),
      "; global BH=",
      format_p(
        srsq_crp$BH_global[1]
      ),
      "."
    )
  )
}


results_summary_ru <- c(
  
  results_summary_ru,
  
  "",
  
  "ВОЗРАСТ И ПОЛ",
  
  paste0(
    "Несмотря на различия BP и BC по возрасту и полу, внутри группы сепсиса ",
    "ни возраст, ни пол не были связаны с пятигенным score, SRSq, SRS или CTS."
  )
)


if (nrow(primary_adjusted_result) == 1) {
  
  results_summary_ru <- c(
    
    results_summary_ru,
    
    paste0(
      "После одновременной поправки на возраст и пол различие пятигенного ",
      "score между BP и BC сохранялось: adjusted effect=",
      format_num(
        primary_adjusted_result$
          condition_effect_BP_vs_BC[1]
      ),
      "; p=",
      format_p(
        primary_adjusted_result$p_value[1]
      ),
      "."
    )
  )
}


results_summary_ru <- c(
  
  results_summary_ru,
  
  "",
  
  "СОПОСТАВЛЕНИЕ С ОПУБЛИКОВАННЫМИ SIGNATURES",
  
  paste0(
    "Все необходимые гены опубликованных signatures присутствовали ",
    "в RNA-seq матрице."
  ),
  
  paste0(
    gene_direction_summary$direction_concordant[1],
    " из ",
    gene_direction_summary$genes_evaluated[1],
    " исследованных biomarker genes имели ожидаемое направление изменения."
  )
)


if (nrow(lifts_primary) == 1) {
  
  results_summary_ru <- c(
    
    results_summary_ru,
    
    paste0(
      "Наиболее сильную внешнюю согласованность с нашей пятигенной панелью ",
      "показал LIFTS-like score: rho=",
      format_num(
        lifts_primary$Spearman_rho[1]
      ),
      "; p=",
      format_p(
        lifts_primary$p_value[1]
      ),
      "."
    )
  )
}


if (nrow(lifts_srsq) == 1) {
  
  results_summary_ru <- c(
    
    results_summary_ru,
    
    paste0(
      "Связь LIFTS-like с SRSq: rho=",
      format_num(
        lifts_srsq$Spearman_rho[1]
      ),
      "; p=",
      format_p(
        lifts_srsq$p_value[1]
      ),
      "."
    )
  )
}


results_summary_ru <- c(
  
  results_summary_ru,
  
  "",
  
  "ОБЩАЯ ИНТЕРПРЕТАЦИЯ",
  
  paste0(
    "Полученные данные поддерживают иерархическую организацию системного ",
    "ответа при сепсисе: SRS отражает доминирующую ось иммунной дисфункции, ",
    "CTS разделяет ее на биологически различимые состояния, а пятигенный ",
    "score представляет эту структуру в компактной форме."
  ),
  
  paste0(
    "Согласованность с независимо разработанными опубликованными ",
    "транскриптомными signatures указывает на воспроизводимость общей ",
    "myeloid-adaptive host-response оси."
  ),
  
  "",
  
  "ОГРАНИЧЕНИЯ",
  
  paste0(
    "Пятигенная панель была выбрана в текущей когорте и поэтому пока ",
    "не является независимо валидированным диагностическим тестом."
  ),
  
  paste0(
    "Контрольная группа представлена здоровыми добровольцами, а не ",
    "пациентами с неинфекционным системным воспалением или SIRS."
  )
)


writeLines(
  results_summary_ru,
  file.path(
    text_dir,
    "139_results_evidence_summary_RU.txt"
  )
)


# ==============================================================================
# 20. SAVE CSV TABLES
# ==============================================================================

write.csv(
  key_results,
  file.path(
    tables_dir,
    "139_key_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  cohort_summary,
  file.path(
    tables_dir,
    "139_cohort_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  srs_score_summary,
  file.path(
    tables_dir,
    "139_primary_score_by_SRS.csv"
  ),
  row.names = FALSE
)


write.csv(
  cts_score_summary,
  file.path(
    tables_dir,
    "139_primary_score_by_CTS.csv"
  ),
  row.names = FALSE
)


write.csv(
  integrated_state_summary,
  file.path(
    tables_dir,
    "139_primary_score_by_integrated_CTS_SRS.csv"
  ),
  row.names = FALSE
)


write.csv(
  clinical_fdr,
  file.path(
    tables_dir,
    "139_clinical_associations_global_FDR.csv"
  ),
  row.names = FALSE
)


write.csv(
  external_auc,
  file.path(
    tables_dir,
    "139_external_signature_AUC_contextual.csv"
  ),
  row.names = FALSE
)


write.csv(
  signature_primary_correlations,
  file.path(
    tables_dir,
    "139_signature_correlations_primary.csv"
  ),
  row.names = FALSE
)


write.csv(
  signature_srsq_correlations,
  file.path(
    tables_dir,
    "139_signature_correlations_SRSq.csv"
  ),
  row.names = FALSE
)


write.csv(
  signature_cts_associations,
  file.path(
    tables_dir,
    "139_signature_CTS_effects.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 21. PUBLICATION WORKBOOK
# ==============================================================================

wb <- openxlsx::createWorkbook()


run_info <- tibble(
  
  parameter = c(
    "script",
    "run_date",
    "study_scope",
    "BP_n",
    "BC_n",
    "SRS1_n",
    "SRS2_n",
    "CTS1_n",
    "CTS2_n",
    "CTS3_n",
    "primary_panel",
    "analysis_status",
    "interpretation_limit"
  ),
  
  value = c(
    script_name,
    as.character(run_date),
    "Blood RNA-seq only",
    "35",
    "10",
    "28",
    "7",
    "14",
    "6",
    "15",
    "CD177 + HK3 + IRAK3 - CARD11 - IKZF2",
    "Publication evidence consolidation; no new feature selection",
    paste0(
      "Internal candidate signature; healthy controls; ",
      "no independent diagnostic validation"
    )
  )
)


sheet_list <- list(
  
  "00_run_info" =
    run_info,
  
  "01_key_results" =
    key_results,
  
  "02_cohort" =
    cohort_summary,
  
  "03_demographic_balance" =
    demographic_balance136b,
  
  "04_score_by_SRS" =
    srs_score_summary,
  
  "05_SRS_tests" =
    bind_rows(
      srs_test_summary,
      srsq_summary
    ),
  
  "06_score_by_CTS" =
    cts_score_summary,
  
  "07_CTS_tests" =
    bind_rows(
      cts_test_summary,
      cts_pairwise
    ),
  
  "08_integrated_CTS_SRS" =
    integrated_state_summary,
  
  "09_integrated_test" =
    integrated_test_summary,
  
  "10_clinical_FDR" =
    clinical_fdr,
  
  "11_clinical_ranked" =
    clinical_ranked,
  
  "12_age_sex_score" =
    adjusted_score136b,
  
  "13_age_sex_genes" =
    primary_gene_adjusted,
  
  "14_signature_status" =
    signature_status137,
  
  "15_signature_AUC" =
    external_auc,
  
  "16_signature_vs_primary" =
    signature_primary_correlations,
  
  "17_signature_vs_SRSq" =
    signature_srsq_correlations,
  
  "18_signature_CTS" =
    signature_cts_associations,
  
  "19_signature_SRS" =
    signature_srs_associations,
  
  "20_gene_direction_audit" =
    gene_audit137,
  
  "21_figure_manifest" =
    figure_manifest138,
  
  "22_color_palette" =
    palette138
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)


for (sheet_name in names(sheet_list)) {
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  
  openxlsx::writeData(
    wb,
    sheet_name,
    sheet_list[[sheet_name]],
    headerStyle = header_style
  )
  
  
  openxlsx::freezePane(
    wb,
    sheet_name,
    firstRow = TRUE
  )
  
  
  openxlsx::setColWidths(
    wb,
    sheet_name,
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


workbook_file <- file.path(
  tables_dir,
  "139_blood_endotype_results_package.xlsx"
)


openxlsx::saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)


# ==============================================================================
# 22. INPUT MANIFEST
# ==============================================================================

input_manifest <- tibble(
  
  input = c(
    "Script135_scores",
    "Script136b_clinical_tests",
    "Script136b_demographics",
    "Script136b_demographic_balance",
    "Script136b_adjusted_score",
    "Script136b_adjusted_genes",
    "Script137_AUC",
    "Script137_correlations",
    "Script137_endotypes",
    "Script137_gene_audit",
    "Script137_signature_status",
    "Script138_figure_manifest",
    "Script138_palette"
  ),
  
  path = required_files
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
    "139_input_file_manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 23. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "139_sessionInfo.txt"
  )
)


# ==============================================================================
# 24. FINAL REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 139 completed successfully.\n")
cat("====================================================================\n\n")


cat("FINAL COHORT:\n")
cat("- BP n=35\n")
cat("- BC n=10\n")
cat("- SRS1 n=28\n")
cat("- SRS2 n=7\n")
cat("- CTS1 n=14\n")
cat("- CTS2 n=6\n")
cat("- CTS3 n=15\n\n")


cat("Primary five-gene score:\n")

cat(
  "- SRS1 vs SRS2: p=",
  format_p(
    srs_test_summary$p_value[1]
  ),
  "\n",
  sep = ""
)


cat(
  "- Score vs SRSq: rho=",
  format_num(
    srsq_summary$Spearman_rho[1]
  ),
  "; p=",
  format_p(
    srsq_summary$p_value[1]
  ),
  "\n",
  sep = ""
)


cat(
  "- CTS effect: epsilon^2=",
  format_num(
    cts_test_summary$epsilon_squared[1]
  ),
  "; p=",
  format_p(
    cts_test_summary$p_value[1]
  ),
  "\n",
  sep = ""
)


cat(
  "- Integrated CTS/SRS effect: epsilon^2=",
  format_num(
    integrated_test_summary$epsilon_squared[1]
  ),
  "; p=",
  format_p(
    integrated_test_summary$p_value[1]
  ),
  "\n\n",
  sep = ""
)


cat("Clinical associations:\n")

cat(
  "- Total tests: ",
  nrow(
    clinical_tests136b
  ),
  "\n",
  sep = ""
)

cat(
  "- Nominal p<0.05: ",
  nrow(
    clinical_nominal
  ),
  "\n",
  sep = ""
)

cat(
  "- Global BH FDR<0.05: ",
  nrow(
    clinical_fdr
  ),
  "\n\n",
  sep = ""
)


cat("External signature benchmarking:\n")

cat(
  "- Signatures complete: ",
  sum(
    signature_status137$complete,
    na.rm = TRUE
  ),
  "/",
  nrow(
    signature_status137
  ),
  "\n",
  sep = ""
)


cat(
  "- Expected gene direction: ",
  gene_direction_summary$direction_concordant[1],
  "/",
  gene_direction_summary$genes_evaluated[1],
  "\n\n",
  sep = ""
)


if (nrow(lifts_primary) == 1) {
  
  cat(
    "- LIFTS-like vs primary score: rho=",
    format_num(
      lifts_primary$Spearman_rho[1]
    ),
    "; p=",
    format_p(
      lifts_primary$p_value[1]
    ),
    "\n",
    sep = ""
  )
}


if (nrow(lifts_srsq) == 1) {
  
  cat(
    "- LIFTS-like vs SRSq: rho=",
    format_num(
      lifts_srsq$Spearman_rho[1]
    ),
    "; p=",
    format_p(
      lifts_srsq$p_value[1]
    ),
    "\n\n",
    sep = ""
  )
}


cat("Main output workbook:\n")

cat(
  normalizePath(
    workbook_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat("Text outputs:\n")

cat(
  "1) ",
  normalizePath(
    file.path(
      text_dir,
      "139_results_evidence_summary_EN.txt"
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
      text_dir,
      "139_results_evidence_summary_RU.txt"
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
      text_dir,
      "139_figure_legends_EN.txt"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat("IMPORTANT:\n")
cat("- No new feature selection.\n")
cat("- No new prediction model.\n")
cat("- Primary five-gene panel remains frozen.\n")
cat("- Final CTS source remains BP-only: 14/6/15.\n")
cat("- AUC remains contextual and supplementary.\n")
cat("- No urine.\n")
cat("- No lncRNA.\n\n")


cat("Next planned script:\n")
cat("140_blood_endotype_manuscript_results.R\n\n")

cat("Done.\n")