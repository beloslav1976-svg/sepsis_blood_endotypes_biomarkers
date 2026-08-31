################################################################################
# Script 150
# FINAL VISUAL-POLISH VERSION
#
# Main Figure 5
# External association of the frozen five-gene host-response score
# with organ-dysfunction severity
#
# Dataset: GSE185263
# Project: Sepsis_DESeq2
#
#
# FIGURE STRUCTURE
# ----------------
# A. Primary five-gene score vs continuous 24-h SOFA
# B. Secondary SOFA >=2 vs SOFA 0-1 comparison
# C. Individual frozen component genes vs continuous SOFA
# D. Age-, sex-, and location-adjusted SOFA coefficient
# E. Geographic sensitivity within GSE185263
#
#
# STATUS
# ------
# FINAL publication-packaging / visual-polish script.
#
#
# IMPORTANT
# ---------
# This script does NOT:
#
#   - download or reprocess GSE185263;
#   - normalize raw counts;
#   - remap Ensembl identifiers;
#   - recalculate the five-gene score;
#   - recalculate z-scores;
#   - select or substitute genes;
#   - refit signature coefficients;
#   - optimize a cutoff;
#   - reverse score direction;
#   - rerun Spearman tests;
#   - rerun Wilcoxon tests;
#   - rerun ROC analyses;
#   - rerun BH correction;
#   - rerun the adjusted linear model;
#   - rerun location-specific associations;
#   - perform new hypothesis testing.
#
#
# VISUALIZATION-ONLY CALCULATIONS
# -------------------------------
#
# Panel A:
#   descriptive linear trend line only.
#
# Panel D:
#   95% CI = frozen beta +/- 1.96 * frozen SE
#   solely for graphical presentation.
#
#
# FINAL VISUAL CHANGES
# --------------------
#
# Panel A:
#   shorter title.
#
# Panel C:
#   shorter title/subtitle;
#   individual q-value labels removed;
#   rho values retained.
#
# Panel D:
#   compact title/subtitle;
#   narrow panel width.
#
# Panel E:
#   shorter title/subtitle;
#   individual q-value labels removed;
#   rho values retained;
#   interpretation focused on directional geographic consistency.
#
# Layout:
#
#   AAAAAA BBBBBB
#   CCCCC  DD EEEEE
#
# i.e. lower-row widths approximately 5 : 2 : 5.
#
#
# FROZEN PANEL
# ------------
# UP:
#   CD177
#   HK3
#   IRAK3
#
# DOWN:
#   CARD11
#   IKZF2
#
#
# FROZEN SCORE
# ------------
# mean[z(CD177), z(HK3), z(IRAK3)] -
# mean[z(CARD11), z(IKZF2)]
#
#
# PRIMARY STANDARDIZATION
# -----------------------
# Gene-wise z-standardization parameters were estimated across ALL 348
# available sepsis samples, independently of SOFA availability.
#
#
# PRIMARY EXTERNAL ENDPOINT
# -------------------------
#
# Five-gene score vs continuous 24-h SOFA:
#
# n = 345
# Spearman rho = 0.31149631009461
# P = 3.36887411623066e-09
#
#
# SECONDARY SOFA ANALYSIS
# -----------------------
#
# SOFA >=2 vs SOFA 0-1
#
# n = 207 vs 138
#
# median:
#   SOFA >=2 = 0.442576323149592
#   SOFA 0-1 = -0.481157996747498
#
# P = 8.95261154099103e-08
# BH = 1.79052230819821e-07
#
# AUC = 0.669887278582931
# 95% CI = 0.612671243297642 - 0.727103313868219
#
#
# COMPONENT GENES vs SOFA
# -----------------------
#
# CD177   rho = +0.257134974578282
# HK3     rho = +0.345021311654618
# IRAK3   rho = +0.308239409184258
# CARD11  rho = -0.215317239074883
# IKZF2   rho = -0.181133513354905
#
# All five:
#   directionally concordant
#   BH-significant within the five-gene family.
#
#
# ADJUSTED MODEL
# --------------
#
# five_gene_score ~ SOFA + age + sex + collection_location
#
# n = 345
#
# SOFA:
# beta = 0.116609263940537
# SE   = 0.0260606000130327
# P    = 1.04867638525246e-05
#
#
# LOCATION-SPECIFIC SOFA ASSOCIATIONS
# -----------------------------------
#
# Australia:
#   n=84
#   rho=0.26015301561854
#
# Colombia:
#   n=67
#   rho=0.0383826200917412
#
# Netherlands:
#   n=104
#   rho=0.279704076367344
#
# Toronto:
#   n=79
#   rho=0.328421558769912
#
# Vancouver:
#   n=11
#   rho=0.25434031378315
#
# All five location-specific estimates are positive.
#
#
# IMPORTANT LOCATION INTERPRETATION
# ---------------------------------
#
# Geographic analyses are sensitivity analyses WITHIN GSE185263.
#
# They are NOT independent validation cohorts.
#
# Script 142b also calculated a descriptive fixed-effect Fisher-z synthesis.
# It is audited/exported but NOT presented as random-effects meta-analysis.
#
#
# OUTPUT
# ------
#
# results/blood_endotypes_biomarkers/
#   150_Figure5_GSE185263_severity_validation/
#
################################################################################


cat("====================================================================\n")
cat("Running Script 150\n")
cat("FINAL VISUAL-POLISH Main Figure 5\n")
cat("GSE185263 organ-dysfunction severity\n")
cat("====================================================================\n\n")


# =============================================================================
# 1. PROJECT DIRECTORY
# =============================================================================

project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))

if (!dir.exists(project_dir)) {
  project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))
}

if (!dir.exists(project_dir)) {
  stop(
    "Sepsis_DESeq2 project directory not found."
  )
}

setwd(
  project_dir
)


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
  "tidyr",
  "ggplot2",
  "patchwork",
  "openxlsx",
  "stringr",
  "forcats",
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
  library(tidyr)
  
  library(ggplot2)
  library(patchwork)
  
  library(openxlsx)
  
  library(stringr)
  library(forcats)
  library(scales)
  
})


# =============================================================================
# 3. HELPERS
# =============================================================================

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
          
          cat(
            "  ",
            h,
            "\n",
            sep = ""
          )
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


format_p <- function(x) {
  
  x <- as.numeric(
    x
  )
  
  
  if (!is.finite(x)) {
    
    return(
      "NA"
    )
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
  
  
  sprintf(
    "%.4f",
    x
  )
}


# =============================================================================
# 4. LOCATE FROZEN SCRIPT 142b WORKBOOK
# =============================================================================

source_workbook <- find_project_file(
  
  candidates = c(
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "142b_external_validation_GSE185263",
      "tables",
      "142b_GSE185263_external_validation.xlsx"
    ),
    
    file.path(
      "results",
      "blood_endotypes_biomarkers",
      "142b_external_validation_GSE185263",
      "142b_GSE185263_external_validation.xlsx"
    )
  ),
  
  recursive_pattern =
    "^142b_GSE185263_external_validation\\.xlsx$",
  
  description =
    "Script 142b GSE185263 external-validation workbook"
)


if (
  is.na(
    source_workbook
  ) ||
  !file.exists(
    source_workbook
  )
) {
  
  stop(
    paste0(
      "Could not locate frozen Script 142b workbook:\n",
      "142b_GSE185263_external_validation.xlsx"
    )
  )
}


cat("\n====================================================================\n")
cat("SOURCE WORKBOOK\n")
cat("====================================================================\n")


cat(
  normalizePath(
    source_workbook,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


# =============================================================================
# 5. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "150_Figure5_GSE185263_severity_validation"
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


logs_dir <- file.path(
  output_dir,
  "logs"
)


for (d in c(
  output_dir,
  figures_dir,
  tables_dir,
  text_dir,
  logs_dir
)) {
  
  dir.create(
    d,
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
# 6. WORKBOOK SHEET AUDIT
# =============================================================================

required_sheets <- c(
  "00_run_info",
  "01_PRIMARY_SOFA",
  "02_secondary",
  "03_gene_SOFA",
  "04_adjusted_model",
  "05_adjusted_summary",
  "06_location_SOFA",
  "07_location_pooled",
  "08_scaling",
  "09_SOFA_summary",
  "10_SOFA_groups",
  "11_mortality",
  "12_site_summary",
  "13_location_summary",
  "14_gene_coverage",
  "15_sample_scores"
)


available_sheets <- openxlsx::getSheetNames(
  source_workbook
)


cat("\nAvailable Script 142b sheets:\n")

print(
  available_sheets
)


missing_sheets <- setdiff(
  required_sheets,
  available_sheets
)


if (length(missing_sheets) > 0) {
  
  stop(
    paste0(
      "Required Script 142b sheet(s) missing:\n",
      paste(
        missing_sheets,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 7. LOAD FROZEN TABLES
# =============================================================================

run_info <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "00_run_info"
)


primary_sofa <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "01_PRIMARY_SOFA"
)


secondary_results <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "02_secondary"
)


gene_sofa <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "03_gene_SOFA"
)


adjusted_model <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "04_adjusted_model"
)


adjusted_summary <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "05_adjusted_summary"
)


location_sofa <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "06_location_SOFA"
)


location_pooled <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "07_location_pooled"
)


scaling <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "08_scaling"
)


sofa_summary <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "09_SOFA_summary"
)


sofa_groups <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "10_SOFA_groups"
)


mortality_summary <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "11_mortality"
)


site_summary <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "12_site_summary"
)


location_summary <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "13_location_summary"
)


gene_coverage <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "14_gene_coverage"
)


sample_scores <- openxlsx::read.xlsx(
  source_workbook,
  sheet = "15_sample_scores"
)


# =============================================================================
# 8. REQUIRED COLUMN AUDIT
# =============================================================================

required_primary_columns <- c(
  "analysis",
  "expected_direction",
  "n",
  "rho",
  "p_value"
)


required_secondary_columns <- c(
  "analysis",
  "score",
  "case",
  "control",
  "n_case",
  "n_control",
  "median_case",
  "median_control",
  "median_difference",
  "p_value",
  "AUC",
  "CI_low",
  "CI_high",
  "BH_secondary"
)


required_gene_columns <- c(
  "gene",
  "expected_direction",
  "observed_direction",
  "direction_concordant",
  "n",
  "rho",
  "p_value",
  "BH_five_genes"
)


required_adjusted_columns <- c(
  "term",
  "estimate",
  "SE",
  "t_value",
  "p_value"
)


required_location_columns <- c(
  "collection_location",
  "direction_concordant",
  "n",
  "rho",
  "p_value",
  "BH_location"
)


required_score_columns <- c(
  "geo_accession",
  "sample_title",
  "disease_state",
  "age_numeric",
  "sex_standardized",
  "collection_location",
  "collection_site",
  "mortality",
  "sofa",
  "five_gene_score",
  "five_gene_score_all_reference",
  "five_gene_score_healthy_reference",
  "sofa_group"
)


check_required_columns <- function(
    data,
    required_columns,
    table_name
) {
  
  missing_columns <- setdiff(
    required_columns,
    names(
      data
    )
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        "Missing columns in ",
        table_name,
        ":\n",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
}


check_required_columns(
  primary_sofa,
  required_primary_columns,
  "01_PRIMARY_SOFA"
)


check_required_columns(
  secondary_results,
  required_secondary_columns,
  "02_secondary"
)


check_required_columns(
  gene_sofa,
  required_gene_columns,
  "03_gene_SOFA"
)


check_required_columns(
  adjusted_model,
  required_adjusted_columns,
  "04_adjusted_model"
)


check_required_columns(
  location_sofa,
  required_location_columns,
  "06_location_SOFA"
)


check_required_columns(
  sample_scores,
  required_score_columns,
  "15_sample_scores"
)


# =============================================================================
# 9. NORMALIZE TYPES
# =============================================================================

primary_sofa <- primary_sofa %>%
  
  dplyr::mutate(
    
    n =
      as.numeric(
        n
      ),
    
    rho =
      as.numeric(
        rho
      ),
    
    p_value =
      as.numeric(
        p_value
      )
  )


secondary_results <- secondary_results %>%
  
  dplyr::mutate(
    
    n_case =
      as.numeric(
        n_case
      ),
    
    n_control =
      as.numeric(
        n_control
      ),
    
    median_case =
      as.numeric(
        median_case
      ),
    
    median_control =
      as.numeric(
        median_control
      ),
    
    median_difference =
      as.numeric(
        median_difference
      ),
    
    p_value =
      as.numeric(
        p_value
      ),
    
    AUC =
      as.numeric(
        AUC
      ),
    
    CI_low =
      as.numeric(
        CI_low
      ),
    
    CI_high =
      as.numeric(
        CI_high
      ),
    
    BH_secondary =
      as.numeric(
        BH_secondary
      )
  )


gene_sofa <- gene_sofa %>%
  
  dplyr::mutate(
    
    n =
      as.numeric(
        n
      ),
    
    rho =
      as.numeric(
        rho
      ),
    
    p_value =
      as.numeric(
        p_value
      ),
    
    BH_five_genes =
      as.numeric(
        BH_five_genes
      ),
    
    direction_concordant =
      as.logical(
        direction_concordant
      )
  )


adjusted_model <- adjusted_model %>%
  
  dplyr::mutate(
    
    estimate =
      as.numeric(
        estimate
      ),
    
    SE =
      as.numeric(
        SE
      ),
    
    t_value =
      as.numeric(
        t_value
      ),
    
    p_value =
      as.numeric(
        p_value
      )
  )


location_sofa <- location_sofa %>%
  
  dplyr::mutate(
    
    n =
      as.numeric(
        n
      ),
    
    rho =
      as.numeric(
        rho
      ),
    
    p_value =
      as.numeric(
        p_value
      ),
    
    BH_location =
      as.numeric(
        BH_location
      ),
    
    direction_concordant =
      as.logical(
        direction_concordant
      )
  )


sample_scores <- sample_scores %>%
  
  dplyr::mutate(
    
    disease_state =
      as.character(
        disease_state
      ),
    
    age_numeric =
      as.numeric(
        age_numeric
      ),
    
    sofa =
      as.numeric(
        sofa
      ),
    
    five_gene_score =
      as.numeric(
        five_gene_score
      ),
    
    five_gene_score_all_reference =
      as.numeric(
        five_gene_score_all_reference
      ),
    
    five_gene_score_healthy_reference =
      as.numeric(
        five_gene_score_healthy_reference
      ),
    
    sofa_group =
      as.character(
        sofa_group
      )
  )


# =============================================================================
# 10. EXTRACT FROZEN RESULTS
# =============================================================================

primary_row <- primary_sofa %>%
  
  dplyr::filter(
    analysis ==
      "PRIMARY_five_gene_score_vs_SOFA"
  )


if (nrow(primary_row) != 1) {
  
  stop(
    "Primary SOFA result not uniquely identified."
  )
}


sofa_binary_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "SOFA_ge2_vs_SOFA_0_1"
  )


mortality_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "Died_vs_Survived"
  )


site_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "ICU_vs_Emergency_Room"
  )


context_row <- secondary_results %>%
  
  dplyr::filter(
    analysis ==
      "Sepsis_vs_healthy_contextual"
  )


if (nrow(sofa_binary_row) != 1) {
  
  stop(
    "SOFA >=2 vs SOFA 0-1 result not uniquely identified."
  )
}


if (nrow(mortality_row) != 1) {
  
  stop(
    "Mortality result not uniquely identified."
  )
}


if (nrow(site_row) != 1) {
  
  stop(
    "Collection-site result not uniquely identified."
  )
}


if (nrow(context_row) != 1) {
  
  stop(
    "Contextual sepsis-vs-healthy result not uniquely identified."
  )
}


adjusted_sofa_row <- adjusted_model %>%
  
  dplyr::filter(
    term ==
      "sofa"
  )


if (nrow(adjusted_sofa_row) != 1) {
  
  stop(
    "Adjusted SOFA coefficient not uniquely identified."
  )
}


# =============================================================================
# 11. ANALYSIS DATASETS
# =============================================================================

sepsis_scores <- sample_scores %>%
  
  dplyr::filter(
    disease_state ==
      "sepsis"
  )


sofa_complete <- sepsis_scores %>%
  
  dplyr::filter(
    is.finite(
      sofa
    ),
    is.finite(
      five_gene_score
    )
  )


if (
  nrow(
    sepsis_scores
  ) != 348
) {
  
  stop(
    paste0(
      "Expected 348 sepsis samples; observed ",
      nrow(
        sepsis_scores
      )
    )
  )
}


if (
  nrow(
    sofa_complete
  ) != 345
) {
  
  stop(
    paste0(
      "Expected 345 SOFA-complete sepsis samples; observed ",
      nrow(
        sofa_complete
      )
    )
  )
}


# =============================================================================
# 12. EXPECTED FROZEN VALUES
# =============================================================================

expected_total_samples <- 392
expected_sepsis_samples <- 348
expected_healthy_samples <- 44
expected_sofa_n <- 345


expected_primary_rho <- 0.31149631009461
expected_primary_p <- 3.36887411623066e-09


expected_sofa_case_n <- 207
expected_sofa_control_n <- 138

expected_sofa_case_median <- 0.442576323149592
expected_sofa_control_median <- -0.481157996747498

expected_sofa_p <- 8.95261154099103e-08
expected_sofa_q <- 1.79052230819821e-07

expected_sofa_auc <- 0.669887278582931
expected_sofa_ci_low <- 0.612671243297642
expected_sofa_ci_high <- 0.727103313868219


expected_adjusted_beta <- 0.116609263940537
expected_adjusted_se <- 0.0260606000130327
expected_adjusted_p <- 1.04867638525246e-05


expected_gene_rho <- c(
  
  CD177 =
    0.257134974578282,
  
  HK3 =
    0.345021311654618,
  
  IRAK3 =
    0.308239409184258,
  
  CARD11 =
    -0.215317239074883,
  
  IKZF2 =
    -0.181133513354905
)


expected_location_rho <- c(
  
  australia =
    0.26015301561854,
  
  colombia =
    0.0383826200917412,
  
  netherlands =
    0.279704076367344,
  
  toronto =
    0.328421558769912,
  
  vancouver =
    0.25434031378315
)


expected_location_n <- c(
  
  australia =
    84,
  
  colombia =
    67,
  
  netherlands =
    104,
  
  toronto =
    79,
  
  vancouver =
    11
)


expected_pooled_rho <- 0.240783921963788
expected_pooled_ci_low <- 0.136847674716661
expected_pooled_ci_high <- 0.339476360533898


# =============================================================================
# 13. NUMERICAL AUDIT
# =============================================================================

observed_total_samples <- nrow(
  sample_scores
)


observed_sepsis_samples <- sum(
  sample_scores$disease_state ==
    "sepsis",
  na.rm = TRUE
)


observed_healthy_samples <- sum(
  sample_scores$disease_state ==
    "healthy",
  na.rm = TRUE
)


if (
  observed_total_samples !=
  expected_total_samples
) {
  
  stop(
    "Total GSE185263 sample-count audit failed."
  )
}


if (
  observed_sepsis_samples !=
  expected_sepsis_samples
) {
  
  stop(
    "GSE185263 sepsis sample-count audit failed."
  )
}


if (
  observed_healthy_samples !=
  expected_healthy_samples
) {
  
  stop(
    "GSE185263 healthy sample-count audit failed."
  )
}


if (
  primary_row$n[1] !=
  expected_sofa_n
) {
  
  stop(
    "Primary SOFA n mismatch."
  )
}


if (
  abs(
    primary_row$rho[1] -
    expected_primary_rho
  ) >
  1e-12
) {
  
  stop(
    "Primary score-SOFA rho mismatch."
  )
}


if (
  abs(
    primary_row$p_value[1] -
    expected_primary_p
  ) >
  1e-14
) {
  
  stop(
    "Primary score-SOFA P mismatch."
  )
}


sofa_binary_checks <- c(
  
  sofa_binary_row$n_case[1] ==
    expected_sofa_case_n,
  
  sofa_binary_row$n_control[1] ==
    expected_sofa_control_n,
  
  abs(
    sofa_binary_row$median_case[1] -
      expected_sofa_case_median
  ) <
    1e-12,
  
  abs(
    sofa_binary_row$median_control[1] -
      expected_sofa_control_median
  ) <
    1e-12,
  
  abs(
    sofa_binary_row$p_value[1] -
      expected_sofa_p
  ) <
    1e-14,
  
  abs(
    sofa_binary_row$BH_secondary[1] -
      expected_sofa_q
  ) <
    1e-14,
  
  abs(
    sofa_binary_row$AUC[1] -
      expected_sofa_auc
  ) <
    1e-12,
  
  abs(
    sofa_binary_row$CI_low[1] -
      expected_sofa_ci_low
  ) <
    1e-12,
  
  abs(
    sofa_binary_row$CI_high[1] -
      expected_sofa_ci_high
  ) <
    1e-12
)


if (!all(
  sofa_binary_checks
)) {
  
  stop(
    "Secondary SOFA binary audit failed."
  )
}


if (
  nrow(
    gene_sofa
  ) != 5
) {
  
  stop(
    "Expected five component-gene SOFA results."
  )
}


for (
  gene_name in names(
    expected_gene_rho
  )
) {
  
  gene_row <- gene_sofa %>%
    
    dplyr::filter(
      gene ==
        gene_name
    )
  
  
  if (
    nrow(
      gene_row
    ) != 1
  ) {
    
    stop(
      paste0(
        "Missing or duplicated gene result: ",
        gene_name
      )
    )
  }
  
  
  if (
    abs(
      gene_row$rho[1] -
      expected_gene_rho[gene_name]
    ) >
    1e-12
  ) {
    
    stop(
      paste0(
        "SOFA rho mismatch for ",
        gene_name
      )
    )
  }
}


n_gene_direction_concordant <- sum(
  gene_sofa$direction_concordant,
  na.rm = TRUE
)


n_gene_BH_significant <- sum(
  gene_sofa$BH_five_genes <
    0.05,
  na.rm = TRUE
)


if (
  n_gene_direction_concordant != 5
) {
  
  stop(
    "Not all five component genes retain prespecified SOFA direction."
  )
}


if (
  n_gene_BH_significant != 5
) {
  
  stop(
    paste0(
      "Expected 5/5 component genes significant after BH; observed ",
      n_gene_BH_significant,
      "/5."
    )
  )
}


adjusted_n_row <- adjusted_summary %>%
  
  dplyr::filter(
    metric ==
      "n"
  )


if (
  nrow(
    adjusted_n_row
  ) != 1
) {
  
  stop(
    "Adjusted model n not uniquely identified."
  )
}


if (
  as.numeric(
    adjusted_n_row$value[1]
  ) != 345
) {
  
  stop(
    "Adjusted model n mismatch."
  )
}


if (
  abs(
    adjusted_sofa_row$estimate[1] -
    expected_adjusted_beta
  ) >
  1e-12
) {
  
  stop(
    "Adjusted SOFA beta mismatch."
  )
}


if (
  abs(
    adjusted_sofa_row$SE[1] -
    expected_adjusted_se
  ) >
  1e-12
) {
  
  stop(
    "Adjusted SOFA SE mismatch."
  )
}


if (
  abs(
    adjusted_sofa_row$p_value[1] -
    expected_adjusted_p
  ) >
  1e-14
) {
  
  stop(
    "Adjusted SOFA P mismatch."
  )
}


if (
  nrow(
    location_sofa
  ) != 5
) {
  
  stop(
    paste0(
      "Expected 5 eligible locations; observed ",
      nrow(
        location_sofa
      )
    )
  )
}


for (
  location_name in names(
    expected_location_rho
  )
) {
  
  location_row <- location_sofa %>%
    
    dplyr::filter(
      collection_location ==
        location_name
    )
  
  
  if (
    nrow(
      location_row
    ) != 1
  ) {
    
    stop(
      paste0(
        "Missing or duplicated location: ",
        location_name
      )
    )
  }
  
  
  if (
    location_row$n[1] !=
    expected_location_n[location_name]
  ) {
    
    stop(
      paste0(
        "Location n mismatch for ",
        location_name
      )
    )
  }
  
  
  if (
    abs(
      location_row$rho[1] -
      expected_location_rho[location_name]
    ) >
    1e-12
  ) {
    
    stop(
      paste0(
        "Location rho mismatch for ",
        location_name
      )
    )
  }
}


n_positive_locations <- sum(
  location_sofa$rho > 0,
  na.rm = TRUE
)


if (
  n_positive_locations != 5
) {
  
  stop(
    "Expected all five location-specific correlations to be positive."
  )
}


if (
  nrow(
    location_pooled
  ) != 1
) {
  
  stop(
    "Location pooled result not uniquely identified."
  )
}


if (
  abs(
    as.numeric(
      location_pooled$pooled_rho[1]
    ) -
    expected_pooled_rho
  ) >
  1e-12
) {
  
  stop(
    "Descriptive fixed-effect pooled rho mismatch."
  )
}


if (
  abs(
    as.numeric(
      location_pooled$CI_low[1]
    ) -
    expected_pooled_ci_low
  ) >
  1e-12
) {
  
  stop(
    "Descriptive fixed-effect pooled CI low mismatch."
  )
}


if (
  abs(
    as.numeric(
      location_pooled$CI_high[1]
    ) -
    expected_pooled_ci_high
  ) >
  1e-12
) {
  
  stop(
    "Descriptive fixed-effect pooled CI high mismatch."
  )
}


audit_table <- data.frame(
  
  metric = c(
    "Total GSE185263 samples",
    "Sepsis samples",
    "Healthy samples",
    "SOFA-complete sepsis samples",
    "Primary score-SOFA rho",
    "Primary score-SOFA P",
    "SOFA >=2 n",
    "SOFA 0-1 n",
    "SOFA >=2 median",
    "SOFA 0-1 median",
    "SOFA binary P",
    "SOFA binary BH",
    "SOFA binary AUC",
    "SOFA binary AUC CI low",
    "SOFA binary AUC CI high",
    "Directionally concordant component genes",
    "BH-significant component genes",
    "Adjusted model n",
    "Adjusted SOFA beta",
    "Adjusted SOFA SE",
    "Adjusted SOFA P",
    "Eligible geographic locations",
    "Positive location-specific rho estimates",
    "Descriptive fixed-effect pooled rho",
    "Descriptive fixed-effect pooled CI low",
    "Descriptive fixed-effect pooled CI high"
  ),
  
  observed = c(
    observed_total_samples,
    observed_sepsis_samples,
    observed_healthy_samples,
    nrow(
      sofa_complete
    ),
    primary_row$rho[1],
    primary_row$p_value[1],
    sofa_binary_row$n_case[1],
    sofa_binary_row$n_control[1],
    sofa_binary_row$median_case[1],
    sofa_binary_row$median_control[1],
    sofa_binary_row$p_value[1],
    sofa_binary_row$BH_secondary[1],
    sofa_binary_row$AUC[1],
    sofa_binary_row$CI_low[1],
    sofa_binary_row$CI_high[1],
    n_gene_direction_concordant,
    n_gene_BH_significant,
    as.numeric(
      adjusted_n_row$value[1]
    ),
    adjusted_sofa_row$estimate[1],
    adjusted_sofa_row$SE[1],
    adjusted_sofa_row$p_value[1],
    nrow(
      location_sofa
    ),
    n_positive_locations,
    as.numeric(
      location_pooled$pooled_rho[1]
    ),
    as.numeric(
      location_pooled$CI_low[1]
    ),
    as.numeric(
      location_pooled$CI_high[1]
    )
  ),
  
  expected = c(
    expected_total_samples,
    expected_sepsis_samples,
    expected_healthy_samples,
    expected_sofa_n,
    expected_primary_rho,
    expected_primary_p,
    expected_sofa_case_n,
    expected_sofa_control_n,
    expected_sofa_case_median,
    expected_sofa_control_median,
    expected_sofa_p,
    expected_sofa_q,
    expected_sofa_auc,
    expected_sofa_ci_low,
    expected_sofa_ci_high,
    5,
    5,
    345,
    expected_adjusted_beta,
    expected_adjusted_se,
    expected_adjusted_p,
    5,
    5,
    expected_pooled_rho,
    expected_pooled_ci_low,
    expected_pooled_ci_high
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


cat("\nNumerical audit passed successfully.\n")


# =============================================================================
# 14. COLORS
# =============================================================================

col_primary <- "#D55E00"

col_low_sofa <- "#56B4E9"
col_high_sofa <- "#D55E00"

col_up_gene <- "#D55E00"
col_down_gene <- "#0072B2"

col_location <- "#0072B2"

col_adjusted <- "#009E73"

col_reference <- "#555555"


# =============================================================================
# 15. PANEL A — PRIMARY CONTINUOUS SOFA
# =============================================================================

panel_A_annotation <- paste0(
  
  "Spearman \u03c1 = ",
  sprintf(
    "%.3f",
    primary_row$rho[1]
  ),
  
  "\nP = ",
  format_p(
    primary_row$p_value[1]
  )
)


p_A <- ggplot2::ggplot(
  
  sofa_complete,
  
  ggplot2::aes(
    x = sofa,
    y = five_gene_score
  )
  
) +
  
  ggplot2::geom_jitter(
    width = 0.10,
    height = 0,
    size = 1.8,
    alpha = 0.58,
    color = col_primary
  ) +
  
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.80,
    color = "#4D4D4D"
  ) +
  
  ggplot2::annotate(
    "label",
    x = Inf,
    y = -Inf,
    label = panel_A_annotation,
    hjust = 1.05,
    vjust = -0.15,
    size = 3.3,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = seq(
      0,
      16,
      by = 2
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    
    tag =
      "A",
    
    title =
      "Five-gene score tracks continuous SOFA severity",
    
    subtitle =
      "GSE185263 sepsis samples with 24-h SOFA; n=345",
    
    x =
      "SOFA score at 24 h",
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 16. PANEL B — SECONDARY SOFA THRESHOLD
# =============================================================================

panel_B_data <- sofa_complete %>%
  
  dplyr::mutate(
    
    SOFA_group =
      dplyr::case_when(
        
        sofa >= 2 ~
          "SOFA \u22652",
        
        sofa <= 1 ~
          "SOFA 0\u20131",
        
        TRUE ~
          NA_character_
      ),
    
    SOFA_group =
      factor(
        SOFA_group,
        levels = c(
          "SOFA 0\u20131",
          "SOFA \u22652"
        )
      )
  )


panel_B_annotation <- paste0(
  
  "Median: ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_control[1]
  ),
  " vs ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_case[1]
  ),
  
  "\nP = ",
  format_p(
    sofa_binary_row$p_value[1]
  ),
  
  "\nAUC = ",
  sprintf(
    "%.3f",
    sofa_binary_row$AUC[1]
  ),
  
  " [",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_low[1]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_high[1]
  ),
  "]"
)


sofa_group_colors <- c(
  
  "SOFA 0\u20131" =
    col_low_sofa,
  
  "SOFA \u22652" =
    col_high_sofa
)


p_B <- ggplot2::ggplot(
  
  panel_B_data,
  
  ggplot2::aes(
    x = SOFA_group,
    y = five_gene_score,
    fill = SOFA_group,
    color = SOFA_group
  )
  
) +
  
  ggplot2::geom_boxplot(
    width = 0.58,
    alpha = 0.40,
    outlier.shape = NA,
    linewidth = 0.60
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.7,
    alpha = 0.60
  ) +
  
  ggplot2::scale_fill_manual(
    values = sofa_group_colors
  ) +
  
  ggplot2::scale_color_manual(
    values = sofa_group_colors
  ) +
  
  ggplot2::annotate(
    "label",
    x = 1.02,
    y = Inf,
    label = panel_B_annotation,
    hjust = 0,
    vjust = 1.10,
    size = 3.0,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    legend.position = "none"
  ) +
  
  ggplot2::labs(
    
    tag =
      "B",
    
    title =
      "Secondary SOFA threshold analysis",
    
    subtitle =
      paste0(
        "SOFA 0\u20131 (n=",
        sofa_binary_row$n_control[1],
        ") vs SOFA \u22652 (n=",
        sofa_binary_row$n_case[1],
        ")"
      ),
    
    x =
      NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 17. PANEL C — COMPONENT GENES
# =============================================================================

gene_order <- c(
  "HK3",
  "IRAK3",
  "CD177",
  "CARD11",
  "IKZF2"
)


gene_sofa_plot <- gene_sofa %>%
  
  dplyr::mutate(
    
    gene =
      factor(
        gene,
        levels =
          rev(
            gene_order
          )
      ),
    
    Program =
      dplyr::if_else(
        expected_direction ==
          "POSITIVE",
        "Increased component",
        "Decreased component"
      )
  )


gene_program_colors <- c(
  
  "Increased component" =
    col_up_gene,
  
  "Decreased component" =
    col_down_gene
)


p_C <- ggplot2::ggplot(
  
  gene_sofa_plot,
  
  ggplot2::aes(
    x = rho,
    y = gene
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.55,
    color = "grey45"
  ) +
  
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = rho,
      yend = gene,
      color = Program
    ),
    linewidth = 1.0
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      color = Program
    ),
    size = 4.1
  ) +
  
  ggplot2::geom_text(
    
    data =
      gene_sofa_plot %>%
      dplyr::filter(
        rho >= 0
      ),
    
    ggplot2::aes(
      label =
        sprintf(
          "%.3f",
          rho
        )
    ),
    
    hjust = -0.35,
    size = 3.1
  ) +
  
  ggplot2::geom_text(
    
    data =
      gene_sofa_plot %>%
      dplyr::filter(
        rho < 0
      ),
    
    ggplot2::aes(
      label =
        sprintf(
          "%.3f",
          rho
        )
    ),
    
    hjust = 1.35,
    size = 3.1
  ) +
  
  ggplot2::scale_color_manual(
    values = gene_program_colors
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      -0.42,
      0.42
    ),
    breaks = seq(
      -0.4,
      0.4,
      by = 0.2
    )
  ) +
  
  theme_publication(
    9.5
  ) +
  
  ggplot2::theme(
    
    axis.title.y =
      ggplot2::element_blank(),
    
    legend.position =
      "bottom",
    
    legend.text =
      ggplot2::element_text(
        size = 7.5
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "C",
    
    title =
      "Component genes retain concordant SOFA associations",
    
    subtitle =
      "All five genes remained significant after BH correction",
    
    x =
      "Spearman \u03c1 with SOFA",
    
    color =
      NULL
  )


# =============================================================================
# 18. PANEL D — ADJUSTED SOFA ASSOCIATION
# =============================================================================

adjusted_beta <- adjusted_sofa_row$estimate[1]

adjusted_se <- adjusted_sofa_row$SE[1]


adjusted_ci_low <-
  adjusted_beta -
  1.96 *
  adjusted_se


adjusted_ci_high <-
  adjusted_beta +
  1.96 *
  adjusted_se


adjusted_plot <- data.frame(
  
  term =
    "SOFA",
  
  beta =
    adjusted_beta,
  
  CI_low =
    adjusted_ci_low,
  
  CI_high =
    adjusted_ci_high,
  
  stringsAsFactors = FALSE
)


panel_D_annotation <- paste0(
  
  "\u03b2 = ",
  sprintf(
    "%.3f",
    adjusted_beta
  ),
  
  "\n95% CI ",
  sprintf(
    "%.3f",
    adjusted_ci_low
  ),
  "\u2013",
  sprintf(
    "%.3f",
    adjusted_ci_high
  ),
  
  "\nP = ",
  format_p(
    adjusted_sofa_row$p_value[1]
  )
)


p_D <- ggplot2::ggplot(
  
  adjusted_plot,
  
  ggplot2::aes(
    x = beta,
    y = term
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey55"
  ) +
  
  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    width = 0.16,
    linewidth = 1.0,
    color = col_adjusted
  ) +
  
  ggplot2::geom_point(
    size = 4.5,
    color = col_adjusted
  ) +
  
  ggplot2::annotate(
    "label",
    x = Inf,
    y = 1,
    label = panel_D_annotation,
    hjust = 1.04,
    vjust = -0.45,
    size = 2.8,
    linewidth = 0.25,
    fill = "white"
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      -0.02,
      0.22
    ),
    breaks = seq(
      0,
      0.20,
      by = 0.05
    )
  ) +
  
  theme_publication(
    9.0
  ) +
  
  ggplot2::theme(
    
    axis.title.y =
      ggplot2::element_blank(),
    
    axis.text.y =
      ggplot2::element_text(
        face = "bold"
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 7.8,
        color = "grey25"
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "D",
    
    title =
      "Adjusted SOFA association",
    
    subtitle =
      "Adjusted for age, sex, and location; n=345",
    
    x =
      "Adjusted SOFA coefficient"
  )


# =============================================================================
# 19. PANEL E — GEOGRAPHIC SENSITIVITY
# =============================================================================

location_label_map <- c(
  
  australia =
    "Australia",
  
  colombia =
    "Colombia",
  
  netherlands =
    "Netherlands",
  
  toronto =
    "Toronto",
  
  vancouver =
    "Vancouver"
)


location_plot <- location_sofa %>%
  
  dplyr::mutate(
    
    display_location =
      unname(
        location_label_map[
          collection_location
        ]
      ),
    
    display_label =
      paste0(
        display_location,
        " (n=",
        n,
        ")"
      ),
    
    display_label =
      forcats::fct_reorder(
        display_label,
        rho
      )
  )


p_E <- ggplot2::ggplot(
  
  location_plot,
  
  ggplot2::aes(
    x = rho,
    y = display_label
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept =
      primary_row$rho[1],
    linetype = "dashed",
    linewidth = 0.65,
    color = col_reference
  ) +
  
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = rho,
      yend = display_label
    ),
    linewidth = 0.90,
    color = "#B8B8B8"
  ) +
  
  ggplot2::geom_point(
    size = 4.0,
    color = col_location
  ) +
  
  ggplot2::geom_text(
    
    data =
      location_plot %>%
      dplyr::filter(
        rho < 0.36
      ),
    
    ggplot2::aes(
      label =
        sprintf(
          "%.3f",
          rho
        )
    ),
    
    hjust = -0.35,
    size = 3.0
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      0.40
    ),
    breaks = seq(
      0,
      0.4,
      by = 0.1
    ),
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.04
        )
      )
  ) +
  
  theme_publication(
    9.5
  ) +
  
  ggplot2::theme(
    
    axis.title.y =
      ggplot2::element_blank(),
    
    axis.text.y =
      ggplot2::element_text(
        size = 8.3
      )
  ) +
  
  ggplot2::labs(
    
    tag =
      "E",
    
    title =
      "Geographic sensitivity",
    
    subtitle =
      paste0(
        "All five estimates positive; dashed line = overall \u03c1 = ",
        sprintf(
          "%.3f",
          primary_row$rho[1]
        )
      ),
    
    x =
      "Location-specific Spearman \u03c1"
  )


# =============================================================================
# 20. ASSEMBLE FINAL FIGURE 5
# =============================================================================
#
# Final layout:
#
# AAAAAABBBBBB
# CCCCCDDEEEEE
#
# Lower-row relative widths:
#
# C = 5 units
# D = 2 units
# E = 5 units
#
# =============================================================================

figure5_design <- "
AAAAAABBBBBB
CCCCCDDEEEEE
"


figure5 <- p_A +
  p_B +
  p_C +
  p_D +
  p_E +
  
  patchwork::plot_layout(
    design =
      figure5_design,
    heights = c(
      1.04,
      1.00
    )
  )


# =============================================================================
# 21. EXPORT FINAL FIGURE
# =============================================================================

figure_png <- file.path(
  figures_dir,
  "150_Figure5_GSE185263_severity_validation.png"
)


figure_pdf <- file.path(
  figures_dir,
  "150_Figure5_GSE185263_severity_validation.pdf"
)


figure_tiff <- file.path(
  figures_dir,
  "150_Figure5_GSE185263_severity_validation.tiff"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure5,
  width = 14.8,
  height = 9.3,
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure5,
  width = 14.8,
  height = 9.3,
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
  plot = figure5,
  width = 14.8,
  height = 9.3,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 22. EXPORT INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  
  A_score_vs_SOFA =
    p_A,
  
  B_score_by_SOFA_group =
    p_B,
  
  C_component_gene_SOFA =
    p_C,
  
  D_adjusted_SOFA =
    p_D,
  
  E_location_SOFA =
    p_E
)


for (
  panel_name in names(
    individual_panels
  )
) {
  
  panel_width <- if (
    panel_name ==
    "D_adjusted_SOFA"
  ) {
    4.6
  } else {
    6.7
  }
  
  
  ggplot2::ggsave(
    
    filename =
      file.path(
        figures_dir,
        paste0(
          "150_panel_",
          panel_name,
          ".png"
        )
      ),
    
    plot =
      individual_panels[[panel_name]],
    
    width =
      panel_width,
    
    height =
      5.2,
    
    dpi =
      600,
    
    bg =
      "white"
  )
}


# =============================================================================
# 23. SOURCE WORKBOOK
# =============================================================================

run_info_150 <- data.frame(
  
  item = c(
    "script",
    "figure",
    "status",
    "analysis_mode",
    "source_script",
    "source_dataset",
    "total_samples",
    "sepsis_samples",
    "SOFA_complete_sepsis",
    "primary_score_standardization",
    "primary_panel",
    "primary_score",
    "primary_endpoint",
    "new_hypothesis_testing",
    "feature_selection",
    "coefficient_refitting",
    "cutoff_optimization",
    "score_direction_flipping",
    "Panel_C_individual_q_labels",
    "Panel_D_CI",
    "Panel_E_individual_q_labels",
    "location_interpretation",
    "pooled_location_interpretation",
    "final_layout"
  ),
  
  value = c(
    "150_build_Figure5_GSE185263_severity_validation.R",
    "Main Figure 5",
    "FINAL visual-polish version",
    "publication packaging of frozen Script 142b results",
    "142b_external_validation_GSE185263.R",
    "GSE185263",
    "392",
    "348",
    "345",
    paste0(
      "gene-wise z parameters estimated across all 348 sepsis samples ",
      "independently of SOFA availability"
    ),
    "CD177; HK3; IRAK3; CARD11; IKZF2",
    "mean z(CD177,HK3,IRAK3) - mean z(CARD11,IKZF2)",
    "Spearman five-gene score versus continuous 24-h SOFA in sepsis",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "Removed from final main figure",
    "95% CI calculated as frozen beta +/- 1.96*SE for visualization only",
    "Removed from final main figure",
    "within-cohort geographic sensitivity; NOT independent validation cohorts",
    paste0(
      "Script 142b descriptive fixed-effect Fisher-z result audited/exported ",
      "but not used as independent-cohort meta-analysis"
    ),
    "AAAAAABBBBBB / CCCCCDDEEEEE"
  ),
  
  stringsAsFactors = FALSE
)


source_output_workbook <- file.path(
  tables_dir,
  "150_Figure5_source_data.xlsx"
)


openxlsx::write.xlsx(
  
  list(
    
    Run_info =
      run_info_150,
    
    Numerical_audit =
      audit_table,
    
    Primary_SOFA =
      primary_row,
    
    SOFA_binary =
      sofa_binary_row,
    
    SOFA_sample_data =
      sofa_complete,
    
    Component_gene_SOFA =
      gene_sofa,
    
    Adjusted_model =
      adjusted_model,
    
    Adjusted_model_summary =
      adjusted_summary,
    
    Adjusted_SOFA_display =
      adjusted_plot,
    
    Location_SOFA =
      location_sofa,
    
    Location_pooled_descriptive =
      location_pooled,
    
    Scaling_sensitivity =
      scaling,
    
    SOFA_summary =
      sofa_summary,
    
    SOFA_groups =
      sofa_groups,
    
    Mortality_context =
      mortality_summary,
    
    Site_context =
      site_summary,
    
    Location_summary =
      location_summary,
    
    Frozen_gene_coverage =
      gene_coverage
    
  ),
  
  source_output_workbook,
  
  overwrite = TRUE
)


# =============================================================================
# 24. PROVENANCE
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
    "Primary continuous five-gene score-SOFA association",
    "Secondary SOFA >=2 vs SOFA 0-1 comparison",
    "Individual frozen component genes vs SOFA",
    "Covariate-adjusted SOFA association",
    "Location-specific score-SOFA sensitivity within GSE185263"
  ),
  
  Upstream_source = c(
    "142b sheets 01 + 15",
    "142b sheets 02 + 10 + 15",
    "142b sheet 03",
    "142b sheets 04 + 05",
    "142b sheet 06"
  ),
  
  New_hypothesis_testing = rep(
    "NO",
    5
  ),
  
  Visualization_only_calculation = c(
    "Descriptive LM trend line",
    "NO",
    "NO",
    "95% CI from frozen beta +/- 1.96*SE",
    "NO"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  provenance,
  file.path(
    tables_dir,
    "150_Figure5_provenance.csv"
  ),
  row.names = FALSE
)


write.csv(
  audit_table,
  file.path(
    tables_dir,
    "150_Figure5_numerical_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 25. FINAL FIGURE LEGEND — ENGLISH
# =============================================================================

caption_en <- paste0(
  
  "Figure 5. The five-gene host-response score tracks organ-dysfunction ",
  "severity in an independent sepsis RNA-seq cohort. ",
  
  "(A) Association between the frozen five-gene host-response score and ",
  "continuous Sequential Organ Failure Assessment (SOFA) score at 24 h in ",
  "GSE185263. Primary gene-wise standardization parameters were estimated ",
  "across all 348 sepsis samples independently of SOFA availability. Among ",
  "345 sepsis samples with available SOFA measurements, the score was ",
  "positively associated with organ-dysfunction severity (Spearman rho=",
  sprintf(
    "%.3f",
    primary_row$rho[1]
  ),
  ", P=",
  format(
    primary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "). The fitted line is shown as a descriptive visual guide; inferential ",
  "statistics are based on Spearman rank correlation. ",
  
  "(B) Secondary comparison between patients with SOFA 0-1 (n=",
  sofa_binary_row$n_control[1],
  ") and SOFA >=2 (n=",
  sofa_binary_row$n_case[1],
  "). Median five-gene scores were ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_control[1]
  ),
  " and ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_case[1]
  ),
  ", respectively (P=",
  format(
    sofa_binary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; BH-adjusted P=",
  format(
    sofa_binary_row$BH_secondary[1],
    scientific = TRUE,
    digits = 3
  ),
  "), with a fixed-direction AUC of ",
  sprintf(
    "%.3f",
    sofa_binary_row$AUC[1]
  ),
  " (95% CI ",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_low[1]
  ),
  "-",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_high[1]
  ),
  "). ",
  
  "(C) Spearman associations between continuous SOFA and each frozen ",
  "signature component. CD177, HK3, and IRAK3 were positively associated ",
  "with SOFA, whereas CARD11 and IKZF2 were negatively associated, preserving ",
  "all five prespecified directions. All five associations remained significant ",
  "after Benjamini-Hochberg correction within the component-gene family. ",
  
  "(D) Association between SOFA and the five-gene score after adjustment for ",
  "age, sex, and collection location (n=345). The adjusted SOFA coefficient ",
  "was beta=",
  sprintf(
    "%.3f",
    adjusted_beta
  ),
  " score units per SOFA point (SE=",
  sprintf(
    "%.3f",
    adjusted_se
  ),
  "; P=",
  format(
    adjusted_sofa_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "). The displayed 95% confidence interval was calculated from the frozen ",
  "coefficient and standard error for visualization. ",
  
  "(E) Geographic sensitivity of the score-SOFA relationship across eligible ",
  "recruitment locations within GSE185263. All five location-specific Spearman ",
  "estimates were positive; the dashed line denotes the overall primary ",
  "correlation. These location-stratified analyses represent within-cohort ",
  "sensitivity analyses and should not be interpreted as independent external ",
  "validation cohorts. Together, these findings support external replication ",
  "of the association between the frozen molecular host-response score and ",
  "organ-dysfunction severity rather than validation of a calibrated clinical ",
  "prognostic or diagnostic assay."
)


writeLines(
  caption_en,
  file.path(
    text_dir,
    "150_Figure5_caption_EN.txt"
  )
)


# =============================================================================
# 26. FINAL FIGURE LEGEND — RUSSIAN
# =============================================================================

caption_ru <- paste0(
  
  "Рисунок 5. Пятигенный host-response score отражает тяжесть органной ",
  "дисфункции в независимой RNA-seq когорте пациентов с сепсисом. ",
  
  "(A) Связь замороженного пятигенного host-response score с непрерывным ",
  "SOFA через 24 часа в GSE185263. Параметры primary gene-wise стандартизации ",
  "были оценены по всем 348 образцам сепсиса независимо от наличия SOFA. ",
  "Среди 345 пациентов с доступным SOFA score положительно коррелировал с ",
  "тяжестью органной дисфункции (Spearman rho=",
  sprintf(
    "%.3f",
    primary_row$rho[1]
  ),
  ", P=",
  format(
    primary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "). Линия является описательным визуальным трендом; статистический вывод ",
  "основан на ранговой корреляции Spearman. ",
  
  "(B) Вторичное сравнение пациентов с SOFA 0-1 (n=",
  sofa_binary_row$n_control[1],
  ") и SOFA >=2 (n=",
  sofa_binary_row$n_case[1],
  "). Медианы score составляли ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_control[1]
  ),
  " и ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_case[1]
  ),
  " соответственно (P=",
  format(
    sofa_binary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; BH-adjusted P=",
  format(
    sofa_binary_row$BH_secondary[1],
    scientific = TRUE,
    digits = 3
  ),
  "), fixed-direction AUC = ",
  sprintf(
    "%.3f",
    sofa_binary_row$AUC[1]
  ),
  " (95% ДИ ",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_low[1]
  ),
  "-",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_high[1]
  ),
  "). ",
  
  "(C) Корреляции Spearman между непрерывным SOFA и каждым из пяти ",
  "замороженных компонентов. CD177, HK3 и IRAK3 имели положительные, а ",
  "CARD11 и IKZF2 отрицательные ассоциации, сохраняя все пять заранее ",
  "определенных направлений. Все ассоциации сохраняли значимость после ",
  "Benjamini-Hochberg correction внутри component-gene family. ",
  
  "(D) Связь SOFA с пятигенным score после коррекции по возрасту, полу и ",
  "месту набора (n=345). Скорректированный коэффициент SOFA составлял beta=",
  sprintf(
    "%.3f",
    adjusted_beta
  ),
  " единиц score на один балл SOFA (SE=",
  sprintf(
    "%.3f",
    adjusted_se
  ),
  "; P=",
  format(
    adjusted_sofa_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "). Отображаемый 95% ДИ рассчитан из замороженных коэффициента и SE ",
  "исключительно для визуализации. ",
  
  "(E) Географический sensitivity analysis связи score-SOFA по eligible ",
  "местам набора внутри GSE185263. Во всех пяти location-specific analyses ",
  "Spearman rho был положительным; пунктирная линия обозначает общую primary ",
  "correlation. Эти location-stratified analyses являются sensitivity analyses ",
  "внутри одной когорты и не должны трактоваться как независимые external ",
  "validation cohorts. Совокупность результатов подтверждает внешнюю ",
  "репликацию связи molecular host-response score с тяжестью органной ",
  "дисфункции, но не является валидацией калиброванного клинического ",
  "прогностического или диагностического теста."
)


writeLines(
  caption_ru,
  file.path(
    text_dir,
    "150_Figure5_caption_RU.txt"
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
    "150_Figure5_caption_EN_RU.txt"
  )
)


# =============================================================================
# 27. RESULTS SECTION 3.9 — ENGLISH
# =============================================================================

results_39_en <- paste0(
  
  "The association of the frozen five-gene host-response signature with ",
  "organ-dysfunction severity was evaluated in the independent whole-blood ",
  "RNA-seq dataset GSE185263. The dataset comprised 348 sepsis samples and ",
  "44 healthy controls; primary gene-wise standardization parameters were ",
  "estimated across all sepsis samples independently of SOFA availability. ",
  "Among 345 sepsis samples with available 24-h SOFA measurements, the ",
  "five-gene score was positively associated with continuous SOFA ",
  "(Spearman rho=",
  sprintf(
    "%.3f",
    primary_row$rho[1]
  ),
  ", P=",
  format(
    primary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; Fig. 5A). ",
  
  "In a prespecified secondary analysis, patients with SOFA >=2 had a higher ",
  "score than those with SOFA 0-1 (median ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_case[1]
  ),
  " versus ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_control[1]
  ),
  "; P=",
  format(
    sofa_binary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; BH-adjusted P=",
  format(
    sofa_binary_row$BH_secondary[1],
    scientific = TRUE,
    digits = 3
  ),
  "), with a fixed-direction AUC of ",
  sprintf(
    "%.3f",
    sofa_binary_row$AUC[1]
  ),
  " (95% CI ",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_low[1]
  ),
  "-",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_high[1]
  ),
  ") (Fig. 5B). ",
  
  "All five frozen component genes retained their prespecified direction of ",
  "association with SOFA, and all five remained significant after correction ",
  "within the component-gene family (Fig. 5C). The score-SOFA relationship ",
  "also persisted after adjustment for age, sex, and collection location ",
  "(beta=",
  sprintf(
    "%.3f",
    adjusted_sofa_row$estimate[1]
  ),
  " score units per SOFA point, SE=",
  sprintf(
    "%.3f",
    adjusted_sofa_row$SE[1]
  ),
  ", P=",
  format(
    adjusted_sofa_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; Fig. 5D). ",
  
  "Location-stratified analyses showed positive score-SOFA associations in ",
  "all five eligible geographic groups (Fig. 5E), although the magnitude and ",
  "statistical precision varied by location. These analyses were treated as ",
  "within-cohort geographic sensitivity analyses rather than independent ",
  "validation cohorts. Secondary associations with mortality and collection ",
  "site were retained as contextual analyses and were not interpreted as ",
  "validated prognostic or triage performance. Overall, GSE185263 independently ",
  "replicated the relationship between the frozen five-gene host-response score ",
  "and organ-dysfunction severity."
)


writeLines(
  results_39_en,
  file.path(
    text_dir,
    "150_Results_3.9_EN.txt"
  )
)


# =============================================================================
# 28. RESULTS SECTION 3.9 — RUSSIAN
# =============================================================================

results_39_ru <- paste0(
  
  "Связь замороженной пятигенной host-response сигнатуры с тяжестью органной ",
  "дисфункции была оценена в независимом whole-blood RNA-seq наборе ",
  "GSE185263. Набор включал 348 образцов сепсиса и 44 здоровых контроля; ",
  "параметры primary gene-wise стандартизации оценивались по всем образцам ",
  "сепсиса независимо от наличия SOFA. Среди 345 пациентов с доступным ",
  "24-часовым SOFA пятигенный score положительно коррелировал с непрерывным ",
  "SOFA (Spearman rho=",
  sprintf(
    "%.3f",
    primary_row$rho[1]
  ),
  ", P=",
  format(
    primary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; рис. 5A). ",
  
  "Во вторичном заранее определенном анализе пациенты с SOFA >=2 имели более ",
  "высокий score, чем пациенты с SOFA 0-1 (медианы ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_case[1]
  ),
  " против ",
  sprintf(
    "%.3f",
    sofa_binary_row$median_control[1]
  ),
  "; P=",
  format(
    sofa_binary_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; BH-adjusted P=",
  format(
    sofa_binary_row$BH_secondary[1],
    scientific = TRUE,
    digits = 3
  ),
  "), fixed-direction AUC = ",
  sprintf(
    "%.3f",
    sofa_binary_row$AUC[1]
  ),
  " (95% ДИ ",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_low[1]
  ),
  "-",
  sprintf(
    "%.3f",
    sofa_binary_row$CI_high[1]
  ),
  ") (рис. 5B). ",
  
  "Все пять компонентов сохранили заранее определенное направление связи с ",
  "SOFA, и все пять оставались статистически значимыми после коррекции внутри ",
  "component-gene family (рис. 5C). Связь score-SOFA также сохранялась после ",
  "коррекции по возрасту, полу и месту набора (beta=",
  sprintf(
    "%.3f",
    adjusted_sofa_row$estimate[1]
  ),
  " единиц score на один балл SOFA, SE=",
  sprintf(
    "%.3f",
    adjusted_sofa_row$SE[1]
  ),
  ", P=",
  format(
    adjusted_sofa_row$p_value[1],
    scientific = TRUE,
    digits = 3
  ),
  "; рис. 5D). ",
  
  "Location-stratified analyses показали положительную связь score-SOFA во ",
  "всех пяти eligible geographic groups (рис. 5E), хотя величина эффектов и ",
  "статистическая точность различались. Эти анализы рассматривались как ",
  "географические sensitivity analyses внутри одной когорты, а не как ",
  "независимые validation cohorts. Вторичные ассоциации с mortality и местом ",
  "сбора сохранялись как contextual analyses и не трактовались как ",
  "валидированная прогностическая или triage performance. В целом GSE185263 ",
  "независимо реплицировал связь замороженного пятигенного host-response score ",
  "с тяжестью органной дисфункции."
)


writeLines(
  results_39_ru,
  file.path(
    text_dir,
    "150_Results_3.9_RU.txt"
  )
)


# =============================================================================
# 29. FIGURE-TO-RESULTS PLACEMENT
# =============================================================================

placement <- data.frame(
  
  Results_section = rep(
    "3.9",
    5
  ),
  
  Panel = c(
    "5A",
    "5B",
    "5C",
    "5D",
    "5E"
  ),
  
  Content = c(
    "Primary continuous five-gene score-SOFA association",
    "Secondary SOFA >=2 vs SOFA 0-1 analysis",
    "Individual frozen component genes vs SOFA",
    "Covariate-adjusted SOFA association",
    "Geographic sensitivity within GSE185263"
  ),
  
  Recommended_first_citation = c(
    "Sentence reporting primary Spearman score-SOFA association",
    "Sentence reporting secondary SOFA threshold result and AUC",
    "Sentence reporting 5/5 directionally concordant and BH-significant genes",
    "Sentence reporting age/sex/location-adjusted SOFA coefficient",
    "Sentence reporting all-positive location-specific associations"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  placement,
  file.path(
    tables_dir,
    "150_Figure5_Results_placement.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 30. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "150_sessionInfo.txt"
  )
)


# =============================================================================
# 31. FINAL CONSOLE REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 150 completed successfully.\n")
cat("FINAL VISUAL-POLISH Figure 5 generated.\n")
cat("====================================================================\n\n")


cat("GSE185263 COHORT\n")
cat("---------------\n")


cat(
  "Total samples = ",
  observed_total_samples,
  "\n",
  sep = ""
)


cat(
  "Sepsis = ",
  observed_sepsis_samples,
  "\n",
  sep = ""
)


cat(
  "Healthy = ",
  observed_healthy_samples,
  "\n",
  sep = ""
)


cat(
  "Sepsis with 24-h SOFA = ",
  nrow(
    sofa_complete
  ),
  "\n",
  sep = ""
)


cat("\nPRIMARY EXTERNAL SEVERITY ENDPOINT\n")
cat("----------------------------------\n")


cat(
  "n = ",
  primary_row$n[1],
  "\n",
  sep = ""
)


cat(
  "Spearman rho = ",
  sprintf(
    "%.6f",
    primary_row$rho[1]
  ),
  "\n",
  sep = ""
)


cat(
  "P = ",
  format(
    primary_row$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat("\nSECONDARY SOFA THRESHOLD ANALYSIS\n")
cat("---------------------------------\n")


cat(
  "n = ",
  sofa_binary_row$n_case[1],
  " vs ",
  sofa_binary_row$n_control[1],
  "\n",
  sep = ""
)


cat(
  "Median = ",
  sprintf(
    "%.6f",
    sofa_binary_row$median_case[1]
  ),
  " vs ",
  sprintf(
    "%.6f",
    sofa_binary_row$median_control[1]
  ),
  "\n",
  sep = ""
)


cat(
  "P = ",
  format(
    sofa_binary_row$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat(
  "BH = ",
  format(
    sofa_binary_row$BH_secondary[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat(
  "AUC = ",
  sprintf(
    "%.6f",
    sofa_binary_row$AUC[1]
  ),
  "\n",
  sep = ""
)


cat(
  "95% CI = [",
  sprintf(
    "%.6f",
    sofa_binary_row$CI_low[1]
  ),
  ", ",
  sprintf(
    "%.6f",
    sofa_binary_row$CI_high[1]
  ),
  "]\n",
  sep = ""
)


cat("\nCOMPONENT-GENE SOFA ASSOCIATIONS\n")
cat("--------------------------------\n")


print(
  gene_sofa %>%
    
    dplyr::select(
      gene,
      expected_direction,
      observed_direction,
      n,
      rho,
      p_value,
      BH_five_genes
    ),
  row.names = FALSE
)


cat(
  "\nDirectionally concordant genes = ",
  n_gene_direction_concordant,
  "/5\n",
  sep = ""
)


cat(
  "BH-significant component genes = ",
  n_gene_BH_significant,
  "/5\n",
  sep = ""
)


cat("\nADJUSTED SOFA ASSOCIATION\n")
cat("-------------------------\n")


cat(
  "SOFA beta = ",
  sprintf(
    "%.6f",
    adjusted_sofa_row$estimate[1]
  ),
  "\n",
  sep = ""
)


cat(
  "SE = ",
  sprintf(
    "%.6f",
    adjusted_sofa_row$SE[1]
  ),
  "\n",
  sep = ""
)


cat(
  "Display 95% CI = [",
  sprintf(
    "%.6f",
    adjusted_ci_low
  ),
  ", ",
  sprintf(
    "%.6f",
    adjusted_ci_high
  ),
  "]\n",
  sep = ""
)


cat(
  "P = ",
  format(
    adjusted_sofa_row$p_value[1],
    scientific = TRUE,
    digits = 6
  ),
  "\n",
  sep = ""
)


cat("\nGEOGRAPHIC SENSITIVITY\n")
cat("----------------------\n")


print(
  location_sofa %>%
    
    dplyr::select(
      collection_location,
      n,
      rho,
      p_value,
      BH_location,
      direction_concordant
    ),
  row.names = FALSE
)


cat(
  "\nPositive location-specific rho estimates = ",
  n_positive_locations,
  "/5\n",
  sep = ""
)


cat("\nFINAL VISUAL STATUS\n")
cat("-------------------\n")


cat(
  "Panel C individual q labels: REMOVED\n"
)


cat(
  "Panel E individual q labels: REMOVED\n"
)


cat(
  "Panel D compact width: YES\n"
)


cat(
  "Final lower-row width ratio: 5 : 2 : 5\n"
)


cat(
  "Global figure title inside image: NO\n"
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
    source_output_workbook,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nFIGURE LEGEND\n")
cat("-------------\n")


cat(
  normalizePath(
    file.path(
      text_dir,
      "150_Figure5_caption_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nRESULTS 3.9\n")
cat("-----------\n")


cat(
  normalizePath(
    file.path(
      text_dir,
      "150_Results_3.9_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n"
)


cat("\nFIGURE PLACEMENT\n")
cat("----------------\n")


cat(
  "Figure 5A-E -> Results Section 3.9\n"
)


cat(
  "Physical placement -> end of Section 3.9 / end of Results\n"
)


cat("\nINTERPRETATION GUARDRAILS\n")
cat("-------------------------\n")


cat(
  "- Primary endpoint = continuous SOFA association.\n"
)


cat(
  "- SOFA >=2 analysis = SECONDARY.\n"
)


cat(
  "- Mortality analysis = SECONDARY / supplementary-contextual.\n"
)


cat(
  "- Geographic groups are NOT independent validation cohorts.\n"
)


cat(
  "- Fixed-effect Fisher-z synthesis is DESCRIPTIVE only.\n"
)


cat(
  "- Do NOT describe the score as a validated prognostic or diagnostic assay.\n"
)


cat(
  "- Figure 5 completes the Results.\n"
)


cat("\n")
cat("If visual inspection passes, freeze as:\n")
cat("FIGURE 5 — FROZEN\n")
cat("Analytics frozen | source data frozen | composition frozen\n")
cat("\nDone.\n")