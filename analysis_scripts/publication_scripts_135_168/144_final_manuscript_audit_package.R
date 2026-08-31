# ==============================================================================
# Script 144
# Final manuscript audit and publication package
# Blood transcriptomic endotypes + frozen five-gene host-response signature
#
# Project: Sepsis_DESeq2
#
# PURPOSE
# -------
# This script performs a FINAL PUBLICATION AUDIT of the blood-only manuscript.
#
# It does NOT:
#   - select new genes
#   - refit the five-gene score
#   - optimize diagnostic thresholds
#   - rerun exploratory biomarker searches
#   - create new endotypes
#   - combine clinically different external endpoints into one meta-analysis
#
# It DOES:
#   1. integrate the locked numerical evidence from Scripts 126, 135, 136,
#      136b, 137, 141, 142b and 143;
#   2. classify claims as:
#        DISCOVERY
#        INTERNAL CHARACTERIZATION
#        EXTERNAL REPLICATION
#        EXTERNAL SEVERITY VALIDATION
#        SECONDARY / EXPLORATORY
#   3. define manuscript-safe wording and prohibited overclaims;
#   4. generate a master manuscript key-number table;
#   5. generate a manuscript-ready Results section;
#   6. generate Methods analysis/provenance map;
#   7. generate Abstract-ready numerical results;
#   8. generate final figure hierarchy;
#   9. generate a reviewer-risk / overclaim audit;
#  10. export a master Excel workbook.
#
#
# FROZEN FIVE-GENE PANEL
# ----------------------
# UP:
#   CD177
#   HK3
#   IRAK3
#
# DOWN:
#   CARD11
#   IKZF2
#
# SCORE:
#   mean[z(CD177), z(HK3), z(IRAK3)] -
#   mean[z(CARD11), z(IKZF2)]
#
#
# CENTRAL MANUSCRIPT INTERPRETATION
# ---------------------------------
# The five-gene signature is interpreted as a compact molecular readout of
# a myeloid-adaptive host-response axis associated with transcriptomic
# endotypes and organ-dysfunction severity.
#
# It is NOT interpreted as a clinically calibrated diagnostic or prognostic
# assay.
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

script_name <- "144_final_manuscript_audit_package.R"

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
cat("Running Script 144\n")
cat("Final manuscript audit and publication package\n")
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
  "144_final_manuscript_audit"
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


for (
  dir_path in c(
    output_dir,
    tables_dir,
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


five_gene_formula <- paste0(
  "mean[z(CD177), z(HK3), z(IRAK3)] - ",
  "mean[z(CARD11), z(IKZF2)]"
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


fmt_p <- function(
    x,
    digits = 3
) {
  
  if (
    length(x) == 0 ||
    is.na(x)
  ) {
    
    return(
      "NA"
    )
  }
  
  
  if (
    x < 0.001
  ) {
    
    return(
      format(
        x,
        scientific = TRUE,
        digits = digits
      )
    )
  }
  
  
  formatC(
    x,
    format = "f",
    digits = digits
  )
}


fmt_num <- function(
    x,
    digits = 3
) {
  
  if (
    length(x) == 0 ||
    is.na(x)
  ) {
    
    return(
      "NA"
    )
  }
  
  
  formatC(
    x,
    format = "f",
    digits = digits
  )
}


safe_one_row <- function(
    df,
    filter_expression,
    label
) {
  
  out <- df %>%
    dplyr::filter(
      {{ filter_expression }}
    )
  
  
  if (
    nrow(
      out
    ) != 1
  ) {
    
    stop(
      paste0(
        "Expected exactly one row for: ",
        label,
        ". Found ",
        nrow(
          out
        ),
        "."
      )
    )
  }
  
  
  out
}


write_text <- function(
    lines,
    filename
) {
  
  path <- file.path(
    text_dir,
    filename
  )
  
  
  writeLines(
    lines,
    con = path,
    useBytes = TRUE
  )
  
  
  return(
    path
  )
}


# ==============================================================================
# 5. REQUIRED SCRIPT 143 OUTPUTS
# ==============================================================================

dir_143 <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "143_multicohort_integration"
)


if (!dir.exists(
  dir_143
)) {
  
  stop(
    paste0(
      "Script 143 output directory not found:\n",
      dir_143
    )
  )
}


file_143_key_numbers <- file.path(
  dir_143,
  "tables",
  "143_manuscript_key_numbers.csv"
)


file_143_evidence <- file.path(
  dir_143,
  "tables",
  "143_multicohort_evidence_summary.csv"
)


file_143_direction <- file.path(
  dir_143,
  "tables",
  "143_cross_cohort_gene_direction_matrix.csv"
)


file_143_endotypes <- file.path(
  dir_143,
  "tables",
  "143_discovery_integrated_endotype_summary.csv"
)


file_143_benchmark <- file.path(
  dir_143,
  "tables",
  "143_published_signature_benchmark_summary.csv"
)


file_143_gse154918 <- file.path(
  dir_143,
  "tables",
  "143_GSE154918_group_summary.csv"
)


file_143_location_meta <- file.path(
  dir_143,
  "tables",
  "143_GSE185263_location_random_effects_meta.csv"
)


file_143_external_auc <- file.path(
  dir_143,
  "tables",
  "143_external_binary_endpoint_summary.csv"
)


required_143_files <- c(
  file_143_key_numbers,
  file_143_evidence,
  file_143_direction,
  file_143_endotypes,
  file_143_benchmark,
  file_143_gse154918,
  file_143_location_meta,
  file_143_external_auc
)


missing_143 <- required_143_files[
  !file.exists(
    required_143_files
  )
]


if (
  length(
    missing_143
  ) > 0
) {
  
  stop(
    paste0(
      "Required Script 143 output(s) missing:\n",
      paste(
        missing_143,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# 6. DIRECT SCRIPT 141 / 142b OUTPUTS
# ==============================================================================

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


file_141_comparisons_candidates <- list.files(
  dir_141,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)


file_141_comparisons <- NA_character_


for (
  candidate in file_141_comparisons_candidates
) {
  
  header <- tryCatch(
    
    data.table::fread(
      candidate,
      nrows = 0,
      data.table = FALSE,
      check.names = FALSE
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  
  if (
    !is.null(
      header
    ) &&
    all(
      c(
        "comparison",
        "auc_fixed_direction"
      ) %in%
      names(
        header
      )
    )
  ) {
    
    file_141_comparisons <- candidate
    
    break
  }
}


if (
  is.na(
    file_141_comparisons
  )
) {
  
  stop(
    "Unable to locate Script 141 comparison CSV."
  )
}


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


required_direct_files <- c(
  file_141_comparisons,
  file_142_primary,
  file_142_secondary,
  file_142_gene_sofa,
  file_142_adjusted
)


missing_direct <- required_direct_files[
  !file.exists(
    required_direct_files
  )
]


if (
  length(
    missing_direct
  ) > 0
) {
  
  stop(
    paste0(
      "Required direct upstream output(s) missing:\n",
      paste(
        missing_direct,
        collapse = "\n"
      )
    )
  )
}


cat(
  "Required integrated and direct upstream outputs found.\n\n"
)


# ==============================================================================
# 7. READ TABLES
# ==============================================================================

key_143 <- read_csv_fast(
  file_143_key_numbers
)


evidence_143 <- read_csv_fast(
  file_143_evidence
)


direction_143 <- read_csv_fast(
  file_143_direction
)


endotypes_143 <- read_csv_fast(
  file_143_endotypes
)


benchmark_143 <- read_csv_fast(
  file_143_benchmark
)


gse154918_groups <- read_csv_fast(
  file_143_gse154918
)


location_meta <- read_csv_fast(
  file_143_location_meta
)


external_auc <- read_csv_fast(
  file_143_external_auc
)


comparisons_141 <- read_csv_fast(
  file_141_comparisons
)


primary_142 <- read_csv_fast(
  file_142_primary
)


secondary_142 <- read_csv_fast(
  file_142_secondary
)


gene_sofa_142 <- read_csv_fast(
  file_142_gene_sofa
)


adjusted_142_all <- read_csv_fast(
  file_142_adjusted
)


cat(
  "Upstream numerical evidence loaded successfully.\n\n"
)


# ==============================================================================
# 8. EXTRACT CORE NUMERICAL RESULTS
# ==============================================================================

row_141_primary <- comparisons_141 %>%
  
  dplyr::filter(
    comparison ==
      "Sepsis_or_shock_vs_uncomplicated"
  )


row_141_shock <- comparisons_141 %>%
  
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
    "GSE154918 primary comparison not uniquely identified."
  )
}


if (
  nrow(
    row_141_shock
  ) != 1
) {
  
  stop(
    "GSE154918 Shock vs uncomplicated comparison not uniquely identified."
  )
}


if (
  nrow(
    primary_142
  ) != 1
) {
  
  stop(
    "Script 142b primary table should contain exactly one row."
  )
}


adjusted_142 <- adjusted_142_all %>%
  
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
    "Adjusted SOFA coefficient not uniquely identified."
  )
}


secondary_sofa <- secondary_142 %>%
  
  dplyr::filter(
    analysis ==
      "SOFA_ge2_vs_SOFA_0_1"
  )


secondary_mortality <- secondary_142 %>%
  
  dplyr::filter(
    analysis ==
      "Died_vs_Survived"
  )


secondary_site <- secondary_142 %>%
  
  dplyr::filter(
    analysis ==
      "ICU_vs_Emergency_Room"
  )


secondary_disease <- secondary_142 %>%
  
  dplyr::filter(
    analysis ==
      "Sepsis_vs_healthy_contextual"
  )


if (
  any(
    c(
      nrow(
        secondary_sofa
      ),
      nrow(
        secondary_mortality
      ),
      nrow(
        secondary_site
      ),
      nrow(
        secondary_disease
      )
    ) != 1
  )
) {
  
  stop(
    "One or more Script 142b secondary analyses were not uniquely identified."
  )
}


# ==============================================================================
# 9. LOCKED UPSTREAM PROVENANCE
#
# These values are not recomputed here.
# They come from previously completed scripts.
# ==============================================================================

locked_provenance <- tibble::tibble(
  
  item = c(
    
    "Candidate genes in primary biology-guided pool",
    
    "Candidate combinations searched without forced DCAF17",
    
    "Candidate combinations searched with forced DCAF17",
    
    "Primary frozen five-gene panel",
    
    "Primary score formula",
    
    "Blood BP samples",
    
    "Healthy blood controls",
    
    "Robust core blood DEGs",
    
    "Robust core UP genes",
    
    "Robust core DOWN genes",
    
    "SRS1 BP",
    
    "SRS2 BP",
    
    "CTS1 BP",
    
    "CTS2 BP",
    
    "CTS3 BP"
  ),
  
  value = c(
    
    "13",
    
    "5432",
    
    "2707",
    
    paste(
      five_genes,
      collapse = "; "
    ),
    
    five_gene_formula,
    
    "35",
    
    "10",
    
    "1796",
    
    "1133",
    
    "663",
    
    "28",
    
    "7",
    
    "14",
    
    "6",
    
    "15"
  ),
  
  source = c(
    
    "Script 126",
    
    "Script 126",
    
    "Script 126",
    
    "Script 126",
    
    "Scripts 126/135",
    
    "Discovery cohort",
    
    "Discovery cohort",
    
    "Robust-core blood DE analysis",
    
    "Robust-core blood DE analysis",
    
    "Robust-core blood DE analysis",
    
    "Script 135",
    
    "Script 135",
    
    "Script 135",
    
    "Script 135",
    
    "Script 135"
  )
)


# ==============================================================================
# 10. LOCKED INTERNAL DISCOVERY / CHARACTERIZATION RESULTS
# ==============================================================================

internal_results <- tibble::tibble(
  
  result_id = c(
    "INT01",
    "INT02",
    "INT03",
    "INT04",
    "INT05",
    "INT06",
    "INT07",
    "INT08"
  ),
  
  analysis = c(
    
    "Primary five-gene apparent BP-vs-BC discrimination",
    
    "Primary five-gene repeated five-fold CV",
    
    "Primary score vs SRSq",
    
    "Primary score across CTS",
    
    "Integrated CTS/SRS hierarchy",
    
    "Age/sex-adjusted BP-vs-BC score association",
    
    "Primary score vs CRP",
    
    "SRSq vs CRP"
  ),
  
  estimate = c(
    
    1.000,
    
    1.000,
    
    0.765,
    
    0.661,
    
    0.695,
    
    3.88,
    
    0.574,
    
    0.526
  ),
  
  estimate_type = c(
    
    "Apparent AUC",
    
    "Mean repeated-CV AUC",
    
    "Spearman rho",
    
    "Kruskal-Wallis epsilon2",
    
    "Kruskal-Wallis epsilon2",
    
    "Adjusted score difference",
    
    "Spearman rho",
    
    "Spearman rho"
  ),
  
  p_value = c(
    
    NA_real_,
    
    NA_real_,
    
    8.742e-08,
    
    9.437e-06,
    
    1.91e-05,
    
    3.73e-11,
    
    3.09e-04,
    
    1.17e-03
  ),
  
  adjusted_p = c(
    
    NA_real_,
    
    NA_real_,
    
    NA_real_,
    
    NA_real_,
    
    NA_real_,
    
    NA_real_,
    
    0.0185,
    
    0.0352
  ),
  
  evidence_stage = c(
    
    "Internal characterization",
    
    "Internal resampling",
    
    "Internal characterization",
    
    "Internal characterization",
    
    "Internal characterization",
    
    "Internal sensitivity",
    
    "Exploratory clinical association",
    
    "Exploratory clinical association"
  )
)


# ==============================================================================
# 11. EXTERNAL VALIDATION TABLE
# ==============================================================================

external_results <- tibble::tibble(
  
  result_id = c(
    "EXT01",
    "EXT02",
    "EXT03",
    "EXT04",
    "EXT05",
    "EXT06",
    "EXT07",
    "EXT08",
    "EXT09",
    "EXT10"
  ),
  
  cohort = c(
    
    "GSE154918",
    
    "GSE154918",
    
    "GSE154918",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263 geographic subcohorts"
  ),
  
  analysis = c(
    
    "Frozen five-gene directional concordance",
    
    "Sepsis/septic shock vs uncomplicated infection",
    
    "Septic shock vs uncomplicated infection",
    
    "Primary five-gene score vs continuous SOFA",
    
    "Frozen component genes vs SOFA",
    
    "Age/sex/location-adjusted SOFA coefficient",
    
    "SOFA >=2 vs SOFA 0-1",
    
    "Died vs Survived",
    
    "Sepsis vs healthy",
    
    "Random-effects geographic score-SOFA synthesis"
  ),
  
  estimate = c(
    
    5,
    
    row_141_primary$auc_fixed_direction,
    
    row_141_shock$auc_fixed_direction,
    
    primary_142$rho,
    
    sum(
      gene_sofa_142$direction_concordant,
      na.rm = TRUE
    ),
    
    adjusted_142$estimate,
    
    secondary_sofa$AUC,
    
    secondary_mortality$AUC,
    
    secondary_disease$AUC,
    
    location_meta$random_rho
  ),
  
  estimate_type = c(
    
    "Concordant genes / 5",
    
    "AUC",
    
    "AUC",
    
    "Spearman rho",
    
    "Concordant genes / 5",
    
    "Adjusted beta per SOFA point",
    
    "AUC",
    
    "AUC",
    
    "AUC",
    
    "Random-effects pooled rho"
  ),
  
  CI_low = c(
    
    NA_real_,
    
    row_141_primary$auc_ci_low,
    
    row_141_shock$auc_ci_low,
    
    NA_real_,
    
    NA_real_,
    
    NA_real_,
    
    secondary_sofa$CI_low,
    
    secondary_mortality$CI_low,
    
    secondary_disease$CI_low,
    
    location_meta$random_CI_low
  ),
  
  CI_high = c(
    
    NA_real_,
    
    row_141_primary$auc_ci_high,
    
    row_141_shock$auc_ci_high,
    
    NA_real_,
    
    NA_real_,
    
    NA_real_,
    
    secondary_sofa$CI_high,
    
    secondary_mortality$CI_high,
    
    secondary_disease$CI_high,
    
    location_meta$random_CI_high
  ),
  
  p_value = c(
    
    NA_real_,
    
    row_141_primary$p_value,
    
    row_141_shock$p_value,
    
    primary_142$p_value,
    
    NA_real_,
    
    adjusted_142$p_value,
    
    secondary_sofa$p_value,
    
    secondary_mortality$p_value,
    
    secondary_disease$p_value,
    
    location_meta$random_p
  ),
  
  evidence_stage = c(
    
    "External replication",
    
    "Prespecified external primary",
    
    "External secondary",
    
    "Prespecified external primary",
    
    "External component-level replication",
    
    "External sensitivity",
    
    "External secondary",
    
    "External secondary",
    
    "External contextual",
    
    "External geographic sensitivity"
  )
)


# ==============================================================================
# 12. BENCHMARKING SUMMARY
# ==============================================================================

published_only <- benchmark_143 %>%
  
  dplyr::filter(
    category ==
      "Published comparator"
  )


benchmark_summary <- tibble::tibble(
  
  metric = c(
    
    "Number of published comparator signatures",
    
    "Minimum SRSq rho among published comparators",
    
    "Maximum SRSq rho among published comparators",
    
    "Minimum CTS epsilon2 among published comparators",
    
    "Maximum CTS epsilon2 among published comparators",
    
    "Primary five-gene SRSq rho",
    
    "Primary five-gene CTS epsilon2"
  ),
  
  value = c(
    
    nrow(
      published_only
    ),
    
    min(
      published_only$SRSq_rho
    ),
    
    max(
      published_only$SRSq_rho
    ),
    
    min(
      published_only$CTS_epsilon2
    ),
    
    max(
      published_only$CTS_epsilon2
    ),
    
    benchmark_143$SRSq_rho[
      benchmark_143$signature ==
        "Primary five-gene"
    ],
    
    benchmark_143$CTS_epsilon2[
      benchmark_143$signature ==
        "Primary five-gene"
    ]
  )
)


# ==============================================================================
# 13. MASTER CLAIM AUDIT
# ==============================================================================

claim_audit <- tibble::tibble(
  
  claim_id = c(
    
    "C01",
    "C02",
    "C03",
    "C04",
    "C05",
    "C06",
    "C07",
    "C08",
    "C09",
    "C10",
    "C11",
    "C12",
    "C13",
    "C14",
    "C15"
  ),
  
  manuscript_domain = c(
    
    "Blood transcriptomic response",
    
    "Panel derivation",
    
    "Internal discrimination",
    
    "SRS",
    
    "CTS",
    
    "Integrated endotypes",
    
    "Clinical association",
    
    "Published signatures",
    
    "GSE154918",
    
    "GSE154918",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "GSE185263",
    
    "Overall interpretation"
  ),
  
  evidence_stage = c(
    
    "Discovery",
    
    "Discovery / derivation",
    
    "Internal characterization",
    
    "Internal characterization",
    
    "Internal characterization",
    
    "Internal characterization",
    
    "Exploratory clinical",
    
    "Internal benchmarking",
    
    "External replication",
    
    "External secondary",
    
    "Prespecified external primary",
    
    "External component replication",
    
    "External sensitivity",
    
    "External geographic sensitivity",
    
    "Integrated multicohort"
  ),
  
  supported_claim = c(
    
    paste0(
      "Sepsis blood shows a robust transcriptional program characterized ",
      "by concurrent myeloid activation and adaptive immune suppression."
    ),
    
    paste0(
      "The five-gene panel was obtained from a biologically constrained ",
      "candidate pool representing the dominant myeloid and adaptive ",
      "components of the robust blood transcriptional response."
    ),
    
    paste0(
      "The five-gene score showed complete apparent and internally ",
      "resampled separation of sepsis from healthy controls in the ",
      "discovery dataset."
    ),
    
    paste0(
      "The five-gene score closely tracks the SRS immune-dysfunction axis."
    ),
    
    paste0(
      "The five-gene score differs substantially across consensus ",
      "transcriptomic subtypes."
    ),
    
    paste0(
      "Integrated CTS/SRS groups define an ordered host-response hierarchy."
    ),
    
    paste0(
      "The score is associated with systemic inflammatory activity, ",
      "including CRP, in exploratory clinical analyses."
    ),
    
    paste0(
      "Previously published transcriptomic signatures converge on a ",
      "similar host-response axis in the discovery cohort."
    ),
    
    paste0(
      "In GSE154918, all five frozen genes reproduced the prespecified ",
      "direction and the score increased across clinical states toward ",
      "septic shock."
    ),
    
    paste0(
      "The prespecified sepsis/septic-shock versus uncomplicated-infection ",
      "comparison showed only modest discrimination."
    ),
    
    paste0(
      "In GSE185263, the frozen five-gene score independently tracked ",
      "continuous organ-dysfunction severity measured by SOFA."
    ),
    
    paste0(
      "All five component genes independently reproduced their expected ",
      "direction of association with SOFA."
    ),
    
    paste0(
      "The score-SOFA association persisted after adjustment for age, ",
      "sex and geographic collection location."
    ),
    
    paste0(
      "The score-SOFA association was directionally consistent across ",
      "five geographic subcohorts, with no detectable between-location ",
      "heterogeneity."
    ),
    
    paste0(
      "Across discovery and external cohorts, the five-gene score is best ",
      "interpreted as a compact molecular readout of a reproducible ",
      "myeloid-adaptive host-response axis."
    )
  ),
  
  manuscript_safe = c(
    
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),
  
  prohibited_overclaim = c(
    
    "Do not claim that the robust core is universal to all sepsis populations.",
    
    "Do not call the panel endotype-selected or SRS/CTS-selected.",
    
    paste0(
      "Do not call AUC=1 independent validation or clinical diagnostic ",
      "validation."
    ),
    
    "Do not claim that the score is identical to the published SRS classifier.",
    
    "Do not claim CTS causality or treatment-response prediction.",
    
    "Do not call the four integrated groups a continuous variable.",
    
    paste0(
      "Do not claim independent prognostic validation from exploratory ",
      "clinical associations."
    ),
    
    paste0(
      "Do not describe published comparator calculations in the same ",
      "dataset as independent validation."
    ),
    
    paste0(
      "Do not claim that GSE154918 validates a calibrated diagnostic ",
      "threshold."
    ),
    
    paste0(
      "Do not hide or reinterpret the nonsignificant prespecified primary ",
      "comparison."
    ),
    
    paste0(
      "Do not describe rho=0.311 as strong diagnostic discrimination; ",
      "it is an independent severity association."
    ),
    
    "Do not claim mechanistic causality for individual genes.",
    
    paste0(
      "Do not claim that demographic/geographic adjustment removes all ",
      "possible confounding."
    ),
    
    paste0(
      "Do not claim five geographic locations are five independent ",
      "external validation cohorts."
    ),
    
    paste0(
      "Do not claim clinical assay validation, generalizable diagnostic ",
      "accuracy, or validated prognosis."
    )
  )
)


# ==============================================================================
# 14. OVERCLAIM / REVIEWER-RISK AUDIT
# ==============================================================================

reviewer_risk <- tibble::tibble(
  
  risk_id = paste0(
    "R",
    sprintf(
      "%02d",
      1:12
    )
  ),
  
  topic = c(
    
    "Same-cohort panel derivation",
    
    "Apparent AUC=1",
    
    "Small healthy-control group",
    
    "Healthy controls rather than infected controls",
    
    "Batch/chip structure",
    
    "Age/sex imbalance",
    
    "SRS/CTS circularity",
    
    "External score standardization",
    
    "GSE154918 primary result",
    
    "Mortality association",
    
    "Geographic meta-analysis",
    
    "Clinical translation"
  ),
  
  risk = c(
    
    paste0(
      "The panel was selected using the discovery cohort; apparent and ",
      "internal-CV discrimination may be optimistic."
    ),
    
    paste0(
      "AUC=1 can be overinterpreted as diagnostic perfection."
    ),
    
    "Only 10 healthy blood controls were available in the discovery cohort.",
    
    paste0(
      "BP-vs-BC discrimination does not establish specificity against ",
      "non-septic infection or sterile inflammation."
    ),
    
    paste0(
      "Condition and sequencing-chip structure may not be fully separable ",
      "in a single multivariable model."
    ),
    
    paste0(
      "Discovery sepsis and healthy groups differed in age and sex."
    ),
    
    paste0(
      "SRS and CTS were evaluated in the same discovery expression dataset."
    ),
    
    paste0(
      "External scores were standardized within external cohorts rather ",
      "than using a clinically locked absolute calibration."
    ),
    
    paste0(
      "The prespecified GSE154918 sepsis/shock vs uncomplicated-infection ",
      "comparison was not statistically significant."
    ),
    
    paste0(
      "Mortality analysis in GSE185263 was secondary and unadjusted for ",
      "full clinical prognostic covariates."
    ),
    
    paste0(
      "Five geographic strata are subcohorts of one study, not five ",
      "independent studies."
    ),
    
    paste0(
      "Transcriptomic replication does not equal validation of a ",
      "laboratory qPCR/ddPCR clinical assay."
    )
  ),
  
  mitigation = c(
    
    paste0(
      "Describe panel as biology-guided candidate signature; emphasize ",
      "external frozen-panel replication."
    ),
    
    paste0(
      "Report as apparent/internal resampled performance and place ",
      "diagnostic ROC emphasis in supplementary material."
    ),
    
    "Report exact sample numbers and confidence intervals.",
    
    paste0(
      "Use GSE154918 uncomplicated infection as the more clinically ",
      "relevant specificity challenge."
    ),
    
    paste0(
      "Retain prior chip-adjusted DE sensitivity as complementary evidence; ",
      "avoid unstable fully adjusted models."
    ),
    
    paste0(
      "Report age/sex-adjusted score effect from Script 136b."
    ),
    
    paste0(
      "State explicitly that SRS/CTS labels were NOT used for panel ",
      "feature selection."
    ),
    
    paste0(
      "Describe external evidence as transcriptomic replication rather ",
      "than calibration validation."
    ),
    
    paste0(
      "Report the result transparently: AUC=0.656, p=0.107; emphasize ",
      "directional replication and shock/severity gradient."
    ),
    
    paste0(
      "Describe as supportive secondary association only."
    ),
    
    paste0(
      "Use wording 'geographic sensitivity analysis within GSE185263'."
    ),
    
    paste0(
      "Require future prospective qPCR/ddPCR analytical and clinical ",
      "validation before clinical-use claims."
    )
  )
)


# ==============================================================================
# 15. FINAL FIGURE HIERARCHY
# ==============================================================================

figure_plan <- tibble::tibble(
  
  figure = c(
    
    "Figure 1",
    
    "Figure 2",
    
    "Figure 3",
    
    "Figure 4",
    
    "Figure 5",
    
    "Supplementary Figure S1",
    
    "Supplementary Figure S2",
    
    "Supplementary Figure S3",
    
    "Supplementary Figure S4"
  ),
  
  primary_message = c(
    
    paste0(
      "Sepsis blood is characterized by coordinated myeloid activation ",
      "and adaptive immune suppression."
    ),
    
    paste0(
      "SRS and CTS reveal hierarchical host-response organization and ",
      "the five-gene score captures this axis."
    ),
    
    paste0(
      "The five-gene score aligns with previously published sepsis ",
      "transcriptomic signatures and selected clinical inflammatory features."
    ),
    
    paste0(
      "The frozen five-gene signature reproduces gene direction and a ",
      "clinical-state gradient in GSE154918."
    ),
    
    paste0(
      "The frozen five-gene score independently tracks SOFA severity in ",
      "GSE185263 and is geographically reproducible."
    ),
    
    "Discovery demographic and batch sensitivity analyses.",
    
    "Internal BP-vs-BC ROC and repeated cross-validation.",
    
    "Full comparator-signature benchmarking.",
    
    "External binary AUC summary."
  ),
  
  recommended_panels = c(
    
    paste0(
      "Volcano / robust-core summary; representative genes; pathway or ",
      "functional summary."
    ),
    
    paste0(
      "SRS; CTS; integrated CTS/SRS hierarchy; five-gene heatmap or ",
      "direction schematic."
    ),
    
    paste0(
      "Published-signature convergence; CRP/clinical association."
    ),
    
    paste0(
      "GSE154918 clinical-state boxplot; five-gene direction audit."
    ),
    
    paste0(
      "GSE185263 score-SOFA plot; component gene-SOFA effects; ",
      "geographic forest plot."
    ),
    
    "Age/sex sensitivity; chip/batch sensitivity.",
    
    "Internal ROC / resampling distributions.",
    
    "All published comparator analyses.",
    
    "Script 143 external binary endpoint forest plot."
  ),
  
  source = c(
    
    "Upstream blood DE / enrichment scripts",
    
    "Scripts 135 and 143",
    
    "Scripts 136/137/143",
    
    "Scripts 141/143",
    
    "Scripts 142b/143",
    
    "Scripts 136b and prior blood sensitivity analyses",
    
    "Script 135",
    
    "Script 137",
    
    "Script 143"
  )
)


# ==============================================================================
# 16. FINAL MANUSCRIPT SECTION MAP
# ==============================================================================

results_structure <- tibble::tibble(
  
  section = c(
    
    "3.1",
    
    "3.2",
    
    "3.3",
    
    "3.4",
    
    "3.5",
    
    "3.6",
    
    "3.7",
    
    "3.8",
    
    "3.9"
  ),
  
  title = c(
    
    "Blood transcriptomic remodeling in sepsis",
    
    paste0(
      "A robust blood transcriptional signature characterizes ",
      "sepsis-associated host-response dysregulation"
    ),
    
    paste0(
      "Established transcriptomic endotypes reveal hierarchical ",
      "host-response organization"
    ),
    
    paste0(
      "A biology-guided five-gene signature captures the ",
      "myeloid-adaptive host-response axis"
    ),
    
    paste0(
      "The five-gene score recapitulates SRS and consensus ",
      "transcriptomic subtype states"
    ),
    
    paste0(
      "Clinical associations and previously published signatures support ",
      "a shared host-response continuum"
    ),
    
    paste0(
      "Independent GSE154918 replication preserves gene direction and ",
      "reveals a clinical-state gradient"
    ),
    
    paste0(
      "The frozen five-gene score independently tracks organ-dysfunction ",
      "severity in GSE185263"
    ),
    
    paste0(
      "The score-SOFA relationship is reproducible across geographic ",
      "subcohorts"
    )
  ),
  
  evidence_class = c(
    
    "Discovery",
    
    "Discovery",
    
    "Internal characterization",
    
    "Derivation",
    
    "Internal characterization",
    
    "Exploratory / benchmarking",
    
    "External replication",
    
    "Prespecified external validation",
    
    "External geographic sensitivity"
  )
)


# ==============================================================================
# 17. METHODS ANALYSIS / PROVENANCE MAP
# ==============================================================================

methods_map <- tibble::tibble(
  
  methods_subsection = c(
    
    "Blood differential expression",
    
    "Robust-core definition",
    
    "Transcriptomic endotype assignment",
    
    "Candidate-gene pool definition",
    
    "Compact-panel search",
    
    "Five-gene score calculation",
    
    "Internal discrimination",
    
    "Clinical association analyses",
    
    "Demographic sensitivity",
    
    "Published-signature benchmarking",
    
    "GSE154918 external replication",
    
    "GSE185263 external severity validation",
    
    "Geographic sensitivity analysis"
  ),
  
  data_used = c(
    
    "BP35 vs BC10",
    
    "Simple and chip-adjusted BP-vs-BC DE results",
    
    "BP35",
    
    "Robust blood host-response genes",
    
    "Discovery BP/BC expression",
    
    "Discovery and external normalized expression",
    
    "BP35 vs BC10",
    
    "BP35",
    
    "BP34/35 plus BC10 depending on variable completeness",
    
    "Discovery blood expression",
    
    "GSE154918 baseline whole-blood RNA-seq",
    
    "GSE185263 whole-blood RNA-seq",
    
    "Five GSE185263 geographic strata"
  ),
  
  principal_method = c(
    
    "DESeq2",
    
    "Intersection/concordance of primary DE models",
    
    "Published SRS framework and BP-only CTS framework",
    
    "Biology-guided myeloid/adaptive candidate definition",
    
    "Exhaustive compact-panel combination search",
    
    five_gene_formula,
    
    "Apparent ROC + repeated 5-fold CV",
    
    "Spearman / Wilcoxon / Kruskal-Wallis with BH correction",
    
    "Linear models with age and sex",
    
    "Predefined published transcriptomic score implementations",
    
    paste0(
      "Frozen genes and frozen score; primary Seps_P+Shock_P vs Inf1_P; ",
      "no refitting"
    ),
    
    paste0(
      "Frozen score; primary Spearman correlation with continuous ",
      "24-h SOFA"
    ),
    
    paste0(
      "Location-specific Spearman correlations followed by Fisher-z ",
      "fixed/random-effects synthesis"
    )
  ),
  
  leakage_control = c(
    
    "Not applicable",
    
    "Not applicable",
    
    "Endotype labels not used to select five-gene panel",
    
    "Defined before SRS/CTS validation of panel",
    
    "Discovery-stage optimization; acknowledged as such",
    
    "Gene identities and directions frozen before external testing",
    
    "Internal only; not treated as independent validation",
    
    "Exploratory family clearly labeled",
    
    "Sensitivity analysis only",
    
    "No new feature selection",
    
    "No feature selection, refitting, threshold optimization, or direction flipping",
    
    paste0(
      "Primary endpoint declared after metadata audit and before ",
      "five-gene expression analysis"
    ),
    
    "Within-cohort sensitivity only; not treated as five independent studies"
  )
)


# ==============================================================================
# 18. ABSTRACT KEY-NUMBER SET
# ==============================================================================

abstract_numbers <- tibble::tibble(
  
  order = 1:8,
  
  result = c(
    
    "Robust core blood DEGs",
    
    "Five-gene score vs SRSq",
    
    "CTS effect size",
    
    "GSE154918 directional concordance",
    
    "GSE154918 primary AUC",
    
    "GSE185263 primary score-SOFA association",
    
    "GSE185263 adjusted SOFA coefficient",
    
    "GSE185263 geographic random-effects association"
  ),
  
  value = c(
    
    "1,796",
    
    paste0(
      "rho=",
      fmt_num(
        0.765
      ),
      "; p=",
      fmt_p(
        8.742e-08
      )
    ),
    
    paste0(
      "epsilon²=",
      fmt_num(
        0.661
      ),
      "; p=",
      fmt_p(
        9.437e-06
      )
    ),
    
    "5/5 genes",
    
    paste0(
      "AUC=",
      fmt_num(
        row_141_primary$auc_fixed_direction
      ),
      "; 95% CI ",
      fmt_num(
        row_141_primary$auc_ci_low
      ),
      "-",
      fmt_num(
        row_141_primary$auc_ci_high
      ),
      "; p=",
      fmt_p(
        row_141_primary$p_value
      )
    ),
    
    paste0(
      "rho=",
      fmt_num(
        primary_142$rho
      ),
      "; p=",
      fmt_p(
        primary_142$p_value
      ),
      "; n=",
      primary_142$n
    ),
    
    paste0(
      "beta=",
      fmt_num(
        adjusted_142$estimate
      ),
      "; p=",
      fmt_p(
        adjusted_142$p_value
      )
    ),
    
    paste0(
      "rho=",
      fmt_num(
        location_meta$random_rho
      ),
      "; 95% CI ",
      fmt_num(
        location_meta$random_CI_low
      ),
      "-",
      fmt_num(
        location_meta$random_CI_high
      ),
      "; I²=",
      fmt_num(
        location_meta$I2_percent,
        digits = 1
      ),
      "%"
    )
  ),
  
  priority = c(
    
    "HIGH",
    
    "HIGH",
    
    "HIGH",
    
    "HIGH",
    
    "REPORT HONESTLY",
    
    "VERY HIGH",
    
    "HIGH",
    
    "HIGH"
  )
)


# ==============================================================================
# 19. MANUSCRIPT-READY RESULTS SECTION - ENGLISH
# ==============================================================================

results_en <- c(
  
  "3. RESULTS",
  
  "",
  
  "3.1. Blood transcriptomic remodeling in sepsis",
  
  paste0(
    "Whole-blood transcriptomic profiling identified extensive ",
    "sepsis-associated transcriptional remodeling. The simple BP-versus-BC ",
    "model identified 2,659 differentially expressed genes, whereas the ",
    "chip-adjusted model identified 4,125 genes. The two models showed ",
    "strong concordance of gene-level effects, and 1,796 genes formed a ",
    "robust core signature, including 1,133 genes upregulated and 663 genes ",
    "downregulated in sepsis. Representative upregulated genes included ",
    "CD177, IL1R2, MMP9, S100A12, HK3 and IRAK3, whereas CARD11, IKZF2, ",
    "P2RY10 and ST6GAL1 were among the downregulated genes. Collectively, ",
    "these changes were consistent with enhanced myeloid/neutrophil ",
    "activation accompanied by suppression of adaptive lymphoid programs."
  ),
  
  "",
  
  paste0(
    "3.2. A robust blood transcriptional signature characterizes ",
    "sepsis-associated host-response dysregulation"
  ),
  
  paste0(
    "The robust blood signature therefore defined two opposing components ",
    "of the systemic host response: a sepsis-upregulated myeloid/innate ",
    "program and a sepsis-downregulated adaptive/lymphoid program. This ",
    "dual organization was used as the biological basis for subsequent ",
    "endotype characterization and compact biomarker development."
  ),
  
  "",
  
  paste0(
    "3.3. Established transcriptomic endotypes reveal hierarchical ",
    "host-response organization"
  ),
  
  paste0(
    "Among 35 patients with sepsis, SRS classification assigned 28 ",
    "patients to SRS1 and seven to SRS2. Consensus transcriptomic ",
    "subtyping assigned 14 patients to CTS1, six to CTS2 and 15 to CTS3. ",
    "The two frameworks were related but non-equivalent: CTS1 and CTS2 ",
    "occurred exclusively within SRS1, whereas CTS3 contained both SRS1 ",
    "and SRS2 samples. This pattern indicated a hierarchical organization ",
    "of the blood host response rather than a single binary endotype axis."
  ),
  
  "",
  
  paste0(
    "3.4. A biology-guided five-gene signature captures the ",
    "myeloid-adaptive host-response axis"
  ),
  
  paste0(
    "Candidate genes were not selected from the full transcriptome solely ",
    "according to statistical rank. A biologically constrained candidate ",
    "pool was defined to represent the two dominant components of the ",
    "robust blood response. The myeloid-associated pool comprised CD177, ",
    "HK3, IRAK3, PFKFB3, S100A12 and MMP9, whereas the adaptive-associated ",
    "pool comprised CARD11, IKZF2, NR1D2, P2RY10, RPS6 and ST6GAL1; ",
    "DCAF17 was additionally retained as a leading single-gene candidate. ",
    "An exhaustive search evaluated 5,432 permitted compact blood-panel ",
    "combinations, together with 2,707 combinations in which DCAF17 was ",
    "forced into the panel. The primary five-gene signature comprised ",
    "CD177, HK3 and IRAK3 as the upregulated component and CARD11 and IKZF2 ",
    "as the downregulated component. SRS and CTS assignments were not used ",
    "as feature-selection criteria."
  ),
  
  "",
  
  paste0(
    "The composite score was calculated as the mean standardized ",
    "expression of CD177, HK3 and IRAK3 minus the mean standardized ",
    "expression of CARD11 and IKZF2. In the discovery cohort, the score ",
    "showed complete apparent separation of sepsis and healthy controls. ",
    "Because the panel was derived in the same cohort, this result was ",
    "treated as internal characterization rather than independent ",
    "diagnostic validation."
  ),
  
  "",
  
  paste0(
    "3.5. The five-gene score recapitulates SRS and consensus ",
    "transcriptomic subtype states"
  ),
  
  paste0(
    "The five-gene score was strongly associated with the quantitative ",
    "SRS axis (Spearman rho = 0.765, p = 8.74 × 10^-8) and differed ",
    "substantially across CTS classes (Kruskal-Wallis p = 9.44 × 10^-6; ",
    "epsilon² = 0.661). Integrated CTS/SRS groups formed an ordered ",
    "host-response hierarchy, progressing from CTS1/SRS1 through ",
    "CTS2/SRS1 and CTS3/SRS1 to CTS3/SRS2. The corresponding integrated ",
    "group effect size was epsilon² = 0.695. Thus, the compact score ",
    "recapitulated independently derived transcriptomic states despite ",
    "SRS and CTS labels not being used during panel selection."
  ),
  
  "",
  
  paste0(
    "3.6. Clinical associations and previously published signatures ",
    "support a shared host-response continuum"
  ),
  
  paste0(
    "Exploratory clinical analyses identified two associations that ",
    "remained significant after global false-discovery-rate correction: ",
    "the five-gene score correlated with C-reactive protein ",
    "(rho = 0.574, p = 3.09 × 10^-4, BH-adjusted p = 0.0185), and the ",
    "quantitative SRS axis correlated with C-reactive protein ",
    "(rho = 0.526, p = 0.00117, BH-adjusted p = 0.0352). Associations with ",
    "creatinine, lactate and several other clinical variables were nominal ",
    "and were therefore treated as exploratory. Because the sepsis and ",
    "healthy-control groups differed in age and sex distribution, ",
    "demographic sensitivity analyses were performed. The five-gene score ",
    "remained strongly associated with sepsis after simultaneous adjustment ",
    "for age and sex (adjusted difference 3.88 score units; SE 0.43; ",
    "p = 3.73 × 10^-11)."
  ),
  
  "",
  
  paste0(
    "Five previously published transcriptomic comparators were evaluated ",
    "together with the two current-study five-gene implementations. ",
    "Published signatures showed SRSq correlations ranging from ",
    fmt_num(
      min(
        published_only$SRSq_rho
      )
    ),
    " to ",
    fmt_num(
      max(
        published_only$SRSq_rho
      )
    ),
    " and CTS effect sizes ranging from ",
    fmt_num(
      min(
        published_only$CTS_epsilon2
      )
    ),
    " to ",
    fmt_num(
      max(
        published_only$CTS_epsilon2
      )
    ),
    ". The primary five-gene score showed an SRSq correlation of 0.765 ",
    "and a CTS effect size of 0.661, supporting convergence of multiple ",
    "transcriptomic signatures on a related host-response continuum."
  ),
  
  "",
  
  paste0(
    "3.7. Independent GSE154918 replication preserves gene direction ",
    "and reveals a clinical-state gradient"
  ),
  
  paste0(
    "The frozen panel was next evaluated without feature reselection, ",
    "coefficient refitting, threshold optimization or post hoc direction ",
    "reversal in the independent GSE154918 whole-blood RNA-seq cohort. ",
    "All five component genes reproduced their prespecified direction of ",
    "change. Across baseline clinical groups, the median five-gene score ",
    "increased from healthy participants to uncomplicated infection, sepsis ",
    "and septic shock. The prespecified primary comparison of sepsis or ",
    "septic shock with uncomplicated infection showed modest discrimination ",
    "(AUC = ",
    fmt_num(
      row_141_primary$auc_fixed_direction
    ),
    ", 95% CI ",
    fmt_num(
      row_141_primary$auc_ci_low
    ),
    "-",
    fmt_num(
      row_141_primary$auc_ci_high
    ),
    "; p = ",
    fmt_p(
      row_141_primary$p_value
    ),
    "). In contrast, septic shock was more clearly separated from ",
    "uncomplicated infection (AUC = ",
    fmt_num(
      row_141_shock$auc_fixed_direction
    ),
    "; p = ",
    fmt_p(
      row_141_shock$p_value
    ),
    "). These findings supported directional and clinical-state ",
    "replication but did not establish a calibrated diagnostic threshold."
  ),
  
  "",
  
  paste0(
    "3.8. The frozen five-gene score independently tracks ",
    "organ-dysfunction severity in GSE185263"
  ),
  
  paste0(
    "The larger independent GSE185263 cohort comprised 348 sepsis samples ",
    "and 44 healthy controls, with SOFA available for 345 patients with ",
    "sepsis. The prespecified primary external endpoint was met: the frozen ",
    "five-gene score correlated positively with continuous 24-hour SOFA ",
    "(Spearman rho = ",
    fmt_num(
      primary_142$rho
    ),
    ", p = ",
    fmt_p(
      primary_142$p_value
    ),
    "; n = ",
    primary_142$n,
    "). All five component genes reproduced their expected direction of ",
    "association with SOFA, and each remained significant after gene-level ",
    "false-discovery-rate correction."
  ),
  
  "",
  
  paste0(
    "The association persisted after adjustment for age, sex and geographic ",
    "collection location (beta = ",
    fmt_num(
      adjusted_142$estimate
    ),
    " score units per SOFA point; SE = ",
    fmt_num(
      adjusted_142$SE
    ),
    "; p = ",
    fmt_p(
      adjusted_142$p_value
    ),
    "). In secondary analyses, patients with SOFA >=2 had higher scores ",
    "than those with SOFA 0-1 (AUC = ",
    fmt_num(
      secondary_sofa$AUC
    ),
    ", 95% CI ",
    fmt_num(
      secondary_sofa$CI_low
    ),
    "-",
    fmt_num(
      secondary_sofa$CI_high
    ),
    "; p = ",
    fmt_p(
      secondary_sofa$p_value
    ),
    "). Nonsurvivors also showed higher scores than survivors ",
    "(AUC = ",
    fmt_num(
      secondary_mortality$AUC
    ),
    "; p = ",
    fmt_p(
      secondary_mortality$p_value
    ),
    "), although mortality was a secondary analysis and was not interpreted ",
    "as independent prognostic validation."
  ),
  
  "",
  
  paste0(
    "3.9. The score-SOFA relationship is reproducible across ",
    "geographic subcohorts"
  ),
  
  paste0(
    "The direction of the score-SOFA relationship was positive in all five ",
    "geographic subcohorts represented in GSE185263. Random-effects ",
    "Fisher-z synthesis yielded a pooled Spearman rho of ",
    fmt_num(
      location_meta$random_rho
    ),
    " (95% CI ",
    fmt_num(
      location_meta$random_CI_low
    ),
    "-",
    fmt_num(
      location_meta$random_CI_high
    ),
    "; p = ",
    fmt_p(
      location_meta$random_p
    ),
    "). There was no detectable between-location heterogeneity ",
    "(Q = ",
    fmt_num(
      location_meta$Q
    ),
    ", df = ",
    location_meta$Q_df,
    ", heterogeneity p = ",
    fmt_p(
      location_meta$Q_p
    ),
    "; I² = ",
    fmt_num(
      location_meta$I2_percent,
      digits = 1
    ),
    "%; tau² = ",
    fmt_num(
      location_meta$tau2_fisher_z
    ),
    "). Because these geographic strata originated from a single public ",
    "study, this analysis was interpreted as a geographic sensitivity ",
    "analysis rather than as five independent external validation cohorts."
  )
)


# ==============================================================================
# 20. RUSSIAN RESULTS SUMMARY
# ==============================================================================

results_ru <- c(
  
  "КЛЮЧЕВАЯ ЛОГИКА РАЗДЕЛА RESULTS",
  
  "",
  
  paste0(
    "1. В крови при сепсисе выявлено устойчивое транскриптомное ",
    "перестроение: robust core составил 1 796 генов, включая 1 133 ",
    "повышенных и 663 сниженных."
  ),
  
  paste0(
    "2. Основная биологическая ось включает усиление myeloid/neutrophil ",
    "программ и подавление adaptive/lymphoid программ."
  ),
  
  paste0(
    "3. Пятигенная панель CD177, HK3, IRAK3, CARD11 и IKZF2 была ",
    "получена не из SRS/CTS-меток, а из биологически ограниченного ",
    "набора кандидатов."
  ),
  
  paste0(
    "4. Пятигенный score тесно связан с SRSq ",
    "(rho = 0.765; p = 8.74e-08) и CTS ",
    "(epsilon² = 0.661; p = 9.44e-06)."
  ),
  
  paste0(
    "5. В GSE154918 все 5/5 генов воспроизвели ожидаемое направление. ",
    "Заранее заданное сравнение sepsis/shock против uncomplicated ",
    "infection показало только умеренное разделение ",
    "(AUC = ",
    fmt_num(
      row_141_primary$auc_fixed_direction
    ),
    "; p = ",
    fmt_p(
      row_141_primary$p_value
    ),
    "), однако septic shock отличался существенно лучше ",
    "(AUC = ",
    fmt_num(
      row_141_shock$auc_fixed_direction
    ),
    "; p = ",
    fmt_p(
      row_141_shock$p_value
    ),
    ")."
  ),
  
  paste0(
    "6. В GSE185263 основная заранее заданная внешняя конечная точка ",
    "подтверждена: score коррелировал с SOFA ",
    "(rho = ",
    fmt_num(
      primary_142$rho
    ),
    "; p = ",
    fmt_p(
      primary_142$p_value
    ),
    "; n = ",
    primary_142$n,
    ")."
  ),
  
  paste0(
    "7. Все 5/5 генов в GSE185263 показали ожидаемое направление связи ",
    "с SOFA и сохранили значимость после FDR-коррекции."
  ),
  
  paste0(
    "8. После поправки на возраст, пол и географическое место набора ",
    "связь с SOFA сохранялась: beta = ",
    fmt_num(
      adjusted_142$estimate
    ),
    "; p = ",
    fmt_p(
      adjusted_142$p_value
    ),
    "."
  ),
  
  paste0(
    "9. Random-effects анализ пяти географических подкогорт GSE185263 ",
    "дал pooled rho = ",
    fmt_num(
      location_meta$random_rho
    ),
    " (95% ДИ ",
    fmt_num(
      location_meta$random_CI_low
    ),
    "-",
    fmt_num(
      location_meta$random_CI_high
    ),
    "), I² = ",
    fmt_num(
      location_meta$I2_percent,
      digits = 1
    ),
    "%."
  ),
  
  "",
  
  paste0(
    "Итоговая интерпретация: панель следует позиционировать как ",
    "компактный молекулярный readout воспроизводимой myeloid-adaptive ",
    "host-response axis, связанной с транскриптомными эндотипами и ",
    "тяжестью органной дисфункции, а не как уже валидированный ",
    "клинический диагностический или прогностический тест."
  )
)


# ==============================================================================
# 21. ABSTRACT-READY RESULTS PARAGRAPH
# ==============================================================================

abstract_results_en <- c(
  
  paste0(
    "Blood transcriptomic analysis identified a robust core of 1,796 ",
    "sepsis-associated genes characterized by increased myeloid/neutrophil ",
    "activity and reduced adaptive-immune programs. A biology-guided ",
    "five-gene score comprising CD177, HK3 and IRAK3 as the upregulated ",
    "component and CARD11 and IKZF2 as the downregulated component was ",
    "strongly associated with the SRS axis (rho = 0.765, p = 8.74 × 10^-8) ",
    "and consensus transcriptomic subtypes (epsilon² = 0.661, ",
    "p = 9.44 × 10^-6)."
  ),
  
  paste0(
    "In the independent GSE154918 cohort, all five genes reproduced the ",
    "prespecified direction, although the prespecified sepsis/septic-shock ",
    "versus uncomplicated-infection comparison showed only modest ",
    "discrimination (AUC = ",
    fmt_num(
      row_141_primary$auc_fixed_direction
    ),
    ", 95% CI ",
    fmt_num(
      row_141_primary$auc_ci_low
    ),
    "-",
    fmt_num(
      row_141_primary$auc_ci_high
    ),
    ")."
  ),
  
  paste0(
    "In GSE185263, the frozen score independently correlated with ",
    "continuous 24-hour SOFA (rho = ",
    fmt_num(
      primary_142$rho
    ),
    ", p = ",
    fmt_p(
      primary_142$p_value
    ),
    "; n = ",
    primary_142$n,
    "), with all five component genes reproducing the expected direction. ",
    "The association persisted after adjustment for age, sex and ",
    "collection location (beta = ",
    fmt_num(
      adjusted_142$estimate
    ),
    ", p = ",
    fmt_p(
      adjusted_142$p_value
    ),
    "). Geographic sensitivity analysis yielded a random-effects pooled ",
    "rho of ",
    fmt_num(
      location_meta$random_rho
    ),
    " (95% CI ",
    fmt_num(
      location_meta$random_CI_low
    ),
    "-",
    fmt_num(
      location_meta$random_CI_high
    ),
    "; I² = ",
    fmt_num(
      location_meta$I2_percent,
      digits = 1
    ),
    "%)."
  )
)


# ==============================================================================
# 22. FINAL TITLE OPTIONS
# ==============================================================================

title_options <- tibble::tibble(
  
  rank = 1:4,
  
  title = c(
    
    paste0(
      "Blood Transcriptomic Endotypes and a Five-Gene Host-Response ",
      "Signature for Molecular Stratification of Sepsis"
    ),
    
    paste0(
      "A Five-Gene Myeloid-Adaptive Host-Response Signature Captures ",
      "Blood Transcriptomic Endotypes and Organ-Dysfunction Severity in Sepsis"
    ),
    
    paste0(
      "Blood Transcriptomic Endotypes and a Reproducible Five-Gene ",
      "Host-Response Axis in Sepsis"
    ),
    
    paste0(
      "A Compact Five-Gene Host-Response Score Links Sepsis Blood ",
      "Endotypes to Organ-Dysfunction Severity"
    )
  ),
  
  comment = c(
    
    "Balanced and conservative; suitable default manuscript title.",
    
    "Strongest biological framing; highlights external severity result.",
    
    "Most concise multicohort framing.",
    
    "More mechanistic/interpretive; slightly stronger wording."
  )
)


# ==============================================================================
# 23. FINAL CENTRAL CLAIM
# ==============================================================================

central_claim <- c(
  
  "PRIMARY MANUSCRIPT CLAIM",
  
  "",
  
  paste0(
    "A compact five-gene myeloid-adaptive host-response score ",
    "recapitulates established blood transcriptomic endotypes and ",
    "reproducibly tracks host-response severity across independent ",
    "whole-blood RNA-seq cohorts."
  ),
  
  "",
  
  "SUPPORTED SECONDARY INTERPRETATION",
  
  "",
  
  paste0(
    "The score is associated with inflammatory activity and selected ",
    "severity-related clinical features, but these analyses do not ",
    "constitute independent prognostic validation."
  ),
  
  "",
  
  "NOT SUPPORTED",
  
  "",
  
  "- Clinically validated sepsis diagnostic test.",
  
  "- Validated prognostic mortality score.",
  
  "- Proven causal endotype mechanism.",
  
  "- Precalibrated clinical threshold.",
  
  "- Universal superiority over existing sepsis transcriptomic signatures."
)


# ==============================================================================
# 24. EXPORT CSV TABLES
# ==============================================================================

write.csv(
  locked_provenance,
  file.path(
    tables_dir,
    "144_locked_analysis_provenance.csv"
  ),
  row.names = FALSE
)


write.csv(
  internal_results,
  file.path(
    tables_dir,
    "144_internal_results_audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  external_results,
  file.path(
    tables_dir,
    "144_external_results_audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  claim_audit,
  file.path(
    tables_dir,
    "144_master_claim_audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  reviewer_risk,
  file.path(
    tables_dir,
    "144_reviewer_risk_audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  figure_plan,
  file.path(
    tables_dir,
    "144_final_figure_plan.csv"
  ),
  row.names = FALSE
)


write.csv(
  results_structure,
  file.path(
    tables_dir,
    "144_final_results_structure.csv"
  ),
  row.names = FALSE
)


write.csv(
  methods_map,
  file.path(
    tables_dir,
    "144_methods_analysis_provenance_map.csv"
  ),
  row.names = FALSE
)


write.csv(
  abstract_numbers,
  file.path(
    tables_dir,
    "144_abstract_key_numbers.csv"
  ),
  row.names = FALSE
)


write.csv(
  title_options,
  file.path(
    tables_dir,
    "144_title_options.csv"
  ),
  row.names = FALSE
)


write.csv(
  benchmark_summary,
  file.path(
    tables_dir,
    "144_benchmark_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 25. WRITE TEXT FILES
# ==============================================================================

results_en_file <- write_text(
  results_en,
  "144_final_results_section_EN.txt"
)


results_ru_file <- write_text(
  results_ru,
  "144_final_results_summary_RU.txt"
)


abstract_file <- write_text(
  abstract_results_en,
  "144_abstract_results_paragraph_EN.txt"
)


central_claim_file <- write_text(
  central_claim,
  "144_central_claim_and_boundaries.txt"
)


# ==============================================================================
# 26. MASTER EXCEL WORKBOOK
# ==============================================================================

run_info <- tibble::tibble(
  
  parameter = c(
    
    "script",
    
    "run_date",
    
    "project",
    
    "purpose",
    
    "frozen_panel",
    
    "score_formula",
    
    "new_feature_selection",
    
    "new_external_model_fitting",
    
    "new_cutoff_optimization",
    
    "new_endotype_creation",
    
    "new_cross_dataset_meta_analysis"
  ),
  
  value = c(
    
    script_name,
    
    as.character(
      run_date
    ),
    
    "Sepsis_DESeq2",
    
    "Final manuscript evidence audit and publication package",
    
    paste(
      five_genes,
      collapse = "; "
    ),
    
    five_gene_formula,
    
    "NO",
    
    "NO",
    
    "NO",
    
    "NO",
    
    "NO"
  )
)


workbook_file <- file.path(
  tables_dir,
  "144_final_manuscript_audit_package.xlsx"
)


wb <- openxlsx::createWorkbook()


sheet_list <- list(
  
  "00_run_info" =
    run_info,
  
  "01_provenance" =
    locked_provenance,
  
  "02_internal_results" =
    internal_results,
  
  "03_external_results" =
    external_results,
  
  "04_claim_audit" =
    claim_audit,
  
  "05_reviewer_risks" =
    reviewer_risk,
  
  "06_results_structure" =
    results_structure,
  
  "07_methods_map" =
    methods_map,
  
  "08_figure_plan" =
    figure_plan,
  
  "09_abstract_numbers" =
    abstract_numbers,
  
  "10_title_options" =
    title_options,
  
  "11_benchmark_summary" =
    benchmark_summary,
  
  "12_direction_matrix" =
    direction_143,
  
  "13_endotype_summary" =
    endotypes_143,
  
  "14_GSE154918_groups" =
    gse154918_groups,
  
  "15_GSE185263_gene_SOFA" =
    gene_sofa_142,
  
  "16_external_AUCs" =
    external_auc,
  
  "17_geographic_meta" =
    location_meta
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)


wrap_style <- openxlsx::createStyle(
  wrapText = TRUE,
  valign = "top"
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
  
  
  current_data <- sheet_list[[sheet_name]]
  
  
  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = current_data,
    headerStyle = header_style
  )
  
  
  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  
  if (
    nrow(
      current_data
    ) > 0
  ) {
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = wrap_style,
      rows = 2:(
        nrow(
          current_data
        ) + 1
      ),
      cols = seq_len(
        ncol(
          current_data
        )
      ),
      gridExpand = TRUE
    )
  }
  
  
  openxlsx::setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_len(
      max(
        1,
        ncol(
          current_data
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
# 27. INPUT MANIFEST
# ==============================================================================

manifest <- tibble::tibble(
  
  item = c(
    
    "Script143 key numbers",
    
    "Script143 evidence summary",
    
    "Script143 direction matrix",
    
    "Script143 endotype summary",
    
    "Script143 benchmark summary",
    
    "Script143 GSE154918 groups",
    
    "Script143 location meta-analysis",
    
    "Script143 external AUC summary",
    
    "Script141 comparisons",
    
    "Script142b primary SOFA",
    
    "Script142b secondary",
    
    "Script142b gene-SOFA",
    
    "Script142b adjusted model"
  ),
  
  path = c(
    
    file_143_key_numbers,
    
    file_143_evidence,
    
    file_143_direction,
    
    file_143_endotypes,
    
    file_143_benchmark,
    
    file_143_gse154918,
    
    file_143_location_meta,
    
    file_143_external_auc,
    
    file_141_comparisons,
    
    file_142_primary,
    
    file_142_secondary,
    
    file_142_gene_sofa,
    
    file_142_adjusted
  )
)


manifest_info <- file.info(
  manifest$path
)


manifest$file_size_bytes <- manifest_info$size


manifest$modified_time <- as.character(
  manifest_info$mtime
)


manifest$md5 <- unname(
  tools::md5sum(
    manifest$path
  )
)


write.csv(
  manifest,
  file.path(
    logs_dir,
    "144_input_manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 28. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "144_sessionInfo.txt"
  )
)


# ==============================================================================
# 29. FINAL CONSOLE REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 144 completed successfully.\n")
cat("====================================================================\n\n")


cat(
  "FINAL CENTRAL CLAIM:\n"
)


cat(
  paste0(
    "A compact five-gene myeloid-adaptive host-response score ",
    "recapitulates established blood transcriptomic endotypes and ",
    "reproducibly tracks host-response severity across independent ",
    "whole-blood RNA-seq cohorts.\n\n"
  )
)


cat(
  "DISCOVERY / INTERNAL:\n"
)


cat(
  "Robust core DEGs = 1796\n"
)


cat(
  "SRSq rho = 0.765; p = 8.742e-08\n"
)


cat(
  "CTS epsilon2 = 0.661; p = 9.437e-06\n"
)


cat(
  "Age/sex-adjusted BP-vs-BC score difference = 3.88; p = 3.73e-11\n\n"
)


cat(
  "GSE154918:\n"
)


cat(
  "Expected gene directions = 5/5\n"
)


cat(
  "Primary AUC = ",
  fmt_num(
    row_141_primary$auc_fixed_direction
  ),
  " [",
  fmt_num(
    row_141_primary$auc_ci_low
  ),
  ", ",
  fmt_num(
    row_141_primary$auc_ci_high
  ),
  "]; p = ",
  fmt_p(
    row_141_primary$p_value
  ),
  "\n",
  sep = ""
)


cat(
  "Shock vs uncomplicated AUC = ",
  fmt_num(
    row_141_shock$auc_fixed_direction
  ),
  "; p = ",
  fmt_p(
    row_141_shock$p_value
  ),
  "\n\n",
  sep = ""
)


cat(
  "GSE185263:\n"
)


cat(
  "Primary score-SOFA rho = ",
  fmt_num(
    primary_142$rho
  ),
  "; p = ",
  fmt_p(
    primary_142$p_value
  ),
  "; n = ",
  primary_142$n,
  "\n",
  sep = ""
)


cat(
  "Expected gene-SOFA directions = ",
  sum(
    gene_sofa_142$direction_concordant,
    na.rm = TRUE
  ),
  "/5\n",
  sep = ""
)


cat(
  "Adjusted SOFA beta = ",
  fmt_num(
    adjusted_142$estimate
  ),
  "; p = ",
  fmt_p(
    adjusted_142$p_value
  ),
  "\n",
  sep = ""
)


cat(
  "Geographic random-effects rho = ",
  fmt_num(
    location_meta$random_rho
  ),
  " [",
  fmt_num(
    location_meta$random_CI_low
  ),
  ", ",
  fmt_num(
    location_meta$random_CI_high
  ),
  "]\n",
  sep = ""
)


cat(
  "I2 = ",
  fmt_num(
    location_meta$I2_percent,
    digits = 1
  ),
  "%; heterogeneity p = ",
  fmt_p(
    location_meta$Q_p
  ),
  "\n\n",
  sep = ""
)


cat(
  "MASTER WORKBOOK:\n"
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
  "FINAL RESULTS SECTION:\n"
)


cat(
  normalizePath(
    results_en_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "ABSTRACT RESULTS PARAGRAPH:\n"
)


cat(
  normalizePath(
    abstract_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "CLAIM BOUNDARIES:\n"
)


cat(
  normalizePath(
    central_claim_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "MANUSCRIPT AUDIT SUMMARY:\n"
)


cat(
  "- No new feature selection.\n"
)


cat(
  "- No external model refitting.\n"
)


cat(
  "- No cutoff optimization.\n"
)


cat(
  "- No new endotype creation.\n"
)


cat(
  "- Internal AUC=1 retained as internal characterization only.\n"
)


cat(
  "- GSE154918 nonsignificant prespecified primary result retained transparently.\n"
)


cat(
  "- GSE185263 score-SOFA association retained as primary external severity evidence.\n"
)


cat(
  "- Geographic strata treated as sensitivity analysis within one study.\n"
)


cat(
  "- Clinical diagnostic/prognostic validation claims explicitly prohibited.\n\n"
)


cat(
  "Done.\n"
)