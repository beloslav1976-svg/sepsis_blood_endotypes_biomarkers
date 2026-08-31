################################################################################
# Script 165 FINAL v2
#
# Supplementary Table S10
# External evaluation of the frozen five-gene host-response score in GSE185263
#
# PURPOSE
# -------
# Package the previously frozen GSE185263 external-evaluation results into a
# publication-ready Supplementary Table S10.
#
# FROZEN PRIMARY PANEL
# --------------------
# UP components:
#   CD177, HK3, IRAK3
#
# DOWN components:
#   CARD11, IKZF2
#
# Five-gene score:
#   mean z(CD177, HK3, IRAK3) - mean z(CARD11, IKZF2)
#
# PRIMARY EXTERNAL ENDPOINT
# -------------------------
# Spearman association between the frozen five-gene score and continuous
# 24-h SOFA among patients with sepsis.
#
# IMPORTANT
# ---------
# This script performs NO new inferential statistical analysis.
#
# It does NOT:
#   - reselect genes
#   - refit coefficients
#   - optimize a cutoff
#   - reverse score direction
#   - recompute Spearman correlations
#   - recompute Wilcoxon tests
#   - recompute ROC curves or AUCs
#   - recompute confidence intervals
#   - recompute BH-adjusted P values
#   - recompute adjusted regression models
#   - reconstruct endotypes
#
# All inferential quantities are copied from frozen Scripts 142b, 143, and 150.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 165 FINAL v2\n")
cat("Supplementary Table S10\n")
cat("GSE185263 external severity evaluation\n")
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
# 2. REQUIRED PACKAGES
# =============================================================================

required_packages <- c(
  "readxl",
  "openxlsx",
  "dplyr",
  "stringr"
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
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(stringr)
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "165_TableS10_GSE185263_external_evaluation"
)


tables_dir <- file.path(
  output_dir,
  "tables"
)


audit_dir <- file.path(
  output_dir,
  "audit"
)


text_dir <- file.path(
  output_dir,
  "text"
)


dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  text_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


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
# 4. FROZEN SOURCE DIRECTORIES
# =============================================================================

source142_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "142b_external_validation_GSE185263",
  "tables"
)


source143_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "143_multicohort_integration",
  "tables"
)


source150_dir <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "150_Figure5_GSE185263_severity_validation",
  "tables"
)


if (!dir.exists(source142_dir)) {
  stop("Frozen Script 142b table directory not found.")
}


if (!dir.exists(source150_dir)) {
  stop("Frozen Main Figure 5 source directory not found.")
}


# =============================================================================
# 5. CANONICAL FROZEN SCRIPT 142b FILES
# =============================================================================

f_primary <- file.path(
  source142_dir,
  "142b_PRIMARY_score_vs_SOFA.csv"
)


f_secondary <- file.path(
  source142_dir,
  "142b_secondary_score_associations.csv"
)


f_gene_sofa <- file.path(
  source142_dir,
  "142b_component_gene_SOFA_associations.csv"
)


f_adjusted_model <- file.path(
  source142_dir,
  "142b_age_sex_location_adjusted_model.csv"
)


f_adjusted_summary <- file.path(
  source142_dir,
  "142b_adjusted_model_summary.csv"
)


f_location_sofa <- file.path(
  source142_dir,
  "142b_location_specific_SOFA_correlations.csv"
)


f_location_pooled <- file.path(
  source142_dir,
  "142b_cross_location_pooled_correlation.csv"
)


f_scaling <- file.path(
  source142_dir,
  "142b_scaling_sensitivity.csv"
)


f_sofa_summary <- file.path(
  source142_dir,
  "142b_SOFA_summary.csv"
)


f_sofa_groups <- file.path(
  source142_dir,
  "142b_SOFA_group_summary.csv"
)


f_mortality <- file.path(
  source142_dir,
  "142b_mortality_summary.csv"
)


f_site_summary <- file.path(
  source142_dir,
  "142b_collection_site_summary.csv"
)


f_location_summary <- file.path(
  source142_dir,
  "142b_collection_location_summary.csv"
)


f_gene_coverage <- file.path(
  source142_dir,
  "142b_frozen_gene_coverage.csv"
)


f_sample_scores <- file.path(
  source142_dir,
  "142b_GSE185263_sample_scores.csv"
)


f_workbook142 <- file.path(
  source142_dir,
  "142b_GSE185263_external_validation.xlsx"
)


required142_files <- c(
  f_primary,
  f_secondary,
  f_gene_sofa,
  f_adjusted_model,
  f_adjusted_summary,
  f_location_sofa,
  f_location_pooled,
  f_scaling,
  f_sofa_summary,
  f_sofa_groups,
  f_mortality,
  f_site_summary,
  f_location_summary,
  f_gene_coverage,
  f_sample_scores,
  f_workbook142
)


missing142 <- required142_files[
  !file.exists(
    required142_files
  )
]


if (length(missing142) > 0) {
  
  stop(
    paste0(
      "Missing frozen Script 142b file(s):\n",
      paste(
        missing142,
        collapse = "\n"
      )
    )
  )
}


cat("\n## CANONICAL FROZEN SCRIPT 142b SOURCES\n\n")

print(
  normalizePath(
    required142_files,
    winslash = "\\",
    mustWork = TRUE
  )
)


# =============================================================================
# 6. OPTIONAL / FROZEN SCRIPT 143 LOCATION-CI FILES
# =============================================================================

f_location_ci <- file.path(
  source143_dir,
  "143_GSE185263_location_correlations_with_CI.csv"
)


f_location_meta <- file.path(
  source143_dir,
  "143_GSE185263_location_random_effects_meta.csv"
)


location_ci_available <-
  file.exists(f_location_ci) &&
  file.exists(f_location_meta)


cat(
  "\nFrozen Script 143 location-CI sensitivity available: ",
  location_ci_available,
  "\n",
  sep = ""
)


# =============================================================================
# 7. FROZEN MAIN FIGURE 5 AUDIT SOURCES
# =============================================================================

f_fig5_audit <- file.path(
  source150_dir,
  "150_Figure5_numerical_audit.csv"
)


f_fig5_workbook <- file.path(
  source150_dir,
  "150_Figure5_source_data.xlsx"
)


if (!file.exists(f_fig5_audit)) {
  stop("Frozen Main Figure 5 numerical audit not found.")
}


if (!file.exists(f_fig5_workbook)) {
  stop("Frozen Main Figure 5 source workbook not found.")
}


cat(
  "\nMain Figure 5 frozen numerical audit available: TRUE\n"
)


# =============================================================================
# 8. READ CANONICAL CSV TABLES
# =============================================================================

read_csv_frozen <- function(path) {
  
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


primary <- read_csv_frozen(
  f_primary
)


secondary <- read_csv_frozen(
  f_secondary
)


gene_sofa <- read_csv_frozen(
  f_gene_sofa
)


adjusted_model <- read_csv_frozen(
  f_adjusted_model
)


adjusted_summary <- read_csv_frozen(
  f_adjusted_summary
)


location_sofa <- read_csv_frozen(
  f_location_sofa
)


location_pooled <- read_csv_frozen(
  f_location_pooled
)


scaling <- read_csv_frozen(
  f_scaling
)


sofa_summary <- read_csv_frozen(
  f_sofa_summary
)


sofa_groups <- read_csv_frozen(
  f_sofa_groups
)


mortality <- read_csv_frozen(
  f_mortality
)


site_summary <- read_csv_frozen(
  f_site_summary
)


location_summary <- read_csv_frozen(
  f_location_summary
)


gene_coverage <- read_csv_frozen(
  f_gene_coverage
)


sample_scores <- read_csv_frozen(
  f_sample_scores
)


run_info <- readxl::read_excel(
  f_workbook142,
  sheet = "00_run_info"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


# =============================================================================
# 9. EXACT KNOWN SCHEMA AUDIT
# =============================================================================

expected_secondary_cols <- c(
  "analysis",
  "score",
  "case",
  "control",
  "n_case",
  "n_control",
  "median_case",
  "median_control",
  "median_difference",
  "W",
  "p_value",
  "AUC",
  "CI_low",
  "CI_high",
  "BH_secondary"
)


expected_gene_cols <- c(
  "gene",
  "expected_direction",
  "observed_direction",
  "direction_concordant",
  "n",
  "rho",
  "p_value",
  "BH_five_genes"
)


expected_adjusted_cols <- c(
  "term",
  "estimate",
  "SE",
  "t_value",
  "p_value"
)


expected_gene_coverage_cols <- c(
  "gene",
  "expected_sofa_direction",
  "present"
)


expected_sample_cols <- c(
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
  "sofa_group",
  "CD177_logCPM",
  "CD177_z_sepsis",
  "HK3_logCPM",
  "HK3_z_sepsis",
  "IRAK3_logCPM",
  "IRAK3_z_sepsis",
  "CARD11_logCPM",
  "CARD11_z_sepsis",
  "IKZF2_logCPM",
  "IKZF2_z_sepsis"
)


schema_audit <- data.frame(
  
  object = c(
    "secondary associations",
    "component-gene SOFA",
    "adjusted model",
    "gene coverage",
    "sample scores"
  ),
  
  schema_match = c(
    identical(
      names(secondary),
      expected_secondary_cols
    ),
    identical(
      names(gene_sofa),
      expected_gene_cols
    ),
    identical(
      names(adjusted_model),
      expected_adjusted_cols
    ),
    identical(
      names(gene_coverage),
      expected_gene_coverage_cols
    ),
    identical(
      names(sample_scores),
      expected_sample_cols
    )
  ),
  
  stringsAsFactors = FALSE
)


cat("\n## EXACT FROZEN SCHEMA AUDIT\n\n")

print(
  schema_audit,
  row.names = FALSE
)


if (!all(schema_audit$schema_match)) {
  
  stop(
    "At least one canonical frozen Script 142b schema does not match expectation."
  )
}


# =============================================================================
# 10. READ SCRIPT 142b WORKBOOK COPIES
# =============================================================================

workbook_sheet_map <- list(
  
  "Primary SOFA" =
    "01_PRIMARY_SOFA",
  
  "Secondary associations" =
    "02_secondary",
  
  "Component-gene SOFA" =
    "03_gene_SOFA",
  
  "Adjusted model" =
    "04_adjusted_model",
  
  "Adjusted summary" =
    "05_adjusted_summary",
  
  "Location SOFA" =
    "06_location_SOFA",
  
  "Location pooled" =
    "07_location_pooled",
  
  "Scaling sensitivity" =
    "08_scaling",
  
  "SOFA summary" =
    "09_SOFA_summary",
  
  "SOFA groups" =
    "10_SOFA_groups",
  
  "Mortality summary" =
    "11_mortality",
  
  "Site summary" =
    "12_site_summary",
  
  "Location summary" =
    "13_location_summary",
  
  "Gene coverage" =
    "14_gene_coverage",
  
  "Sample scores" =
    "15_sample_scores"
)


csv_table_map <- list(
  
  "Primary SOFA" =
    primary,
  
  "Secondary associations" =
    secondary,
  
  "Component-gene SOFA" =
    gene_sofa,
  
  "Adjusted model" =
    adjusted_model,
  
  "Adjusted summary" =
    adjusted_summary,
  
  "Location SOFA" =
    location_sofa,
  
  "Location pooled" =
    location_pooled,
  
  "Scaling sensitivity" =
    scaling,
  
  "SOFA summary" =
    sofa_summary,
  
  "SOFA groups" =
    sofa_groups,
  
  "Mortality summary" =
    mortality,
  
  "Site summary" =
    site_summary,
  
  "Location summary" =
    location_summary,
  
  "Gene coverage" =
    gene_coverage,
  
  "Sample scores" =
    sample_scores
)


# =============================================================================
# 11. ROBUST CSV-vs-WORKBOOK TABLE COMPARISON
# =============================================================================

compare_vectors <- function(x, y, tolerance = 1e-8) {
  
  if (length(x) != length(y)) {
    return(FALSE)
  }
  
  
  x_num <- suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
  
  
  y_num <- suppressWarnings(
    as.numeric(
      as.character(y)
    )
  )
  
  
  x_nonmissing <- !is.na(x)
  y_nonmissing <- !is.na(y)
  
  
  x_numeric_fraction <- if (
    sum(x_nonmissing) > 0
  ) {
    mean(
      is.finite(
        x_num[x_nonmissing]
      )
    )
  } else {
    1
  }
  
  
  y_numeric_fraction <- if (
    sum(y_nonmissing) > 0
  ) {
    mean(
      is.finite(
        y_num[y_nonmissing]
      )
    )
  } else {
    1
  }
  
  
  if (
    x_numeric_fraction == 1 &&
    y_numeric_fraction == 1
  ) {
    
    same_na <- identical(
      is.na(x_num),
      is.na(y_num)
    )
    
    
    if (!same_na) {
      return(FALSE)
    }
    
    
    good <- is.finite(x_num) &
      is.finite(y_num)
    
    
    if (!any(good)) {
      return(TRUE)
    }
    
    
    return(
      all(
        abs(
          x_num[good] -
            y_num[good]
        ) <=
          tolerance
      )
    )
  }
  
  
  x_chr <- trimws(
    as.character(x)
  )
  
  
  y_chr <- trimws(
    as.character(y)
  )
  
  
  x_chr[is.na(x_chr)] <- "<NA>"
  y_chr[is.na(y_chr)] <- "<NA>"
  
  
  identical(
    x_chr,
    y_chr
  )
}


compare_tables <- function(
    csv_df,
    workbook_df,
    tolerance = 1e-8
) {
  
  dimensions_match <-
    identical(
      dim(csv_df),
      dim(workbook_df)
    )
  
  
  columns_match <-
    identical(
      names(csv_df),
      names(workbook_df)
    )
  
  
  if (
    !dimensions_match ||
    !columns_match
  ) {
    
    return(
      list(
        dimensions_match = dimensions_match,
        columns_match = columns_match,
        content_match = FALSE,
        overall_match = FALSE
      )
    )
  }
  
  
  column_matches <- vapply(
    
    names(csv_df),
    
    function(column_name) {
      
      compare_vectors(
        csv_df[[column_name]],
        workbook_df[[column_name]],
        tolerance = tolerance
      )
    },
    
    logical(1)
  )
  
  
  content_match <- all(
    column_matches
  )
  
  
  list(
    
    dimensions_match =
      dimensions_match,
    
    columns_match =
      columns_match,
    
    content_match =
      content_match,
    
    overall_match =
      dimensions_match &&
      columns_match &&
      content_match
  )
}


equivalence_rows <- list()


for (
  source_name in names(
    workbook_sheet_map
  )
) {
  
  workbook_copy <- readxl::read_excel(
    
    f_workbook142,
    
    sheet =
      workbook_sheet_map[[source_name]]
    
  ) %>%
    
    as.data.frame(
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  
  
  comparison_result <- compare_tables(
    
    csv_df =
      csv_table_map[[source_name]],
    
    workbook_df =
      workbook_copy
  )
  
  
  equivalence_rows[[source_name]] <- data.frame(
    
    source =
      source_name,
    
    dimensions_match =
      comparison_result$dimensions_match,
    
    columns_match =
      comparison_result$columns_match,
    
    content_match =
      comparison_result$content_match,
    
    overall_match =
      comparison_result$overall_match,
    
    stringsAsFactors = FALSE
  )
}


equivalence_audit <- dplyr::bind_rows(
  equivalence_rows
)


cat("\n## CSV-vs-WORKBOOK EQUIVALENCE AUDIT\n\n")

print(
  equivalence_audit,
  row.names = FALSE
)


if (!all(equivalence_audit$overall_match)) {
  
  stop(
    "At least one canonical Script 142b CSV table differs from its frozen workbook copy."
  )
}


cat(
  "\nAll 15 canonical Script 142b CSV tables are equivalent to their workbook copies.\n"
)


# =============================================================================
# 12. RUN-INFO EXTRACTION
# =============================================================================

if (
  !all(
    c(
      "parameter",
      "value"
    ) %in%
    names(run_info)
  )
) {
  
  stop(
    "Unexpected Script 142b run-info schema."
  )
}


run_value <- function(parameter_name) {
  
  hit <- run_info[
    run_info$parameter ==
      parameter_name,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(hit) != 1) {
    
    stop(
      paste0(
        "Run-info parameter not uniquely found: ",
        parameter_name
      )
    )
  }
  
  
  as.character(
    hit$value[1]
  )
}


total_samples <- as.numeric(
  run_value(
    "total_samples"
  )
)


sepsis_samples <- as.numeric(
  run_value(
    "sepsis_samples"
  )
)


healthy_samples <- as.numeric(
  run_value(
    "healthy_samples"
  )
)


sofa_available <- as.numeric(
  run_value(
    "SOFA_available_sepsis"
  )
)


mortality_available <- as.numeric(
  run_value(
    "mortality_available_sepsis"
  )
)


primary_panel <- run_value(
  "primary_panel"
)


primary_score_definition <- run_value(
  "primary_score"
)


primary_endpoint <- run_value(
  "primary_endpoint"
)


guardrail_parameters <- c(
  "feature_selection_external",
  "coefficient_refitting_external",
  "cutoff_optimization_external",
  "score_direction_flipping_external",
  "endotype_reconstruction"
)


guardrail_values <- vapply(
  
  guardrail_parameters,
  
  run_value,
  
  character(1)
)


run_guardrails <- data.frame(
  
  item = c(
    "GEO_accession",
    "total_samples",
    "sepsis_samples",
    "healthy_samples",
    "SOFA_available_sepsis",
    "mortality_available_sepsis",
    "primary_panel",
    "primary_score",
    "primary_endpoint",
    guardrail_parameters
  ),
  
  value = c(
    run_value(
      "GEO_accession"
    ),
    total_samples,
    sepsis_samples,
    healthy_samples,
    sofa_available,
    mortality_available,
    primary_panel,
    primary_score_definition,
    primary_endpoint,
    guardrail_values
  ),
  
  stringsAsFactors = FALSE
)


cat("\n## SCRIPT 142b RUN-INFO GUARDRAILS\n\n")

print(
  run_guardrails,
  row.names = FALSE
)


if (
  !all(
    toupper(
      guardrail_values
    ) ==
    "NO"
  )
) {
  
  stop(
    "At least one frozen external-analysis guardrail is not NO."
  )
}


# =============================================================================
# 13. SAMPLE-LEVEL COHORT AUDIT
# =============================================================================

sample_total_observed <- nrow(
  sample_scores
)


sample_sepsis_observed <- sum(
  tolower(
    sample_scores$disease_state
  ) ==
    "sepsis",
  na.rm = TRUE
)


sample_healthy_observed <- sum(
  tolower(
    sample_scores$disease_state
  ) ==
    "healthy",
  na.rm = TRUE
)


sepsis_rows <- tolower(
  sample_scores$disease_state
) ==
  "sepsis"


sample_sofa_observed <- sum(
  sepsis_rows &
    !is.na(
      sample_scores$sofa
    )
)


sample_mortality_observed <- sum(
  sepsis_rows &
    !is.na(
      sample_scores$mortality
    ) &
    trimws(
      sample_scores$mortality
    ) != ""
)


cohort_audit <- data.frame(
  
  metric = c(
    "Total samples",
    "Sepsis samples",
    "Healthy samples",
    "SOFA-complete sepsis samples",
    "Mortality-complete sepsis samples"
  ),
  
  expected = c(
    total_samples,
    sepsis_samples,
    healthy_samples,
    sofa_available,
    mortality_available
  ),
  
  observed = c(
    sample_total_observed,
    sample_sepsis_observed,
    sample_healthy_observed,
    sample_sofa_observed,
    sample_mortality_observed
  ),
  
  stringsAsFactors = FALSE
)


cohort_audit$pass <-
  cohort_audit$expected ==
  cohort_audit$observed


cat("\n## COHORT AUDIT\n\n")

print(
  cohort_audit,
  row.names = FALSE
)


if (!all(cohort_audit$pass)) {
  
  stop(
    "Sample-level cohort counts do not match frozen Script 142b run-info anchors."
  )
}


# =============================================================================
# 14. PRIMARY SOFA COLUMN DETECTION
# =============================================================================

find_one_column <- function(
    df,
    candidates
) {
  
  hits <- candidates[
    candidates %in%
      names(df)
  ]
  
  
  if (length(hits) != 1) {
    
    stop(
      paste0(
        "Could not uniquely identify column. Candidates: ",
        paste(
          candidates,
          collapse = ", "
        ),
        ". Available columns: ",
        paste(
          names(df),
          collapse = ", "
        )
      )
    )
  }
  
  
  hits[1]
}


primary_n_col <- find_one_column(
  primary,
  c(
    "n",
    "N"
  )
)


primary_rho_col <- find_one_column(
  primary,
  c(
    "rho",
    "spearman_rho",
    "Spearman_rho"
  )
)


primary_p_col <- find_one_column(
  primary,
  c(
    "p_value",
    "P",
    "p"
  )
)


primary_n <- as.numeric(
  primary[[primary_n_col]][1]
)


primary_rho <- as.numeric(
  primary[[primary_rho_col]][1]
)


primary_p <- as.numeric(
  primary[[primary_p_col]][1]
)


# =============================================================================
# 15. SECONDARY ANALYSIS ROWS
# =============================================================================

get_secondary_row <- function(
    analysis_name
) {
  
  out <- secondary[
    secondary$analysis ==
      analysis_name,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(out) != 1) {
    
    stop(
      paste0(
        "Secondary analysis row not uniquely found: ",
        analysis_name
      )
    )
  }
  
  
  out
}


sofa_binary <- get_secondary_row(
  "SOFA_ge2_vs_SOFA_0_1"
)


mortality_test <- get_secondary_row(
  "Died_vs_Survived"
)


icu_test <- get_secondary_row(
  "ICU_vs_Emergency_Room"
)


sepsis_healthy <- get_secondary_row(
  "Sepsis_vs_healthy_contextual"
)


# =============================================================================
# 16. COMPONENT-GENE SOFA AUDIT
# =============================================================================

expected_gene_order <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "CARD11",
  "IKZF2"
)


if (
  !setequal(
    gene_sofa$gene,
    expected_gene_order
  )
) {
  
  stop(
    "Frozen component-gene table does not contain exactly the five primary genes."
  )
}


n_gene_direction_concordant <- sum(
  as.logical(
    gene_sofa$direction_concordant
  ),
  na.rm = TRUE
)


n_gene_BH_significant <- sum(
  gene_sofa$BH_five_genes <
    0.05,
  na.rm = TRUE
)


cat("\n## COMPONENT-GENE SOFA REPLICATION\n\n")

cat(
  "Directionally concordant genes = ",
  n_gene_direction_concordant,
  "/5\n",
  sep = ""
)


cat(
  "BH-adjusted P < 0.05 = ",
  n_gene_BH_significant,
  "/5\n",
  sep = ""
)


if (
  n_gene_direction_concordant !=
  5
) {
  
  stop(
    "Expected 5/5 directional concordance was not reproduced."
  )
}


if (
  n_gene_BH_significant !=
  5
) {
  
  stop(
    "Expected 5/5 BH-significant component genes were not reproduced."
  )
}


# =============================================================================
# 17. ADJUSTED MODEL SOFA ROW
# =============================================================================

adjusted_sofa_row <- adjusted_model[
  tolower(
    adjusted_model$term
  ) ==
    "sofa",
  ,
  drop = FALSE
]


if (nrow(adjusted_sofa_row) != 1) {
  
  stop(
    "SOFA term not uniquely found in frozen adjusted model."
  )
}


adjusted_beta <- as.numeric(
  adjusted_sofa_row$estimate[1]
)


adjusted_se <- as.numeric(
  adjusted_sofa_row$SE[1]
)


adjusted_p <- as.numeric(
  adjusted_sofa_row$p_value[1]
)


# =============================================================================
# 18. EXTRACT FROZEN ADJUSTED MODEL R-SQUARED
# =============================================================================

extract_numeric_value_from_matching_row <- function(
    df,
    pattern
) {
  
  row_text <- apply(
    df,
    1,
    function(x) {
      paste(
        x,
        collapse = " "
      )
    }
  )
  
  
  hits <- which(
    grepl(
      pattern,
      row_text,
      ignore.case = TRUE
    )
  )
  
  
  if (length(hits) == 0) {
    return(NA_real_)
  }
  
  
  one_row <- df[
    hits[1],
    ,
    drop = FALSE
  ]
  
  
  preferred_columns <- c(
    "value",
    "estimate",
    "R2",
    "r2",
    "R_squared",
    "r_squared"
  )
  
  
  preferred_columns <- preferred_columns[
    preferred_columns %in%
      names(one_row)
  ]
  
  
  if (length(preferred_columns) > 0) {
    
    for (
      column_name in preferred_columns
    ) {
      
      one_value <- suppressWarnings(
        as.numeric(
          as.character(
            one_row[[column_name]][1]
          )
        )
      )
      
      
      if (is.finite(one_value)) {
        return(one_value)
      }
    }
  }
  
  
  for (
    column_name in names(one_row)
  ) {
    
    one_value <- suppressWarnings(
      as.numeric(
        as.character(
          one_row[[column_name]][1]
        )
      )
    )
    
    
    if (is.finite(one_value)) {
      return(one_value)
    }
  }
  
  
  NA_real_
}


adjusted_r2 <- extract_numeric_value_from_matching_row(
  
  adjusted_summary,
  
  "R.?2|R.?squared|R_squared|r_squared"
)


# =============================================================================
# 19. FROZEN DISPLAY CI FOR ADJUSTED SOFA FROM MAIN FIGURE 5
# =============================================================================

adjusted_sofa_display <- readxl::read_excel(
  f_fig5_workbook,
  sheet = "Adjusted_SOFA_display"
) %>%
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


if (
  !identical(
    names(adjusted_sofa_display),
    c(
      "term",
      "beta",
      "CI_low",
      "CI_high"
    )
  )
) {
  
  stop(
    "Unexpected frozen Figure 5 Adjusted_SOFA_display schema."
  )
}


if (
  nrow(adjusted_sofa_display) != 1 ||
  toupper(
    adjusted_sofa_display$term[1]
  ) !=
  "SOFA"
) {
  
  stop(
    "Unexpected frozen Figure 5 adjusted-SOFA display row."
  )
}


adjusted_ci_low <- as.numeric(
  adjusted_sofa_display$CI_low[1]
)


adjusted_ci_high <- as.numeric(
  adjusted_sofa_display$CI_high[1]
)


# =============================================================================
# 20. SCRIPT 143 LOCATION-SENSITIVITY TABLES
# =============================================================================

if (location_ci_available) {
  
  location_ci <- read_csv_frozen(
    f_location_ci
  )
  
  
  location_meta_all <- read_csv_frozen(
    f_location_meta
  )
  
  
  required_location_ci_cols <- c(
    "collection_location",
    "direction_concordant",
    "n",
    "rho",
    "p_value",
    "BH_location",
    "CI_low",
    "CI_high"
  )
  
  
  if (
    !identical(
      names(location_ci),
      required_location_ci_cols
    )
  ) {
    
    stop(
      "Unexpected frozen Script 143 location-correlation schema."
    )
  }
  
  
  required_meta_cols <- c(
    "k",
    "total_n",
    "fixed_rho",
    "fixed_CI_low",
    "fixed_CI_high",
    "fixed_p"
  )
  
  
  if (
    !all(
      required_meta_cols %in%
      names(location_meta_all)
    )
  ) {
    
    stop(
      "Required fixed-effect columns not found in frozen Script 143 location meta table."
    )
  }
  
  
  location_pooled_fixed <- location_meta_all %>%
    
    dplyr::select(
      dplyr::all_of(
        required_meta_cols
      )
    )
  
  
  if (nrow(location_pooled_fixed) != 1) {
    
    stop(
      "Unexpected number of frozen fixed-effect location-pooling rows."
    )
  }
  
  
  location_join_audit <- merge(
    
    location_sofa,
    
    location_ci[
      ,
      c(
        "collection_location",
        "n",
        "rho",
        "p_value",
        "BH_location"
      ),
      drop = FALSE
    ],
    
    by = "collection_location",
    
    suffixes = c(
      "_142b",
      "_143"
    ),
    
    all = TRUE
  )
  
  
  location_join_audit$n_match <-
    location_join_audit$n_142b ==
    location_join_audit$n_143
  
  
  location_join_audit$rho_match <-
    abs(
      location_join_audit$rho_142b -
        location_join_audit$rho_143
    ) <
    1e-10
  
  
  location_join_audit$p_match <-
    abs(
      location_join_audit$p_value_142b -
        location_join_audit$p_value_143
    ) <
    1e-10
  
  
  location_join_audit$BH_match <-
    abs(
      location_join_audit$BH_location_142b -
        location_join_audit$BH_location_143
    ) <
    1e-10
  
  
  location_join_audit$overall_match <-
    location_join_audit$n_match &
    location_join_audit$rho_match &
    location_join_audit$p_match &
    location_join_audit$BH_match
  
  
  if (
    !all(
      location_join_audit$overall_match
    )
  ) {
    
    stop(
      "Script 142b and Script 143 location-specific correlations do not match."
    )
  }
  
  
} else {
  
  location_ci <- data.frame()
  
  location_pooled_fixed <- data.frame()
  
  location_join_audit <- data.frame()
}


# =============================================================================
# 21. FROZEN NUMERICAL ANCHOR AUDIT
# =============================================================================

assert_close <- function(
    observed,
    expected,
    tolerance,
    label
) {
  
  absolute_difference <- abs(
    observed -
      expected
  )
  
  
  pass <- is.finite(observed) &&
    absolute_difference <=
    tolerance
  
  
  data.frame(
    
    metric =
      label,
    
    observed =
      observed,
    
    expected =
      expected,
    
    tolerance =
      tolerance,
    
    absolute_difference =
      absolute_difference,
    
    pass =
      pass,
    
    stringsAsFactors = FALSE
  )
}


anchor_audit <- dplyr::bind_rows(
  
  assert_close(
    total_samples,
    392,
    0,
    "Total samples"
  ),
  
  assert_close(
    sepsis_samples,
    348,
    0,
    "Sepsis samples"
  ),
  
  assert_close(
    healthy_samples,
    44,
    0,
    "Healthy samples"
  ),
  
  assert_close(
    sofa_available,
    345,
    0,
    "SOFA-complete sepsis samples"
  ),
  
  assert_close(
    primary_n,
    345,
    0,
    "Primary SOFA n"
  ),
  
  assert_close(
    primary_rho,
    0.3114963,
    1e-6,
    "Primary score-SOFA rho"
  ),
  
  assert_close(
    primary_p,
    3.368874e-09,
    1e-12,
    "Primary score-SOFA P"
  ),
  
  assert_close(
    as.numeric(
      sofa_binary$n_case[1]
    ),
    207,
    0,
    "SOFA >=2 n"
  ),
  
  assert_close(
    as.numeric(
      sofa_binary$n_control[1]
    ),
    138,
    0,
    "SOFA 0-1 n"
  ),
  
  assert_close(
    as.numeric(
      sofa_binary$AUC[1]
    ),
    0.6698873,
    1e-6,
    "SOFA binary AUC"
  ),
  
  assert_close(
    as.numeric(
      sofa_binary$p_value[1]
    ),
    8.952612e-08,
    1e-11,
    "SOFA binary P"
  ),
  
  assert_close(
    as.numeric(
      mortality_test$AUC[1]
    ),
    0.6272644,
    1e-6,
    "Mortality AUC"
  ),
  
  assert_close(
    as.numeric(
      mortality_test$p_value[1]
    ),
    0.003447767,
    1e-8,
    "Mortality P"
  ),
  
  assert_close(
    as.numeric(
      icu_test$AUC[1]
    ),
    0.6438658,
    1e-6,
    "ICU vs Emergency Room AUC"
  ),
  
  assert_close(
    as.numeric(
      icu_test$p_value[1]
    ),
    8.173562e-05,
    1e-9,
    "ICU vs Emergency Room P"
  ),
  
  assert_close(
    as.numeric(
      sepsis_healthy$AUC[1]
    ),
    0.9494514,
    1e-6,
    "Contextual sepsis vs healthy AUC"
  ),
  
  assert_close(
    adjusted_beta,
    0.1166093,
    1e-6,
    "Adjusted SOFA beta"
  ),
  
  assert_close(
    adjusted_se,
    0.0260606,
    1e-6,
    "Adjusted SOFA SE"
  ),
  
  assert_close(
    adjusted_p,
    1.048676e-05,
    1e-9,
    "Adjusted SOFA P"
  )
)


if (is.finite(adjusted_r2)) {
  
  anchor_audit <- dplyr::bind_rows(
    
    anchor_audit,
    
    assert_close(
      adjusted_r2,
      0.2610886,
      1e-4,
      "Adjusted model R-squared"
    )
  )
}


if (location_ci_available) {
  
  anchor_audit <- dplyr::bind_rows(
    
    anchor_audit,
    
    assert_close(
      as.numeric(
        location_pooled_fixed$fixed_rho[1]
      ),
      0.2407839,
      1e-6,
      "Descriptive fixed-effect location pooled rho"
    )
  )
}


cat("\n## FROZEN NUMERICAL ANCHOR AUDIT\n\n")

print(
  anchor_audit,
  row.names = FALSE
)


if (!all(anchor_audit$pass)) {
  
  stop(
    "At least one frozen GSE185263 numerical anchor failed."
  )
}


# =============================================================================
# 22. MAIN FIGURE 5 CONSISTENCY AUDIT
# =============================================================================

fig5_audit <- read_csv_frozen(
  f_fig5_audit
)


if (
  !all(
    c(
      "metric",
      "observed"
    ) %in%
    names(fig5_audit)
  )
) {
  
  stop(
    "Unexpected Main Figure 5 numerical-audit schema."
  )
}


fig5_value <- function(
    metric_name
) {
  
  hit <- fig5_audit[
    fig5_audit$metric ==
      metric_name,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(hit) != 1) {
    
    stop(
      paste0(
        "Main Figure 5 metric not uniquely found: ",
        metric_name
      )
    )
  }
  
  
  as.numeric(
    hit$observed[1]
  )
}


fig5_consistency <- data.frame(
  
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
    "Adjusted SOFA SE"
  ),
  
  tableS10_value = c(
    total_samples,
    sepsis_samples,
    healthy_samples,
    sofa_available,
    primary_rho,
    primary_p,
    as.numeric(
      sofa_binary$n_case[1]
    ),
    as.numeric(
      sofa_binary$n_control[1]
    ),
    as.numeric(
      sofa_binary$median_case[1]
    ),
    as.numeric(
      sofa_binary$median_control[1]
    ),
    as.numeric(
      sofa_binary$p_value[1]
    ),
    as.numeric(
      sofa_binary$BH_secondary[1]
    ),
    as.numeric(
      sofa_binary$AUC[1]
    ),
    as.numeric(
      sofa_binary$CI_low[1]
    ),
    as.numeric(
      sofa_binary$CI_high[1]
    ),
    n_gene_direction_concordant,
    n_gene_BH_significant,
    sofa_available,
    adjusted_beta,
    adjusted_se
  ),
  
  stringsAsFactors = FALSE
)


fig5_consistency$figure5_value <- vapply(
  
  fig5_consistency$metric,
  
  fig5_value,
  
  numeric(1)
)


fig5_consistency$absolute_difference <-
  abs(
    fig5_consistency$tableS10_value -
      fig5_consistency$figure5_value
  )


fig5_consistency$pass <-
  fig5_consistency$absolute_difference <
  1e-10


cat("\n## MAIN FIGURE 5 CONSISTENCY AUDIT\n\n")

print(
  fig5_consistency,
  row.names = FALSE
)


if (!all(fig5_consistency$pass)) {
  
  stop(
    "Supplementary Table S10 is not numerically consistent with frozen Main Figure 5."
  )
}


cat(
  "\nTable S10 is numerically consistent with frozen Main Figure 5.\n"
)


# =============================================================================
# 23. PUBLICATION-FACING EXTERNAL SUMMARY
# =============================================================================

external_summary <- data.frame(
  
  analysis_role = c(
    
    "Primary external severity endpoint",
    
    "Secondary severity analysis",
    
    "Secondary outcome-associated analysis",
    
    "Secondary collection-site analysis",
    
    "Contextual disease-state comparison",
    
    "Covariate-adjusted severity analysis",
    
    if (
      location_ci_available
    ) {
      "Descriptive location sensitivity"
    } else {
      NULL
    }
  ),
  
  analysis = c(
    
    "Five-gene score versus continuous 24-h SOFA",
    
    "SOFA >=2 versus SOFA 0-1",
    
    "In-hospital death versus survival",
    
    "ICU versus Emergency Room",
    
    "Sepsis versus healthy controls",
    
    "Five-gene score ~ SOFA + age + sex + collection location",
    
    if (
      location_ci_available
    ) {
      "Fixed-effect descriptive pooling of five within-dataset location-specific correlations"
    } else {
      NULL
    }
  ),
  
  population = c(
    
    "Sepsis with SOFA available",
    
    "Sepsis with SOFA available",
    
    "Sepsis with mortality available",
    
    "Sepsis",
    
    "Sepsis and healthy controls",
    
    "Sepsis with complete model covariates",
    
    if (
      location_ci_available
    ) {
      "Five collection locations within GSE185263"
    } else {
      NULL
    }
  ),
  
  n = c(
    
    primary_n,
    
    primary_n,
    
    sum(
      as.numeric(
        mortality_test$n_case[1]
      ),
      as.numeric(
        mortality_test$n_control[1]
      )
    ),
    
    sum(
      as.numeric(
        icu_test$n_case[1]
      ),
      as.numeric(
        icu_test$n_control[1]
      )
    ),
    
    sum(
      as.numeric(
        sepsis_healthy$n_case[1]
      ),
      as.numeric(
        sepsis_healthy$n_control[1]
      )
    ),
    
    sofa_available,
    
    if (
      location_ci_available
    ) {
      as.numeric(
        location_pooled_fixed$total_n[1]
      )
    } else {
      NULL
    }
  ),
  
  n_case = c(
    
    NA,
    
    as.numeric(
      sofa_binary$n_case[1]
    ),
    
    as.numeric(
      mortality_test$n_case[1]
    ),
    
    as.numeric(
      icu_test$n_case[1]
    ),
    
    as.numeric(
      sepsis_healthy$n_case[1]
    ),
    
    NA,
    
    if (
      location_ci_available
    ) {
      NA
    } else {
      NULL
    }
  ),
  
  n_control = c(
    
    NA,
    
    as.numeric(
      sofa_binary$n_control[1]
    ),
    
    as.numeric(
      mortality_test$n_control[1]
    ),
    
    as.numeric(
      icu_test$n_control[1]
    ),
    
    as.numeric(
      sepsis_healthy$n_control[1]
    ),
    
    NA,
    
    if (
      location_ci_available
    ) {
      NA
    } else {
      NULL
    }
  ),
  
  median_case = c(
    
    NA,
    
    as.numeric(
      sofa_binary$median_case[1]
    ),
    
    as.numeric(
      mortality_test$median_case[1]
    ),
    
    as.numeric(
      icu_test$median_case[1]
    ),
    
    as.numeric(
      sepsis_healthy$median_case[1]
    ),
    
    NA,
    
    if (
      location_ci_available
    ) {
      NA
    } else {
      NULL
    }
  ),
  
  median_control = c(
    
    NA,
    
    as.numeric(
      sofa_binary$median_control[1]
    ),
    
    as.numeric(
      mortality_test$median_control[1]
    ),
    
    as.numeric(
      icu_test$median_control[1]
    ),
    
    as.numeric(
      sepsis_healthy$median_control[1]
    ),
    
    NA,
    
    if (
      location_ci_available
    ) {
      NA
    } else {
      NULL
    }
  ),
  
  statistic = c(
    
    "Spearman rho",
    
    "AUC",
    
    "AUC",
    
    "AUC",
    
    "AUC",
    
    "Regression beta per one SOFA point",
    
    if (
      location_ci_available
    ) {
      "Descriptive fixed-effect pooled rho"
    } else {
      NULL
    }
  ),
  
  estimate = c(
    
    primary_rho,
    
    as.numeric(
      sofa_binary$AUC[1]
    ),
    
    as.numeric(
      mortality_test$AUC[1]
    ),
    
    as.numeric(
      icu_test$AUC[1]
    ),
    
    as.numeric(
      sepsis_healthy$AUC[1]
    ),
    
    adjusted_beta,
    
    if (
      location_ci_available
    ) {
      as.numeric(
        location_pooled_fixed$fixed_rho[1]
      )
    } else {
      NULL
    }
  ),
  
  SE = c(
    
    NA,
    
    NA,
    
    NA,
    
    NA,
    
    NA,
    
    adjusted_se,
    
    if (
      location_ci_available
    ) {
      NA
    } else {
      NULL
    }
  ),
  
  CI_low = c(
    
    NA,
    
    as.numeric(
      sofa_binary$CI_low[1]
    ),
    
    as.numeric(
      mortality_test$CI_low[1]
    ),
    
    as.numeric(
      icu_test$CI_low[1]
    ),
    
    as.numeric(
      sepsis_healthy$CI_low[1]
    ),
    
    adjusted_ci_low,
    
    if (
      location_ci_available
    ) {
      as.numeric(
        location_pooled_fixed$fixed_CI_low[1]
      )
    } else {
      NULL
    }
  ),
  
  CI_high = c(
    
    NA,
    
    as.numeric(
      sofa_binary$CI_high[1]
    ),
    
    as.numeric(
      mortality_test$CI_high[1]
    ),
    
    as.numeric(
      icu_test$CI_high[1]
    ),
    
    as.numeric(
      sepsis_healthy$CI_high[1]
    ),
    
    adjusted_ci_high,
    
    if (
      location_ci_available
    ) {
      as.numeric(
        location_pooled_fixed$fixed_CI_high[1]
      )
    } else {
      NULL
    }
  ),
  
  p_value = c(
    
    primary_p,
    
    as.numeric(
      sofa_binary$p_value[1]
    ),
    
    as.numeric(
      mortality_test$p_value[1]
    ),
    
    as.numeric(
      icu_test$p_value[1]
    ),
    
    as.numeric(
      sepsis_healthy$p_value[1]
    ),
    
    adjusted_p,
    
    if (
      location_ci_available
    ) {
      as.numeric(
        location_pooled_fixed$fixed_p[1]
      )
    } else {
      NULL
    }
  ),
  
  BH_adjusted_p = c(
    
    NA,
    
    as.numeric(
      sofa_binary$BH_secondary[1]
    ),
    
    as.numeric(
      mortality_test$BH_secondary[1]
    ),
    
    as.numeric(
      icu_test$BH_secondary[1]
    ),
    
    as.numeric(
      sepsis_healthy$BH_secondary[1]
    ),
    
    NA,
    
    if (
      location_ci_available
    ) {
      NA
    } else {
      NULL
    }
  ),
  
  interpretation = c(
    
    "Prespecified primary external endpoint; positive association with organ-dysfunction severity.",
    
    "Secondary binary severity analysis; moderate discrimination between SOFA strata.",
    
    "Secondary association with in-hospital mortality; not a validated prognostic endpoint.",
    
    "Secondary collection-site association; ICU and Emergency Room groups are not independent cohorts.",
    
    "Contextual comparison only; healthy-control discrimination is not the primary clinically relevant external endpoint.",
    
    "SOFA association persisted after adjustment for age, sex, and collection location.",
    
    if (
      location_ci_available
    ) {
      "Descriptive sensitivity analysis across locations within one GEO dataset; locations must not be presented as independent validation cohorts."
    } else {
      NULL
    }
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 24. LOCATION-SENSITIVITY INTERPRETATION
# =============================================================================

if (location_ci_available) {
  
  all_location_directions_positive <-
    all(
      location_ci$rho >
        0
    )
  
  
  n_locations <- nrow(
    location_ci
  )
  
  
  n_location_BH_sig <- sum(
    location_ci$BH_location <
      0.05,
    na.rm = TRUE
  )
  
  
} else {
  
  all_location_directions_positive <- NA
  
  n_locations <- NA
  
  n_location_BH_sig <- NA
}


# =============================================================================
# 25. README TABLE
# =============================================================================

readme <- data.frame(
  
  Item = c(
    
    "Table",
    
    "Dataset",
    
    "Dataset size",
    
    "Primary external population",
    
    "SOFA-complete sepsis samples",
    
    "Primary frozen panel",
    
    "Score definition",
    
    "Score standardization",
    
    "Primary endpoint",
    
    "Primary result",
    
    "Component-gene replication",
    
    "Secondary analyses",
    
    "Adjusted model",
    
    "Location sensitivity",
    
    "Contextual healthy-control comparison",
    
    "Feature selection in external dataset",
    
    "Coefficient refitting in external dataset",
    
    "Cutoff optimization in external dataset",
    
    "Direction flipping in external dataset",
    
    "Interpretive boundary",
    
    "Source provenance"
  ),
  
  Details = c(
    
    "Supplementary Table S10. External evaluation of the frozen five-gene host-response score in GSE185263.",
    
    "GSE185263 whole-blood RNA-seq dataset.",
    
    paste0(
      total_samples,
      " samples: ",
      sepsis_samples,
      " sepsis and ",
      healthy_samples,
      " healthy."
    ),
    
    paste0(
      sepsis_samples,
      " sepsis samples."
    ),
    
    as.character(
      sofa_available
    ),
    
    primary_panel,
    
    primary_score_definition,
    
    paste0(
      "For the primary severity analysis, gene-wise z-score parameters were estimated across all ",
      sepsis_samples,
      " sepsis samples independently of SOFA availability."
    ),
    
    primary_endpoint,
    
    paste0(
      "n=",
      primary_n,
      "; Spearman rho=",
      format(
        primary_rho,
        digits = 7
      ),
      "; P=",
      format(
        primary_p,
        scientific = TRUE,
        digits = 7
      ),
      "."
    ),
    
    paste0(
      n_gene_direction_concordant,
      "/5 genes showed the expected direction and ",
      n_gene_BH_significant,
      "/5 remained significant after BH correction across the five component genes."
    ),
    
    "Frozen secondary analyses include SOFA >=2 versus 0-1, in-hospital mortality, ICU versus Emergency Room, and contextual sepsis-versus-healthy discrimination.",
    
    paste0(
      "Linear model: five-gene score ~ SOFA + age + sex + collection location. Frozen SOFA beta=",
      format(
        adjusted_beta,
        digits = 7
      ),
      ", SE=",
      format(
        adjusted_se,
        digits = 7
      ),
      ", P=",
      format(
        adjusted_p,
        scientific = TRUE,
        digits = 7
      ),
      if (
        is.finite(
          adjusted_r2
        )
      ) {
        paste0(
          ", model R-squared=",
          format(
            adjusted_r2,
            digits = 7
          )
        )
      } else {
        ""
      },
      "."
    ),
    
    if (
      location_ci_available
    ) {
      paste0(
        n_locations,
        " collection locations were examined within GSE185263. All location-specific rho estimates were positive; ",
        n_location_BH_sig,
        "/",
        n_locations,
        " were BH-significant. These are within-dataset sensitivity analyses, not independent validation cohorts."
      )
    } else {
      "Location-specific analyses are within-dataset sensitivity analyses, not independent validation cohorts."
    },
    
    "Sepsis-versus-healthy discrimination is contextual only and must not be presented as the primary clinically relevant validation endpoint.",
    
    "NO",
    
    "NO",
    
    "NO",
    
    "NO",
    
    "The five-gene signature is interpreted as a molecular readout of host-response state and severity, not as a clinically validated diagnostic or prognostic assay.",
    
    "All inferential quantities are copied from frozen Scripts 142b, 143, and 150; Script 165 performs packaging and numerical provenance auditing only."
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 26. INTERPRETATION TABLE
# =============================================================================

interpretation_table <- data.frame(
  
  Topic = c(
    
    "Primary finding",
    
    "Binary SOFA analysis",
    
    "Component genes",
    
    "Adjusted association",
    
    "Mortality",
    
    "Collection site",
    
    "Healthy controls",
    
    "Location-specific analyses",
    
    "Pooled location estimate",
    
    "Clinical interpretation"
  ),
  
  Interpretation = c(
    
    "The frozen five-gene score was positively associated with continuous 24-h SOFA in the independent GSE185263 dataset.",
    
    "Higher scores were observed in patients with SOFA >=2 than in those with SOFA 0-1; this analysis is secondary to the prespecified continuous-SOFA endpoint.",
    
    "All five component genes associated with SOFA in their prespecified direction and all five remained significant after BH correction across the five genes.",
    
    "The association between the score and SOFA persisted after adjustment for age, sex, and collection location.",
    
    "Higher scores among patients who died support an outcome-related association but do not constitute validation of a prognostic assay.",
    
    "Higher scores in ICU than Emergency Room samples are interpreted as additional severity/context information rather than independent cohort validation.",
    
    "Sepsis-versus-healthy discrimination is contextual because the main clinical question concerns variation in host-response severity within sepsis.",
    
    "Collection locations are strata within a single GEO dataset and must not be described as independent external validation cohorts.",
    
    "The fixed-effect pooled location estimate is descriptive only; no random-effects or heterogeneity result is required for the manuscript's primary claim.",
    
    "GSE185263 supports external molecular replication of the association between the frozen five-gene host-response score and organ-dysfunction severity; it does not establish a calibrated clinical diagnostic or prognostic test."
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# 27. CREATE SUBMISSION WORKBOOK
# =============================================================================

output_workbook <- file.path(
  tables_dir,
  "165_TableS10_GSE185263_external_evaluation.xlsx"
)


wb <- openxlsx::createWorkbook()


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


scientific_style <- openxlsx::createStyle(
  
  numFmt = "0.000E+00"
)


decimal_style <- openxlsx::createStyle(
  
  numFmt = "0.000000"
)


write_publication_sheet <- function(
    workbook,
    sheet_name,
    data,
    freeze_first_row = TRUE
) {
  
  openxlsx::addWorksheet(
    workbook,
    sheet_name
  )
  
  
  openxlsx::writeData(
    workbook,
    sheet_name,
    data,
    headerStyle = header_style
  )
  
  
  if (
    nrow(data) >
    0 &&
    ncol(data) >
    0
  ) {
    
    openxlsx::addStyle(
      workbook,
      sheet_name,
      wrap_style,
      rows = seq_len(
        nrow(data)
      ) + 1,
      cols = seq_len(
        ncol(data)
      ),
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  
  if (freeze_first_row) {
    
    openxlsx::freezePane(
      workbook,
      sheet_name,
      firstActiveRow = 2
    )
  }
  
  
  openxlsx::setColWidths(
    workbook,
    sheet_name,
    cols = seq_len(
      ncol(data)
    ),
    widths = "auto"
  )
}


write_publication_sheet(
  wb,
  "S10_ReadMe",
  readme
)


write_publication_sheet(
  wb,
  "External_summary",
  external_summary
)


write_publication_sheet(
  wb,
  "Primary_SOFA",
  primary
)


write_publication_sheet(
  wb,
  "Secondary_associations",
  secondary
)


write_publication_sheet(
  wb,
  "Component_gene_SOFA",
  gene_sofa
)


write_publication_sheet(
  wb,
  "Adjusted_model",
  adjusted_model
)


write_publication_sheet(
  wb,
  "Adjusted_summary",
  adjusted_summary
)


write_publication_sheet(
  wb,
  "Location_SOFA",
  location_sofa
)


write_publication_sheet(
  wb,
  "Location_pooled",
  location_pooled
)


if (location_ci_available) {
  
  write_publication_sheet(
    wb,
    "Location_SOFA_CI",
    location_ci
  )
  
  
  write_publication_sheet(
    wb,
    "Location_pooled_fixed",
    location_pooled_fixed
  )
}


write_publication_sheet(
  wb,
  "Scaling_sensitivity",
  scaling
)


write_publication_sheet(
  wb,
  "SOFA_summary",
  sofa_summary
)


write_publication_sheet(
  wb,
  "SOFA_groups",
  sofa_groups
)


write_publication_sheet(
  wb,
  "Mortality_context",
  mortality
)


write_publication_sheet(
  wb,
  "Site_summary",
  site_summary
)


write_publication_sheet(
  wb,
  "Location_summary",
  location_summary
)


write_publication_sheet(
  wb,
  "Frozen_gene_coverage",
  gene_coverage
)


write_publication_sheet(
  wb,
  "Sample_scores",
  sample_scores
)


write_publication_sheet(
  wb,
  "Interpretation",
  interpretation_table
)


openxlsx::setColWidths(
  wb,
  "S10_ReadMe",
  cols = 1,
  widths = 32
)


openxlsx::setColWidths(
  wb,
  "S10_ReadMe",
  cols = 2,
  widths = 90
)


openxlsx::setColWidths(
  wb,
  "Interpretation",
  cols = 1,
  widths = 30
)


openxlsx::setColWidths(
  wb,
  "Interpretation",
  cols = 2,
  widths = 95
)


openxlsx::saveWorkbook(
  wb,
  output_workbook,
  overwrite = TRUE
)


# =============================================================================
# 28. INTERNAL AUDIT WORKBOOK
# =============================================================================

audit_workbook <- file.path(
  audit_dir,
  "165_INTERNAL_AUDIT_TableS10_GSE185263.xlsx"
)


audit_wb <- openxlsx::createWorkbook()


write_publication_sheet(
  audit_wb,
  "Schema_audit",
  schema_audit
)


write_publication_sheet(
  audit_wb,
  "CSV_XLSX_equivalence",
  equivalence_audit
)


write_publication_sheet(
  audit_wb,
  "Run_guardrails",
  run_guardrails
)


write_publication_sheet(
  audit_wb,
  "Cohort_audit",
  cohort_audit
)


write_publication_sheet(
  audit_wb,
  "Numerical_anchors",
  anchor_audit
)


write_publication_sheet(
  audit_wb,
  "Figure5_consistency",
  fig5_consistency
)


if (location_ci_available) {
  
  write_publication_sheet(
    audit_wb,
    "Location_142b_vs_143",
    location_join_audit
  )
}


write_publication_sheet(
  audit_wb,
  "Figure5_frozen_audit",
  fig5_audit
)


openxlsx::saveWorkbook(
  audit_wb,
  audit_workbook,
  overwrite = TRUE
)


# =============================================================================
# 29. TABLE TITLE AND NOTE
# =============================================================================

table_title_note <- c(
  
  "Supplementary Table S10. External evaluation of the frozen five-gene host-response score in GSE185263.",
  
  "",
  
  paste0(
    "The prespecified primary external endpoint was the Spearman association between the frozen five-gene score and continuous 24-h SOFA among patients with sepsis. The score was defined as ",
    primary_score_definition,
    "."
  ),
  
  "",
  
  paste0(
    "Gene-wise z-score parameters for the primary severity analysis were estimated across all ",
    sepsis_samples,
    " sepsis samples independently of SOFA availability; ",
    sofa_available,
    " sepsis samples had available 24-h SOFA values."
  ),
  
  "",
  
  "Secondary analyses included comparison of SOFA >=2 versus SOFA 0-1, in-hospital mortality, ICU versus Emergency Room collection site, and contextual sepsis-versus-healthy discrimination.",
  
  "",
  
  "The covariate-adjusted model included SOFA, age, sex, and collection location.",
  
  "",
  
  "Component-gene P values were adjusted across the five constituent genes using the Benjamini-Hochberg procedure. Secondary score-comparison P values retain the frozen multiplicity adjustment generated by Script 142b.",
  
  "",
  
  "Collection-location analyses represent sensitivity analyses within a single GEO dataset and must not be interpreted as independent external validation cohorts. The fixed-effect pooled location estimate is descriptive only.",
  
  "",
  
  "No feature reselection, coefficient refitting, cutoff optimization, score-direction reversal, or endotype reconstruction was performed in the external dataset.",
  
  "",
  
  "The sepsis-versus-healthy analysis is contextual and is not the primary clinically relevant external validation endpoint.",
  
  "",
  
  "A numerically reported P value of 0 in the frozen contextual sepsis-versus-healthy table reflects computational precision and should not be interpreted as a literally zero probability.",
  
  "",
  
  "The five-gene signature is interpreted as a molecular readout of host-response state and organ-dysfunction severity rather than as a clinically validated diagnostic or prognostic assay."
)


title_note_file <- file.path(
  text_dir,
  "165_TableS10_title_and_note_EN.txt"
)


writeLines(
  table_title_note,
  title_note_file
)


# =============================================================================
# 30. PROPOSED RESULTS 3.9
# =============================================================================

r2_sentence <- if (
  is.finite(
    adjusted_r2
  )
) {
  
  paste0(
    " The complete adjusted model had R-squared = ",
    format(
      adjusted_r2,
      digits = 3
    ),
    "."
  )
  
} else {
  
  ""
}


location_sentence <- if (
  location_ci_available
) {
  
  paste0(
    " Across the five collection locations, all location-specific correlations were positive, although their individual statistical precision varied. A descriptive fixed-effect summary yielded rho = ",
    format(
      as.numeric(
        location_pooled_fixed$fixed_rho[1]
      ),
      digits = 3
    ),
    " (95% CI ",
    format(
      as.numeric(
        location_pooled_fixed$fixed_CI_low[1]
      ),
      digits = 3
    ),
    " to ",
    format(
      as.numeric(
        location_pooled_fixed$fixed_CI_high[1]
      ),
      digits = 3
    ),
    "; P = ",
    format(
      as.numeric(
        location_pooled_fixed$fixed_p[1]
      ),
      scientific = TRUE,
      digits = 3
    ),
    "). These locations represent strata within a single dataset rather than independent validation cohorts."
  )
  
} else {
  
  ""
}


results_3_9 <- paste0(
  
  "### 3.9 External replication of the five-gene host-response score with organ-dysfunction severity in GSE185263\n\n",
  
  "The association between the frozen five-gene host-response score and organ-dysfunction severity was evaluated in the independent GSE185263 whole-blood RNA-seq dataset. Of 348 sepsis samples, 345 had available 24-h SOFA measurements. The prespecified primary analysis showed a positive association between the five-gene score and continuous SOFA (Spearman rho = ",
  
  format(
    primary_rho,
    digits = 3
  ),
  
  ", P = ",
  
  format(
    primary_p,
    scientific = TRUE,
    digits = 3
  ),
  
  "; Fig. 5 and Supplementary Table S10). ",
  
  "Consistent with this continuous association, patients with SOFA >=2 (n = ",
  
  as.numeric(
    sofa_binary$n_case[1]
  ),
  
  ") had higher scores than those with SOFA 0-1 (n = ",
  
  as.numeric(
    sofa_binary$n_control[1]
  ),
  
  "; median ",
  
  format(
    as.numeric(
      sofa_binary$median_case[1]
    ),
    digits = 3
  ),
  
  " versus ",
  
  format(
    as.numeric(
      sofa_binary$median_control[1]
    ),
    digits = 3
  ),
  
  "; P = ",
  
  format(
    as.numeric(
      sofa_binary$p_value[1]
    ),
    scientific = TRUE,
    digits = 3
  ),
  
  ", BH-adjusted P = ",
  
  format(
    as.numeric(
      sofa_binary$BH_secondary[1]
    ),
    scientific = TRUE,
    digits = 3
  ),
  
  "), with an AUC of ",
  
  format(
    as.numeric(
      sofa_binary$AUC[1]
    ),
    digits = 3
  ),
  
  " (95% CI ",
  
  format(
    as.numeric(
      sofa_binary$CI_low[1]
    ),
    digits = 3
  ),
  
  " to ",
  
  format(
    as.numeric(
      sofa_binary$CI_high[1]
    ),
    digits = 3
  ),
  
  "). All five component genes were associated with SOFA in their prespecified direction, and all five remained significant after Benjamini-Hochberg correction across the component genes (Supplementary Table S10). ",
  
  "The score-SOFA association also persisted in a linear model adjusted for age, sex, and collection location (beta = ",
  
  format(
    adjusted_beta,
    digits = 4
  ),
  
  " per SOFA point, SE = ",
  
  format(
    adjusted_se,
    digits = 3
  ),
  
  ", P = ",
  
  format(
    adjusted_p,
    scientific = TRUE,
    digits = 3
  ),
  
  ").",
  
  r2_sentence,
  
  " Secondary analyses showed higher scores among patients who died than among survivors (n = ",
  
  as.numeric(
    mortality_test$n_case[1]
  ),
  
  " versus ",
  
  as.numeric(
    mortality_test$n_control[1]
  ),
  
  "; AUC = ",
  
  format(
    as.numeric(
      mortality_test$AUC[1]
    ),
    digits = 3
  ),
  
  ", P = ",
  
  format(
    as.numeric(
      mortality_test$p_value[1]
    ),
    scientific = TRUE,
    digits = 3
  ),
  
  ") and among ICU versus Emergency Room samples (n = ",
  
  as.numeric(
    icu_test$n_case[1]
  ),
  
  " versus ",
  
  as.numeric(
    icu_test$n_control[1]
  ),
  
  "; AUC = ",
  
  format(
    as.numeric(
      icu_test$AUC[1]
    ),
    digits = 3
  ),
  
  ", P = ",
  
  format(
    as.numeric(
      icu_test$p_value[1]
    ),
    scientific = TRUE,
    digits = 3
  ),
  
  ").",
  
  location_sentence,
  
  " Contextual discrimination between sepsis and healthy controls was high (AUC = ",
  
  format(
    as.numeric(
      sepsis_healthy$AUC[1]
    ),
    digits = 3
  ),
  
  "), but this comparison was not the primary clinically relevant external endpoint. Together, these findings independently reproduce the association between the frozen five-gene host-response score and organ-dysfunction severity, while supporting its interpretation as a molecular host-response index rather than a calibrated diagnostic or prognostic assay."
)


results_file <- file.path(
  text_dir,
  "165_proposed_Results_3.9_GSE185263_external_severity_EN.txt"
)


writeLines(
  results_3_9,
  results_file
)


# =============================================================================
# 31. FINAL REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 165 FINAL v2 completed successfully.\n")
cat("====================================================================\n\n")


cat("## GSE185263 COHORT\n\n")

cat(
  "Total samples = ",
  total_samples,
  "\n",
  sep = ""
)


cat(
  "Sepsis = ",
  sepsis_samples,
  "\n",
  sep = ""
)


cat(
  "Healthy = ",
  healthy_samples,
  "\n",
  sep = ""
)


cat(
  "SOFA-complete sepsis = ",
  sofa_available,
  "\n",
  sep = ""
)


cat(
  "Mortality-complete sepsis = ",
  mortality_available,
  "\n"
)


cat("\n## PRIMARY EXTERNAL SEVERITY ENDPOINT\n\n")

cat(
  "n = ",
  primary_n,
  "\n",
  sep = ""
)


cat(
  "Spearman rho = ",
  primary_rho,
  "\n",
  sep = ""
)


cat(
  "P = ",
  primary_p,
  "\n",
  sep = ""
)


cat("\n## SOFA >=2 vs SOFA 0-1\n\n")

cat(
  "n = ",
  sofa_binary$n_case[1],
  " vs ",
  sofa_binary$n_control[1],
  "\n",
  sep = ""
)


cat(
  "Median score = ",
  sofa_binary$median_case[1],
  " vs ",
  sofa_binary$median_control[1],
  "\n",
  sep = ""
)


cat(
  "P = ",
  sofa_binary$p_value[1],
  "\n",
  sep = ""
)


cat(
  "BH q = ",
  sofa_binary$BH_secondary[1],
  "\n",
  sep = ""
)


cat(
  "AUC = ",
  sofa_binary$AUC[1],
  "\n",
  sep = ""
)


cat(
  "95% CI = ",
  sofa_binary$CI_low[1],
  " to ",
  sofa_binary$CI_high[1],
  "\n",
  sep = ""
)


cat("\n## COMPONENT-GENE SOFA REPLICATION\n\n")

cat(
  "Expected-direction concordance = ",
  n_gene_direction_concordant,
  "/5\n",
  sep = ""
)


cat(
  "BH-adjusted P < 0.05 = ",
  n_gene_BH_significant,
  "/5\n",
  sep = ""
)


cat("\n## ADJUSTED SOFA ASSOCIATION\n\n")

cat(
  "Beta = ",
  adjusted_beta,
  "\n",
  sep = ""
)


cat(
  "SE = ",
  adjusted_se,
  "\n",
  sep = ""
)


cat(
  "P = ",
  adjusted_p,
  "\n",
  sep = ""
)


if (is.finite(adjusted_r2)) {
  
  cat(
    "Model R-squared = ",
    adjusted_r2,
    "\n",
    sep = ""
  )
}


cat("\n## SECONDARY CONTEXT\n\n")

cat(
  "Mortality AUC = ",
  mortality_test$AUC[1],
  "; P = ",
  mortality_test$p_value[1],
  "\n",
  sep = ""
)


cat(
  "ICU vs Emergency Room AUC = ",
  icu_test$AUC[1],
  "; P = ",
  icu_test$p_value[1],
  "\n",
  sep = ""
)


cat(
  "Contextual sepsis vs healthy AUC = ",
  sepsis_healthy$AUC[1],
  "\n",
  sep = ""
)


if (location_ci_available) {
  
  cat("\n## LOCATION SENSITIVITY\n\n")
  
  print(
    location_ci,
    row.names = FALSE
  )
  
  
  cat(
    "\nAll five location-specific rho estimates positive = ",
    all_location_directions_positive,
    "\n",
    sep = ""
  )
  
  
  cat(
    "Descriptive fixed-effect pooled rho = ",
    location_pooled_fixed$fixed_rho[1],
    "\n",
    sep = ""
  )
  
  
  cat(
    "95% CI = ",
    location_pooled_fixed$fixed_CI_low[1],
    " to ",
    location_pooled_fixed$fixed_CI_high[1],
    "\n",
    sep = ""
  )
  
  
  cat(
    "P = ",
    location_pooled_fixed$fixed_p[1],
    "\n",
    sep = ""
  )
}


cat("\n## OUTPUT FILES\n\n")


cat(
  "Supplementary Table S10:\n"
)


cat(
  normalizePath(
    output_workbook,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n"
)


cat(
  "Internal audit:\n"
)


cat(
  normalizePath(
    audit_workbook,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n"
)


cat(
  "Table title/note:\n"
)


cat(
  normalizePath(
    title_note_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n\n"
)


cat(
  "Proposed Results 3.9:\n"
)


cat(
  normalizePath(
    results_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


cat("\n## REPORTING GUARDRAILS\n\n")

cat(
  "- All inferential statistics come from frozen Scripts 142b, 143, and 150.\n"
)

cat(
  "- Script 165 performs no new inferential statistical analysis.\n"
)

cat(
  "- The prespecified primary endpoint is continuous 24-h SOFA in sepsis.\n"
)

cat(
  "- Gene-wise z parameters were estimated across all 348 sepsis samples independently of SOFA availability.\n"
)

cat(
  "- All five component genes show the expected SOFA direction and all five are BH-significant.\n"
)

cat(
  "- The adjusted SOFA association persists after age, sex, and collection-location adjustment.\n"
)

cat(
  "- Mortality and ICU-versus-Emergency-Room analyses are secondary.\n"
)

cat(
  "- Sepsis-versus-healthy discrimination is contextual, not the primary clinically relevant validation endpoint.\n"
)

cat(
  "- The five collection locations are sensitivity strata within one dataset, not five independent validation cohorts.\n"
)

cat(
  "- The fixed-effect pooled location estimate is descriptive only.\n"
)

cat(
  "- No random-effects location meta-analysis is required for the main manuscript claim.\n"
)

cat(
  "- GSE185263 supports external replication of an organ-dysfunction-severity association, not calibration of a clinical diagnostic or prognostic assay.\n"
)


cat("\nDone.\n")