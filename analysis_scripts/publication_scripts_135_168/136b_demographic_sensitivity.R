# ==============================================================================
# Script 136b
# Demographic integration and age/sex sensitivity analysis
#
# Project: Sepsis_DESeq2
#
# PURPOSE
# 1. Integrate age and sex from Supplementary Table S1
# 2. Quantify demographic imbalance between BP and BC
# 3. Test age/sex associations with SRS, SRSq, CTS and the 5-gene score
# 4. Recalculate global clinical BH FDR including age and sex
# 5. Test whether the primary five-gene score remains associated with
#    sepsis after adjustment for age and sex
# 6. Test individual primary genes after age/sex adjustment
#
# IMPORTANT
# - blood only
# - no new feature selection
# - no urine
# - no lncRNA
# - no mortality prediction model
# - no diagnostic validation claim
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

script_name <- "136b_demographic_sensitivity.R"
run_date <- Sys.time()

cat("\n")
cat("====================================================================\n")
cat("Running Script 136b\n")
cat("Age/sex integration and demographic sensitivity analysis\n")
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
  "readxl",
  "openxlsx",
  "ggplot2",
  "tibble"
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
      "Missing packages: ",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(tibble)
})

cat("Packages loaded successfully.\n\n")


# ==============================================================================
# 2. HELPER FUNCTIONS
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
    length(x[["statistic"]]) > 0
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
    length(x[["p.value"]]) > 0
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
    length(x[["estimate"]]) > 0
  ) {
    
    return(
      as.numeric(
        x[["estimate"]][1]
      )
    )
  }
  
  return(NA_real_)
}


make_spearman_row <- function(
    x,
    y,
    framework,
    clinical_variable,
    clinical_label
) {
  
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
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = "continuous_correlation",
        test = "Spearman",
        n = length(x2),
        statistic = NA_real_,
        effect = NA_real_,
        effect_name = "Spearman_rho",
        p_value = NA_real_,
        group_summary = "",
        note = "Insufficient observations"
      )
    )
  }
  
  ht <- suppressWarnings(
    stats::cor.test(
      x = x2,
      y = y2,
      method = "spearman",
      exact = FALSE
    )
  )
  
  return(
    tibble(
      framework = framework,
      clinical_variable = clinical_variable,
      clinical_label = clinical_label,
      test_family = "continuous_correlation",
      test = "Spearman",
      n = length(x2),
      statistic = safe_htest_statistic(ht),
      effect = safe_htest_estimate(ht),
      effect_name = "Spearman_rho",
      p_value = safe_htest_pvalue(ht),
      group_summary = "",
      note = "Demographic extension from Script 136b"
    )
  )
}


make_wilcoxon_row <- function(
    x,
    group,
    framework,
    clinical_variable,
    clinical_label,
    test_family
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
  
  if (nlevels(g2) != 2) {
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = test_family,
        test = "Wilcoxon_rank_sum",
        n = length(x2),
        statistic = NA_real_,
        effect = NA_real_,
        effect_name = NA_character_,
        p_value = NA_real_,
        group_summary = "",
        note = "Expected exactly two groups"
      )
    )
  }
  
  lev <- levels(g2)
  
  x1 <- x2[
    g2 == lev[1]
  ]
  
  x2b <- x2[
    g2 == lev[2]
  ]
  
  ht <- stats::wilcox.test(
    x = x1,
    y = x2b,
    exact = FALSE,
    paired = FALSE
  )
  
  med1 <- median(
    x1,
    na.rm = TRUE
  )
  
  med2 <- median(
    x2b,
    na.rm = TRUE
  )
  
  return(
    tibble(
      framework = framework,
      clinical_variable = clinical_variable,
      clinical_label = clinical_label,
      test_family = test_family,
      test = "Wilcoxon_rank_sum",
      n = length(x2),
      statistic = safe_htest_statistic(ht),
      effect = med2 - med1,
      effect_name = paste0(
        "median_difference_",
        lev[2],
        "_minus_",
        lev[1]
      ),
      p_value = safe_htest_pvalue(ht),
      group_summary = paste0(
        lev[1],
        ": n=",
        length(x1),
        ", median=",
        signif(med1, 4),
        " | ",
        lev[2],
        ": n=",
        length(x2b),
        ", median=",
        signif(med2, 4)
      ),
      note = "Demographic extension from Script 136b"
    )
  )
}


make_kruskal_row <- function(
    x,
    group,
    framework,
    clinical_variable,
    clinical_label
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
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = "continuous_by_endotype",
        test = "Kruskal_Wallis",
        n = length(x2),
        statistic = NA_real_,
        effect = NA_real_,
        effect_name = "epsilon_squared",
        p_value = NA_real_,
        group_summary = "",
        note = "Insufficient groups"
      )
    )
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
  
  epsilon2 <- max(
    0,
    epsilon2
  )
  
  return(
    tibble(
      framework = framework,
      clinical_variable = clinical_variable,
      clinical_label = clinical_label,
      test_family = "continuous_by_endotype",
      test = "Kruskal_Wallis",
      n = n_total,
      statistic = H,
      effect = epsilon2,
      effect_name = "epsilon_squared",
      p_value = safe_htest_pvalue(ht),
      group_summary = "",
      note = "Demographic extension from Script 136b"
    )
  )
}


make_fisher_row <- function(
    molecular_group,
    clinical_group,
    framework,
    clinical_variable,
    clinical_label
) {
  
  keep <- complete.cases(
    molecular_group,
    clinical_group
  )
  
  g1 <- droplevels(
    factor(
      molecular_group[keep]
    )
  )
  
  g2 <- droplevels(
    factor(
      clinical_group[keep]
    )
  )
  
  if (
    length(g1) < 6 ||
    nlevels(g1) < 2 ||
    nlevels(g2) < 2
  ) {
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = "categorical_endotype",
        test = "Fisher_exact",
        n = length(g1),
        statistic = NA_real_,
        effect = NA_real_,
        effect_name = NA_character_,
        p_value = NA_real_,
        group_summary = "",
        note = "Insufficient groups"
      )
    )
  }
  
  tab <- table(
    g1,
    g2
  )
  
  ht <- stats::fisher.test(
    tab
  )
  
  return(
    tibble(
      framework = framework,
      clinical_variable = clinical_variable,
      clinical_label = clinical_label,
      test_family = "categorical_endotype",
      test = "Fisher_exact",
      n = sum(tab),
      statistic = NA_real_,
      effect = NA_real_,
      effect_name = NA_character_,
      p_value = safe_htest_pvalue(ht),
      group_summary = paste(
        capture.output(
          print(tab)
        ),
        collapse = " "
      ),
      note = "Demographic extension from Script 136b"
    )
  )
}


# ==============================================================================
# 3. INPUT FILES
# ==============================================================================

file_name_S1 <-
  "Table_S1_Deidentified_RNAseq_Sample_Annotation_REVISED.xlsx"


user_home <- Sys.getenv(
  "USERPROFILE"
)


S1_candidates <- c(
  
  file.path(
    "data",
    file_name_S1
  ),
  
  file.path(
    getwd(),
    file_name_S1
  ),
  
  file.path(
    user_home,
    "Downloads",
    file_name_S1
  ),
  
)


S1_file <- S1_candidates[
  file.exists(
    S1_candidates
  )
]


if (length(S1_file) == 0) {
  
  stop(
    paste0(
      "Supplementary Table S1 not found.\n",
      "Copy the file to:\n",
      file.path(
        getwd(),
        "data",
        file_name_S1
      )
    )
  )
}


S1_file <- S1_file[1]


script136_tests_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136_clinical_associations",
  "tables",
  "136_all_clinical_association_tests.csv"
)


script136_bp_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136_clinical_associations",
  "tables",
  "136_blood_molecular_clinical_merged.csv"
)


script135_scores_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "135_validation",
  "tables",
  "135_blood_scores_with_final_endotypes.csv"
)


script135_expression_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "135_validation",
  "tables",
  "135_selected_gene_expression_logCPM.csv"
)


required_files <- c(
  S1_file,
  script136_tests_file,
  script136_bp_file,
  script135_scores_file,
  script135_expression_file
)


if (any(!file.exists(required_files))) {
  
  missing_files <- required_files[
    !file.exists(
      required_files
    )
  ]
  
  stop(
    paste0(
      "Missing required files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


cat("Input files:\n")

cat(
  "Table S1: ",
  S1_file,
  "\n",
  sep = ""
)

cat(
  "Script 136 tests: ",
  script136_tests_file,
  "\n",
  sep = ""
)

cat(
  "Script 136 BP table: ",
  script136_bp_file,
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
  "Script 135 expression: ",
  script135_expression_file,
  "\n\n",
  sep = ""
)


# ==============================================================================
# 4. OUTPUT DIRECTORIES
# ==============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136b_demographic_sensitivity"
)


tables_dir <- file.path(
  output_dir,
  "tables"
)


figures_dir <- file.path(
  output_dir,
  "figures"
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
# 5. READ PARTICIPANT METADATA
# ==============================================================================

participant_metadata <- readxl::read_excel(
  S1_file,
  sheet = "Participant_Metadata"
)


required_demo_columns <- c(
  "participant_id",
  "cohort",
  "sex",
  "age_years",
  "blood_rnaseq_sample_id"
)


missing_demo_columns <- setdiff(
  required_demo_columns,
  names(
    participant_metadata
  )
)


if (length(missing_demo_columns) > 0) {
  
  stop(
    paste0(
      "Missing Participant_Metadata columns: ",
      paste(
        missing_demo_columns,
        collapse = ", "
      )
    )
  )
}


demographics <- participant_metadata %>%
  transmute(
    
    participant_id =
      as.character(
        participant_id
      ),
    
    sample_id =
      clean_sample_id(
        blood_rnaseq_sample_id
      ),
    
    cohort =
      as.character(
        cohort
      ),
    
    sex =
      factor(
        as.character(
          sex
        ),
        levels = c(
          "Female",
          "Male"
        )
      ),
    
    age_years =
      suppressWarnings(
        as.numeric(
          age_years
        )
      )
    
  ) %>%
  filter(
    grepl(
      "^(BP|BC)[0-9]+$",
      sample_id
    )
  )


if (nrow(demographics) != 45) {
  
  stop(
    paste0(
      "Expected 45 blood participants; observed ",
      nrow(demographics)
    )
  )
}


demographics$condition <- dplyr::case_when(
  
  grepl(
    "^BP[0-9]+$",
    demographics$sample_id
  ) ~ "BP",
  
  grepl(
    "^BC[0-9]+$",
    demographics$sample_id
  ) ~ "BC",
  
  TRUE ~ NA_character_
)


if (any(is.na(demographics$condition))) {
  
  stop(
    "Unable to assign BP/BC condition from sample_id in demographics."
  )
}


demographics$condition <- factor(
  demographics$condition,
  levels = c(
    "BC",
    "BP"
  )
)


cat("Demographic table validation:\n")

cat(
  "Blood participants: ",
  nrow(demographics),
  "\n",
  sep = ""
)

cat(
  "BP: ",
  sum(
    demographics$condition == "BP"
  ),
  "\n",
  sep = ""
)

cat(
  "BC: ",
  sum(
    demographics$condition == "BC"
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 6. DEMOGRAPHIC SUMMARY
# ==============================================================================

demographic_summary <- demographics %>%
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


cat("Demographic summary:\n")

print(
  tibble::as_tibble(
    demographic_summary
  ),
  n = Inf
)

cat("\n")


# ==============================================================================
# 7. BP vs BC DEMOGRAPHIC BALANCE
# ==============================================================================

age_bp <- demographics$age_years[
  demographics$condition == "BP"
]


age_bc <- demographics$age_years[
  demographics$condition == "BC"
]


age_balance_ht <- stats::wilcox.test(
  x = age_bp,
  y = age_bc,
  exact = FALSE,
  paired = FALSE
)


sex_table <- table(
  demographics$condition,
  demographics$sex
)


sex_balance_ht <- stats::fisher.test(
  sex_table
)


demographic_balance <- tibble(
  
  variable = c(
    "Age",
    "Sex"
  ),
  
  comparison = c(
    "BP_vs_BC",
    "BP_vs_BC"
  ),
  
  test = c(
    "Wilcoxon_rank_sum",
    "Fisher_exact"
  ),
  
  statistic = c(
    safe_htest_statistic(
      age_balance_ht
    ),
    NA_real_
  ),
  
  p_value = c(
    safe_htest_pvalue(
      age_balance_ht
    ),
    safe_htest_pvalue(
      sex_balance_ht
    )
  ),
  
  interpretation = c(
    "Age distribution differs between BP and BC",
    "Sex distribution differs between BP and BC"
  )
)


cat("BP vs BC demographic balance:\n")

print(
  tibble::as_tibble(
    demographic_balance
  ),
  n = Inf
)

cat("\n")


cat("Sex table:\n")

print(
  sex_table
)

cat("\n")


# ==============================================================================
# 8. READ SCRIPT 136 BP TABLE
# ==============================================================================

bp136 <- read.csv(
  script136_bp_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


bp136$sample_id <- clean_sample_id(
  bp136$sample_id
)


bp_demo <- bp136 %>%
  left_join(
    
    demographics %>%
      select(
        sample_id,
        age_years,
        sex
      ),
    
    by = "sample_id"
  )


if (nrow(bp_demo) != 35) {
  
  stop(
    "Expected 35 BP rows after demographic merge."
  )
}


cat("Sepsis demographic availability:\n")

cat(
  "Age available: ",
  sum(
    !is.na(
      bp_demo$age_years
    )
  ),
  "/35\n",
  sep = ""
)

cat(
  "Sex available: ",
  sum(
    !is.na(
      bp_demo$sex
    )
  ),
  "/35\n\n",
  sep = ""
)


# ==============================================================================
# 9. AGE ASSOCIATIONS WITHIN SEPSIS
# ==============================================================================

new_tests <- list()


new_tests[[1]] <- make_spearman_row(
  
  x =
    bp_demo$primary_5gene_score,
  
  y =
    bp_demo$age_years,
  
  framework =
    "Primary_5gene_score",
  
  clinical_variable =
    "age_years",
  
  clinical_label =
    "Age"
)


new_tests[[2]] <- make_spearman_row(
  
  x =
    bp_demo$SRSq,
  
  y =
    bp_demo$age_years,
  
  framework =
    "SRSq",
  
  clinical_variable =
    "age_years",
  
  clinical_label =
    "Age"
)


new_tests[[3]] <- make_wilcoxon_row(
  
  x =
    bp_demo$age_years,
  
  group =
    bp_demo$SRS,
  
  framework =
    "SRS_class",
  
  clinical_variable =
    "age_years",
  
  clinical_label =
    "Age",
  
  test_family =
    "continuous_by_endotype"
)


new_tests[[4]] <- make_kruskal_row(
  
  x =
    bp_demo$age_years,
  
  group =
    bp_demo$CTS,
  
  framework =
    "CTS_class",
  
  clinical_variable =
    "age_years",
  
  clinical_label =
    "Age"
)


# ==============================================================================
# 10. SEX ASSOCIATIONS WITHIN SEPSIS
# ==============================================================================

new_tests[[5]] <- make_fisher_row(
  
  molecular_group =
    bp_demo$SRS,
  
  clinical_group =
    bp_demo$sex,
  
  framework =
    "SRS_class",
  
  clinical_variable =
    "sex",
  
  clinical_label =
    "Sex"
)


new_tests[[6]] <- make_fisher_row(
  
  molecular_group =
    bp_demo$CTS,
  
  clinical_group =
    bp_demo$sex,
  
  framework =
    "CTS_class",
  
  clinical_variable =
    "sex",
  
  clinical_label =
    "Sex"
)


new_tests[[7]] <- make_wilcoxon_row(
  
  x =
    bp_demo$primary_5gene_score,
  
  group =
    bp_demo$sex,
  
  framework =
    "Primary_5gene_score",
  
  clinical_variable =
    "sex",
  
  clinical_label =
    "Sex",
  
  test_family =
    "molecular_score_by_clinical_group"
)


new_tests[[8]] <- make_wilcoxon_row(
  
  x =
    bp_demo$SRSq,
  
  group =
    bp_demo$sex,
  
  framework =
    "SRSq",
  
  clinical_variable =
    "sex",
  
  clinical_label =
    "Sex",
  
  test_family =
    "molecular_score_by_clinical_group"
)


demographic_association_tests <- bind_rows(
  new_tests
)


cat("Age/sex associations within sepsis:\n")


print(
  demographic_association_tests %>%
    select(
      framework,
      clinical_label,
      test,
      n,
      effect,
      p_value
    ) %>%
    tibble::as_tibble(),
  n = Inf
)


cat("\n")


# ==============================================================================
# 11. RECALCULATE GLOBAL FDR INCLUDING AGE / SEX
# ==============================================================================

old_tests <- read.csv(
  script136_tests_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


base_columns <- c(
  "framework",
  "clinical_variable",
  "clinical_label",
  "test_family",
  "test",
  "n",
  "statistic",
  "effect",
  "effect_name",
  "p_value",
  "group_summary",
  "note"
)


missing_base_columns <- setdiff(
  base_columns,
  names(
    old_tests
  )
)


if (length(missing_base_columns) > 0) {
  
  stop(
    paste0(
      "Missing Script 136 columns: ",
      paste(
        missing_base_columns,
        collapse = ", "
      )
    )
  )
}


old_tests_base <- old_tests %>%
  select(
    all_of(
      base_columns
    )
  ) %>%
  filter(
    !clinical_variable %in%
      c(
        "age_years",
        "sex",
        "age_numeric",
        "sex_std"
      )
  )


all_tests_updated <- bind_rows(
  tibble::as_tibble(
    old_tests_base
  ),
  demographic_association_tests
)


all_tests_updated$test_id <- paste(
  all_tests_updated$framework,
  all_tests_updated$clinical_variable,
  all_tests_updated$test,
  sep = "__"
)


# ==============================================================================
# 12. GLOBAL BH FDR
# ==============================================================================

all_tests_updated$BH_global <- NA_real_


valid_idx <- which(
  !is.na(
    all_tests_updated$p_value
  ) &
    is.finite(
      all_tests_updated$p_value
    )
)


if (length(valid_idx) > 0) {
  
  all_tests_updated$BH_global[
    valid_idx
  ] <- stats::p.adjust(
    all_tests_updated$p_value[
      valid_idx
    ],
    method = "BH"
  )
}


# ==============================================================================
# 13. BH WITHIN FRAMEWORK
# ==============================================================================

all_tests_updated$BH_within_framework <- NA_real_


for (
  framework_name in unique(
    all_tests_updated$framework
  )
) {
  
  idx <- which(
    all_tests_updated$framework ==
      framework_name &
      !is.na(
        all_tests_updated$p_value
      ) &
      is.finite(
        all_tests_updated$p_value
      )
  )
  
  if (length(idx) > 0) {
    
    all_tests_updated$BH_within_framework[
      idx
    ] <- stats::p.adjust(
      all_tests_updated$p_value[
        idx
      ],
      method = "BH"
    )
  }
}


# ==============================================================================
# 14. BH WITHIN TEST FAMILY
# ==============================================================================

all_tests_updated$BH_within_test_family <- NA_real_


for (
  family_name in unique(
    all_tests_updated$test_family
  )
) {
  
  idx <- which(
    all_tests_updated$test_family ==
      family_name &
      !is.na(
        all_tests_updated$p_value
      ) &
      is.finite(
        all_tests_updated$p_value
      )
  )
  
  if (length(idx) > 0) {
    
    all_tests_updated$BH_within_test_family[
      idx
    ] <- stats::p.adjust(
      all_tests_updated$p_value[
        idx
      ],
      method = "BH"
    )
  }
}


# ==============================================================================
# 15. SIGNIFICANCE LABELS
# ==============================================================================

all_tests_updated$significance_global <- dplyr::case_when(
  
  !is.na(
    all_tests_updated$BH_global
  ) &
    all_tests_updated$BH_global < 0.05 ~
    "global_FDR_lt_0.05",
  
  !is.na(
    all_tests_updated$p_value
  ) &
    all_tests_updated$p_value < 0.05 ~
    "nominal_p_lt_0.05",
  
  !is.na(
    all_tests_updated$p_value
  ) &
    all_tests_updated$p_value < 0.10 ~
    "trend_p_lt_0.10",
  
  TRUE ~
    "not_significant"
)


all_tests_updated <- all_tests_updated %>%
  arrange(
    BH_global,
    p_value
  )


# ==============================================================================
# 16. SAFE CONSOLE PRINT
# ==============================================================================

cat("Updated top clinical associations:\n")


updated_top_table <- all_tests_updated %>%
  dplyr::select(
    framework,
    clinical_label,
    test,
    n,
    effect,
    p_value,
    BH_global,
    significance_global
  ) %>%
  head(25) %>%
  tibble::as_tibble()


print(
  updated_top_table,
  n = 25
)


cat("\n")


# ==============================================================================
# 17. READ SCRIPT 135 BLOOD SCORES
# ==============================================================================

blood_scores <- read.csv(
  script135_scores_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


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
      "Expected 45 blood samples in Script 135 score file; observed ",
      nrow(blood_scores)
    )
  )
}


if (!"primary_5gene_score" %in% names(blood_scores)) {
  
  stop(
    "primary_5gene_score not found in Script 135 blood score table."
  )
}


# ==============================================================================
# 18. READ SCRIPT 135 SELECTED GENE EXPRESSION
# ==============================================================================

blood_expression <- read.csv(
  script135_expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


blood_expression$sample_id <- clean_sample_id(
  blood_expression$sample_id
)


primary_genes <- c(
  "CD177",
  "HK3",
  "IRAK3",
  "CARD11",
  "IKZF2"
)


missing_primary_genes <- setdiff(
  primary_genes,
  names(
    blood_expression
  )
)


if (length(missing_primary_genes) > 0) {
  
  stop(
    paste0(
      "Missing primary genes in Script 135 expression table: ",
      paste(
        missing_primary_genes,
        collapse = ", "
      )
    )
  )
}


expression_subset <- blood_expression %>%
  select(
    sample_id,
    all_of(
      primary_genes
    )
  )


# ==============================================================================
# 19. BUILD FINAL ALL-BLOOD TABLE FOR ADJUSTED MODELS
# ==============================================================================

blood_all <- blood_scores %>%
  select(
    sample_id,
    primary_5gene_score
  ) %>%
  left_join(
    expression_subset,
    by = "sample_id"
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


if (nrow(blood_all) != 45) {
  
  stop(
    paste0(
      "Expected 45 blood samples after score/expression/demographic merge; observed ",
      nrow(blood_all)
    )
  )
}


# ==============================================================================
# 20. ASSIGN CONDITION DIRECTLY FROM sample_id
# ==============================================================================

blood_all$condition_binary <- dplyr::case_when(
  
  grepl(
    "^BP[0-9]+$",
    blood_all$sample_id
  ) ~ 1L,
  
  grepl(
    "^BC[0-9]+$",
    blood_all$sample_id
  ) ~ 0L,
  
  TRUE ~ NA_integer_
)


blood_all$condition <- dplyr::case_when(
  
  blood_all$condition_binary == 1L ~ "BP",
  
  blood_all$condition_binary == 0L ~ "BC",
  
  TRUE ~ NA_character_
)


blood_all$condition <- factor(
  blood_all$condition,
  levels = c(
    "BC",
    "BP"
  )
)


if (
  any(
    is.na(
      blood_all$condition_binary
    )
  )
) {
  
  stop(
    "Unable to assign BP/BC condition from sample_id."
  )
}


cat("Blood condition distribution for adjusted models:\n")

print(
  table(
    blood_all$condition,
    useNA = "ifany"
  )
)

cat("\n")


cat(
  "Age available in all-blood adjusted dataset: ",
  sum(
    !is.na(
      blood_all$age_years
    )
  ),
  "/45\n",
  sep = ""
)


cat(
  "Sex available in all-blood adjusted dataset: ",
  sum(
    !is.na(
      blood_all$sex
    )
  ),
  "/45\n\n",
  sep = ""
)


# ==============================================================================
# 21. AGE/SEX-ADJUSTED PRIMARY FIVE-GENE SCORE
# ==============================================================================

score_model_data <- blood_all %>%
  filter(
    complete.cases(
      primary_5gene_score,
      condition_binary,
      age_years,
      sex
    )
  )


cat(
  "Samples available for adjusted five-gene model: ",
  nrow(
    score_model_data
  ),
  "\n\n",
  sep = ""
)


if (nrow(score_model_data) < 40) {
  
  warning(
    paste0(
      "Only ",
      nrow(score_model_data),
      " samples available for adjusted score model."
    )
  )
}


score_model <- stats::lm(
  primary_5gene_score ~
    condition_binary +
    age_years +
    sex,
  data = score_model_data
)


score_coeff <- summary(
  score_model
)$coefficients


if (!"condition_binary" %in% rownames(score_coeff)) {
  
  stop(
    "condition_binary coefficient not found in adjusted five-gene model."
  )
}


condition_row <- score_coeff[
  "condition_binary",
  ,
  drop = FALSE
]


score_adjusted_result <- tibble(
  
  model =
    "primary_5gene_score ~ condition + age + sex",
  
  n =
    nrow(
      score_model_data
    ),
  
  condition_effect_BP_vs_BC =
    condition_row[
      1,
      "Estimate"
    ],
  
  standard_error =
    condition_row[
      1,
      "Std. Error"
    ],
  
  t_value =
    condition_row[
      1,
      "t value"
    ],
  
  p_value =
    condition_row[
      1,
      "Pr(>|t|)"
    ],
  
  adjusted_R2 =
    summary(
      score_model
    )$adj.r.squared
)


cat("Age/sex-adjusted primary score model:\n")

print(
  score_adjusted_result,
  n = Inf
)

cat("\n")


# ==============================================================================
# 22. FULL PRIMARY SCORE MODEL COEFFICIENTS
# ==============================================================================

score_full_coefficients <- tibble(
  
  term =
    rownames(
      score_coeff
    ),
  
  estimate =
    score_coeff[
      ,
      "Estimate"
    ],
  
  standard_error =
    score_coeff[
      ,
      "Std. Error"
    ],
  
  t_value =
    score_coeff[
      ,
      "t value"
    ],
  
  p_value =
    score_coeff[
      ,
      "Pr(>|t|)"
    ]
)


cat("Full adjusted primary-score model coefficients:\n")

print(
  score_full_coefficients,
  n = Inf
)

cat("\n")


# ==============================================================================
# 23. AGE/SEX-ADJUSTED INDIVIDUAL PRIMARY GENES
# ==============================================================================

gene_adjusted_results <- vector(
  "list",
  length(
    primary_genes
  )
)


for (i in seq_along(
  primary_genes
)) {
  
  gene <- primary_genes[i]
  
  
  model_data <- tibble(
    
    expression =
      blood_all[[gene]],
    
    condition_binary =
      blood_all$condition_binary,
    
    age_years =
      blood_all$age_years,
    
    sex =
      blood_all$sex
    
  ) %>%
    filter(
      complete.cases(
        expression,
        condition_binary,
        age_years,
        sex
      )
    )
  
  
  fit <- stats::lm(
    expression ~
      condition_binary +
      age_years +
      sex,
    data = model_data
  )
  
  
  coef_table <- summary(
    fit
  )$coefficients
  
  
  if (!"condition_binary" %in% rownames(coef_table)) {
    
    gene_adjusted_results[[i]] <- tibble(
      
      gene =
        gene,
      
      n =
        nrow(
          model_data
        ),
      
      adjusted_condition_logCPM_difference =
        NA_real_,
      
      standard_error =
        NA_real_,
      
      t_value =
        NA_real_,
      
      p_value =
        NA_real_
    )
    
    next
  }
  
  
  condition_coef <- coef_table[
    "condition_binary",
    ,
    drop = FALSE
  ]
  
  
  gene_adjusted_results[[i]] <- tibble(
    
    gene =
      gene,
    
    n =
      nrow(
        model_data
      ),
    
    adjusted_condition_logCPM_difference =
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
      ]
  )
}


gene_adjusted_results <- bind_rows(
  gene_adjusted_results
)


gene_adjusted_results$BH_primary_genes <-
  stats::p.adjust(
    gene_adjusted_results$p_value,
    method = "BH"
  )


gene_adjusted_results <- gene_adjusted_results %>%
  arrange(
    BH_primary_genes,
    p_value
  )


cat("Age/sex-adjusted primary gene effects:\n")

print(
  gene_adjusted_results,
  n = Inf
)

cat("\n")


# ==============================================================================
# 24. CHECK WHETHER SCORE EFFECT REMAINS STRONGLY SIGNIFICANT
# ==============================================================================

score_adjustment_interpretation <- dplyr::case_when(
  
  score_adjusted_result$p_value < 0.001 ~
    "Strong_sepsis_effect_after_age_sex_adjustment",
  
  score_adjusted_result$p_value < 0.05 ~
    "Significant_sepsis_effect_after_age_sex_adjustment",
  
  TRUE ~
    "No_significant_adjusted_sepsis_effect"
)


score_adjusted_result$interpretation <-
  score_adjustment_interpretation


# ==============================================================================
# 25. FIGURE A — AGE DISTRIBUTION BP vs BC
# ==============================================================================

p_age <- ggplot(
  demographics,
  aes(
    x = condition,
    y = age_years
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 2
  ) +
  
  labs(
    title =
      "Age distribution in the blood RNA-seq cohort",
    
    x =
      NULL,
    
    y =
      "Age, years"
  ) +
  
  theme_bw(
    base_size = 12
  )


ggsave(
  filename = file.path(
    figures_dir,
    "136b_Figure_A_age_BP_vs_BC.png"
  ),
  plot = p_age,
  width = 6,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 26. FIGURE B — PRIMARY SCORE vs AGE WITHIN SEPSIS
# ==============================================================================

bp_age_plot <- bp_demo %>%
  filter(
    !is.na(
      age_years
    ),
    !is.na(
      primary_5gene_score
    )
  )


p_score_age <- ggplot(
  bp_age_plot,
  aes(
    x = age_years,
    y = primary_5gene_score
  )
) +
  
  geom_point(
    size = 2.4
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  labs(
    title =
      "Five-gene host-response score versus age in sepsis",
    
    x =
      "Age, years",
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  
  theme_bw(
    base_size = 12
  )


ggsave(
  filename = file.path(
    figures_dir,
    "136b_Figure_B_primary_score_vs_age_BP.png"
  ),
  plot = p_score_age,
  width = 6,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 27. FIGURE C — PRIMARY SCORE BY SEX WITHIN SEPSIS
# ==============================================================================

bp_sex_plot <- bp_demo %>%
  filter(
    !is.na(
      sex
    ),
    !is.na(
      primary_5gene_score
    )
  )


p_score_sex <- ggplot(
  bp_sex_plot,
  aes(
    x = sex,
    y = primary_5gene_score
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 2.2
  ) +
  
  labs(
    title =
      "Five-gene host-response score by sex in sepsis",
    
    x =
      NULL,
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  
  theme_bw(
    base_size = 12
  )


ggsave(
  filename = file.path(
    figures_dir,
    "136b_Figure_C_primary_score_by_sex_BP.png"
  ),
  plot = p_score_sex,
  width = 6,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 28. FIGURE D — PRIMARY SCORE BP vs BC
# ==============================================================================

p_adjusted_score <- ggplot(
  blood_all,
  aes(
    x = condition,
    y = primary_5gene_score
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 2.1
  ) +
  
  labs(
    title =
      "Five-gene host-response score in sepsis and healthy controls",
    
    subtitle =
      "Age and sex are evaluated as covariates in the adjusted model",
    
    x =
      NULL,
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  
  theme_bw(
    base_size = 12
  )


ggsave(
  filename = file.path(
    figures_dir,
    "136b_Figure_D_primary_score_BP_vs_BC.png"
  ),
  plot = p_adjusted_score,
  width = 6,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 29. SAVE TABLES
# ==============================================================================

write.csv(
  demographics,
  file.path(
    tables_dir,
    "136b_blood_demographics.csv"
  ),
  row.names = FALSE
)


write.csv(
  demographic_summary,
  file.path(
    tables_dir,
    "136b_demographic_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  demographic_balance,
  file.path(
    tables_dir,
    "136b_BP_BC_demographic_balance.csv"
  ),
  row.names = FALSE
)


write.csv(
  demographic_association_tests,
  file.path(
    tables_dir,
    "136b_age_sex_endotype_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  all_tests_updated,
  file.path(
    tables_dir,
    "136b_all_clinical_tests_updated_FDR.csv"
  ),
  row.names = FALSE
)


write.csv(
  blood_all,
  file.path(
    tables_dir,
    "136b_all_blood_score_expression_demographics.csv"
  ),
  row.names = FALSE
)


write.csv(
  score_adjusted_result,
  file.path(
    tables_dir,
    "136b_primary_score_age_sex_adjusted.csv"
  ),
  row.names = FALSE
)


write.csv(
  score_full_coefficients,
  file.path(
    tables_dir,
    "136b_primary_score_age_sex_model_coefficients.csv"
  ),
  row.names = FALSE
)


write.csv(
  gene_adjusted_results,
  file.path(
    tables_dir,
    "136b_primary_genes_age_sex_adjusted.csv"
  ),
  row.names = FALSE
)


write.csv(
  bp_demo,
  file.path(
    tables_dir,
    "136b_BP_molecular_clinical_demographics.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 30. PUBLICATION WORKBOOK
# ==============================================================================

wb <- openxlsx::createWorkbook()


sheet_list <- list(
  
  "00_demographic_summary" =
    demographic_summary,
  
  "01_demographic_balance" =
    demographic_balance,
  
  "02_BP_demographics" =
    bp_demo,
  
  "03_age_sex_endotypes" =
    demographic_association_tests,
  
  "04_updated_all_tests" =
    all_tests_updated,
  
  "05_all_blood_adjusted" =
    blood_all,
  
  "06_adjusted_5gene_score" =
    score_adjusted_result,
  
  "07_score_model_coeff" =
    score_full_coefficients,
  
  "08_adjusted_primary_genes" =
    gene_adjusted_results,
  
  "09_all_blood_demographics" =
    demographics
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
    sheet_name,
    sheet_list[[sheet_name]]
  )
}


workbook_file <- file.path(
  tables_dir,
  "136b_demographic_sensitivity.xlsx"
)


openxlsx::saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)


# ==============================================================================
# 31. INPUT MANIFEST
# ==============================================================================

input_manifest <- tibble(
  
  input = c(
    "Supplementary_Table_S1",
    "Script136_tests",
    "Script136_BP_table",
    "Script135_scores",
    "Script135_expression"
  ),
  
  path = c(
    S1_file,
    script136_tests_file,
    script136_bp_file,
    script135_scores_file,
    script135_expression_file
  )
)


file_info <- file.info(
  input_manifest$path
)


input_manifest$file_size_bytes <-
  file_info$size


input_manifest$modified_time <-
  as.character(
    file_info$mtime
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
    "136b_input_file_manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 32. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "136b_sessionInfo.txt"
  )
)


# ==============================================================================
# 33. FINAL COUNTS
# ==============================================================================

n_global_significant <- sum(
  !is.na(
    all_tests_updated$BH_global
  ) &
    all_tests_updated$BH_global < 0.05
)


n_nominal_significant <- sum(
  !is.na(
    all_tests_updated$p_value
  ) &
    all_tests_updated$p_value < 0.05
)


# ==============================================================================
# 34. FINAL REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 136b completed successfully.\n")
cat("====================================================================\n\n")


cat("DEMOGRAPHIC SUMMARY:\n")

print(
  tibble::as_tibble(
    demographic_summary
  ),
  n = Inf
)

cat("\n")


cat("BP vs BC demographic tests:\n")

print(
  tibble::as_tibble(
    demographic_balance
  ),
  n = Inf
)

cat("\n")


cat("Age/sex associations within sepsis:\n")

print(
  demographic_association_tests %>%
    select(
      framework,
      clinical_label,
      test,
      n,
      effect,
      p_value
    ) %>%
    tibble::as_tibble(),
  n = Inf
)

cat("\n")


cat("Blood condition distribution for adjusted models:\n")

print(
  table(
    blood_all$condition,
    useNA = "ifany"
  )
)

cat("\n")


cat("Age/sex-adjusted primary score model:\n")

print(
  tibble::as_tibble(
    score_adjusted_result
  ),
  n = Inf
)

cat("\n")


cat("Full adjusted primary-score model coefficients:\n")

print(
  tibble::as_tibble(
    score_full_coefficients
  ),
  n = Inf
)

cat("\n")


cat("Age/sex-adjusted primary gene effects:\n")

print(
  tibble::as_tibble(
    gene_adjusted_results
  ),
  n = Inf
)

cat("\n")


cat(
  "Updated main clinical tests: ",
  nrow(
    all_tests_updated
  ),
  "\n",
  sep = ""
)


cat(
  "Updated nominal p<0.05: ",
  n_nominal_significant,
  "\n",
  sep = ""
)


cat(
  "Updated global BH FDR<0.05: ",
  n_global_significant,
  "\n\n",
  sep = ""
)


cat("Top updated clinical associations:\n")

print(
  updated_top_table,
  n = 25
)

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
      "136b_primary_score_age_sex_adjusted.csv"
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
      "136b_primary_genes_age_sex_adjusted.csv"
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
      "136b_all_clinical_tests_updated_FDR.csv"
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
      "136b_BP_BC_demographic_balance.csv"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat("IMPORTANT:\n")
cat("- Age and sex are taken from Participant_Metadata in Table S1.\n")
cat("- BP/BC demographic imbalance is explicitly quantified.\n")
cat("- BP/BC class for adjusted models is derived directly from sample_id.\n")
cat("- Primary five-gene score is adjusted for age and sex.\n")
cat("- Each primary gene is independently adjusted for age and sex.\n")
cat("- Age and sex are tested against SRS, SRSq and CTS within sepsis.\n")
cat("- Global clinical FDR is recalculated after adding age/sex tests.\n")
cat("- No new biomarker feature selection was performed.\n")
cat("- No diagnostic/prognostic validation claim should be made.\n")
cat("- Blood only; no urine; no lncRNA.\n\n")


cat("Done.\n")
