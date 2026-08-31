# ==============================================================================
# Script 136
# Blood-only endotype and five-gene signature clinical associations
#
# Project: Sepsis_DESeq2
#
# INPUT:
# 1) Script 135 final BP molecular/endotype table
# 2) Clinical metadata
#
# COHORT:
# BP only, n = 35
#
# MOLECULAR VARIABLES:
# - primary_5gene_score
# - SRSq
# - SRS class
# - CTS class
#
# STATISTICS:
# - Spearman correlations
# - Wilcoxon rank-sum
# - Kruskal-Wallis
# - Fisher exact test
# - global Benjamini-Hochberg FDR
# - framework-specific BH FDR
# - test-family-specific BH FDR
#
# IMPORTANT:
# - exploratory clinical associations only
# - no new feature selection
# - no mortality prediction model
# - no ventilation prediction model
# - no urine
# - no lncRNA
# ==============================================================================


# ==============================================================================
# 0. GENERAL SETTINGS
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

script_id <- "136"
script_name <- "136_blood_endotype_clinical_associations.R"
run_date <- Sys.time()

cat("\n")
cat("====================================================================\n")
cat("Running Script 136\n")
cat("Blood-only endotype and five-gene clinical associations\n")
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
# 1. REQUIRED PACKAGES
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "readxl",
  "openxlsx",
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
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(tibble)
})

cat("Required packages loaded successfully.\n\n")


# ==============================================================================
# 2. BASIC HELPERS
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


normalize_column_name <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  x <- gsub(
    "[^a-z0-9]",
    "",
    x
  )
  
  return(x)
}


detect_exact_column <- function(
    df,
    candidates,
    required = FALSE
) {
  
  original_names <- names(df)
  
  normalized_names <- normalize_column_name(
    original_names
  )
  
  normalized_candidates <- normalize_column_name(
    candidates
  )
  
  for (candidate in normalized_candidates) {
    
    hit <- which(
      normalized_names == candidate
    )
    
    if (length(hit) >= 1) {
      return(
        original_names[hit[1]]
      )
    }
  }
  
  if (required) {
    
    stop(
      paste0(
        "Required column not detected.\n",
        "Candidates: ",
        paste(
          candidates,
          collapse = ", "
        ),
        "\nAvailable columns:\n",
        paste(
          original_names,
          collapse = ", "
        )
      )
    )
  }
  
  return(NA_character_)
}


find_existing_file <- function(paths) {
  
  for (p in paths) {
    
    if (file.exists(p)) {
      return(p)
    }
  }
  
  return(NA_character_)
}


get_column_or_na <- function(
    df,
    column_name
) {
  
  if (is.na(column_name)) {
    return(
      rep(
        NA,
        nrow(df)
      )
    )
  }
  
  if (!column_name %in% names(df)) {
    return(
      rep(
        NA,
        nrow(df)
      )
    )
  }
  
  return(df[[column_name]])
}


parse_numeric_safe <- function(x) {
  
  if (is.numeric(x)) {
    return(
      as.numeric(x)
    )
  }
  
  x <- as.character(x)
  
  result <- vapply(
    x,
    function(value) {
      
      if (
        is.na(value) ||
        trimws(value) == ""
      ) {
        return(NA_real_)
      }
      
      value <- trimws(value)
      
      value <- gsub(
        ",",
        ".",
        value,
        fixed = TRUE
      )
      
      value <- gsub(
        "−",
        "-",
        value,
        fixed = TRUE
      )
      
      hit <- regexpr(
        "-?[0-9]+(?:\\.[0-9]+)?",
        value,
        perl = TRUE
      )
      
      if (hit[1] == -1) {
        return(NA_real_)
      }
      
      number_text <- regmatches(
        value,
        hit
      )
      
      return(
        suppressWarnings(
          as.numeric(
            number_text
          )
        )
      )
    },
    numeric(1)
  )
  
  return(result)
}


# ==============================================================================
# 3. SAFE EXTRACTION FROM htest OBJECTS
# ==============================================================================

safe_htest_statistic <- function(object) {
  
  if (
    is.list(object) &&
    "statistic" %in% names(object) &&
    length(object[["statistic"]]) >= 1
  ) {
    
    return(
      as.numeric(
        object[["statistic"]][1]
      )
    )
  }
  
  return(NA_real_)
}


safe_htest_pvalue <- function(object) {
  
  if (
    is.list(object) &&
    "p.value" %in% names(object) &&
    length(object[["p.value"]]) >= 1
  ) {
    
    return(
      as.numeric(
        object[["p.value"]][1]
      )
    )
  }
  
  return(NA_real_)
}


safe_htest_estimate <- function(object) {
  
  if (
    is.list(object) &&
    "estimate" %in% names(object) &&
    length(object[["estimate"]]) >= 1
  ) {
    
    return(
      as.numeric(
        object[["estimate"]][1]
      )
    )
  }
  
  return(NA_real_)
}


# ==============================================================================
# 4. CLINICAL NORMALIZATION HELPERS
# ==============================================================================

normalize_outcome_text <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  result <- rep(
    NA_character_,
    length(x)
  )
  
  died_flag <- grepl(
    "died|dead|death|deceased|умер|летал|сконч",
    x,
    perl = TRUE
  )
  
  survived_flag <- grepl(
    "discharg|surviv|выпис|жив",
    x,
    perl = TRUE
  )
  
  result[died_flag] <- "Died"
  
  result[
    survived_flag &
      is.na(result)
  ] <- "Discharged"
  
  return(result)
}


normalize_mortality_numeric <- function(x) {
  
  x_numeric <- suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
  
  result <- rep(
    NA_character_,
    length(x_numeric)
  )
  
  result[
    !is.na(x_numeric) &
      x_numeric == 1
  ] <- "Died"
  
  result[
    !is.na(x_numeric) &
      x_numeric == 0
  ] <- "Discharged"
  
  return(result)
}


normalize_yes_no <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  result <- rep(
    NA_character_,
    length(x)
  )
  
  negative <- (
    x %in%
      c(
        "0",
        "no",
        "n",
        "false",
        "нет"
      )
  ) |
    grepl(
      "without|not ventilated|не провод|не было",
      x,
      perl = TRUE
    )
  
  positive <- (
    x %in%
      c(
        "1",
        "yes",
        "y",
        "true",
        "да"
      )
  ) |
    grepl(
      "mechanical ventilation|ventilat|ivl|ивл",
      x,
      perl = TRUE
    )
  
  result[negative] <- "No"
  
  result[
    positive &
      is.na(result)
  ] <- "Yes"
  
  return(result)
}


normalize_culture <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  result <- rep(
    NA_character_,
    length(x)
  )
  
  negative <- (
    x %in%
      c(
        "negative",
        "neg",
        "0",
        "no",
        "нет"
      )
  ) |
    grepl(
      "negative|отриц",
      x,
      perl = TRUE
    )
  
  positive <- (
    x %in%
      c(
        "positive",
        "pos",
        "1",
        "yes",
        "да"
      )
  ) |
    grepl(
      "positive|полож",
      x,
      perl = TRUE
    )
  
  result[negative] <- "Negative"
  
  result[
    positive &
      is.na(result)
  ] <- "Positive"
  
  return(result)
}


normalize_sex <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  result <- rep(
    NA_character_,
    length(x)
  )
  
  male <- (
    x %in%
      c(
        "m",
        "male",
        "1",
        "м",
        "муж"
      )
  ) |
    grepl(
      "male|муж",
      x,
      perl = TRUE
    )
  
  female <- (
    x %in%
      c(
        "f",
        "female",
        "0",
        "ж",
        "жен"
      )
  ) |
    grepl(
      "female|жен",
      x,
      perl = TRUE
    )
  
  result[male] <- "Male"
  
  result[
    female &
      is.na(result)
  ] <- "Female"
  
  return(result)
}


# ==============================================================================
# 5. STATISTICAL HELPERS
# ==============================================================================

epsilon_squared_kw <- function(
    x,
    group
) {
  
  keep <- complete.cases(
    x,
    group
  )
  
  keep <- keep &
    is.finite(x)
  
  x2 <- x[keep]
  
  group2 <- droplevels(
    factor(
      group[keep]
    )
  )
  
  if (
    length(x2) < 5 ||
    nlevels(group2) < 2
  ) {
    return(NA_real_)
  }
  
  ht <- stats::kruskal.test(
    x2 ~ group2
  )
  
  H <- safe_htest_statistic(
    ht
  )
  
  if (!is.finite(H)) {
    return(NA_real_)
  }
  
  n <- length(x2)
  k <- nlevels(group2)
  
  epsilon2 <- (
    H - k + 1
  ) / (
    n - k
  )
  
  return(
    max(
      0,
      epsilon2
    )
  )
}


make_group_summary <- function(
    x,
    group
) {
  
  tmp <- tibble(
    x = x,
    group = as.character(group)
  ) %>%
    filter(
      !is.na(x),
      is.finite(x),
      !is.na(group),
      group != ""
    ) %>%
    group_by(
      group
    ) %>%
    summarise(
      n = n(),
      
      median = median(
        x,
        na.rm = TRUE
      ),
      
      IQR = IQR(
        x,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  if (nrow(tmp) == 0) {
    return("")
  }
  
  pieces <- paste0(
    tmp$group,
    ": n=",
    tmp$n,
    ", median=",
    signif(
      tmp$median,
      4
    ),
    ", IQR=",
    signif(
      tmp$IQR,
      4
    )
  )
  
  return(
    paste(
      pieces,
      collapse = " | "
    )
  )
}


safe_spearman <- function(
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
  
  keep <- keep &
    is.finite(x) &
    is.finite(y)
  
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
        note = "Insufficient valid observations"
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
      note = ""
    )
  )
}


safe_group_continuous_test <- function(
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
  
  keep <- keep &
    is.finite(x)
  
  x2 <- x[keep]
  
  group2 <- droplevels(
    factor(
      group[keep]
    )
  )
  
  n_total <- length(x2)
  
  if (
    n_total < 6 ||
    nlevels(group2) < 2
  ) {
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = test_family,
        test = NA_character_,
        n = n_total,
        statistic = NA_real_,
        effect = NA_real_,
        effect_name = NA_character_,
        p_value = NA_real_,
        group_summary = make_group_summary(
          x2,
          group2
        ),
        note = "Insufficient groups"
      )
    )
  }
  
  group_counts <- table(
    group2
  )
  
  if (any(group_counts < 2)) {
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = test_family,
        test = NA_character_,
        n = n_total,
        statistic = NA_real_,
        effect = NA_real_,
        effect_name = NA_character_,
        p_value = NA_real_,
        group_summary = make_group_summary(
          x2,
          group2
        ),
        note = "At least one group has fewer than 2 observations"
      )
    )
  }
  
  # --------------------------------------------------------------------------
  # TWO GROUPS — explicit vector-based Wilcoxon test
  # --------------------------------------------------------------------------
  
  if (nlevels(group2) == 2) {
    
    group_levels <- levels(
      group2
    )
    
    group1_name <- group_levels[1]
    group2_name <- group_levels[2]
    
    x_group1 <- x2[
      group2 == group1_name
    ]
    
    x_group2 <- x2[
      group2 == group2_name
    ]
    
    ht <- stats::wilcox.test(
      x = x_group1,
      y = x_group2,
      exact = FALSE,
      paired = FALSE,
      conf.int = FALSE
    )
    
    median_group1 <- median(
      x_group1,
      na.rm = TRUE
    )
    
    median_group2 <- median(
      x_group2,
      na.rm = TRUE
    )
    
    effect_value <- (
      median_group2 -
        median_group1
    )
    
    effect_name_value <- paste0(
      "median_difference_",
      group2_name,
      "_minus_",
      group1_name
    )
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = test_family,
        test = "Wilcoxon_rank_sum",
        n = n_total,
        statistic = safe_htest_statistic(ht),
        effect = effect_value,
        effect_name = effect_name_value,
        p_value = safe_htest_pvalue(ht),
        group_summary = make_group_summary(
          x2,
          group2
        ),
        note = ""
      )
    )
  }
  
  # --------------------------------------------------------------------------
  # THREE OR MORE GROUPS — Kruskal-Wallis
  # --------------------------------------------------------------------------
  
  ht <- stats::kruskal.test(
    x2 ~ group2
  )
  
  epsilon2 <- epsilon_squared_kw(
    x = x2,
    group = group2
  )
  
  return(
    tibble(
      framework = framework,
      clinical_variable = clinical_variable,
      clinical_label = clinical_label,
      test_family = test_family,
      test = "Kruskal_Wallis",
      n = n_total,
      statistic = safe_htest_statistic(ht),
      effect = epsilon2,
      effect_name = "epsilon_squared",
      p_value = safe_htest_pvalue(ht),
      group_summary = make_group_summary(
        x2,
        group2
      ),
      note = ""
    )
  )
}


safe_fisher_test <- function(
    molecular_group,
    clinical_group,
    framework,
    clinical_variable,
    clinical_label
) {
  
  keep <- !is.na(
    molecular_group
  ) &
    !is.na(
      clinical_group
    )
  
  molecular_group2 <- droplevels(
    factor(
      molecular_group[keep]
    )
  )
  
  clinical_group2 <- droplevels(
    factor(
      clinical_group[keep]
    )
  )
  
  n_total <- length(
    molecular_group2
  )
  
  if (
    n_total < 6 ||
    nlevels(molecular_group2) < 2 ||
    nlevels(clinical_group2) < 2
  ) {
    
    return(
      tibble(
        framework = framework,
        clinical_variable = clinical_variable,
        clinical_label = clinical_label,
        test_family = "categorical_endotype",
        test = "Fisher_exact",
        n = n_total,
        statistic = NA_real_,
        effect = NA_real_,
        effect_name = NA_character_,
        p_value = NA_real_,
        group_summary = "",
        note = "Insufficient groups"
      )
    )
  }
  
  contingency <- table(
    molecular_group2,
    clinical_group2
  )
  
  ht <- stats::fisher.test(
    contingency
  )
  
  summary_text <- paste(
    capture.output(
      print(
        contingency
      )
    ),
    collapse = " "
  )
  
  return(
    tibble(
      framework = framework,
      clinical_variable = clinical_variable,
      clinical_label = clinical_label,
      test_family = "categorical_endotype",
      test = "Fisher_exact",
      n = n_total,
      statistic = NA_real_,
      effect = NA_real_,
      effect_name = NA_character_,
      p_value = safe_htest_pvalue(ht),
      group_summary = summary_text,
      note = ""
    )
  )
}


# ==============================================================================
# 6. INPUT FILES
# ==============================================================================

script135_file <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "135_validation",
  "tables",
  "135_sepsis_blood_scores_with_SRS_CTS.csv"
)


clinical_file <- find_existing_file(
  c(
    file.path(
      "data",
      "metadata_clinical_annotated_biomarkers.xlsx"
    ),
    
    file.path(
      "data",
      "metadata_clinical_annotated.csv"
    ),
    
    file.path(
      "data",
      "metadata_clinical.xlsx"
    )
  )
)


input_check <- tibble(
  
  input = c(
    "Script135_BP_scores",
    "clinical_metadata"
  ),
  
  path = c(
    script135_file,
    clinical_file
  ),
  
  exists = c(
    file.exists(
      script135_file
    ),
    
    !is.na(clinical_file) &&
      file.exists(
        clinical_file
      )
  )
)


cat("Input file check:\n")
print(input_check)
cat("\n")


if (any(!input_check$exists)) {
  
  stop(
    "Required input file missing."
  )
}


cat(
  "All required input files identified successfully.\n\n"
)


# ==============================================================================
# 7. OUTPUT DIRECTORIES
# ==============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "136_clinical_associations"
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
# 8. READ SCRIPT 135 RESULTS
# ==============================================================================

cat("Reading Script 135 output:\n")
cat(
  script135_file,
  "\n\n"
)


molecular <- read.csv(
  script135_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


molecular$sample_id <- clean_sample_id(
  molecular$sample_id
)


molecular <- molecular %>%
  filter(
    grepl(
      "^BP[0-9]+$",
      sample_id
    )
  )


if (nrow(molecular) != 35) {
  
  stop(
    paste0(
      "Expected 35 BP samples from Script 135; observed ",
      nrow(molecular)
    )
  )
}


required_135_columns <- c(
  "sample_id",
  "SRS",
  "SRSq",
  "CTS",
  "primary_5gene_score",
  "myeloid_UP_score",
  "adaptive_suppression_score"
)


missing_135_columns <- setdiff(
  required_135_columns,
  names(molecular)
)


if (length(missing_135_columns) > 0) {
  
  stop(
    paste0(
      "Missing Script 135 columns: ",
      paste(
        missing_135_columns,
        collapse = ", "
      )
    )
  )
}


if (!"mNN_outlier_final" %in% names(molecular)) {
  
  molecular$mNN_outlier_final <-
    molecular$sample_id %in%
    c(
      "BP27",
      "BP31",
      "BP26",
      "BP10"
    )
}


if (!"CTS_SRS_group" %in% names(molecular)) {
  
  molecular$CTS_SRS_group <- paste(
    molecular$CTS,
    molecular$SRS,
    sep = "/"
  )
}


cat("Script 135 SRS distribution:\n")

print(
  table(
    molecular$SRS
  )
)


cat("\nScript 135 CTS distribution:\n")

print(
  table(
    molecular$CTS
  )
)

cat("\n")


observed_srs <- table(
  factor(
    molecular$SRS,
    levels = c(
      "SRS1",
      "SRS2"
    )
  )
)


if (!identical(
  as.integer(observed_srs),
  c(
    28L,
    7L
  )
)) {
  
  stop(
    "Script 135 SRS structure does not match expected 28/7."
  )
}


observed_cts <- table(
  factor(
    molecular$CTS,
    levels = c(
      "CTS1",
      "CTS2",
      "CTS3"
    )
  )
)


if (!identical(
  as.integer(observed_cts),
  c(
    14L,
    6L,
    15L
  )
)) {
  
  stop(
    "Script 135 CTS structure does not match expected 14/6/15."
  )
}


cat(
  "Script 135 molecular input validation PASSED.\n\n"
)


# ==============================================================================
# 9. READ CLINICAL METADATA
# ==============================================================================

cat("Reading clinical metadata:\n")
cat(
  clinical_file,
  "\n\n"
)


clinical_extension <- tolower(
  tools::file_ext(
    clinical_file
  )
)


selected_clinical_sheet <- "CSV"


if (
  clinical_extension %in%
  c(
    "xlsx",
    "xls"
  )
) {
  
  clinical_sheets <- readxl::excel_sheets(
    clinical_file
  )
  
  
  candidate_sheet_data <- vector(
    "list",
    length(
      clinical_sheets
    )
  )
  
  
  candidate_sheet_scores <- numeric(
    length(
      clinical_sheets
    )
  )
  
  
  for (i in seq_along(
    clinical_sheets
  )) {
    
    tmp <- suppressMessages(
      readxl::read_excel(
        clinical_file,
        sheet = clinical_sheets[i]
      )
    )
    
    
    candidate_sheet_data[[i]] <- tmp
    
    
    sample_col_test <- detect_exact_column(
      tmp,
      c(
        "sample_id",
        "Sample_ID",
        "sample"
      ),
      required = FALSE
    )
    
    
    score <- 0
    
    
    if (!is.na(sample_col_test)) {
      score <- score + 10
    }
    
    
    useful_candidates <- c(
      "clinical_outcome",
      "mortality_28d",
      "ivl_status",
      "creatinine_value",
      "crp_mg_l",
      "procalcitonin_ng_ml",
      "lactate_mmol_l",
      "wbc_10e9_l",
      "platelets_10e9_l"
    )
    
    
    for (candidate in useful_candidates) {
      
      candidate_hit <- detect_exact_column(
        tmp,
        candidate,
        required = FALSE
      )
      
      
      if (!is.na(candidate_hit)) {
        score <- score + 1
      }
    }
    
    
    candidate_sheet_scores[i] <- score
  }
  
  
  best_sheet_index <- which.max(
    candidate_sheet_scores
  )
  
  
  clinical_raw <- candidate_sheet_data[[best_sheet_index]]
  
  
  selected_clinical_sheet <-
    clinical_sheets[best_sheet_index]
  
  
} else {
  
  clinical_raw <- read.csv(
    clinical_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


cat("Selected clinical worksheet/source:\n")
cat(
  selected_clinical_sheet,
  "\n\n"
)


cat(
  "Clinical metadata dimensions: ",
  nrow(clinical_raw),
  " rows x ",
  ncol(clinical_raw),
  " columns\n\n",
  sep = ""
)


# ==============================================================================
# 10. STRICT CLINICAL COLUMN MAPPING
# ==============================================================================

column_mapping <- list(
  
  sample_id = detect_exact_column(
    clinical_raw,
    c(
      "sample_id",
      "Sample_ID",
      "sample"
    ),
    required = TRUE
  ),
  
  
  patient_id = detect_exact_column(
    clinical_raw,
    c(
      "patient_id",
      "Patient_ID",
      "patient"
    ),
    required = FALSE
  ),
  
  
  outcome = detect_exact_column(
    clinical_raw,
    c(
      "clinical_outcome",
      "outcome",
      "outcome_status"
    ),
    required = FALSE
  ),
  
  
  mortality = detect_exact_column(
    clinical_raw,
    c(
      "mortality_28d",
      "mortality",
      "mortality28d",
      "death_28d"
    ),
    required = FALSE
  ),
  
  
  ventilation = detect_exact_column(
    clinical_raw,
    c(
      "ivl_status",
      "IVL_status",
      "ivl",
      "mechanical_ventilation",
      "ventilation_status",
      "ventilation"
    ),
    required = FALSE
  ),
  
  
  culture = detect_exact_column(
    clinical_raw,
    c(
      "culture_status",
      "Culture_status",
      "culture"
    ),
    required = FALSE
  ),
  
  
  creatinine = detect_exact_column(
    clinical_raw,
    c(
      "creatinine_value",
      "creatinine",
      "serum_creatinine",
      "creatinine_umol_l"
    ),
    required = FALSE
  ),
  
  
  crp = detect_exact_column(
    clinical_raw,
    c(
      "crp_mg_l",
      "CRP",
      "crp",
      "crp_value"
    ),
    required = FALSE
  ),
  
  
  pct = detect_exact_column(
    clinical_raw,
    c(
      "procalcitonin_ng_ml",
      "procalcitonin",
      "PCT",
      "pct"
    ),
    required = FALSE
  ),
  
  
  lactate = detect_exact_column(
    clinical_raw,
    c(
      "lactate_mmol_l",
      "lactate",
      "Lactate"
    ),
    required = FALSE
  ),
  
  
  wbc = detect_exact_column(
    clinical_raw,
    c(
      "wbc_10e9_l",
      "WBC",
      "wbc"
    ),
    required = FALSE
  ),
  
  
  lymphocytes = detect_exact_column(
    clinical_raw,
    c(
      "lymphocytes_10e9_l",
      "lymphocyte_10e9_l",
      "lymphocytes",
      "lymphocyte",
      "LYM",
      "lymphocyte_count"
    ),
    required = FALSE
  ),
  
  
  platelets = detect_exact_column(
    clinical_raw,
    c(
      "platelets_10e9_l",
      "platelets",
      "platelet",
      "PLT"
    ),
    required = FALSE
  ),
  
  
  albumin = detect_exact_column(
    clinical_raw,
    c(
      "albumin_g_l",
      "albumin",
      "Albumin"
    ),
    required = FALSE
  ),
  
  
  alt = detect_exact_column(
    clinical_raw,
    c(
      "alt_u_l",
      "ALT",
      "alt"
    ),
    required = FALSE
  ),
  
  
  ast = detect_exact_column(
    clinical_raw,
    c(
      "ast_u_l",
      "AST",
      "ast"
    ),
    required = FALSE
  ),
  
  
  sofa = detect_exact_column(
    clinical_raw,
    c(
      "sofa_score",
      "SOFA_score",
      "SOFA",
      "sofa"
    ),
    required = FALSE
  ),
  
  
  age = detect_exact_column(
    clinical_raw,
    c(
      "age_years",
      "age",
      "Age"
    ),
    required = FALSE
  ),
  
  
  sex = detect_exact_column(
    clinical_raw,
    c(
      "sex",
      "Sex",
      "gender"
    ),
    required = FALSE
  )
)


mapping_table <- tibble(
  
  standardized_variable =
    names(
      column_mapping
    ),
  
  detected_column =
    unname(
      unlist(
        column_mapping
      )
    )
)


cat("Clinical column mapping:\n")

print(
  mapping_table,
  n = Inf
)

cat("\n")


# ==============================================================================
# 11. STANDARDIZE CLINICAL DATA
# ==============================================================================

clinical_std <- tibble(
  
  sample_id = clean_sample_id(
    clinical_raw[[column_mapping$sample_id]]
  ),
  
  
  patient_id = as.character(
    get_column_or_na(
      clinical_raw,
      column_mapping$patient_id
    )
  ),
  
  
  outcome_raw = as.character(
    get_column_or_na(
      clinical_raw,
      column_mapping$outcome
    )
  ),
  
  
  mortality_raw = as.character(
    get_column_or_na(
      clinical_raw,
      column_mapping$mortality
    )
  ),
  
  
  ventilation_raw = as.character(
    get_column_or_na(
      clinical_raw,
      column_mapping$ventilation
    )
  ),
  
  
  culture_raw = as.character(
    get_column_or_na(
      clinical_raw,
      column_mapping$culture
    )
  ),
  
  
  creatinine_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$creatinine
    )
  ),
  
  
  crp_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$crp
    )
  ),
  
  
  pct_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$pct
    )
  ),
  
  
  lactate_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$lactate
    )
  ),
  
  
  wbc_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$wbc
    )
  ),
  
  
  lymphocytes_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$lymphocytes
    )
  ),
  
  
  platelets_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$platelets
    )
  ),
  
  
  albumin_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$albumin
    )
  ),
  
  
  alt_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$alt
    )
  ),
  
  
  ast_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$ast
    )
  ),
  
  
  sofa_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$sofa
    )
  ),
  
  
  age_numeric = parse_numeric_safe(
    get_column_or_na(
      clinical_raw,
      column_mapping$age
    )
  ),
  
  
  sex_raw = as.character(
    get_column_or_na(
      clinical_raw,
      column_mapping$sex
    )
  )
)


clinical_std$outcome_status <-
  normalize_outcome_text(
    clinical_std$outcome_raw
  )


mortality_fallback <-
  normalize_mortality_numeric(
    clinical_std$mortality_raw
  )


missing_outcome <- is.na(
  clinical_std$outcome_status
)


clinical_std$outcome_status[
  missing_outcome
] <- mortality_fallback[
  missing_outcome
]


clinical_std$ventilation_status <-
  normalize_yes_no(
    clinical_std$ventilation_raw
  )


clinical_std$culture_status <-
  normalize_culture(
    clinical_std$culture_raw
  )


clinical_std$sex_std <-
  normalize_sex(
    clinical_std$sex_raw
  )


clinical_std$creatinine_group_110 <-
  ifelse(
    is.na(
      clinical_std$creatinine_numeric
    ),
    NA_character_,
    ifelse(
      clinical_std$creatinine_numeric > 110,
      "High_gt110",
      "Low_le110"
    )
  )


clinical_std <- clinical_std %>%
  filter(
    grepl(
      "^BP[0-9]+$",
      sample_id
    )
  )


if (
  anyDuplicated(
    clinical_std$sample_id
  ) > 0
) {
  
  duplicated_ids <- unique(
    clinical_std$sample_id[
      duplicated(
        clinical_std$sample_id
      )
    ]
  )
  
  stop(
    paste0(
      "Duplicated BP sample IDs in clinical metadata: ",
      paste(
        duplicated_ids,
        collapse = ", "
      )
    )
  )
}


cat(
  "BP clinical rows detected: ",
  nrow(clinical_std),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 12. MERGE MOLECULAR + CLINICAL DATA
# ==============================================================================

bp_merged <- molecular %>%
  left_join(
    clinical_std,
    by = "sample_id"
  )


if (nrow(bp_merged) != 35) {
  
  stop(
    "Merged BP table does not contain exactly 35 samples."
  )
}


cat(
  "Merged molecular-clinical table: ",
  nrow(bp_merged),
  " BP samples.\n\n",
  sep = ""
)


# ==============================================================================
# 13. CLINICAL COHORT CHECKS
# ==============================================================================

cat("Clinical cohort checks:\n\n")


cat("Outcome:\n")

print(
  table(
    bp_merged$outcome_status,
    useNA = "ifany"
  )
)


cat("\nMechanical ventilation:\n")

print(
  table(
    bp_merged$ventilation_status,
    useNA = "ifany"
  )
)


cat("\nCulture:\n")

print(
  table(
    bp_merged$culture_status,
    useNA = "ifany"
  )
)


cat("\nCreatinine >110 umol/L:\n")

print(
  table(
    bp_merged$creatinine_group_110,
    useNA = "ifany"
  )
)


cat("\n")


# ==============================================================================
# 14. CLINICAL VARIABLE DEFINITIONS
# ==============================================================================

continuous_variables <- tibble(
  
  variable = c(
    "crp_numeric",
    "pct_numeric",
    "lactate_numeric",
    "creatinine_numeric",
    "wbc_numeric",
    "lymphocytes_numeric",
    "platelets_numeric",
    "albumin_numeric",
    "alt_numeric",
    "ast_numeric",
    "sofa_numeric",
    "age_numeric"
  ),
  
  
  label = c(
    "CRP",
    "Procalcitonin",
    "Lactate",
    "Creatinine",
    "WBC",
    "Lymphocytes",
    "Platelets",
    "Albumin",
    "ALT",
    "AST",
    "SOFA",
    "Age"
  ),
  
  
  priority = c(
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "secondary",
    "secondary",
    "secondary",
    "secondary",
    "secondary"
  )
)


categorical_variables <- tibble(
  
  variable = c(
    "outcome_status",
    "ventilation_status",
    "culture_status",
    "creatinine_group_110",
    "sex_std"
  ),
  
  
  label = c(
    "Outcome",
    "Mechanical ventilation",
    "Culture status",
    "Creatinine >110",
    "Sex"
  ),
  
  
  priority = c(
    "primary",
    "primary",
    "secondary",
    "secondary",
    "secondary"
  )
)


# ==============================================================================
# 15. VARIABLE AVAILABILITY
# ==============================================================================

continuous_availability_list <- vector(
  "list",
  nrow(
    continuous_variables
  )
)


for (i in seq_len(
  nrow(
    continuous_variables
  )
)) {
  
  variable_name <-
    continuous_variables$variable[i]
  
  
  values <- bp_merged[[variable_name]]
  
  
  valid_values <- values[
    !is.na(values) &
      is.finite(values)
  ]
  
  
  continuous_availability_list[[i]] <-
    tibble(
      
      variable =
        variable_name,
      
      label =
        continuous_variables$label[i],
      
      priority =
        continuous_variables$priority[i],
      
      n_available =
        length(
          valid_values
        ),
      
      n_unique =
        length(
          unique(
            valid_values
          )
        ),
      
      available_for_analysis =
        length(valid_values) >= 8 &&
        length(
          unique(
            valid_values
          )
        ) >= 3
    )
}


availability_continuous <- bind_rows(
  continuous_availability_list
)


categorical_availability_list <- vector(
  "list",
  nrow(
    categorical_variables
  )
)


for (i in seq_len(
  nrow(
    categorical_variables
  )
)) {
  
  variable_name <-
    categorical_variables$variable[i]
  
  
  values <- bp_merged[[variable_name]]
  
  
  valid_values <- values[
    !is.na(values)
  ]
  
  
  categorical_availability_list[[i]] <-
    tibble(
      
      variable =
        variable_name,
      
      label =
        categorical_variables$label[i],
      
      priority =
        categorical_variables$priority[i],
      
      n_available =
        length(
          valid_values
        ),
      
      n_levels =
        length(
          unique(
            valid_values
          )
        ),
      
      available_for_analysis =
        length(valid_values) >= 6 &&
        length(
          unique(
            valid_values
          )
        ) >= 2
    )
}


availability_categorical <- bind_rows(
  categorical_availability_list
)


cat("Continuous clinical availability:\n")

print(
  availability_continuous,
  n = Inf
)


cat("\nCategorical clinical availability:\n")

print(
  availability_categorical,
  n = Inf
)


cat("\n")


# ==============================================================================
# 16. CONTINUOUS CLINICAL ASSOCIATIONS
# ==============================================================================

continuous_results <- list()

result_counter <- 1L


for (i in seq_len(
  nrow(
    availability_continuous
  )
)) {
  
  if (
    !availability_continuous$
    available_for_analysis[i]
  ) {
    next
  }
  
  
  variable_name <-
    availability_continuous$variable[i]
  
  
  variable_label <-
    availability_continuous$label[i]
  
  
  clinical_values <-
    bp_merged[[variable_name]]
  
  
  cat(
    "Testing continuous variable: ",
    variable_label,
    "\n",
    sep = ""
  )
  
  
  continuous_results[[result_counter]] <-
    safe_spearman(
      x =
        bp_merged$primary_5gene_score,
      
      y =
        clinical_values,
      
      framework =
        "Primary_5gene_score",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label
    )
  
  
  result_counter <-
    result_counter + 1L
  
  
  continuous_results[[result_counter]] <-
    safe_spearman(
      x =
        bp_merged$SRSq,
      
      y =
        clinical_values,
      
      framework =
        "SRSq",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label
    )
  
  
  result_counter <-
    result_counter + 1L
  
  
  continuous_results[[result_counter]] <-
    safe_group_continuous_test(
      x =
        clinical_values,
      
      group =
        bp_merged$SRS,
      
      framework =
        "SRS_class",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label,
      
      test_family =
        "continuous_by_endotype"
    )
  
  
  result_counter <-
    result_counter + 1L
  
  
  continuous_results[[result_counter]] <-
    safe_group_continuous_test(
      x =
        clinical_values,
      
      group =
        bp_merged$CTS,
      
      framework =
        "CTS_class",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label,
      
      test_family =
        "continuous_by_endotype"
    )
  
  
  result_counter <-
    result_counter + 1L
}


continuous_results_table <- bind_rows(
  continuous_results
)


cat(
  "\nContinuous clinical association tests completed: ",
  nrow(
    continuous_results_table
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 17. CATEGORICAL CLINICAL ASSOCIATIONS
# ==============================================================================

categorical_results <- list()

result_counter <- 1L


for (i in seq_len(
  nrow(
    availability_categorical
  )
)) {
  
  if (
    !availability_categorical$
    available_for_analysis[i]
  ) {
    next
  }
  
  
  variable_name <-
    availability_categorical$variable[i]
  
  
  variable_label <-
    availability_categorical$label[i]
  
  
  clinical_group <-
    bp_merged[[variable_name]]
  
  
  cat(
    "Testing categorical variable: ",
    variable_label,
    "\n",
    sep = ""
  )
  
  
  categorical_results[[result_counter]] <-
    safe_fisher_test(
      molecular_group =
        bp_merged$SRS,
      
      clinical_group =
        clinical_group,
      
      framework =
        "SRS_class",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label
    )
  
  
  result_counter <-
    result_counter + 1L
  
  
  categorical_results[[result_counter]] <-
    safe_fisher_test(
      molecular_group =
        bp_merged$CTS,
      
      clinical_group =
        clinical_group,
      
      framework =
        "CTS_class",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label
    )
  
  
  result_counter <-
    result_counter + 1L
  
  
  categorical_results[[result_counter]] <-
    safe_group_continuous_test(
      x =
        bp_merged$primary_5gene_score,
      
      group =
        clinical_group,
      
      framework =
        "Primary_5gene_score",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label,
      
      test_family =
        "molecular_score_by_clinical_group"
    )
  
  
  result_counter <-
    result_counter + 1L
  
  
  categorical_results[[result_counter]] <-
    safe_group_continuous_test(
      x =
        bp_merged$SRSq,
      
      group =
        clinical_group,
      
      framework =
        "SRSq",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label,
      
      test_family =
        "molecular_score_by_clinical_group"
    )
  
  
  result_counter <-
    result_counter + 1L
}


categorical_results_table <- bind_rows(
  categorical_results
)


cat(
  "\nCategorical clinical association tests completed: ",
  nrow(
    categorical_results_table
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 18. COMBINE MAIN TESTS
# ==============================================================================

all_tests <- bind_rows(
  continuous_results_table,
  categorical_results_table
)


all_tests$test_id <- paste(
  all_tests$framework,
  all_tests$clinical_variable,
  all_tests$test,
  sep = "__"
)


# ==============================================================================
# 19. GLOBAL BH FDR
# ==============================================================================

all_tests$BH_global <- NA_real_


valid_p <- which(
  !is.na(
    all_tests$p_value
  ) &
    is.finite(
      all_tests$p_value
    )
)


if (length(valid_p) > 0) {
  
  all_tests$BH_global[
    valid_p
  ] <- stats::p.adjust(
    all_tests$p_value[
      valid_p
    ],
    method = "BH"
  )
}


# ==============================================================================
# 20. BH WITHIN FRAMEWORK
# ==============================================================================

all_tests$BH_within_framework <- NA_real_


for (
  framework_name in unique(
    all_tests$framework
  )
) {
  
  idx <- which(
    all_tests$framework ==
      framework_name &
      !is.na(
        all_tests$p_value
      ) &
      is.finite(
        all_tests$p_value
      )
  )
  
  
  if (length(idx) > 0) {
    
    all_tests$BH_within_framework[
      idx
    ] <- stats::p.adjust(
      all_tests$p_value[
        idx
      ],
      method = "BH"
    )
  }
}


# ==============================================================================
# 21. BH WITHIN TEST FAMILY
# ==============================================================================

all_tests$BH_within_test_family <- NA_real_


for (
  family_name in unique(
    all_tests$test_family
  )
) {
  
  idx <- which(
    all_tests$test_family ==
      family_name &
      !is.na(
        all_tests$p_value
      ) &
      is.finite(
        all_tests$p_value
      )
  )
  
  
  if (length(idx) > 0) {
    
    all_tests$BH_within_test_family[
      idx
    ] <- stats::p.adjust(
      all_tests$p_value[
        idx
      ],
      method = "BH"
    )
  }
}


# ==============================================================================
# 22. SIGNIFICANCE LABELS
# ==============================================================================

all_tests$significance_global <- dplyr::case_when(
  
  !is.na(
    all_tests$BH_global
  ) &
    all_tests$BH_global < 0.05 ~
    "global_FDR_lt_0.05",
  
  !is.na(
    all_tests$p_value
  ) &
    all_tests$p_value < 0.05 ~
    "nominal_p_lt_0.05",
  
  !is.na(
    all_tests$p_value
  ) &
    all_tests$p_value < 0.10 ~
    "trend_p_lt_0.10",
  
  TRUE ~
    "not_significant"
)


all_tests <- all_tests %>%
  arrange(
    BH_global,
    p_value
  )


cat("Top clinical associations:\n")


print(
  all_tests %>%
    select(
      framework,
      clinical_label,
      test,
      n,
      effect,
      effect_name,
      p_value,
      BH_global,
      significance_global
    ) %>%
    head(25),
  n = 25
)


cat("\n")


# ==============================================================================
# 23. PRIMARY-PRIORITY ASSOCIATIONS
# ==============================================================================

primary_variable_names <- c(
  
  continuous_variables$variable[
    continuous_variables$priority ==
      "primary"
  ],
  
  categorical_variables$variable[
    categorical_variables$priority ==
      "primary"
  ]
)


primary_association_table <- all_tests %>%
  filter(
    clinical_variable %in%
      primary_variable_names
  ) %>%
  arrange(
    BH_global,
    p_value
  )


# ==============================================================================
# 24. SRSq SENSITIVITY WITHOUT mNN OUTLIERS
# ==============================================================================

bp_no_mnn <- bp_merged %>%
  filter(
    !mNN_outlier_final
  )


if (nrow(bp_no_mnn) != 31) {
  
  warning(
    paste0(
      "Expected n=31 after mNN exclusion; observed ",
      nrow(bp_no_mnn)
    )
  )
}


srsq_sensitivity_list <- list()

result_counter <- 1L


for (i in seq_len(
  nrow(
    availability_continuous
  )
)) {
  
  if (
    !availability_continuous$
    available_for_analysis[i]
  ) {
    next
  }
  
  
  variable_name <-
    availability_continuous$variable[i]
  
  
  variable_label <-
    availability_continuous$label[i]
  
  
  clinical_values <-
    bp_no_mnn[[variable_name]]
  
  
  srsq_sensitivity_list[[result_counter]] <-
    safe_spearman(
      x =
        bp_no_mnn$SRSq,
      
      y =
        clinical_values,
      
      framework =
        "SRSq_without_mNN",
      
      clinical_variable =
        variable_name,
      
      clinical_label =
        variable_label
    )
  
  
  result_counter <-
    result_counter + 1L
}


srsq_sensitivity <- bind_rows(
  srsq_sensitivity_list
)


srsq_sensitivity$BH_sensitivity <- NA_real_


valid_sensitivity <- which(
  !is.na(
    srsq_sensitivity$p_value
  ) &
    is.finite(
      srsq_sensitivity$p_value
    )
)


if (length(valid_sensitivity) > 0) {
  
  srsq_sensitivity$BH_sensitivity[
    valid_sensitivity
  ] <- stats::p.adjust(
    srsq_sensitivity$p_value[
      valid_sensitivity
    ],
    method = "BH"
  )
}


srsq_sensitivity <- srsq_sensitivity %>%
  arrange(
    BH_sensitivity,
    p_value
  )


cat("SRSq sensitivity without mNN outliers:\n")


print(
  srsq_sensitivity %>%
    select(
      clinical_label,
      n,
      effect,
      p_value,
      BH_sensitivity
    ),
  n = Inf
)


cat("\n")


# ==============================================================================
# 25. CONSISTENCY WITH PREVIOUS SRS RESULTS
# ==============================================================================

previous_consistency <- tibble(
  clinical_variable =
    character(),
  
  observed_rho =
    numeric(),
  
  expected_rho_reference =
    numeric(),
  
  absolute_difference =
    numeric(),
  
  status =
    character()
)


crp_row <- all_tests %>%
  filter(
    framework == "SRSq",
    clinical_variable == "crp_numeric",
    test == "Spearman"
  )


if (nrow(crp_row) == 1) {
  
  crp_observed <-
    crp_row$effect[1]
  
  
  previous_consistency <- bind_rows(
    
    previous_consistency,
    
    tibble(
      clinical_variable =
        "CRP",
      
      observed_rho =
        crp_observed,
      
      expected_rho_reference =
        0.526,
      
      absolute_difference =
        abs(
          crp_observed -
            0.526
        ),
      
      status =
        ifelse(
          abs(
            crp_observed -
              0.526
          ) <= 0.05,
          "consistent",
          "review_mapping"
        )
    )
  )
}


creatinine_row <- all_tests %>%
  filter(
    framework == "SRSq",
    clinical_variable == "creatinine_numeric",
    test == "Spearman"
  )


if (nrow(creatinine_row) == 1) {
  
  creatinine_observed <-
    creatinine_row$effect[1]
  
  
  previous_consistency <- bind_rows(
    
    previous_consistency,
    
    tibble(
      clinical_variable =
        "Creatinine",
      
      observed_rho =
        creatinine_observed,
      
      expected_rho_reference =
        0.382,
      
      absolute_difference =
        abs(
          creatinine_observed -
            0.382
        ),
      
      status =
        ifelse(
          abs(
            creatinine_observed -
              0.382
          ) <= 0.05,
          "consistent",
          "review_mapping"
        )
    )
  )
}


cat(
  "Consistency with previous SRS analysis:\n"
)


print(
  previous_consistency
)


cat("\n")


# ==============================================================================
# 26. CLINICAL RATES BY SRS AND CTS
# ==============================================================================

rate_variables <- c(
  "outcome_status",
  "ventilation_status",
  "culture_status",
  "creatinine_group_110"
)


rate_labels <- c(
  outcome_status =
    "Outcome",
  
  ventilation_status =
    "Mechanical ventilation",
  
  culture_status =
    "Culture status",
  
  creatinine_group_110 =
    "Creatinine >110"
)


clinical_rates_list <- list()

rate_counter <- 1L


for (variable_name in rate_variables) {
  
  clinical_values <-
    bp_merged[[variable_name]]
  
  
  for (
    framework_name in c(
      "SRS",
      "CTS"
    )
  ) {
    
    molecular_values <-
      bp_merged[[framework_name]]
    
    
    tmp <- tibble(
      
      clinical_group =
        clinical_values,
      
      molecular_group =
        molecular_values
      
    ) %>%
      filter(
        !is.na(
          clinical_group
        ),
        !is.na(
          molecular_group
        )
      ) %>%
      count(
        clinical_group,
        molecular_group,
        name = "n"
      ) %>%
      group_by(
        clinical_group
      ) %>%
      mutate(
        clinical_group_total =
          sum(n),
        
        percent_within_clinical_group =
          100 *
          n /
          clinical_group_total
      ) %>%
      ungroup() %>%
      mutate(
        clinical_variable =
          variable_name,
        
        clinical_label =
          unname(
            rate_labels[
              variable_name
            ]
          ),
        
        framework =
          framework_name
      )
    
    
    clinical_rates_list[[rate_counter]] <-
      tmp
    
    
    rate_counter <-
      rate_counter + 1L
  }
}


clinical_rates <- bind_rows(
  clinical_rates_list
)


# ==============================================================================
# 27. FIGURE A — PRIMARY SCORE BY OUTCOME
# ==============================================================================

outcome_plot_data <- bp_merged %>%
  filter(
    !is.na(
      outcome_status
    )
  )


if (
  length(
    unique(
      outcome_plot_data$outcome_status
    )
  ) >= 2
) {
  
  p_outcome <- ggplot(
    outcome_plot_data,
    aes(
      x =
        outcome_status,
      
      y =
        primary_5gene_score
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
        "Five-gene host-response score by clinical outcome",
      
      x =
        "Clinical outcome",
      
      y =
        "Myeloid-adaptive balance score"
    ) +
    
    theme_bw(
      base_size = 12
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "136_Figure_A_primary_score_by_outcome.png"
    ),
    plot = p_outcome,
    width = 6.5,
    height = 5.5,
    dpi = 300
  )
}


# ==============================================================================
# 28. FIGURE B — PRIMARY SCORE BY VENTILATION
# ==============================================================================

ivl_plot_data <- bp_merged %>%
  filter(
    !is.na(
      ventilation_status
    )
  )


if (
  length(
    unique(
      ivl_plot_data$ventilation_status
    )
  ) >= 2
) {
  
  p_ivl <- ggplot(
    ivl_plot_data,
    aes(
      x =
        ventilation_status,
      
      y =
        primary_5gene_score
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
        "Five-gene host-response score by mechanical ventilation",
      
      x =
        "Mechanical ventilation",
      
      y =
        "Myeloid-adaptive balance score"
    ) +
    
    theme_bw(
      base_size = 12
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "136_Figure_B_primary_score_by_ventilation.png"
    ),
    plot = p_ivl,
    width = 6.5,
    height = 5.5,
    dpi = 300
  )
}


# ==============================================================================
# 29. FIGURE C — SCORE / SRSq CONTINUOUS CLINICAL CORRELATIONS
# ==============================================================================

correlation_plot_data <- all_tests %>%
  filter(
    framework %in%
      c(
        "Primary_5gene_score",
        "SRSq"
      ),
    
    test == "Spearman",
    
    !is.na(
      effect
    )
  )


if (nrow(correlation_plot_data) > 0) {
  
  p_corr <- ggplot(
    correlation_plot_data,
    aes(
      x =
        effect,
      
      y =
        reorder(
          clinical_label,
          effect
        )
    )
  ) +
    
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    
    geom_point(
      size = 3
    ) +
    
    facet_wrap(
      ~ framework
    ) +
    
    labs(
      title =
        "Clinical correlates of blood host-response scores",
      
      x =
        "Spearman rho",
      
      y =
        NULL
    ) +
    
    theme_bw(
      base_size = 11
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "136_Figure_C_score_SRSq_clinical_correlations.png"
    ),
    plot = p_corr,
    width = 9,
    height = 6.5,
    dpi = 300
  )
}


# ==============================================================================
# 30. FIGURE D — TOP ASSOCIATIONS
# ==============================================================================

overview_plot_data <- all_tests %>%
  filter(
    !is.na(
      p_value
    )
  ) %>%
  mutate(
    minus_log10_p =
      -log10(
        pmax(
          p_value,
          1e-300
        )
      ),
    
    association_name =
      paste(
        framework,
        clinical_label,
        sep = " — "
      )
  ) %>%
  arrange(
    p_value
  ) %>%
  slice_head(
    n = 25
  )


if (nrow(overview_plot_data) > 0) {
  
  p_overview <- ggplot(
    overview_plot_data,
    aes(
      x =
        minus_log10_p,
      
      y =
        reorder(
          association_name,
          minus_log10_p
        )
    )
  ) +
    
    geom_point(
      size = 2.8
    ) +
    
    geom_vline(
      xintercept =
        -log10(
          0.05
        ),
      linetype = 2
    ) +
    
    labs(
      title =
        "Top exploratory blood endotype-clinical associations",
      
      x =
        "-log10(nominal p-value)",
      
      y =
        NULL
    ) +
    
    theme_bw(
      base_size = 10
    )
  
  
  ggsave(
    filename = file.path(
      figures_dir,
      "136_Figure_D_top_clinical_associations.png"
    ),
    plot = p_overview,
    width = 10,
    height = 8,
    dpi = 300
  )
}


# ==============================================================================
# 31. SAVE CSV TABLES
# ==============================================================================

write.csv(
  mapping_table,
  file.path(
    tables_dir,
    "136_clinical_column_mapping.csv"
  ),
  row.names = FALSE
)


write.csv(
  availability_continuous,
  file.path(
    tables_dir,
    "136_continuous_clinical_availability.csv"
  ),
  row.names = FALSE
)


write.csv(
  availability_categorical,
  file.path(
    tables_dir,
    "136_categorical_clinical_availability.csv"
  ),
  row.names = FALSE
)


write.csv(
  bp_merged,
  file.path(
    tables_dir,
    "136_blood_molecular_clinical_merged.csv"
  ),
  row.names = FALSE
)


write.csv(
  all_tests,
  file.path(
    tables_dir,
    "136_all_clinical_association_tests.csv"
  ),
  row.names = FALSE
)


write.csv(
  primary_association_table,
  file.path(
    tables_dir,
    "136_primary_clinical_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  srsq_sensitivity,
  file.path(
    tables_dir,
    "136_SRSq_sensitivity_without_mNN.csv"
  ),
  row.names = FALSE
)


write.csv(
  previous_consistency,
  file.path(
    tables_dir,
    "136_previous_SRS_consistency_check.csv"
  ),
  row.names = FALSE
)


write.csv(
  clinical_rates,
  file.path(
    tables_dir,
    "136_clinical_rates_by_SRS_CTS.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 32. PUBLICATION WORKBOOK
# ==============================================================================

wb <- openxlsx::createWorkbook()


run_info <- tibble(
  
  parameter = c(
    "script",
    "run_date",
    "analysis_scope",
    "n_BP",
    "molecular_source",
    "clinical_source",
    "selected_clinical_sheet",
    "FDR_rule",
    "interpretation"
  ),
  
  
  value = c(
    script_name,
    as.character(
      run_date
    ),
    "Blood sepsis only",
    "35",
    script135_file,
    clinical_file,
    selected_clinical_sheet,
    "Global Benjamini-Hochberg FDR across all main association tests",
    "Exploratory clinical associations; no prediction claim"
  )
)


sheet_data <- list(
  
  "00_run_info" =
    run_info,
  
  "01_column_mapping" =
    mapping_table,
  
  "02_cont_availability" =
    availability_continuous,
  
  "03_cat_availability" =
    availability_categorical,
  
  "04_BP_merged" =
    bp_merged,
  
  "05_all_tests" =
    all_tests,
  
  "06_primary_ranked" =
    primary_association_table,
  
  "07_continuous_tests" =
    continuous_results_table,
  
  "08_categorical_tests" =
    categorical_results_table,
  
  "09_SRSq_sensitivity" =
    srsq_sensitivity,
  
  "10_SRS_consistency" =
    previous_consistency,
  
  "11_clinical_rates" =
    clinical_rates
)


for (
  sheet_name in names(
    sheet_data
  )
) {
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  openxlsx::writeData(
    wb,
    sheet_name,
    sheet_data[[sheet_name]]
  )
}


workbook_file <- file.path(
  tables_dir,
  "136_blood_endotype_clinical_associations.xlsx"
)


openxlsx::saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)


# ==============================================================================
# 33. SUMMARY TABLES
# ==============================================================================

global_significant <- all_tests %>%
  filter(
    !is.na(
      BH_global
    ),
    BH_global < 0.05
  ) %>%
  arrange(
    BH_global,
    p_value
  )


nominal_significant <- all_tests %>%
  filter(
    !is.na(
      p_value
    ),
    p_value < 0.05
  ) %>%
  arrange(
    p_value
  )


format_result_line <- function(row) {
  
  effect_text <- ""
  
  if (
    length(
      row$effect
    ) == 1 &&
    !is.na(
      row$effect
    )
  ) {
    
    effect_text <- paste0(
      "; ",
      row$effect_name,
      "=",
      signif(
        row$effect,
        4
      )
    )
  }
  
  return(
    paste0(
      row$framework,
      " vs ",
      row$clinical_label,
      ": ",
      row$test,
      ", n=",
      row$n,
      effect_text,
      ", p=",
      signif(
        row$p_value,
        4
      ),
      ", BH_global=",
      signif(
        row$BH_global,
        4
      )
    )
  )
}


# ==============================================================================
# 34. RUSSIAN SUMMARY
# ==============================================================================

summary_ru <- c(
  
  "SCRIPT 136 — BLOOD ENDOTYPE / CLINICAL ASSOCIATIONS",
  
  "====================================================================",
  
  "",
  
  "АНАЛИЗ",
  
  "Только blood sepsis cohort, BP n=35.",
  
  paste0(
    "Molecular source: ",
    script135_file
  ),
  
  paste0(
    "Clinical source: ",
    clinical_file
  ),
  
  paste0(
    "Clinical worksheet: ",
    selected_clinical_sheet
  ),
  
  "",
  
  "МОЛЕКУЛЯРНЫЕ ОСИ",
  
  "- primary five-gene host-response score",
  
  "- SRSq",
  
  "- SRS1 vs SRS2",
  
  "- CTS1 / CTS2 / CTS3",
  
  "",
  
  "СТАТИСТИКА",
  
  "Основной multiple-testing control: global Benjamini-Hochberg FDR.",
  
  "",
  
  paste0(
    "Main association tests: ",
    nrow(
      all_tests
    )
  ),
  
  paste0(
    "Nominal p<0.05: ",
    nrow(
      nominal_significant
    )
  ),
  
  paste0(
    "Global BH FDR<0.05: ",
    nrow(
      global_significant
    )
  ),
  
  "",
  
  "GLOBAL FDR-SIGNIFICANT ASSOCIATIONS:"
)


if (nrow(global_significant) == 0) {
  
  summary_ru <- c(
    summary_ru,
    "None."
  )
  
} else {
  
  for (
    i in seq_len(
      min(
        20,
        nrow(
          global_significant
        )
      )
    )
  ) {
    
    summary_ru <- c(
      summary_ru,
      format_result_line(
        global_significant[
          i,
        ]
      )
    )
  }
}


summary_ru <- c(
  summary_ru,
  "",
  "TOP NOMINAL ASSOCIATIONS:"
)


if (nrow(nominal_significant) == 0) {
  
  summary_ru <- c(
    summary_ru,
    "None."
  )
  
} else {
  
  for (
    i in seq_len(
      min(
        20,
        nrow(
          nominal_significant
        )
      )
    )
  ) {
    
    summary_ru <- c(
      summary_ru,
      format_result_line(
        nominal_significant[
          i,
        ]
      )
    )
  }
}


summary_ru <- c(
  
  summary_ru,
  
  "",
  
  "INTERPRETATION",
  
  paste0(
    "Clinical associations are exploratory and hypothesis-generating. ",
    "No mortality or ventilation prediction model was constructed."
  ),
  
  "",
  
  paste0(
    "The primary five-gene score remains a candidate endotype-informed ",
    "molecular signature rather than a clinically validated prognostic assay."
  )
)


summary_ru_file <- file.path(
  text_dir,
  "136_summary_RU.txt"
)


writeLines(
  summary_ru,
  summary_ru_file
)


# ==============================================================================
# 35. ENGLISH SUMMARY
# ==============================================================================

summary_en <- c(
  
  "SCRIPT 136 — BLOOD ENDOTYPE / CLINICAL ASSOCIATIONS",
  
  "====================================================================",
  
  "",
  
  "Blood sepsis cohort only, BP n=35.",
  
  "",
  
  paste0(
    "Main association tests: ",
    nrow(
      all_tests
    )
  ),
  
  paste0(
    "Nominal p<0.05: ",
    nrow(
      nominal_significant
    )
  ),
  
  paste0(
    "Global BH FDR<0.05: ",
    nrow(
      global_significant
    )
  ),
  
  "",
  
  paste0(
    "Clinical associations are exploratory and hypothesis-generating; ",
    "no mortality or ventilation prediction model was developed."
  ),
  
  "",
  
  paste0(
    "The primary five-gene score remains a candidate endotype-informed ",
    "molecular signature requiring independent validation."
  )
)


summary_en_file <- file.path(
  text_dir,
  "136_summary_EN.txt"
)


writeLines(
  summary_en,
  summary_en_file
)


# ==============================================================================
# 36. REPRODUCIBILITY MANIFEST
# ==============================================================================

input_manifest <- tibble(
  
  input = c(
    "Script135_BP_scores",
    "clinical_metadata"
  ),
  
  path = c(
    script135_file,
    clinical_file
  )
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
    "136_input_file_manifest.csv"
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
    "136_sessionInfo.txt"
  )
)


# ==============================================================================
# 38. FINAL CONSOLE REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 136 completed successfully.\n")
cat("====================================================================\n\n")


cat("Clinical metadata source:\n")

cat(
  normalizePath(
    clinical_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat("Clinical worksheet:\n")

cat(
  selected_clinical_sheet,
  "\n\n"
)


cat("Main output folder:\n")

cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat("Key numbers:\n")


cat(
  "Main association tests: ",
  nrow(
    all_tests
  ),
  "\n",
  sep = ""
)


cat(
  "Nominal p<0.05: ",
  nrow(
    nominal_significant
  ),
  "\n",
  sep = ""
)


cat(
  "Global BH FDR<0.05: ",
  nrow(
    global_significant
  ),
  "\n\n",
  sep = ""
)


cat("Main outputs:\n")


cat(
  "1) ",
  normalizePath(
    workbook_file,
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
      "136_all_clinical_association_tests.csv"
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
      "136_primary_clinical_associations.csv"
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
      "136_SRSq_sensitivity_without_mNN.csv"
    ),
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat(
  "5) ",
  normalizePath(
    summary_ru_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat("Open first in workbook:\n")
cat("  01_column_mapping\n")
cat("  02_cont_availability\n")
cat("  03_cat_availability\n")
cat("  05_all_tests\n")
cat("  06_primary_ranked\n")
cat("  09_SRSq_sensitivity\n")
cat("  10_SRS_consistency\n")
cat("  11_clinical_rates\n\n")


cat("IMPORTANT:\n")
cat("- BP n=35 only.\n")
cat("- Primary five-gene panel frozen from Script 135.\n")
cat("- No new feature selection.\n")
cat("- No mortality prediction model.\n")
cat("- No ventilation prediction model.\n")
cat("- Main multiplicity correction = BH_global.\n")
cat("- Nominal associations are exploratory.\n")
cat("- No urine.\n")
cat("- No lncRNA.\n\n")


cat("Done.\n")