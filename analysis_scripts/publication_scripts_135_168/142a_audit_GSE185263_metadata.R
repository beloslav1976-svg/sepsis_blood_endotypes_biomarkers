# ==============================================================================
# Script 142a
# Audit of GSE185263 metadata and raw-count matrix structure
#
# Project: Sepsis_DESeq2
#
# PURPOSE
# Before external validation of the frozen five-gene host-response signature,
# inspect the independent GSE185263 cohort WITHOUT examining expression of the
# five frozen genes and WITHOUT selecting an endpoint based on expression data.
#
# This script:
#   1. downloads GEO phenotype metadata;
#   2. parses GEO characteristics_ch1 fields into structured variables;
#   3. audits disease state, age, sex, collection site/location,
#      SOFA, mortality and any other available clinical fields;
#   4. audits sample-name prefixes / possible cohort structure;
#   5. detects possible repeated/time-point samples;
#   6. searches metadata ONLY for published endotype labels
#      NPS, INF, IHD, IFN, ADA;
#   7. downloads the published raw-count file;
#   8. reads ONLY its header for sample-column mapping;
#   9. verifies how many of the 392 GEO samples map to the count matrix;
#  10. generates an Excel/CSV/text audit package.
#
# IMPORTANT
#   - NO five-gene expression analysis
#   - NO CD177/HK3/IRAK3/CARD11/IKZF2 lookup
#   - NO score calculation
#   - NO feature selection
#   - NO AUC
#   - NO outcome model
#   - NO endpoint selection from expression results
#
# The purpose is to determine the correct PREDECLARED design for Script 142b.
# ==============================================================================


# ==============================================================================
# 0. GLOBAL SETTINGS
# ==============================================================================

options(
  stringsAsFactors = FALSE
)


project_dir <- Sys.getenv("SEPSIS_PROJECT_DIR", unset = path.expand("~/Sepsis_DESeq2"))

gse_id <- "GSE185263"

script_name <- "142a_audit_GSE185263_metadata.R"

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
cat("Running Script 142a\n")
cat("GSE185263 metadata and cohort audit\n")
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
  "142a_GSE185263_metadata_audit"
)


raw_dir <- file.path(
  output_dir,
  "raw_download"
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
    raw_dir,
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


if (length(missing_cran) > 0) {
  
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
  "Biobase"
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
  
  library(openxlsx)
  
  library(GEOquery)
  
  library(Biobase)
})


cat(
  "Required packages loaded successfully.\n\n"
)


# ==============================================================================
# 3. HELPER FUNCTIONS
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


clean_variable_name <- function(x) {
  
  x <- as.character(
    x
  )
  
  x <- trimws(
    x
  )
  
  x <- tolower(
    x
  )
  
  x <- gsub(
    "[^a-z0-9]+",
    "_",
    x
  )
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  return(
    x
  )
}


missing_like <- function(x) {
  
  x_chr <- trimws(
    as.character(
      x
    )
  )
  
  
  x_lower <- tolower(
    x_chr
  )
  
  
  is.na(x) |
    x_chr == "" |
    x_lower %in%
    c(
      "na",
      "n/a",
      "nan",
      "null",
      "none",
      "not available",
      "not applicable",
      "unknown",
      "missing"
    )
}


safe_examples <- function(
    x,
    n = 8,
    max_chars = 500
) {
  
  x_chr <- as.character(
    x
  )
  
  
  x_chr <- x_chr[
    !missing_like(
      x_chr
    )
  ]
  
  
  x_chr <- unique(
    x_chr
  )
  
  
  if (length(x_chr) == 0) {
    
    return(
      ""
    )
  }
  
  
  x_chr <- head(
    x_chr,
    n
  )
  
  
  out <- paste(
    x_chr,
    collapse = " | "
  )
  
  
  if (
    nchar(
      out
    ) > max_chars
  ) {
    
    out <- paste0(
      substr(
        out,
        1,
        max_chars
      ),
      "..."
    )
  }
  
  
  return(
    out
  )
}


audit_dataframe_columns <- function(df) {
  
  result_list <- lapply(
    names(
      df
    ),
    function(column_name) {
      
      x <- df[[column_name]]
      
      
      missing_flag <- missing_like(
        x
      )
      
      
      x_nonmissing <- x[
        !missing_flag
      ]
      
      
      unique_values <- unique(
        as.character(
          x_nonmissing
        )
      )
      
      
      numeric_candidate <- suppressWarnings(
        as.numeric(
          as.character(
            x_nonmissing
          )
        )
      )
      
      
      numeric_fraction <- if (
        length(
          x_nonmissing
        ) == 0
      ) {
        
        NA_real_
        
      } else {
        
        mean(
          is.finite(
            numeric_candidate
          )
        )
      }
      
      
      tibble::tibble(
        
        variable =
          column_name,
        
        n_total =
          length(
            x
          ),
        
        n_nonmissing =
          sum(
            !missing_flag
          ),
        
        n_missing =
          sum(
            missing_flag
          ),
        
        percent_available =
          round(
            100 *
              sum(
                !missing_flag
              ) /
              length(
                x
              ),
            1
          ),
        
        n_unique =
          length(
            unique_values
          ),
        
        numeric_fraction =
          numeric_fraction,
        
        examples =
          safe_examples(
            x
          )
      )
    }
  )
  
  
  return(
    dplyr::bind_rows(
      result_list
    )
  )
}


extract_sample_prefix <- function(x) {
  
  x <- as.character(
    x
  )
  
  
  prefix <- stringr::str_extract(
    x,
    "^[A-Za-z]+"
  )
  
  
  prefix <- tolower(
    prefix
  )
  
  
  return(
    prefix
  )
}


extract_timepoint_suffix <- function(x) {
  
  x <- as.character(
    x
  )
  
  
  suffix <- stringr::str_extract(
    x,
    "(?i)(T[0-9]+[A-Za-z]?|W[0-9]+)$"
  )
  
  
  return(
    suffix
  )
}


remove_timepoint_suffix <- function(x) {
  
  x <- as.character(
    x
  )
  
  
  out <- stringr::str_remove(
    x,
    "(?i)(T[0-9]+[A-Za-z]?|W[0-9]+)$"
  )
  
  
  return(
    out
  )
}


make_low_cardinality_distribution <- function(
    df,
    variable_name,
    max_unique = 30
) {
  
  x <- df[[variable_name]]
  
  
  x_chr <- as.character(
    x
  )
  
  
  x_chr[missing_like(
    x_chr
  )] <- NA_character_
  
  
  n_unique <- length(
    unique(
      x_chr[
        !is.na(
          x_chr
        )
      ]
    )
  )
  
  
  if (
    n_unique == 0 ||
    n_unique > max_unique
  ) {
    
    return(
      NULL
    )
  }
  
  
  out <- tibble::tibble(
    
    value =
      x_chr
  ) %>%
    
    dplyr::count(
      value,
      name = "n",
      .drop = FALSE
    ) %>%
    
    dplyr::mutate(
      
      variable =
        variable_name,
      
      value =
        ifelse(
          is.na(
            value
          ),
          "<MISSING>",
          value
        )
    ) %>%
    
    dplyr::select(
      variable,
      value,
      n
    )
  
  
  return(
    out
  )
}


numeric_summary_if_possible <- function(
    df,
    variable_name,
    min_numeric_fraction = 0.80
) {
  
  x <- df[[variable_name]]
  
  
  missing_flag <- missing_like(
    x
  )
  
  
  x_nonmissing <- x[
    !missing_flag
  ]
  
  
  if (
    length(
      x_nonmissing
    ) < 5
  ) {
    
    return(
      NULL
    )
  }
  
  
  x_numeric <- suppressWarnings(
    as.numeric(
      as.character(
        x_nonmissing
      )
    )
  )
  
  
  numeric_fraction <- mean(
    is.finite(
      x_numeric
    )
  )
  
  
  if (
    !is.finite(
      numeric_fraction
    ) ||
    numeric_fraction <
    min_numeric_fraction
  ) {
    
    return(
      NULL
    )
  }
  
  
  x_numeric <- x_numeric[
    is.finite(
      x_numeric
    )
  ]
  
  
  return(
    tibble::tibble(
      
      variable =
        variable_name,
      
      n =
        length(
          x_numeric
        ),
      
      min =
        min(
          x_numeric
        ),
      
      q1 =
        stats::quantile(
          x_numeric,
          0.25
        ),
      
      median =
        stats::median(
          x_numeric
        ),
      
      mean =
        mean(
          x_numeric
        ),
      
      q3 =
        stats::quantile(
          x_numeric,
          0.75
        ),
      
      max =
        max(
          x_numeric
        )
    )
  )
}


# ==============================================================================
# 4. DOWNLOAD GEO SERIES METADATA
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
    FALSE,
  
  check.names =
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


if (
  nrow(
    pheno_raw
  ) != 392
) {
  
  warning(
    paste0(
      "Expected 392 GEO samples from the Series record; observed ",
      nrow(
        pheno_raw
      ),
      "."
    )
  )
}


# ==============================================================================
# 5. CORE SAMPLE IDENTIFIERS
# ==============================================================================

if (
  "geo_accession" %in%
  names(
    pheno_raw
  )
) {
  
  geo_accession_vector <- as.character(
    pheno_raw[["geo_accession"]]
  )
  
} else {
  
  geo_accession_vector <- as.character(
    pheno_raw[[".geo_rowname"]]
  )
}


if (
  "title" %in%
  names(
    pheno_raw
  )
) {
  
  sample_title_vector <- as.character(
    pheno_raw[["title"]]
  )
  
} else {
  
  sample_title_vector <-
    geo_accession_vector
}


if (
  "source_name_ch1" %in%
  names(
    pheno_raw
  )
) {
  
  source_name_vector <- as.character(
    pheno_raw[["source_name_ch1"]]
  )
  
} else {
  
  source_name_vector <- rep(
    NA_character_,
    nrow(
      pheno_raw
    )
  )
}


metadata_core <- tibble::tibble(
  
  geo_accession =
    geo_accession_vector,
  
  sample_title =
    sample_title_vector,
  
  source_name =
    source_name_vector,
  
  sample_prefix =
    extract_sample_prefix(
      sample_title_vector
    ),
  
  timepoint_suffix =
    extract_timepoint_suffix(
      sample_title_vector
    ),
  
  participant_id_guess =
    remove_timepoint_suffix(
      sample_title_vector
    )
)


# ==============================================================================
# 6. AUDIT SAMPLE-NAME / PREFIX STRUCTURE
# ==============================================================================

sample_prefix_summary <- metadata_core %>%
  
  dplyr::count(
    sample_prefix,
    name = "n_samples"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      n_samples
    )
  )


cat(
  "Sample-title prefix distribution:\n"
)


print(
  sample_prefix_summary,
  n = Inf
)


cat("\n")


timepoint_summary <- metadata_core %>%
  
  dplyr::count(
    timepoint_suffix,
    name = "n_samples",
    .drop = FALSE
  ) %>%
  
  dplyr::mutate(
    
    timepoint_suffix =
      ifelse(
        is.na(
          timepoint_suffix
        ),
        "<NO_SUFFIX>",
        timepoint_suffix
      )
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      n_samples
    )
  )


cat(
  "Detected sample-title timepoint suffixes:\n"
)


print(
  timepoint_summary,
  n = Inf
)


cat("\n")


repeated_participant_guesses <- metadata_core %>%
  
  dplyr::count(
    participant_id_guess,
    name = "n_samples"
  ) %>%
  
  dplyr::filter(
    n_samples > 1
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      n_samples
    ),
    participant_id_guess
  )


repeated_sample_details <- metadata_core %>%
  
  dplyr::semi_join(
    
    repeated_participant_guesses,
    
    by =
      "participant_id_guess"
  ) %>%
  
  dplyr::arrange(
    participant_id_guess,
    sample_title
  )


cat(
  "Possible repeated/time-point participants based on sample-title suffix parsing: ",
  nrow(
    repeated_participant_guesses
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 7. IDENTIFY GEO CHARACTERISTICS COLUMNS
# ==============================================================================

characteristic_cols <- names(
  pheno_raw
)[
  grepl(
    "^characteristics_ch1",
    names(
      pheno_raw
    ),
    ignore.case = TRUE
  )
]


cat(
  "GEO characteristics columns detected: ",
  length(
    characteristic_cols
  ),
  "\n",
  sep = ""
)


if (
  length(
    characteristic_cols
  ) > 0
) {
  
  print(
    characteristic_cols
  )
}


cat("\n")


if (
  length(
    characteristic_cols
  ) == 0
) {
  
  stop(
    "No characteristics_ch1 columns found in GEO metadata."
  )
}


# ==============================================================================
# 8. PARSE CHARACTERISTICS INTO LONG FORMAT
# ==============================================================================

characteristics_source <- pheno_raw %>%
  
  dplyr::mutate(
    
    geo_accession =
      geo_accession_vector,
    
    sample_title =
      sample_title_vector
  ) %>%
  
  dplyr::select(
    
    geo_accession,
    
    sample_title,
    
    dplyr::all_of(
      characteristic_cols
    )
  )


characteristics_long <- characteristics_source %>%
  
  tidyr::pivot_longer(
    
    cols =
      dplyr::all_of(
        characteristic_cols
      ),
    
    names_to =
      "characteristic_field",
    
    values_to =
      "characteristic_raw"
  ) %>%
  
  dplyr::mutate(
    
    characteristic_raw =
      trimws(
        as.character(
          characteristic_raw
        )
      ),
    
    has_colon =
      grepl(
        ":",
        characteristic_raw,
        fixed = TRUE
      ),
    
    characteristic_key =
      ifelse(
        
        has_colon,
        
        trimws(
          sub(
            ":.*$",
            "",
            characteristic_raw
          )
        ),
        
        characteristic_field
      ),
    
    characteristic_value =
      ifelse(
        
        has_colon,
        
        trimws(
          sub(
            "^[^:]+:[[:space:]]*",
            "",
            characteristic_raw
          )
        ),
        
        characteristic_raw
      ),
    
    characteristic_key_clean =
      clean_variable_name(
        characteristic_key
      )
  ) %>%
  
  dplyr::filter(
    !missing_like(
      characteristic_value
    )
  )


cat(
  "Parsed non-missing GEO characteristic entries: ",
  nrow(
    characteristics_long
  ),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 9. CHARACTERISTIC KEY AUDIT
# ==============================================================================

characteristic_key_summary <- characteristics_long %>%
  
  dplyr::group_by(
    characteristic_key_clean
  ) %>%
  
  dplyr::summarise(
    
    original_key_examples =
      safe_examples(
        characteristic_key,
        n = 5
      ),
    
    n_samples =
      dplyr::n_distinct(
        geo_accession
      ),
    
    n_unique_values =
      dplyr::n_distinct(
        characteristic_value
      ),
    
    value_examples =
      safe_examples(
        characteristic_value,
        n = 8
      ),
    
    .groups =
      "drop"
  ) %>%
  
  dplyr::arrange(
    characteristic_key_clean
  )


cat(
  "Parsed clinical/characteristic variables:\n"
)


print(
  characteristic_key_summary,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 10. CONVERT CHARACTERISTICS TO WIDE SAMPLE METADATA
# ==============================================================================

characteristics_wide <- characteristics_long %>%
  
  dplyr::group_by(
    geo_accession,
    characteristic_key_clean
  ) %>%
  
  dplyr::summarise(
    
    characteristic_value =
      paste(
        unique(
          characteristic_value
        ),
        collapse = " | "
      ),
    
    .groups =
      "drop"
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from =
      characteristic_key_clean,
    
    values_from =
      characteristic_value
  )


metadata_parsed <- metadata_core %>%
  
  dplyr::left_join(
    
    characteristics_wide,
    
    by =
      "geo_accession"
  )


cat(
  "Structured metadata table: ",
  nrow(
    metadata_parsed
  ),
  " samples x ",
  ncol(
    metadata_parsed
  ),
  " variables\n\n",
  sep = ""
)


# ==============================================================================
# 11. COMPLETE VARIABLE AUDIT
# ==============================================================================

metadata_variable_audit <- audit_dataframe_columns(
  metadata_parsed
)


candidate_pattern <- paste0(
  "disease|sepsis|severity|shock|infection|",
  "outcome|mortality|death|surviv|",
  "sofa|apache|organ|icu|emergency|",
  "site|location|cohort|center|centre|",
  "age|sex|gender|",
  "endotype|subtype|cluster|",
  "nps|ihd|ifn|ada|inf|",
  "time|visit|day|admission"
)


candidate_metadata_variables <- metadata_variable_audit %>%
  
  dplyr::filter(
    
    grepl(
      candidate_pattern,
      variable,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::arrange(
    variable
  )


cat(
  "Candidate clinical / severity / outcome / endotype variables:\n"
)


print(
  candidate_metadata_variables,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 12. LOW-CARDINALITY DISTRIBUTIONS
# ==============================================================================

distribution_list <- list()


distribution_counter <- 1L


for (
  variable_name in candidate_metadata_variables$variable
) {
  
  distribution_table <- make_low_cardinality_distribution(
    
    metadata_parsed,
    
    variable_name = variable_name,
    
    max_unique = 30
  )
  
  
  if (!is.null(
    distribution_table
  )) {
    
    distribution_list[[distribution_counter]] <-
      distribution_table
    
    
    distribution_counter <-
      distribution_counter +
      1L
  }
}


if (
  length(
    distribution_list
  ) > 0
) {
  
  candidate_variable_distributions <- dplyr::bind_rows(
    distribution_list
  )
  
} else {
  
  candidate_variable_distributions <- tibble::tibble(
    
    variable =
      character(),
    
    value =
      character(),
    
    n =
      integer()
  )
}


cat(
  "Low-cardinality candidate-variable distributions:\n"
)


if (
  nrow(
    candidate_variable_distributions
  ) > 0
) {
  
  print(
    candidate_variable_distributions,
    n = Inf,
    width = Inf
  )
  
} else {
  
  cat(
    "None detected.\n"
  )
}


cat("\n")


# ==============================================================================
# 13. NUMERIC VARIABLE AUDIT
# ==============================================================================

numeric_summary_list <- list()


numeric_counter <- 1L


for (
  variable_name in names(
    metadata_parsed
  )
) {
  
  numeric_table <- numeric_summary_if_possible(
    
    metadata_parsed,
    
    variable_name =
      variable_name,
    
    min_numeric_fraction =
      0.80
  )
  
  
  if (!is.null(
    numeric_table
  )) {
    
    numeric_summary_list[[numeric_counter]] <-
      numeric_table
    
    
    numeric_counter <-
      numeric_counter +
      1L
  }
}


if (
  length(
    numeric_summary_list
  ) > 0
) {
  
  numeric_variable_summary <- dplyr::bind_rows(
    numeric_summary_list
  )
  
} else {
  
  numeric_variable_summary <- tibble::tibble(
    
    variable =
      character(),
    
    n =
      integer(),
    
    min =
      numeric(),
    
    q1 =
      numeric(),
    
    median =
      numeric(),
    
    mean =
      numeric(),
    
    q3 =
      numeric(),
    
    max =
      numeric()
  )
}


cat(
  "Numeric metadata variables:\n"
)


print(
  numeric_variable_summary,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 14. SPECIFIC PUBLISHED ENDOTYPE-LABEL METADATA AUDIT
#
# IMPORTANT:
# Search is performed ONLY in metadata.
# No transcriptomic endotype is inferred or reconstructed here.
# ==============================================================================

published_endotype_labels <- c(
  "NPS",
  "INF",
  "IHD",
  "IFN",
  "ADA"
)


metadata_character_strings <- characteristics_long %>%
  
  dplyr::select(
    geo_accession,
    sample_title,
    characteristic_key_clean,
    characteristic_value
  )


endotype_hit_list <- list()


endotype_counter <- 1L


for (
  endotype_label in published_endotype_labels
) {
  
  pattern <- paste0(
    "(^|[^A-Za-z])",
    endotype_label,
    "([^A-Za-z]|$)"
  )
  
  
  hits <- metadata_character_strings %>%
    
    dplyr::filter(
      
      grepl(
        pattern,
        characteristic_value,
        ignore.case = TRUE,
        perl = TRUE
      )
    ) %>%
    
    dplyr::mutate(
      published_endotype_label =
        endotype_label
    )
  
  
  if (
    nrow(
      hits
    ) > 0
  ) {
    
    endotype_hit_list[[endotype_counter]] <-
      hits
    
    
    endotype_counter <-
      endotype_counter +
      1L
  }
}


if (
  length(
    endotype_hit_list
  ) > 0
) {
  
  published_endotype_metadata_hits <- dplyr::bind_rows(
    endotype_hit_list
  ) %>%
    
    dplyr::select(
      published_endotype_label,
      geo_accession,
      sample_title,
      characteristic_key_clean,
      characteristic_value
    )
  
} else {
  
  published_endotype_metadata_hits <- tibble::tibble(
    
    published_endotype_label =
      character(),
    
    geo_accession =
      character(),
    
    sample_title =
      character(),
    
    characteristic_key_clean =
      character(),
    
    characteristic_value =
      character()
  )
}


cat(
  "Published endotype labels found directly in GEO metadata: ",
  nrow(
    published_endotype_metadata_hits
  ),
  " entries\n\n",
  sep = ""
)


if (
  nrow(
    published_endotype_metadata_hits
  ) > 0
) {
  
  print(
    published_endotype_metadata_hits,
    n = Inf,
    width = Inf
  )
  
  
  cat("\n")
}


# ==============================================================================
# 15. DISEASE-STATE / SITE / SOFA / MORTALITY KEY SEARCH
# ==============================================================================

priority_variable_patterns <- tibble::tibble(
  
  concept = c(
    "Disease state",
    "Age",
    "Sex",
    "Collection location",
    "Collection site",
    "SOFA",
    "Mortality",
    "Outcome",
    "Endotype",
    "Severity"
  ),
  
  pattern = c(
    "disease.*state",
    "^age$|age",
    "^sex$|gender",
    "collection.*location|location",
    "collection.*site|site",
    "sofa",
    "mortality|death|surviv",
    "outcome",
    "endotype|subtype|cluster",
    "severity|severe"
  )
)


priority_variable_matches <- lapply(
  seq_len(
    nrow(
      priority_variable_patterns
    )
  ),
  function(i) {
    
    concept_i <-
      priority_variable_patterns$concept[i]
    
    
    pattern_i <-
      priority_variable_patterns$pattern[i]
    
    
    matches <- metadata_variable_audit %>%
      
      dplyr::filter(
        
        grepl(
          pattern_i,
          variable,
          ignore.case = TRUE
        )
      )
    
    
    if (
      nrow(
        matches
      ) == 0
    ) {
      
      return(
        tibble::tibble(
          
          concept =
            concept_i,
          
          variable =
            NA_character_,
          
          n_nonmissing =
            NA_integer_,
          
          percent_available =
            NA_real_,
          
          n_unique =
            NA_integer_,
          
          examples =
            ""
        )
      )
    }
    
    
    matches %>%
      
      dplyr::mutate(
        concept =
          concept_i
      ) %>%
      
      dplyr::select(
        concept,
        variable,
        n_nonmissing,
        percent_available,
        n_unique,
        examples
      )
  }
)


priority_variable_matches <- dplyr::bind_rows(
  priority_variable_matches
)


cat(
  "Priority clinical-variable availability:\n"
)


print(
  priority_variable_matches,
  n = Inf,
  width = Inf
)


cat("\n")


# ==============================================================================
# 16. SAVE GEO METADATA TABLES
# ==============================================================================

write.csv(
  
  pheno_raw,
  
  file.path(
    tables_dir,
    "142a_GSE185263_GEO_pData_raw.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  metadata_core,
  
  file.path(
    tables_dir,
    "142a_sample_core_identifiers.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  characteristics_long,
  
  file.path(
    tables_dir,
    "142a_characteristics_long.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  characteristic_key_summary,
  
  file.path(
    tables_dir,
    "142a_characteristic_key_summary.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  metadata_parsed,
  
  file.path(
    tables_dir,
    "142a_GSE185263_metadata_parsed.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  metadata_variable_audit,
  
  file.path(
    tables_dir,
    "142a_metadata_variable_audit.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  candidate_metadata_variables,
  
  file.path(
    tables_dir,
    "142a_candidate_clinical_variables.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  candidate_variable_distributions,
  
  file.path(
    tables_dir,
    "142a_candidate_variable_distributions.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  numeric_variable_summary,
  
  file.path(
    tables_dir,
    "142a_numeric_variable_summary.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  priority_variable_matches,
  
  file.path(
    tables_dir,
    "142a_priority_variable_matches.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  published_endotype_metadata_hits,
  
  file.path(
    tables_dir,
    "142a_published_endotype_metadata_hits.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  sample_prefix_summary,
  
  file.path(
    tables_dir,
    "142a_sample_prefix_summary.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  timepoint_summary,
  
  file.path(
    tables_dir,
    "142a_timepoint_suffix_summary.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  repeated_participant_guesses,
  
  file.path(
    tables_dir,
    "142a_possible_repeated_participants.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  repeated_sample_details,
  
  file.path(
    tables_dir,
    "142a_possible_repeated_sample_details.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 17. DOWNLOAD RAW-COUNT MATRIX
#
# IMPORTANT:
# We inspect ONLY the matrix HEADER / column names here.
# No gene rows are analyzed in Script 142a.
# ==============================================================================

raw_counts_filename <-
  "GSE185263_raw_counts.csv.gz"


raw_counts_file <- file.path(
  raw_dir,
  raw_counts_filename
)


if (
  !file.exists(
    raw_counts_file
  )
) {
  
  cat(
    "Downloading GSE185263 supplementary files...\n"
  )
  
  
  GEOquery::getGEOSuppFiles(
    
    gse_id,
    
    makeDirectory =
      TRUE,
    
    baseDir =
      raw_dir
  )
  
  
  candidate_count_files <- list.files(
    
    raw_dir,
    
    pattern =
      "GSE185263_raw_counts\\.csv\\.gz$",
    
    recursive =
      TRUE,
    
    full.names =
      TRUE
  )
  
  
  if (
    length(
      candidate_count_files
    ) == 0
  ) {
    
    stop(
      paste0(
        "Could not locate ",
        raw_counts_filename,
        " after GEO download."
      )
    )
  }
  
  
  candidate_file <-
    candidate_count_files[1]
  
  
  candidate_norm <- normalizePath(
    candidate_file,
    winslash = "/",
    mustWork = FALSE
  )
  
  
  target_norm <- normalizePath(
    raw_counts_file,
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
        raw_counts_file,
      
      overwrite =
        TRUE
    )
    
    
    if (!isTRUE(
      copy_ok
    )) {
      
      stop(
        "Unable to copy raw-count matrix into audit raw_download folder."
      )
    }
  }
}


if (
  !file.exists(
    raw_counts_file
  )
) {
  
  stop(
    paste0(
      "Raw-count matrix not found: ",
      raw_counts_file
    )
  )
}


cat(
  "Raw-count matrix:\n"
)


cat(
  normalizePath(
    raw_counts_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)


cat(
  "File size: ",
  round(
    file.info(
      raw_counts_file
    )$size /
      1024^2,
    2
  ),
  " MB\n\n",
  sep = ""
)


# ==============================================================================
# 18. READ ONLY RAW-COUNT MATRIX HEADER
# ==============================================================================

cat(
  "Reading raw-count matrix header only...\n"
)


counts_header <- tryCatch(
  
  {
    
    data.table::fread(
      
      raw_counts_file,
      
      nrows =
        0,
      
      data.table =
        FALSE,
      
      check.names =
        FALSE
    )
  },
  
  error = function(e) {
    
    message(
      paste0(
        "fread(nrows=0) failed: ",
        e$message,
        "\nReading one row only for column-name audit."
      )
    )
    
    
    data.table::fread(
      
      raw_counts_file,
      
      nrows =
        1,
      
      data.table =
        FALSE,
      
      check.names =
        FALSE
    )
  }
)


count_column_names <- names(
  counts_header
)


cat(
  "Raw-count file contains ",
  length(
    count_column_names
  ),
  " columns.\n\n",
  sep = ""
)


cat(
  "First count-matrix column names:\n"
)


print(
  head(
    count_column_names,
    20
  )
)


cat("\n")


# ==============================================================================
# 19. MAP COUNT-MATRIX COLUMNS TO GEO SAMPLES
# ==============================================================================

sample_lookup <- metadata_core %>%
  
  dplyr::mutate(
    
    accession_key =
      normalize_key(
        geo_accession
      ),
    
    title_key =
      normalize_key(
        sample_title
      )
  )


title_to_accession <- stats::setNames(
  sample_lookup$geo_accession,
  sample_lookup$title_key
)


accession_to_accession <- stats::setNames(
  sample_lookup$geo_accession,
  sample_lookup$accession_key
)


count_column_mapping <- tibble::tibble(
  
  count_column =
    count_column_names,
  
  count_column_key =
    normalize_key(
      count_column_names
    )
) %>%
  
  dplyr::mutate(
    
    geo_from_title =
      unname(
        title_to_accession[
          count_column_key
        ]
      ),
    
    geo_from_accession =
      unname(
        accession_to_accession[
          count_column_key
        ]
      ),
    
    geo_accession =
      dplyr::coalesce(
        geo_from_title,
        geo_from_accession
      )
  ) %>%
  
  dplyr::left_join(
    
    metadata_core %>%
      
      dplyr::select(
        geo_accession,
        sample_title,
        sample_prefix,
        timepoint_suffix,
        participant_id_guess
      ),
    
    by =
      "geo_accession"
  )


mapped_count_columns <- count_column_mapping %>%
  
  dplyr::filter(
    !is.na(
      geo_accession
    )
  )


unmapped_count_columns <- count_column_mapping %>%
  
  dplyr::filter(
    is.na(
      geo_accession
    )
  )


geo_samples_missing_from_counts <- metadata_core %>%
  
  dplyr::filter(
    !geo_accession %in%
      mapped_count_columns$geo_accession
  )


cat(
  "Count-matrix columns mapped to GEO samples: ",
  nrow(
    mapped_count_columns
  ),
  "\n",
  sep = ""
)


cat(
  "Non-sample / unmapped count-matrix columns: ",
  nrow(
    unmapped_count_columns
  ),
  "\n",
  sep = ""
)


cat(
  "GEO samples missing from count matrix: ",
  nrow(
    geo_samples_missing_from_counts
  ),
  "\n\n",
  sep = ""
)


if (
  nrow(
    mapped_count_columns
  ) != 392
) {
  
  warning(
    paste0(
      "Expected 392 GEO sample columns to map to raw counts; mapped ",
      nrow(
        mapped_count_columns
      ),
      ". Inspect mapping before validation."
    )
  )
}


cat(
  "Unmapped count-matrix columns:\n"
)


print(
  unmapped_count_columns %>%
    
    dplyr::select(
      count_column
    ),
  n = Inf
)


cat("\n")


if (
  nrow(
    geo_samples_missing_from_counts
  ) > 0
) {
  
  cat(
    "GEO samples not mapped to count matrix:\n"
  )
  
  
  print(
    geo_samples_missing_from_counts,
    n = Inf
  )
  
  
  cat("\n")
}


write.csv(
  
  count_column_mapping,
  
  file.path(
    tables_dir,
    "142a_raw_count_column_mapping.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  unmapped_count_columns,
  
  file.path(
    tables_dir,
    "142a_unmapped_raw_count_columns.csv"
  ),
  
  row.names =
    FALSE
)


write.csv(
  
  geo_samples_missing_from_counts,
  
  file.path(
    tables_dir,
    "142a_GEO_samples_missing_from_raw_counts.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 20. COUNT-MATRIX STRUCTURE SUMMARY
# ==============================================================================

count_matrix_structure <- tibble::tibble(
  
  metric = c(
    "GEO_samples",
    "count_file_columns_total",
    "mapped_GEO_sample_columns",
    "unmapped_count_columns",
    "GEO_samples_missing_from_counts",
    "raw_count_file_size_MB"
  ),
  
  value = c(
    nrow(
      metadata_core
    ),
    length(
      count_column_names
    ),
    nrow(
      mapped_count_columns
    ),
    nrow(
      unmapped_count_columns
    ),
    nrow(
      geo_samples_missing_from_counts
    ),
    round(
      file.info(
        raw_counts_file
      )$size /
        1024^2,
      2
    )
  )
)


cat(
  "Raw-count matrix structure audit:\n"
)


print(
  count_matrix_structure,
  n = Inf
)


cat("\n")


write.csv(
  
  count_matrix_structure,
  
  file.path(
    tables_dir,
    "142a_raw_count_matrix_structure.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 21. BUILD AUDIT DECISION TABLE
# ==============================================================================

find_availability <- function(pattern) {
  
  matched <- metadata_variable_audit %>%
    
    dplyr::filter(
      
      grepl(
        pattern,
        variable,
        ignore.case = TRUE
      )
    )
  
  
  if (
    nrow(
      matched
    ) == 0
  ) {
    
    return(
      list(
        available = FALSE,
        variables = "",
        max_n = 0L
      )
    )
  }
  
  
  return(
    list(
      
      available =
        TRUE,
      
      variables =
        paste(
          matched$variable,
          collapse = "; "
        ),
      
      max_n =
        max(
          matched$n_nonmissing,
          na.rm = TRUE
        )
    )
  )
}


disease_avail <- find_availability(
  "disease.*state|diagnos|sepsis"
)


sofa_avail <- find_availability(
  "sofa"
)


mortality_avail <- find_availability(
  "mortality|death|surviv"
)


age_avail <- find_availability(
  "^age$|age"
)


sex_avail <- find_availability(
  "^sex$|gender"
)


site_avail <- find_availability(
  "collection.*site|site"
)


location_avail <- find_availability(
  "collection.*location|location"
)


severity_avail <- find_availability(
  "severity|severe"
)


endotype_avail <- find_availability(
  "endotype|subtype|cluster"
)


audit_decision_table <- tibble::tibble(
  
  domain = c(
    "Disease state / diagnosis",
    "SOFA",
    "Mortality",
    "Age",
    "Sex",
    "Collection site",
    "Collection location",
    "Severity field",
    "Explicit endotype field",
    "Published endotype labels in metadata",
    "Raw-count sample mapping"
  ),
  
  available = c(
    disease_avail$available,
    sofa_avail$available,
    mortality_avail$available,
    age_avail$available,
    sex_avail$available,
    site_avail$available,
    location_avail$available,
    severity_avail$available,
    endotype_avail$available,
    nrow(
      published_endotype_metadata_hits
    ) > 0,
    nrow(
      mapped_count_columns
    ) > 0
  ),
  
  available_n_or_entries = c(
    disease_avail$max_n,
    sofa_avail$max_n,
    mortality_avail$max_n,
    age_avail$max_n,
    sex_avail$max_n,
    site_avail$max_n,
    location_avail$max_n,
    severity_avail$max_n,
    endotype_avail$max_n,
    nrow(
      published_endotype_metadata_hits
    ),
    nrow(
      mapped_count_columns
    )
  ),
  
  variables_or_note = c(
    disease_avail$variables,
    sofa_avail$variables,
    mortality_avail$variables,
    age_avail$variables,
    sex_avail$variables,
    site_avail$variables,
    location_avail$variables,
    severity_avail$variables,
    endotype_avail$variables,
    paste(
      unique(
        published_endotype_metadata_hits$published_endotype_label
      ),
      collapse = "; "
    ),
    paste0(
      nrow(
        mapped_count_columns
      ),
      " GEO sample columns mapped"
    )
  )
)


cat(
  "AUDIT DECISION TABLE:\n"
)


print(
  audit_decision_table,
  n = Inf,
  width = Inf
)


cat("\n")


write.csv(
  
  audit_decision_table,
  
  file.path(
    tables_dir,
    "142a_audit_decision_table.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 22. EXCEL WORKBOOK
# ==============================================================================

run_info <- tibble::tibble(
  
  parameter = c(
    "script",
    "run_date",
    "GEO_accession",
    "purpose",
    "expression_values_analyzed",
    "five_gene_panel_tested",
    "score_calculated",
    "endpoint_selected",
    "expected_GEO_samples",
    "raw_count_file"
  ),
  
  value = c(
    script_name,
    as.character(
      run_date
    ),
    gse_id,
    "Metadata / cohort / raw-count-header audit before external validation",
    "NO",
    "NO",
    "NO",
    "NO",
    "392",
    raw_counts_filename
  )
)


workbook_file <- file.path(
  tables_dir,
  "142a_GSE185263_metadata_audit.xlsx"
)


wb <- openxlsx::createWorkbook()


sheet_list <- list(
  
  "00_run_info" =
    run_info,
  
  "01_audit_decision" =
    audit_decision_table,
  
  "02_priority_variables" =
    priority_variable_matches,
  
  "03_candidate_variables" =
    candidate_metadata_variables,
  
  "04_variable_distributions" =
    candidate_variable_distributions,
  
  "05_numeric_summary" =
    numeric_variable_summary,
  
  "06_characteristic_keys" =
    characteristic_key_summary,
  
  "07_metadata_parsed" =
    metadata_parsed,
  
  "08_metadata_raw" =
    pheno_raw,
  
  "09_sample_prefixes" =
    sample_prefix_summary,
  
  "10_timepoint_suffixes" =
    timepoint_summary,
  
  "11_repeated_participants" =
    repeated_participant_guesses,
  
  "12_repeated_samples" =
    repeated_sample_details,
  
  "13_endotype_metadata_hits" =
    published_endotype_metadata_hits,
  
  "14_count_structure" =
    count_matrix_structure,
  
  "15_count_column_mapping" =
    count_column_mapping,
  
  "16_missing_count_samples" =
    geo_samples_missing_from_counts
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
# 23. AUTOMATED TEXT AUDIT SUMMARY
# ==============================================================================

audit_summary <- c(
  
  "GSE185263 METADATA AUDIT",
  
  "====================================================================",
  
  "",
  
  paste0(
    "Audit date: ",
    run_date
  ),
  
  "",
  
  paste0(
    "GEO samples: ",
    nrow(
      metadata_core
    )
  ),
  
  paste0(
    "Structured metadata variables: ",
    ncol(
      metadata_parsed
    )
  ),
  
  paste0(
    "Parsed GEO characteristic keys: ",
    nrow(
      characteristic_key_summary
    )
  ),
  
  paste0(
    "Count-matrix GEO sample columns mapped: ",
    nrow(
      mapped_count_columns
    )
  ),
  
  paste0(
    "GEO samples missing from raw-count matrix: ",
    nrow(
      geo_samples_missing_from_counts
    )
  ),
  
  "",
  
  "SAMPLE PREFIXES:",
  
  paste(
    paste0(
      sample_prefix_summary$sample_prefix,
      "=",
      sample_prefix_summary$n_samples
    ),
    collapse = "; "
  ),
  
  "",
  
  "TIMEPOINT SUFFIXES:",
  
  paste(
    paste0(
      timepoint_summary$timepoint_suffix,
      "=",
      timepoint_summary$n_samples
    ),
    collapse = "; "
  ),
  
  "",
  
  paste0(
    "Possible repeated/time-point participant IDs: ",
    nrow(
      repeated_participant_guesses
    )
  ),
  
  "",
  
  "PRIORITY CLINICAL VARIABLES:",
  
  paste(
    apply(
      priority_variable_matches,
      1,
      function(row_values) {
        
        paste0(
          row_values[["concept"]],
          ": ",
          ifelse(
            is.na(
              row_values[["variable"]]
            ),
            "NOT FOUND",
            paste0(
              row_values[["variable"]],
              " (n=",
              row_values[["n_nonmissing"]],
              ")"
            )
          )
        )
      }
    ),
    collapse = "\n"
  ),
  
  "",
  
  paste0(
    "Published NPS/INF/IHD/IFN/ADA label entries found directly ",
    "in GEO metadata: ",
    nrow(
      published_endotype_metadata_hits
    )
  ),
  
  "",
  
  "IMPORTANT:",
  
  "- No expression values were analyzed.",
  
  "- The five frozen genes were not queried.",
  
  "- No five-gene score was calculated.",
  
  "- No validation endpoint was selected from expression results.",
  
  "",
  
  paste0(
    "The next step is to review these metadata audit results and ",
    "freeze the GSE185263 validation design before Script 142b."
  )
)


audit_summary_file <- file.path(
  text_dir,
  "142a_GSE185263_metadata_audit_summary.txt"
)


writeLines(
  
  audit_summary,
  
  con =
    audit_summary_file,
  
  useBytes =
    TRUE
)


# ==============================================================================
# 24. INPUT / DOWNLOAD MANIFEST
# ==============================================================================

manifest <- tibble::tibble(
  
  item = c(
    "GSE185263_raw_counts"
  ),
  
  path = c(
    raw_counts_file
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
    "142a_input_download_manifest.csv"
  ),
  
  row.names =
    FALSE
)


# ==============================================================================
# 25. SESSION INFO
# ==============================================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      logs_dir,
      "142a_sessionInfo.txt"
    )
)


# ==============================================================================
# 26. FINAL REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 142a completed successfully.\n")
cat("====================================================================\n\n")


cat(
  "GSE185263 GEO SAMPLES:\n"
)


cat(
  nrow(
    metadata_core
  ),
  "\n\n"
)


cat(
  "SAMPLE PREFIX DISTRIBUTION:\n"
)


print(
  sample_prefix_summary,
  n = Inf
)


cat("\n")


cat(
  "TIMEPOINT SUFFIX DISTRIBUTION:\n"
)


print(
  timepoint_summary,
  n = Inf
)


cat("\n")


cat(
  "POSSIBLE REPEATED/TIMEPOINT PARTICIPANTS:\n"
)


cat(
  nrow(
    repeated_participant_guesses
  ),
  "\n\n"
)


cat(
  "PRIORITY CLINICAL VARIABLE AVAILABILITY:\n"
)


print(
  priority_variable_matches,
  n = Inf,
  width = Inf
)


cat("\n")


cat(
  "PUBLISHED ENDOTYPE LABELS DIRECTLY FOUND IN METADATA:\n"
)


if (
  nrow(
    published_endotype_metadata_hits
  ) == 0
) {
  
  cat(
    "None detected directly in GEO metadata.\n\n"
  )
  
} else {
  
  print(
    published_endotype_metadata_hits,
    n = Inf,
    width = Inf
  )
  
  
  cat("\n")
}


cat(
  "RAW COUNT MATRIX MAPPING:\n"
)


cat(
  "Total count-file columns = ",
  length(
    count_column_names
  ),
  "\n",
  sep = ""
)


cat(
  "Mapped GEO sample columns = ",
  nrow(
    mapped_count_columns
  ),
  "\n",
  sep = ""
)


cat(
  "Unmapped/non-sample columns = ",
  nrow(
    unmapped_count_columns
  ),
  "\n",
  sep = ""
)


cat(
  "GEO samples missing from counts = ",
  nrow(
    geo_samples_missing_from_counts
  ),
  "\n\n",
  sep = ""
)


cat(
  "AUDIT DECISION TABLE:\n"
)


print(
  audit_decision_table,
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
  "Audit summary:\n"
)


cat(
  normalizePath(
    audit_summary_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


cat(
  "IMPORTANT:\n"
)


cat(
  "- No five-gene expression analysis performed.\n"
)


cat(
  "- No five-gene score calculated.\n"
)


cat(
  "- No AUC calculated.\n"
)


cat(
  "- No validation endpoint selected from expression results.\n"
)


cat(
  "- Script 142b must be designed only AFTER reviewing this metadata audit.\n\n"
)


cat(
  "Done.\n"
)