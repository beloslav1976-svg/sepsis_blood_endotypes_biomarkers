# ==============================================================================
# Script 135
# Blood-only endotype-informed biomarker validation
#
# Project: Sepsis_DESeq2
#
# Manuscript:
# Blood Transcriptomic Endotypes and Endotype-Informed Candidate Biomarkers
# for Molecular Stratification of Sepsis
#
# Version: corrected CTS source = BP-only
# Date: 2026-08-17
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

script_id <- "135"
script_name <- "135_blood_endotype_biomarker_validation.R"
run_date <- Sys.time()

cat("\n")
cat("====================================================================\n")
cat("Running Script 135\n")
cat("Blood-only endotype-informed biomarker validation\n")
cat("====================================================================\n\n")

cat("Project directory:\n")
cat(normalizePath(getwd(), winslash = "/", mustWork = FALSE), "\n\n")

cat("Run date:\n")
cat(as.character(run_date), "\n\n")


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_packages <- c(
  "edgeR",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readxl",
  "openxlsx",
  "pROC",
  "glmnet",
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
      "Missing packages:\n",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(pROC)
  library(glmnet)
  library(tibble)
})

cat("Required packages loaded successfully.\n\n")


# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

clean_sample_id <- function(x) {
  
  x <- toupper(as.character(x))
  x <- trimws(x)
  x <- gsub("[^A-Z0-9]", "", x)
  
  return(x)
}


normalize_srs_label <- function(x) {
  
  x <- toupper(trimws(as.character(x)))
  
  x[x == "1"] <- "SRS1"
  x[x == "2"] <- "SRS2"
  x[x == "3"] <- "SRS3"
  
  return(x)
}


normalize_cts_label <- function(x) {
  
  x <- toupper(trimws(as.character(x)))
  
  x[x == "1"] <- "CTS1"
  x[x == "2"] <- "CTS2"
  x[x == "3"] <- "CTS3"
  
  x[x == "CTS 1"] <- "CTS1"
  x[x == "CTS 2"] <- "CTS2"
  x[x == "CTS 3"] <- "CTS3"
  
  return(x)
}


find_existing_file <- function(
    preferred_paths,
    fallback_pattern = NULL,
    search_root = "."
) {
  
  for (p in preferred_paths) {
    
    if (file.exists(p)) {
      return(p)
    }
  }
  
  if (!is.null(fallback_pattern)) {
    
    candidates <- list.files(
      path = search_root,
      pattern = fallback_pattern,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    candidates <- candidates[file.exists(candidates)]
    
    if (length(candidates) == 1) {
      return(candidates[1])
    }
    
    if (length(candidates) > 1) {
      
      cat("\nMultiple candidate files found:\n")
      print(candidates)
      
      stop(
        paste0(
          "More than one file matched: ",
          fallback_pattern
        )
      )
    }
  }
  
  return(NA_character_)
}


find_column <- function(
    df,
    candidates,
    required = TRUE
) {
  
  nm <- names(df)
  
  nm_lower <- tolower(trimws(nm))
  candidates_lower <- tolower(trimws(candidates))
  
  exact_hits <- which(
    nm_lower %in% candidates_lower
  )
  
  if (length(exact_hits) > 0) {
    return(nm[exact_hits[1]])
  }
  
  hits <- integer(0)
  
  for (candidate in candidates_lower) {
    
    h <- grep(
      candidate,
      nm_lower,
      fixed = TRUE
    )
    
    hits <- unique(
      c(
        hits,
        h
      )
    )
  }
  
  if (length(hits) == 1) {
    return(nm[hits])
  }
  
  if (required) {
    
    stop(
      paste0(
        "Required column could not be identified.\n",
        "Candidates: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns: ",
        paste(nm, collapse = ", ")
      )
    )
  }
  
  return(NA_character_)
}


find_excel_sheet_with_columns <- function(
    file,
    sample_candidates,
    class_candidates,
    score_candidates = NULL
) {
  
  sheets <- readxl::excel_sheets(file)
  
  for (sheet_name in sheets) {
    
    tmp <- suppressMessages(
      readxl::read_excel(
        file,
        sheet = sheet_name
      )
    )
    
    sample_col <- tryCatch(
      find_column(
        tmp,
        sample_candidates,
        required = TRUE
      ),
      error = function(e) NA_character_
    )
    
    class_col <- tryCatch(
      find_column(
        tmp,
        class_candidates,
        required = TRUE
      ),
      error = function(e) NA_character_
    )
    
    score_col <- NA_character_
    
    if (!is.null(score_candidates)) {
      
      score_col <- tryCatch(
        find_column(
          tmp,
          score_candidates,
          required = TRUE
        ),
        error = function(e) NA_character_
      )
    }
    
    base_ok <- !is.na(sample_col) &&
      !is.na(class_col)
    
    score_ok <- TRUE
    
    if (!is.null(score_candidates)) {
      score_ok <- !is.na(score_col)
    }
    
    if (base_ok && score_ok) {
      
      return(
        list(
          data = tmp,
          sheet = sheet_name,
          sample_col = sample_col,
          class_col = class_col,
          score_col = score_col
        )
      )
    }
  }
  
  stop(
    paste0(
      "Could not identify suitable worksheet in:\n",
      file
    )
  )
}


safe_z <- function(x) {
  
  s <- sd(
    x,
    na.rm = TRUE
  )
  
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  
  return(
    (
      x - mean(
        x,
        na.rm = TRUE
      )
    ) / s
  )
}


row_zscore <- function(mat) {
  
  z <- t(
    apply(
      mat,
      1,
      safe_z
    )
  )
  
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  
  return(z)
}


signed_panel_score <- function(
    z_matrix,
    up_genes,
    down_genes
) {
  
  up_present <- intersect(
    up_genes,
    rownames(z_matrix)
  )
  
  down_present <- intersect(
    down_genes,
    rownames(z_matrix)
  )
  
  if (length(up_present) == 0) {
    stop("No UP genes available.")
  }
  
  if (length(down_present) == 0) {
    stop("No DOWN genes available.")
  }
  
  up_score <- colMeans(
    z_matrix[
      up_present,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  down_score <- colMeans(
    z_matrix[
      down_present,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  return(
    up_score - down_score
  )
}


calc_auc <- function(
    truth,
    score
) {
  
  keep <- complete.cases(
    truth,
    score
  )
  
  keep <- keep &
    is.finite(score)
  
  truth2 <- truth[keep]
  score2 <- score[keep]
  
  if (length(unique(truth2)) != 2) {
    return(NA_real_)
  }
  
  roc_obj <- suppressWarnings(
    pROC::roc(
      response = truth2,
      predictor = score2,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
  )
  
  return(
    as.numeric(
      pROC::auc(
        roc_obj
      )
    )
  )
}


calc_auc_with_ci <- function(
    truth,
    score
) {
  
  keep <- complete.cases(
    truth,
    score
  )
  
  keep <- keep &
    is.finite(score)
  
  truth2 <- truth[keep]
  score2 <- score[keep]
  
  if (length(unique(truth2)) != 2) {
    
    return(
      tibble(
        auc = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    )
  }
  
  roc_obj <- suppressWarnings(
    pROC::roc(
      response = truth2,
      predictor = score2,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
  )
  
  auc_value <- as.numeric(
    pROC::auc(
      roc_obj
    )
  )
  
  ci_value <- suppressWarnings(
    pROC::ci.auc(
      roc_obj,
      method = "delong"
    )
  )
  
  return(
    tibble(
      auc = auc_value,
      ci_low = as.numeric(ci_value[1]),
      ci_high = as.numeric(ci_value[3])
    )
  )
}


make_stratified_foldid <- function(
    y,
    k = 5
) {
  
  foldid <- integer(
    length(y)
  )
  
  for (cl in sort(unique(y))) {
    
    idx <- which(
      y == cl
    )
    
    assignments <- rep(
      seq_len(k),
      length.out = length(idx)
    )
    
    assignments <- sample(
      assignments,
      size = length(assignments),
      replace = FALSE
    )
    
    foldid[idx] <- assignments
  }
  
  return(foldid)
}


standardize_train_test <- function(
    train_expr,
    test_expr,
    genes
) {
  
  genes_present <- intersect(
    genes,
    rownames(train_expr)
  )
  
  if (length(genes_present) == 0) {
    stop("No requested genes present.")
  }
  
  train_x <- train_expr[
    genes_present,
    ,
    drop = FALSE
  ]
  
  test_x <- test_expr[
    genes_present,
    ,
    drop = FALSE
  ]
  
  train_mean <- rowMeans(
    train_x,
    na.rm = TRUE
  )
  
  train_sd <- apply(
    train_x,
    1,
    sd,
    na.rm = TRUE
  )
  
  train_sd[
    !is.finite(train_sd)
  ] <- 1
  
  train_sd[
    train_sd == 0
  ] <- 1
  
  train_z <- sweep(
    train_x,
    1,
    train_mean,
    "-"
  )
  
  train_z <- sweep(
    train_z,
    1,
    train_sd,
    "/"
  )
  
  test_z <- sweep(
    test_x,
    1,
    train_mean,
    "-"
  )
  
  test_z <- sweep(
    test_z,
    1,
    train_sd,
    "/"
  )
  
  return(
    list(
      train = train_z,
      test = test_z
    )
  )
}


ridge_predict_outer_fold <- function(
    train_z,
    test_z,
    y_train
) {
  
  x_train <- t(train_z)
  x_test <- t(test_z)
  
  if (length(unique(y_train)) != 2) {
    return(rep(NA_real_, nrow(x_test)))
  }
  
  k_inner <- min(
    5,
    sum(y_train == 0),
    sum(y_train == 1)
  )
  
  if (k_inner < 3) {
    return(rep(NA_real_, nrow(x_test)))
  }
  
  inner_foldid <- make_stratified_foldid(
    y = y_train,
    k = k_inner
  )
  
  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      alpha = 0,
      foldid = inner_foldid,
      type.measure = "deviance",
      standardize = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(rep(NA_real_, nrow(x_test)))
  }
  
  prediction <- predict(
    fit,
    newx = x_test,
    s = "lambda.1se",
    type = "response"
  )
  
  return(
    as.numeric(
      prediction
    )
  )
}


run_repeated_cv <- function(
    expr,
    truth,
    model_name,
    genes,
    model_type,
    up_genes = NULL,
    down_genes = NULL,
    n_repeats = 100,
    n_folds = 5,
    seed = 20260817
) {
  
  if (!model_type %in% c(
    "signed_score",
    "ridge"
  )) {
    stop("Unknown model_type.")
  }
  
  set.seed(seed)
  
  sample_ids <- colnames(expr)
  
  auc_list <- vector(
    "list",
    n_repeats
  )
  
  prediction_list <- vector(
    "list",
    n_repeats
  )
  
  for (r in seq_len(n_repeats)) {
    
    foldid <- make_stratified_foldid(
      y = truth,
      k = n_folds
    )
    
    prediction <- rep(
      NA_real_,
      length(truth)
    )
    
    for (fold_id in seq_len(n_folds)) {
      
      test_idx <- which(
        foldid == fold_id
      )
      
      train_idx <- which(
        foldid != fold_id
      )
      
      train_expr <- expr[
        ,
        train_idx,
        drop = FALSE
      ]
      
      test_expr <- expr[
        ,
        test_idx,
        drop = FALSE
      ]
      
      standardized <- standardize_train_test(
        train_expr = train_expr,
        test_expr = test_expr,
        genes = genes
      )
      
      if (model_type == "signed_score") {
        
        prediction[test_idx] <- signed_panel_score(
          z_matrix = standardized$test,
          up_genes = up_genes,
          down_genes = down_genes
        )
        
      } else {
        
        prediction[test_idx] <- ridge_predict_outer_fold(
          train_z = standardized$train,
          test_z = standardized$test,
          y_train = truth[train_idx]
        )
      }
    }
    
    auc_value <- calc_auc(
      truth = truth,
      score = prediction
    )
    
    auc_list[[r]] <- tibble(
      model = model_name,
      repeat_id = r,
      auc = auc_value
    )
    
    prediction_list[[r]] <- tibble(
      model = model_name,
      repeat_id = r,
      sample_id = sample_ids,
      truth = truth,
      prediction = prediction
    )
  }
  
  auc_table <- bind_rows(
    auc_list
  )
  
  prediction_table <- bind_rows(
    prediction_list
  )
  
  sample_prediction_table <- prediction_table %>%
    group_by(
      model,
      sample_id,
      truth
    ) %>%
    summarise(
      mean_prediction = mean(
        prediction,
        na.rm = TRUE
      ),
      median_prediction = median(
        prediction,
        na.rm = TRUE
      ),
      n_predictions = sum(
        is.finite(
          prediction
        )
      ),
      .groups = "drop"
    )
  
  pooled_auc <- calc_auc(
    truth = sample_prediction_table$truth,
    score = sample_prediction_table$mean_prediction
  )
  
  auc_summary <- tibble(
    model = model_name,
    
    n_repeats = sum(
      is.finite(
        auc_table$auc
      )
    ),
    
    median_repeat_auc = median(
      auc_table$auc,
      na.rm = TRUE
    ),
    
    mean_repeat_auc = mean(
      auc_table$auc,
      na.rm = TRUE
    ),
    
    repeat_auc_q025 = as.numeric(
      quantile(
        auc_table$auc,
        probs = 0.025,
        na.rm = TRUE
      )
    ),
    
    repeat_auc_q975 = as.numeric(
      quantile(
        auc_table$auc,
        probs = 0.975,
        na.rm = TRUE
      )
    ),
    
    pooled_sample_level_auc = pooled_auc
  )
  
  return(
    list(
      auc_table = auc_table,
      auc_summary = auc_summary,
      predictions = prediction_table,
      sample_predictions = sample_prediction_table
    )
  )
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
  group2 <- factor(
    group[keep]
  )
  
  if (nlevels(group2) < 2) {
    return(NA_real_)
  }
  
  kw <- kruskal.test(
    x2 ~ group2
  )
  
  H <- as.numeric(
    kw$statistic
  )
  
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


# ==============================================================================
# 3. INPUT FILES
# ==============================================================================

counts_file <- find_existing_file(
  preferred_paths = c(
    "data/counts_all.csv",
    "data/processed/counts_all.csv"
  ),
  fallback_pattern = "^counts_all\\.csv$",
  search_root = "."
)


srs_file <- find_existing_file(
  preferred_paths = c(
    file.path(
      "results",
      "sepstratifier",
      "publication_summary",
      "43_SRS_by_each_sample.xlsx"
    )
  ),
  fallback_pattern = "43_SRS_by_each_sample\\.xlsx$",
  search_root = "results"
)


# IMPORTANT CORRECTION:
# Use BP-only formal CTS classification.
#
# DO NOT use:
# 107_CTS_predictions_BP_BC.xlsx
#
# because that file contains a different 18/5/12 classification.
#
# The final manuscript CTS classification is from BP-only:
# CTS1=14, CTS2=6, CTS3=15.

cts_file <- find_existing_file(
  preferred_paths = c(
    file.path(
      "results",
      "cts_consensus",
      "tables",
      "107_CTS_predictions_BP_only.xlsx"
    )
  ),
  fallback_pattern = "107_CTS_predictions_BP_only\\.xlsx$",
  search_root = "results"
)


input_check <- tibble(
  input = c(
    "counts_file",
    "FINAL_SRS_file",
    "FINAL_CTS_BP_only_file"
  ),
  path = c(
    counts_file,
    srs_file,
    cts_file
  )
)


input_check$exists <- !is.na(
  input_check$path
) &
  file.exists(
    input_check$path
  )


cat("Input file check:\n")
print(input_check)
cat("\n")


if (any(!input_check$exists)) {
  
  stop(
    paste0(
      "Required input file missing:\n",
      paste(
        input_check$input[
          !input_check$exists
        ],
        collapse = ", "
      )
    )
  )
}


cat("All required input files identified successfully.\n\n")


# ==============================================================================
# 4. OUTPUT DIRECTORIES
# ==============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "135_validation"
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
# 5. INPUT MANIFEST
# ==============================================================================

input_manifest <- input_check

fi <- file.info(
  input_manifest$path
)

input_manifest$file_size_bytes <- fi$size
input_manifest$modified_time <- as.character(fi$mtime)

input_manifest$md5 <- unname(
  tools::md5sum(
    input_manifest$path
  )
)


write.csv(
  input_manifest,
  file.path(
    logs_dir,
    "135_input_file_manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 6. PRE-SPECIFIED PANELS
# ==============================================================================

primary_up <- c(
  "CD177",
  "HK3",
  "IRAK3"
)

primary_down <- c(
  "CARD11",
  "IKZF2"
)

primary_genes <- c(
  primary_up,
  primary_down
)


dcaf_up <- c(
  "CD177",
  "HK3",
  "IRAK3"
)

dcaf_down <- c(
  "CARD11",
  "DCAF17"
)

dcaf_genes <- c(
  dcaf_up,
  dcaf_down
)


septicyte_genes <- c(
  "CEACAM4",
  "LAMP1",
  "PLA2G7",
  "PLAC8"
)


all_candidate_genes <- unique(
  c(
    primary_genes,
    dcaf_genes,
    septicyte_genes
  )
)


panel_definition <- bind_rows(
  
  tibble(
    panel = "Primary_5_gene",
    gene = primary_up,
    expected_direction = "UP_in_sepsis",
    biological_arm = "myeloid_neutrophil_activation"
  ),
  
  tibble(
    panel = "Primary_5_gene",
    gene = primary_down,
    expected_direction = "DOWN_in_sepsis",
    biological_arm = "adaptive_T_cell_program"
  ),
  
  tibble(
    panel = "DCAF17_5_gene",
    gene = dcaf_up,
    expected_direction = "UP_in_sepsis",
    biological_arm = "myeloid_neutrophil_activation"
  ),
  
  tibble(
    panel = "DCAF17_5_gene",
    gene = dcaf_down,
    expected_direction = "DOWN_in_sepsis",
    biological_arm = "supporting_down_arm"
  ),
  
  tibble(
    panel = "SeptiCyte_related_4_gene",
    gene = septicyte_genes,
    expected_direction = "not_imposed",
    biological_arm = "published_external_benchmark"
  )
)


cat("Pre-specified panels:\n")
print(panel_definition)
cat("\n")


# ==============================================================================
# 7. READ BLOOD COUNTS
# ==============================================================================

cat("Reading counts:\n")
cat(counts_file, "\n\n")


counts_df <- read.csv(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


gene_col <- find_column(
  counts_df,
  candidates = c(
    "Gene",
    "gene",
    "GeneSymbol",
    "gene_symbol",
    "symbol"
  ),
  required = TRUE
)


cat("Detected gene column:\n")
cat(gene_col, "\n\n")


clean_column_names <- clean_sample_id(
  names(counts_df)
)


blood_column_flag <- grepl(
  "^(BP|BC)[0-9]+$",
  clean_column_names
)


blood_columns <- names(counts_df)[
  blood_column_flag
]


blood_sample_ids <- clean_column_names[
  blood_column_flag
]


if (length(blood_columns) != 45) {
  
  cat("Detected blood columns:\n")
  print(blood_columns)
  
  stop(
    paste0(
      "Expected 45 blood samples, detected ",
      length(blood_columns)
    )
  )
}


gene_symbols <- as.character(
  counts_df[[gene_col]]
)


valid_gene <- !is.na(gene_symbols) &
  trimws(gene_symbols) != ""


counts_sub <- counts_df[
  valid_gene,
  blood_columns,
  drop = FALSE
]


gene_symbols <- gene_symbols[
  valid_gene
]


blood_counts <- data.matrix(
  counts_sub
)


rownames(blood_counts) <- gene_symbols
colnames(blood_counts) <- blood_sample_ids


if (anyDuplicated(colnames(blood_counts)) > 0) {
  stop("Duplicated blood sample IDs detected.")
}


if (anyDuplicated(rownames(blood_counts)) > 0) {
  
  cat(
    "Duplicated gene symbols detected. ",
    "Collapsing by summed raw counts.\n\n"
  )
  
  blood_counts <- rowsum(
    blood_counts,
    group = rownames(blood_counts),
    reorder = FALSE
  )
}


ercc_flag <- grepl(
  "^ERCC[-_]",
  rownames(blood_counts),
  ignore.case = TRUE
)


if (any(ercc_flag)) {
  
  cat(
    "Removing ",
    sum(ercc_flag),
    " technical ERCC spike-in rows.\n\n",
    sep = ""
  )
  
  blood_counts <- blood_counts[
    !ercc_flag,
    ,
    drop = FALSE
  ]
}


if (any(!is.finite(blood_counts))) {
  stop("Non-finite values detected in raw counts.")
}


if (any(blood_counts < 0)) {
  stop("Negative raw counts detected.")
}


blood_group <- ifelse(
  grepl(
    "^BP",
    colnames(blood_counts)
  ),
  "BP",
  "BC"
)


blood_group <- factor(
  blood_group,
  levels = c(
    "BC",
    "BP"
  )
)


sample_metadata <- tibble(
  sample_id = colnames(blood_counts),
  condition = as.character(blood_group),
  truth_sepsis = ifelse(
    blood_group == "BP",
    1L,
    0L
  )
)


cat("Blood sample distribution:\n")
print(table(blood_group))
cat("\n")


if (sum(blood_group == "BP") != 35) {
  stop("Expected BP=35.")
}


if (sum(blood_group == "BC") != 10) {
  stop("Expected BC=10.")
}


cat(
  "Blood dataset validation PASSED: BP=35, BC=10.\n\n"
)


# ==============================================================================
# 8. GENE AVAILABILITY
# ==============================================================================

gene_availability <- tibble(
  gene = all_candidate_genes,
  detected = all_candidate_genes %in%
    rownames(blood_counts)
)


cat("Candidate gene availability:\n")
print(gene_availability)
cat("\n")


missing_primary <- setdiff(
  primary_genes,
  rownames(blood_counts)
)


if (length(missing_primary) > 0) {
  
  stop(
    paste0(
      "Primary genes missing: ",
      paste(
        missing_primary,
        collapse = ", "
      )
    )
  )
}


missing_dcaf <- setdiff(
  dcaf_genes,
  rownames(blood_counts)
)


if (length(missing_dcaf) > 0) {
  
  stop(
    paste0(
      "DCAF17 panel genes missing: ",
      paste(
        missing_dcaf,
        collapse = ", "
      )
    )
  )
}


missing_septicyte <- setdiff(
  septicyte_genes,
  rownames(blood_counts)
)


if (length(missing_septicyte) > 0) {
  
  warning(
    paste0(
      "SeptiCyte-related genes missing: ",
      paste(
        missing_septicyte,
        collapse = ", "
      )
    )
  )
}


# ==============================================================================
# 9. TMM NORMALIZATION
# ==============================================================================

dge <- edgeR::DGEList(
  counts = blood_counts,
  group = blood_group
)


# Current edgeR name.
dge <- edgeR::normLibSizes(
  dge,
  method = "TMM"
)


logcpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 1
)


cat(
  "TMM-adjusted logCPM created: ",
  nrow(logcpm),
  " genes x ",
  ncol(logcpm),
  " samples.\n\n",
  sep = ""
)


# ==============================================================================
# 10. PANEL SCORES
# ==============================================================================

selected_genes_detected <- intersect(
  all_candidate_genes,
  rownames(logcpm)
)


panel_expr <- logcpm[
  selected_genes_detected,
  ,
  drop = FALSE
]


panel_z <- row_zscore(
  panel_expr
)


primary_score <- signed_panel_score(
  z_matrix = panel_z,
  up_genes = primary_up,
  down_genes = primary_down
)


dcaf_score <- signed_panel_score(
  z_matrix = panel_z,
  up_genes = dcaf_up,
  down_genes = dcaf_down
)


dcaf17_single_score <- -as.numeric(
  panel_z[
    "DCAF17",
    ,
    drop = TRUE
  ]
)


myeloid_up_score <- colMeans(
  panel_z[
    primary_up,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)


adaptive_suppression_score <- -colMeans(
  panel_z[
    primary_down,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)


score_table <- sample_metadata %>%
  mutate(
    primary_5gene_score = primary_score,
    myeloid_UP_score = myeloid_up_score,
    adaptive_suppression_score = adaptive_suppression_score,
    DCAF17_5gene_score = dcaf_score,
    DCAF17_single_score = dcaf17_single_score
  )


cat("Descriptive panel scores calculated.\n\n")


# ==============================================================================
# 11. READ FINAL SRS
# ==============================================================================

cat("Reading FINAL SRS:\n")
cat(srs_file, "\n\n")


srs_sheet <- find_excel_sheet_with_columns(
  file = srs_file,
  sample_candidates = c(
    "sample_id",
    "sample",
    "barcode_id",
    "barcode"
  ),
  class_candidates = c(
    "SRS",
    "SRS_class",
    "srs_class",
    "classification"
  ),
  score_candidates = c(
    "SRSq",
    "srsq",
    "SRS_q"
  )
)


srs_raw <- srs_sheet$data
srs_sample_col <- srs_sheet$sample_col
srs_class_col <- srs_sheet$class_col
srsq_col <- srs_sheet$score_col


cat("Selected SRS worksheet:\n")
cat(srs_sheet$sheet, "\n\n")


mnn_col <- find_column(
  srs_raw,
  candidates = c(
    "mNN_outlier",
    "mnn_outlier",
    "outlier"
  ),
  required = FALSE
)


srs_table <- tibble(
  sample_id = clean_sample_id(
    srs_raw[[srs_sample_col]]
  ),
  SRS = normalize_srs_label(
    srs_raw[[srs_class_col]]
  ),
  SRSq = suppressWarnings(
    as.numeric(
      srs_raw[[srsq_col]]
    )
  )
)


if (!is.na(mnn_col)) {
  
  srs_table$mNN_outlier_from_file <- as.character(
    srs_raw[[mnn_col]]
  )
  
} else {
  
  srs_table$mNN_outlier_from_file <- NA_character_
}


srs_table <- srs_table %>%
  filter(
    grepl(
      "^(BP|BC)[0-9]+$",
      sample_id
    )
  ) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )


cat("SRS distribution:\n")

print(
  srs_table %>%
    mutate(
      condition = ifelse(
        grepl("^BP", sample_id),
        "BP",
        "BC"
      )
    ) %>%
    count(
      condition,
      SRS
    )
)

cat("\n")


# ==============================================================================
# 12. HARD CHECK FINAL SRS
# ==============================================================================

bp_srs <- table(
  factor(
    srs_table$SRS[
      grepl(
        "^BP",
        srs_table$sample_id
      )
    ],
    levels = c(
      "SRS1",
      "SRS2",
      "SRS3"
    )
  )
)


bc_srs <- table(
  factor(
    srs_table$SRS[
      grepl(
        "^BC",
        srs_table$sample_id
      )
    ],
    levels = c(
      "SRS1",
      "SRS2",
      "SRS3"
    )
  )
)


observed_bp_srs <- setNames(
  as.integer(bp_srs),
  c(
    "SRS1",
    "SRS2",
    "SRS3"
  )
)


observed_bc_srs <- setNames(
  as.integer(bc_srs),
  c(
    "SRS1",
    "SRS2",
    "SRS3"
  )
)


expected_bp_srs <- c(
  SRS1 = 28L,
  SRS2 = 7L,
  SRS3 = 0L
)


expected_bc_srs <- c(
  SRS1 = 0L,
  SRS2 = 3L,
  SRS3 = 7L
)


if (!identical(
  observed_bp_srs,
  expected_bp_srs
)) {
  
  print(observed_bp_srs)
  
  stop(
    "FINAL SRS validation FAILED for BP."
  )
}


if (!identical(
  observed_bc_srs,
  expected_bc_srs
)) {
  
  print(observed_bc_srs)
  
  stop(
    "FINAL SRS validation FAILED for BC."
  )
}


cat("FINAL SRS validation PASSED.\n")
cat("BP: SRS1=28, SRS2=7, SRS3=0.\n")
cat("BC: SRS1=0, SRS2=3, SRS3=7.\n\n")


# ==============================================================================
# 13. mNN OUTLIERS
# ==============================================================================

known_mnn_outliers <- c(
  "BP27",
  "BP31",
  "BP26",
  "BP10"
)


srs_table <- srs_table %>%
  mutate(
    mNN_outlier_final =
      sample_id %in%
      known_mnn_outliers
  )


cat("Predefined mNN outliers:\n")

print(
  srs_table %>%
    filter(
      mNN_outlier_final
    ) %>%
    select(
      sample_id,
      SRS,
      SRSq
    )
)

cat("\n")


# ==============================================================================
# 14. READ FINAL CTS — BP ONLY
# ==============================================================================

cat("Reading FINAL BP-only CTS file:\n")
cat(cts_file, "\n\n")


cts_sheet <- find_excel_sheet_with_columns(
  file = cts_file,
  sample_candidates = c(
    "sample_id",
    "sample",
    "Sample",
    "barcode_id"
  ),
  class_candidates = c(
    "CTS",
    "cts",
    "CTS_class",
    "prediction"
  ),
  score_candidates = NULL
)


cts_raw <- cts_sheet$data

cts_sample_col <- cts_sheet$sample_col
cts_class_col <- cts_sheet$class_col


cat("Selected CTS worksheet:\n")
cat(cts_sheet$sheet, "\n\n")


cts_table <- tibble(
  sample_id = clean_sample_id(
    cts_raw[[cts_sample_col]]
  ),
  CTS = normalize_cts_label(
    cts_raw[[cts_class_col]]
  )
)


cts_table <- cts_table %>%
  filter(
    grepl(
      "^BP[0-9]+$",
      sample_id
    )
  ) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )


cat("CTS distribution from BP-only source:\n")

print(
  cts_table %>%
    count(
      CTS
    )
)

cat("\n")


if (nrow(cts_table) != 35) {
  
  stop(
    paste0(
      "Expected 35 BP CTS assignments, observed ",
      nrow(cts_table)
    )
  )
}


# ==============================================================================
# 15. HARD CHECK FINAL CTS
# ==============================================================================

bp_cts <- table(
  factor(
    cts_table$CTS,
    levels = c(
      "CTS1",
      "CTS2",
      "CTS3"
    )
  )
)


observed_bp_cts <- setNames(
  as.integer(bp_cts),
  c(
    "CTS1",
    "CTS2",
    "CTS3"
  )
)


expected_bp_cts <- c(
  CTS1 = 14L,
  CTS2 = 6L,
  CTS3 = 15L
)


if (!identical(
  observed_bp_cts,
  expected_bp_cts
)) {
  
  cat("Observed BP-only CTS:\n")
  print(observed_bp_cts)
  
  stop(
    "FINAL BP-only CTS validation FAILED."
  )
}


cat("FINAL BP-only CTS validation PASSED.\n")
cat("CTS1=14, CTS2=6, CTS3=15.\n\n")


# ==============================================================================
# 16. INTEGRATE SCORES + SRS + CTS
# ==============================================================================

blood_integrated <- score_table %>%
  left_join(
    srs_table,
    by = "sample_id"
  ) %>%
  left_join(
    cts_table,
    by = "sample_id"
  )


bp_integrated <- blood_integrated %>%
  filter(
    condition == "BP"
  )


if (nrow(bp_integrated) != 35) {
  stop("BP integration count is not 35.")
}


if (sum(!is.na(bp_integrated$SRS)) != 35) {
  stop("Incomplete BP SRS assignments.")
}


if (sum(!is.na(bp_integrated$CTS)) != 35) {
  stop("Incomplete BP CTS assignments.")
}


cat(
  "Scores + SRS + BP-only CTS integrated successfully.\n\n"
)


# ==============================================================================
# 17. HARD CHECK CTS x SRS
# ==============================================================================

cts_srs_table <- table(
  factor(
    bp_integrated$CTS,
    levels = c(
      "CTS1",
      "CTS2",
      "CTS3"
    )
  ),
  factor(
    bp_integrated$SRS,
    levels = c(
      "SRS1",
      "SRS2"
    )
  )
)


dimnames(cts_srs_table) <- list(
  CTS = c(
    "CTS1",
    "CTS2",
    "CTS3"
  ),
  SRS = c(
    "SRS1",
    "SRS2"
  )
)


expected_cts_srs <- matrix(
  c(
    14, 0,
    6, 0,
    8, 7
  ),
  nrow = 3,
  byrow = TRUE
)


dimnames(expected_cts_srs) <- list(
  CTS = c(
    "CTS1",
    "CTS2",
    "CTS3"
  ),
  SRS = c(
    "SRS1",
    "SRS2"
  )
)


cat("FINAL CTS x SRS table:\n")
print(cts_srs_table)
cat("\n")


if (!all(
  cts_srs_table ==
  expected_cts_srs
)) {
  
  cat("Expected CTS x SRS:\n")
  print(expected_cts_srs)
  
  stop(
    "FINAL CTS x SRS validation FAILED."
  )
}


cat(
  "FINAL CTS x SRS validation PASSED.\n\n"
)


# ==============================================================================
# 18. APPARENT BP vs BC DISCRIMINATION
# ==============================================================================

truth <- blood_integrated$truth_sepsis


primary_bp_bc_test <- wilcox.test(
  primary_5gene_score ~ condition,
  data = blood_integrated,
  exact = FALSE
)


dcaf_bp_bc_test <- wilcox.test(
  DCAF17_5gene_score ~ condition,
  data = blood_integrated,
  exact = FALSE
)


primary_auc <- calc_auc_with_ci(
  truth = truth,
  score = blood_integrated$primary_5gene_score
)


dcaf_auc <- calc_auc_with_ci(
  truth = truth,
  score = blood_integrated$DCAF17_5gene_score
)


dcaf17_auc <- calc_auc_with_ci(
  truth = truth,
  score = blood_integrated$DCAF17_single_score
)


descriptive_auc <- bind_rows(
  
  primary_auc %>%
    mutate(
      model = "Primary_5_gene_signed_score",
      comparison = "BP_vs_BC",
      wilcoxon_p =
        primary_bp_bc_test$p.value
    ),
  
  dcaf_auc %>%
    mutate(
      model = "DCAF17_5_gene_signed_score",
      comparison = "BP_vs_BC",
      wilcoxon_p =
        dcaf_bp_bc_test$p.value
    ),
  
  dcaf17_auc %>%
    mutate(
      model = "DCAF17_single_oriented",
      comparison = "BP_vs_BC",
      wilcoxon_p = NA_real_
    )
) %>%
  select(
    model,
    comparison,
    auc,
    ci_low,
    ci_high,
    wilcoxon_p
  )


cat("Descriptive apparent BP vs BC performance:\n")
print(descriptive_auc)
cat("\n")


# ==============================================================================
# 19. INDIVIDUAL GENE PERFORMANCE
# ==============================================================================

candidate_gene_results <- vector(
  "list",
  length(
    selected_genes_detected
  )
)


for (i in seq_along(
  selected_genes_detected
)) {
  
  gene <- selected_genes_detected[i]
  
  gene_values <- as.numeric(
    logcpm[
      gene,
      ,
      drop = TRUE
    ]
  )
  
  median_bp <- median(
    gene_values[
      blood_group == "BP"
    ],
    na.rm = TRUE
  )
  
  median_bc <- median(
    gene_values[
      blood_group == "BC"
    ],
    na.rm = TRUE
  )
  
  orientation <- ifelse(
    median_bp >= median_bc,
    1,
    -1
  )
  
  oriented_values <- gene_values *
    orientation
  
  auc_info <- calc_auc_with_ci(
    truth = truth,
    score = oriented_values
  )
  
  wt <- wilcox.test(
    gene_values ~ blood_group,
    exact = FALSE
  )
  
  candidate_gene_results[[i]] <- tibble(
    gene = gene,
    median_logCPM_BP = median_bp,
    median_logCPM_BC = median_bc,
    
    direction_BP_vs_BC = ifelse(
      median_bp >= median_bc,
      "UP_in_BP",
      "DOWN_in_BP"
    ),
    
    apparent_auc_oriented =
      auc_info$auc,
    
    auc_ci_low =
      auc_info$ci_low,
    
    auc_ci_high =
      auc_info$ci_high,
    
    wilcoxon_p =
      wt$p.value
  )
}


candidate_gene_performance <- bind_rows(
  candidate_gene_results
)


cat("Individual candidate-gene performance:\n")
print(candidate_gene_performance)
cat("\n")


# ==============================================================================
# 20. REPEATED CROSS-VALIDATION
# ==============================================================================

cv_repeats <- 100L
cv_folds <- 5L
cv_seed <- 20260817L


cat(
  "Running ",
  cv_repeats,
  " repeats of ",
  cv_folds,
  "-fold stratified CV...\n\n",
  sep = ""
)


cv_primary_signed <- run_repeated_cv(
  expr = logcpm,
  truth = truth,
  model_name =
    "Primary_5_gene_signed_score",
  genes = primary_genes,
  model_type = "signed_score",
  up_genes = primary_up,
  down_genes = primary_down,
  n_repeats = cv_repeats,
  n_folds = cv_folds,
  seed = cv_seed
)


cv_dcaf_signed <- run_repeated_cv(
  expr = logcpm,
  truth = truth,
  model_name =
    "DCAF17_5_gene_signed_score",
  genes = dcaf_genes,
  model_type = "signed_score",
  up_genes = dcaf_up,
  down_genes = dcaf_down,
  n_repeats = cv_repeats,
  n_folds = cv_folds,
  seed = cv_seed + 1L
)


cv_primary_ridge <- run_repeated_cv(
  expr = logcpm,
  truth = truth,
  model_name =
    "Primary_5_gene_ridge",
  genes = primary_genes,
  model_type = "ridge",
  n_repeats = cv_repeats,
  n_folds = cv_folds,
  seed = cv_seed + 2L
)


cv_dcaf_ridge <- run_repeated_cv(
  expr = logcpm,
  truth = truth,
  model_name =
    "DCAF17_5_gene_ridge",
  genes = dcaf_genes,
  model_type = "ridge",
  n_repeats = cv_repeats,
  n_folds = cv_folds,
  seed = cv_seed + 3L
)


cv_list <- list(
  cv_primary_signed,
  cv_dcaf_signed,
  cv_primary_ridge,
  cv_dcaf_ridge
)


if (length(
  missing_septicyte
) == 0) {
  
  cat(
    "Running SeptiCyte-related ",
    "four-gene ridge benchmark.\n\n"
  )
  
  cv_septicyte_ridge <- run_repeated_cv(
    expr = logcpm,
    truth = truth,
    model_name =
      "SeptiCyte_related_4_gene_ridge",
    genes = septicyte_genes,
    model_type = "ridge",
    n_repeats = cv_repeats,
    n_folds = cv_folds,
    seed = cv_seed + 4L
  )
  
  cv_list <- c(
    cv_list,
    list(
      cv_septicyte_ridge
    )
  )
}


cv_auc_summary <- bind_rows(
  lapply(
    cv_list,
    function(x) x$auc_summary
  )
)


cv_auc_by_repeat <- bind_rows(
  lapply(
    cv_list,
    function(x) x$auc_table
  )
)


cv_all_predictions <- bind_rows(
  lapply(
    cv_list,
    function(x) x$predictions
  )
)


cv_sample_predictions <- bind_rows(
  lapply(
    cv_list,
    function(x) x$sample_predictions
  )
)


cat("Cross-validation summary:\n")
print(cv_auc_summary)
cat("\n")


# ==============================================================================
# 21. PRIMARY SCORE vs SRS
# ==============================================================================

srs_sepsis <- bp_integrated %>%
  filter(
    SRS %in% c(
      "SRS1",
      "SRS2"
    )
  )


srs_wilcox <- wilcox.test(
  primary_5gene_score ~ SRS,
  data = srs_sepsis,
  exact = FALSE
)


srs_summary <- srs_sepsis %>%
  group_by(
    SRS
  ) %>%
  summarise(
    n = n(),
    
    median_score = median(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    IQR_score = IQR(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    mean_score = mean(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    sd_score = sd(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


srs_auc <- calc_auc(
  truth = ifelse(
    srs_sepsis$SRS == "SRS1",
    1L,
    0L
  ),
  score =
    srs_sepsis$primary_5gene_score
)


srs_test_summary <- tibble(
  comparison =
    "Primary_5gene_score_SRS1_vs_SRS2",
  
  test =
    "Wilcoxon_rank_sum",
  
  n =
    nrow(srs_sepsis),
  
  p_value =
    srs_wilcox$p.value,
  
  apparent_auc_SRS1_vs_SRS2 =
    srs_auc
)


cat("Primary 5-gene score by SRS:\n")
print(srs_summary)
print(srs_test_summary)
cat("\n")


# ==============================================================================
# 22. PRIMARY SCORE vs SRSq
# ==============================================================================

srsq_test <- cor.test(
  bp_integrated$primary_5gene_score,
  bp_integrated$SRSq,
  method = "spearman",
  exact = FALSE
)


primary_vs_srsq <- tibble(
  comparison =
    "Primary_5gene_score_vs_SRSq",
  
  n =
    sum(
      complete.cases(
        bp_integrated$primary_5gene_score,
        bp_integrated$SRSq
      )
    ),
  
  spearman_rho =
    as.numeric(
      srsq_test$estimate
    ),
  
  p_value =
    srsq_test$p.value
)


cat("Primary 5-gene score vs SRSq:\n")
print(primary_vs_srsq)
cat("\n")


# ==============================================================================
# 23. SRSq SENSITIVITY WITHOUT mNN OUTLIERS
# ==============================================================================

bp_no_mnn <- bp_integrated %>%
  filter(
    !sample_id %in%
      known_mnn_outliers
  )


if (nrow(bp_no_mnn) != 31) {
  
  warning(
    paste0(
      "Expected n=31 after mNN exclusion, observed ",
      nrow(bp_no_mnn)
    )
  )
}


srsq_no_mnn_test <- cor.test(
  bp_no_mnn$primary_5gene_score,
  bp_no_mnn$SRSq,
  method = "spearman",
  exact = FALSE
)


primary_vs_srsq_no_mnn <- tibble(
  comparison =
    "Primary_5gene_score_vs_SRSq_without_mNN_outliers",
  
  n =
    sum(
      complete.cases(
        bp_no_mnn$primary_5gene_score,
        bp_no_mnn$SRSq
      )
    ),
  
  spearman_rho =
    as.numeric(
      srsq_no_mnn_test$estimate
    ),
  
  p_value =
    srsq_no_mnn_test$p.value
)


cat("SRSq sensitivity without mNN outliers:\n")
print(primary_vs_srsq_no_mnn)
cat("\n")


# ==============================================================================
# 24. PRIMARY SCORE vs CTS
# ==============================================================================

cts_kw <- kruskal.test(
  primary_5gene_score ~ CTS,
  data = bp_integrated
)


cts_epsilon2 <- epsilon_squared_kw(
  x =
    bp_integrated$primary_5gene_score,
  
  group =
    bp_integrated$CTS
)


cts_summary <- bp_integrated %>%
  group_by(
    CTS
  ) %>%
  summarise(
    n = n(),
    
    median_score = median(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    IQR_score = IQR(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    mean_score = mean(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    sd_score = sd(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


cts_pairwise <- pairwise.wilcox.test(
  x =
    bp_integrated$primary_5gene_score,
  
  g =
    bp_integrated$CTS,
  
  p.adjust.method =
    "BH",
  
  exact =
    FALSE
)


cts_pairwise_table <- as.data.frame(
  as.table(
    cts_pairwise$p.value
  )
)


cts_pairwise_table <- cts_pairwise_table %>%
  filter(
    !is.na(Freq)
  ) %>%
  rename(
    group1 = Var1,
    group2 = Var2,
    BH_adjusted_p = Freq
  )


cts_test_summary <- tibble(
  comparison =
    "Primary_5gene_score_by_CTS",
  
  test =
    "Kruskal_Wallis",
  
  n =
    nrow(bp_integrated),
  
  p_value =
    cts_kw$p.value,
  
  epsilon_squared =
    cts_epsilon2
)


cat("Primary score by CTS:\n")
print(cts_summary)
print(cts_test_summary)

cat("Pairwise CTS comparisons:\n")
print(cts_pairwise_table)
cat("\n")


# ==============================================================================
# 25. INTEGRATED CTS/SRS STATES
# ==============================================================================

bp_integrated <- bp_integrated %>%
  mutate(
    CTS_SRS_group = paste(
      CTS,
      SRS,
      sep = "/"
    )
  )


integrated_summary <- bp_integrated %>%
  group_by(
    CTS_SRS_group
  ) %>%
  summarise(
    n = n(),
    
    median_primary_score = median(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    mean_primary_score = mean(
      primary_5gene_score,
      na.rm = TRUE
    ),
    
    median_myeloid_UP_score = median(
      myeloid_UP_score,
      na.rm = TRUE
    ),
    
    median_adaptive_suppression_score = median(
      adaptive_suppression_score,
      na.rm = TRUE
    ),
    
    median_SRSq = median(
      SRSq,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      mean_primary_score
    )
  )


integrated_kw <- kruskal.test(
  primary_5gene_score ~ CTS_SRS_group,
  data = bp_integrated
)


integrated_test_summary <- tibble(
  comparison =
    "Primary_5gene_score_by_integrated_CTS_SRS",
  
  test =
    "Kruskal_Wallis",
  
  n =
    nrow(bp_integrated),
  
  p_value =
    integrated_kw$p.value,
  
  epsilon_squared =
    epsilon_squared_kw(
      x =
        bp_integrated$primary_5gene_score,
      
      group =
        bp_integrated$CTS_SRS_group
    )
)


cat("Integrated CTS/SRS states:\n")
print(integrated_summary)
print(integrated_test_summary)
cat("\n")


# ==============================================================================
# 26. PRIMARY PANEL vs DCAF17 PANEL
# ==============================================================================

panel_correlation <- cor.test(
  bp_integrated$primary_5gene_score,
  bp_integrated$DCAF17_5gene_score,
  method = "spearman",
  exact = FALSE
)


panel_comparison <- tibble(
  comparison =
    "Primary_5gene_vs_DCAF17_5gene",
  
  n =
    nrow(bp_integrated),
  
  spearman_rho =
    as.numeric(
      panel_correlation$estimate
    ),
  
  p_value =
    panel_correlation$p.value
)


cat("Primary vs DCAF17 panel:\n")
print(panel_comparison)
cat("\n")


# ==============================================================================
# 27. INDIVIDUAL PRIMARY GENES vs SRSq
# ==============================================================================

primary_gene_srsq_list <- vector(
  "list",
  length(
    primary_genes
  )
)


for (i in seq_along(
  primary_genes
)) {
  
  gene <- primary_genes[i]
  
  gene_expression <- as.numeric(
    logcpm[
      gene,
      bp_integrated$sample_id,
      drop = TRUE
    ]
  )
  
  gene_test <- cor.test(
    gene_expression,
    bp_integrated$SRSq,
    method = "spearman",
    exact = FALSE
  )
  
  primary_gene_srsq_list[[i]] <- tibble(
    gene = gene,
    
    n =
      sum(
        complete.cases(
          gene_expression,
          bp_integrated$SRSq
        )
      ),
    
    spearman_rho =
      as.numeric(
        gene_test$estimate
      ),
    
    p_value =
      gene_test$p.value
  )
}


primary_gene_srsq <- bind_rows(
  primary_gene_srsq_list
)


primary_gene_srsq$BH_adjusted_p <- p.adjust(
  primary_gene_srsq$p_value,
  method = "BH"
)


cat("Individual primary genes vs SRSq:\n")
print(primary_gene_srsq)
cat("\n")


# ==============================================================================
# 28. SELECTED EXPRESSION TABLE
# ==============================================================================

selected_expression <- as.data.frame(
  t(
    logcpm[
      selected_genes_detected,
      ,
      drop = FALSE
    ]
  )
)


selected_expression$sample_id <- rownames(
  selected_expression
)


selected_expression <- selected_expression %>%
  relocate(
    sample_id
  ) %>%
  left_join(
    blood_integrated,
    by = "sample_id"
  )


# ==============================================================================
# 29. FIGURES
# ==============================================================================

p1 <- ggplot(
  blood_integrated,
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
    size = 2.2
  ) +
  labs(
    title =
      "Primary five-gene blood host-response score",
    
    subtitle =
      "CD177 + HK3 + IRAK3 versus CARD11 + IKZF2",
    
    x = NULL,
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  theme_bw(
    base_size = 12
  )


ggsave(
  file.path(
    figures_dir,
    "135_Figure_A_primary_score_BP_vs_BC.png"
  ),
  p1,
  width = 6.5,
  height = 5.5,
  dpi = 300
)


p2 <- ggplot(
  bp_integrated,
  aes(
    x = SRS,
    y = primary_5gene_score
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.2
  ) +
  labs(
    title =
      "Five-gene score across blood SRS classes",
    
    x = "SRS",
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  theme_bw(
    base_size = 12
  )


ggsave(
  file.path(
    figures_dir,
    "135_Figure_B_primary_score_by_SRS.png"
  ),
  p2,
  width = 6.5,
  height = 5.5,
  dpi = 300
)


p3 <- ggplot(
  bp_integrated,
  aes(
    x = SRSq,
    y = primary_5gene_score
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE
  ) +
  labs(
    title =
      "Five-gene host-response score and SRSq",
    
    x = "SRSq",
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  theme_bw(
    base_size = 12
  )


ggsave(
  file.path(
    figures_dir,
    "135_Figure_C_primary_score_vs_SRSq.png"
  ),
  p3,
  width = 6.5,
  height = 5.5,
  dpi = 300
)


p4 <- ggplot(
  bp_integrated,
  aes(
    x = CTS,
    y = primary_5gene_score
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.2
  ) +
  labs(
    title =
      "Five-gene score across Consensus Transcriptomic Subtypes",
    
    x = "CTS",
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  theme_bw(
    base_size = 12
  )


ggsave(
  file.path(
    figures_dir,
    "135_Figure_D_primary_score_by_CTS.png"
  ),
  p4,
  width = 6.5,
  height = 5.5,
  dpi = 300
)


p5 <- ggplot(
  bp_integrated,
  aes(
    x = CTS_SRS_group,
    y = primary_5gene_score
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.2
  ) +
  labs(
    title =
      "Five-gene score across integrated CTS/SRS states",
    
    x =
      "CTS/SRS state",
    
    y =
      "Myeloid-adaptive balance score"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )


ggsave(
  file.path(
    figures_dir,
    "135_Figure_E_primary_score_integrated_CTS_SRS.png"
  ),
  p5,
  width = 8,
  height = 5.5,
  dpi = 300
)


p6 <- ggplot(
  cv_auc_by_repeat,
  aes(
    x = model,
    y = auc
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.25,
    size = 1
  ) +
  coord_cartesian(
    ylim = c(
      0.5,
      1
    )
  ) +
  labs(
    title =
      "Repeated stratified cross-validation",
    
    subtitle =
      paste0(
        cv_repeats,
        " repeats of ",
        cv_folds,
        "-fold CV; internal resampling only"
      ),
    
    x = NULL,
    
    y = "AUROC"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )


ggsave(
  file.path(
    figures_dir,
    "135_Figure_F_cross_validation_AUC.png"
  ),
  p6,
  width = 9,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 30. SAVE CSV TABLES
# ==============================================================================

write.csv(
  blood_integrated,
  file.path(
    tables_dir,
    "135_blood_scores_with_final_endotypes.csv"
  ),
  row.names = FALSE
)


write.csv(
  bp_integrated,
  file.path(
    tables_dir,
    "135_sepsis_blood_scores_with_SRS_CTS.csv"
  ),
  row.names = FALSE
)


write.csv(
  descriptive_auc,
  file.path(
    tables_dir,
    "135_descriptive_apparent_AUC.csv"
  ),
  row.names = FALSE
)


write.csv(
  candidate_gene_performance,
  file.path(
    tables_dir,
    "135_candidate_gene_descriptive_performance.csv"
  ),
  row.names = FALSE
)


write.csv(
  cv_auc_summary,
  file.path(
    tables_dir,
    "135_cross_validation_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  cv_auc_by_repeat,
  file.path(
    tables_dir,
    "135_cross_validation_AUC_by_repeat.csv"
  ),
  row.names = FALSE
)


write.csv(
  cv_all_predictions,
  file.path(
    tables_dir,
    "135_cross_validation_all_predictions.csv"
  ),
  row.names = FALSE
)


write.csv(
  cv_sample_predictions,
  file.path(
    tables_dir,
    "135_cross_validation_sample_predictions.csv"
  ),
  row.names = FALSE
)


write.csv(
  primary_gene_srsq,
  file.path(
    tables_dir,
    "135_primary_gene_vs_SRSq.csv"
  ),
  row.names = FALSE
)


write.csv(
  selected_expression,
  file.path(
    tables_dir,
    "135_selected_gene_expression_logCPM.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 31. EXCEL WORKBOOK
# ==============================================================================

wb <- openxlsx::createWorkbook()


run_info <- tibble(
  parameter = c(
    "script",
    "run_date",
    "analysis_scope",
    "BP_n",
    "BC_n",
    "CTS_source",
    "CV_repeats",
    "CV_folds",
    "primary_UP",
    "primary_DOWN",
    "DCAF17_UP",
    "DCAF17_DOWN",
    "SeptiCyte_related",
    "critical_limitation"
  ),
  
  value = c(
    script_name,
    as.character(run_date),
    "Blood only",
    "35",
    "10",
    cts_file,
    as.character(cv_repeats),
    as.character(cv_folds),
    paste(primary_up, collapse = "; "),
    paste(primary_down, collapse = "; "),
    paste(dcaf_up, collapse = "; "),
    paste(dcaf_down, collapse = "; "),
    paste(septicyte_genes, collapse = "; "),
    paste0(
      "Primary five-gene panel was previously selected ",
      "using this cohort; CV is internal."
    )
  )
)


sheet_list <- list(
  "00_run_info" =
    run_info,
  
  "01_input_manifest" =
    input_manifest,
  
  "02_gene_panels" =
    panel_definition,
  
  "03_gene_availability" =
    gene_availability,
  
  "04_blood_scores" =
    blood_integrated,
  
  "05_BP_endotype_scores" =
    bp_integrated,
  
  "06_FINAL_SRS" =
    srs_table,
  
  "07_FINAL_CTS_BP_only" =
    cts_table,
  
  "09_apparent_AUC" =
    descriptive_auc,
  
  "10_gene_performance" =
    candidate_gene_performance,
  
  "11_CV_summary" =
    cv_auc_summary,
  
  "12_CV_AUC_repeats" =
    cv_auc_by_repeat,
  
  "13_CV_predictions" =
    cv_sample_predictions,
  
  "14_primary_by_SRS" =
    srs_summary,
  
  "15_primary_SRS_test" =
    srs_test_summary,
  
  "16_primary_vs_SRSq" =
    bind_rows(
      primary_vs_srsq,
      primary_vs_srsq_no_mnn
    ),
  
  "17_primary_by_CTS" =
    cts_summary,
  
  "18_primary_CTS_test" =
    cts_test_summary,
  
  "19_CTS_pairwise" =
    cts_pairwise_table,
  
  "20_integrated_CTS_SRS" =
    integrated_summary,
  
  "21_integrated_test" =
    integrated_test_summary,
  
  "22_panel_comparison" =
    panel_comparison,
  
  "23_primary_gene_SRSq" =
    primary_gene_srsq,
  
  "24_selected_expression" =
    selected_expression
)


for (sheet_name in names(
  sheet_list
)) {
  
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


openxlsx::addWorksheet(
  wb,
  "08_FINAL_CTSxSRS"
)


openxlsx::writeData(
  wb,
  "08_FINAL_CTSxSRS",
  as.data.frame.matrix(
    cts_srs_table
  ),
  rowNames = TRUE
)


workbook_file <- file.path(
  tables_dir,
  "135_blood_endotype_biomarker_validation.xlsx"
)


openxlsx::saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)


# ==============================================================================
# 32. SUMMARY TEXT
# ==============================================================================

primary_signed_cv <- cv_auc_summary %>%
  filter(
    model ==
      "Primary_5_gene_signed_score"
  )


primary_ridge_cv <- cv_auc_summary %>%
  filter(
    model ==
      "Primary_5_gene_ridge"
  )


summary_ru <- c(
  "SCRIPT 135 — BLOOD-ONLY ENDOTYPE-INFORMED BIOMARKER VALIDATION",
  "====================================================================",
  "",
  "Только кровь.",
  "BP n=35.",
  "BC n=10.",
  "",
  "FINAL SRS:",
  "BP: SRS1=28; SRS2=7; SRS3=0.",
  "BC: SRS1=0; SRS2=3; SRS3=7.",
  "",
  "FINAL BP-ONLY CTS:",
  "CTS1=14; CTS2=6; CTS3=15.",
  "",
  "CTS x SRS:",
  "CTS1/SRS1=14",
  "CTS2/SRS1=6",
  "CTS3/SRS1=8",
  "CTS3/SRS2=7",
  "",
  "PRIMARY FIVE-GENE PANEL:",
  "UP: CD177, HK3, IRAK3",
  "DOWN: CARD11, IKZF2",
  "",
  paste0(
    "Primary score vs SRS: p=",
    signif(
      srs_wilcox$p.value,
      5
    ),
    "; AUC=",
    round(
      srs_auc,
      3
    )
  ),
  "",
  paste0(
    "Primary score vs SRSq: rho=",
    round(
      primary_vs_srsq$spearman_rho,
      3
    ),
    "; p=",
    signif(
      primary_vs_srsq$p_value,
      5
    )
  ),
  "",
  paste0(
    "Without mNN outliers: rho=",
    round(
      primary_vs_srsq_no_mnn$spearman_rho,
      3
    ),
    "; p=",
    signif(
      primary_vs_srsq_no_mnn$p_value,
      5
    )
  ),
  "",
  paste0(
    "Primary score by CTS: Kruskal-Wallis p=",
    signif(
      cts_kw$p.value,
      5
    ),
    "; epsilon-squared=",
    round(
      cts_epsilon2,
      3
    )
  ),
  "",
  paste0(
    "Primary signed-score repeated-CV median AUC=",
    round(
      primary_signed_cv$median_repeat_auc,
      3
    ),
    "; empirical range ",
    round(
      primary_signed_cv$repeat_auc_q025,
      3
    ),
    "–",
    round(
      primary_signed_cv$repeat_auc_q975,
      3
    )
  ),
  "",
  paste0(
    "Primary ridge repeated-CV median AUC=",
    round(
      primary_ridge_cv$median_repeat_auc,
      3
    )
  ),
  "",
  "IMPORTANT:",
  paste0(
    "Primary five-gene panel was originally selected ",
    "using this same cohort."
  ),
  paste0(
    "Repeated CV therefore measures internal stability ",
    "and is NOT independent validation."
  ),
  "",
  "Recommended terminology:",
  "candidate endotype-informed five-gene blood signature"
)


summary_ru_file <- file.path(
  text_dir,
  "135_summary_RU.txt"
)


writeLines(
  summary_ru,
  summary_ru_file
)


summary_en <- c(
  "SCRIPT 135 — BLOOD-ONLY ENDOTYPE-INFORMED BIOMARKER VALIDATION",
  "====================================================================",
  "",
  "Blood only.",
  "BP n=35.",
  "BC n=10.",
  "",
  "Final SRS:",
  "BP: SRS1=28; SRS2=7; SRS3=0.",
  "BC: SRS1=0; SRS2=3; SRS3=7.",
  "",
  "Final BP-only CTS:",
  "CTS1=14; CTS2=6; CTS3=15.",
  "",
  "CTS x SRS:",
  "CTS1/SRS1=14",
  "CTS2/SRS1=6",
  "CTS3/SRS1=8",
  "CTS3/SRS2=7",
  "",
  "Recommended terminology:",
  "candidate endotype-informed five-gene blood signature",
  "",
  paste0(
    "Critical limitation: the five-gene panel was ",
    "previously selected using this same cohort."
  )
)


summary_en_file <- file.path(
  text_dir,
  "135_summary_EN.txt"
)


writeLines(
  summary_en,
  summary_en_file
)


# ==============================================================================
# 33. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    logs_dir,
    "135_sessionInfo.txt"
  )
)


run_parameters <- tibble(
  parameter = c(
    "script",
    "run_date",
    "BP_n",
    "BC_n",
    "CTS_file",
    "CV_repeats",
    "CV_folds",
    "CV_seed"
  ),
  
  value = c(
    script_name,
    as.character(run_date),
    "35",
    "10",
    cts_file,
    as.character(cv_repeats),
    as.character(cv_folds),
    as.character(cv_seed)
  )
)


write.csv(
  run_parameters,
  file.path(
    logs_dir,
    "135_run_parameters.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 34. FINAL CONSOLE REPORT
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 135 completed successfully.\n")
cat("====================================================================\n\n")


cat("FINAL ENDOTYPE CHECKS:\n")
cat("SRS BP: SRS1=28, SRS2=7, SRS3=0 — PASSED\n")
cat("SRS BC: SRS1=0, SRS2=3, SRS3=7 — PASSED\n")
cat("CTS BP-only: CTS1=14, CTS2=6, CTS3=15 — PASSED\n")
cat("CTS x SRS — PASSED\n\n")


cat("CTS source used:\n")
cat(
  normalizePath(
    cts_file,
    winslash = "/",
    mustWork = FALSE
  ),
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
      "135_cross_validation_summary.csv"
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
      "135_sepsis_blood_scores_with_SRS_CTS.csv"
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
    summary_ru_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "5) ",
  normalizePath(
    summary_en_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat("Open first in workbook:\n")
cat("  08_FINAL_CTSxSRS\n")
cat("  11_CV_summary\n")
cat("  14_primary_by_SRS\n")
cat("  16_primary_vs_SRSq\n")
cat("  17_primary_by_CTS\n")
cat("  19_CTS_pairwise\n")
cat("  20_integrated_CTS_SRS\n")
cat("  23_primary_gene_SRSq\n\n")


cat("IMPORTANT:\n")
cat("- Blood only.\n")
cat("- CTS source = FINAL BP-only classification.\n")
cat("- Do NOT use 107_CTS_predictions_BP_BC.xlsx for final CTS.\n")
cat("- No urine.\n")
cat("- No lncRNA.\n")
cat("- CV is internal, not external validation.\n")
cat("- BP vs BC = sepsis vs healthy-control discrimination.\n\n")

cat("Done.\n")