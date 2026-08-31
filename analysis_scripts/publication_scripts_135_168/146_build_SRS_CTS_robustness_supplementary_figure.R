################################################################################
# Script 146
# Supplementary publication figure:
# Robustness of blood transcriptomic endotype assignments
#
# Project: Sepsis_DESeq2
#
# PURPOSE
# -------
# Build a publication-ready supplementary figure summarizing already completed
# robustness analyses from:
#
#   Script 40b:
#     SRS batch and k sensitivity in all blood samples (BP + BC)
#
#   Script 40c:
#     BP-only SRS cohort-composition and within-BP batch sensitivity
#
#   Script 107d:
#     CTS random-forest reproducibility, BP-only preprocessing sensitivity,
#     batch-adjustment sensitivity, and CTS x sequencing-batch association
#
# IMPORTANT
# ---------
#   - NO SRS classification is recalculated.
#   - NO CTS classification is recalculated.
#   - NO random forest is rerun.
#   - NO batch correction is rerun.
#   - NO inferential statistics are recalculated.
#   - All numerical results are read from finalized output tables.
#
# Figure layout
# -------------
#   A. SRS categorical agreement across k
#   B. SRSq correlation across k
#   C. SRS sensitivity to within-BP batch removal
#   D. CTS stability across 50 random-forest refits
#   E. CTS modal agreement across preprocessing strategies
#   F. CTS x sequencing-batch distribution
#
# Output
# ------
# results/blood_endotypes_biomarkers/146_SRS_CTS_robustness_figure/
#
################################################################################


cat("====================================================================\n")
cat("Running Script 146\n")
cat("SRS / CTS robustness supplementary publication figure\n")
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
# 2. INPUT FILES
# =============================================================================

file_40b <- file.path(
  "results",
  "sepstratifier",
  "40b_batch_k_sensitivity",
  "tables",
  "40b_SRS_batch_k_sensitivity_audit.xlsx"
)

file_40c <- file.path(
  "results",
  "sepstratifier",
  "40c_BP_only_batch_sensitivity",
  "tables",
  "40c_BP_only_SRS_batch_sensitivity_audit.xlsx"
)

file_107d <- file.path(
  "results",
  "cts_consensus",
  "107d_reproducibility_batch_sensitivity",
  "tables",
  "107d_CTS_reproducibility_batch_sensitivity_audit.xlsx"
)


required_files <- c(
  file_40b,
  file_40c,
  file_107d
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


cat("\nInput files:\n")

for (f in required_files) {
  
  print(
    normalizePath(
      f,
      winslash = "\\",
      mustWork = TRUE
    )
  )
}


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

out_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "146_SRS_CTS_robustness_figure"
)

figure_dir <- file.path(
  out_dir,
  "figures"
)

table_dir <- file.path(
  out_dir,
  "tables"
)

text_dir <- file.path(
  out_dir,
  "text"
)


for (d in c(
  out_dir,
  figure_dir,
  table_dir,
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
# 4. PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "openxlsx",
  "readr",
  "scales"
)


for (pkg in required_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    stop(
      "Required package not installed: ",
      pkg
    )
  }
}


suppressPackageStartupMessages({
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(readr)
  library(scales)
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
    "146_package_versions.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file = file.path(
    text_dir,
    "146_sessionInfo.txt"
  )
)


# =============================================================================
# 5. PUBLICATION COLORS
# =============================================================================

col_original <- "#2166AC"

col_adjusted <- "#B2182B"

col_neutral <- "#636363"


col_preprocessing <- c(
  "Primary" = "#2166AC",
  "BP-only VST" = "#67A9CF",
  "BP-only + batch removal" = "#B2182B"
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
      
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = base_size + 2
      ),
      
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      
      legend.title = ggplot2::element_blank(),
      
      legend.position = "top",
      
      plot.margin = ggplot2::margin(
        7,
        7,
        7,
        7
      )
    )
}


# -----------------------------------------------------------------------------
# Read workbook sheet by required columns rather than hard-coded sheet name.
# -----------------------------------------------------------------------------

read_sheet_with_columns <- function(
    workbook,
    required_columns,
    preferred_pattern = NULL
) {
  
  sheets <- openxlsx::getSheetNames(
    workbook
  )
  
  
  if (!is.null(
    preferred_pattern
  )) {
    
    preferred <- sheets[
      grepl(
        preferred_pattern,
        sheets,
        ignore.case = TRUE
      )
    ]
    
    other <- setdiff(
      sheets,
      preferred
    )
    
    sheets <- c(
      preferred,
      other
    )
  }
  
  
  for (sheet in sheets) {
    
    dat <- tryCatch(
      
      openxlsx::read.xlsx(
        workbook,
        sheet = sheet
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    if (is.null(dat)) {
      next
    }
    
    
    names(dat) <- clean_names_simple(
      names(dat)
    )
    
    
    if (all(
      required_columns %in%
      names(dat)
    )) {
      
      return(
        list(
          data = dat,
          sheet = sheet
        )
      )
    }
  }
  
  
  stop(
    "Could not find a sheet in:\n",
    workbook,
    "\ncontaining required columns:\n",
    paste(
      required_columns,
      collapse = ", "
    )
  )
}


# -----------------------------------------------------------------------------
# Read CTS x batch sheet
# -----------------------------------------------------------------------------

read_cts_batch_sheet <- function(
    workbook
) {
  
  sheets <- openxlsx::getSheetNames(
    workbook
  )
  
  
  preferred <- sheets[
    grepl(
      "CTS.*batch|batch.*CTS",
      sheets,
      ignore.case = TRUE
    )
  ]
  
  
  sheets <- c(
    preferred,
    setdiff(
      sheets,
      preferred
    )
  )
  
  
  for (sheet in sheets) {
    
    dat <- tryCatch(
      
      openxlsx::read.xlsx(
        workbook,
        sheet = sheet
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    if (is.null(dat)) {
      next
    }
    
    
    names(dat) <- clean_names_simple(
      names(dat)
    )
    
    
    if (
      "CTS" %in%
      names(dat) &&
      any(
        grepl(
          "^chip",
          names(dat),
          ignore.case = TRUE
        )
      )
    ) {
      
      return(
        list(
          data = dat,
          sheet = sheet
        )
      )
    }
  }
  
  
  stop(
    "CTS x batch table was not found in ",
    workbook
  )
}


# =============================================================================
# 7. READ SCRIPT 40b DATA
# =============================================================================

res_40b_global <- read_sheet_with_columns(
  workbook = file_40b,
  required_columns = c(
    "preprocessing",
    "k",
    "SRS_agreement_percent",
    "SRSq_spearman_vs_primary"
  ),
  preferred_pattern = "global|sensitivity"
)


srs_40b <- res_40b_global$data


cat("\n40b global sensitivity sheet:\n")
print(res_40b_global$sheet)


srs_40b_original <- srs_40b %>%
  
  dplyr::filter(
    grepl(
      "Original",
      preprocessing,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::mutate(
    
    k = as.numeric(
      k
    ),
    
    SRS_agreement_percent =
      as.numeric(
        SRS_agreement_percent
      ),
    
    SRSq_spearman_vs_primary =
      as.numeric(
        SRSq_spearman_vs_primary
      )
  ) %>%
  
  dplyr::arrange(
    k
  )


if (nrow(srs_40b_original) == 0) {
  
  stop(
    "Could not identify Original_logCPM results in Script 40b."
  )
}


# =============================================================================
# 8. READ SCRIPT 40c DATA
# =============================================================================

res_40c_samek <- read_sheet_with_columns(
  workbook = file_40c,
  required_columns = c(
    "k",
    "SRS_agreement_percent",
    "SRSq_spearman"
  ),
  preferred_pattern = "same|original.*batch"
)


srs_40c_samek <- res_40c_samek$data %>%
  
  dplyr::mutate(
    
    k = as.numeric(
      k
    ),
    
    SRS_agreement_percent =
      as.numeric(
        SRS_agreement_percent
      ),
    
    SRSq_spearman =
      as.numeric(
        SRSq_spearman
      )
  ) %>%
  
  dplyr::arrange(
    k
  )


cat("\n40c same-k sheet:\n")
print(res_40c_samek$sheet)


row_40c_k20 <- srs_40c_samek %>%
  dplyr::filter(
    k == 20
  )


if (nrow(row_40c_k20) != 1) {
  
  stop(
    "Expected exactly one k=20 row in 40c same-k comparison."
  )
}


k20_agreement <- as.numeric(
  row_40c_k20$SRS_agreement_percent
)

k20_rho <- as.numeric(
  row_40c_k20$SRSq_spearman
)


# =============================================================================
# 9. READ SCRIPT 107d DATA
# =============================================================================

res_107d_stability <- read_sheet_with_columns(
  workbook = file_107d,
  required_columns = c(
    "preprocessing",
    "sample_id",
    "modal_CTS",
    "modal_percent",
    "modal_matches_original"
  ),
  preferred_pattern = "stability"
)


cts_stability <- res_107d_stability$data


res_107d_modal <- read_sheet_with_columns(
  workbook = file_107d,
  required_columns = c(
    "comparison",
    "agreement_n",
    "n_total",
    "agreement_percent"
  ),
  preferred_pattern = "modal.*compare|comparison"
)


cts_modal_compare <- res_107d_modal$data


res_107d_batch_assoc <- read_sheet_with_columns(
  workbook = file_107d,
  required_columns = c(
    "test",
    "p_value",
    "cramers_V",
    "n"
  ),
  preferred_pattern = "batch.*association"
)


cts_batch_assoc <- res_107d_batch_assoc$data


res_107d_batch_table <- read_cts_batch_sheet(
  workbook = file_107d
)


cts_batch_wide <- res_107d_batch_table$data


cat("\n107d sample stability sheet:\n")
print(res_107d_stability$sheet)

cat("\n107d modal comparison sheet:\n")
print(res_107d_modal$sheet)

cat("\n107d batch association sheet:\n")
print(res_107d_batch_assoc$sheet)

cat("\n107d CTS x batch sheet:\n")
print(res_107d_batch_table$sheet)


# =============================================================================
# 10. PREPARE CTS SEED-STABILITY DATA
# =============================================================================

cts_primary_stability <- cts_stability %>%
  
  dplyr::filter(
    preprocessing ==
      "Primary_BP_from_BP_BC_VST"
  ) %>%
  
  dplyr::mutate(
    
    modal_percent =
      as.numeric(
        modal_percent
      ),
    
    sample_id =
      as.character(
        sample_id
      )
  ) %>%
  
  dplyr::arrange(
    modal_percent,
    sample_id
  ) %>%
  
  dplyr::mutate(
    sample_index =
      dplyr::row_number()
  )


if (nrow(cts_primary_stability) != 35) {
  
  stop(
    "Expected 35 primary CTS stability rows; found ",
    nrow(cts_primary_stability)
  )
}


n_cts_100 <- sum(
  cts_primary_stability$modal_percent == 100
)

n_cts_ge95 <- sum(
  cts_primary_stability$modal_percent >= 95
)

min_cts_stability <- min(
  cts_primary_stability$modal_percent
)


# =============================================================================
# 11. PREPARE CTS PREPROCESSING COMPARISON
# =============================================================================

cts_preprocessing_plot <- cts_modal_compare %>%
  
  dplyr::filter(
    comparison %in% c(
      "Primary modal vs original Script 107",
      "BP-only modal vs original Script 107",
      "Batch-adjusted modal vs original Script 107"
    )
  ) %>%
  
  dplyr::mutate(
    
    preprocessing = dplyr::case_when(
      
      comparison ==
        "Primary modal vs original Script 107" ~
        "Primary",
      
      comparison ==
        "BP-only modal vs original Script 107" ~
        "BP-only VST",
      
      comparison ==
        "Batch-adjusted modal vs original Script 107" ~
        "BP-only + batch removal",
      
      TRUE ~ comparison
    ),
    
    preprocessing = factor(
      preprocessing,
      levels = c(
        "Primary",
        "BP-only VST",
        "BP-only + batch removal"
      )
    ),
    
    agreement_percent =
      as.numeric(
        agreement_percent
      ),
    
    agreement_n =
      as.numeric(
        agreement_n
      ),
    
    n_total =
      as.numeric(
        n_total
      )
  )


if (nrow(cts_preprocessing_plot) != 3) {
  
  stop(
    "Expected 3 CTS preprocessing comparison rows; found ",
    nrow(cts_preprocessing_plot)
  )
}


# =============================================================================
# 12. PREPARE CTS x BATCH TABLE
# =============================================================================

cts_batch_long <- cts_batch_wide %>%
  
  tidyr::pivot_longer(
    cols = -CTS,
    names_to = "batch",
    values_to = "n"
  ) %>%
  
  dplyr::mutate(
    
    CTS = factor(
      CTS,
      levels = c(
        "CTS1",
        "CTS2",
        "CTS3"
      )
    ),
    
    n = as.numeric(
      n
    )
  )


batch_p <- as.numeric(
  cts_batch_assoc$p_value[1]
)

batch_v <- as.numeric(
  cts_batch_assoc$cramers_V[1]
)


# =============================================================================
# 13. NUMERICAL QC
# =============================================================================

cat("\nNumerical audit:\n")

cat(
  "Minimum original-k SRS agreement: ",
  round(
    min(
      srs_40b_original$SRS_agreement_percent,
      na.rm = TRUE
    ),
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Minimum original-k SRSq rho: ",
  round(
    min(
      srs_40b_original$SRSq_spearman_vs_primary,
      na.rm = TRUE
    ),
    3
  ),
  "\n",
  sep = ""
)


cat(
  "BP-only batch sensitivity at k=20: agreement ",
  round(
    k20_agreement,
    1
  ),
  "%; rho ",
  round(
    k20_rho,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "CTS 100% seed-stable: ",
  n_cts_100,
  "/35\n",
  sep = ""
)


cat(
  "CTS >=95% seed-stable: ",
  n_cts_ge95,
  "/35\n",
  sep = ""
)


cat(
  "Minimum CTS stability: ",
  min_cts_stability,
  "%\n",
  sep = ""
)


cat(
  "CTS x batch Fisher p: ",
  signif(
    batch_p,
    6
  ),
  "\n",
  sep = ""
)


cat(
  "CTS x batch Cramer's V: ",
  round(
    batch_v,
    4
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 14. PANEL A
# SRS CATEGORICAL AGREEMENT ACROSS k
# =============================================================================

p_A <- ggplot2::ggplot(
  srs_40b_original,
  ggplot2::aes(
    x = k,
    y = SRS_agreement_percent
  )
) +
  
  ggplot2::annotate(
    "rect",
    xmin = 9,
    xmax = 14,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08,
    fill = col_original
  ) +
  
  ggplot2::geom_hline(
    yintercept = 95,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  
  ggplot2::geom_vline(
    xintercept = 20,
    linetype = "dotted",
    linewidth = 0.55,
    color = col_neutral
  ) +
  
  ggplot2::geom_line(
    linewidth = 0.8,
    color = col_original
  ) +
  
  ggplot2::geom_point(
    size = 2.7,
    color = col_original
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = srs_40b_original$k
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(
      94,
      101
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    tag = "A",
    
    title =
      "SRS classification is stable across k",
    
    subtitle =
      "Original TMM-logCPM; shaded region indicates the 20–30% k window",
    
    x =
      "Mutual-nearest-neighbour parameter (k)",
    
    y =
      "Agreement with primary SRS (%)"
  )


# =============================================================================
# 15. PANEL B
# SRSq CORRELATION ACROSS k
# =============================================================================

p_B <- ggplot2::ggplot(
  srs_40b_original,
  ggplot2::aes(
    x = k,
    y = SRSq_spearman_vs_primary
  )
) +
  
  ggplot2::annotate(
    "rect",
    xmin = 9,
    xmax = 14,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08,
    fill = col_original
  ) +
  
  ggplot2::geom_vline(
    xintercept = 20,
    linetype = "dotted",
    linewidth = 0.55,
    color = col_neutral
  ) +
  
  ggplot2::geom_line(
    linewidth = 0.8,
    color = col_original
  ) +
  
  ggplot2::geom_point(
    size = 2.7,
    color = col_original
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = srs_40b_original$k
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(
      0.92,
      1.01
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    tag = "B",
    
    title =
      "Quantitative SRSq is robust across k",
    
    subtitle =
      "Minimum Spearman correlation with the primary analysis was 0.938",
    
    x =
      "Mutual-nearest-neighbour parameter (k)",
    
    y =
      "Spearman ρ vs primary SRSq"
  )


# =============================================================================
# 16. PANEL C
# BP-ONLY ORIGINAL vs BATCH-ADJUSTED SRS
# =============================================================================

p_C <- ggplot2::ggplot(
  srs_40c_samek,
  ggplot2::aes(
    x = k,
    y = SRS_agreement_percent
  )
) +
  
  ggplot2::annotate(
    "rect",
    xmin = 7,
    xmax = 10,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.07,
    fill = col_adjusted
  ) +
  
  ggplot2::geom_vline(
    xintercept = 20,
    linetype = "dotted",
    linewidth = 0.55,
    color = col_neutral
  ) +
  
  ggplot2::geom_line(
    linewidth = 0.8,
    color = col_adjusted
  ) +
  
  ggplot2::geom_point(
    size = 2.7,
    color = col_adjusted
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = srs_40c_samek$k
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(
      74,
      91
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    tag = "C",
    
    title =
      "SRS sensitivity to within-BP batch removal",
    
    subtitle =
      paste0(
        "At k=20: ",
        round(
          k20_agreement,
          1
        ),
        "% categorical agreement; SRSq ρ=",
        round(
          k20_rho,
          3
        )
      ),
    
    x =
      "Mutual-nearest-neighbour parameter (k)",
    
    y =
      "Original vs batch-adjusted SRS agreement (%)"
  )


# =============================================================================
# 17. PANEL D
# CTS REPRODUCIBILITY ACROSS 50 RANDOM-FOREST REFITS
# =============================================================================

p_D <- ggplot2::ggplot(
  cts_primary_stability,
  ggplot2::aes(
    x = sample_index,
    y = modal_percent
  )
) +
  
  ggplot2::geom_hline(
    yintercept = 95,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey60"
  ) +
  
  ggplot2::geom_point(
    size = 2.9,
    color = col_original,
    alpha = 0.95
  ) +
  
  ggplot2::scale_y_continuous(
    breaks = c(
      95,
      96,
      97,
      98,
      99,
      100
    )
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(
      97.5,
      100.35
    )
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_blank(),
    
    axis.ticks.x =
      ggplot2::element_blank()
  ) +
  
  ggplot2::labs(
    tag = "D",
    
    title =
      "CTS is reproducible across random-forest refits",
    
    subtitle =
      paste0(
        n_cts_100,
        "/35 identical in all 50 runs; ",
        "minimum modal stability ",
        round(
          min_cts_stability,
          0
        ),
        "%"
      ),
    
    x =
      "Sepsis blood samples",
    
    y =
      "Modal CTS frequency across 50 seeds (%)"
  )


# =============================================================================
# 18. PANEL E
# CTS PREPROCESSING SENSITIVITY
# =============================================================================

p_E <- ggplot2::ggplot(
  cts_preprocessing_plot,
  ggplot2::aes(
    x = preprocessing,
    y = agreement_percent,
    fill = preprocessing
  )
) +
  
  ggplot2::geom_col(
    width = 0.64
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        agreement_n,
        "/",
        n_total,
        "\n",
        round(
          agreement_percent,
          1
        ),
        "%"
      )
    ),
    vjust = -0.32,
    size = 3.4,
    fontface = "bold"
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_preprocessing
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(
      0,
      112
    ),
    clip = "off"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    axis.text.x =
      ggplot2::element_text(
        angle = 14,
        hjust = 1,
        size = 9
      )
  ) +
  
  ggplot2::labs(
    tag = "E",
    
    title =
      "CTS is stable to BP-only variance stabilization",
    
    subtitle =
      "Batch removal was retained as a conservative sensitivity analysis",
    
    x = NULL,
    
    y =
      "Modal agreement with original CTS (%)"
  )


# =============================================================================
# 19. PANEL F
# CTS x SEQUENCING BATCH
# =============================================================================

p_F <- ggplot2::ggplot(
  cts_batch_long,
  ggplot2::aes(
    x = batch,
    y = CTS,
    fill = n
  )
) +
  
  ggplot2::geom_tile(
    linewidth = 1,
    color = "white"
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = n
    ),
    fontface = "bold",
    size = 3.6
  ) +
  
  ggplot2::scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2166AC"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_text(
        angle = 28,
        hjust = 1,
        size = 8.5
      ),
    
    legend.position =
      "right"
  ) +
  
  ggplot2::labs(
    tag = "F",
    
    title =
      "CTS distribution across sequencing batches",
    
    subtitle =
      paste0(
        "Monte Carlo Fisher p=",
        format_p(
          batch_p
        ),
        "; Cramér's V=",
        round(
          batch_v,
          3
        )
      ),
    
    x =
      "Sequencing batch",
    
    y =
      "Consensus Transcriptomic Subtype",
    
    fill =
      "Patients"
  )


# =============================================================================
# 20. ASSEMBLE SUPPLEMENTARY FIGURE
# =============================================================================

row1 <- (
  p_A |
    p_B |
    p_C
) +
  
  patchwork::plot_layout(
    widths = c(
      1,
      1,
      1
    )
  )


row2 <- (
  p_D |
    p_E |
    p_F
) +
  
  patchwork::plot_layout(
    widths = c(
      1.05,
      0.95,
      1.10
    )
  )


main_figure <- (
  row1 /
    row2
) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1.02
    )
  ) +
  
  patchwork::plot_annotation(
    
    title =
      "Robustness of blood transcriptomic endotype assignments",
    
    subtitle =
      paste0(
        "Sensitivity of SRS and CTS to parameter choice, cohort composition, ",
        "stochastic refitting, and batch preprocessing"
      ),
    
    theme = ggplot2::theme(
      
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 16
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 11,
        margin = ggplot2::margin(
          b = 9
        )
      )
    )
  )


# =============================================================================
# 21. EXPORT MAIN FIGURE
# =============================================================================

main_png <- file.path(
  figure_dir,
  "146_Supplementary_Figure_Sx_SRS_CTS_robustness.png"
)

main_pdf <- file.path(
  figure_dir,
  "146_Supplementary_Figure_Sx_SRS_CTS_robustness.pdf"
)

main_tiff <- file.path(
  figure_dir,
  "146_Supplementary_Figure_Sx_SRS_CTS_robustness.tiff"
)


ggplot2::ggsave(
  filename = main_png,
  plot = main_figure,
  width = 15,
  height = 9.4,
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = main_pdf,
  plot = main_figure,
  width = 15,
  height = 9.4,
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = main_tiff,
  plot = main_figure,
  width = 15,
  height = 9.4,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 22. EXPORT INDIVIDUAL PANELS
# =============================================================================

individual_panels <- list(
  A_SRS_k_agreement = p_A,
  B_SRSq_k_robustness = p_B,
  C_BP_batch_sensitivity = p_C,
  D_CTS_seed_reproducibility = p_D,
  E_CTS_preprocessing = p_E,
  F_CTS_batch_distribution = p_F
)


for (nm in names(
  individual_panels
)) {
  
  panel_plot <- individual_panels[[nm]]
  
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "146_panel_",
        nm,
        ".png"
      )
    ),
    plot = panel_plot,
    width = 6,
    height = 5,
    dpi = 600,
    bg = "white"
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "146_panel_",
        nm,
        ".pdf"
      )
    ),
    plot = panel_plot,
    width = 6,
    height = 5,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}


# =============================================================================
# 23. SAVE FIGURE DATA
# =============================================================================

master_xlsx <- file.path(
  table_dir,
  "146_SRS_CTS_robustness_figure_data.xlsx"
)


openxlsx::write.xlsx(
  
  list(
    
    Run_info = data.frame(
      
      item = c(
        "script",
        "source_40b",
        "source_40c",
        "source_107d",
        "new_classification_performed",
        "new_random_forest_performed",
        "new_batch_correction_performed",
        "new_statistical_tests_performed"
      ),
      
      value = c(
        "146_build_SRS_CTS_robustness_supplementary_figure.R",
        file_40b,
        file_40c,
        file_107d,
        "NO",
        "NO",
        "NO",
        "NO"
      ),
      
      stringsAsFactors = FALSE
    ),
    
    SRS_40b_original_k =
      srs_40b_original,
    
    SRS_40c_same_k =
      srs_40c_samek,
    
    CTS_seed_stability =
      cts_primary_stability,
    
    CTS_preprocessing =
      cts_preprocessing_plot,
    
    CTS_by_batch =
      cts_batch_long,
    
    CTS_batch_association =
      cts_batch_assoc
    
  ),
  
  master_xlsx,
  
  overwrite = TRUE
)


# =============================================================================
# 24. FIGURE CAPTION — ENGLISH
# =============================================================================

caption_en <- paste0(
  
  "Supplementary Figure Sx. Robustness of blood transcriptomic endotype ",
  "assignments. ",
  
  "(A) Agreement of categorical Sepsis Response Signature (SRS) assignments ",
  "with the primary k=20 analysis across alternative mutual-nearest-neighbour ",
  "parameter values using the original TMM-normalized logCPM input. ",
  "The shaded region indicates the approximate 20–30% sample-size-based ",
  "k window. ",
  
  "(B) Spearman correlation of quantitative SRSq values with the primary ",
  "analysis across alternative k values. ",
  
  "(C) Categorical agreement between original and conservatively ",
  "batch-adjusted BP-only SRS analyses at matched k values. At k=20, ",
  "categorical agreement was ",
  round(
    k20_agreement,
    1
  ),
  "% and the SRSq Spearman correlation was ",
  round(
    k20_rho,
    3
  ),
  ". Batch removal was evaluated as a sensitivity analysis and did not ",
  "replace the primary classification. ",
  
  "(D) Sample-level CTS reproducibility across 50 repeated random-forest ",
  "refits using the primary expression input. Thirty-four of 35 patients ",
  "retained the same CTS assignment in all 50 runs, and the minimum ",
  "individual modal classification frequency was ",
  round(
    min_cts_stability,
    0
  ),
  "%. ",
  
  "(E) Modal CTS agreement with the original Script 107 assignments under ",
  "the primary preprocessing, BP-only variance stabilization, and exploratory ",
  "within-BP batch removal. BP-only variance stabilization preserved all ",
  "35 assignments, whereas batch removal yielded 29/35 concordant modal ",
  "assignments. ",
  
  "(F) Distribution of original CTS assignments across sequencing batches. ",
  "No statistically significant CTS-by-batch association was detected by ",
  "Monte Carlo Fisher exact testing (p=",
  format_p(
    batch_p
  ),
  "; Cramér's V=",
  round(
    batch_v,
    3
  ),
  ")."
)


writeLines(
  caption_en,
  file.path(
    text_dir,
    "146_Supplementary_Figure_Sx_caption_EN.txt"
  )
)


# =============================================================================
# 25. FIGURE CAPTION — RUSSIAN
# =============================================================================

caption_ru <- paste0(
  
  "Дополнительный рисунок Sx. Устойчивость транскриптомной классификации ",
  "эндотипов крови. ",
  
  "(A) Совпадение категориальных классов Sepsis Response Signature (SRS) ",
  "с основным анализом k=20 при альтернативных значениях параметра ",
  "mutual-nearest-neighbour на исходных TMM-нормализованных logCPM. ",
  "Затенённая область соответствует приблизительному диапазону k, равному ",
  "20–30% размера выборки. ",
  
  "(B) Корреляция Спирмена количественных значений SRSq с основным анализом ",
  "при различных значениях k. ",
  
  "(C) Совпадение SRS между исходным и консервативным batch-adjusted ",
  "BP-only анализом при одинаковых значениях k. При k=20 категориальное ",
  "совпадение составило ",
  round(
    k20_agreement,
    1
  ),
  "%, а корреляция SRSq — ρ=",
  round(
    k20_rho,
    3
  ),
  ". Удаление batch-эффекта использовалось только как sensitivity analysis. ",
  
  "(D) Воспроизводимость CTS для отдельных пациентов при 50 повторных ",
  "обучениях random-forest classifier на основном наборе экспрессии. ",
  "У 34 из 35 пациентов классификация была идентична во всех 50 запусках, ",
  "а минимальная индивидуальная modal stability составила ",
  round(
    min_cts_stability,
    0
  ),
  "%. ",
  
  "(E) Совпадение modal CTS с исходным Script 107 при основном preprocessing, ",
  "BP-only variance stabilization и exploratory within-BP batch removal. ",
  "BP-only VST сохранил исходную классификацию у всех 35 пациентов, тогда ",
  "как после удаления batch-эффекта совпадение составило 29/35. ",
  
  "(F) Распределение исходных CTS по sequencing batches. Статистически ",
  "значимой связи между CTS и batch по Monte Carlo Fisher exact test ",
  "не выявлено (p=",
  format_p(
    batch_p
  ),
  "; Cramér's V=",
  round(
    batch_v,
    3
  ),
  ")."
)


writeLines(
  caption_ru,
  file.path(
    text_dir,
    "146_Supplementary_Figure_Sx_caption_RU.txt"
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
    "146_Supplementary_Figure_Sx_caption_EN_RU.txt"
  )
)


# =============================================================================
# 26. MANUSCRIPT-READY NUMERICAL SUMMARY
# =============================================================================

summary_lines <- c(
  
  "SCRIPT 146 — SRS / CTS ROBUSTNESS FIGURE",
  "====================================================================",
  "",
  
  "SRS — Script 40b:",
  
  paste0(
    "Minimum categorical agreement across original-logCPM k sensitivity: ",
    round(
      min(
        srs_40b_original$SRS_agreement_percent,
        na.rm = TRUE
      ),
      2
    ),
    "%"
  ),
  
  paste0(
    "Minimum SRSq Spearman rho across original-logCPM k sensitivity: ",
    round(
      min(
        srs_40b_original$SRSq_spearman_vs_primary,
        na.rm = TRUE
      ),
      3
    )
  ),
  
  "",
  
  "SRS — Script 40c:",
  
  paste0(
    "BP-only original vs batch-adjusted at k=20: agreement ",
    round(
      k20_agreement,
      1
    ),
    "%; SRSq rho ",
    round(
      k20_rho,
      3
    )
  ),
  
  "",
  
  "CTS — Script 107d:",
  
  paste0(
    "Primary CTS samples with 100% stability: ",
    n_cts_100,
    "/35"
  ),
  
  paste0(
    "Primary CTS samples with >=95% stability: ",
    n_cts_ge95,
    "/35"
  ),
  
  paste0(
    "Minimum individual modal stability: ",
    round(
      min_cts_stability,
      1
    ),
    "%"
  ),
  
  paste0(
    "CTS x sequencing batch Fisher p: ",
    signif(
      batch_p,
      5
    )
  ),
  
  paste0(
    "CTS x sequencing batch Cramer's V: ",
    round(
      batch_v,
      3
    )
  ),
  
  "",
  
  "INTERPRETATION:",
  
  paste0(
    "Primary SRS and CTS assignments were highly robust to parameter choice, ",
    "cohort composition, and stochastic classifier refitting. Explicit ",
    "within-cohort batch removal produced moderate reclassification and was ",
    "therefore retained as a conservative sensitivity analysis rather than ",
    "used to redefine primary endotype assignments."
  )
)


summary_file <- file.path(
  text_dir,
  "146_SRS_CTS_robustness_summary.txt"
)


writeLines(
  summary_lines,
  summary_file
)


# =============================================================================
# 27. SAVE FINAL NUMERICAL AUDIT
# =============================================================================

final_audit <- data.frame(
  
  item = c(
    "minimum_SRS_agreement_across_original_k_percent",
    "minimum_SRSq_rho_across_original_k",
    "BP_only_batch_adjusted_k20_SRS_agreement_percent",
    "BP_only_batch_adjusted_k20_SRSq_rho",
    "CTS_100_percent_seed_stable_n",
    "CTS_ge95_percent_seed_stable_n",
    "CTS_minimum_modal_stability_percent",
    "CTS_batch_Fisher_p",
    "CTS_batch_Cramers_V"
  ),
  
  value = c(
    
    min(
      srs_40b_original$SRS_agreement_percent,
      na.rm = TRUE
    ),
    
    min(
      srs_40b_original$SRSq_spearman_vs_primary,
      na.rm = TRUE
    ),
    
    k20_agreement,
    
    k20_rho,
    
    n_cts_100,
    
    n_cts_ge95,
    
    min_cts_stability,
    
    batch_p,
    
    batch_v
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  final_audit,
  file.path(
    table_dir,
    "146_final_numerical_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 28. FINAL CONSOLE OUTPUT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 146 completed successfully.\n")
cat("====================================================================\n\n")


cat("SCRIPT 40b — ORIGINAL k SENSITIVITY:\n")

print(
  srs_40b_original %>%
    dplyr::select(
      k,
      SRS_agreement_percent,
      SRSq_spearman_vs_primary
    )
)


cat("\nSCRIPT 40c — BP-ONLY SAME-k ORIGINAL vs BATCH-ADJUSTED:\n")

print(
  srs_40c_samek %>%
    dplyr::select(
      k,
      SRS_agreement_percent,
      SRSq_spearman
    )
)


cat("\nCTS PRIMARY SEED STABILITY:\n")

cat(
  "100% stable: ",
  n_cts_100,
  "/35\n",
  sep = ""
)

cat(
  ">=95% stable: ",
  n_cts_ge95,
  "/35\n",
  sep = ""
)

cat(
  "Minimum modal stability: ",
  min_cts_stability,
  "%\n",
  sep = ""
)


cat("\nCTS PREPROCESSING AGREEMENT:\n")

print(
  cts_preprocessing_plot %>%
    dplyr::select(
      preprocessing,
      agreement_n,
      n_total,
      agreement_percent
    )
)


cat("\nCTS x BATCH ASSOCIATION:\n")

cat(
  "Fisher p = ",
  signif(
    batch_p,
    6
  ),
  "\n",
  sep = ""
)

cat(
  "Cramer's V = ",
  round(
    batch_v,
    4
  ),
  "\n",
  sep = ""
)


cat("\nMAIN FIGURE:\n")

print(
  normalizePath(
    main_png,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nMASTER WORKBOOK:\n")

print(
  normalizePath(
    master_xlsx,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nNUMERICAL AUDIT:\n")

print(
  normalizePath(
    file.path(
      table_dir,
      "146_final_numerical_audit.csv"
    ),
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nCAPTION:\n")

print(
  normalizePath(
    file.path(
      text_dir,
      "146_Supplementary_Figure_Sx_caption_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nDone.\n")