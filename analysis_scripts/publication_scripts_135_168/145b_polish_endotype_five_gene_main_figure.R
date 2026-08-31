################################################################################
# Script 145b
# Final visual polishing of Figure 2
#
# Project: Sepsis_DESeq2
#
# Purpose:
#   Rebuild the main endotype/five-gene figure using finalized Script 145
#   data WITHOUT recalculating SRS, CTS, five-gene score, or statistics.
#
# Final visual changes:
#   - strict A-F panel labeling
#   - heatmap is a single panel E
#   - compact title for integrated CTS/SRS panel
#   - compact title for panel E: "Five-gene expression gradient"
#   - one descriptive linear trend in panel F across all 35 patients
#   - no CI band in panel F
#   - fill=SRS applied only to points in panel F to avoid ggplot warnings
#   - improved panel proportions
#   - no redundant legend in panel A
#   - publication-oriented heatmap labels
#
# IMPORTANT:
#   - NO statistical tests are recalculated.
#   - NO SRS classifications are recalculated.
#   - NO CTS classifications are recalculated.
#   - NO five-gene score is recalculated.
#   - All inferential statistics are inherited from finalized Script 145.
#   - The linear regression line in panel F is descriptive only.
#   - Statistical inference in panel F remains based on Spearman correlation.
#
# Input:
#   results/blood_endotypes_biomarkers/145_main_endotype_figure/
#     tables/145_endotype_five_gene_main_figure_data.xlsx
#     tables/145_main_figure_statistics.csv
#
# Output:
#   results/blood_endotypes_biomarkers/145b_polished_main_figure/
#
################################################################################


cat("====================================================================\n")
cat("Running Script 145b\n")
cat("Final polishing of Figure 2\n")
cat("====================================================================\n\n")


# =============================================================================
# 1. PROJECT
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
# 2. INPUT / OUTPUT
# =============================================================================

input_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "145_main_endotype_figure"
)

input_xlsx <- file.path(
  input_dir,
  "tables",
  "145_endotype_five_gene_main_figure_data.xlsx"
)

input_stats <- file.path(
  input_dir,
  "tables",
  "145_main_figure_statistics.csv"
)


out_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "145b_polished_main_figure"
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


required_files <- c(
  input_xlsx,
  input_stats
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "Missing required Script 145 outputs:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


cat("\nInput workbook:\n")

print(
  normalizePath(
    input_xlsx,
    winslash = "\\",
    mustWork = TRUE
  )
)


cat("\nInput statistics:\n")

print(
  normalizePath(
    input_stats,
    winslash = "\\",
    mustWork = TRUE
  )
)


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
    "145b_package_versions.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file = file.path(
    text_dir,
    "145b_sessionInfo.txt"
  )
)


# =============================================================================
# 4. COLORS
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
  low = "#2166AC",
  mid = "#F7F7F7",
  high = "#B2182B"
)


# =============================================================================
# 5. HELPERS
# =============================================================================

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
      
      legend.position = "top",
      
      legend.title = ggplot2::element_blank(),
      
      plot.margin = ggplot2::margin(
        6,
        6,
        6,
        6
      )
    )
}


# =============================================================================
# 6. READ FINALIZED SCRIPT 145 DATA
# =============================================================================

blood <- openxlsx::read.xlsx(
  input_xlsx,
  sheet = "Final_BP_data"
)

group_summary <- openxlsx::read.xlsx(
  input_xlsx,
  sheet = "Endotype_group_summary"
)

cross_tab <- openxlsx::read.xlsx(
  input_xlsx,
  sheet = "CTS_by_SRS"
)

heat_long <- openxlsx::read.xlsx(
  input_xlsx,
  sheet = "Five_gene_heatmap_z"
)

stats_summary <- readr::read_csv(
  input_stats,
  show_col_types = FALSE
) %>%
  as.data.frame()


cat("\nLoaded finalized data:\n")
cat("Blood rows: ", nrow(blood), "\n", sep = "")
cat("Heatmap rows: ", nrow(heat_long), "\n", sep = "")


# =============================================================================
# 7. BASIC INPUT QC
# =============================================================================

required_blood_columns <- c(
  "sample_id",
  "five_gene_score",
  "SRS",
  "SRSq",
  "CTS",
  "integrated_group"
)

missing_blood_columns <- setdiff(
  required_blood_columns,
  names(blood)
)

if (length(missing_blood_columns) > 0) {
  
  stop(
    "Missing required columns in Final_BP_data:\n",
    paste(
      missing_blood_columns,
      collapse = ", "
    )
  )
}


required_heat_columns <- c(
  "gene",
  "sample_id",
  "z_expression"
)

missing_heat_columns <- setdiff(
  required_heat_columns,
  names(heat_long)
)

if (length(missing_heat_columns) > 0) {
  
  stop(
    "Missing required columns in Five_gene_heatmap_z:\n",
    paste(
      missing_heat_columns,
      collapse = ", "
    )
  )
}


if (nrow(blood) != 35) {
  
  stop(
    "Expected 35 BP samples; found ",
    nrow(blood),
    "."
  )
}


# =============================================================================
# 8. RESTORE FACTOR LEVELS
# =============================================================================

blood <- blood %>%
  
  dplyr::mutate(
    
    sample_id = as.character(
      sample_id
    ),
    
    five_gene_score = as.numeric(
      five_gene_score
    ),
    
    SRSq = as.numeric(
      SRSq
    ),
    
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
    ),
    
    integrated_group = factor(
      integrated_group,
      levels = c(
        "CTS1/SRS1",
        "CTS2/SRS1",
        "CTS3/SRS1",
        "CTS3/SRS2"
      )
    )
  )


cross_tab <- cross_tab %>%
  
  dplyr::mutate(
    
    n = as.numeric(
      n
    ),
    
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
# 9. ENDOTYPE DISTRIBUTION QC
# =============================================================================

expected_srs <- c(
  SRS1 = 28,
  SRS2 = 7
)

expected_cts <- c(
  CTS1 = 14,
  CTS2 = 6,
  CTS3 = 15
)

expected_integrated <- c(
  "CTS1/SRS1" = 14,
  "CTS2/SRS1" = 6,
  "CTS3/SRS1" = 8,
  "CTS3/SRS2" = 7
)


observed_srs <- table(
  blood$SRS
)

observed_cts <- table(
  blood$CTS
)

observed_integrated <- table(
  blood$integrated_group
)


cat("\nSRS distribution:\n")
print(observed_srs)

cat("\nCTS distribution:\n")
print(observed_cts)

cat("\nIntegrated CTS/SRS distribution:\n")
print(observed_integrated)


if (!all(
  as.integer(
    observed_srs[
      names(expected_srs)
    ]
  ) ==
  as.integer(expected_srs)
)) {
  
  stop(
    "SRS distribution QC failed."
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
    "CTS distribution QC failed."
  )
}


if (!all(
  as.integer(
    observed_integrated[
      names(expected_integrated)
    ]
  ) ==
  as.integer(expected_integrated)
)) {
  
  stop(
    "Integrated CTS/SRS distribution QC failed."
  )
}


cat("\nEndotype distribution QC: PASS\n")


# =============================================================================
# 10. SAMPLE ORDER FOR HEATMAP
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


heat_long <- heat_long %>%
  
  dplyr::mutate(
    
    sample_id = factor(
      as.character(sample_id),
      levels = sample_order
    ),
    
    gene = factor(
      as.character(gene),
      levels = rev(
        c(
          "CD177",
          "HK3",
          "IRAK3",
          "CARD11",
          "IKZF2"
        )
      )
    ),
    
    z_expression = as.numeric(
      z_expression
    )
  )


# =============================================================================
# 11. EXTRACT FINAL STATISTICS
# =============================================================================

get_stat_row <- function(
    pattern
) {
  
  hit <- stats_summary[
    grepl(
      pattern,
      stats_summary$analysis,
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(hit) == 0) {
    
    stop(
      "Could not find statistic for pattern: ",
      pattern
    )
  }
  
  hit[1, , drop = FALSE]
}


stat_srs <- get_stat_row(
  "score by SRS"
)

stat_cts <- get_stat_row(
  "score by CTS"
)

stat_integrated <- get_stat_row(
  "integrated CTS"
)

stat_srsq <- get_stat_row(
  "score vs SRSq"
)


p_srs <- as.numeric(
  stat_srs$p_value
)

p_cts <- as.numeric(
  stat_cts$p_value
)

p_integrated <- as.numeric(
  stat_integrated$p_value
)

p_srsq <- as.numeric(
  stat_srsq$p_value
)


eps_cts <- as.numeric(
  stat_cts$effect_size
)

eps_integrated <- as.numeric(
  stat_integrated$effect_size
)

rho_srsq <- as.numeric(
  stat_srsq$effect_size
)


cat("\nFinalized Script 145 statistics retained:\n")

cat(
  "SRS p = ",
  signif(
    p_srs,
    6
  ),
  "\n",
  sep = ""
)

cat(
  "CTS p = ",
  signif(
    p_cts,
    6
  ),
  "; epsilon^2 = ",
  round(
    eps_cts,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Integrated p = ",
  signif(
    p_integrated,
    6
  ),
  "; epsilon^2 = ",
  round(
    eps_integrated,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "SRSq rho = ",
  round(
    rho_srsq,
    4
  ),
  "; p = ",
  signif(
    p_srsq,
    6
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 12. PANEL A
# SRS x CTS
# =============================================================================

p_A <- ggplot2::ggplot(
  cross_tab,
  ggplot2::aes(
    x = CTS,
    y = SRS,
    fill = n
  )
) +
  
  ggplot2::geom_tile(
    linewidth = 1.1,
    color = "white"
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        "n=",
        n
      )
    ),
    size = 4,
    fontface = "bold"
  ) +
  
  ggplot2::scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2166AC",
    guide = "none"
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::labs(
    tag = "A",
    
    title =
      "Hierarchical overlap of SRS and CTS",
    
    subtitle =
      "CTS1 and CTS2 are confined to SRS1, whereas CTS3 spans SRS1 and SRS2",
    
    x =
      "Consensus Transcriptomic Subtype",
    
    y =
      "Sepsis Response Signature"
  )


# =============================================================================
# 13. PANEL B
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
    alpha = 0.23,
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
    width = 0.09,
    size = 2,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_srs
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    legend.position = "none"
  ) +
  
  ggplot2::labs(
    tag = "B",
    
    title =
      "Five-gene score separates SRS states",
    
    subtitle =
      paste0(
        "Wilcoxon p=",
        format_p(
          p_srs
        )
      ),
    
    x = NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 14. PANEL C
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
    alpha = 0.23,
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
    width = 0.09,
    size = 2,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_cts
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    legend.position = "none"
  ) +
  
  ggplot2::labs(
    tag = "C",
    
    title =
      "Five-gene score spans CTS classes",
    
    subtitle =
      paste0(
        "Kruskal-Wallis p=",
        format_p(
          p_cts
        ),
        "; ε²=",
        round(
          eps_cts,
          3
        )
      ),
    
    x = NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 15. PANEL D
# INTEGRATED CTS/SRS
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
    width = 0.82,
    alpha = 0.20,
    trim = FALSE,
    color = NA
  ) +
  
  ggplot2::geom_boxplot(
    width = 0.40,
    outlier.shape = NA,
    alpha = 0.74,
    linewidth = 0.45
  ) +
  
  ggplot2::geom_jitter(
    width = 0.075,
    size = 2,
    alpha = 0.85
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_integrated
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    legend.position = "none",
    
    axis.text.x = ggplot2::element_text(
      angle = 15,
      hjust = 1,
      size = 9
    )
  ) +
  
  ggplot2::labs(
    tag = "D",
    
    title =
      "Integrated CTS/SRS states form a continuum",
    
    subtitle =
      paste0(
        "Kruskal-Wallis p=",
        format_p(
          p_integrated
        ),
        "; ε²=",
        round(
          eps_integrated,
          3
        )
      ),
    
    x = NULL,
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 16. PANEL E
# FIVE-GENE HEATMAP
# =============================================================================

annotation_long <- blood %>%
  
  dplyr::filter(
    sample_id %in%
      sample_order
  ) %>%
  
  dplyr::select(
    sample_id,
    SRS,
    CTS
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
      as.character(sample_id),
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


p_E_ann <- ggplot2::ggplot(
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
  
  ggplot2::theme_minimal(
    base_size = 8
  ) +
  
  ggplot2::theme(
    
    axis.title =
      ggplot2::element_blank(),
    
    axis.text.x =
      ggplot2::element_blank(),
    
    axis.text.y =
      ggplot2::element_text(
        size = 8
      ),
    
    axis.ticks =
      ggplot2::element_blank(),
    
    panel.grid =
      ggplot2::element_blank(),
    
    legend.position =
      "top",
    
    legend.key.size =
      grid::unit(
        0.28,
        "cm"
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 7.5
      ),
    
    plot.margin =
      ggplot2::margin(
        0,
        4,
        0,
        4
      )
  )


p_E_heat <- ggplot2::ggplot(
  heat_long,
  ggplot2::aes(
    x = sample_id,
    y = gene,
    fill = z_expression
  )
) +
  
  ggplot2::geom_tile() +
  
  ggplot2::scale_fill_gradient2(
    
    low =
      col_expression["low"],
    
    mid =
      col_expression["mid"],
    
    high =
      col_expression["high"],
    
    midpoint = 0,
    
    limits = c(
      -2.5,
      2.5
    ),
    
    oob =
      scales::squish
  ) +
  
  ggplot2::theme_classic(
    base_size = 9
  ) +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_blank(),
    
    axis.ticks.x =
      ggplot2::element_blank(),
    
    axis.text.y =
      ggplot2::element_text(
        size = 9
      ),
    
    axis.title =
      ggplot2::element_blank(),
    
    legend.position =
      "right",
    
    legend.title =
      ggplot2::element_text(
        size = 8
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 7
      ),
    
    plot.margin =
      ggplot2::margin(
        0,
        4,
        4,
        4
      )
  ) +
  
  ggplot2::labs(
    fill =
      "Gene-wise\nz-score"
  )


p_E_core <- (
  p_E_ann /
    p_E_heat
) +
  
  patchwork::plot_layout(
    heights = c(
      0.22,
      1
    ),
    guides = "collect"
  )


# -----------------------------------------------------------------------------
# Dedicated compact panel-E header
# -----------------------------------------------------------------------------

p_E_header <- ggplot2::ggplot() +
  
  ggplot2::annotate(
    "text",
    x = 0,
    y = 1,
    label = "E",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 5
  ) +
  
  ggplot2::annotate(
    "text",
    x = 0.08,
    y = 1,
    label = "Five-gene expression gradient",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  ggplot2::xlim(
    0,
    1
  ) +
  
  ggplot2::ylim(
    0,
    1
  ) +
  
  ggplot2::theme_void()


p_E <- (
  p_E_header /
    p_E_core
) +
  
  patchwork::plot_layout(
    heights = c(
      0.12,
      1
    )
  )


# =============================================================================
# 17. PANEL F
# FIVE-GENE SCORE vs SRSq
# =============================================================================

# IMPORTANT:
#
# Global ggplot mapping contains only x and y.
#
# fill = SRS is assigned ONLY to geom_point().
#
# Therefore geom_smooth() receives no fill grouping aesthetic and will not
# generate the "aesthetics were dropped" warning.
#
# group = 1 explicitly forces ONE descriptive regression line across all
# 35 patients.
#
# Inferential statistic remains Spearman rho inherited from Script 145.


p_F <- ggplot2::ggplot(
  blood,
  ggplot2::aes(
    x = SRSq,
    y = five_gene_score
  )
) +
  
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = 1
    ),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.75,
    color = "grey35"
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      fill = SRS
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.4,
    alpha = 0.90
  ) +
  
  ggplot2::scale_fill_manual(
    values = col_srs
  ) +
  
  theme_publication(
    10
  ) +
  
  ggplot2::theme(
    legend.position = "top"
  ) +
  
  ggplot2::labs(
    tag = "F",
    
    title =
      "Five-gene score tracks quantitative SRSq",
    
    subtitle =
      paste0(
        "Spearman ρ=",
        round(
          rho_srsq,
          3
        ),
        "; p=",
        format_p(
          p_srsq
        )
      ),
    
    x =
      "Quantitative SRS score (SRSq)",
    
    y =
      "Five-gene host-response score"
  )


# =============================================================================
# 18. FINAL LAYOUT
# =============================================================================

row1 <- (
  p_A |
    p_B |
    p_C
) +
  
  patchwork::plot_layout(
    widths = c(
      1.05,
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
      1.18,
      1.48,
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
      1.22
    )
  ) +
  
  patchwork::plot_annotation(
    
    title =
      "Blood transcriptomic endotypes converge on a five-gene host-response axis",
    
    subtitle =
      paste0(
        "Primary signature: CD177 + HK3 + IRAK3 − CARD11 − IKZF2; ",
        "n=35 patients with sepsis"
      ),
    
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
      )
    )
  )


# =============================================================================
# 19. EXPORT MAIN FIGURE
# =============================================================================

main_png <- file.path(
  figure_dir,
  "145b_Figure2_endotype_five_gene_host_response.png"
)

main_pdf <- file.path(
  figure_dir,
  "145b_Figure2_endotype_five_gene_host_response.pdf"
)

main_tiff <- file.path(
  figure_dir,
  "145b_Figure2_endotype_five_gene_host_response.tiff"
)


ggplot2::ggsave(
  filename = main_png,
  plot = main_figure,
  width = 15,
  height = 9.6,
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = main_pdf,
  plot = main_figure,
  width = 15,
  height = 9.6,
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = main_tiff,
  plot = main_figure,
  width = 15,
  height = 9.6,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 20. SAVE INDIVIDUAL POLISHED PANELS
# =============================================================================

individual_panels <- list(
  A_SRS_CTS_overlap = p_A,
  B_score_by_SRS = p_B,
  C_score_by_CTS = p_C,
  D_integrated_continuum = p_D,
  E_expression_gradient = p_E,
  F_score_vs_SRSq = p_F
)


for (nm in names(
  individual_panels
)) {
  
  panel_plot <- individual_panels[[nm]]
  
  if (nm == "E_expression_gradient") {
    
    panel_width <- 8
    panel_height <- 5
    
  } else {
    
    panel_width <- 6
    panel_height <- 5
  }
  
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "145b_panel_",
        nm,
        ".png"
      )
    ),
    plot = panel_plot,
    width = panel_width,
    height = panel_height,
    dpi = 600,
    bg = "white"
  )
  
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "145b_panel_",
        nm,
        ".pdf"
      )
    ),
    plot = panel_plot,
    width = panel_width,
    height = panel_height,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}


# =============================================================================
# 21. FINAL FIGURE CAPTION
# =============================================================================

caption_en <- paste0(
  
  "Figure 2. Blood transcriptomic endotypes converge on a five-gene ",
  "host-response axis. ",
  
  "(A) Cross-classification of 35 sepsis blood transcriptomes by ",
  "Sepsis Response Signature (SRS) and Consensus Transcriptomic Subtype ",
  "(CTS). CTS1 and CTS2 were confined to SRS1, whereas CTS3 included both ",
  "SRS1 and SRS2 samples. ",
  
  "(B) Five-gene host-response score across SRS classes. ",
  
  "(C) Five-gene host-response score across CTS classes. ",
  
  "(D) Integrated CTS/SRS states form an ordered molecular continuum from ",
  "CTS1/SRS1 through CTS2/SRS1 and CTS3/SRS1 to CTS3/SRS2. ",
  
  "(E) Gene-wise standardized expression of CD177, HK3, IRAK3, CARD11, ",
  "and IKZF2, with patients ordered from highest to lowest composite score; ",
  "annotation tracks indicate SRS and CTS assignments. ",
  
  "(F) Association between the five-gene host-response score and quantitative ",
  "SRSq. The regression line is shown for descriptive visualization only, ",
  "whereas statistical inference is based on Spearman rank correlation. ",
  
  "The five-gene signature was developed independently of SRS and CTS ",
  "assignments."
)


caption_ru <- paste0(
  
  "Рисунок 2. Транскриптомные эндотипы крови сходятся на общей ",
  "пятигенной оси ответа хозяина. ",
  
  "(A) Совместная классификация 35 транскриптомов крови пациентов с сепсисом ",
  "по Sepsis Response Signature (SRS) и Consensus Transcriptomic Subtype ",
  "(CTS). CTS1 и CTS2 встречались только внутри SRS1, тогда как CTS3 ",
  "включал образцы как SRS1, так и SRS2. ",
  
  "(B) Пятигенный показатель ответа хозяина в классах SRS. ",
  
  "(C) Пятигенный показатель в классах CTS. ",
  
  "(D) Интегрированные состояния CTS/SRS формируют упорядоченный ",
  "молекулярный континуум CTS1/SRS1 → CTS2/SRS1 → CTS3/SRS1 → CTS3/SRS2. ",
  
  "(E) Стандартизированная по генам экспрессия CD177, HK3, IRAK3, CARD11 ",
  "и IKZF2; пациенты расположены от максимального к минимальному значению ",
  "композитного пятигенного показателя. Верхние полосы показывают классы ",
  "SRS и CTS. ",
  
  "(F) Связь пятигенного показателя с количественным SRSq. Прямая ",
  "регрессии приведена только для визуализации тенденции; статистический ",
  "вывод основан на ранговой корреляции Спирмена. ",
  
  "SRS и CTS не использовались при отборе генов для пятигенной сигнатуры."
)


writeLines(
  caption_en,
  file.path(
    text_dir,
    "145b_Figure2_caption_EN.txt"
  )
)


writeLines(
  caption_ru,
  file.path(
    text_dir,
    "145b_Figure2_caption_RU.txt"
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
    "145b_Figure2_caption_EN_RU.txt"
  )
)


# =============================================================================
# 22. SAVE POLISHING PROVENANCE
# =============================================================================

polish_info <- data.frame(
  
  item = c(
    "script",
    "source_workbook",
    "source_statistics",
    "statistics_recalculated",
    "SRS_reassigned",
    "CTS_reassigned",
    "five_gene_score_recalculated",
    "panel_E_title",
    "panel_F_global_mapping",
    "panel_F_visual_trend",
    "panel_F_inference"
  ),
  
  value = c(
    "145b_polish_endotype_five_gene_main_figure.R",
    input_xlsx,
    input_stats,
    "NO",
    "NO",
    "NO",
    "NO",
    "Five-gene expression gradient",
    "x=SRSq; y=five_gene_score; fill applied only to geom_point",
    "One linear regression line across all 35 patients; visualization only",
    "Spearman rank correlation inherited from Script 145"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  polish_info,
  file.path(
    table_dir,
    "145b_polishing_provenance.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 23. SAVE FINAL NUMERICAL AUDIT
# =============================================================================

final_audit <- data.frame(
  
  item = c(
    "n_BP",
    "SRS1_n",
    "SRS2_n",
    "CTS1_n",
    "CTS2_n",
    "CTS3_n",
    "CTS1_SRS1_n",
    "CTS2_SRS1_n",
    "CTS3_SRS1_n",
    "CTS3_SRS2_n",
    "SRS_p",
    "CTS_p",
    "CTS_epsilon_squared",
    "Integrated_p",
    "Integrated_epsilon_squared",
    "SRSq_Spearman_rho",
    "SRSq_p"
  ),
  
  value = c(
    nrow(blood),
    as.integer(observed_srs["SRS1"]),
    as.integer(observed_srs["SRS2"]),
    as.integer(observed_cts["CTS1"]),
    as.integer(observed_cts["CTS2"]),
    as.integer(observed_cts["CTS3"]),
    as.integer(observed_integrated["CTS1/SRS1"]),
    as.integer(observed_integrated["CTS2/SRS1"]),
    as.integer(observed_integrated["CTS3/SRS1"]),
    as.integer(observed_integrated["CTS3/SRS2"]),
    p_srs,
    p_cts,
    eps_cts,
    p_integrated,
    eps_integrated,
    rho_srsq,
    p_srsq
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  final_audit,
  file.path(
    table_dir,
    "145b_final_figure_numerical_audit.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 24. FINAL CONSOLE OUTPUT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 145b completed successfully.\n")
cat("====================================================================\n\n")


cat("SOURCE DATA:\n")

print(
  normalizePath(
    input_xlsx,
    winslash = "\\",
    mustWork = TRUE
  )
)


cat("\nENDOTYPE DISTRIBUTIONS RETAINED:\n")

cat("SRS:\n")
print(observed_srs)

cat("\nCTS:\n")
print(observed_cts)

cat("\nIntegrated CTS/SRS:\n")
print(observed_integrated)


cat("\nFINALIZED STATISTICS RETAINED:\n")

cat(
  "SRS p = ",
  signif(
    p_srs,
    6
  ),
  "\n",
  sep = ""
)

cat(
  "CTS p = ",
  signif(
    p_cts,
    6
  ),
  "; epsilon^2 = ",
  round(
    eps_cts,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Integrated p = ",
  signif(
    p_integrated,
    6
  ),
  "; epsilon^2 = ",
  round(
    eps_integrated,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "SRSq rho = ",
  round(
    rho_srsq,
    4
  ),
  "; p = ",
  signif(
    p_srsq,
    6
  ),
  "\n",
  sep = ""
)


cat("\nPANEL F:\n")
cat(
  "One descriptive regression line across all 35 patients.\n"
)

cat(
  "Spearman correlation retained for statistical inference.\n"
)


cat("\nMAIN PNG:\n")

print(
  normalizePath(
    main_png,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nMAIN PDF:\n")

print(
  normalizePath(
    main_pdf,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nMAIN TIFF:\n")

print(
  normalizePath(
    main_tiff,
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nCAPTION:\n")

print(
  normalizePath(
    file.path(
      text_dir,
      "145b_Figure2_caption_EN.txt"
    ),
    winslash = "\\",
    mustWork = FALSE
  )
)


cat("\nDone.\n")