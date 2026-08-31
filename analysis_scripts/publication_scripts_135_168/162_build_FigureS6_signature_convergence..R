################################################################################
# Script 162
# FINAL v2
#
# Supplementary Figure S6
#
# Convergence of study-derived and published transcriptomic signatures
# on SRSq and CTS endotype structure
#
# Project:
#   Sepsis_DESeq2
#
#
# PURPOSE
# -------
#
# Generate Supplementary Figure S6 exclusively from frozen
# Supplementary Table S8.
#
#
# PANELS
# ------
#
# A. Spearman correlation with continuous SRSq
#
# B. Kruskal-Wallis epsilon-squared across CTS classes
#
#
# IMPORTANT
# ---------
#
# No statistical analysis is rerun.
#
# The figure demonstrates cross-signature biological convergence,
# not superiority of one signature over another.
#
################################################################################


cat("====================================================================\n")
cat("Running Script 162 FINAL v2\n")
cat("Supplementary Figure S6\n")
cat("Transcriptomic-signature convergence\n")
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
# 2. PACKAGES
# =============================================================================

required_packages <- c(
  "dplyr",
  "ggplot2",
  "readxl",
  "openxlsx",
  "gridExtra"
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
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(gridExtra)
})


# =============================================================================
# 3. OUTPUT DIRECTORIES
# =============================================================================

output_dir <- file.path(
  "results",
  "blood_endotypes_biomarkers",
  "162_FigureS6_signature_convergence"
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


audit_dir <- file.path(
  output_dir,
  "audit"
)


for (
  one_dir in c(
    output_dir,
    figures_dir,
    tables_dir,
    text_dir,
    audit_dir
  )
) {
  
  dir.create(
    one_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# =============================================================================
# 4. INPUT — FROZEN TABLE S8
# =============================================================================

tableS8_file <- file.path(
  project_dir,
  "results",
  "blood_endotypes_biomarkers",
  "161_TableS8_published_signature_benchmarking",
  "tables",
  "161_TableS8_published_signature_benchmarking.xlsx"
)


if (!file.exists(tableS8_file)) {
  
  stop(
    paste0(
      "Frozen Supplementary Table S8 not found:\n",
      tableS8_file
    )
  )
}


benchmark <- readxl::read_excel(
  tableS8_file,
  sheet = "Benchmark_summary"
) %>%
  
  as.data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat("\nFrozen Table S8:\n")

cat(
  normalizePath(
    tableS8_file,
    winslash = "\\",
    mustWork = TRUE
  ),
  "\n"
)


cat(
  "Benchmark dimensions = ",
  nrow(benchmark),
  " x ",
  ncol(benchmark),
  "\n",
  sep = ""
)


# =============================================================================
# 5. REQUIRED COLUMN AUDIT
# =============================================================================

required_columns <- c(
  "Signature",
  "Signature_type",
  "n",
  "SRSq_Spearman_rho",
  "SRSq_nominal_P",
  "SRSq_BH_adjusted_P",
  "CTS_epsilon_squared",
  "CTS_nominal_P",
  "CTS_BH_adjusted_P",
  "SRSq_rank_descriptive",
  "CTS_rank_descriptive"
)


missing_columns <- setdiff(
  required_columns,
  names(benchmark)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "Missing required Table S8 columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


if (nrow(benchmark) != 7) {
  
  stop(
    paste0(
      "Expected seven frozen signatures; observed ",
      nrow(benchmark),
      "."
    )
  )
}


# =============================================================================
# 6. NUMERIC STANDARDIZATION
# =============================================================================

benchmark <- benchmark %>%
  
  dplyr::mutate(
    
    n =
      suppressWarnings(
        as.numeric(n)
      ),
    
    SRSq_Spearman_rho =
      suppressWarnings(
        as.numeric(SRSq_Spearman_rho)
      ),
    
    SRSq_nominal_P =
      suppressWarnings(
        as.numeric(SRSq_nominal_P)
      ),
    
    SRSq_BH_adjusted_P =
      suppressWarnings(
        as.numeric(SRSq_BH_adjusted_P)
      ),
    
    CTS_epsilon_squared =
      suppressWarnings(
        as.numeric(CTS_epsilon_squared)
      ),
    
    CTS_nominal_P =
      suppressWarnings(
        as.numeric(CTS_nominal_P)
      ),
    
    CTS_BH_adjusted_P =
      suppressWarnings(
        as.numeric(CTS_BH_adjusted_P)
      ),
    
    SRSq_rank_descriptive =
      suppressWarnings(
        as.numeric(SRSq_rank_descriptive)
      ),
    
    CTS_rank_descriptive =
      suppressWarnings(
        as.numeric(CTS_rank_descriptive)
      )
  )


# =============================================================================
# 7. FROZEN ANCHORS
# =============================================================================

expected_signatures <- c(
  "LIFTS-like",
  "DCAF17 five-gene alternative",
  "Primary five-gene",
  "FAIM3:PLAC8-related",
  "MetaScore-like",
  "RAPID-related (PLAC8-PLA2G7 contrast)",
  "SeptiCyte LAB-like"
)


anchor <- data.frame(
  
  Signature =
    expected_signatures,
  
  expected_rho = c(
    0.8518207,
    0.7943978,
    0.7649860,
    0.7162465,
    0.6848739,
    0.4946779,
    0.4624650
  ),
  
  expected_epsilon2 = c(
    0.6978685,
    0.6096457,
    0.6606774,
    0.6492474,
    0.5292361,
    0.4361125,
    0.3281250
  ),
  
  stringsAsFactors = FALSE
)


if (
  !setequal(
    benchmark$Signature,
    expected_signatures
  )
) {
  
  stop(
    "Frozen Table S8 signature set differs from expected set."
  )
}


anchor_audit <- anchor %>%
  
  dplyr::left_join(
    benchmark,
    by = "Signature"
  ) %>%
  
  dplyr::mutate(
    
    rho_difference =
      abs(
        SRSq_Spearman_rho -
          expected_rho
      ),
    
    epsilon2_difference =
      abs(
        CTS_epsilon_squared -
          expected_epsilon2
      ),
    
    rho_pass =
      rho_difference <
      1e-6,
    
    epsilon2_pass =
      epsilon2_difference <
      1e-6,
    
    overall_pass =
      rho_pass &
      epsilon2_pass
  )


if (
  !all(
    anchor_audit$overall_pass
  )
) {
  
  stop(
    "Frozen Figure S6 anchor audit failed."
  )
}


if (
  !all(
    benchmark$n ==
    35
  )
) {
  
  stop(
    "Expected n=35 for all seven benchmark signatures."
  )
}


if (
  sum(
    benchmark$SRSq_BH_adjusted_P <
    0.05
  ) !=
  7
) {
  
  stop(
    "Expected 7/7 SRSq associations with frozen BH <0.05."
  )
}


if (
  sum(
    benchmark$CTS_BH_adjusted_P <
    0.05
  ) !=
  7
) {
  
  stop(
    "Expected 7/7 CTS associations with frozen BH <0.05."
  )
}


cat(
  "\nFrozen numerical audit passed: 7/7 signatures.\n"
)


# =============================================================================
# 8. FIGURE GROUP
# =============================================================================

benchmark <- benchmark %>%
  
  dplyr::mutate(
    
    Figure_group =
      dplyr::case_when(
        
        Signature ==
          "Primary five-gene" ~
          "Primary five-gene",
        
        Signature ==
          "DCAF17 five-gene alternative" ~
          "Alternative five-gene",
        
        TRUE ~
          "Published comparator"
      )
  )


# =============================================================================
# 9. SHORT DISPLAY LABELS
# =============================================================================
#
# Only graphical labels are shortened.
# Frozen manuscript names remain unchanged in Table S8.
#
# =============================================================================

benchmark <- benchmark %>%
  
  dplyr::mutate(
    
    Signature_display =
      dplyr::case_when(
        
        Signature ==
          "RAPID-related (PLAC8-PLA2G7 contrast)" ~
          "RAPID-related\n(PLAC8-PLA2G7 contrast)",
        
        TRUE ~
          Signature
      )
  )


# =============================================================================
# 10. FIXED DISPLAY ORDER
# =============================================================================

signature_order <- benchmark %>%
  
  dplyr::arrange(
    SRSq_rank_descriptive
  ) %>%
  
  dplyr::pull(
    Signature_display
  )


benchmark$Signature_display <- factor(
  benchmark$Signature_display,
  levels = rev(
    signature_order
  )
)


cat("\nFIXED DISPLAY ORDER\n")
cat("-------------------\n")

print(
  signature_order
)


# =============================================================================
# 11. COLORS AND SHAPES
# =============================================================================

group_colors <- c(
  
  "Primary five-gene" =
    "#54278F",
  
  "Alternative five-gene" =
    "#E08214",
  
  "Published comparator" =
    "#2B8CBE"
)


group_shapes <- c(
  
  "Primary five-gene" =
    17,
  
  "Alternative five-gene" =
    15,
  
  "Published comparator" =
    16
)


# =============================================================================
# 12. COMMON THEME
# =============================================================================

benchmark_theme <- ggplot2::theme_bw(
  base_size = 11
) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 9.5
      ),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    panel.grid.major.y =
      ggplot2::element_blank(),
    
    axis.text.y =
      ggplot2::element_text(
        size = 9.3,
        lineheight = 0.95
      ),
    
    legend.position =
      "bottom",
    
    legend.title =
      ggplot2::element_blank(),
    
    plot.margin =
      ggplot2::margin(
        8,
        14,
        8,
        8
      )
  )


# =============================================================================
# 13. PANEL A — SRSq
# =============================================================================

panel_A <- ggplot2::ggplot(
  
  benchmark,
  
  ggplot2::aes(
    x = SRSq_Spearman_rho,
    y = Signature_display,
    color = Figure_group,
    shape = Figure_group
  )
  
) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend = SRSq_Spearman_rho,
      yend = Signature_display
    ),
    
    linewidth = 0.7,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_point(
    size = 4
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        sprintf(
          "%.3f",
          SRSq_Spearman_rho
        )
    ),
    
    hjust = -0.25,
    size = 3.4,
    color = "#333333",
    show.legend = FALSE
  ) +
  
  ggplot2::scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  
  ggplot2::scale_shape_manual(
    values = group_shapes,
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    
    limits = c(
      0,
      1.0
    ),
    
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  ggplot2::labs(
    
    title =
      "A  Alignment with continuous SRSq",
    
    subtitle =
      "Spearman correlation with SRSq; frozen BH-adjusted P < 0.05 for all 7",
    
    x =
      "Spearman rho with SRSq",
    
    y =
      NULL,
    
    color =
      NULL,
    
    shape =
      NULL
  ) +
  
  benchmark_theme


# =============================================================================
# 14. PANEL B — CTS
# =============================================================================

panel_B <- ggplot2::ggplot(
  
  benchmark,
  
  ggplot2::aes(
    x = CTS_epsilon_squared,
    y = Signature_display,
    color = Figure_group,
    shape = Figure_group
  )
  
) +
  
  ggplot2::geom_segment(
    
    ggplot2::aes(
      x = 0,
      xend = CTS_epsilon_squared,
      yend = Signature_display
    ),
    
    linewidth = 0.7,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_point(
    size = 4
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label =
        sprintf(
          "%.3f",
          CTS_epsilon_squared
        )
    ),
    
    hjust = -0.25,
    size = 3.4,
    color = "#333333",
    show.legend = FALSE
  ) +
  
  ggplot2::scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  
  ggplot2::scale_shape_manual(
    values = group_shapes,
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    
    limits = c(
      0,
      0.8
    ),
    
    breaks = seq(
      0,
      0.8,
      by = 0.2
    ),
    
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  ggplot2::labs(
    
    title =
      "B  Alignment with CTS structure",
    
    subtitle =
      "CTS epsilon-squared; frozen BH-adjusted P < 0.05 for all 7",
    
    x =
      "CTS epsilon-squared",
    
    y =
      NULL,
    
    color =
      NULL,
    
    shape =
      NULL
  ) +
  
  benchmark_theme


# =============================================================================
# 15. COMBINED FIGURE
# =============================================================================

figure_S6 <- gridExtra::arrangeGrob(
  
  panel_A,
  panel_B,
  
  ncol = 2,
  
  widths = c(
    1,
    1
  )
)


# =============================================================================
# 16. SAVE FIGURE
# =============================================================================

figure_pdf <- file.path(
  figures_dir,
  "162_FigureS6_signature_convergence.pdf"
)


figure_png <- file.path(
  figures_dir,
  "162_FigureS6_signature_convergence.png"
)


figure_tiff <- file.path(
  figures_dir,
  "162_FigureS6_signature_convergence.tiff"
)


ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure_S6,
  width = 15,
  height = 7.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_png,
  plot = figure_S6,
  width = 15,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure_S6,
  width = 15,
  height = 7.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# =============================================================================
# 17. SOURCE DATA
# =============================================================================

source_data_file <- file.path(
  tables_dir,
  "162_FigureS6_source_data.xlsx"
)


source_table <- benchmark %>%
  
  dplyr::mutate(
    Signature_display =
      as.character(
        Signature_display
      )
  ) %>%
  
  dplyr::select(
    Signature,
    Signature_display,
    Figure_group,
    n,
    SRSq_Spearman_rho,
    SRSq_nominal_P,
    SRSq_BH_adjusted_P,
    CTS_epsilon_squared,
    CTS_nominal_P,
    CTS_BH_adjusted_P,
    SRSq_rank_descriptive,
    CTS_rank_descriptive
  )


wb <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wb,
  "FigureS6_data"
)


openxlsx::writeData(
  wb,
  "FigureS6_data",
  source_table,
  withFilter = TRUE
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE
)


openxlsx::addStyle(
  wb,
  "FigureS6_data",
  header_style,
  rows = 1,
  cols = seq_len(
    ncol(source_table)
  ),
  gridExpand = TRUE
)


openxlsx::freezePane(
  wb,
  "FigureS6_data",
  firstActiveRow = 2
)


openxlsx::setColWidths(
  wb,
  "FigureS6_data",
  cols = seq_len(
    ncol(source_table)
  ),
  widths = "auto"
)


openxlsx::saveWorkbook(
  wb,
  source_data_file,
  overwrite = TRUE
)


# =============================================================================
# 18. FIGURE LEGEND
# =============================================================================

figure_legend <- paste0(
  
  "Supplementary Figure S6. Convergence of study-derived and published ",
  "blood transcriptomic signatures on SRS and CTS endotype structure. ",
  
  "(A) Spearman correlations between each frozen signature score and the ",
  "continuous SRSq output. ",
  
  "(B) Epsilon-squared effect sizes from frozen Kruskal-Wallis analyses ",
  "across CTS classes. ",
  
  "Signatures are shown in the same descriptive order in both panels, ",
  "defined by their frozen SRSq correlation. Purple and orange symbols ",
  "identify the primary and alternative study-derived five-gene signatures, ",
  "respectively; blue circles identify published-signature RNA-seq analogues. ",
  
  "All seven SRSq associations and all seven CTS associations retained ",
  "Benjamini-Hochberg-adjusted P values below 0.05 in the frozen Script 137 ",
  "benchmarking analysis. The LIFTS-like score showed the largest observed ",
  "alignment with both SRSq (rho=",
  sprintf(
    "%.3f",
    max(
      benchmark$SRSq_Spearman_rho
    )
  ),
  ") and CTS structure (epsilon-squared=",
  sprintf(
    "%.3f",
    max(
      benchmark$CTS_epsilon_squared
    )
  ),
  "). Differences in effect magnitude are descriptive and were not tested ",
  "as formal pairwise differences between signatures. Published comparators ",
  "represent the RNA-seq gene-based implementations used in this study; ",
  "SeptiCyte LAB-like does not represent the proprietary clinical assay ",
  "output. The figure illustrates cross-signature biological convergence ",
  "rather than diagnostic superiority or independent external validation."
)


legend_file <- file.path(
  text_dir,
  "162_FigureS6_legend_EN.txt"
)


writeLines(
  figure_legend,
  legend_file
)


# =============================================================================
# 19. RESULTS 3.7 TEXT
# =============================================================================

primary_row <- benchmark[
  benchmark$Signature ==
    "Primary five-gene",
  ,
  drop = FALSE
]


dcaf_row <- benchmark[
  benchmark$Signature ==
    "DCAF17 five-gene alternative",
  ,
  drop = FALSE
]


lifts_row <- benchmark[
  benchmark$Signature ==
    "LIFTS-like",
  ,
  drop = FALSE
]


faim_row <- benchmark[
  benchmark$Signature ==
    "FAIM3:PLAC8-related",
  ,
  drop = FALSE
]


results_3_7 <- paste0(
  
  "Published and study-derived blood transcriptomic signatures showed broad ",
  "convergence on the same molecular endotype axis (Supplementary Fig. S6 ",
  "and Supplementary Table S8). All seven benchmark implementations were ",
  "positively correlated with continuous SRSq and differed significantly ",
  "across CTS classes after the multiple-testing corrections retained from ",
  "the frozen benchmarking analysis. The largest observed SRSq correlation ",
  "was obtained for the LIFTS-like score (Spearman rho=",
  sprintf(
    "%.3f",
    lifts_row$SRSq_Spearman_rho
  ),
  "), followed by the DCAF17 alternative five-gene score (rho=",
  sprintf(
    "%.3f",
    dcaf_row$SRSq_Spearman_rho
  ),
  "), the primary five-gene score (rho=",
  sprintf(
    "%.3f",
    primary_row$SRSq_Spearman_rho
  ),
  "), and the FAIM3:PLAC8-related score (rho=",
  sprintf(
    "%.3f",
    faim_row$SRSq_Spearman_rho
  ),
  "). CTS-associated effect sizes showed a broadly similar pattern, with ",
  "epsilon-squared=",
  sprintf(
    "%.3f",
    lifts_row$CTS_epsilon_squared
  ),
  " for LIFTS-like, ",
  sprintf(
    "%.3f",
    primary_row$CTS_epsilon_squared
  ),
  " for the primary five-gene score, ",
  sprintf(
    "%.3f",
    faim_row$CTS_epsilon_squared
  ),
  " for FAIM3:PLAC8-related, and ",
  sprintf(
    "%.3f",
    dcaf_row$CTS_epsilon_squared
  ),
  " for the DCAF17 alternative. Thus, the primary five-gene signature did ",
  "not uniquely define the molecular continuum; several independently ",
  "developed transcriptomic signatures converged on the same underlying ",
  "host-response structure."
)


results_file <- file.path(
  text_dir,
  "162_proposed_Results_3.7_signature_convergence_EN.txt"
)


writeLines(
  results_3_7,
  results_file
)


# =============================================================================
# 20. AUDIT
# =============================================================================

audit_file <- file.path(
  audit_dir,
  "162_FigureS6_audit.xlsx"
)


audit_summary <- data.frame(
  
  Metric = c(
    "Number of signatures",
    "Samples per signature",
    "SRSq BH-significant signatures",
    "CTS BH-significant signatures",
    "Maximum SRSq rho",
    "Maximum CTS epsilon-squared"
  ),
  
  Value = c(
    nrow(benchmark),
    paste(
      unique(benchmark$n),
      collapse = ", "
    ),
    sum(
      benchmark$SRSq_BH_adjusted_P <
        0.05
    ),
    sum(
      benchmark$CTS_BH_adjusted_P <
        0.05
    ),
    max(
      benchmark$SRSq_Spearman_rho
    ),
    max(
      benchmark$CTS_epsilon_squared
    )
  ),
  
  stringsAsFactors = FALSE
)


wb_audit <- openxlsx::createWorkbook()


openxlsx::addWorksheet(
  wb_audit,
  "Audit_summary"
)


openxlsx::writeData(
  wb_audit,
  "Audit_summary",
  audit_summary
)


openxlsx::addWorksheet(
  wb_audit,
  "Frozen_anchor_audit"
)


openxlsx::writeData(
  wb_audit,
  "Frozen_anchor_audit",
  anchor_audit
)


openxlsx::saveWorkbook(
  wb_audit,
  audit_file,
  overwrite = TRUE
)


capture.output(
  sessionInfo(),
  file = file.path(
    audit_dir,
    "162_sessionInfo.txt"
  )
)


# =============================================================================
# 21. FINAL REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("Script 162 FINAL v2 completed successfully.\n")
cat("====================================================================\n\n")


cat("KEY CONVERGENCE RESULTS\n")
cat("-----------------------\n")


cat(
  "Highest SRSq rho: LIFTS-like = ",
  max(
    benchmark$SRSq_Spearman_rho
  ),
  "\n",
  sep = ""
)


cat(
  "Highest CTS epsilon-squared: LIFTS-like = ",
  max(
    benchmark$CTS_epsilon_squared
  ),
  "\n",
  sep = ""
)


cat(
  "SRSq BH-significant = ",
  sum(
    benchmark$SRSq_BH_adjusted_P <
      0.05
  ),
  "/7\n",
  sep = ""
)


cat(
  "CTS BH-significant = ",
  sum(
    benchmark$CTS_BH_adjusted_P <
      0.05
  ),
  "/7\n",
  sep = ""
)


cat("\nOUTPUT FILES\n")
cat("------------\n")


cat(
  "Figure S6 PDF:\n  ",
  normalizePath(
    figure_pdf,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S6 PNG:\n  ",
  normalizePath(
    figure_png,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)


cat(
  "Figure S6 TIFF:\n  ",
  normalizePath(
    figure_tiff,
    winslash = "\\",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)


cat("\nREPORTING GUARDRAILS\n")
cat("--------------------\n")


cat(
  "- Figure S6 uses frozen Table S8 only.\n"
)


cat(
  "- No statistical result is recalculated.\n"
)


cat(
  "- All seven SRSq and CTS associations retain frozen BH <0.05.\n"
)


cat(
  "- Display ordering is descriptive, not a formal performance ranking.\n"
)


cat(
  "- Effect-size differences are not superiority tests.\n"
)


cat(
  "- Benchmarking demonstrates cross-signature biological convergence.\n"
)


cat(
  "- This is not independent external validation.\n"
)


cat("\nDone.\n")