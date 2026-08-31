# Running the Public Scripts

The publication-facing R scripts in `publication_scripts_135_168/` were adapted
for public release by removing personal local paths. They use the
`SEPSIS_PROJECT_DIR` environment variable to locate the analysis project root.

## Expected Project Layout

Set `SEPSIS_PROJECT_DIR` to a directory that contains the required analysis
inputs and results folders used by the scripts, for example:

```text
Sepsis_DESeq2/
  data/
  results/
    blood_endotypes_biomarkers/
```

If `SEPSIS_PROJECT_DIR` is not set, the scripts fall back to
`~/Sepsis_DESeq2`.

## Example

```r
Sys.setenv(SEPSIS_PROJECT_DIR = "/path/to/Sepsis_DESeq2")
source("analysis_scripts/publication_scripts_135_168/168_audit_and_freeze_complete_main_manuscript_package.R")
```

The staged repository includes manuscript-facing outputs and freeze records, but
does not include raw primary human sequencing data.
