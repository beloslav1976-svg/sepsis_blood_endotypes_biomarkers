# R Script Portability Audit

Date: 2026-08-31

## Scope

Staged R scripts in:

`analysis_scripts/publication_scripts_135_168/`

## Actions

- Replaced personal local project-root paths with `SEPSIS_PROJECT_DIR`.
- Removed staged script lookups that depended on the user's Desktop.
- Added `analysis_scripts/RUNNING_PUBLIC_SCRIPTS.md`.
- Added `config/project_paths.example.R`.
- Added `config/project_paths.R` to `.gitignore` so private local paths can be
  used without being committed.

## Verification

- Personal Windows user-path matches: 0
- Personal POSIX-style user-path matches: 0
- Desktop-path matches in staged R scripts: 0
- DOCX files in repository staging area: 0
- Files larger than 100 MiB: 0
- R parse check: 39 staged `.R` files parsed; parse errors: 0

## Remaining Limitations

- Rscript was located at `C:\Program Files\R\R-4.6.0\bin\Rscript.exe` and used
  for a non-executing syntax parse of the staged R scripts.
- R emitted startup locale warnings and text-encoding warnings for symbols such
  as `≤`, `≥`, and en dash in character strings. These did not produce parse
  errors.
- Rerunning the scripts requires the expected `Sepsis_DESeq2` input/results
  directory layout.
- Raw or controlled primary human sequencing data are not included in this
  repository staging package.
- Three double-dot script names remain as provenance artifacts unless cleaned
  after final author review.
