# Sepsis Blood Endotypes and Biomarkers

Publication-facing staging package for the manuscript on blood transcriptomic
endotypes and biology-guided host-response biomarkers in sepsis.

## Contents

- `analysis_scripts/publication_scripts_135_168/` - publication-facing R
  scripts used to build, audit, and freeze the blood endotype/biomarker package.
- `publication_package/main_figures_tables/` - main Table 1 and Figures 1-5 in
  publication formats.
- `publication_package/supplementary_files/` - Supplementary Tables S1-S10 and
  Supplementary Figures S1-S8 in publication formats.
- `freeze_records/supplementary_package_167/` - frozen supplementary package
  audit, provenance, freeze lock, and MD5 manifest.
- `freeze_records/main_package_168/` - frozen main manuscript package audit,
  numerical guardrails, provenance, freeze lock, and MD5 manifest.
- `docs/` - repository preparation notes, methods text, release-readiness
  audits, and GitHub/Zenodo citation templates.

## Frozen Analytical Scope

This repository staging package is blood-only. Urine transcriptomic analyses,
lncRNA-focused analyses, and paraspeckle analyses are outside the scope of this
manuscript package.

The primary study-derived five-gene panel is:

`CD177`, `HK3`, `IRAK3`, `CARD11`, and `IKZF2`.

The five-gene score should be described as a biology-guided molecular
host-response readout, not as a clinically calibrated diagnostic assay.

## Provenance

The publication-facing freeze records are defined by:

- Script 167: complete supplementary package freeze.
- Script 168: complete main manuscript package freeze.

The freeze records include MD5 manifests, audit summaries, source-script
provenance files, cross-reference maps, and reporting guardrails.

The staged R-script range is 135-168. This range covers the blood-only
endotype/biomarker validation, clinical associations, external evaluations,
publication figures/tables, and final package freeze audits. Earlier raw-import
and primary DESeq2 scripts are not staged here pending final public-data and
privacy decisions.

For rerunning scripts, set `SEPSIS_PROJECT_DIR` to the local analysis project
root. See `analysis_scripts/RUNNING_PUBLIC_SCRIPTS.md` and
`config/project_paths.example.R`.

## Current Release Status

This repository has been published on GitHub:

`https://github.com/beloslav1976-svg/sepsis_blood_endotypes_biomarkers`

It is a public repository package pending final release metadata. The manuscript
DOCX is intentionally excluded from this repository; the peer-reviewed article
should cite the final GitHub repository and Zenodo DOI after release. See
`docs/RELEASE_READINESS_AUDIT_2026-08-31.md` before tagging or archiving.

Publication metadata templates are provided in `docs/CITATION.cff.template`,
`docs/ZENODO_METADATA.template.json`, and `docs/GITHUB_RELEASE_DRAFT.md`.
Root `CITATION.cff`, `.zenodo.json`, and `LICENSE` have been prepared for the
initial Zenodo archival release. The article DOI is not yet available and is not
included in the current Zenodo metadata.

First-push instructions are provided in `docs/GITHUB_PUBLISH_RUNBOOK.md`.
Zenodo archival instructions are provided in `docs/ZENODO_RELEASE_RUNBOOK.md`.

## Data Availability Status

Primary human sequencing data accession details and final Zenodo DOI are still
pending final verification. Do not treat this repository as the final archived
release until the final accession statement and Zenodo DOI are inserted.

## Repository Release Checklist

- Verify that no direct identifiers are present in public tables.
- Replace manuscript placeholders for GitHub URL and Zenodo DOI.
- Add a final data availability statement aligned with the journal submission.
- MIT license file added for the repository release.
- Finalize the Zenodo DOI in documentation after archival.
- Follow `docs/GITHUB_PUBLISH_RUNBOOK.md` for the first GitHub push.
- Follow `docs/ZENODO_RELEASE_RUNBOOK.md` for Zenodo archival.
- Tag the final release and archive it with Zenodo.
