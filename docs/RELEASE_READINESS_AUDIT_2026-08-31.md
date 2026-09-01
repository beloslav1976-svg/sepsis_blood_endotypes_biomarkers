# Release Readiness Audit

Date: 2026-09-01

## Current Repository State

- Local repository branch: `main`
- Clean public-history source: `github_repo_public_release_clean` repository
  initialized from the sanitized public-ready file tree.
- Remote repository:
  `https://github.com/beloslav1976-svg/sepsis_blood_endotypes_biomarkers`
- Manifest entries after Zenodo DOI update: 135 files
- Manuscript DOCX: intentionally excluded from the public-ready repository
- Files larger than 100 MiB: none detected in the current working tree
- Absolute local path prefixes in freeze/provenance artifacts: replaced with
  `SEPSIS_PROJECT_DIR` placeholders for public release packaging

## Release Status

Current status: **published public GitHub repository with Zenodo archival DOI
issued for `v0.1.0`**.

The current repository is suitable as the public GitHub and Zenodo software
archive base. The article DOI is not yet known and is intentionally not included
in the current Zenodo metadata.

Zenodo version DOI: `10.5281/zenodo.22228576`

Zenodo concept DOI: `10.5281/zenodo.22228575`

The public release should be created from `github_repo_public_release_clean`,
not from any intermediate staging history that contained local provenance paths.

## Remaining Blockers Before Public Release

1. Final accession-level data availability wording for primary human sequencing
   data is still pending.
2. The article DOI should be added later when known.
3. The staged R scripts use `SEPSIS_PROJECT_DIR` instead of personal local
   paths, but rerunning them still requires the expected input data layout.

## Privacy and Workbook Review

A deep workbook privacy screen covered all 13 staged `.xlsx` files and 228,246
populated text cells. It detected no email, phone-like, or explicit personal-code
findings. The remaining seven automated findings were:

- two freeze-manifest timestamps;
- five participant-identifier/de-identification documentation labels in Table 1
  and Supplementary Table S1.

No direct-contact patterns were detected by the automated scan. This remains an
automated screen, not a substitute for the final study-team privacy and
data-use-agreement review.

See `docs/DEEP_WORKBOOK_PRIVACY_SCAN_2026-08-31.md` for the workbook-level
details.

## Portability Review

The staged publication-facing R scripts preserve the analysis provenance from
Scripts 135-168. Personal local paths were removed from the staged copies. The
scripts now use the `SEPSIS_PROJECT_DIR` environment variable, with
`~/Sepsis_DESeq2` as a fallback.

After Rscript became available locally, all 39 staged `.R` files parsed
successfully with zero parse errors.

Residual portability requirement: public users will need the expected
`Sepsis_DESeq2` input/results layout and access to any controlled primary human
sequencing data that are not included in the repository.

See `docs/R_SCRIPT_PORTABILITY_AUDIT_2026-08-31.md` for the script-level
portability audit.

Three staged R-script files contain double-dot names:

- `143_multicohort_five_gene_evidence_integration..R`
- `162_build_FigureS6_signature_convergence..R`
- `165_build_TableS10_GSE185263_external_evaluation..R`

These should be treated as provenance artifacts unless a clean public script
set is prepared.

## Recommended Public Release Path

1. Keep the manuscript DOCX outside the public repository; the journal article
   should cite the released repository and Zenodo DOI.
2. Confirm data availability wording.
3. Add article DOI when available.
4. Create a new GitHub/Zenodo version if any repository content changes are
   needed after submission.

See `docs/GITHUB_PUBLISH_RUNBOOK.md` for the first-push command sequence and
pre-push checks.
See `docs/ZENODO_RELEASE_RUNBOOK.md` for Zenodo archival steps.
