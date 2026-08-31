# Release Readiness Audit

Date: 2026-08-31

## Current Repository State

- Local repository branch: `main`
- Clean public-history source: sanitized public-ready staging tree derived from
  base commit `ea138ee`
- Remote repository: not configured
- Manifest entries after metadata/template update: 130 files
- Manuscript DOCX: intentionally excluded from the public-ready repository
- Files larger than 100 MiB: none detected in the current working tree
- Absolute local path prefixes in freeze/provenance artifacts: replaced with
  `SEPSIS_PROJECT_DIR` placeholders for public release packaging

## Release Status

Current status: **public-ready staging repository, pending final release
metadata**.

The current repository is suitable as a clean local base for a public GitHub
repository after the final license, data availability, citation, and Zenodo
metadata are confirmed. It should not yet be treated as the final public release
because the final GitHub URL, Zenodo DOI, article DOI, and license are still
placeholders.

The public release should be created from the sanitized file tree, not from any
intermediate history that contained local provenance paths.

## Remaining Blockers Before Public Release

1. Final GitHub URL and Zenodo DOI are not yet available.
2. License terms for code, tables, and figure exports are not yet confirmed.
3. Final accession-level data availability wording for primary human sequencing
   data is still pending.
4. `docs/CITATION.cff.template` and `docs/ZENODO_METADATA.template.json` must
   be finalized after author, repository, DOI, date, and license metadata are
   confirmed.
5. The staged R scripts use `SEPSIS_PROJECT_DIR` instead of personal local
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
2. Confirm license and data availability wording.
3. Finalize `CITATION.cff` and `.zenodo.json` from the provided templates.
4. Regenerate `MANIFEST.csv` after the final public file set is chosen.
5. Create the public GitHub repository, push `main`, create a release tag, and
   archive the tagged release with Zenodo.
