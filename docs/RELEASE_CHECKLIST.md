# Release Checklist

## Before First Public GitHub Push

- Confirm repository name and visibility:
  `https://github.com/beloslav1976-svg/sepsis_blood_endotypes_biomarkers`
- MIT license terms selected for the repository release.
- Keep the working manuscript DOCX out of the public GitHub repository; the
  journal article will cite the repository after release.
- Use `github_repo_public_release_clean/sepsis_blood_endotypes_biomarkers` for
  the first public push, not earlier staging repositories.
- Decide whether staged R scripts should remain as provenance snapshots or be
  converted into a portable public analysis workflow.
- Test `SEPSIS_PROJECT_DIR` execution on a clean local clone before public
  release.
- Confirm final data availability wording for primary human sequencing data.
- Confirm whether raw count matrices can be deposited publicly, under controlled
  access, or only described by accession.
- Verify no direct identifiers are present in all staged spreadsheets.
- Decide whether earlier import/DESeq2 scripts before Script 135 should be
  included after privacy review.
- Root `.zenodo.json` prepared for `v0.1.0`; article DOI is not yet included.

## Zenodo Archival Status

- GitHub URL for manuscript:
  `https://github.com/beloslav1976-svg/sepsis_blood_endotypes_biomarkers`
- Zenodo DOI for manuscript:
  `https://doi.org/10.5281/zenodo.22228576`
- Concept DOI for all versions:
  `https://doi.org/10.5281/zenodo.22228575`
- GitHub release tag created: `v0.1.0`.
- Release archived on Zenodo.
- Save the Zenodo DOI and release tag in the submission files.

## Current Clean Package Commit Message

`Initial public release clean package`
