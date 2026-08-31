# Release Checklist

## Before First Public GitHub Push

- Confirm repository name and visibility.
- Confirm license terms for code, tables, and figure exports.
- Keep the working manuscript DOCX out of the public GitHub repository; the
  journal article will cite the repository after release.
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

## Before Zenodo Archival

- Insert final GitHub URL in the manuscript.
- Insert final Zenodo DOI in the manuscript.
- Regenerate `MANIFEST.csv` after the final file set is frozen.
- Create a GitHub release tag.
- Archive the release on Zenodo.
- Save the Zenodo DOI and final release tag in the submission files.

## Suggested Initial Commit Message

`Stage blood endotype biomarker publication package`
