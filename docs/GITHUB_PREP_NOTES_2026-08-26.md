# GitHub Preparation Notes

Date: 2026-08-26

## Source Locations

- Manuscript and publication files:
  local manuscript/article working folder
- Frozen endotype/biomarker results:
  local `Sepsis_DESeq2` analysis results folder

## Actions Completed

- Created a working manuscript copy:
  `Статья Эндотипы Сепсис_CriticalCare_working.docx`
- Replaced the pending Ion Torrent software provenance paragraph with the
  verified Torrent Suite and ampliSeqRNA plugin versions supplied by the user.
- Removed the editorial placeholder requesting insertion of the software
  versions before submission.
- Staged main figures, main table, supplementary figures, supplementary tables,
  and freeze records from Scripts 167 and 168.
- Audited the referenced Table S2 workbook as data only. The workbook contains
  two visible sheets: `S2_ReadMe` and `Complete_DE`. The raw worksheet XML shows
  12,394 rows and 14 columns in `Complete_DE`, despite the workbook dimension
  metadata advertising `A1`.
- Added publication-facing R scripts 135-168 to preserve the build and freeze
  provenance for the blood endotype/biomarker package.
- Added a staged workbook identifier-label scan. The automated screen found
  label-level references to de-identified participant IDs, but no obvious direct
  identifier labels in the scanned rows.
- Removed the working manuscript DOCX from the GitHub staging repository. The
  manuscript file remains in the local article folder and is not intended for
  public repository release.
- Removed personal local paths from staged R-script copies and added
  `SEPSIS_PROJECT_DIR` as the public project-root setting.

## Release Cautions

- The staging package is not yet a public GitHub release.
- The manuscript still needs final journal-format cleanup beyond the methods
  paragraph fixed here.
- The working manuscript DOCX is intentionally not included in the public
  repository staging package.
- The final data accession statement, GitHub URL, Zenodo DOI, and license remain
  to be inserted after repository archival decisions are finalized.
