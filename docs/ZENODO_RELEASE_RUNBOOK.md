# Zenodo Release Runbook

This repository is published at:

`https://github.com/beloslav1976-svg/sepsis_blood_endotypes_biomarkers`

Use this runbook to archive a tagged GitHub release with Zenodo and obtain a
DOI for the manuscript data availability statement.

## Current Status

Zenodo archival is not yet complete.

Do not create the final GitHub release until the root `.zenodo.json`, license,
and final author metadata are confirmed.

## Required Before Release

- Confirm final creator list in Zenodo format: `Family, Given`.
- Confirm creator affiliations and ORCID identifiers, if available.
- Confirm repository license terms for code, tables, figures, and freeze
  records.
- Add the selected `LICENSE` file.
- Copy `docs/ZENODO_METADATA.template.json` to `.zenodo.json` and replace all
  placeholders.
- Optionally copy `docs/CITATION.cff.template` to `CITATION.cff`; when both
  files exist, Zenodo uses `.zenodo.json` for GitHub archival metadata.
- Regenerate `MANIFEST.csv` after final file changes.

## Enable GitHub Integration In Zenodo

1. Sign in to Zenodo.
2. Open the profile menu and select GitHub.
3. Click `Sync now`.
4. Find `beloslav1976-svg/sepsis_blood_endotypes_biomarkers`.
5. Enable the repository toggle.

After this is enabled, new GitHub releases are ingested by Zenodo.

## Create The GitHub Release

Create the release only after `.zenodo.json`, `LICENSE`, and final metadata are
committed and pushed.

```powershell
git tag -a v0.1.0 -m "Sepsis blood endotypes and biomarkers publication package"
git push origin v0.1.0
gh release create v0.1.0 --title "Sepsis blood endotypes and biomarkers publication package" --notes-file docs/GITHUB_RELEASE_DRAFT.md
```

## After Zenodo Processing

1. Wait for Zenodo to finish processing the GitHub release.
2. Open the created Zenodo record and copy the DOI.
3. Update:
   - `README.md`;
   - `docs/DATA_AVAILABILITY_DRAFT.md`;
   - `docs/RELEASE_READINESS_AUDIT_2026-08-31.md`;
   - manuscript data availability statement.
4. Commit and push the DOI update.

## Notes

- Zenodo metadata should be final before archival.
- If the Zenodo release fails, check the Zenodo release error panel first; most
  failures at this stage are metadata-format or license-related.
- Do not archive earlier local staging repositories.
