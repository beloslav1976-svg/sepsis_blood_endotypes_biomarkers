# GitHub Publish Runbook

Use this directory as the public repository source:

`github_repo_public_release_clean/sepsis_blood_endotypes_biomarkers`

Do not push earlier staging directories that may contain intermediate Git
history or local provenance paths.

## Pre-Push Checks

Run these checks from the repository root before the first public push:

```powershell
git status --short
git log --oneline --decorate -3
git remote -v
git grep -n --fixed-strings 'PRIVATE_USER_PATH_TOKEN' HEAD
git grep -n --fixed-strings 'PRIVATE_PROJECT_PATH_TOKEN' HEAD
```

Expected results:

- `git status --short` is empty;
- history starts from a clean public-release commit;
- no remote is configured until the final GitHub repository is selected;
- local-path searches using any known private path prefix return no matches.

## First Push

After the GitHub repository is created and the final visibility is confirmed:

```powershell
git remote add origin https://github.com/beloslav1976-svg/sepsis_blood_endotypes_biomarkers.git
git push -u origin main
```

The repository has already been published at:

`https://github.com/beloslav1976-svg/sepsis_blood_endotypes_biomarkers`

If a release tag is ready:

```powershell
git tag -a v0.1.0 -m "Sepsis blood endotypes and biomarkers publication package"
git push origin v0.1.0
```

## Metadata Finalization

Before tagging the final public release:

- copy `docs/CITATION.cff.template` to `CITATION.cff` and replace all
  placeholders;
- copy `docs/ZENODO_METADATA.template.json` to `.zenodo.json` and replace all
  placeholders;
- add the selected `LICENSE` file;
- update the manuscript data availability statement with the final GitHub URL
  and Zenodo DOI;
- regenerate `MANIFEST.csv` after any final file changes.
