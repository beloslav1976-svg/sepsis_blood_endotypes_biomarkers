# Deep Staged Workbook Privacy Scan

Scope: all populated string cells in staged `.xlsx` files.
This automated screen searches for obvious direct-contact patterns, calendar-date strings, and identifier-related labels.
It does not replace final manual privacy review.

Workbooks scanned: 13
String cells scanned: 228246
Findings: 7

## Workbook Summary

| Workbook | Sheet | Rows | Columns | String cells scanned |
|---|---:|---:|---:|---:|
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `00_Freeze_status` | 10 | 30 | 20 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `01_Audit_summary` | 12 | 30 | 48 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `02_Main_Table1` | 2 | 30 | 14 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `03_Main_Figures1_5` | 6 | 30 | 54 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `04_Source_provenance` | 7 | 30 | 29 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `05_Caption_audit` | 6 | 30 | 25 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `06_MD5_manifest` | 22 | 30 | 157 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `07_Numerical_guardrails` | 20 | 30 | 80 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `08_Story_map` | 7 | 30 | 21 |
| `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `09_Reporting_guardrails` | 19 | 30 | 38 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `00_Freeze_status` | 11 | 30 | 22 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `01_Audit_summary` | 16 | 30 | 64 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `02_Tables_S1_S10` | 11 | 30 | 92 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `03_Figures_S1_S8` | 9 | 30 | 86 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `04_MD5_freeze_manifest` | 43 | 30 | 304 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `05_Table_content_notes` | 11 | 30 | 36 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `06_Caption_notes` | 9 | 30 | 30 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `07_Workbook_sheet_audit` | 11 | 30 | 48 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `08_Numbering_audit` | 3 | 30 | 10 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `09_Source_provenance` | 19 | 30 | 77 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `10_Cross_reference_map` | 10 | 30 | 30 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `11_Reporting_guardrails` | 17 | 30 | 34 |
| `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `12_Terminology` | 11 | 30 | 22 |
| `publication_package\main_figures_tables\151_Table1_discovery_blood_cohort.xlsx` | `Table1` | 27 | 4 | 70 |
| `publication_package\main_figures_tables\151_Table1_discovery_blood_cohort.xlsx` | `Clinical_summary_source` | 11 | 6 | 26 |
| `publication_package\main_figures_tables\151_Table1_discovery_blood_cohort.xlsx` | `Missingness` | 11 | 5 | 15 |
| `publication_package\main_figures_tables\151_Table1_discovery_blood_cohort.xlsx` | `Validation_audit` | 25 | 2 | 26 |
| `publication_package\main_figures_tables\151_Table1_discovery_blood_cohort.xlsx` | `Source_cross_audit` | 8 | 5 | 19 |
| `publication_package\main_figures_tables\151_Table1_discovery_blood_cohort.xlsx` | `Raw_source_audit` | 8 | 2 | 16 |
| `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `S1_ReadMe` | 11 | 2 | 22 |
| `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `Participant_Metadata` | 46 | 18 | 338 |
| `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `Blood_Sample_Metadata` | 46 | 6 | 276 |
| `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `Data_Dictionary` | 23 | 4 | 92 |
| `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `Validation_Audit` | 25 | 2 | 26 |
| `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `Source_Cross_Audit` | 8 | 5 | 19 |
| `publication_package\supplementary_files\152_TableS2_complete_blood_differential_expression.xlsx` | `S2_ReadMe` | 11 | 2 | 22 |
| `publication_package\supplementary_files\152_TableS2_complete_blood_differential_expression.xlsx` | `Complete_DE` | 12401 | 14 | 39010 |
| `publication_package\supplementary_files\152_TableS3_robust_core_and_functional_enrichment.xlsx` | `S3_ReadMe` | 12 | 2 | 24 |
| `publication_package\supplementary_files\152_TableS3_robust_core_and_functional_enrichment.xlsx` | `Robust_core_full` | 1797 | 14 | 7198 |
| `publication_package\supplementary_files\152_TableS3_robust_core_and_functional_enrichment.xlsx` | `Robust_core_no_sex` | 1797 | 14 | 7198 |
| `publication_package\supplementary_files\152_TableS3_robust_core_and_functional_enrichment.xlsx` | `GO_BP_ORA` | 541 | 14 | 3794 |
| `publication_package\supplementary_files\152_TableS3_robust_core_and_functional_enrichment.xlsx` | `KEGG_ORA` | 64 | 16 | 583 |
| `publication_package\supplementary_files\152_TableS3_robust_core_and_functional_enrichment.xlsx` | `WikiPathways_ORA` | 54 | 14 | 385 |
| `publication_package\supplementary_files\152_TableS3_robust_core_and_functional_enrichment.xlsx` | `Hallmark_GSEA` | 51 | 16 | 416 |
| `publication_package\supplementary_files\152b_TableS2_FINAL_SUBMISSION_biological_targets.xlsx` | `S2_ReadMe` | 16 | 2 | 32 |
| `publication_package\supplementary_files\152b_TableS2_FINAL_SUBMISSION_biological_targets.xlsx` | `Complete_DE` | 12394 | 14 | 38989 |
| `publication_package\supplementary_files\155_TableS4_SRS_CTS_assignments_and_robustness.xlsx` | `S4_ReadMe` | 14 | 2 | 28 |
| `publication_package\supplementary_files\155_TableS4_SRS_CTS_assignments_and_robustness.xlsx` | `SRS_all_blood` | 46 | 7 | 232 |
| `publication_package\supplementary_files\155_TableS4_SRS_CTS_assignments_and_robustness.xlsx` | `Sepsis_SRS_CTS` | 36 | 7 | 182 |
| `publication_package\supplementary_files\155_TableS4_SRS_CTS_assignments_and_robustness.xlsx` | `SRS_summary` | 5 | 4 | 12 |
| `publication_package\supplementary_files\155_TableS4_SRS_CTS_assignments_and_robustness.xlsx` | `CTS_summary` | 4 | 3 | 6 |
| `publication_package\supplementary_files\155_TableS4_SRS_CTS_assignments_and_robustness.xlsx` | `CTSxSRS` | 4 | 4 | 7 |
| `publication_package\supplementary_files\155_TableS4_SRS_CTS_assignments_and_robustness.xlsx` | `Robustness_summary` | 10 | 4 | 40 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `S5_ReadMe` | 14 | 2 | 28 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `Candidate_pool` | 14 | 7 | 98 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `Panel_size_summary` | 9 | 10 | 18 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `Manuscript_panels` | 3 | 9 | 15 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `Manuscript_panel_genes` | 11 | 5 | 55 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `Recommended_blood` | 3 | 42 | 72 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `All_eligible_panels` | 5433 | 42 | 81522 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `DCAF17_forced_panels` | 2708 | 42 | 40647 |
| `publication_package\supplementary_files\156_TableS5_candidate_gene_pool_and_exhaustive_panel_screening.xlsx` | `Selection_provenance` | 14 | 2 | 28 |
| `publication_package\supplementary_files\157_TableS6_repeated_internal_cross_validation.xlsx` | `S6_ReadMe` | 12 | 2 | 24 |
| `publication_package\supplementary_files\157_TableS6_repeated_internal_cross_validation.xlsx` | `CV_summary` | 6 | 9 | 24 |
| `publication_package\supplementary_files\157_TableS6_repeated_internal_cross_validation.xlsx` | `CV_source_complete` | 6 | 7 | 12 |
| `publication_package\supplementary_files\157_TableS6_repeated_internal_cross_validation.xlsx` | `Score_definitions` | 15 | 4 | 60 |
| `publication_package\supplementary_files\157_TableS6_repeated_internal_cross_validation.xlsx` | `CV_method` | 13 | 2 | 26 |
| `publication_package\supplementary_files\157_TableS6_repeated_internal_cross_validation.xlsx` | `Interpretation` | 6 | 2 | 12 |
| `publication_package\supplementary_files\159_TableS7_complete_exploratory_clinical_associations.xlsx` | `S7_ReadMe` | 14 | 2 | 28 |
| `publication_package\supplementary_files\159_TableS7_complete_exploratory_clinical_associations.xlsx` | `Complete_60_tests` | 61 | 20 | 657 |
| `publication_package\supplementary_files\159_TableS7_complete_exploratory_clinical_associations.xlsx` | `Global_BH_significant` | 3 | 20 | 40 |
| `publication_package\supplementary_files\159_TableS7_complete_exploratory_clinical_associations.xlsx` | `Continuous_associations` | 41 | 20 | 443 |
| `publication_package\supplementary_files\159_TableS7_complete_exploratory_clinical_associations.xlsx` | `Categorical_associations` | 41 | 20 | 455 |
| `publication_package\supplementary_files\159_TableS7_complete_exploratory_clinical_associations.xlsx` | `Test_family_summary` | 9 | 3 | 19 |
| `publication_package\supplementary_files\159_TableS7_complete_exploratory_clinical_associations.xlsx` | `Multiplicity_summary` | 10 | 2 | 11 |
| `publication_package\supplementary_files\161_TableS8_published_signature_benchmarking.xlsx` | `S8_ReadMe` | 15 | 2 | 30 |
| `publication_package\supplementary_files\161_TableS8_published_signature_benchmarking.xlsx` | `Benchmark_summary` | 8 | 15 | 50 |
| `publication_package\supplementary_files\161_TableS8_published_signature_benchmarking.xlsx` | `Significance_summary` | 3 | 3 | 5 |
| `publication_package\supplementary_files\161_TableS8_published_signature_benchmarking.xlsx` | `Interpretation` | 7 | 2 | 14 |
| `publication_package\supplementary_files\161_TableS8_published_signature_benchmarking.xlsx` | `Reference_registry` | 8 | 4 | 32 |
| `publication_package\supplementary_files\161_TableS8_published_signature_benchmarking.xlsx` | `Gene_completeness` | 8 | 6 | 27 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `S10_ReadMe` | 22 | 2 | 44 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `External_summary` | 8 | 16 | 51 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Primary_SOFA` | 2 | 5 | 7 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Secondary_associations` | 5 | 15 | 31 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Component_gene_SOFA` | 6 | 8 | 23 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Adjusted_model` | 9 | 5 | 13 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Adjusted_summary` | 5 | 2 | 6 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Location_SOFA` | 6 | 6 | 11 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Location_pooled` | 2 | 6 | 6 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Location_SOFA_CI` | 6 | 8 | 13 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Location_pooled_fixed` | 2 | 6 | 6 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Scaling_sensitivity` | 3 | 4 | 6 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `SOFA_summary` | 2 | 7 | 7 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `SOFA_groups` | 3 | 5 | 7 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Mortality_context` | 3 | 6 | 8 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Site_summary` | 3 | 2 | 4 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Location_summary` | 6 | 2 | 7 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Frozen_gene_coverage` | 6 | 3 | 13 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Sample_scores` | 393 | 23 | 3065 |
| `publication_package\supplementary_files\165_TableS10_GSE185263_external_evaluation.xlsx` | `Interpretation` | 11 | 2 | 22 |

## Findings

| Kind | Workbook | Sheet | Row | Column | Value |
|---|---|---:|---:|---:|---|
| `calendar_date` | `freeze_records\main_package_168\168_complete_main_manuscript_package_freeze\tables\168_COMPLETE_Main_Manuscript_Package_Manifest.xlsx` | `00_Freeze_status` | 3 | 2 | 2026-08-21 23:38:27.325787 |
| `calendar_date` | `freeze_records\supplementary_package_167\167_complete_supplementary_package_freeze\tables\167_COMPLETE_Supplementary_Package_Manifest.xlsx` | `00_Freeze_status` | 3 | 2 | 2026-08-21 23:27:54.298627 |
| `identifier_label` | `publication_package\main_figures_tables\151_Table1_discovery_blood_cohort.xlsx` | `Validation_audit` | 8 | 1 | Unique participant IDs |
| `identifier_label` | `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `S1_ReadMe` | 2 | 2 | Supplementary Table S1. De-identified participant- and blood RNA-seq sample-level metadata for the discovery cohort |
| `identifier_label` | `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `S1_ReadMe` | 10 | 2 | Names, raw patient IDs, dates, birth years, original study codes, hospital identifiers, pathogen names, and free-text infection sites are not included. |
| `identifier_label` | `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `Data_Dictionary` | 2 | 2 | De-identified participant identifier. |
| `identifier_label` | `publication_package\supplementary_files\151_TableS1_deidentified_blood_metadata.xlsx` | `Validation_Audit` | 8 | 1 | Unique participant IDs |

## Interpretation Guidance

- Automated conclusion for this run: no `email`, `phone_like`, or `personal_code_label` findings were detected.
- `identifier_label` and `personal_code_label` findings may be benign when they occur in readme/data-dictionary text documenting de-identification.
- `calendar_date`, `email`, and `phone_like` findings should be manually reviewed before public release.
- De-identified participant IDs can still be sensitive depending on the data-use agreement and should be approved by the study team before publication.
