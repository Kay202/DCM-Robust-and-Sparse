# NARPS real-data workflow

## Deposited group-level input

`res2_corrected.rds` contains de-identified posterior means, posterior covariance matrices, parameter names, subject identifiers, and within-subject PEB objects for the 48 analyzed participants. `included_48_participants.csv` contains the task version, gender, and age used in the group design.

The sample consists of the first 48 available participants in ascending subject-ID order after excluding `sub-020` and `sub-021`, whose NIfTI files were not detected during the original local processing stage. `verify_sample_selection.R` checks this rule against `participants_full.tsv`.

## Group-level analysis

From the repository root:

```bash
Rscript real_data/real_data_analysis.R
Rscript real_data/task_version_sensitivity.R
Rscript real_data/make_real_data_tables.R
```

The first two commands write to `real_data/results_regenerated/`. They do not overwrite the deposited `real_data/results/` folder. Set `RSDCM_REAL_RESULTS_DIR` to choose another output directory.

The table script uses regenerated results when they exist. Otherwise, it uses the deposited results. It writes `real_data_tables_generated.tex` into the selected result directory.

## First-level reproduction

Raw NARPS imaging data are available from OpenNeuro accession `ds001734` and are not included here. The first-level workflow requires MATLAB and SPM25. Set these environment variables or pass the equivalent paths directly to the MATLAB functions:

- `NARPS_ROOT`: local NARPS dataset and derivative root.
- `NARPS_REVISION_ROOT`: directory containing corrected run-level DCM and within-subject PEB folders.
- `SPM_DIR`: SPM25 installation directory.
- `AAL_ATLAS`: local AAL V4 atlas image used to build the four ROI masks.

The main order is:

1. `first_level_glm_multiple.m`
2. `make_aal_rois_v4.m`
3. `extract_vois_all_multiple.m`
4. `audit_first_level_dcm_inputs.m`
5. `multiple_dcm_runs.m`
6. `run_within_subject_peb_multiple.m`
7. `export_corrected_pebs.R`
8. `extract_first_level_revision_diagnostics.m`
9. `summarize_first_level_revision_diagnostics.R`

The raw imaging data, run-level DCM `.mat` files, and within-subject PEB `.mat` files are not redistributed. Their derived, de-identified group-level export and the complete diagnostic CSV files are included.

## Parameter scale

Off-diagonal intrinsic and modulatory parameters are reported in Hz. SPM stores diagonal intrinsic parameters as dimensionless log scalings. The physical self-inhibition shown in the connectivity figure is

```text
a_ii in Hz = -0.5 * exp(internal diagonal parameter).
```

PIPs for diagonal parameters quantify evidence for departure from SPM's default self-inhibition of -0.5 Hz.

