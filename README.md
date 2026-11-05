# Robust and Sparse Group Dynamic Causal Modeling

This repository contains the reference R implementation and reproducibility materials for **“Robust and Sparse Group Dynamic Causal Modeling via Student-t Parametric Empirical Bayes and Nonlocal Priors”** by Godfred Arhin and Nilotpal Sanyal.

The proposed group-level method combines a Student-t likelihood with a product-moment spike-and-slab prior. It preserves first-level posterior covariance information, uses generalized expectation-maximization with conditional maximization for mode fitting, and computes posterior inclusion probabilities through parallel-tempering posterior integration.

`rsDCM.R` is a standalone reference implementation. This repository is not an installable R package.

## Repository contents

- `rsDCM.R`: proposed robust and sparse group-DCM estimator.
- `simulation/`: simulation driver, SPM-PEB reference port, manuscript-output script, empirical covariance input, and deposited results.
- `real_data/`: NARPS first-level and group-level analysis scripts, corrected group-level input, participant metadata, diagnostics, and deposited results.
- `capture_session_info.R`: writes R and package-version information for configured analysis environments.

## Software requirements

The R workflows use R 4.5.x and the packages `Rcpp`, `Matrix`, `MASS`, `expm`, `R.matlab`, `snowfall`, `snow`, and `Rmpi`. The optional Windows parallel path in the SPM port uses `future` and `future.apply`. Compiling the parallel-tempering routine requires a working C++ toolchain.

The first-level neuroimaging workflow requires MATLAB and SPM25. SPM must be obtained separately. MPI is required only for the parallel full simulation workflow.

## Recreate simulation tables and figures

From the repository root, run:

```bash
Rscript simulation/simulation_manuscript_outputs.R
```

This reads the deposited CSV files and writes the regenerated LaTeX tables and figures to `simulation/results/`. It does not refit any model.

The complete simulation can be run stage by stage with `simulation/simulation_study.R`. See [simulation/README.md](simulation/README.md). Quick runs are written to `simulation/results_quick/` and cannot overwrite the deposited full results.

## Recreate the group-level NARPS analysis

The file `real_data/res2_corrected.rds` contains de-identified subject-level posterior summaries for the 48 analyzed participants. From the repository root, run:

```bash
Rscript real_data/real_data_analysis.R
Rscript real_data/task_version_sensitivity.R
Rscript real_data/make_real_data_tables.R
```

Regenerated outputs are written to `real_data/results_regenerated/`. Deposited results under `real_data/results/` are never overwritten by default. See [real_data/README.md](real_data/README.md) for the first-level workflow and data requirements.

## Data availability

The empirical analysis uses the openly shared NARPS mixed-gambles dataset, OpenNeuro accession `ds001734`. Raw imaging data are not redistributed here. The repository contains only de-identified participant metadata, derived posterior summaries, covariance matrices, diagnostics, and aggregate results needed to reproduce the reported group analysis.

The analyzed sample consists of the first 48 available participants in ascending subject-ID order after excluding `sub-020` and `sub-021`, whose NIfTI files were not detected during the original local processing stage. Run the following command to verify the selection from the deposited metadata:

```bash
Rscript real_data/verify_sample_selection.R
```

## Reproducibility notes

- The final simulation diagnostic files report convergence for every retained fit.
- Diagonal DCM `A` parameters are stored internally as log scalings. Figures reporting physical self-inhibition use the transformation \(a_{ii}^{\mathrm{Hz}}=-0.5\exp(\widetilde a_{ii})\).
- The Student-t robustness weights apply to coordinates defined by the selected pre-whitening transformation.

## Citation

Please cite the accompanying manuscript and this repository. Machine-readable citation metadata are provided in `CITATION.cff`.

## License and attribution

The authors' reference implementation is distributed under GPL-2.0-or-later. See `LICENSE`. The SPM-related R port is derived from the SPM computational framework. See `NOTICE` and the source-file headers. SPM itself is not included.

