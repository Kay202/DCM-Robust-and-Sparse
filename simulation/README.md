# Simulation workflow

## Deposited results

The `results/` folder contains the complete raw, coefficient-level, summary, diagnostic, tuning, validation, and manuscript-output files from the final simulation study. All retained fits in the final diagnostic files are marked as converged.

## Requirements

Install the packages listed in the root README. Full parallel runs require an MPI installation compatible with `Rmpi`, `snow`, and `snowfall`. The environment variable `RSDCM_CORES` controls the number of workers.

## Recreate manuscript outputs without refitting

From the repository root:

```bash
Rscript simulation/simulation_manuscript_outputs.R
```

This recreates `simulation_tables_generated.tex`, `fig_simulation_main.png`, and `fig_simulation_sparse.png` in `simulation/results/`. The generated TeX also contains the empirical-covariance sensitivity table.

## Quick smoke test

```bash
RSDCM_CORES=4 Rscript simulation/simulation_study.R tune quick
RSDCM_CORES=4 Rscript simulation/simulation_study.R main quick
```

Quick outputs are written to `simulation/results_quick/`.

## Full analysis order

Run stages separately so that any failed fit can be examined before continuing:

```bash
RSDCM_CORES=39 Rscript simulation/simulation_study.R tune
RSDCM_CORES=39 Rscript simulation/simulation_study.R main
RSDCM_CORES=39 Rscript simulation/simulation_study.R moderate
RSDCM_CORES=39 Rscript simulation/simulation_study.R breadth
RSDCM_CORES=39 Rscript simulation/simulation_study.R sparse
RSDCM_CORES=39 Rscript simulation/simulation_study.R empirical
RSDCM_CORES=39 Rscript simulation/simulation_study.R whitening
RSDCM_CORES=39 Rscript simulation/simulation_study.R validation
RSDCM_CORES=39 Rscript simulation/simulation_study.R nu
RSDCM_CORES=39 Rscript simulation/simulation_study.R fairness
Rscript simulation/simulation_study.R report
Rscript simulation/simulation_manuscript_outputs.R
```

Each stage writes a separate `sim_run_settings_<stage>.csv`. Set `RSDCM_RESULTS_DIR` to direct a new run to another folder.

## Principal inputs

- `../rsDCM.R`: proposed estimator.
- `spm_dcm_peb_port.R`: standard SPM-PEB reference implementation in R.
- `exact_core_final_2_optimized.R` and `helpers_final_optimized.R`: dependencies of the SPM port.
- `narps_covariances.rds`: corrected empirical covariance matrices used in the NARPS covariance sensitivity experiment.

