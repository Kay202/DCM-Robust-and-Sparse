# Task-version sensitivity analysis for corrected NARPS DCM results

script_directory <- function() {
  sourced_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile) || !length(frame$ofile)) "" else as.character(frame$ofile[1])
  }, character(1))
  sourced_files <- sourced_files[nzchar(sourced_files)]
  if (length(sourced_files)) return(dirname(normalizePath(tail(sourced_files, 1L))))

  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", all_args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  normalizePath(getwd())
}
script_dir <- script_directory()
repository_dir <- normalizePath(file.path(script_dir, ".."))
deposited_results_dir <- file.path(script_dir, "results")
results_dir <- Sys.getenv(
  "RSDCM_REAL_RESULTS_DIR",
  unset = file.path(script_dir, "results_regenerated")
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

source(
  file.path(repository_dir, "rsDCM.R"),
  chdir = TRUE
)
source(file.path(script_dir, "spm_dcm_peb_port.R"), chdir = TRUE)

res2 <- readRDS(
  file.path(script_dir, "res2_corrected.rds")
)
parts <- read.csv(
  file.path(
    script_dir, "included_48_participants.csv"
  ),
  stringsAsFactors = FALSE
)

eta <- res2$eta
Cp <- res2$Cp_list
pnames <- res2$parameter_names
N <- nrow(eta)
p <- ncol(eta)

sid <- sub(".*(sub-[0-9]+).*", "\\1", res2$subject_names)
idx <- match(sid, parts$participant_id)

if (anyNA(idx))
  stop("Some analyzed participants could not be matched.")

meta <- parts[idx, ]

task_range <- as.integer(meta$task_version == "equalRange")
gender_m <- as.integer(meta$gender == "M")

# Center binary covariates so that the intercept remains the sample average.
task_range_c <- task_range - mean(task_range)
gender_m_c <- gender_m - mean(gender_m)
age_z <- as.numeric(scale(meta$age))

X_task <- cbind(
  intercept = 1,
  task_equalRange = task_range_c,
  genderM = gender_m_c,
  age = age_z
)

V_list <- lapply(seq_len(p), function(k) {
  Vk <- matrix(0, p, p)
  Vk[k, k] <- 1
  Vk
})

set.seed(20260718)

fit_task <- fit_group_dcm_t_MOM(
  eta_theta_y = eta,
  C_theta_y_list = Cp,
  X_G = X_task,
  V_list = V_list,
  nu = 3L,
  pi = 0.5,
  tau0 = 0.05,
  tau1 = 1,
  estimate_sigma2 = FALSE,
  sigma2_fixed = 1,
  whitening_method = "precision_chol",
  n_starts = 3L,
  start_seed = 20260718L,
  max_iter = 2000L,
  tol = 1e-6,
  pip_mcmc_iter = 30000L,
  pip_mcmc_burn = 10000L,
  pip_mcmc_thin = 5L,
  pip_mcmc_chains = 4L,
  pip_mcmc_temperatures = 16L,
  pip_mcmc_min_inverse_temp = 0.05,
  pip_mcmc_swap_every = 2L,
  verbose = TRUE
)

rownames(fit_task$beta_mat) <- pnames
rownames(fit_task$inclusion) <- pnames

# Matching SPM-PEB sensitivity model
pebs <- lapply(res2$PEB_list, .normalize_mat_peb)

M <- list(
  X = X_task,
  Xnames = as.list(colnames(X_task)),
  Q = "single",
  alpha = 1,
  beta = 16,
  maxit = 64
)

PEB_task <- spm_dcm_peb_run(
  spm_dcm_peb_prepare(P = pebs, M = M, field = seq_len(p)),
  verbose = FALSE
)$PEB

peb_block <- function(k) {
  ii <- ((k - 1L) * p + 1L):(k * p)
  list(
    estimate = as.numeric(PEB_task$Ep)[ii],
    sd = sqrt(pmax(diag(PEB_task$Cp[ii, ii, drop = FALSE]), 0))
  )
}

peb_intercept <- peb_block(1L)
peb_task_effect <- peb_block(2L)

primary_results_dir <- if (file.exists(file.path(results_dir, "real_data_fits.rds"))) {
  results_dir
} else {
  deposited_results_dir
}
primary <- readRDS(file.path(primary_results_dir, "real_data_fits.rds"))$fit_int

task_results <- data.frame(
  parameter = pnames,
  primary_intercept_estimate = primary$beta_mat[, "intercept"],
  primary_intercept_PIP = primary$inclusion[, "intercept"],
  adjusted_intercept_estimate = fit_task$beta_mat[, "intercept"],
  adjusted_intercept_PIP = fit_task$inclusion[, "intercept"],
  task_equalRange_estimate = fit_task$beta_mat[, "task_equalRange"],
  task_equalRange_PIP = fit_task$inclusion[, "task_equalRange"],
  SPM_adjusted_intercept = peb_intercept$estimate,
  SPM_adjusted_intercept_SD = peb_intercept$sd,
  SPM_task_equalRange = peb_task_effect$estimate,
  SPM_task_equalRange_SD = peb_task_effect$sd
)

write.csv(
  task_results,
  file.path(results_dir, "task_version_sensitivity.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    fit_task = fit_task,
    PEB_task = PEB_task,
    design = X_task,
    metadata = meta
  ),
  file.path(results_dir, "task_version_sensitivity_fits.rds"),
  compress = "xz"
)

cat("\nTask-version composition:\n")
print(table(meta$task_version))

cat("\nFit diagnostics:\n")
print(fit_task[c(
  "converged", "iterations", "objective_monotone",
  "q_monotone", "pip_max_chain_range"
)])
cat("Maximum PIP MCSE:",
    max(fit_task$inclusion_mcse, na.rm = TRUE), "\n")
print(fit_task$multi_start_summary)

cat("\nAdjusted intercept effects with PIP > 0.5:\n")
print(
  task_results[
    task_results$adjusted_intercept_PIP > 0.5,
    c(
      "parameter", "adjusted_intercept_estimate",
      "adjusted_intercept_PIP"
    )
  ]
)

cat("\nTask-version effects with PIP > 0.5:\n")
print(
  task_results[
    task_results$task_equalRange_PIP > 0.5,
    c(
      "parameter", "task_equalRange_estimate",
      "task_equalRange_PIP"
    )
  ]
)
