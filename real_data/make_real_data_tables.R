# Generate the real-data LaTeX tables from deposited or regenerated results.

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

regenerated_dir <- file.path(script_dir, "results_regenerated")
deposited_dir <- file.path(script_dir, "results")
result_dir <- Sys.getenv(
  "RSDCM_REAL_RESULTS_DIR",
  unset = if (file.exists(file.path(regenerated_dir, "real_data_fits.rds"))) {
    regenerated_dir
  } else {
    deposited_dir
  }
)

required <- c(
  "real_data_fits.rds", "method_comparison.csv", "covariate_comparison.csv"
)
missing <- required[!file.exists(file.path(result_dir, required))]
if (length(missing)) stop("Missing result file(s): ", paste(missing, collapse = ", "))

fits <- readRDS(file.path(result_dir, "real_data_fits.rds"))
intercept <- read.csv(file.path(result_dir, "method_comparison.csv"), check.names = FALSE)
covariates <- read.csv(file.path(result_dir, "covariate_comparison.csv"), check.names = FALSE)

sd_matrix <- function(fit) {
  if (!is.null(fit$beta_sd_approx)) {
    out <- as.matrix(fit$beta_sd_approx)
    if (all(dim(out) == dim(fit$beta_mat))) {
      dimnames(out) <- dimnames(fit$beta_mat)
      return(out)
    }
    return(matrix(as.numeric(out), nrow = nrow(fit$beta_mat),
                  ncol = ncol(fit$beta_mat), dimnames = dimnames(fit$beta_mat)))
  }
  if (is.null(fit$beta_cov_approx)) stop("Conditional Laplace covariance is unavailable.")
  matrix(
    sqrt(pmax(diag(fit$beta_cov_approx), 0)),
    nrow = nrow(fit$beta_mat), ncol = ncol(fit$beta_mat),
    dimnames = dimnames(fit$beta_mat)
  )
}

sd_int <- sd_matrix(fits$fit_int)[, "intercept"]
sd_cov <- sd_matrix(fits$fit_cov)

fmt <- function(x, digits = 3L) formatC(x, digits = digits, format = "f")
signed <- function(x) sprintf("%+.3f", x)
pip_cell <- function(x) {
  value <- fmt(x)
  if (is.finite(x) && x > 0.5) paste0("\\textbf{", value, "}") else value
}
pm_cell <- function(estimate, standard_error) {
  paste0("$", signed(estimate), "\\!\\pm\\!", fmt(standard_error), "$")
}
row <- function(x) paste0(paste(x, collapse = " & "), " \\\\")

table_intercept <- c(
  "\\begin{table}[tbp]",
  "\\centering",
  "\\small",
  "\\caption{Group-average (intercept) estimates of all $22$ DCM parameters under standard \\emph{SPM-PEB} and the \\emph{Proposed} (robust $+$ sparse) model. $A(i,j)$ denotes the intrinsic connection from region $j$ to region $i$, and $B(i,j,k)$ denotes input $k$ modulating that connection (region indices $1$--$4$ are vmPFC, vStr, amygdala, and aIns). \\emph{SPM-PEB} reports a posterior mean $\\pm$ SD, with \\textsuperscript{$\\ast$} marking $|\\mathrm{E}p|>3\\,\\mathrm{SD}$. The \\emph{Proposed} model reports a coefficient estimate $\\pm$ conditional Laplace SD and a PIP, with bold PIPs exceeding $0.5$.}",
  "\\label{tab:realparams}",
  "\\begin{tabular}{l r@{\\,$\\pm$\\,}l r@{\\,$\\pm$\\,}l c}",
  "\\toprule",
  row(c("Parameter", "\\multicolumn{2}{c}{\\emph{SPM-PEB} ($\\mathrm{E}p\\pm$SD)}", "\\multicolumn{2}{c}{\\emph{Proposed} ($\\widehat{\\beta}\\pm$SD)}", "\\emph{Proposed} PIP")),
  "\\midrule"
)
for (i in seq_len(nrow(intercept))) {
  spm_estimate <- signed(intercept$PEB_Ep[i])
  if (isTRUE(intercept$PEB_strong[i])) {
    spm_estimate <- paste0(spm_estimate, "\\textsuperscript{$\\ast$}")
  }
  table_intercept <- c(table_intercept, row(c(
    intercept$parameter[i], spm_estimate, fmt(intercept$PEB_SD[i]),
    signed(intercept$Robust_est[i]), fmt(sd_int[i]), pip_cell(intercept$Robust_PIP[i])
  )))
}
table_intercept <- c(table_intercept, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")

table_covariates <- c(
  "\\begin{table}[tbp]",
  "\\centering",
  "\\small",
  "\\caption{Gender and age covariate effects on all $22$ NARPS DCM parameters under standard \\emph{SPM-PEB} and the \\emph{Proposed} (robust $+$ sparse) model. Gender is coded as M$-$F and age is measured per standard deviation. \\emph{SPM-PEB} reports a posterior mean $\\pm$ SD. The \\emph{Proposed} model reports a coefficient estimate $\\pm$ conditional Laplace SD and a PIP, with bold PIPs exceeding $0.5$. Effects on diagonal $A$ parameters are on SPM's internal log-scaling scale. Parameter naming follows Table~\\ref{tab:realparams}.}",
  "\\label{tab:realcov}",
  "\\setlength{\\tabcolsep}{2.5pt}",
  "\\begin{tabular}{l ccc ccc}",
  "\\toprule",
  row(c("", "\\multicolumn{3}{c}{Gender (M$-$F)}", "\\multicolumn{3}{c}{Age (per SD)}")),
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  row(c("Parameter", "PEB\\,$\\pm$\\,SD", "\\emph{Proposed} $\\widehat{\\beta}\\pm$SD", "PIP", "PEB\\,$\\pm$\\,SD", "\\emph{Proposed} $\\widehat{\\beta}\\pm$SD", "PIP")),
  "\\midrule"
)
for (i in seq_len(nrow(covariates))) {
  table_covariates <- c(table_covariates, row(c(
    covariates$parameter[i],
    pm_cell(covariates$gender_PEB[i], covariates$gender_PEB_SD[i]),
    pm_cell(covariates$gender_Robust[i], sd_cov[i, "genderM"]),
    pip_cell(covariates$gender_PIP[i]),
    pm_cell(covariates$age_PEB[i], covariates$age_PEB_SD[i]),
    pm_cell(covariates$age_Robust[i], sd_cov[i, "age"]),
    pip_cell(covariates$age_PIP[i])
  )))
}
table_covariates <- c(
  table_covariates, "\\bottomrule", "\\end{tabular}", "\\end{table}", ""
)

output_file <- file.path(result_dir, "real_data_tables_generated.tex")
writeLines(c(table_intercept, table_covariates), output_file)
message("Wrote: ", output_file)

