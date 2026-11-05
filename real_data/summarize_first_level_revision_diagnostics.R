# Summarize the five CSV files produced by
# extract_first_level_revision_diagnostics.m.
# Run from the new 13_revision_diagnostics_YYYYMMDD_HHMMSS folder:
# source("/full/path/to/summarize_first_level_revision_diagnostics.R")

required_files <- c(
  "subject_run_condition_counts.csv",
  "subject_run_design_ranks.csv",
  "run_level_dcm_fit_diagnostics.csv",
  "run_level_variance_changes.csv",
  "within_subject_peb_variance_changes.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing files: ", paste(missing_files, collapse = ", "))
}

counts <- read.csv("subject_run_condition_counts.csv", check.names = FALSE)
ranks <- read.csv("subject_run_design_ranks.csv", check.names = FALSE)
fits <- read.csv("run_level_dcm_fit_diagnostics.csv", check.names = FALSE)
run_var <- read.csv("run_level_variance_changes.csv", check.names = FALSE)
peb_var <- read.csv("within_subject_peb_variance_changes.csv", check.names = FALSE)

for (object_name in c("counts", "ranks", "fits", "run_var", "peb_var")) {
  object <- get(object_name)
  object$subject <- as.character(object$subject)
  object$task_version <- as.character(object$task_version)
  object$task_version[is.na(object$task_version) | object$task_version == ""] <- "not_merged"
  assign(object_name, object)
}

# Fill absent conditions with zero within each task version. This is important
# when a condition is missing entirely from one run's SPM.Sess(r).U structure.
count_groups <- split(counts, counts$task_version)
counts_complete <- do.call(rbind, lapply(count_groups, function(d) {
  grid <- expand.grid(
    subject = unique(d$subject),
    run = sort(unique(d$run)),
    condition = sort(unique(d$condition)),
    stringsAsFactors = FALSE
  )
  grid$task_version <- unique(d$task_version)[1]
  out <- merge(
    grid,
    d[c("subject", "task_version", "run", "condition", "n_events")],
    by = c("subject", "task_version", "run", "condition"),
    all.x = TRUE,
    sort = FALSE
  )
  out$n_events[is.na(out$n_events)] <- 0
  out
}))
rownames(counts_complete) <- NULL
write.csv(counts_complete, "subject_run_condition_counts_complete.csv", row.names = FALSE)

summarize_numeric <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(n = 0, min = NA, q1 = NA, median = NA, q3 = NA, max = NA))
  }
  c(
    n = length(x),
    min = min(x),
    q1 = unname(quantile(x, 0.25)),
    median = median(x),
    q3 = unname(quantile(x, 0.75)),
    max = max(x)
  )
}

split_apply <- function(data, grouping, function_for_group) {
  keys <- interaction(data[grouping], drop = TRUE, lex.order = TRUE)
  output <- lapply(split(data, keys), function_for_group)
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  result
}

count_summary <- split_apply(
  counts_complete,
  c("task_version", "run", "condition"),
  function(d) {
    s <- summarize_numeric(d$n_events)
    data.frame(
      task_version = d$task_version[1],
      run = d$run[1],
      condition = d$condition[1],
      n_subjects = length(unique(d$subject)),
      subjects_with_zero_events = sum(d$n_events == 0),
      min_events = s["min"],
      median_events = s["median"],
      max_events = s["max"],
      total_events = sum(d$n_events)
    )
  }
)
write.csv(count_summary, "revision_condition_count_summary.csv", row.names = FALSE)

rank_summary <- split_apply(
  ranks,
  "task_version",
  function(d) {
    data.frame(
      task_version = d$task_version[1],
      subject_runs = nrow(d),
      dcm_input_rank_min = min(d$dcm_input_rank),
      dcm_input_rank_median = median(d$dcm_input_rank),
      dcm_input_rank_max = max(d$dcm_input_rank),
      dcm_input_full_rank_fraction = mean(d$dcm_input_rank == d$n_nonzero_dcm_input_columns),
      glm_task_rank_min = min(d$glm_task_rank),
      glm_task_rank_median = median(d$glm_task_rank),
      glm_task_rank_max = max(d$glm_task_rank),
      glm_task_full_rank_fraction = mean(d$glm_task_rank == d$n_glm_task_columns)
    )
  }
)
write.csv(rank_summary, "revision_rank_summary.csv", row.names = FALSE)

parameter_group <- function(parameter) {
  ifelse(
    grepl("^A\\(", parameter),
    "intrinsic_A",
    ifelse(grepl(",1\\)$", parameter), "gain_modulatory_B", "loss_modulatory_B")
  )
}

run_var$parameter_group <- parameter_group(run_var$parameter)
peb_var$parameter_group <- parameter_group(peb_var$parameter)
variance_all <- rbind(run_var, peb_var)

variance_summary <- split_apply(
  variance_all,
  c("source", "task_version", "parameter"),
  function(d) {
    x <- d$posterior_to_prior_ratio
    x <- x[is.finite(x)]
    data.frame(
      source = d$source[1],
      task_version = d$task_version[1],
      parameter_group = d$parameter_group[1],
      parameter = d$parameter[1],
      n = length(x),
      median_posterior_to_prior_ratio = if (length(x)) median(x) else NA,
      q1_posterior_to_prior_ratio = if (length(x)) unname(quantile(x, 0.25)) else NA,
      q3_posterior_to_prior_ratio = if (length(x)) unname(quantile(x, 0.75)) else NA,
      median_variance_reduction = if (length(x)) median(1 - x) else NA,
      fraction_with_variance_reduction = if (length(x)) mean(x < 1) else NA
    )
  }
)
write.csv(variance_summary, "revision_parameter_variance_summary.csv", row.names = FALSE)

fit_summary <- split_apply(
  fits,
  "task_version",
  function(d) {
    f <- summarize_numeric(d$free_energy)
    p <- summarize_numeric(d$variance_explained_percent)
    r <- summarize_numeric(d$residual_rms)
    data.frame(
      task_version = d$task_version[1],
      subject_runs = nrow(d),
      free_energy_median = f["median"],
      free_energy_q1 = f["q1"],
      free_energy_q3 = f["q3"],
      variance_explained_median = p["median"],
      variance_explained_q1 = p["q1"],
      variance_explained_q3 = p["q3"],
      variance_explained_min = p["min"],
      residual_rms_median = r["median"],
      covariance_finite_fraction = mean(d$posterior_covariance_finite),
      covariance_psd_fraction = mean(d$posterior_covariance_positive_semidefinite)
    )
  }
)
write.csv(fit_summary, "revision_first_level_fit_summary.csv", row.names = FALSE)

# A compact text report for drafting the manuscript. Complete subject-run
# results remain in the CSV files and should be deposited in the repository.
sink("revision_diagnostics_report.txt")
cat("TASK-VERSION AND CONDITION COUNTS\n")
print(count_summary, row.names = FALSE)
cat("\nINPUT AND DESIGN RANKS\n")
print(rank_summary, row.names = FALSE)
cat("\nPRIOR-TO-POSTERIOR VARIANCE RATIOS FOR MODULATORY PARAMETERS\n")
print(
  variance_summary[grepl("modulatory", variance_summary$parameter_group), ],
  row.names = FALSE
)
cat("\nFIRST-LEVEL DCM FIT DIAGNOSTICS\n")
print(fit_summary, row.names = FALSE)
sink()

cat(
  "Wrote complete counts and four manuscript-summary files in:\n",
  normalizePath("."), "\n",
  "Please inspect revision_diagnostics_report.txt before drafting numerical text.\n",
  sep = ""
)
