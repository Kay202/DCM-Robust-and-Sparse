# =============================================================================
# Simulation study for robust and sparse group DCM
#
# This script locates the repository files relative to its own location.
# Full-run outputs are written to simulation/results. Quick-run outputs are
# written to simulation/results_quick so they cannot overwrite full results.
# Set RSDCM_RESULTS_DIR to use a different output directory.
#
# Examples:
#   Run "tune" before any other individual stage.
#   Rscript simulation_study.R tune
#   Rscript simulation_study.R main
#   Rscript simulation_study.R all
#   Rscript simulation_study.R main quick
#   RSDCM_CORES=39 Rscript simulation_study.R all

# Run in this order
# RSDCM_CORES=39 Rscript simulation_study.R tune
# RSDCM_CORES=39 Rscript simulation_study.R main
# RSDCM_CORES=39 Rscript simulation_study.R moderate
# RSDCM_CORES=39 Rscript simulation_study.R breadth
# RSDCM_CORES=39 Rscript simulation_study.R sparse
# RSDCM_CORES=39 Rscript simulation_study.R whitening
# RSDCM_CORES=39 Rscript simulation_study.R validation
# RSDCM_CORES=39 Rscript simulation_study.R nu
# RSDCM_CORES=39 Rscript simulation_study.R fairness
# RSDCM_CORES=35 Rscript simulation_study.R empirical
# Rscript simulation_study.R report
# =============================================================================

library(snowfall)

args <- commandArgs(trailingOnly = TRUE)
STAGE <- if (length(args)) args[1] else "all"
QUICK <- "quick" %in% args

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
SCRIPT_DIR <- script_directory()
REPOSITORY_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
default_results <- file.path(SCRIPT_DIR, if (QUICK) "results_quick" else "results")
RESULTS_DIR <- Sys.getenv("RSDCM_RESULTS_DIR", unset = default_results)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

METHOD_FILE <- file.path(REPOSITORY_DIR, "rsDCM.R")
SPM_FILE <- file.path(SCRIPT_DIR, "spm_dcm_peb_port.R")
NARPS_COVARIANCE_FILE <- file.path(SCRIPT_DIR, "narps_covariances.rds")

source(METHOD_FILE)
source(SPM_FILE, chdir = TRUE)
setwd(RESULTS_DIR)

N_REP <- if (QUICK) 3L else 50L #100L
N_TUNE <- if (QUICK) 2L else 5L #20L
N_CORE <- as.integer(Sys.getenv("RSDCM_CORES", unset = "39"))

N_STARTS <- 3L #if (QUICK) 2L else 3L
PIP_ITER <- if (QUICK) 1200L else 5000L #10000L
PIP_BURN <- if (QUICK) 300L else 1500L  #3000L
PIP_THIN <- 5L
PIP_CHAINS <- 4L
PIP_TEMPERATURES <- 8L

SEED <- list(
  main_truth = 100000L,
  main_design = 150000L,
  main_data = 500000L,
  tune_truth = 200000L,
  tune_data = 600000L,
  sparse_truth = 300000L,
  sparse_data = 700000L,
  breadth_truth = 400000L,
  breadth_design = 450000L,
  breadth_data = 800000L,
  validation_truth = 900000L,
  validation_data = 950000L,
  mcmc = 1200000L
)

diag_basis <- function(p) {
  lapply(seq_len(p), function(j) {
    V <- matrix(0, p, p)
    V[j, j] <- 1
    V
  })
}

parallel_map <- function(tasks, FUN, ...) {
  sfClusterApplyLB(tasks, FUN, ...)
}

make_tasks <- function(n_condition, n_rep = N_REP) {
  unlist(lapply(seq_len(n_condition), function(i) {
    lapply(seq_len(n_rep), function(rep) list(condition = i, rep = rep))
  }), recursive = FALSE)
}

# ---- Data-generating process -------------------------------------------------

make_truth <- function(p, n_regressor, seed) {
  set.seed(seed)
  gain <- c(p - 5L, p - 4L)
  loss <- c(p - 3L, p - 2L, p - 1L, p)
  beta <- matrix(0, p, n_regressor)

  inclusion_uniform <- runif(1000)
  signs <- sample(c(-1, 1), 1000, replace = TRUE)
  magnitudes <- runif(1000, 0.3, 0.6)
  intrinsic <- which(inclusion_uniform[seq_len(p - 6L)] < 0.25)
  if (length(intrinsic) < min(3L, p - 6L)) {
    intrinsic <- sort(order(inclusion_uniform[seq_len(p - 6L)])[seq_len(min(3L, p - 6L))])
  }
  beta[intrinsic, 1] <- signs[intrinsic] * magnitudes[intrinsic]
  beta[gain, 1] <- c(0.30, 0.25)
  beta[loss, 1] <- -c(0.28, 0.10, 0.30, 0.08)

  if (n_regressor >= 2L) beta[c(gain[1], loss[1]), 2] <- c(0.35, 0.30)
  if (n_regressor >= 3L) beta[loss[3], 3] <- -0.25

  list(beta = beta, gain = gain, loss = loss)
}

make_sparse_truth <- function(p, active_fraction, seed) {
  set.seed(seed)
  beta <- matrix(0, p, 1)
  active <- sample(seq_len(p), max(1L, round(active_fraction * p)))
  beta[active, 1] <- sample(c(-1, 1), length(active), replace = TRUE) *
    runif(length(active), 0.4, 0.7)
  list(beta = beta, gain = 1:2, loss = 3:6)
}

make_covariance <- function(p, regime) {
  if (regime == "small_diagonal") return(diag(0.05, p))
  if (regime == "moderate_diagonal") return(diag(0.20, p))
  if (regime == "large_diagonal") return(diag(0.60, p))

  Z <- matrix(rnorm(p * p), p, p)
  S <- tcrossprod(Z) / p
  D <- diag(1 / sqrt(diag(S)))
  0.20 * D %*% S %*% D
}

simulate_data <- function(
    N,
    beta_true,
    X,
    gain,
    loss,
    between_subject_sd = 0.15,
    covariance_regime = "moderate_diagonal",
    heterogeneous_covariance = FALSE,
    covariance_pool = NULL,
    contamination_fraction = 0.08,
    contamination_magnitude = 6,
    contamination_geometry = "cell",
    contamination_subgroup = NULL,
    seed
) {
  set.seed(seed)
  p <- nrow(beta_true)

  if (is.null(covariance_pool)) {
    common_covariance <- make_covariance(p, covariance_regime)
    C_list <- lapply(seq_len(N), function(i) {
      if (heterogeneous_covariance) make_covariance(p, covariance_regime) else common_covariance
    })
  } else {
    selected <- sample(seq_along(covariance_pool), N, replace = TRUE)
    C_list <- lapply(selected, function(i) covariance_pool[[i]])
  }

  mean_matrix <- X %*% t(beta_true)
  random_effect <- matrix(rnorm(N * p, sd = between_subject_sd), N, p)
  first_level_error <- t(vapply(seq_len(N), function(i) {
    as.numeric(crossprod(chol(C_list[[i]]), rnorm(p)))
  }, numeric(p)))
  eta <- mean_matrix + random_effect + first_level_error

  mask <- matrix(FALSE, N, p)
  eligible <- if (is.null(contamination_subgroup)) seq_len(N) else which(contamination_subgroup)

  if (contamination_fraction > 0) {
    if (contamination_geometry == "cell") {
      mask[eligible, ] <- matrix(
        runif(length(eligible) * p) < contamination_fraction,
        nrow = length(eligible)
      )
      signs <- matrix(sample(c(-1, 1), N * p, replace = TRUE), N, p)
      eta[mask] <- eta[mask] + contamination_magnitude * signs[mask]
    } else {
      n_contaminated <- max(1L, round(contamination_fraction * length(eligible)))
      subjects <- sample(eligible, n_contaminated)
      mask[subjects, ] <- TRUE
      if (contamination_geometry == "subject") {
        signs <- matrix(sample(c(-1, 1), N * p, replace = TRUE), N, p)
        eta[mask] <- eta[mask] + contamination_magnitude * signs[mask]
      } else {
        signs <- sample(c(-1, 1), length(subjects), replace = TRUE)
        for (i in seq_along(subjects)) {
          eta[subjects[i], ] <- eta[subjects[i], ] +
            contamination_magnitude * signs[i]
        }
      }
    }
  }

  list(
    eta = eta,
    C_list = C_list,
    X = X,
    beta_true = beta_true,
    active = beta_true != 0,
    gain = gain,
    loss = loss,
    contamination_mask = mask
  )
}

# ---- Metrics ----------------------------------------------------------------

roc_auc <- function(score, active) {
  active <- as.logical(active)
  n1 <- sum(active)
  n0 <- sum(!active)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[active]) - n1 * (n1 + 1) / 2) /
    (n1 * n0)
}

pr_auc <- function(score, active) {
  active <- as.logical(active)
  if (!any(active)) return(NA_real_)
  order_score <- order(score, decreasing = TRUE)
  y <- active[order_score]
  tp <- cumsum(y)
  fp <- cumsum(!y)
  sum(tp / (tp + fp) * diff(c(0, tp / sum(active))))
}

sign_confidence <- function(estimate, standard_error) {
  1 - 2 * pnorm(-abs(estimate) / pmax(standard_error, 1e-12))
}

evaluate_fit <- function(fit, data, columns = NULL, threshold = 0.5) {
  truth <- data$beta_true
  active <- data$active
  estimate <- fit$estimate
  score <- fit$score
  standard_error <- fit$standard_error

  if (!is.null(columns)) {
    truth <- truth[, columns, drop = FALSE]
    active <- active[, columns, drop = FALSE]
    estimate <- estimate[, columns, drop = FALSE]
    score <- score[, columns, drop = FALSE]
    standard_error <- standard_error[, columns, drop = FALSE]
  }

  truth <- as.vector(truth)
  active <- as.vector(active)
  estimate <- as.vector(estimate)
  score <- as.vector(score)
  selected <- score > threshold

  TP <- sum(selected & active)
  FP <- sum(selected & !active)
  FN <- sum(!selected & active)
  TN <- sum(!selected & !active)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  data.frame(
    rmse = sqrt(mean((estimate - truth)^2)),
    rmse_active = if (any(active)) sqrt(mean((estimate[active] - truth[active])^2)) else NA,
    rmse_null = if (any(!active)) sqrt(mean(estimate[!active]^2)) else NA,
    pr_auc = pr_auc(score, active),
    roc_auc = roc_auc(score, active),
    tpr_at_05 = TP / max(1, TP + FN),
    fpr_at_05 = FP / max(1, FP + TN),
    fdr_at_05 = if (TP + FP) FP / (TP + FP) else 0,
    mcc_at_05 = if (denominator > 0) (TP * TN - FP * FN) / denominator else 0,
    n_selected_at_05 = sum(selected),
    fdr_at_095 = {
      strong <- score > 0.95
      if (any(strong)) mean(!active[strong]) else NA
    },
    n_selected_at_095 = sum(score > 0.95),
    runtime = fit$runtime
  )
}

score_rows <- function(fit, data, columns = NULL) {
  truth <- data$beta_true
  active <- data$active
  estimate <- fit$estimate
  standard_error <- fit$standard_error
  score <- fit$score

  if (!is.null(columns)) {
    truth <- truth[, columns, drop = FALSE]
    active <- active[, columns, drop = FALSE]
    estimate <- estimate[, columns, drop = FALSE]
    standard_error <- standard_error[, columns, drop = FALSE]
    score <- score[, columns, drop = FALSE]
  }

  data.frame(
    coefficient = seq_along(score),
    active = as.vector(active),
    truth = as.vector(truth),
    estimate = as.vector(estimate),
    standard_error = as.vector(standard_error),
    score = as.vector(score)
  )
}

contrast_probabilities <- function(fit, data, seed, n_draw = 4000L) {
  set.seed(seed)
  p <- ncol(data$eta)
  covariance <- fit$covariance[seq_len(p), seq_len(p), drop = FALSE]
  covariance <- covariance + 1e-10 * diag(p)
  draws <- matrix(fit$estimate[, 1], n_draw, p, byrow = TRUE) +
    matrix(rnorm(n_draw * p), n_draw, p) %*% chol(covariance)
  c(
    probability_gain_positive = mean(rowMeans(draws[, data$gain, drop = FALSE]) > 0),
    probability_loss_negative = mean(rowMeans(draws[, data$loss, drop = FALSE]) < 0)
  )
}

# ---- Estimators -------------------------------------------------------------
# R3.1, R3.2, R3.3: use the identified GEM/ECM fitter, integrated PIPs, and
# convergence diagnostics. R2.6 and R2.7: use multiple starts and return the
# conditional Laplace standard errors supplied by the revised fitter.

fit_pmom <- function(data, V_list, nu, pi, fit_seed, whitening_method = "precision_chol") {
  start_time <- Sys.time()
  fit <- fit_group_dcm_t_MOM(
    eta_theta_y = data$eta,
    C_theta_y_list = data$C_list,
    X_G = data$X,
    V_list = V_list,
    nu = nu,
    pi = pi,
    tau0 = 0.05,
    tau1 = 1,
    estimate_sigma2 = FALSE,
    sigma2_fixed = 1,
    whitening_method = whitening_method,
    n_starts = N_STARTS,
    start_seed = fit_seed,
    pip_mcmc_iter = PIP_ITER,
    pip_mcmc_burn = PIP_BURN,
    pip_mcmc_thin = PIP_THIN,
    pip_mcmc_chains = PIP_CHAINS,
    pip_mcmc_temperatures = PIP_TEMPERATURES,
    pip_mcmc_min_inverse_temp = 0.05,
    pip_mcmc_swap_every = 2L,
    pip_mcmc_large_jump_prob = 0.10,
    pip_mcmc_seed = fit_seed + 1L,
    max_iter = 2000L, #400L
    tol = 5e-4,
    gradient_tol = 1e-3,
    verbose = FALSE
  )

  objective_range <- diff(range(fit$multi_start_summary$objective))
  list(
    estimate = fit$beta_posterior_mean,
    posterior_mode = fit$beta_mat,
    score = fit$inclusion,
    standard_error = fit$beta_sd_approx,
    covariance = fit$beta_cov_approx,
    weights = fit$weights,
    runtime = as.numeric(Sys.time() - start_time, units = "secs"),
    diagnostics = data.frame(
      converged = fit$converged,
      stopping_reason = fit$stopping_reason,
      iterations = fit$iterations,
      final_parameter_change = fit$final_parameter_change,
      final_objective_change = fit$final_objective_change,
      alpha_optim_convergence = fit$alpha_optim_convergence,
      objective_monotone = fit$objective_monotone,
      q_monotone = fit$q_monotone,
      beta_score_max = fit$beta_score_max,
      best_start = fit$best_start,
      start_objective_range = objective_range,
      pip_mcse_max = max(fit$inclusion_mcse),
      pip_chain_range = fit$pip_max_chain_range,
      laplace_regularized = fit$laplace_regularized,
      whitening_method = fit$whitening_method
    ),
    raw_fit = fit
  )
}

ridge_objective <- function(beta, alpha, data, V_list, nu, ridge, whitening_method) {
  Sigma_b <- build_Sigma_b(alpha, V_list, ncol(data$eta))
  wd <- whiten_data(
    data$eta, data$C_list, data$X, Sigma_b,
    whitening_method = whitening_method
  )
  residual <- as.numeric(wd$y_star - wd$X_star %*% beta)
  wd$log_det_B + sum(dt(residual, df = nu, log = TRUE)) -
    0.5 * ridge * sum(beta^2)
}

fit_ridge <- function(data, V_list, nu, fit_seed, ridge = 1e-3,
                      whitening_method = "precision_chol",
                      max_iter = 2000L, tol = 5e-4,
                      alpha_max_iter = 200L) {
  start_time <- Sys.time()
  p <- ncol(data$eta)
  n_regressor <- ncol(data$X)
  d <- p * n_regressor
  K <- length(V_list)
  X_big <- kronecker(data$X, diag(p))
  y_big <- as.vector(t(data$eta))
  beta_base <- as.numeric(solve(crossprod(X_big) + 1e-6 * diag(d), crossprod(X_big, y_big)))

  set.seed(fit_seed)
  starts <- lapply(seq_len(N_STARTS), function(start) {
    beta <- if (start == 1L) beta_base else beta_base + rnorm(d, sd = 0.05)
    alpha <- if (K) rep(0.1, K) * exp(rnorm(K, sd = if (start == 1L) 0 else 0.25)) else numeric(0)
    objective_trace <- ridge_objective(beta, alpha, data, V_list, nu, ridge, whitening_method)
    converged <- FALSE
    stopping_reason <- "maximum iterations reached"
    parameter_change <- Inf
    objective_change <- Inf
    alpha_optim_convergence <- if (K > 0L) NA_integer_ else 0L

    for (iteration in seq_len(max_iter)) {
      beta_old <- beta
      alpha_old <- alpha
      objective_old <- tail(objective_trace, 1)

      Sigma_b <- build_Sigma_b(alpha, V_list, p)
      wd <- whiten_data(data$eta, data$C_list, data$X, Sigma_b,
                        whitening_method = whitening_method)
      residual <- as.numeric(wd$y_star - wd$X_star %*% beta)
      weights <- (nu + 1) / (nu + residual^2)
      beta <- as.numeric(solve(
        crossprod(wd$X_star, weights * wd$X_star) + ridge * diag(d),
        crossprod(wd$X_star, weights * wd$y_star)
      ))

      if (K) {
        alpha_objective <- function(log_alpha) {
          Sigma_try <- build_Sigma_b(exp(log_alpha), V_list, p)

          wd_try <- whiten_data(
            data$eta,
            data$C_list,
            data$X,
            Sigma_try,
            whitening_method = whitening_method
          )

          residual_try <- as.numeric(
            wd_try$y_star - wd_try$X_star %*% beta
          )

          -(
            wd_try$log_det_B -
              0.5 * sum(weights * residual_try^2) -
              0.5 * ridge * sum(beta^2)
          )
        }

        alpha_fit <- optim(
          log(pmax(alpha, 1e-8)),
          alpha_objective,
          method = "L-BFGS-B",
          lower = rep(log(1e-8), K),
          upper = rep(log(100), K),
          control = list(maxit = 500L)
        )

        if (alpha_fit$convergence == 0L && is.finite(alpha_fit$value)) {
          alpha <- exp(alpha_fit$par)
          alpha_optim_convergence <- 0L
        } else {
          alpha_retry <- nlminb(
            start = log(pmax(alpha, 1e-8)),
            objective = alpha_objective,
            lower = rep(log(1e-8), K),
            upper = rep(log(100), K),
            control = list(
              iter.max = 500L,
              eval.max = 1000L,
              rel.tol = 1e-8
            )
          )

          alpha <- exp(alpha_retry$par)
          alpha_optim_convergence <- alpha_retry$convergence
        }
      }

      objective_new <- ridge_objective(beta, alpha, data, V_list, nu, ridge, whitening_method)
      if (objective_new < objective_old - 1e-10) {
        beta <- beta_old
        alpha <- alpha_old
        stopping_reason <- "objective decreased; final iteration rolled back"
        break
      }
      objective_trace <- c(objective_trace, objective_new)

      parameter_change <- max(
        max(abs(beta - beta_old) / (1 + abs(beta_old))),
        if (K) max(abs(log(alpha / alpha_old))) else 0
      )
      objective_change <- abs(objective_new - objective_old) / (1 + abs(objective_old))
      alpha_stationary <- K == 0L || alpha_optim_convergence == 0L
      if (parameter_change < tol && objective_change < tol && alpha_stationary) {
        converged <- TRUE
        stopping_reason <- "parameter, objective, and covariance-optimizer tolerances satisfied"
        break
      }
    }

    Sigma_b <- build_Sigma_b(alpha, V_list, p)
    wd <- whiten_data(data$eta, data$C_list, data$X, Sigma_b,
                      whitening_method = whitening_method)
    residual <- as.numeric(wd$y_star - wd$X_star %*% beta)
    weights <- (nu + 1) / (nu + residual^2)
    score_vector <- as.numeric(crossprod(wd$X_star, weights * residual)) - ridge * beta

    list(
      beta = beta,
      alpha = alpha,
      weights = weights,
      wd = wd,
      objective = tail(objective_trace, 1),
      objective_trace = objective_trace,
      converged = converged,
      stopping_reason = stopping_reason,
      iterations = iteration,
      final_parameter_change = parameter_change,
      final_objective_change = objective_change,
      alpha_optim_convergence = alpha_optim_convergence,
      beta_score_max = max(abs(score_vector))
    )
  })

  objectives <- vapply(starts, `[[`, numeric(1), "objective")
  best_start <- which.max(objectives)
  best <- starts[[best_start]]
  precision <- crossprod(best$wd$X_star, best$weights * best$wd$X_star) +
    ridge * diag(d)
  covariance <- solve(precision)
  standard_error <- matrix(sqrt(diag(covariance)), p, n_regressor)
  estimate <- matrix(best$beta, p, n_regressor)

  list(
    estimate = estimate,
    score = sign_confidence(estimate, standard_error),
    standard_error = standard_error,
    covariance = covariance,
    weights = best$weights,
    alpha = best$alpha,
    whitened_y = best$wd$y_star,
    whitened_X = best$wd$X_star,
    runtime = as.numeric(Sys.time() - start_time, units = "secs"),
    diagnostics = data.frame(
      converged = best$converged,
      stopping_reason = best$stopping_reason,
      iterations = best$iterations,
      final_parameter_change = best$final_parameter_change,
      final_objective_change = best$final_objective_change,
      alpha_optim_convergence = best$alpha_optim_convergence,
      objective_monotone = all(diff(best$objective_trace) >= -sqrt(.Machine$double.eps)),
      q_monotone = NA,
      beta_score_max = best$beta_score_max,
      best_start = best_start,
      start_objective_range = diff(range(objectives)),
      pip_mcse_max = NA,
      pip_chain_range = NA,
      laplace_regularized = FALSE,
      whitening_method = whitening_method
    )
  )
}

fit_spm_peb <- function(data, prior_variance = 1, max_iterations = 32L) {
  start_time <- Sys.time()
  N <- nrow(data$eta)
  p <- ncol(data$eta)
  n_regressor <- ncol(data$X)
  prior_mean <- matrix(0, p, 1)
  prior_covariance <- diag(prior_variance, p)

  subject_fits <- lapply(seq_len(N), function(i) {
    list(
      Ep = matrix(data$eta[i, ], ncol = 1),
      Cp = data$C_list[[i]],
      F = 0,
      M = list(pE = prior_mean, pC = prior_covariance)
    )
  })
  model <- list(
    X = data$X,
    Xnames = as.list(colnames(data$X)),
    Q = "single",
    alpha = 1,
    beta = 16,
    maxit = max_iterations
  )
  peb <- spm_dcm_peb_run(
    spm_dcm_peb_prepare(P = subject_fits, M = model, field = seq_len(p)),
    verbose = FALSE
  )$PEB

  estimate <- matrix(as.numeric(peb$Ep), ncol = n_regressor)
  covariance <- as.matrix(peb$Cp)
  standard_error <- matrix(sqrt(diag(covariance)), ncol = n_regressor)
  prior_estimate <- as.numeric(peb$M$pE)
  prior_covariance <- as.matrix(peb$M$pC)

  bmr <- vapply(seq_along(peb$Ep), function(j) {
    reduced_mean <- prior_estimate
    reduced_mean[j] <- 0
    reduced_covariance <- prior_covariance
    reduced_covariance[j, ] <- 0
    reduced_covariance[, j] <- 0
    reduced_covariance[j, j] <- 1e-8
    1 / (1 + exp(spm_log_evidence_reduce(
      as.numeric(peb$Ep), covariance, prior_estimate, prior_covariance,
      reduced_mean, reduced_covariance
    )$F))
  }, numeric(1))

  list(
    estimate = estimate,
    score = matrix(bmr, ncol = n_regressor),
    standard_error = standard_error,
    covariance = covariance,
    runtime = as.numeric(Sys.time() - start_time, units = "secs"),
    diagnostics = data.frame(
      converged = TRUE,
      stopping_reason = "SPM-PEB completed",
      iterations = NA,
      final_parameter_change = NA_real_,
      final_objective_change = NA_real_,
      alpha_optim_convergence = NA_integer_,
      objective_monotone = NA,
      q_monotone = NA,
      beta_score_max = NA,
      best_start = NA,
      start_objective_range = NA,
      pip_mcse_max = NA,
      pip_chain_range = NA,
      laplace_regularized = FALSE,
      whitening_method = NA
    )
  )
}

# R2.8 and R3.6: local Gaussian spike-and-slab control.
# Conditional on the covariance
# estimate from the robust ridge fit, this model has conjugate Gaussian,
# gamma, and Bernoulli updates. The local slab variance is three times the
# pMOM scale variance, matching the variance of the first-order pMOM slab.
local_responsibility <- function(beta, tau0, tau1, pi) {
  log_spike <- log1p(-pi) + dnorm(beta, sd = tau0, log = TRUE)
  log_slab <- log(pi) + dnorm(beta, sd = tau1, log = TRUE)
  plogis(log_slab - log_spike)
}

fit_local_slab <- function(data, robust_fit, nu, pi, fit_seed) {
  start_time <- Sys.time()
  y <- robust_fit$whitened_y
  X <- robust_fit$whitened_X
  d <- ncol(X)
  tau0 <- 0.05
  tau1 <- sqrt(3)
  n_keep <- floor((PIP_ITER - PIP_BURN) / PIP_THIN)
  chain_pip <- matrix(NA_real_, PIP_CHAINS, d)
  chain_beta <- vector("list", PIP_CHAINS)
  chain_responsibility <- vector("list", PIP_CHAINS)

  for (chain in seq_len(PIP_CHAINS)) {
    set.seed(fit_seed + chain)
    beta <- as.vector(robust_fit$estimate)
    gamma <- rbinom(d, 1, if (chain == 1L) pi else runif(1, 0.1, 0.9))
    retained_beta <- matrix(NA_real_, n_keep, d)
    retained_responsibility <- matrix(NA_real_, n_keep, d)
    keep <- 0L

    for (iteration in seq_len(PIP_ITER)) {
      residual <- y - as.numeric(X %*% beta)
      lambda <- rgamma(
        length(y),
        shape = (nu + 1) / 2,
        rate = (nu + residual^2) / 2
      )
      prior_variance <- ifelse(gamma == 1L, tau1^2, tau0^2)
      precision <- crossprod(X, lambda * X) + diag(1 / prior_variance, d)
      precision_chol <- chol(precision)
      posterior_mean <- backsolve(
        precision_chol,
        forwardsolve(t(precision_chol), crossprod(X, lambda * y))
      )
      beta <- posterior_mean + backsolve(precision_chol, rnorm(d))
      responsibility <- local_responsibility(beta, tau0, tau1, pi)
      gamma <- rbinom(d, 1, responsibility)

      if (iteration > PIP_BURN &&
          (iteration - PIP_BURN) %% PIP_THIN == 0L) {
        keep <- keep + 1L
        retained_beta[keep, ] <- beta
        retained_responsibility[keep, ] <- responsibility
      }
    }

    chain_beta[[chain]] <- retained_beta
    chain_responsibility[[chain]] <- retained_responsibility
    chain_pip[chain, ] <- colMeans(retained_responsibility)
  }

  beta_draws <- do.call(rbind, chain_beta)
  responsibility_draws <- do.call(rbind, chain_responsibility)
  batch_means <- do.call(rbind, lapply(chain_responsibility, function(draws) {
    batch_size <- max(10L, floor(sqrt(nrow(draws))))
    n_batch <- floor(nrow(draws) / batch_size)
    draws <- draws[seq_len(n_batch * batch_size), , drop = FALSE]
    t(vapply(seq_len(n_batch), function(batch) {
      index <- ((batch - 1L) * batch_size + 1L):(batch * batch_size)
      colMeans(draws[index, , drop = FALSE])
    }, numeric(d)))
  }))

  p <- ncol(data$eta)
  n_regressor <- ncol(data$X)
  covariance <- cov(beta_draws)
  list(
    estimate = matrix(colMeans(beta_draws), p, n_regressor),
    score = matrix(colMeans(responsibility_draws), p, n_regressor),
    standard_error = matrix(apply(beta_draws, 2, sd), p, n_regressor),
    covariance = covariance,
    weights = robust_fit$weights,
    runtime = robust_fit$runtime + as.numeric(Sys.time() - start_time, units = "secs"),
    diagnostics = data.frame(
      converged = robust_fit$diagnostics$converged,
      stopping_reason = robust_fit$diagnostics$stopping_reason,
      iterations = robust_fit$diagnostics$iterations,
      final_parameter_change = robust_fit$diagnostics$final_parameter_change,
      final_objective_change = robust_fit$diagnostics$final_objective_change,
      alpha_optim_convergence = robust_fit$diagnostics$alpha_optim_convergence,
      objective_monotone = robust_fit$diagnostics$objective_monotone,
      q_monotone = NA,
      beta_score_max = NA,
      best_start = robust_fit$diagnostics$best_start,
      start_objective_range = robust_fit$diagnostics$start_objective_range,
      pip_mcse_max = max(apply(batch_means, 2, sd) / sqrt(nrow(batch_means))),
      pip_chain_range = max(apply(chain_pip, 2, function(x) diff(range(x)))),
      laplace_regularized = FALSE,
      whitening_method = robust_fit$diagnostics$whitening_method
    )
  )
}
# ---- Arm settings and common scenario runner --------------------------------

default_parameters <- data.frame(
  arm = c("proposed", "robust_only", "sparse_only", "local_slab"),
  nu = c(2, 2, 1e6, 2),
  pi = c(0.5, NA, 0.5, 0.5)
)

read_best_parameters <- function() {
  read.csv("sim_best_parameters.csv", stringsAsFactors = FALSE)
}

shared_parameters <- function() {
  best <- read_best_parameters()
  proposed <- best[best$arm == "proposed", ]
  data.frame(
    arm = default_parameters$arm,
    nu = c(proposed$nu, proposed$nu, 1e6, proposed$nu),
    pi = c(proposed$pi, NA, proposed$pi, proposed$pi)
  )
}

arm_parameters <- function(arm, parameter_table) {
  parameter_table[parameter_table$arm == arm, , drop = FALSE]
}

fit_all_arms <- function(data, V_list, parameter_table, fit_seed,
                         whitening_method = "precision_chol") {
  proposed_parameter <- arm_parameters("proposed", parameter_table)
  robust_parameter <- arm_parameters("robust_only", parameter_table)
  sparse_parameter <- arm_parameters("sparse_only", parameter_table)
  local_parameter <- arm_parameters("local_slab", parameter_table)

  proposed_fit <- fit_pmom(
    data, V_list, proposed_parameter$nu, proposed_parameter$pi,
    fit_seed + 10L, whitening_method
  )
  robust_fit <- fit_ridge(
    data, V_list, robust_parameter$nu,
    fit_seed + 20L, whitening_method = whitening_method
  )
  sparse_fit <- fit_pmom(
    data, V_list, 1e6, sparse_parameter$pi,
    fit_seed + 30L, whitening_method
  )
  local_robust_fit <- if (local_parameter$nu == robust_parameter$nu) {
    robust_fit
  } else {
    fit_ridge(
      data, V_list, local_parameter$nu,
      fit_seed + 35L, whitening_method = whitening_method
    )
  }
  local_fit <- fit_local_slab(
    data, local_robust_fit, local_parameter$nu, local_parameter$pi,
    fit_seed + 40L
  )

  list(
    proposed = proposed_fit,
    robust_only = robust_fit,
    sparse_only = sparse_fit,
    local_slab = local_fit,
    spm_peb = fit_spm_peb(data)
  )
}

condition_row <- function(condition, arm, rep, parameter_table) {
  # R1.1--R1.6 and R1.8: retain every condition-defining value in each result
  # row so tables, figures, captions, and nu comparisons can be cross-checked.
  parameters <- arm_parameters(arm, parameter_table)
  data.frame(
    axis = condition$axis,
    condition = condition$label,
    arm = arm,
    rep = rep,
    N = condition$N,
    p = condition$p,
    n_regressor = condition$n_regressor,
    effect_scale = condition$effect_scale,
    active_fraction = if (is.null(condition$active_fraction)) NA else condition$active_fraction,
    contamination_fraction = condition$contamination_fraction,
    contamination_magnitude = condition$contamination_magnitude,
    contamination_geometry = condition$contamination_geometry,
    covariance_regime = condition$covariance_regime,
    nu = if (nrow(parameters)) parameters$nu else NA,
    pi = if (nrow(parameters)) parameters$pi else NA
  )
}

run_condition <- function(task, conditions, parameter_table,
                          truth_seed, design_seed, data_seed) {
  condition <- conditions[[task$condition]]
  p <- condition$p
  n_regressor <- condition$n_regressor
  truth <- if (condition$sparse_truth) {
    make_sparse_truth(p, condition$active_fraction, truth_seed + task$rep)
  } else {
    make_truth(p, n_regressor, seed = truth_seed + task$rep)
  }
  truth$beta <- condition$effect_scale * truth$beta

  set.seed(design_seed + task$rep)
  X <- matrix(1, condition$N, 1, dimnames = list(NULL, "intercept"))
  if (n_regressor == 3L) {
    X <- cbind(
      intercept = 1,
      genderM = rbinom(condition$N, 1, 0.5),
      age = as.numeric(scale(rnorm(condition$N)))
    )
  }
  subgroup <- if (condition$subgroup_contamination) X[, "genderM"] == 1 else NULL

  data <- simulate_data(
    N = condition$N,
    beta_true = truth$beta,
    X = X,
    gain = truth$gain,
    loss = truth$loss,
    between_subject_sd = condition$between_subject_sd,
    covariance_regime = condition$covariance_regime,
    heterogeneous_covariance = condition$heterogeneous_covariance,
    covariance_pool = condition$covariance_pool,
    contamination_fraction = condition$contamination_fraction,
    contamination_magnitude = condition$contamination_magnitude,
    contamination_geometry = condition$contamination_geometry,
    contamination_subgroup = subgroup,
    seed = data_seed + task$rep
  )

  V_list <- if (condition$single_variance_component) {
    list(diag(p))
  } else {
    diag_basis(p)
  }
  fits <- fit_all_arms(
    data,
    V_list,
    parameter_table,
    SEED$mcmc + 10000L * task$condition + 100L * task$rep,
    condition$whitening_method
  )
  evaluation_columns <- if (condition$evaluate_covariates) 2:n_regressor else NULL

  metrics <- list()
  scores <- list()
  diagnostics <- list()
  for (arm in names(fits)) {
    metadata <- condition_row(condition, arm, task$rep, parameter_table)
    contrast <- if (n_regressor == 1L) {
      contrast_probabilities(
        fits[[arm]], data,
        SEED$mcmc + 900000L + 10000L * task$condition + 100L * task$rep
      )
    } else {
      c(probability_gain_positive = NA, probability_loss_negative = NA)
    }
    metrics[[arm]] <- cbind(
      metadata,
      evaluate_fit(fits[[arm]], data, evaluation_columns),
      probability_gain_positive = unname(contrast[1]),
      probability_loss_negative = unname(contrast[2])
    )
    scores[[arm]] <- cbind(metadata, score_rows(fits[[arm]], data, evaluation_columns))
    diagnostics[[arm]] <- cbind(metadata, fits[[arm]]$diagnostics)
  }

  list(
    metrics = do.call(rbind, metrics),
    scores = do.call(rbind, scores),
    diagnostics = do.call(rbind, diagnostics)
  )
}

run_scenarios <- function(conditions, prefix, parameter_table = shared_parameters(),
                          truth_seed = SEED$main_truth,
                          design_seed = SEED$main_design,
                          data_seed = SEED$main_data) {
  tasks <- make_tasks(length(conditions))
  cat(sprintf("[%s] %d conditions x %d replicates on %d cores\n",
              prefix, length(conditions), N_REP, N_CORE))
  results <- parallel_map(
    tasks,
    run_condition,
    conditions = conditions,
    parameter_table = parameter_table,
    truth_seed = truth_seed,
    design_seed = design_seed,
    data_seed = data_seed
  )

  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))

  write.csv(metrics, paste0(prefix, "_raw.csv"), row.names = FALSE)
  write.csv(scores, paste0(prefix, "_scores.csv"), row.names = FALSE)
  write.csv(diagnostics, paste0(prefix, "_diagnostics.csv"), row.names = FALSE)

  failed <- diagnostics[
    !is.na(diagnostics$converged) & !diagnostics$converged,
  ]

  if (nrow(failed) > 0L) {
    write.csv(
      failed,
      paste0(prefix, "_nonconverged.csv"),
      row.names = FALSE
    )
    stop(
      nrow(failed),
      " fits did not converge in ", prefix,
      ". See ", prefix, "_nonconverged.csv."
    )
  }

  write.csv(summarize_metrics(metrics), paste0(prefix, "_summary.csv"), row.names = FALSE)
  invisible(list(metrics = metrics, scores = scores, diagnostics = diagnostics))
}

summarize_metrics <- function(data) {
  data$nu_group <- ifelse(is.na(data$nu), -1, data$nu)
  data$pi_group <- ifelse(is.na(data$pi), -1, data$pi)
  data$active_fraction_group <- ifelse(is.na(data$active_fraction), -1, data$active_fraction)
  group <- c(
    "axis", "condition", "arm", "N", "p", "n_regressor",
    "effect_scale", "active_fraction_group", "contamination_fraction", "contamination_magnitude",
    "contamination_geometry", "covariance_regime", "nu_group", "pi_group"
  )
  measures <- c(
    "rmse", "rmse_active", "rmse_null", "pr_auc", "roc_auc",
    "tpr_at_05", "fpr_at_05", "fdr_at_05", "mcc_at_05",
    "n_selected_at_05", "fdr_at_095", "n_selected_at_095", "runtime",
    "probability_gain_positive", "probability_loss_negative"
  )
  summary <- aggregate(data[measures], data[group], mean, na.rm = TRUE)
  mcse <- aggregate(data[measures], data[group], function(x) {
    x <- x[is.finite(x)]
    if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_
  })
  denominators <- aggregate(data[measures], data[group], function(x) sum(is.finite(x)))
  names(mcse)[match(measures, names(mcse))] <- paste0(measures, "_mcse")
  names(denominators)[match(measures, names(denominators))] <- paste0(measures, "_n")
  summary <- merge(summary, mcse, by = group, sort = FALSE)
  summary <- merge(summary, denominators, by = group, sort = FALSE)
  names(summary)[names(summary) == "nu_group"] <- "nu"
  names(summary)[names(summary) == "pi_group"] <- "pi"
  names(summary)[names(summary) == "active_fraction_group"] <- "active_fraction"
  summary$nu[summary$nu < 0] <- NA
  summary$pi[summary$pi < 0] <- NA
  summary$active_fraction[summary$active_fraction < 0] <- NA
  summary
}

base_condition <- function() {
  list(
    axis = "baseline",
    label = "baseline_phi08_mag6",
    N = 48L,
    p = 16L,
    n_regressor = 1L,
    effect_scale = 1,
    active_fraction = NULL,
    sparse_truth = FALSE,
    contamination_fraction = 0.08,
    contamination_magnitude = 6,
    contamination_geometry = "cell",
    covariance_regime = "moderate_diagonal",
    heterogeneous_covariance = FALSE,
    covariance_pool = NULL,
    between_subject_sd = 0.15,
    subgroup_contamination = FALSE,
    evaluate_covariates = FALSE,
    single_variance_component = FALSE,
    whitening_method = "precision_chol"
  )
}

change_condition <- function(axis, label, ...) {
  modifyList(modifyList(base_condition(), list(axis = axis, label = label)), list(...))
}

# ---- Separate tuning experiment --------------------------------------------
# R2.2: tune on a deliberately perturbed generative family with different N,
# p, contamination magnitudes/fractions, heterogeneous C_n, and correlation.
# R4.6: tune the original robust-only and sparse-only ablations separately for
# the arm-specific sensitivity run; the local-slab control uses the proposed
# method's settings so that only the slab family changes.

tuning_candidates <- function() {
  nu_grid <- c(2, 3, 4, 6, 8)
  pi_grid <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
  rbind(
    transform(expand.grid(nu = nu_grid, pi = pi_grid), arm = "proposed"),
    data.frame(nu = nu_grid, pi = NA, arm = "robust_only"),
    data.frame(nu = 1e6, pi = pi_grid, arm = "sparse_only")
  )
}

tuning_conditions <- function() {
  data.frame(
    contamination_fraction = c(0, 0.06, 0.15, 0.15),
    contamination_magnitude = c(3, 3, 3, 6),
    covariance_regime = c(
      "moderate_diagonal", "correlated", "moderate_diagonal", "correlated"
    ),
    stringsAsFactors = FALSE
  )
}

run_tuning_task <- function(task, candidates, conditions) {
  candidate <- candidates[task$candidate, ]
  condition <- conditions[task$condition, ]
  p <- 18L
  N <- 54L
  truth <- make_truth(p, 1L, seed = SEED$tune_truth + task$rep)
  data <- simulate_data(
    N = N,
    beta_true = truth$beta,
    X = matrix(1, N, 1, dimnames = list(NULL, "intercept")),
    gain = truth$gain,
    loss = truth$loss,
    covariance_regime = condition$covariance_regime,
    heterogeneous_covariance = TRUE,
    contamination_fraction = condition$contamination_fraction,
    contamination_magnitude = condition$contamination_magnitude,
    contamination_geometry = "cell",
    seed = SEED$tune_data + task$rep
  )
  V_list <- diag_basis(p)
  fit_seed <- SEED$mcmc + 500000L + 10000L * task$condition +
    100L * task$rep + task$candidate

  fit <- switch(
    candidate$arm,
    proposed = fit_pmom(data, V_list, candidate$nu, candidate$pi, fit_seed),
    robust_only = fit_ridge(data, V_list, candidate$nu, fit_seed),
    sparse_only = fit_pmom(data, V_list, 1e6, candidate$pi, fit_seed)
  )
  evaluation <- evaluate_fit(fit, data)

  data.frame(
    arm = candidate$arm,
    nu = candidate$nu,
    pi = candidate$pi,
    contamination_fraction = condition$contamination_fraction,
    contamination_magnitude = condition$contamination_magnitude,
    covariance_regime = condition$covariance_regime,
    rep = task$rep,
    pr_auc = evaluation$pr_auc,
    rmse = evaluation$rmse,
    converged = fit$diagnostics$converged,
    stopping_reason = fit$diagnostics$stopping_reason,
    iterations = fit$diagnostics$iterations,
    final_parameter_change = fit$diagnostics$final_parameter_change,
    final_objective_change = fit$diagnostics$final_objective_change,
    alpha_optim_convergence = fit$diagnostics$alpha_optim_convergence,
    beta_score_max = fit$diagnostics$beta_score_max,
    pip_mcse_max = fit$diagnostics$pip_mcse_max
  )
}

normalize_to_unit_interval <- function(x) {
  if (diff(range(x)) == 0) return(rep(1, length(x)))
  (x - min(x)) / diff(range(x))
}

run_tuning <- function() {
  candidates <- tuning_candidates()
  conditions <- tuning_conditions()
  tasks <- unlist(lapply(seq_len(nrow(candidates)), function(candidate) {
    unlist(lapply(seq_len(nrow(conditions)), function(condition) {
      lapply(seq_len(N_TUNE), function(rep) {
        list(candidate = candidate, condition = condition, rep = rep)
      })
    }), recursive = FALSE)
  }), recursive = FALSE)

  cat(sprintf("[tune] %d candidates x %d perturbed conditions x %d replicates\n",
              nrow(candidates), nrow(conditions), N_TUNE))
  raw <- do.call(rbind, parallel_map(
    tasks, run_tuning_task, candidates = candidates, conditions = conditions
  ))
  write.csv(raw, "sim_tuning_raw.csv", row.names = FALSE)
  failed <- raw[!is.na(raw$converged) & !raw$converged, ]

  if (nrow(failed) > 0L) {
    write.csv(failed, "sim_tuning_nonconverged.csv", row.names = FALSE)
    stop(
      nrow(failed),
      " tuning fits did not converge. See sim_tuning_nonconverged.csv."
    )
  }

  raw$pi_group <- ifelse(is.na(raw$pi), -1, raw$pi)
  grid <- aggregate(
    cbind(pr_auc, rmse, converged) ~ arm + nu + pi_group,
    data = raw,
    FUN = mean
  )
  names(grid)[names(grid) == "pi_group"] <- "pi"
  grid$pi[grid$pi < 0] <- NA
  grid$score <- NA_real_
  for (arm in unique(grid$arm)) {
    index <- grid$arm == arm
    grid$score[index] <- 0.5 * normalize_to_unit_interval(grid$pr_auc[index]) +
      0.5 * (1 - normalize_to_unit_interval(grid$rmse[index]))
  }
  grid <- grid[order(grid$arm, -grid$score), ]
  write.csv(grid, "sim_tuning_grid.csv", row.names = FALSE)

  best <- do.call(rbind, lapply(split(grid, grid$arm), function(x) x[which.max(x$score), ]))
  local <- best[best$arm == "proposed", ]
  local$arm <- "local_slab"
  best <- rbind(best, local)
  best <- best[match(default_parameters$arm, best$arm), ]
  write.csv(best[c("arm", "nu", "pi", "pr_auc", "rmse", "score", "converged")],
            "sim_best_parameters.csv", row.names = FALSE)
  print(best)
  invisible(best)
}

# ---- Principal and supplementary simulation conditions ---------------------

main_conditions <- function() {
  conditions <- list(change_condition("baseline", "baseline_phi08_mag6"))

  for (fraction in c(0, 0.10, 0.20)) {
    for (N in c(24L, 48L, 80L)) {
      conditions[[length(conditions) + 1L]] <- change_condition(
        "contamination_sample_size",
        sprintf("phi%02d_N%d", round(100 * fraction), N),
        contamination_fraction = fraction,
        N = N
      )
    }
  }
  for (effect in c(0, 0.5, 1.5, 2)) {
    conditions[[length(conditions) + 1L]] <- change_condition(
      "effect_scale", paste0("effect_", effect), effect_scale = effect
    )
  }
  for (p in c(12L, 22L)) {
    conditions[[length(conditions) + 1L]] <- change_condition(
      "network_size", paste0("p", p), p = p
    )
  }
  for (regime in c("small_diagonal", "large_diagonal", "correlated")) {
    conditions[[length(conditions) + 1L]] <- change_condition(
      "covariance", paste0("covariance_", regime), covariance_regime = regime
    )
  }
  conditions[[length(conditions) + 1L]] <- change_condition(
    "covariate", "covariate_design", n_regressor = 3L,
    evaluate_covariates = TRUE
  )

  if (QUICK) conditions[c(1, 6, 11, 15, 17, length(conditions))] else conditions
}

moderate_contamination_conditions <- function() {
  # R2.1: include mild-to-moderate shifts of 2 and 3 in addition to shift 6.
  c(
    list(change_condition(
      "contamination_magnitude", "clean", contamination_fraction = 0,
      contamination_magnitude = 0
    )),
    lapply(c(2, 3, 6), function(magnitude) {
      change_condition(
        "contamination_magnitude",
        paste0("phi10_mag", magnitude),
        contamination_fraction = 0.10,
        contamination_magnitude = magnitude
      )
    })
  )
}

breadth_conditions <- function() {
  conditions <- list()
  for (geometry in c("cell", "subject", "structured")) {
    conditions[[length(conditions) + 1L]] <- change_condition(
      "contamination_geometry",
      paste0("geometry_", geometry),
      contamination_fraction = 0.15,
      contamination_geometry = geometry
    )
  }
  for (fraction in c(0, 0.10, 0.20)) {
    conditions[[length(conditions) + 1L]] <- change_condition(
      "covariate_contamination",
      paste0("covariate_phi", round(100 * fraction)),
      n_regressor = 3L,
      evaluate_covariates = TRUE,
      contamination_fraction = fraction
    )
  }
  conditions[[length(conditions) + 1L]] <- change_condition(
    "subgroup", "subgroup_clean", n_regressor = 3L,
    evaluate_covariates = TRUE, contamination_fraction = 0
  )
  conditions[[length(conditions) + 1L]] <- change_condition(
    "subgroup", "subgroup_genderM", n_regressor = 3L,
    evaluate_covariates = TRUE, contamination_fraction = 0.30,
    contamination_geometry = "subject", subgroup_contamination = TRUE
  )
  conditions
}

sparse_conditions <- function() {
  conditions <- list()
  for (active_fraction in c(0.05, 0.10, 0.25)) {
    for (fraction in c(0, 0.10)) {
      conditions[[length(conditions) + 1L]] <- change_condition(
        "sparse_p40",
        paste0("s", active_fraction, "_phi", round(100 * fraction)),
        p = 40L,
        active_fraction = active_fraction,
        sparse_truth = TRUE,
        contamination_fraction = fraction,
        single_variance_component = TRUE
      )
    }
  }
  conditions
}

# ---- Empirical-covariance, whitening, and objective validation --------------

run_empirical_covariance <- function() {
  # R2.3: resample first-level covariance matrices from the NARPS analysis.
  covariance_pool <- readRDS(NARPS_COVARIANCE_FILE)
  p <- nrow(covariance_pool[[1]])
  conditions <- list(
    change_condition(
      "empirical_covariance", "stylized_correlated",
      p = p,
      covariance_regime = "correlated",
      contamination_fraction = 0.10,
      contamination_magnitude = 3
    ),
    change_condition(
      "empirical_covariance", "resampled_NARPS",
      p = p,
      covariance_regime = "resampled_NARPS",
      covariance_pool = covariance_pool,
      contamination_fraction = 0.10,
      contamination_magnitude = 3
    )
  )
  run_scenarios(conditions, "sim_empirical_covariance")
}

permute_data <- function(data, permutation) {
  data$eta <- data$eta[, permutation, drop = FALSE]
  data$C_list <- lapply(data$C_list, function(C) C[permutation, permutation, drop = FALSE])
  data$beta_true <- data$beta_true[permutation, , drop = FALSE]
  data$active <- data$active[permutation, , drop = FALSE]
  data
}

run_whitening_task <- function(rep, parameter_table) {
  p <- 16L
  N <- 48L
  truth <- make_truth(p, 1L, seed = SEED$main_truth + rep)
  data <- simulate_data(
    N = N,
    beta_true = truth$beta,
    X = matrix(1, N, 1, dimnames = list(NULL, "intercept")),
    gain = truth$gain,
    loss = truth$loss,
    covariance_regime = "correlated",
    contamination_fraction = 0.10,
    contamination_magnitude = 3,
    seed = SEED$main_data + rep
  )
  proposed <- arm_parameters("proposed", parameter_table)
  permutation <- rev(seq_len(p))

  fit_precision <- fit_pmom(
    data, diag_basis(p), proposed$nu, proposed$pi,
    SEED$mcmc + 700000L + 100L * rep,
    "precision_chol"
  )
  fit_symmetric <- fit_pmom(
    data, diag_basis(p), proposed$nu, proposed$pi,
    SEED$mcmc + 710000L + 100L * rep,
    "symmetric"
  )
  permuted_data <- permute_data(data, permutation)
  fit_permuted <- fit_pmom(
    permuted_data, diag_basis(p), proposed$nu, proposed$pi,
    SEED$mcmc + 720000L + 100L * rep,
    "precision_chol"
  )
  fit_permuted$estimate <- fit_permuted$estimate[order(permutation), , drop = FALSE]
  fit_permuted$score <- fit_permuted$score[order(permutation), , drop = FALSE]

  fits <- list(
    precision_chol = fit_precision,
    symmetric = fit_symmetric,
    reversed_order = fit_permuted
  )
  metrics <- do.call(rbind, lapply(names(fits), function(specification) {
    cbind(
      data.frame(rep = rep, specification = specification),
      evaluate_fit(fits[[specification]], data)
    )
  }))
  differences <- do.call(rbind, lapply(c("symmetric", "reversed_order"), function(specification) {
    data.frame(
      rep = rep,
      specification = specification,
      maximum_absolute_estimate_difference = max(abs(
        fits[[specification]]$estimate - fit_precision$estimate
      )),
      maximum_absolute_pip_difference = max(abs(
        fits[[specification]]$score - fit_precision$score
      ))
    )
  }))
  list(metrics = metrics, differences = differences)
}

run_whitening_sensitivity <- function() {
  # R3.5 and R4.4: compare square roots and reverse the parameter ordering.
  parameter_table <- shared_parameters()
  results <- parallel_map(seq_len(N_REP), run_whitening_task,
                          parameter_table = parameter_table)
  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  differences <- do.call(rbind, lapply(results, `[[`, "differences"))
  write.csv(metrics, "sim_whitening_raw.csv", row.names = FALSE)
  write.csv(
    aggregate(metrics[c("rmse", "pr_auc")],
              list(specification = metrics$specification), mean),
    "sim_whitening_summary.csv",
    row.names = FALSE
  )
  write.csv(differences, "sim_whitening_differences.csv", row.names = FALSE)
}

finite_difference_gradient <- function(x, objective, step = 1e-5) {
  vapply(seq_along(x), function(j) {
    upper <- x
    lower <- x
    upper[j] <- upper[j] + step
    lower[j] <- lower[j] - step
    (objective(upper) - objective(lower)) / (2 * step)
  }, numeric(1))
}

validate_integrated_pip <- function() {
  # R3.1: reproduce the reviewer's scalar example and distinguish the
  # responsibility at the posterior mode from the integrated PIP.
  y <- 1
  observation_sd <- 0.5
  tau0 <- 0.05
  tau1 <- 0.5
  pi <- 0.5
  spike <- function(beta) dnorm(beta, sd = tau0)
  slab <- function(beta) (beta^2 / tau1^2) * dnorm(beta, sd = tau1)
  likelihood <- function(beta) dnorm(y, mean = beta, sd = observation_sd)
  mixture <- function(beta) (1 - pi) * spike(beta) + pi * slab(beta)
  responsibility <- function(beta) pi * slab(beta) / mixture(beta)
  grid <- seq(-3, 3, length.out = 60001L)
  mode <- grid[which.max(likelihood(grid) * mixture(grid))]
  integrated_pip <- integrate(
    function(beta) likelihood(beta) * pi * slab(beta), -Inf, Inf
  )$value / integrate(
    function(beta) likelihood(beta) * mixture(beta), -Inf, Inf
  )$value
  data.frame(
    y = y,
    observation_sd = observation_sd,
    tau0 = tau0,
    tau1 = tau1,
    pi = pi,
    posterior_mode = mode,
    responsibility_at_mode = responsibility(mode),
    integrated_pip = integrated_pip
  )
}

run_validation_task <- function(rep, parameter_table) {
  p <- 12L
  N <- 80L
  between_subject_sd <- 0.15
  true_alpha <- between_subject_sd^2
  truth <- make_truth(p, 1L, seed = SEED$validation_truth + rep)
  data <- simulate_data(
    N = N,
    beta_true = truth$beta,
    X = matrix(1, N, 1, dimnames = list(NULL, "intercept")),
    gain = truth$gain,
    loss = truth$loss,
    between_subject_sd = between_subject_sd,
    covariance_regime = "moderate_diagonal",
    contamination_fraction = 0,
    seed = SEED$validation_data + rep
  )
  V_list <- list(diag(p))
  proposed <- arm_parameters("proposed", parameter_table)
  proposed$nu <- 1e6  
  fit <- fit_pmom(
    data, V_list, proposed$nu, proposed$pi,
    SEED$mcmc + 800000L + 100L * rep
  )
  raw_fit <- fit$raw_fit

  alpha_objective <- function(log_alpha) {
    normalized_t_log_posterior(
      raw_fit$beta_vec,
      raw_fit$sigma2,
      raw_fit$tau0,
      raw_fit$tau1,
      exp(log_alpha),
      data$eta,
      data$C_list,
      data$X,
      V_list,
      proposed$nu,
      proposed$pi,
      raw_fit$a0,
      raw_fit$b0,
      raw_fit$a1,
      raw_fit$b1,
      raw_fit$whitening_method
    )
  }
  beta_objective <- function(beta) {
    normalized_t_log_posterior(
      beta,
      raw_fit$sigma2,
      raw_fit$tau0,
      raw_fit$tau1,
      raw_fit$alpha,
      data$eta,
      data$C_list,
      data$X,
      V_list,
      proposed$nu,
      proposed$pi,
      raw_fit$a0,
      raw_fit$b0,
      raw_fit$a1,
      raw_fit$b1,
      raw_fit$whitening_method
    )
  }
  alpha_gradient <- finite_difference_gradient(log(raw_fit$alpha), alpha_objective)
  beta_gradient <- finite_difference_gradient(raw_fit$beta_vec, beta_objective)

  data.frame(
    rep = rep,
    true_alpha = true_alpha,
    estimated_alpha = raw_fit$alpha,
    alpha_error = raw_fit$alpha - true_alpha,
    absolute_log_alpha_gradient = max(abs(alpha_gradient)),
    absolute_beta_gradient = max(abs(beta_gradient)),
    reported_beta_score = raw_fit$beta_score_max,
    converged = raw_fit$converged,
    objective_monotone = raw_fit$objective_monotone,
    q_monotone = raw_fit$q_monotone
  )
}

run_objective_validation <- function() {
  # R3.2 and R3.3: numerical score checks, covariance recovery, monotonicity,
  # convergence, and denominators.
  parameter_table <- shared_parameters()
  validation_replicates <- if (QUICK) 2L else 50L #100L
  results <- parallel_map(
    seq_len(validation_replicates),
    run_validation_task,
    parameter_table = parameter_table
  )
  raw <- do.call(rbind, results)
  write.csv(raw, "sim_objective_validation_raw.csv", row.names = FALSE)
  summary <- data.frame(
    n = nrow(raw),
    alpha_bias = mean(raw$alpha_error),
    alpha_rmse = sqrt(mean(raw$alpha_error^2)),
    max_log_alpha_gradient = max(raw$absolute_log_alpha_gradient),
    max_beta_gradient = max(raw$absolute_beta_gradient),
    convergence_rate = mean(raw$converged),
    objective_monotonicity_rate = mean(raw$objective_monotone),
    q_monotonicity_rate = mean(raw$q_monotone)
  )
  write.csv(summary, "sim_objective_validation_summary.csv", row.names = FALSE)
  write.csv(validate_integrated_pip(), "sim_integrated_pip_validation.csv",
            row.names = FALSE)
}

# ---- Study drivers ----------------------------------------------------------

run_main <- function() {
  run_scenarios(main_conditions(), "sim_main")
}

run_main_retry <- function() {
  required_files <- c(
    "sim_main_nonconverged.csv",
    "sim_main_raw.csv",
    "sim_main_scores.csv",
    "sim_main_diagnostics.csv",
    "sim_best_parameters.csv"
  )
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) {
    stop("Missing files required for main_retry: ", paste(missing_files, collapse = ", "))
  }

  failed <- read.csv("sim_main_nonconverged.csv", stringsAsFactors = FALSE)
  if (!nrow(failed)) stop("sim_main_nonconverged.csv contains no failed fits.")

  conditions <- main_conditions()
  condition_labels <- vapply(conditions, `[[`, character(1L), "label")
  task_keys <- unique(failed[c("axis", "condition", "rep")])
  condition_index <- match(task_keys$condition, condition_labels)
  if (anyNA(condition_index)) {
    stop(
      "Failed conditions were not found in main_conditions(): ",
      paste(task_keys$condition[is.na(condition_index)], collapse = ", ")
    )
  }

  tasks <- lapply(seq_len(nrow(task_keys)), function(i) {
    list(condition = condition_index[i], rep = as.integer(task_keys$rep[i]))
  })
  cat(sprintf("[sim_main retry] %d failed task(s) on %d cores\n", length(tasks), N_CORE))

  results <- parallel_map(
    tasks,
    run_condition,
    conditions = conditions,
    parameter_table = shared_parameters(),
    truth_seed = SEED$main_truth,
    design_seed = SEED$main_design,
    data_seed = SEED$main_data
  )
  retry_metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  retry_scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  retry_diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))

  row_key <- function(data) {
    paste(data$axis, data$condition, data$arm, data$rep, sep = "\r")
  }
  failed_keys <- unique(row_key(failed))
  retry_metrics <- retry_metrics[row_key(retry_metrics) %in% failed_keys, , drop = FALSE]
  retry_scores <- retry_scores[row_key(retry_scores) %in% failed_keys, , drop = FALSE]
  retry_diagnostics <- retry_diagnostics[
    row_key(retry_diagnostics) %in% failed_keys,
    ,
    drop = FALSE
  ]

  retry_failed <- retry_diagnostics[
    !is.na(retry_diagnostics$converged) & !retry_diagnostics$converged,
    ,
    drop = FALSE
  ]
  write.csv(retry_diagnostics, "sim_main_retry_diagnostics.csv", row.names = FALSE)
  if (nrow(retry_failed)) {
    write.csv(retry_failed, "sim_main_retry_nonconverged.csv", row.names = FALSE)
    stop(
      nrow(retry_failed),
      " retried fits still did not converge. See sim_main_retry_nonconverged.csv."
    )
  }

  replace_failed_rows <- function(old, replacement) {
    old_key <- row_key(old)
    replacement_key <- unique(row_key(replacement))
    combined <- rbind(
      old[!old_key %in% replacement_key, , drop = FALSE],
      replacement
    )
    rownames(combined) <- NULL
    combined
  }

  metrics_old <- read.csv("sim_main_raw.csv", stringsAsFactors = FALSE)
  scores_old <- read.csv("sim_main_scores.csv", stringsAsFactors = FALSE)
  diagnostics_old <- read.csv("sim_main_diagnostics.csv", stringsAsFactors = FALSE)
  metrics <- replace_failed_rows(metrics_old, retry_metrics)
  scores <- replace_failed_rows(scores_old, retry_scores)
  diagnostics <- replace_failed_rows(diagnostics_old, retry_diagnostics)

  if (nrow(metrics) != nrow(metrics_old) ||
      nrow(scores) != nrow(scores_old) ||
      nrow(diagnostics) != nrow(diagnostics_old)) {
    stop("Retry merge changed the number of rows. Existing output files were not overwritten.")
  }
  if (any(!diagnostics$converged, na.rm = TRUE)) {
    stop("Nonconverged fits remain after retry. Existing output files were not overwritten.")
  }

  write.csv(metrics, "sim_main_raw.csv", row.names = FALSE)
  write.csv(scores, "sim_main_scores.csv", row.names = FALSE)
  write.csv(diagnostics, "sim_main_diagnostics.csv", row.names = FALSE)
  write.csv(summarize_metrics(metrics), "sim_main_summary.csv", row.names = FALSE)
  cat("The failed main-simulation fits were replaced and sim_main_summary.csv was regenerated.\n")
}

run_moderate_contamination <- function() {
  run_scenarios(
    moderate_contamination_conditions(),
    "sim_moderate_contamination"
  )
}

run_breadth <- function() {
  run_scenarios(
    breadth_conditions(),
    "sim_breadth",
    truth_seed = SEED$breadth_truth,
    design_seed = SEED$breadth_design,
    data_seed = SEED$breadth_data
  )
}

run_sparse <- function() {
  run_scenarios(
    sparse_conditions(),
    "sim_sparse",
    truth_seed = SEED$sparse_truth,
    data_seed = SEED$sparse_data
  )
}

# run_nu_sensitivity <- function() {
#   conditions <- Filter(
#     function(x) x$axis %in% c("baseline", "contamination_sample_size"),
#     main_conditions()
#   )
#   parameters <- shared_parameters()
#   parameters$nu[parameters$arm %in% c("proposed", "robust_only", "local_slab")] <- 3
#   run_scenarios(conditions, "sim_nu3", parameters)
# }
run_nu_sensitivity <- function() {
  conditions <- Filter(
    function(x) {
      x$axis == "baseline" ||
        x$label %in% c("phi00_N48", "phi20_N48")
    },
    main_conditions()
  )

  parameters <- shared_parameters()
  parameters$nu[
    parameters$arm %in% c("proposed", "robust_only", "local_slab")
  ] <- 3

  run_scenarios(conditions, "sim_nu3", parameters)
}

# run_arm_tuning_sensitivity <- function() {
#   # R4.6: contrast the primary shared-tuning ablations with arm-specific tuning.
#   conditions <- Filter(
#     function(x) x$axis %in% c("baseline", "contamination_sample_size"),
#     main_conditions()
#   )
#   run_scenarios(conditions, "sim_arm_tuned", read_best_parameters())
# }
run_arm_tuning_sensitivity <- function() {
  # R4.6: contrast shared and arm-specific tuning in representative conditions.
  conditions <- Filter(
    function(x) {
      x$axis == "baseline" ||
        x$label %in% c("phi00_N48", "phi20_N48")
    },
    main_conditions()
  )

  run_scenarios(
    conditions,
    "sim_arm_tuned",
    read_best_parameters()
  )
}

# ---- Matched-error-rate and diagnostic summaries ----------------------------

matched_fpr_summary <- function(scores, target_fpr = c(0.05, 0.10)) {
  # R3.6: cross-fitted arm-specific thresholds give matched empirical FPRs and
  # report power, FDR, number of discoveries, and discovery frequency.
  scores$fold <- scores$rep %% 2L
  output <- list()
  output_index <- 0L

  for (arm in unique(scores$arm)) {
    arm_scores <- scores[scores$arm == arm, ]
    for (target in target_fpr) {
      selected_rows <- list()
      for (fold in 0:1) {
        calibration_null <- arm_scores$score[
          arm_scores$fold != fold & !arm_scores$active
        ]
        threshold <- as.numeric(quantile(
          calibration_null,
          probs = 1 - target,
          type = 1,
          names = FALSE
        ))
        evaluation <- arm_scores[arm_scores$fold == fold, ]
        evaluation$selected <- evaluation$score > threshold
        evaluation$threshold <- threshold
        selected_rows[[fold + 1L]] <- evaluation
      }
      evaluation <- do.call(rbind, selected_rows)
      TP <- sum(evaluation$selected & evaluation$active)
      FP <- sum(evaluation$selected & !evaluation$active)
      FN <- sum(!evaluation$selected & evaluation$active)
      TN <- sum(!evaluation$selected & !evaluation$active)
      unit <- interaction(evaluation$condition, evaluation$rep, drop = TRUE)
      number_selected <- tapply(evaluation$selected, unit, sum)

      output_index <- output_index + 1L
      output[[output_index]] <- data.frame(
        arm = arm,
        target_fpr = target,
        observed_fpr = FP / (FP + TN),
        power = TP / (TP + FN),
        fdr = if (TP + FP) FP / (TP + FP) else 0,
        mean_number_selected = mean(number_selected),
        discovery_frequency = mean(number_selected > 0),
        number_of_coefficients = nrow(evaluation),
        number_of_replicate_conditions = length(number_selected)
      )
    }
  }
  do.call(rbind, output)
}

paired_summary <- function(data, measures) {
  # R3.9: average conditions within replicate first, then compute the
  # Monte Carlo standard error across the matched replicate-level averages.
  do.call(rbind, lapply(split(data, data$arm), function(arm_data) {
    do.call(rbind, lapply(measures, function(measure) {
      replicate_mean <- tapply(arm_data[[measure]], arm_data$rep, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) mean(x) else NA_real_
      })
      replicate_mean <- replicate_mean[is.finite(replicate_mean)]
      data.frame(
        arm = arm_data$arm[1],
        measure = measure,
        mean = mean(replicate_mean),
        mcse = sd(replicate_mean) / sqrt(length(replicate_mean)),
        n_replicates = length(replicate_mean)
      )
    }))
  }))
}

condition_table <- function(conditions) {
  do.call(rbind, lapply(conditions, function(x) {
    data.frame(
      axis = x$axis,
      condition = x$label,
      N = x$N,
      p = x$p,
      n_regressor = x$n_regressor,
      effect_scale = x$effect_scale,
      active_fraction = if (is.null(x$active_fraction)) NA else x$active_fraction,
      contamination_fraction = x$contamination_fraction,
      contamination_magnitude = x$contamination_magnitude,
      contamination_geometry = x$contamination_geometry,
      covariance_regime = x$covariance_regime,
      heterogeneous_covariance = x$heterogeneous_covariance,
      whitening_method = x$whitening_method
    )
  }))
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

run_report <- function() {
  main_raw <- read.csv("sim_main_raw.csv", stringsAsFactors = FALSE)
  main_scores <- read.csv("sim_main_scores.csv", stringsAsFactors = FALSE)
  headline_raw <- main_raw[
    main_raw$axis %in% c("baseline", "contamination_sample_size"),
  ]
  headline_scores <- main_scores[
    main_scores$axis %in% c("baseline", "contamination_sample_size"),
  ]
  write.csv(
    paired_summary(
      headline_raw,
      c("rmse", "pr_auc", "mcc_at_05", "fdr_at_095", "runtime")
    ),
    "sim_headline_paired_summary.csv",
    row.names = FALSE
  )
  write.csv(
    matched_fpr_summary(headline_scores),
    "sim_matched_fpr.csv",
    row.names = FALSE
  )

  diagnostics <- read.csv("sim_main_diagnostics.csv", stringsAsFactors = FALSE)
  diagnostics_by_arm <- do.call(rbind, lapply(split(diagnostics, diagnostics$arm), function(data) {
    data.frame(
      arm = data$arm[1],
      n_fits = nrow(data),
      n_failed = sum(!data$converged, na.rm = TRUE),
      n_converged = sum(data$converged, na.rm = TRUE),
      convergence_rate = mean(data$converged, na.rm = TRUE),
      objective_monotonicity_rate = mean(data$objective_monotone, na.rm = TRUE),
      q_monotonicity_rate = mean(data$q_monotone, na.rm = TRUE),
      maximum_beta_score = safe_max(data$beta_score_max),
      maximum_pip_mcse = safe_max(data$pip_mcse_max),
      maximum_between_chain_pip_range = safe_max(data$pip_chain_range),
      laplace_regularization_rate = mean(data$laplace_regularized, na.rm = TRUE)
    )
  }))
  write.csv(diagnostics_by_arm, "sim_fit_diagnostics_summary.csv", row.names = FALSE)
  write.csv(condition_table(main_conditions()), "sim_main_conditions.csv", row.names = FALSE)

  main <- read.csv("sim_main_summary.csv", stringsAsFactors = FALSE)
  moderate <- read.csv("sim_moderate_contamination_summary.csv", stringsAsFactors = FALSE)
  colors <- c(
    proposed = "#1b7837",
    robust_only = "#4393c3",
    sparse_only = "#d6604d",
    local_slab = "#8c510a",
    spm_peb = "#7f7f7f"
  )
  labels <- c(
    proposed = "Proposed",
    robust_only = "Robust-only",
    sparse_only = "Sparse-only",
    local_slab = "Local-slab",
    spm_peb = "SPM PEB"
  )

  plot_panel <- function(data, x, y, xlab, ylab, tag, ylim = NULL) {
    x_values <- sort(unique(data[[x]]))
    y_values <- data[[y]][is.finite(data[[y]])]
    if (is.null(ylim)) ylim <- range(y_values)
    plot(NA, xlim = range(x_values), ylim = ylim, xlab = xlab, ylab = ylab)
    grid()
    for (arm in names(colors)) {
      arm_data <- data[data$arm == arm, ]
      means <- tapply(arm_data[[y]], arm_data[[x]], mean, na.rm = TRUE)
      lines(as.numeric(names(means)), means, type = "b", pch = 19,
            lwd = 2, col = colors[arm])
    }
    mtext(tag, side = 3, adj = 0, font = 2)
  }

  grid_data <- main[main$axis == "contamination_sample_size" & main$N == 48, ]
  magnitude_data <- moderate[moderate$axis == "contamination_magnitude", ]
  png("fig_simulation_revision.png", width = 1800, height = 900, res = 180)
  par(mfrow = c(1, 3), mar = c(4.2, 4.5, 2.2, 1))
  plot_panel(
    grid_data, "contamination_fraction", "pr_auc",
    "contamination fraction", "PR-AUC", "(a)", c(0, 1)
  )
  plot_panel(
    magnitude_data, "contamination_magnitude", "rmse",
    "contamination magnitude", "RMSE", "(b)"
  )
  plot.new()
  legend("center", labels, col = colors, lwd = 2, pch = 19,
         bty = "n", title = "Estimator")
  dev.off()
}

# ---- Dispatch ---------------------------------------------------------------

write_run_settings <- function(output_file = "sim_run_settings.csv") {
  # R3.9: record the exact code settings and disjoint seed blocks used.
  settings <- c(
    script = "simulation_study.R",
    method_source = "rsDCM.R",
    spm_source = "simulation/spm_dcm_peb_port.R",
    r_version = R.version.string,
    stage = STAGE,
    quick = QUICK,
    replicates = N_REP,
    tuning_replicates = N_TUNE,
    cores = N_CORE,
    parallel_backend = "snowfall MPI",
    results_directory = if (QUICK) "simulation/results_quick" else "simulation/results",
    starts = N_STARTS,
    pip_iterations = PIP_ITER,
    pip_burn = PIP_BURN,
    pip_thin = PIP_THIN,
    pip_chains = PIP_CHAINS,
    pip_temperatures = PIP_TEMPERATURES
  )
  settings <- c(settings, setNames(unlist(SEED), paste0("seed_", names(SEED))))
  write.csv(
    data.frame(setting = names(settings), value = as.character(settings)),
    output_file,
    row.names = FALSE
  )
}

cat(sprintf("stage=%s, replicates=%d, tuning replicates=%d, cores=%d\n",
            STAGE, N_REP, N_TUNE, N_CORE))
write_run_settings(paste0("sim_run_settings_", STAGE, ".csv"))

if (STAGE != "report") {
  suppressWarnings(sfInit(parallel = N_CORE > 1L, cpus = N_CORE, type = "MPI"))
  sfExportAll()

  # MPI workers are separate R sessions. Source the Rcpp-enabled method file
  # on every worker so that its compiled native routine is loaded locally.
  worker_ready <- unlist(sfClusterCall(
    function(method_file, spm_file) {
      source(method_file, local = .GlobalEnv)
      source(spm_file, local = .GlobalEnv, chdir = TRUE)

      exists(
        ".get_mom_tempering_cpp",
        envir = .GlobalEnv,
        mode = "function",
        inherits = FALSE
      ) &&
        exists(
          ".spmat",
          envir = .GlobalEnv,
          mode = "function",
          inherits = FALSE
        ) &&
        exists(
          "spm_dcm_peb_run",
          envir = .GlobalEnv,
          mode = "function",
          inherits = FALSE
        )
    },
    METHOD_FILE,
    SPM_FILE
  ), use.names = FALSE)

  if (!all(worker_ready)) {
    stop(
      "Required method and SPM routines were not loaded on every MPI worker: ",
      sum(worker_ready), " of ", length(worker_ready), " workers ready."
    )
  }

  cat(sprintf(
    "Method and SPM routines loaded on all %d MPI workers.\n",
    length(worker_ready)
  ))
}

switch(
  STAGE,
  tune = run_tuning(),
  main = run_main(),
  main_retry = run_main_retry(),
  moderate = run_moderate_contamination(),
  breadth = run_breadth(),
  sparse = run_sparse(),
  empirical = run_empirical_covariance(),
  whitening = run_whitening_sensitivity(),
  validation = run_objective_validation(),
  nu = run_nu_sensitivity(),
  fairness = run_arm_tuning_sensitivity(),
  report = run_report(),
  all = {
    run_tuning()
    run_main()
    run_moderate_contamination()
    run_breadth()
    run_sparse()
    if (file.exists(NARPS_COVARIANCE_FILE)) {
      run_empirical_covariance()
    } else {
      cat("[empirical] skipped because simulation/narps_covariances.rds is unavailable\n")
    }
    run_whitening_sensitivity()
    run_objective_validation()
    run_nu_sensitivity()
    run_arm_tuning_sensitivity()
    run_report()
  },
  stop("Unknown stage: ", STAGE)
)

if (STAGE != "report") {
  sfStop()
  Rmpi::mpi.quit(save = "no")
}
cat("Done.\n")
