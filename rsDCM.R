# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Godfred Arhin and Nilotpal Sanyal

# ---- Implementation of novel Group DCM with t and MOM ---- #

# Rcpp acceleration for the full parallel-tempering sampler. The Student-t
# log density uses its equivalent closed form with constants precomputed once,
# and residuals are stored contiguously by temperature. Compilation is lazy and
# occurs once per R session. Requires Rcpp and a working C++ compiler. The
# empirical-Bayes fitting algorithm and simulated example below are unchanged.
.get_mom_tempering_cpp <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      if (!requireNamespace("Rcpp", quietly = TRUE)) {
        stop("Package 'Rcpp' is required for the accelerated sampler.")
      }
      cpp_code <- paste(c(
        "#include <Rcpp.h>",
        "#include <algorithm>",
        "#include <cmath>",
        "#include <limits>",
        "#include <vector>",
        "using namespace Rcpp;",
        "",
        "static inline double mom_log_mix(const double beta, const double tau0,",
        "                                 const double tau1, const double pi) {",
        "  const double log_f0 = R::dnorm(beta, 0.0, tau0, true);",
        "  const double log_f1 = (beta == 0.0) ? R_NegInf :",
        "    -0.5 * std::log(2.0 * std::acos(-1.0)) - 3.0 * std::log(tau1) +",
        "    2.0 * std::log(std::fabs(beta)) - 0.5 * beta * beta / (tau1 * tau1);",
        "  const double a = std::log1p(-pi) + log_f0;",
        "  const double b = std::log(pi) + log_f1;",
        "  if (a >= b) return a + std::log1p(std::exp(b - a));",
        "  return b + std::log1p(std::exp(a - b));",
        "}",
        "",
        "static inline double mom_resp(const double beta, const double tau0,",
        "                              const double tau1, const double pi) {",
        "  const double log_f0 = R::dnorm(beta, 0.0, tau0, true);",
        "  const double log_f1 = (beta == 0.0) ? R_NegInf :",
        "    -0.5 * std::log(2.0 * std::acos(-1.0)) - 3.0 * std::log(tau1) +",
        "    2.0 * std::log(std::fabs(beta)) - 0.5 * beta * beta / (tau1 * tau1);",
        "  const double d = (std::log(pi) + log_f1) - (std::log1p(-pi) + log_f0);",
        "  if (d >= 0.0) return 1.0 / (1.0 + std::exp(-d));",
        "  const double ed = std::exp(d);",
        "  return ed / (1.0 + ed);",
        "}",
        "",
        "static inline double t_log_density_residual(const double residual,",
        "                                             const double inv_nu_sigma2,",
        "                                             const double log_constant,",
        "                                             const double exponent) {",
        "  const double scaled_square = residual * residual * inv_nu_sigma2;",
        "  return log_constant - exponent * std::log1p(scaled_square);",
        "}",
        "",
        "// [[Rcpp::export]]",
        "List mom_parallel_tempering_cpp(NumericVector beta_init,",
        "                                NumericVector y_star,",
        "                                NumericMatrix X_star,",
        "                                NumericVector inverse_temperatures,",
        "                                NumericVector coordinate_sd,",
        "                                double nu, double sigma2,",
        "                                double tau0, double tau1, double pi,",
        "                                int n_iter, int burn, int thin,",
        "                                int n_chains, int n_temperatures,",
        "                                double large_jump_prob, int swap_every) {",
        "  const int q = beta_init.size();",
        "  const int n = y_star.size();",
        "  const int nt = n_temperatures;",
        "  const int nkeep = (n_iter - burn) / thin;",
        "  const double inv_nu_sigma2 = 1.0 / (nu * sigma2);",
        "  const double log_constant = std::lgamma(0.5 * (nu + 1.0)) -",
        "    std::lgamma(0.5 * nu) - 0.5 * std::log(nu * std::acos(-1.0) * sigma2);",
        "  const double t_exponent = 0.5 * (nu + 1.0);",
        "  Function sample_int = Environment::base_env()[\"sample.int\"];",
        "",
        "  List pip_draws(n_chains);",
        "  NumericMatrix beta_means(n_chains, q);",
        "  NumericMatrix local_acceptance(n_chains, nt);",
        "  NumericMatrix large_jump_acceptance(n_chains, nt);",
        "  NumericMatrix swap_acceptance(n_chains, nt - 1);",
        "  NumericMatrix proposal_scales(n_chains, nt);",
        "",
        "  for (int chain = 0; chain < n_chains; ++chain) {",
        "    NumericMatrix states(nt, q);",
        "    for (int t = 0; t < nt; ++t)",
        "      for (int j = 0; j < q; ++j) states(t, j) = beta_init[j];",
        "",
        "    for (int t = 0; t < nt; ++t) {",
        "      if (chain > 0 || t > 0) {",
        "        const double slab_probability = std::max(",
        "          static_cast<double>(chain) / n_chains,",
        "          static_cast<double>(t) / nt",
        "        );",
        "        std::vector<int> slab(q, 0);",
        "        for (int j = 0; j < q; ++j)",
        "          slab[j] = R::runif(0.0, 1.0) < slab_probability;",
        "        for (int j = 0; j < q; ++j)",
        "          if (slab[j]) states(t, j) += R::rnorm(0.0, tau1);",
        "        for (int j = 0; j < q; ++j)",
        "          if (!slab[j]) states(t, j) += R::rnorm(0.0, tau0);",
        "      }",
        "    }",
        "",
        "    NumericMatrix residuals(n, nt);",
        "    NumericMatrix prior_terms(nt, q);",
        "    NumericVector log_likelihoods(nt), prior_sums(nt), log_posteriors(nt);",
        "    for (int t = 0; t < nt; ++t) {",
        "      for (int i = 0; i < n; ++i) {",
        "        double fitted = 0.0;",
        "        for (int j = 0; j < q; ++j) fitted += X_star(i, j) * states(t, j);",
        "        residuals(i, t) = y_star[i] - fitted;",
        "      }",
        "      long double ll_acc = 0.0L;",
        "      for (int i = 0; i < n; ++i) {",
        "        ll_acc += static_cast<long double>(t_log_density_residual(",
        "          residuals(i, t), inv_nu_sigma2, log_constant, t_exponent));",
        "      }",
        "      log_likelihoods[t] = static_cast<double>(ll_acc);",
        "      double psum = 0.0;",
        "      for (int j = 0; j < q; ++j) {",
        "        prior_terms(t, j) = mom_log_mix(states(t, j), tau0, tau1, pi);",
        "        psum += prior_terms(t, j);",
        "      }",
        "      prior_sums[t] = psum;",
        "      log_posteriors[t] = log_likelihoods[t] + psum;",
        "    }",
        "",
        "    NumericVector proposal_scale(nt);",
        "    std::fill(proposal_scale.begin(), proposal_scale.end(), 2.4);",
        "    std::vector<int> accepted_postburn(nt, 0), proposed_postburn(nt, 0);",
        "    std::vector<int> accepted_window(nt, 0), proposed_window(nt, 0);",
        "    std::vector<int> jump_accepted(nt, 0), jump_proposed(nt, 0);",
        "    std::vector<int> swap_accepted(nt - 1, 0), swap_proposed(nt - 1, 0);",
        "    int swap_phase = 0;",
        "    int kept = 0;",
        "    NumericMatrix pip_chain(nkeep, q);",
        "    NumericVector beta_sum(q);",
        "",
        "    for (int iteration = 1; iteration <= n_iter; ++iteration) {",
        "      for (int t = 0; t < nt; ++t) {",
        "        IntegerVector order = sample_int(q);",
        "        for (int pos = 0; pos < q; ++pos) {",
        "          const int j = order[pos] - 1;",
        "          const bool is_large_jump = R::runif(0.0, 1.0) < large_jump_prob;",
        "          const double temp_scale = 1.0 / std::sqrt(inverse_temperatures[t]);",
        "          const double step_sd = is_large_jump ? tau1 * temp_scale :",
        "            proposal_scale[t] * coordinate_sd[j] * temp_scale;",
        "          const double proposed_value = states(t, j) + R::rnorm(0.0, step_sd);",
        "          const double delta = proposed_value - states(t, j);",
        "          long double proposed_ll_acc = 0.0L;",
        "          for (int i = 0; i < n; ++i) {",
        "            const double proposed_resid = residuals(i, t) - X_star(i, j) * delta;",
        "            proposed_ll_acc += static_cast<long double>(t_log_density_residual(",
        "              proposed_resid, inv_nu_sigma2, log_constant, t_exponent));",
        "          }",
        "          const double proposed_log_likelihood = static_cast<double>(proposed_ll_acc);",
        "          const double proposed_log_prior_j = mom_log_mix(proposed_value, tau0, tau1, pi);",
        "          const double proposed_log_prior_sum = prior_sums[t] - prior_terms(t, j) + proposed_log_prior_j;",
        "          const double proposed_log_posterior = proposed_log_likelihood + proposed_log_prior_sum;",
        "          const double log_ratio = inverse_temperatures[t] *",
        "            (proposed_log_posterior - log_posteriors[t]);",
        "          const bool accepted = std::isfinite(proposed_log_posterior) &&",
        "            std::log(R::runif(0.0, 1.0)) < log_ratio;",
        "          if (accepted) {",
        "            states(t, j) = proposed_value;",
        "            for (int i = 0; i < n; ++i)",
        "              residuals(i, t) -= X_star(i, j) * delta;",
        "            log_likelihoods[t] = proposed_log_likelihood;",
        "            prior_terms(t, j) = proposed_log_prior_j;",
        "            prior_sums[t] = proposed_log_prior_sum;",
        "            log_posteriors[t] = proposed_log_posterior;",
        "          }",
        "",
        "          if (iteration > burn) {",
        "            ++proposed_postburn[t];",
        "            if (accepted) ++accepted_postburn[t];",
        "            if (is_large_jump) {",
        "              ++jump_proposed[t];",
        "              if (accepted) ++jump_accepted[t];",
        "            }",
        "          } else if (!is_large_jump) {",
        "            ++proposed_window[t];",
        "            if (accepted) ++accepted_window[t];",
        "          }",
        "        }",
        "      }",
        "",
        "      if (iteration <= burn && iteration % 50 == 0) {",
        "        for (int t = 0; t < nt; ++t) {",
        "          const double rate = static_cast<double>(accepted_window[t]) /",
        "            std::max(proposed_window[t], 1);",
        "          const double adj = std::max(-0.5, std::min(0.5, rate - 0.44));",
        "          proposal_scale[t] *= std::exp(adj);",
        "          accepted_window[t] = proposed_window[t] = 0;",
        "        }",
        "      }",
        "",
        "      if (iteration % swap_every == 0) {",
        "        swap_phase = 1 - swap_phase;",
        "        const int first_pair = swap_phase == 0 ? 0 : 1;",
        "        for (int lower = first_pair; lower < nt - 1; lower += 2) {",
        "          const int upper = lower + 1;",
        "          const double log_swap_ratio =",
        "            (inverse_temperatures[lower] - inverse_temperatures[upper]) *",
        "            (log_posteriors[upper] - log_posteriors[lower]);",
        "          if (iteration > burn) ++swap_proposed[lower];",
        "          if (std::log(R::runif(0.0, 1.0)) < log_swap_ratio) {",
        "            for (int j = 0; j < q; ++j) std::swap(states(lower, j), states(upper, j));",
        "            for (int i = 0; i < n; ++i) std::swap(residuals(i, lower), residuals(i, upper));",
        "            for (int j = 0; j < q; ++j) std::swap(prior_terms(lower, j), prior_terms(upper, j));",
        "            std::swap(log_likelihoods[lower], log_likelihoods[upper]);",
        "            std::swap(prior_sums[lower], prior_sums[upper]);",
        "            std::swap(log_posteriors[lower], log_posteriors[upper]);",
        "            if (iteration > burn) ++swap_accepted[lower];",
        "          }",
        "        }",
        "      }",
        "",
        "      if (iteration > burn && (iteration - burn) % thin == 0) {",
        "        for (int j = 0; j < q; ++j) {",
        "          pip_chain(kept, j) = mom_resp(states(0, j), tau0, tau1, pi);",
        "          beta_sum[j] += states(0, j);",
        "        }",
        "        ++kept;",
        "      }",
        "    }",
        "",
        "    pip_draws[chain] = pip_chain;",
        "    for (int j = 0; j < q; ++j) beta_means(chain, j) = beta_sum[j] / nkeep;",
        "    for (int t = 0; t < nt; ++t) {",
        "      local_acceptance(chain, t) = proposed_postburn[t] > 0 ?",
        "        static_cast<double>(accepted_postburn[t]) / proposed_postburn[t] : NA_REAL;",
        "      large_jump_acceptance(chain, t) = jump_proposed[t] > 0 ?",
        "        static_cast<double>(jump_accepted[t]) / jump_proposed[t] : NA_REAL;",
        "      proposal_scales(chain, t) = proposal_scale[t];",
        "    }",
        "    for (int t = 0; t < nt - 1; ++t)",
        "      swap_acceptance(chain, t) = swap_proposed[t] > 0 ?",
        "        static_cast<double>(swap_accepted[t]) / swap_proposed[t] : NA_REAL;",
        "  }",
        "",
        "  return List::create(",
        "    _[\"pip_draws\"] = pip_draws,",
        "    _[\"beta_means\"] = beta_means,",
        "    _[\"local_acceptance\"] = local_acceptance,",
        "    _[\"large_jump_acceptance\"] = large_jump_acceptance,",
        "    _[\"swap_acceptance\"] = swap_acceptance,",
        "    _[\"proposal_scales\"] = proposal_scales",
        "  );",
        "}"
      ), collapse = "\n")
      compile_env <- new.env(parent = globalenv())
      Rcpp::sourceCpp(
        code = cpp_code, env = compile_env, rebuild = FALSE,
        showOutput = FALSE, verbose = FALSE
      )
      cached <<- compile_env$mom_parallel_tempering_cpp
    }
    cached
  }
})

# 1. MOM spike-and-slab prior helper

prior_components_MOM <- function(beta, tau0, tau1, pi) {
  # Properly NORMALIZED spike and slab densities. The normalizing constants are
  # essential here: (i) the two components are combined additively, so their
  # relative heights matter, and (ii) tau0 and tau1 are learned, so their
  # tau-dependent normalizers must enter the scale-update objective.
  # base::pi is the constant 3.14159..., written explicitly because the argument
  # `pi` (the inclusion probability) shadows it in this function.
  f0 <- (2 * base::pi * tau0^2)^(-0.5) * exp(-0.5 * (beta^2 / tau0^2))  # spike: N(0, tau0^2)
  f1 <- (2 * base::pi)^(-0.5) * tau1^(-3) * (beta^2) * exp(-0.5 * (beta^2 / tau1^2))  # 1st-order pMOM slab

  m  <- (1 - pi) * f0 + pi * f1

  # First derivatives
  f0_prime <- - (beta / tau0^2) * f0

  # Numerically robust treatment near zero
  eps <- 1e-8
  # sign(0) == 0, which would leave a zero in the denominator; map exact zeros
  # (and sub-eps magnitudes) to +/- eps instead.
  beta_safe <- ifelse(abs(beta) < eps, ifelse(beta >= 0, eps, -eps), beta)

  f1_prime <- f1 * (2 / beta_safe - beta_safe / tau1^2)

  m_prime <- (1 - pi) * f0_prime + pi * f1_prime

  # Second derivatives
  f0_second <- ((beta^2 / tau0^4) - (1 / tau0^2)) * f0

  logf1_second <- -2 / (beta_safe^2) - 1 / tau1^2
  logf1_prime  <-  2 / beta_safe - beta_safe / tau1^2
  f1_second    <- f1 * (logf1_prime^2 + logf1_second)

  m_second <- (1 - pi) * f0_second + pi * f1_second

  # Derivatives of log m
  logm_prime  <- m_prime / (m + 1e-16)
  logm_second <- (m_second * m - m_prime^2) / (m^2 + 1e-16)

  list(
    f0          = f0,
    f1          = f1,
    m           = m,
    logm_prime  = logm_prime,
    logm_second = logm_second
  )
}

# Stable log densities and responsibilities for the Gaussian-spike/pMOM-slab
# mixture. These are used by the common objective and the integrated-PIP step.
log_mixture_MOM <- function(beta, tau0, tau1, pi) {
  log_f0 <- dnorm(beta, mean = 0, sd = tau0, log = TRUE)
  log_f1 <- -0.5 * log(2 * base::pi) - 3 * log(tau1) +
    2 * log(abs(beta)) - 0.5 * beta^2 / tau1^2
  log_a <- log1p(-pi) + log_f0
  log_b <- log(pi) + log_f1
  log_m <- ifelse(
    log_a >= log_b,
    log_a + log1p(exp(log_b - log_a)),
    log_b + log1p(exp(log_a - log_b))
  )
  list(log_m = log_m, responsibility = plogis(log_b - log_a))
}

expected_mixture_log_prior_MOM <- function(
    beta, responsibility, tau0, tau1, pi
) {
  responsibility <- pmin(pmax(as.numeric(responsibility), 0), 1)
  log_f0 <- dnorm(beta, mean = 0, sd = tau0, log = TRUE)
  out <- (1 - responsibility) * (log1p(-pi) + log_f0)
  slab_index <- responsibility > 0
  if (any(slab_index)) {
    beta_slab <- beta[slab_index]
    log_f1 <- -0.5 * log(2 * base::pi) - 3 * log(tau1) +
      2 * log(abs(beta_slab)) - 0.5 * beta_slab^2 / tau1^2
    out[slab_index] <- out[slab_index] +
      responsibility[slab_index] * (log(pi) + log_f1)
  }
  sum(out)
}

log_tau_hyperprior <- function(tau0, tau1, a0, b0, a1, b1) {
  tau0_sq <- tau0^2
  tau1_sq <- tau1^2
  a0 * log(b0) - lgamma(a0) - (a0 + 1) * log(tau0_sq) - b0 / tau0_sq +
    a1 * log(b1) - lgamma(a1) - (a1 + 1) * log(tau1_sq) - b1 / tau1_sq
}


# Parallel-tempering posterior sampling for integrated inclusion probabilities.
# Hyperparameters and covariance quantities are fixed at their empirical-Bayes
# estimates. Each independent ladder contains a cold replica targeting the desired
# posterior and heated replicas that cross separated spike and slab modes more
# readily. Adjacent replicas exchange states by Metropolis swaps. At each retained
# cold-replica draw, the conditional slab responsibility is evaluated and averaged,
# giving the Rao-Blackwellized estimate E{Pr(gamma_j = 1 | beta_j, h) | y, h}.
sample_integrated_pip_MOM <- function(
    beta_init,
    y_star,
    X_star,
    nu,
    sigma2,
    tau0,
    tau1,
    pi,
    weights_at_mode,
    n_iter = 6000L,
    burn = 2000L,
    thin = 2L,
    n_chains = 4L,
    n_temperatures = 8L,
    min_inverse_temp = 0.05,
    swap_every = 2L,
    large_jump_prob = 0.10,
    seed = NULL,
    verbose = FALSE
) {
  beta_init <- as.numeric(beta_init)
  y_star <- as.numeric(y_star)
  X_star <- as.matrix(X_star)
  q <- length(beta_init)

  n_iter <- as.integer(n_iter)
  burn <- as.integer(burn)
  thin <- as.integer(thin)
  n_chains <- as.integer(n_chains)
  n_temperatures <- as.integer(n_temperatures)
  swap_every <- as.integer(swap_every)
  if (n_iter <= burn || burn < 0L || thin < 1L || n_chains < 2L) {
    stop("Use n_iter > burn >= 0, thin >= 1, and at least two MCMC chains.")
  }
  if (n_temperatures < 2L || swap_every < 1L || !is.finite(min_inverse_temp) ||
      min_inverse_temp <= 0 || min_inverse_temp >= 1) {
    stop("Use at least two temperatures, swap_every >= 1, and min_inverse_temp in (0, 1).")
  }
  if (!is.finite(large_jump_prob) || large_jump_prob < 0 || large_jump_prob >= 1) {
    stop("large_jump_prob must lie in [0, 1).")
  }
  if (n_iter - burn < 2L * swap_every) {
    stop("The post-burn period must contain at least two swap intervals.")
  }
  if (!is.finite(pi) || pi <= 0 || pi >= 1 || tau0 <= 0 || tau1 <= 0 || sigma2 <= 0) {
    stop("pi must lie in (0, 1), and tau0, tau1, and sigma2 must be positive.")
  }
  if (!is.null(seed)) set.seed(seed)

  log_mix_components <- function(beta) {
    log_f0 <- dnorm(beta, mean = 0, sd = tau0, log = TRUE)
    log_f1 <- -0.5 * log(2 * base::pi) - 3 * log(tau1) +
      2 * log(abs(beta)) - 0.5 * beta^2 / tau1^2
    log_a <- log1p(-pi) + log_f0
    log_b <- log(pi) + log_f1
    log_m <- ifelse(
      log_a >= log_b,
      log_a + log1p(exp(log_b - log_a)),
      log_b + log1p(exp(log_a - log_b))
    )
    list(log_m = log_m, responsibility = plogis(log_b - log_a))
  }

  log_likelihood_residual <- function(residual) {
    sum(
      dt(residual / sqrt(sigma2), df = nu, log = TRUE) - 0.5 * log(sigma2)
    )
  }

  # Use the local weighted curvature only to scale the random-walk proposal. It
  # does not replace posterior sampling or enter the PIP calculation itself.
  pc_mode <- prior_components_MOM(beta_init, tau0 = tau0, tau1 = tau1, pi = pi)
  prior_precision <- pmax(-pc_mode$logm_second, 0)
  weighted_X <- sweep(X_star, 1L, as.numeric(weights_at_mode), "*")
  proposal_precision <- crossprod(X_star, weighted_X) / sigma2 +
    diag(prior_precision + 1e-8, q)
  eig <- eigen(proposal_precision, symmetric = TRUE)
  proposal_cov <- eig$vectors %*%
    diag(1 / pmax(eig$values, 1e-8), q) %*%
    t(eig$vectors)
  proposal_cov <- 0.5 * (proposal_cov + t(proposal_cov))
  coordinate_sd <- sqrt(pmax(diag(proposal_cov), 1e-12))

  # Geometrically spaced inverse temperatures: 1 is the target posterior and
  # smaller values flatten the posterior to facilitate movement between modes.
  inverse_temperatures <- exp(seq(
    from = 0,
    to = log(min_inverse_temp),
    length.out = n_temperatures
  ))

  n_keep <- floor((n_iter - burn) / thin)
  if (n_keep < 20L) stop("The MCMC settings must retain at least 20 draws per chain.")
  native_sampler <- .get_mom_tempering_cpp()(
    beta_init, y_star, X_star, inverse_temperatures, coordinate_sd,
    nu, sigma2, tau0, tau1, pi, n_iter, burn, thin,
    n_chains, n_temperatures, large_jump_prob, swap_every
  )
  pip_draws <- native_sampler$pip_draws
  beta_means <- native_sampler$beta_means
  local_acceptance <- native_sampler$local_acceptance
  large_jump_acceptance <- native_sampler$large_jump_acceptance
  swap_acceptance <- native_sampler$swap_acceptance
  proposal_scales <- native_sampler$proposal_scales

  pip_array <- do.call(rbind, pip_draws)
  chain_pip <- t(vapply(pip_draws, colMeans, numeric(q)))
  integrated_pip <- colMeans(pip_array)
  max_chain_range <- max(apply(chain_pip, 2L, function(z) diff(range(z))))

  # Batch-means Monte Carlo standard errors for the Rao-Blackwellized PIPs.
  batch_means <- lapply(pip_draws, function(z) {
    batch_size <- max(10L, floor(sqrt(nrow(z))))
    n_batch <- floor(nrow(z) / batch_size)
    z <- z[seq_len(n_batch * batch_size), , drop = FALSE]
    vapply(seq_len(n_batch), function(b) {
      idx <- ((b - 1L) * batch_size + 1L):(b * batch_size)
      colMeans(z[idx, , drop = FALSE])
    }, numeric(q))
  })
  batch_means <- t(do.call(cbind, batch_means))
  pip_mcse <- apply(batch_means, 2L, sd) / sqrt(nrow(batch_means))

  if (verbose) {
    cat("Integrated-PIP cold-chain acceptance rates:",
        paste(sprintf("%.3f", local_acceptance[, 1L]), collapse = ", "), "\n")
    cat("Cold-chain large-jump acceptance rates:",
        paste(sprintf("%.3f", large_jump_acceptance[, 1L]), collapse = ", "), "\n")
    cat("Smallest adjacent-temperature swap rate:",
        sprintf("%.3f", min(swap_acceptance, na.rm = TRUE)), "\n")
    cat("Largest between-chain PIP range:",
        sprintf("%.3f", max_chain_range), "\n")
  }

  list(
    pip = integrated_pip,
    pip_mcse = pip_mcse,
    chain_pip = chain_pip,
    beta_posterior_mean = colMeans(beta_means),
    acceptance = local_acceptance[, 1L],
    tempered_acceptance = local_acceptance,
    large_jump_acceptance = large_jump_acceptance,
    swap_acceptance = swap_acceptance,
    inverse_temperatures = inverse_temperatures,
    proposal_scales = proposal_scales,
    max_chain_pip_range = max_chain_range,
    retained_draws_per_chain = n_keep
  )
}


# 1b. Legacy numerical hyperparameter update retained for compatibility

# The generalized-EM fitter below uses responsibility-weighted closed-form
# conditional updates instead. This helper is retained for external callers.
# Update spike (tau0) and slab (tau1) scales given current beta and hyperparameters.
# We work on the log(tau^2) scale and use 1D optimize() for each variance.
update_tau_scales <- function(beta, pi, tau0, tau1, a0, b0, a1, b1,
                              log_tau2_lower = -20, log_tau2_upper = 5) {
  tiny <- 1e-16

  # Update tau0^2
  f_tau0 <- function(log_tau2) {
    tau2 <- exp(log_tau2)
    # spike density depends on tau0^2 (normalized)
    f0 <- (2 * base::pi * tau2)^(-0.5) * exp(-0.5 * (beta^2 / tau2))
    # slab uses current tau1^2 (normalized)
    f1 <- (2 * base::pi)^(-0.5) * (tau1^2)^(-1.5) * (beta^2) * exp(-0.5 * (beta^2 / (tau1^2)))
    m  <- (1 - pi) * f0 + pi * f1
    loglik <- sum(log(m + tiny))

    # Inv-Gamma(a0, b0) prior on tau0^2
    log_prior <- a0 * log(b0) - lgamma(a0) - (a0 + 1) * log(tau2) - b0 / tau2

    -(loglik + log_prior)
  }
  opt0 <- optimize(f_tau0, lower = log_tau2_lower, upper = log_tau2_upper)
  tau0_new <- sqrt(exp(opt0$minimum))

  # Update tau1^2
  f_tau1 <- function(log_tau2) {
    tau2 <- exp(log_tau2)
    # spike uses the newly updated tau0^2 (normalized)
    f0 <- (2 * base::pi * (tau0_new^2))^(-0.5) * exp(-0.5 * (beta^2 / (tau0_new^2)))
    # slab density depends on tau1^2 (normalized)
    f1 <- (2 * base::pi)^(-0.5) * (tau2)^(-1.5) * (beta^2) * exp(-0.5 * (beta^2 / tau2))
    m  <- (1 - pi) * f0 + pi * f1
    loglik <- sum(log(m + tiny))

    # Inv-Gamma(a1, b1) prior on tau1^2
    log_prior <- a1 * log(b1) - lgamma(a1) - (a1 + 1) * log(tau2) - b1 / tau2

    -(loglik + log_prior)
  }
  opt1 <- optimize(f_tau1, lower = log_tau2_lower, upper = log_tau2_upper)
  tau1_new <- sqrt(exp(opt1$minimum))

  list(tau0 = tau0_new, tau1 = tau1_new)
}


# 2. Utilities for covariance, whitening, and covariance-component estimation

# Build Sigma_b(alpha) = sum_k alpha_k V_k  (V_k: p x p)
build_Sigma_b <- function(alpha, V_list, p) {
  K <- length(V_list)
  if (is.null(alpha) || K == 0) {
    return(matrix(0, nrow = p, ncol = p))
  }
  if (length(alpha) != K) {
    stop("Length of alpha must match length of V_list.")
  }
  Sigma_b <- matrix(0, nrow = p, ncol = p)
  for (k in seq_len(K)) {
    if (alpha[k] != 0) {
      Sigma_b <- Sigma_b + alpha[k] * V_list[[k]]
    }
  }
  Sigma_b
}

# Pre-whiten data for a given alpha. The whitening method is part of the
# coordinate-wise Student-t model because independent t coordinates are not
# invariant to orthogonal rotations.
# For each subject n,
# M_n = C_theta_y_n + Sigma_b
# S_n = solve(M_n)
# B_n^T B_n = S_n = (C_n + Sigma_b)^(-1)
# Then y_star = B(α) vec(eta_theta_y), X_star = B(α) (X_G ⊗ I_p)
whiten_data <- function(eta_theta_y, C_theta_y_list, X_G, Sigma_b,
                        whitening_method = c("precision_chol", "covariance_chol", "symmetric")) {
  whitening_method <- match.arg(whitening_method)
  eta_theta_y <- as.matrix(eta_theta_y)
  X_G   <- as.matrix(X_G)
  N     <- nrow(eta_theta_y)
  p     <- ncol(eta_theta_y)
  r     <- ncol(X_G)

  if (length(C_theta_y_list) != N) {
    stop("C_theta_y_list must be a list of length N with p x p matrices.")
  }

  # Total dimension N*p
  Np <- N * p
  q  <- p * r

  y_star <- numeric(Np)
  X_star <- matrix(0, nrow = Np, ncol = q)

  S_list <- vector("list", N)  # store S_n = (C_n + Sigma_b)^(-1)
  B_list <- vector("list", N)  # store the specified whitening factors
  log_det_B <- 0

  row_start <- 1L
  for (n in seq_len(N)) {
    Cn <- C_theta_y_list[[n]]
    if (!all(dim(Cn) == c(p, p))) {
      stop("Each element of C_theta_y_list must be p x p.")
    }

    M_n <- Cn + Sigma_b
    # Small jitter for numerical stability
    M_n <- M_n + 1e-8 * diag(p)

    S_n <- solve(M_n)
    if (whitening_method == "precision_chol") {
      B_n <- chol(S_n)
    } else if (whitening_method == "covariance_chol") {
      U_n <- chol(M_n)
      B_n <- solve(t(U_n))
    } else {
      ee <- eigen(M_n, symmetric = TRUE)
      if (min(ee$values) <= 0) stop("C_n + Sigma_b must be positive definite.")
      B_n <- ee$vectors %*% diag(1 / sqrt(ee$values), p) %*% t(ee$vectors)
    }

    S_list[[n]] <- S_n
    B_list[[n]] <- B_n
    log_det_B <- log_det_B + as.numeric(determinant(B_n, logarithm = TRUE)$modulus)

    # Indices for this subject's block in y_star and X_star
    idx <- row_start:(row_start + p - 1L)

    # Subject-specific design: X_tilde_n = (x_n^T ⊗ I_p), dimension p x (p*r)
    Xn_tilde <- kronecker(matrix(X_G[n, ], nrow = 1), diag(p))

    # Pre-whiten
    y_star[idx]   <- as.numeric(B_n %*% eta_theta_y[n, ])
    X_star[idx, ] <- B_n %*% Xn_tilde

    row_start <- row_start + p
  }

  list(
    y_star = y_star,  # Np-vector
    X_star = X_star,  # Np x (p*r)
    S_list = S_list,  # list of S_n = (C_n + Sigma_b)^(-1)
    B_list = B_list,
    log_det_B = log_det_B,
    whitening_method = whitening_method
  )
}

# Normalized log posterior used for every empirical-Bayes update. The Jacobian
# log|B(alpha)| is essential when alpha changes because whitening is then part
# of the likelihood rather than a fixed preprocessing operation.
normalized_t_log_posterior <- function(
    beta, sigma2, tau0, tau1, alpha, eta_theta_y, C_theta_y_list,
    X_G, V_list, nu, pi, a0, b0, a1, b1,
    whitening_method = "precision_chol"
) {
  p <- ncol(as.matrix(eta_theta_y))
  Sigma_b <- build_Sigma_b(alpha, V_list, p)
  wd <- whiten_data(
    eta_theta_y, C_theta_y_list, X_G, Sigma_b,
    whitening_method = whitening_method
  )
  residual <- as.numeric(wd$y_star - wd$X_star %*% beta)
  log_likelihood <- wd$log_det_B + sum(
    dt(residual / sqrt(sigma2), df = nu, log = TRUE) - 0.5 * log(sigma2)
  )
  log_prior_beta <- sum(log_mixture_MOM(beta, tau0, tau1, pi)$log_m)
  log_prior_tau <- log_tau_hyperprior(tau0, tau1, a0, b0, a1, b1)
  log_prior_sigma <- -log(sigma2)
  log_likelihood + log_prior_beta + log_prior_tau + log_prior_sigma
}

# Expected complete-data log posterior for the normal-gamma representation of
# the Student-t likelihood. The supplied weights are E(lambda_i | current fit)
# from the E-step and remain fixed throughout the following conditional
# maximization steps. Terms independent of the fitted parameters are omitted.
expected_complete_t_log_posterior <- function(
    beta, sigma2, tau0, tau1, alpha, weights, responsibility, eta_theta_y,
    C_theta_y_list, X_G, V_list, pi, a0, b0, a1, b1,
    whitening_method = "precision_chol"
) {
  p <- ncol(as.matrix(eta_theta_y))
  Sigma_b <- build_Sigma_b(alpha, V_list, p)
  wd <- whiten_data(
    eta_theta_y, C_theta_y_list, X_G, Sigma_b,
    whitening_method = whitening_method
  )
  residual <- as.numeric(wd$y_star - wd$X_star %*% beta)
  if (length(weights) != length(residual)) {
    stop("The E-step weights must match the number of whitened observations.")
  }
  q_likelihood <- wd$log_det_B -
    (length(residual) / 2 + 1) * log(sigma2) -
    sum(weights * residual^2) / (2 * sigma2)
  q_likelihood +
    expected_mixture_log_prior_MOM(
      beta, responsibility, tau0, tau1, pi
    ) +
    log_tau_hyperprior(tau0, tau1, a0, b0, a1, b1)
}

# Legacy Gaussian covariance score retained only for backward compatibility.
# The main fitter no longer uses this function because its Q/Tr update does not
# optimize the normalized Student-t objective.
# Compute the Gaussian gradient and Hessian for alpha using block structure:
# Pi = blockdiag(S_1, ..., S_N), with S_n = (C_n + Sigma_b)^(-1).
# dCov/d alpha_k = I_N ⊗ V_k -> blockdiag(V_k,...,V_k).
# g_k = 0.5 * tr(Pi dCov/d alpha_k) - 0.5 * e^T Pi dCov/d alpha_k Pi e
# H_kl = 0.5 * tr(Pi dCov/d alpha_k Pi dCov/d alpha_l)
# where e = vec(eta_theta_y) - (X_G ⊗ I_p) beta.
# Also returns the separated trace and quadratic terms
# Tr_k = sum_n tr(S_n V_k),   Qd_k = sum_n e_n^T S_n V_k S_n e_n,
# These quantities drove the superseded multiplicative Q/Tr update.
alpha_grad_hess <- function(S_list, V_list, X_G, beta_vec, eta_theta_y) {
  eta_theta_y <- as.matrix(eta_theta_y)
  X_G   <- as.matrix(X_G)

  N <- nrow(eta_theta_y)
  p <- ncol(eta_theta_y)
  r <- ncol(X_G)
  q <- length(beta_vec)
  if (q != p * r) {
    stop("Length of beta_vec must be p * r.")
  }

  K <- length(V_list)
  eta_vec <- as.vector(t(eta_theta_y))  # vec(eta_theta_y) in subject-major order
  X_big   <- kronecker(X_G, diag(p))  # Np x (p*r)
  mu_vec  <- as.numeric(X_big %*% beta_vec)
  e_vec   <- eta_vec - mu_vec

  g <- numeric(K)
  Tr <- numeric(K)  # trace terms      sum_n tr(S_n V_k)
  Qd <- numeric(K)  # quadratic terms  sum_n e_n^T S_n V_k S_n e_n
  H <- matrix(0, nrow = K, ncol = K)

  Np_offset <- 0L
  for (n in seq_len(N)) {
    S_n <- S_list[[n]]
    idx <- (Np_offset + 1L):(Np_offset + p)
    e_n <- e_vec[idx]

    for (k in seq_len(K)) {
      V_k <- V_list[[k]]

      # trace term for gradient: tr(S_n V_k)
      tr_k <- sum(diag(S_n %*% V_k))
      Tr[k] <- Tr[k] + tr_k
      g[k]  <- g[k] + 0.5 * tr_k
    }

    # second part of gradient and Hessian
    for (k in seq_len(K)) {
      V_k <- V_list[[k]]
      SK  <- S_n %*% V_k %*% S_n  # S_n V_k S_n

      # e^T S_n V_k S_n e
      q_k   <- as.numeric(t(e_n) %*% SK %*% e_n)
      Qd[k] <- Qd[k] + q_k
      g[k]  <- g[k] - 0.5 * q_k

      for (l in k:K) {
        V_l <- V_list[[l]]
        # tr(S_n V_k S_n V_l S_n)
        H_kl_add <- 0.5 * sum(diag(SK %*% V_l %*% S_n))
        H[k, l] <- H[k, l] + H_kl_add
        if (l != k) {
          H[l, k] <- H[l, k] + H_kl_add
        }
      }
    }

    Np_offset <- Np_offset + p
  }

  list(g = g, H = H, Tr = Tr, Qd = Qd)
}


# 3. Main fitter: robust & sparse group DCM
# The fitting stage is a generalized EM/ECM algorithm. Each E-step computes
# latent Student-t precision expectations. Conditional maximization steps for
# beta, sigma2, the prior scales, and alpha then increase the same expected
# complete-data log posterior. Integrated PIPs are computed after convergence.

fit_group_dcm_t_MOM <- function(
    eta_theta_y,  # N x p matrix: subject-level means eta_{theta|y}
    C_theta_y_list,  # list of length N: p x p matrices C_{theta|y,n}
    X_G,  # N x r group design matrix
    V_list,  # list of p x p basis matrices V_k
    nu            = 3,  # Student-t df (recommended default; simulation study uses nu=2)
    tau0          = 0.05,  # initial spike SD (small)
    tau1          = 1.0,  # initial slab SD (large)
    pi            = 0.4,  # default prior inclusion probability; analyses set this explicitly
    a0            = 2,  # Inv-Gamma(a0, b0) prior on tau0^2 (spike: small scale)
    b0            = 1e-2,  # prior mode of tau0 ~ sqrt(b0/(a0+1)) ~ 0.058
    a1            = 2,  # Inv-Gamma(a1, b1) prior on tau1^2 (slab: large scale)
    b1            = 1.0,  # prior mode of tau1 ~ sqrt(b1/(a1+1)) ~ 0.577
    min_slab_spike_ratio = 1,  # identifiability guard: enforce tau1 >= ratio * tau0
    max_iter      = 200,
    tol           = 1e-6,
    inner_sweeps  = 3,
    verbose       = TRUE,
    max_beta_step  = 1.0,  # step size cap for each beta coordinate Newton step
    min_alpha      = 1e-8,  # lower bound for variance components
    max_alpha      = 100,  # upper bound to stabilize alpha
    whitening_method = c("precision_chol", "covariance_chol", "symmetric"),
    estimate_sigma2 = FALSE,  # FALSE fixes sigma2 and identifies the covariance scale
    sigma2_fixed   = 1,
    beta_init      = NULL,
    alpha_init     = NULL,
    n_starts       = 3L,
    start_jitter   = 0.25,
    start_seed     = NULL,
    objective_tol  = 1e-10,
    gradient_tol   = 1e-4,
    max_backtrack  = 20L,
    pip_all_starts = FALSE,
    pip_mcmc_iter  = 10000L,  # posterior-sampling iterations per chain
    pip_mcmc_burn  = 3000L,  # burn-in iterations per chain
    pip_mcmc_thin  = 2L,  # retain every pip_mcmc_thin-th post-burn draw
    pip_mcmc_chains = 4L,  # multiple dispersed chains for separated modes
    pip_mcmc_temperatures = 8L,  # replicas in each parallel-tempering ladder
    pip_mcmc_min_inverse_temp = 0.05,  # inverse temperature of hottest replica
    pip_mcmc_swap_every = 2L,  # iterations between adjacent-temperature swaps
    pip_mcmc_large_jump_prob = 0.10,  # coordinate proposals using slab-scale jumps
    pip_mcmc_seed  = NULL  # optional reproducibility seed for PIP sampling
) {
  whitening_method <- match.arg(whitening_method)
  eta_theta_y <- as.matrix(eta_theta_y)
  X_G   <- as.matrix(X_G)
  N     <- nrow(eta_theta_y)
  p     <- ncol(eta_theta_y)
  r     <- ncol(X_G)
  q     <- p * r
  Np    <- N * p

  if (nrow(X_G) != N) {
    stop("Number of rows of X_G must match number of rows (subjects) of eta_theta_y.")
  }
  if (length(C_theta_y_list) != N) {
    stop("C_theta_y_list must be a list of length N.")
  }
  if (!is.finite(nu) || nu <= 0 || !is.finite(pi) || pi <= 0 || pi >= 1) {
    stop("nu must be positive and pi must lie in (0, 1).")
  }
  if (!is.finite(sigma2_fixed) || sigma2_fixed <= 0) {
    stop("sigma2_fixed must be positive.")
  }
  n_starts <- as.integer(n_starts)
  if (n_starts < 1L || start_jitter < 0 || max_backtrack < 1L || gradient_tol <= 0) {
    stop("Use n_starts >= 1, start_jitter >= 0, max_backtrack >= 1, and gradient_tol > 0.")
  }
  if (estimate_sigma2) {
    warning(
      "A free sigma2 can be confounded with covariance-component scale when the supplied ",
      "C matrices and variance basis are proportional. Use only after verifying identification."
    )
  }

  K <- length(V_list)
  if (K > 0L) {
    if (any(vapply(V_list, function(V) !all(dim(V) == c(p, p)), logical(1L)))) {
      stop("Every covariance basis matrix in V_list must be p x p.")
    }
    if (any(vapply(V_list, function(V) max(abs(V - t(V))) > 1e-8, logical(1L)))) {
      stop("Every covariance basis matrix in V_list must be symmetric.")
    }
    V_design <- do.call(cbind, lapply(V_list, as.vector))
    if (qr(V_design)$rank < K) {
      stop("V_list is linearly dependent, so the alpha components are not identifiable.")
    }
  }
  X_big   <- kronecker(X_G, diag(p))
  y_big   <- as.vector(t(eta_theta_y))
  XtX     <- crossprod(X_big)
  Xty     <- crossprod(X_big, y_big)
  beta_base <- if (is.null(beta_init)) {
    as.numeric(solve(XtX + 1e-6 * diag(q), Xty))
  } else {
    if (length(beta_init) != q) stop("beta_init must have length p * r.")
    as.numeric(beta_init)
  }
  alpha_base <- if (K == 0L) {
    numeric(0)
  } else if (is.null(alpha_init)) {
    rep(0.1, K)
  } else {
    if (length(alpha_init) != K) stop("alpha_init must match V_list.")
    pmin(pmax(as.numeric(alpha_init), min_alpha), max_alpha)
  }
  sigma2_base <- if (estimate_sigma2) {
    max(var(y_big - as.numeric(X_big %*% beta_base)), 1e-8)
  } else {
    sigma2_fixed
  }

  q_objective_whitened <- function(beta_value, residual, wd, sigma2_value,
                                   tau0_value, tau1_value, weights_value,
                                   responsibility_value) {
    wd$log_det_B - (length(residual) / 2 + 1) * log(sigma2_value) -
      sum(weights_value * residual^2) / (2 * sigma2_value) +
      expected_mixture_log_prior_MOM(
        beta_value, responsibility_value, tau0_value, tau1_value, pi
      ) +
      log_tau_hyperprior(tau0_value, tau1_value, a0, b0, a1, b1)
  }

  fit_one_start <- function(beta_start, alpha_start, sigma2_start,
                            tau0_start, tau1_start, start_id) {
    beta <- beta_start
    alpha <- alpha_start
    sigma2 <- if (estimate_sigma2) sigma2_start else sigma2_fixed
    tau0_current <- tau0_start
    tau1_current <- max(tau1_start, min_slab_spike_ratio * tau0_current)
    objective_trace <- rep(NA_real_, max_iter + 1L)
    objective_trace[1L] <- normalized_t_log_posterior(
      beta, sigma2, tau0_current, tau1_current, alpha,
      eta_theta_y, C_theta_y_list, X_G, V_list, nu, pi,
      a0, b0, a1, b1, whitening_method
    )
    q_start_trace <- rep(NA_real_, max_iter)
    q_end_trace <- rep(NA_real_, max_iter)
    converged <- FALSE
    stopping_reason <- "maximum iterations reached"
    alpha_optim_convergence <- NA_integer_
    final_parameter_change <- Inf
    final_objective_change <- Inf
    current_beta_score_max <- Inf

    if (verbose) {
      cat(sprintf("Starting generalized EM from initialization %d of %d\n",
                  start_id, n_starts))
    }

    for (iter in seq_len(max_iter)) {
      beta_old <- beta
      alpha_old <- alpha
      sigma2_old <- sigma2
      tau0_old <- tau0_current
      tau1_old <- tau1_current
      objective_old <- objective_trace[iter]

      Sigma_b <- build_Sigma_b(alpha, V_list, p)
      wd <- whiten_data(
        eta_theta_y, C_theta_y_list, X_G, Sigma_b,
        whitening_method = whitening_method
      )
      residual <- as.numeric(wd$y_star - wd$X_star %*% beta)

      # E-step: conditional expectations of the latent normal-gamma precisions
      # and spike-and-slab indicators. These quantities remain fixed throughout
      # this generalized-EM iteration. The indicator responsibilities are
      # optimization quantities and are not reported as integrated PIPs.
      weights_e <- (nu + 1) /
        (nu + residual^2 / (sigma2 + 1e-12))
      responsibility_e <- log_mixture_MOM(
        beta, tau0_current, tau1_current, pi
      )$responsibility
      current_q <- q_objective_whitened(
        beta, residual, wd, sigma2, tau0_current, tau1_current,
        weights_e, responsibility_e
      )
      q_start_trace[iter] <- current_q

      # Conditional M-step for beta. The likelihood gradient is recomputed after
      # each accepted coordinate update using the fixed E-step weights.
      # Backtracking ensures that every update increases the common Q function.
      for (sweep in seq_len(inner_sweeps)) {
        for (j in seq_len(q)) {
          xj <- wd$X_star[, j]
          g_lik <- sum(xj * weights_e * residual) / sigma2
          h_lik <- -sum(weights_e * xj^2) / sigma2
          beta_safe <- if (abs(beta[j]) < 1e-8) {
            if (beta[j] >= 0) 1e-8 else -1e-8
          } else beta[j]
          g_prior <- -(1 - responsibility_e[j]) * beta[j] / tau0_current^2 +
            responsibility_e[j] *
            (2 / beta_safe - beta_safe / tau1_current^2)
          h_prior <- -(1 - responsibility_e[j]) / tau0_current^2 +
            responsibility_e[j] *
            (-2 / beta_safe^2 - 1 / tau1_current^2)
          g_j <- g_lik + g_prior
          h_j <- h_lik + h_prior
          if (!is.finite(h_j) || h_j > -1e-8) h_j <- h_lik - 1e-8
          step_j <- -g_j / h_j
          if (!is.finite(step_j)) next
          step_j <- sign(step_j) * min(abs(step_j), max_beta_step)

          for (backtrack in 0:max_backtrack) {
            step_try <- step_j * 0.5^backtrack
            beta_try <- beta
            beta_try[j] <- beta[j] + step_try
            residual_try <- residual - xj * step_try
            q_try <- q_objective_whitened(
              beta_try, residual_try, wd, sigma2,
              tau0_current, tau1_current, weights_e, responsibility_e
            )
            if (is.finite(q_try) && q_try >= current_q - objective_tol) {
              beta <- beta_try
              residual <- residual_try
              current_q <- q_try
              break
            }
          }
        }
      }

      # Conditional M-step for the optional global scale. The closed-form MAP
      # update follows from the fixed-weight Q function and a Jeffreys prior.
      if (estimate_sigma2) {
        sigma2_try <- max(
          sum(weights_e * residual^2) / (length(residual) + 2),
          exp(-20)
        )
        q_try <- q_objective_whitened(
          beta, residual, wd, sigma2_try,
          tau0_current, tau1_current, weights_e, responsibility_e
        )
        if (is.finite(q_try) && q_try >= current_q - objective_tol) {
          sigma2 <- sigma2_try
          current_q <- q_try
        }
      } else {
        sigma2 <- sigma2_fixed
      }

      # Conditional M-step for the spike and slab scales.
      tau0_sq_try <- (b0 + 0.5 * sum((1 - responsibility_e) * beta^2)) /
        (a0 + 1 + 0.5 * sum(1 - responsibility_e))
      tau1_sq_try <- (b1 + 0.5 * sum(responsibility_e * beta^2)) /
        (a1 + 1 + 1.5 * sum(responsibility_e))
      tau0_try <- sqrt(max(tau0_sq_try, exp(-20)))
      tau1_try <- max(
        sqrt(max(tau1_sq_try, exp(-20))),
        min_slab_spike_ratio * tau0_try
      )
      q_try <- q_objective_whitened(
        beta, residual, wd, sigma2, tau0_try, tau1_try,
        weights_e, responsibility_e
      )
      if (is.finite(q_try) && q_try >= current_q - objective_tol) {
        tau0_current <- tau0_try
        tau1_current <- tau1_try
        current_q <- q_try
      }

      # Conditional M-step for covariance components. The fixed-weight Q
      # function includes log|B(alpha)| and the alpha-dependent residuals.
      if (K > 0L) {
        alpha_objective <- function(log_alpha) {
          -expected_complete_t_log_posterior(
            beta, sigma2, tau0_current, tau1_current, exp(log_alpha),
            weights_e, responsibility_e, eta_theta_y, C_theta_y_list,
            X_G, V_list, pi,
            a0, b0, a1, b1, whitening_method
          )
        }

        alpha_opt <- optim(
          par = log(pmax(alpha, min_alpha)),
          fn = alpha_objective,
          method = "L-BFGS-B",
          lower = rep(log(min_alpha), K),
          upper = rep(log(max_alpha), K),
          control = list(maxit = 200L) #50L
        )

        if (alpha_opt$convergence != 0L || !is.finite(alpha_opt$value)) {
          alpha_retry <- nlminb(
            start = log(pmax(alpha, min_alpha)),
            objective = alpha_objective,
            lower = rep(log(min_alpha), K),
            upper = rep(log(max_alpha), K),
            control = list(
              iter.max = 500L,
              eval.max = 1000L,
              rel.tol = 1e-8
            )
          )
          alpha_opt <- list(
            par = alpha_retry$par,
            value = alpha_retry$objective,
            convergence = alpha_retry$convergence
          )
        }

        alpha_optim_convergence <- alpha_opt$convergence
        alpha_try <- exp(alpha_opt$par)
        q_try <- -alpha_opt$value
        if (is.finite(q_try) && q_try >= current_q - objective_tol) {
          alpha <- alpha_try
          current_q <- q_try
        }
      }

      q_end_trace[iter] <- current_q
      observed_objective <- normalized_t_log_posterior(
        beta, sigma2, tau0_current, tau1_current, alpha,
        eta_theta_y, C_theta_y_list, X_G, V_list, nu, pi,
        a0, b0, a1, b1, whitening_method
      )

      # A valid generalized-EM step cannot reduce the observed log posterior.
      # Roll back and stop if numerical optimization violates that safeguard.
      if (!is.finite(observed_objective) ||
          observed_objective < objective_old - objective_tol) {
        beta <- beta_old
        alpha <- alpha_old
        sigma2 <- sigma2_old
        tau0_current <- tau0_old
        tau1_current <- tau1_old
        objective_trace[iter + 1L] <- objective_old
        stopping_reason <- "observed objective decreased; final GEM iteration rolled back"
        break
      }

      objective_trace[iter + 1L] <- observed_objective
      diff_beta <- max(abs(beta - beta_old) / (1 + abs(beta_old)))
      diff_alpha <- if (K > 0L) {
        max(abs(log(pmax(alpha, min_alpha) / pmax(alpha_old, min_alpha))))
      } else 0
      diff_sigma <- abs(log(sigma2 / sigma2_old))
      diff_tau <- max(abs(log(tau0_current / tau0_old)),
                      abs(log(tau1_current / tau1_old)))
      final_parameter_change <- max(diff_beta, diff_alpha, diff_sigma, diff_tau)
      final_objective_change <- abs(observed_objective - objective_old) /
        (1 + abs(objective_old))

      # Evaluate the coefficient score after all updates, including the new
      # covariance components, so convergence also requires local stationarity.
      Sigma_b_score <- build_Sigma_b(alpha, V_list, p)
      wd_score <- whiten_data(
        eta_theta_y, C_theta_y_list, X_G, Sigma_b_score,
        whitening_method = whitening_method
      )
      residual_score <- as.numeric(wd_score$y_star - wd_score$X_star %*% beta)
      weights_score <- (nu + 1) /
        (nu + residual_score^2 / (sigma2 + 1e-12))
      prior_score <- prior_components_MOM(
        beta, tau0_current, tau1_current, pi
      )$logm_prime
      beta_score_current <- as.numeric(crossprod(
        wd_score$X_star, weights_score * residual_score
      )) / sigma2 + prior_score
      current_beta_score_max <- max(abs(beta_score_current))

      if (verbose && (iter %% 20L == 0L || iter == 1L)) {
        cat(sprintf(
          "Start %d, iter %3d: max relative parameter change = %.3e, relative objective change = %.3e, max|score| = %.3e, sigma2 = %.3f, tau0 = %.4f, tau1 = %.4f\n",
          start_id, iter, final_parameter_change, final_objective_change,
          current_beta_score_max, sigma2, tau0_current, tau1_current
        ))
      }

      alpha_stationary <- K == 0L || isTRUE(alpha_optim_convergence == 0L)
      if (final_parameter_change < tol && final_objective_change < tol &&
          current_beta_score_max < gradient_tol && alpha_stationary) {
        converged <- TRUE
        stopping_reason <- "parameter, objective, and stationarity tolerances satisfied"
        break
      }
    }

    iterations <- iter
    objective_trace <- objective_trace[seq_len(iterations + 1L)]
    q_start_trace <- q_start_trace[seq_len(iterations)]
    q_end_trace <- q_end_trace[seq_len(iterations)]
    Sigma_b <- build_Sigma_b(alpha, V_list, p)
    wd <- whiten_data(
      eta_theta_y, C_theta_y_list, X_G, Sigma_b,
      whitening_method = whitening_method
    )
    residual <- as.numeric(wd$y_star - wd$X_star %*% beta)
    weights_final <- (nu + 1) /
      (nu + residual^2 / (sigma2 + 1e-12))
    pc_final <- prior_components_MOM(beta, tau0_current, tau1_current, pi)
    mode_responsibility <- log_mixture_MOM(
      beta, tau0_current, tau1_current, pi
    )$responsibility
    beta_score <- as.numeric(crossprod(
      wd$X_star, weights_final * residual
    )) / sigma2 + pc_final$logm_prime

    list(
      beta = beta,
      alpha = alpha,
      sigma2 = sigma2,
      tau0 = tau0_current,
      tau1 = tau1_current,
      wd = wd,
      residual = residual,
      weights = weights_final,
      mode_responsibility = mode_responsibility,
      objective = tail(objective_trace, 1L),
      objective_trace = objective_trace,
      objective_monotone = all(diff(objective_trace) >= -sqrt(.Machine$double.eps)),
      q_start_trace = q_start_trace,
      q_end_trace = q_end_trace,
      q_monotone = all(q_end_trace >= q_start_trace - sqrt(.Machine$double.eps)),
      converged = converged,
      stopping_reason = stopping_reason,
      iterations = iterations,
      final_parameter_change = final_parameter_change,
      final_objective_change = final_objective_change,
      beta_score_max = max(abs(beta_score)),
      alpha_optim_convergence = alpha_optim_convergence
    )
  }

  if (!is.null(start_seed)) set.seed(start_seed)
  start_fits <- vector("list", n_starts)
  beta_base_sd <- if (q > 1L) sd(beta_base) else 0
  beta_jitter_sd <- start_jitter * max(beta_base_sd, tau0, 0.1)
  for (start in seq_len(n_starts)) {
    if (start == 1L) {
      beta_start <- beta_base
      alpha_start <- alpha_base
      tau0_start <- tau0
      tau1_start <- tau1
      sigma2_start <- sigma2_base
    } else {
      beta_start <- beta_base + rnorm(q, sd = beta_jitter_sd)
      alpha_start <- if (K > 0L) {
        pmin(pmax(alpha_base * exp(rnorm(K, sd = start_jitter)), min_alpha), max_alpha)
      } else numeric(0)
      tau0_start <- tau0 * exp(rnorm(1L, sd = start_jitter))
      tau1_start <- max(
        tau1 * exp(rnorm(1L, sd = start_jitter)),
        min_slab_spike_ratio * tau0_start
      )
      sigma2_start <- if (estimate_sigma2) {
        sigma2_base * exp(rnorm(1L, sd = start_jitter))
      } else sigma2_fixed
    }
    start_fits[[start]] <- fit_one_start(
      beta_start, alpha_start, sigma2_start,
      tau0_start, tau1_start, start
    )
  }

  start_objectives <- vapply(start_fits, `[[`, numeric(1L), "objective")
  best_start <- which.max(start_objectives)
  best <- start_fits[[best_start]]
  beta <- best$beta
  alpha <- best$alpha
  sigma2 <- best$sigma2
  tau0 <- best$tau0
  tau1 <- best$tau1
  w <- best$weights
  wd_final <- best$wd
  y_star <- wd_final$y_star
  X_star <- wd_final$X_star
  r_star <- best$residual

  start_summary <- data.frame(
    start = seq_len(n_starts),
    objective = start_objectives,
    converged = vapply(start_fits, `[[`, logical(1L), "converged"),
    iterations = vapply(start_fits, `[[`, integer(1L), "iterations"),
    objective_monotone = vapply(start_fits, `[[`, logical(1L), "objective_monotone"),
    q_monotone = vapply(start_fits, `[[`, logical(1L), "q_monotone"),
    beta_score_max = vapply(start_fits, `[[`, numeric(1L), "beta_score_max")
  )
  if (verbose) {
    cat(sprintf("Selected initialization %d with normalized log posterior %.6f\n",
                best_start, best$objective))
    if (!best$converged) warning("The selected initialization did not converge.")
  }

  multi_start_beta <- array(
    NA_real_, dim = c(p, r, n_starts),
    dimnames = list(colnames(eta_theta_y), colnames(X_G), paste0("start", seq_len(n_starts)))
  )
  multi_start_alpha <- matrix(
    NA_real_, nrow = n_starts, ncol = K,
    dimnames = list(paste0("start", seq_len(n_starts)), if (K > 0L) paste0("alpha", seq_len(K)) else NULL)
  )
  multi_start_tau <- matrix(
    NA_real_, nrow = n_starts, ncol = 3L,
    dimnames = list(paste0("start", seq_len(n_starts)), c("sigma2", "tau0", "tau1"))
  )
  for (start in seq_len(n_starts)) {
    multi_start_beta[, , start] <- matrix(
      start_fits[[start]]$beta, nrow = p, ncol = r, byrow = FALSE
    )
    if (K > 0L) multi_start_alpha[start, ] <- start_fits[[start]]$alpha
    multi_start_tau[start, ] <- c(
      start_fits[[start]]$sigma2,
      start_fits[[start]]$tau0,
      start_fits[[start]]$tau1
    )
  }

  # Integrated posterior inclusion probabilities. Conditional on the fitted
  # empirical-Bayes hyperparameters, sample the joint posterior of beta and
  # average Pr(gamma_j = 1 | beta_j, h) over the retained posterior draws.
  pip_fit <- sample_integrated_pip_MOM(
    beta_init = beta,
    y_star = y_star,
    X_star = X_star,
    nu = nu,
    sigma2 = sigma2,
    tau0 = tau0,
    tau1 = tau1,
    pi = pi,
    weights_at_mode = w,
    n_iter = pip_mcmc_iter,
    burn = pip_mcmc_burn,
    thin = pip_mcmc_thin,
    n_chains = pip_mcmc_chains,
    n_temperatures = pip_mcmc_temperatures,
    min_inverse_temp = pip_mcmc_min_inverse_temp,
    swap_every = pip_mcmc_swap_every,
    large_jump_prob = pip_mcmc_large_jump_prob,
    seed = pip_mcmc_seed,
    verbose = verbose
  )
  inclusion_prob_vec <- pip_fit$pip

  # Optional full PIP sensitivity analysis across empirical-Bayes starts. It is
  # disabled by default because each additional start requires another set of
  # parallel-tempering chains.
  multi_start_pip <- array(
    NA_real_, dim = c(p, r, n_starts),
    dimnames = list(colnames(eta_theta_y), colnames(X_G), paste0("start", seq_len(n_starts)))
  )
  multi_start_pip_mcse <- multi_start_pip
  multi_start_pip[, , best_start] <- matrix(
    pip_fit$pip, nrow = p, ncol = r, byrow = FALSE
  )
  multi_start_pip_mcse[, , best_start] <- matrix(
    pip_fit$pip_mcse, nrow = p, ncol = r, byrow = FALSE
  )
  if (isTRUE(pip_all_starts) && n_starts > 1L) {
    for (start in setdiff(seq_len(n_starts), best_start)) {
      start_fit <- start_fits[[start]]
      start_pip <- sample_integrated_pip_MOM(
        beta_init = start_fit$beta,
        y_star = start_fit$wd$y_star,
        X_star = start_fit$wd$X_star,
        nu = nu,
        sigma2 = start_fit$sigma2,
        tau0 = start_fit$tau0,
        tau1 = start_fit$tau1,
        pi = pi,
        weights_at_mode = start_fit$weights,
        n_iter = pip_mcmc_iter,
        burn = pip_mcmc_burn,
        thin = pip_mcmc_thin,
        n_chains = pip_mcmc_chains,
        n_temperatures = pip_mcmc_temperatures,
        min_inverse_temp = pip_mcmc_min_inverse_temp,
        swap_every = pip_mcmc_swap_every,
        large_jump_prob = pip_mcmc_large_jump_prob,
        seed = if (is.null(pip_mcmc_seed)) NULL else pip_mcmc_seed + start,
        verbose = verbose
      )
      multi_start_pip[, , start] <- matrix(
        start_pip$pip, nrow = p, ncol = r, byrow = FALSE
      )
      multi_start_pip_mcse[, , start] <- matrix(
        start_pip$pip_mcse, nrow = p, ncol = r, byrow = FALSE
      )
    }
  }

  # Conditional Laplace approximation for coefficient uncertainty. The
  # empirical-Bayes hyperparameters are held fixed, and any non-positive local
  # curvature is reported and regularized only for inversion.
  weighted_X <- sweep(X_star, 1L, w, "*")
  pc_final <- prior_components_MOM(beta, tau0, tau1, pi)
  beta_precision_raw <- crossprod(X_star, weighted_X) / sigma2 -
    diag(pc_final$logm_second, q)
  beta_precision_raw <- 0.5 * (beta_precision_raw + t(beta_precision_raw))
  beta_eigen <- eigen(beta_precision_raw, symmetric = TRUE)
  laplace_min_eigen <- min(beta_eigen$values)
  beta_eigenvalues <- pmax(beta_eigen$values, 1e-8)
  beta_cov_approx <- beta_eigen$vectors %*%
    diag(1 / beta_eigenvalues, q) %*% t(beta_eigen$vectors)
  beta_cov_approx <- 0.5 * (beta_cov_approx + t(beta_cov_approx))
  beta_sd_approx_vec <- sqrt(pmax(diag(beta_cov_approx), 0))

  # beta is ordered to match kronecker(X_G, diag(p)) (regressor-major,
  # parameter-minor), so the p x r reshape must be column-major (byrow = FALSE).
  beta_mat <- matrix(beta, nrow = p, ncol = r, byrow = FALSE)
  colnames(beta_mat) <- colnames(X_G)
  rownames(beta_mat) <- colnames(eta_theta_y)

  inclusion_mat <- matrix(inclusion_prob_vec, nrow = p, ncol = r, byrow = FALSE)
  colnames(inclusion_mat) <- colnames(X_G)
  rownames(inclusion_mat) <- colnames(eta_theta_y)

  em_responsibility_mat <- matrix(
    best$mode_responsibility, nrow = p, ncol = r, byrow = FALSE
  )
  colnames(em_responsibility_mat) <- colnames(X_G)
  rownames(em_responsibility_mat) <- colnames(eta_theta_y)

  inclusion_mcse_mat <- matrix(pip_fit$pip_mcse, nrow = p, ncol = r, byrow = FALSE)
  colnames(inclusion_mcse_mat) <- colnames(X_G)
  rownames(inclusion_mcse_mat) <- colnames(eta_theta_y)

  beta_posterior_mean_mat <- matrix(
    pip_fit$beta_posterior_mean, nrow = p, ncol = r, byrow = FALSE
  )
  colnames(beta_posterior_mean_mat) <- colnames(X_G)
  rownames(beta_posterior_mean_mat) <- colnames(eta_theta_y)

  beta_sd_approx_mat <- matrix(
    beta_sd_approx_vec, nrow = p, ncol = r, byrow = FALSE
  )
  colnames(beta_sd_approx_mat) <- colnames(X_G)
  rownames(beta_sd_approx_mat) <- colnames(eta_theta_y)

  # Fitted means per subject/parameter (on original scale)
  mu_hat_vec <- as.numeric(kronecker(X_G, diag(p)) %*% beta)
  mu_hat_mat <- t(matrix(mu_hat_vec, nrow = p, ncol = N))
  colnames(mu_hat_mat) <- colnames(eta_theta_y)

  list(
    beta_mat   = beta_mat,  # p x r empirical-Bayes posterior modes
    beta_vec   = beta,  # pr-vector of empirical-Bayes posterior modes
    beta_posterior_mean = beta_posterior_mean_mat,  # MCMC posterior means conditional on fitted hyperparameters
    beta_sd_approx = beta_sd_approx_mat,  # conditional Laplace posterior SDs
    beta_cov_approx = beta_cov_approx,
    laplace_min_eigen = laplace_min_eigen,
    laplace_regularized = laplace_min_eigen < 1e-8,
    inclusion  = inclusion_mat,  # p x r integrated posterior inclusion probabilities
    em_responsibility_at_mode = em_responsibility_mat,  # E-step quantity, not a PIP
    inclusion_mcse = inclusion_mcse_mat,  # batch-means Monte Carlo SEs
    pip_chain  = pip_fit$chain_pip,  # chain-specific integrated-PIP estimates
    pip_acceptance = pip_fit$acceptance,
    pip_tempered_acceptance = pip_fit$tempered_acceptance,
    pip_large_jump_acceptance = pip_fit$large_jump_acceptance,
    pip_swap_acceptance = pip_fit$swap_acceptance,
    pip_inverse_temperatures = pip_fit$inverse_temperatures,
    pip_max_chain_range = pip_fit$max_chain_pip_range,
    pip_retained_draws_per_chain = pip_fit$retained_draws_per_chain,
    sigma2     = sigma2,
    weights    = w,  # Student-t weights on whitened residuals
    alpha      = alpha,  # variance components
    tau0       = tau0,  # learned spike scale
    tau1       = tau1,  # learned slab scale
    objective = best$objective,
    objective_trace = best$objective_trace,
    objective_monotone = best$objective_monotone,
    q_start_trace = best$q_start_trace,
    q_end_trace = best$q_end_trace,
    q_monotone = best$q_monotone,
    converged = best$converged,
    stopping_reason = best$stopping_reason,
    iterations = best$iterations,
    final_parameter_change = best$final_parameter_change,
    final_objective_change = best$final_objective_change,
    beta_score_max = best$beta_score_max,
    alpha_optim_convergence = best$alpha_optim_convergence,
    multi_start_summary = start_summary,
    multi_start_objective_traces = lapply(start_fits, `[[`, "objective_trace"),
    multi_start_q_start_traces = lapply(start_fits, `[[`, "q_start_trace"),
    multi_start_q_end_traces = lapply(start_fits, `[[`, "q_end_trace"),
    multi_start_beta = multi_start_beta,
    multi_start_alpha = multi_start_alpha,
    multi_start_scales = multi_start_tau,
    multi_start_pip = multi_start_pip,
    multi_start_pip_mcse = multi_start_pip_mcse,
    best_start = best_start,
    whitening_method = whitening_method,
    whitening_factors = wd_final$B_list,
    weight_space = "coordinates defined by whitening_method and the input parameter order",
    parameter_order = if (is.null(colnames(eta_theta_y))) seq_len(p) else colnames(eta_theta_y),
    sigma2_estimated = estimate_sigma2,
    scale_identification = if (estimate_sigma2) "free scale; identification must be verified" else "sigma2 fixed",
    algorithm = "generalized EM with conditional maximization and post-fit parallel-tempering PIP integration",
    mu_hat     = mu_hat_mat,  # N x p fitted means
    N          = N,
    p          = p,
    r          = r,
    a0         = a0,
    b0         = b0,
    a1         = a1,
    b1         = b1
  )
}


# 4. Simulated example
# (Runs only when this file is executed directly, e.g. via Rscript;
# skipped when the file is sourced from another script such as the
# simulation study.)

# set.seed(2025)

# N  <- 40  # subjects
# p  <- 6  # DCM parameters per subject
# r  <- 2  # intercept + group

# # Group design: intercept + group indicator
# X_G <- cbind(
#   intercept = 1,
#   group     = rbinom(N, 1, 0.5)
# )

# # True group-level effects (p x r)
# beta_true <- matrix(0, nrow = p, ncol = r)
# colnames(beta_true) <- c("intercept", "group")

# beta_true[1, "intercept"] <-  0.5
# beta_true[2, "intercept"] <- -0.4
# beta_true[3, "group"]     <-  1.2
# beta_true[5, "group"]     <- -0.8

# # Construct mean matrix eta_theta_y_mean (N x p)
# eta_theta_y_mean <- X_G %*% t(beta_true)

# # Heavy-tailed noise via Student-t with df = 4
# nu_true      <- 4
# sigma2_true  <- 0.5
# lambda_true  <- rgamma(N * p, shape = nu_true / 2, rate = nu_true / 2)
# eps          <- rnorm(N * p, sd = sqrt(sigma2_true / lambda_true))
# Eps_mat      <- t(matrix(eps, nrow = p, ncol = N))

# eta_theta_y <- eta_theta_y_mean + Eps_mat

# # Add a couple of strong outliers
# eta_theta_y[5,  2] <- eta_theta_y[5,  2] +  6
# eta_theta_y[17, 4] <- eta_theta_y[17, 4] -  7

# # Simple choice of C_{theta|y,n}: all subjects share a diagonal covariance
# C_theta_y_list <- replicate(N, diag(0.2, p), simplify = FALSE)

# # Per-parameter variance-component basis: V_k has a single 1 in position (k,k),
# # so Sigma_b = sum_k alpha_k V_k is diagonal with entry k equal to alpha_k.
# # (Note: lapply(1:p, diag, nrow=p, ncol=p) would instead build k * I_p, making
# # the components collinear and non-identifiable.)
# V_list <- lapply(1:p, function(k) { Vk <- matrix(0, p, p); Vk[k, k] <- 1; Vk })

# # Fit the robust + sparse group DCM
# fit1 <- fit_group_dcm_t_MOM(
#   eta_theta_y    = eta_theta_y,
#   C_theta_y_list = C_theta_y_list,
#   X_G            = X_G,
#   V_list         = V_list,
#   nu             = 3,
#   tau0           = 0.05,  # initial spike scale (small)
#   tau1           = 1.0,  # initial slab scale (large)
#   pi             = 0.4,
#   max_iter       = 200,
#   tol            = 1e-6,
#   inner_sweeps   = 3,
#   verbose        = TRUE
# )

# cat("\nTrue beta (group-level effects):\n")
# print(round(beta_true, 3))

# cat("\nEstimated beta (beta_mat):\n")
# print(round(fit1$beta_mat, 3))

# cat("\nSquared error in beta:\n")
# print(sum((beta_true - fit1$beta_mat)^2))

# cat("\nPosterior inclusion probabilities (p x r):\n")
# print(round(fit1$inclusion, 3))

# cat("\nEstimated alpha (variance components):\n")
# print(round(fit1$alpha, 3))

# cat("\nEstimated tau0 and tau1 (spike and slab scales):\n")
# print(c(tau0 = fit1$tau0, tau1 = fit1$tau1))


# Sample use (this file defines the estimator; the analysis scripts source it):
#   source("rsDCM.R")
#   fit <- fit_group_dcm_t_MOM(eta_theta_y = eta, X_G = matrix(1, nrow(eta), 1),
#                              C_theta_y_list = Cp, V_list = V_list,
#                              nu = 3, pi = 0.5, tau0 = 0.05, tau1 = 1.0)
