# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Godfred Arhin and Nilotpal Sanyal
# real_data_analysis.R
# Real-data analysis on the NARPS data (48 subjects; 22 DCM parameters over
# vmPFC, ventral striatum, amygdala and anterior insula), run top to bottom:
#   1. robust + sparse group DCM fit (intercept and covariate models) and contrasts
#   2. standard SPM PEB baseline
#   3. method comparison and figures
#   4. covariate-model comparison and figures
#   5. sign tests
#   6. annotated connectivity figure
#   7. LaTeX tables (real_data_tables.tex)
# Robustness nu = 3, sparsity pi = 0.5 (set via NU / PI below).
# The corrected analysis reads the deposited group-level input and participant
# metadata. Regenerated outputs are written separately from deposited results.
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
res2_path <- file.path(script_dir, "res2_corrected.rds")
participants_path <- file.path(script_dir, "included_48_participants.csv")
results_dir <- Sys.getenv(
  "RSDCM_REAL_RESULTS_DIR",
  unset = file.path(script_dir, "results_regenerated")
)

# Engines: the proposed method and the ported SPM PEB baseline.
group_method_file <- file.path(repository_dir, "rsDCM.R")
spm_port_file <- file.path(script_dir, "spm_dcm_peb_port.R")
required_input_files <- c(group_method_file, spm_port_file, res2_path, participants_path)
missing_input_files <- required_input_files[!file.exists(required_input_files)]
if (length(missing_input_files)) {
  stop("Missing required file(s):\n", paste(missing_input_files, collapse = "\n"))
}

source(group_method_file, chdir = TRUE)
source(spm_port_file, chdir = TRUE)
set.seed(20260717)

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
overwrite_results <- tolower(Sys.getenv("RSDCM_OVERWRITE_RESULTS", "false")) %in%
  c("true", "1", "yes")
if (file.exists(file.path(results_dir, "real_data_fits.rds")) && !overwrite_results) {
  stop(
    "Corrected group-analysis results already exist in: ", results_dir,
    "\nMove or rename that directory, or set RSDCM_OVERWRITE_RESULTS=true explicitly."
  )
}
setwd(results_dir)

# prior settings: nu (Student-t robustness), pi (spike-and-slab sparsity).
# Set NU / PI in the workspace before source()-ing to override these defaults.
if (!exists("NU")) NU <- 3L
if (!exists("PI")) PI <- 0.5

# 0b. Corrected data and shared helpers
res2   <- readRDS(res2_path)
eta    <- res2$eta; Cp <- res2$Cp_list; pnames <- res2$parameter_names
N      <- nrow(eta); p <- ncol(eta)
stopifnot(N == 48L, p == 22L, length(Cp) == N, ncol(Cp[[1]]) == p)
if (is.null(res2$export_information) ||
    !identical(res2$export_information$source,
               "corrected run-level DCM and within-subject PEB analysis")) {
  stop("The selected R object is not the validated corrected PEB export.")
}

# This covariance-only object is the input for the pending NARPS covariance
# sensitivity analysis in the simulation workflow.
saveRDS(Cp, "narps_covariances.rds", compress = "xz")

sid   <- sub(".*(sub-[0-9]+).*", "\\1", res2$subject_names)
# Earlier input retained for reference.
# parts <- read.delim("participants.tsv", stringsAsFactors = FALSE)
parts <- read.csv(participants_path, stringsAsFactors = FALSE)
idx   <- match(sid, parts$participant_id)
if (any(is.na(idx))) stop("Unmatched subjects: ", paste(sid[is.na(idx)], collapse = ", "))
meta  <- parts[idx, ]
genderM <- as.integer(meta$gender == "M"); age_z <- as.numeric(scale(meta$age))
cat(sprintf("Matched %d subjects. gender %d M / %d F; age mean %.1f (sd %.1f).\n",
            N, sum(genderM), N - sum(genderM), mean(meta$age), sd(meta$age)))

V_list <- lapply(seq_len(p), function(k) { Vk <- matrix(0, p, p); Vk[k, k] <- 1; Vk })

# parameter name -> A/B layout + gain/loss modulatory indices
parse_idx <- function(nm) {
  if (grepl("^A", nm)) { m <- as.integer(strsplit(gsub("[A()]", "", nm), ",")[[1]]); list("A", m[1], m[2], NA) }
  else { m <- as.integer(strsplit(gsub("[B()]", "", nm), ",")[[1]]); list("B", m[1], m[2], m[3]) }
}
pinfo   <- lapply(pnames, parse_idx)
is_gain <- vapply(pinfo, function(z) z[[1]] == "B" && z[[4]] == 1, logical(1))
is_loss <- vapply(pinfo, function(z) z[[1]] == "B" && z[[4]] == 2, logical(1))
gi <- which(is_gain); li <- which(is_loss)
regs <- c("vmPFC", "vStr", "amyg", "aIns")
vec_to_matrices <- function(v) {
  A <- matrix(NA, 4, 4); Bg <- matrix(0, 4, 4); Bl <- matrix(0, 4, 4)
  for (k in seq_along(v)) { z <- pinfo[[k]]
    if (z[[1]] == "A") A[z[[2]], z[[3]]] <- v[k]
    else if (z[[4]] == 1) Bg[z[[2]], z[[3]]] <- v[k] else Bl[z[[2]], z[[3]]] <- v[k] }
  dimnames(A) <- dimnames(Bg) <- dimnames(Bl) <- list(regs, regs)
  list(A = A, Bgain = Bg, Bloss = Bl)
}
mk_png <- function(f, expr, w = 1100, h = 560) {
  tryCatch({ png(f, width = w, height = h, res = 110); on.exit(dev.off()); eval(expr) },
           error = function(e) cat("(figure", f, "skipped:", conditionMessage(e), ")\n"))
}

# Robust + sparse group DCM fit (nu = NU, pi = PI)
common <- list(C_theta_y_list = Cp, V_list = V_list, nu = NU, pi = PI,
               tau0 = 0.05, tau1 = 1.0, max_iter = 500, tol = 1e-6, verbose = FALSE)
X_int <- matrix(1, N, 1, dimnames = list(NULL, "intercept"))
X_cov <- cbind(intercept = 1, genderM = genderM, age = age_z)

cat("\n[1] Fitting robust+sparse intercept and covariate models ...\n")
fit_int <- do.call(fit_group_dcm_t_MOM, c(list(eta_theta_y = eta, X_G = X_int), common))
fit_cov <- do.call(fit_group_dcm_t_MOM, c(list(eta_theta_y = eta, X_G = X_cov), common))
rownames(fit_int$beta_mat) <- rownames(fit_int$inclusion) <- pnames
rownames(fit_cov$beta_mat) <- rownames(fit_cov$inclusion) <- pnames

# Laplace posterior for the intercept coefficients -> Bayesian contrasts
laplace_cov <- function(fit, X_G) {
  Sigma_b <- build_Sigma_b(fit$alpha, V_list, p); wd <- whiten_data(eta, Cp, X_G, Sigma_b)
  prior_prec <- -prior_components_MOM(as.numeric(fit$beta_vec), fit$tau0, fit$tau1, pi = PI)$logm_second
  solve(crossprod(wd$X_star, fit$weights * wd$X_star) / fit$sigma2 + diag(prior_prec, length(fit$beta_vec)))
}
Cov_beta <- laplace_cov(fit_int, X_int); beta_hat <- as.numeric(fit_int$beta_vec)
S <- 20000
L <- tryCatch(chol(Cov_beta), error = function(e) chol(Cov_beta + 1e-8 * diag(nrow(Cov_beta))))
draws <- matrix(beta_hat, S, length(beta_hat), byrow = TRUE) + matrix(rnorm(S * length(beta_hat)), S) %*% L
contrast_report <- function(sel, fun, label, direction) {
  vals <- fun(draws[, sel, drop = FALSE]); prob <- if (direction == ">0") mean(vals > 0) else mean(vals < 0)
  cat(sprintf("  %-38s P(%s)=%.3f  mean=%+.3f\n", label, direction, prob, mean(vals)))
  c(prob = prob, mean = mean(vals), var = var(vals))
}
int_mats <- vec_to_matrices(fit_int$beta_mat[, "intercept"])
pip_mats <- vec_to_matrices(fit_int$inclusion[, "intercept"])
# SPM stores diagonal A parameters as dimensionless log scalings.
# Convert them to physical self-inhibition in Hz for display.
A_internal <- int_mats$A
A_hz <- A_internal
diag(A_hz) <- -0.5 * exp(diag(A_internal))
cat("--- Intercept intrinsic connectivity A [Hz] ---\n"); print(round(A_hz, 3))
cat("--- Bayesian contrasts (Laplace posterior) ---\n")
h1 <- contrast_report(gi, rowMeans, "H1: gain modulation (pooled) > 0", ">0")
h2 <- contrast_report(li, rowMeans, "H2: loss modulation (pooled) < 0", "<0")
h3 <- contrast_report(c(gi, li),
        function(M) rowMeans(abs(M[, seq_along(li) + length(gi), drop = FALSE])) -
                    rowMeans(abs(M[, seq_along(gi), drop = FALSE])), "H3: |loss| > |gain|", ">0")

w_mat <- matrix(fit_int$weights, nrow = N, ncol = p, byrow = TRUE); subj_w <- rowMeans(w_mat)
cat(sprintf("  subjects with mean Student-t weight < 0.6: %d of %d\n", sum(subj_w < 0.6), N))

saveRDS(list(fit_int = fit_int, fit_cov = fit_cov, Cov_beta = Cov_beta, meta = meta,
             contrasts = rbind(H1 = h1, H2 = h2, H3 = h3)), "real_data_fits.rds")
write.csv(data.frame(parameter = pnames,
            est_intercept = fit_int$beta_mat[, "intercept"], pip_intercept = fit_int$inclusion[, "intercept"],
            est_genderM = fit_cov$beta_mat[, "genderM"], pip_genderM = fit_cov$inclusion[, "genderM"],
            est_age = fit_cov$beta_mat[, "age"], pip_age = fit_cov$inclusion[, "age"]),
          "real_data_estimates.csv", row.names = FALSE)
ij <- expand.grid(i = 1:4, j = 1:4)
write.csv(rbind(
  data.frame(panel = "A",   ij, region_i = regs[ij$i], region_j = regs[ij$j], value = as.vector(A_hz)),
  data.frame(panel = "PIP", ij, region_i = regs[ij$i], region_j = regs[ij$j], value = as.vector(pip_mats$A))),
  "real_A_matrices.csv", row.names = FALSE)
heat <- function(M, main, zlim = NULL) {
  if (is.null(zlim)) zlim <- c(-1, 1) * max(abs(M), na.rm = TRUE) + 1e-9
  image(t(M[nrow(M):1, ]), axes = FALSE, main = main, zlim = zlim, col = colorRampPalette(c("blue", "white", "red"))(64))
  axis(1, at = seq(0, 1, length = ncol(M)), labels = colnames(M), tick = FALSE, cex.axis = 0.8)
  axis(2, at = seq(0, 1, length = nrow(M)), labels = rev(rownames(M)), tick = FALSE, las = 1, cex.axis = 0.8) }
mk_png("real_fig1_intercept_A.png", quote({ par(mfrow = c(1, 2), mar = c(3, 4, 3, 1))
  heat(A_hz, "Intercept A (estimate, Hz)"); heat(pip_mats$A, "Intercept A (PIP)", zlim = c(0, 1)) }), w = 1000, h = 520)
mk_png("real_fig2_subject_weights.png", quote({ par(mar = c(4, 4, 3, 1))
  plot(sort(subj_w), pch = 19, ylim = c(min(0.5, min(subj_w)), max(1.3, max(subj_w))),
       xlab = "subject (sorted)", ylab = "mean Student-t weight",
       main = sprintf("Subject robustness weights (N = %d)", length(subj_w)))
  abline(h = 1.0, lty = 3, col = "grey60"); abline(h = 0.6, lty = 2, col = "firebrick"); grid() }), w = 1000, h = 520)

# Standard SPM PEB baseline
cat("\n[2] Running standard SPM PEB (intercept + covariate) ...\n")
pebs <- lapply(res2$PEB_list, .normalize_mat_peb)
run_peb <- function(X, Xnames) {
  M <- list(X = X, Xnames = as.list(Xnames), Q = "single", alpha = 1, beta = 16, maxit = 64)
  spm_dcm_peb_run(spm_dcm_peb_prepare(P = pebs, M = M, field = seq_len(p)), verbose = FALSE)$PEB
}
PEB_int <- run_peb(matrix(1, N, 1), "Intercept")
PEB_cov <- run_peb(cbind(Intercept = 1, genderM = genderM, age = age_z), c("Intercept", "genderM", "age"))
Ep_int <- matrix(as.numeric(PEB_int$Ep), ncol = ncol(PEB_int$M$X))[, 1]
sd_int <- sqrt(diag(PEB_int$Cp))[seq_len(length(Ep_int))]
Ep_cov <- matrix(as.numeric(PEB_cov$Ep), ncol = ncol(PEB_cov$M$X)); colnames(Ep_cov) <- c("Intercept", "genderM", "age")
cat(sprintf("  SPM PEB strong effects (|Ep|>3SD): %d of %d\n", sum(abs(Ep_int) > 3 * sd_int), p))
saveRDS(list(PEB_int = PEB_int, PEB_cov = PEB_cov), "spm_peb_fits.rds")
write.csv(data.frame(parameter = pnames, Ep_intercept = Ep_int, SD_intercept = sd_int,
                     Ep_genderM = Ep_cov[, "genderM"], Ep_age = Ep_cov[, "age"]),
          "spm_peb_estimates.csv", row.names = FALSE)

# Method comparison (robust vs SPM PEB)
cat("\n[3] Comparing robust+sparse vs SPM PEB (intercept) ...\n")
cmp <- data.frame(parameter = pnames, PEB_Ep = Ep_int, PEB_SD = sd_int,
                  PEB_strong = abs(Ep_int) > 3 * sd_int,
                  Robust_est = fit_int$beta_mat[, "intercept"], Robust_PIP = fit_int$inclusion[, "intercept"])
cmp$shrinkage <- cmp$Robust_est / cmp$PEB_Ep                # robust / PEB (signed ratio)
cmp$same_sign <- sign(cmp$Robust_est) == sign(cmp$PEB_Ep)
is_A <- grepl("^A", pnames)
cat(sprintf("  sign agreement: %d/%d (%.0f%%); PEB strong %d; Robust selected %d\n",
            sum(sign(cmp$Robust_est) == sign(cmp$PEB_Ep)), nrow(cmp),
            100 * mean(sign(cmp$Robust_est) == sign(cmp$PEB_Ep)), sum(cmp$PEB_strong), sum(cmp$Robust_PIP > 0.5)))
write.csv(cmp, "method_comparison.csv", row.names = FALSE)
mk_png("cmp_fig1_estimates.png", quote({ par(mar = c(7, 4, 3, 1))
  bp <- barplot(rbind(PEB = cmp$PEB_Ep, Robust = cmp$Robust_est), beside = TRUE, names.arg = pnames,
                las = 2, cex.names = 0.6, col = c("grey70", "steelblue"), ylab = "group-average estimate",
                main = expression("Standard " * italic("SPM-PEB") * " vs " * italic("Proposed") * " (intercept)"))
  legend("topright", c("SPM-PEB", "Proposed"), fill = c("grey70", "steelblue"), bty = "n", text.font = 3)
  sel <- which(cmp$Robust_PIP > 0.5); if (length(sel)) points(bp[2, sel], cmp$Robust_est[sel], pch = 8, col = "red")
  legend("bottomright", expression(italic("Proposed") * " PIP>0.5"), pch = 8, col = "red", bty = "n") }), w = 1100, h = 620)
mk_png("cmp_fig2_pip_vs_peb.png", quote({ par(mar = c(4, 4, 3, 1))
  plot(abs(cmp$PEB_Ep / cmp$PEB_SD), cmp$Robust_PIP, pch = 19, col = ifelse(is_A, "black", "darkorange"),
       xlab = expression(italic("SPM-PEB") * " evidence  |Ep|/SD"), ylab = expression(italic("Proposed") * " posterior inclusion prob."),
       main = expression("What the " * italic("Proposed") * " model selects vs " * italic("SPM-PEB") * " evidence"))
  abline(h = 0.5, lty = 2); abline(v = 3, lty = 3, col = "grey")
  legend("right", c("A (intrinsic)", "B (modulatory)"), pch = 19, col = c("black", "darkorange"), bty = "n") }), w = 1100, h = 620)

# Covariate-model comparison
cat("\n[4] Comparing covariate effects (gender / age) ...\n")
peb_block <- function(k) { i <- ((k - 1) * p + 1):(k * p); list(Ep = as.numeric(PEB_cov$Ep)[i], SD = sqrt(pmax(diag(PEB_cov$Cp[i, i, drop = FALSE]), 0))) }
rob_block <- function(k, nm) { i <- ((k - 1) * p + 1):(k * p); list(Ep = as.numeric(fit_cov$beta_vec)[i], PIP = fit_cov$inclusion[, nm]) }
g_pb <- peb_block(2); g_rb <- rob_block(2, "genderM"); a_pb <- peb_block(3); a_rb <- rob_block(3, "age")
cat(sprintf("  gender: PEB strong %d, Robust selected %d ; age: PEB strong %d, Robust selected %d\n",
            sum(abs(g_pb$Ep) > 3 * g_pb$SD), sum(g_rb$PIP > 0.5), sum(abs(a_pb$Ep) > 3 * a_pb$SD), sum(a_rb$PIP > 0.5)))
write.csv(data.frame(parameter = pnames,
            gender_PEB = g_pb$Ep, gender_PEB_SD = g_pb$SD, gender_Robust = g_rb$Ep, gender_PIP = g_rb$PIP,
            age_PEB = a_pb$Ep, age_PEB_SD = a_pb$SD, age_Robust = a_rb$Ep, age_PIP = a_rb$PIP),
          "covariate_comparison.csv", row.names = FALSE)
bar_cov <- function(pb, rb, ttl) { bp <- barplot(rbind(PEB = pb$Ep, Robust = rb$Ep), beside = TRUE, names.arg = pnames,
    las = 2, cex.names = 0.6, col = c("grey70", "steelblue"), ylab = "effect", main = ttl)
  legend("topleft", c("SPM-PEB", "Proposed"), fill = c("grey70", "steelblue"), bty = "n", cex = 0.8, text.font = 3)
  sel <- which(rb$PIP > 0.5); if (length(sel)) points(bp[2, sel], rb$Ep[sel], pch = 8, col = "red") }
mk_png("cmp_fig3_gender.png", quote({ par(mar = c(7, 4, 3, 1)); bar_cov(g_pb, g_rb, expression("Gender effect (M vs F): " * italic("SPM-PEB") * " vs " * italic("Proposed"))) }))
mk_png("cmp_fig4_age.png",    quote({ par(mar = c(7, 4, 3, 1)); bar_cov(a_pb, a_rb, expression("Age effect (per SD): " * italic("SPM-PEB") * " vs " * italic("Proposed"))) }))

# Sign tests (concordance + shrinkage)
cat("\n[5] Sign tests ...\n")
concordance <- function(a, b, label, keep = NULL) {
  if (!is.null(keep)) { a <- a[keep]; b <- b[keep] }
  nz <- a != 0 & b != 0; a <- a[nz]; b <- b[nz]; n <- length(a); k <- sum(sign(a) == sign(b))
  if (n == 0L) return(data.frame(comparison = label, n = 0L, agree = NA, agree_rate = NA, p_value = NA))
  data.frame(comparison = label, n = n, agree = k, agree_rate = round(k / n, 3),
             p_value = signif(binom.test(k, n, 0.5, "greater")$p.value, 3)) }
cmp_cred <- cmp$PEB_strong | cmp$Robust_PIP > 0.5
gen_cred <- (abs(g_pb$Ep) > 3 * g_pb$SD) | g_rb$PIP > 0.5
age_cred <- (abs(a_pb$Ep) > 3 * a_pb$SD) | a_rb$PIP > 0.5
tabA <- rbind(
  concordance(cmp$PEB_Ep, cmp$Robust_est, "intercept: all 22"),
  concordance(cmp$PEB_Ep, cmp$Robust_est, "intercept: credible only", keep = cmp_cred),
  concordance(g_pb$Ep, g_rb$Ep, "gender: all 22"),
  concordance(g_pb$Ep, g_rb$Ep, "gender: credible only", keep = gen_cred),
  concordance(a_pb$Ep, a_rb$Ep, "age: all 22"),
  concordance(a_pb$Ep, a_rb$Ep, "age: credible only", keep = age_cred))
print(tabA, row.names = FALSE)
shrink_test <- function(a, b, label) { d <- abs(a) - abs(b); d <- d[d != 0]; n <- length(d); k <- sum(d > 0)
  data.frame(comparison = label, n = n, robust_smaller = k, rate = round(k / n, 3),
             p_value = signif(binom.test(k, n, 0.5, "greater")$p.value, 3)) }
tabB <- rbind(shrink_test(cmp$PEB_Ep, cmp$Robust_est, "intercept"),
              shrink_test(g_pb$Ep, g_rb$Ep, "gender"), shrink_test(a_pb$Ep, a_rb$Ep, "age"))
write.csv(tabA, "sign_test_concordance.csv", row.names = FALSE)
write.csv(tabB, "sign_test_shrinkage.csv", row.names = FALSE)

# Annotated connectivity figure
cat("\n[6] Annotated connectivity figure ...\n")
mat <- read.csv("real_A_matrices.csv", stringsAsFactors = FALSE)
Amat <- matrix(NA, 4, 4); Pmat <- matrix(NA, 4, 4)
for (r in seq_len(nrow(mat))) { if (mat$panel[r] == "A") Amat[mat$i[r], mat$j[r]] <- mat$value[r] else Pmat[mat$i[r], mat$j[r]] <- mat$value[r] }
rdbu <- colorRampPalette(c("#2166AC", "#4393C3", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#D6604D", "#B2182B"))(256)
reds <- colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#67000D"))(256)
draw_hm <- function(M, cols, zlim, main, fmt, diverging) {
  par(mar = c(4, 5, 3, 1)); image(1:4, 1:4, t(M)[, 4:1], col = cols, zlim = zlim, axes = FALSE,
        xlab = "source region (j, from)", ylab = "target region (i, to)", main = main)
  axis(1, 1:4, regs, tick = FALSE); axis(2, 1:4, regs[4:1], tick = FALSE, las = 1); box()
  for (i in 1:4) for (j in 1:4) { dark <- if (diverging) abs(M[i, j]) > 0.6 * zlim[2] else (M[i, j] - zlim[1]) / diff(zlim) > 0.6
    text(j, 4 - i + 1, sprintf(fmt, M[i, j]), col = if (dark) "white" else "black") } }
draw_cb <- function(zlim, cols, lab) { par(mar = c(4, 0.5, 3, 3.5))
  image(1, seq(zlim[1], zlim[2], length.out = length(cols)), matrix(seq_along(cols), 1), col = cols, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1); mtext(lab, side = 4, line = 2.3, cex = 0.85); box() }
vlim <- max(abs(Amat), na.rm = TRUE)
mk_png("real_fig1_A_annotated.png", quote({ layout(matrix(1:4, 1), widths = c(5, 1.1, 5, 1.1))
  draw_hm(Amat, rdbu, c(-vlim, vlim), "Group-average intrinsic connectivity  A  (Hz)", "%+.2f", TRUE); draw_cb(c(-vlim, vlim), rdbu, "effect (Hz)")
  draw_hm(Pmat, reds, c(0, 1), "Posterior inclusion probability (PIP)", "%.2f", FALSE); draw_cb(c(0, 1), reds, "PIP") }), w = 1430, h = 550)

# LaTeX tables are written by make_real_data_tables.R.

cat("\nDone. Wrote: real_data_fits.rds, real_data_estimates.csv, real_A_matrices.csv,\n",
    "spm_peb_fits.rds, spm_peb_estimates.csv, method_comparison.csv, covariate_comparison.csv,\n",
    "sign_test_concordance.csv, sign_test_shrinkage.csv, narps_covariances.rds,\n",
    "and the figures in:\n", results_dir, "\n", sep = "")


# Sample run:
#   Command line:
#     Rscript real_data/real_data_analysis.R
#
#   Interactive R:
#     setwd("real_data")
#     source("real_data_analysis.R")
#
# NU and PI may be set before source() to override their defaults.
