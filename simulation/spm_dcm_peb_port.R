# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Godfred Arhin and Nilotpal Sanyal
################################################################################
# Hierarchical group-PEB port
# Third-level PEB-of-PEBs (group DCM) — merged version.
#
# Combines:
#   * a generic loader that reads both .rds and
#     MATLAB .mat subject PEBs (.mat via R.matlab::readMat), and
#   * a third-level orchestrator with custom Xnames,
#     tidy summary data.frame, missing-file handling, auto dir.create on save.
#
# Design choice: .mat PEBs are normalized to native R-list shape at load time
# (.normalize_mat_peb), so a SINGLE R-native prepare/run path serves both input
# formats (no duplicated _2 code path, no [,,1] indexing scattered downstream).
#
# NOTE: the six robust Student-t / MOM group-fitting functions
# (fit_group_dcm_t_MOM, prior_components_MOM, build_Sigma_b, alpha_grad_hess,
#  update_tau_scales, whiten_data) are intentionally NOT included here.
#
# Requires: R.matlab (only when loading .mat files).
################################################################################

source("exact_core_final_2_optimized.R")
source("helpers_final_optimized.R")

# ---- normalize a MATLAB-struct PEB (from R.matlab::readMat) to native shape ---
# readMat wraps structs as 3-D arrays. The estimator only needs the top-level
# struct and its $M sub-struct unwrapped; all other fields are plain after that.
.normalize_mat_peb <- function(x) {
  if (!is.null(dim(x)) && length(dim(x)) >= 3) x <- x[, , 1]
  if (!is.null(x$M) && !is.null(dim(x$M)) && length(dim(x$M)) >= 3) x$M <- x$M[, , 1]
  x
}

detect_peb_files <- function(peb_dir) {
  files <- list.files(peb_dir,
                      pattern = "PEB_sub.*\\.(rds|mat)$",
                      full.names = TRUE,
                      ignore.case = TRUE)
  
  if (!length(files)) stop("No PEB files found.")
  
  # Extract subject IDs (flexible)
  subjects <- sub(".*sub[-_]?([0-9]+).*", "\\1", basename(files))
  subjects <- as.integer(subjects)
  
  valid <- !is.na(subjects)
  
  list(
    files = files[valid],
    subjects = subjects[valid]
  )
}


# ---- load subject PEBs (.rds or .mat) and normalize to native shape ----------
load_pebs_generic <- function(peb_dir, subjects = NULL) {

  D <- detect_peb_files(peb_dir)

  if (!is.null(subjects)) {
    keep <- D$subjects %in% subjects
    D$files    <- D$files[keep]
    D$subjects <- D$subjects[keep]
  }

  if (!length(D$files)) stop("No matching subjects found.")

  pebs <- lapply(D$files, function(f) {
    if (grepl("\\.rds$", f, ignore.case = TRUE)) {
      readRDS(f)                                   # native R PEB
    } else {
      m <- R.matlab::readMat(f, sparseMatrix = "matrix")
      m$PEB                                         # MATLAB PEB struct
    }
  })

  # unify both formats to native R-list shape for the single prepare/run path
  pebs <- lapply(pebs, .normalize_mat_peb)

  names(pebs) <- paste0("sub", sprintf("%03d", D$subjects))
  list(pebs = pebs, subjects = D$subjects, files = D$files)
}


# === Core PEB estimator (R-native; serves both .rds and normalized .mat) ===
# # ======================================================================
# # Part 1: spm_dcm_peb_prepare
# # - Gathers first-level densities
# # - Selects parameter indices (field -> q)
# # - Projects to rank-reduced subspace (U)
# # - Builds second-level priors (bE, bC), hyperpriors (gE,gC),
# #   precision components Q, design sizes, and initial b,g
# # Returns a list ("prep") used by Part 2
# # ======================================================================
spm_dcm_peb_prepare <- function(P, M = list(), field = c("A","B")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  
  if (is.null(M)) {
    M <- list(X = matrix(1, nrow = length(P), ncol = 1))
  }
  if (is.numeric(M)) {
    M <- list(X = M)
  }
  if (!("X" %in% names(M)) || is.null(M$X)) {
    M$X <- matrix(1, nrow = length(P), ncol = 1)
  }

  # ---- normalize field ---------------------------------------------------
  if (missing(field) || is.null(field)) field <- c("A","B")
  if (is.list(field)) field <- unlist(field, use.names = FALSE)

  if (!(is.character(field) || is.numeric(field))) {
    stop("'field' must be character vector (e.g., c('A','B')), the string 'all', or numeric indices")
  }

  # basic inputs
  Ns   <- length(P)
  DCM <- P[[1]]
  if (is.null(DCM$M) || is.null(DCM$M$pE))
    stop("P[[1]]$M$pE is required")

  pE <- DCM$M$pE
  np_total <- length(spm_vec(pE))

  # ---- resolve 'field' names if needed -----------------------------------
  if (is.character(field) && length(field) == 1L && tolower(field) == "all") {
    if (is.null(names(pE)))
      stop("field='all' requires named blocks in M$pE (e.g., A/B/C/...)")
    field <- names(pE)
  }
  
  maxit <- ifelse("maxit" %in% names(M), M$maxit, 64)

  # # ---- build q (indices in spm_vec(M$pE)) --------------------------------
  # if (is.numeric(field)) {
  #   # numeric indexing directly into spm_vec(pE)
  #   q <- as.integer(field)
  #   if (any(q < 1L | q > np_total)) {
  #     stop(sprintf("numeric 'field' indices must be in 1..%d", np_total))
  #   }
  # } else {
  #   # character block names -> concatenate the ranges that correspond to each block
  #   nm <- names(pE)
  #   if (is.null(nm)) {
  #     stop("character 'field' requires named blocks in M$pE")
  #   }
  #   # confirm all requested names exist
  #   missing_names <- setdiff(field, nm)
  #   if (length(missing_names))
  #     stop(sprintf("Unknown field(s): %s", paste(missing_names, collapse = ", ")))
  # 
  #   # walk through pE in list order to reproduce spm_vec() ordering
  #   q <- integer(0)
  #   offset <- 0L
  #   for (k in seq_along(pE)) {
  #     len_k <- length(spm_vec(pE[[k]]))
  #     idx_k <- seq.int(offset + 1L, offset + len_k)
  #     if (nm[k] %in% field) q <- c(q, idx_k)
  #     offset <- offset + len_k
  #   }
  #   q <- as.integer(q)
  #   if (!length(q)) stop("Requested 'field' produced an empty index set")
  # }
  
  
  if ("bC" %in% names(M) && Ns > 1) {
    q <- spm_find_pC(M$bC, M$bE, field)$i
  } else if (is.numeric(field)) {
    q <- field
  } else {
    q <- spm_find_pC(DCM, field)$i
  }

  # Names
  is_peb_of_pebs <- !is.null(DCM$Pnames)
  if (is_peb_of_pebs) {
    Pstr <- character(0)
    for (i in seq_along(DCM$Xnames)) {
      Pstr <- c(Pstr, paste0(DCM$Xnames[[i]], ": ", DCM$Pnames))
    }
  } else {
    # Let your helper produce labels; if it fails, fall back to generic
    Pstr <- tryCatch(spm_fieldindices(DCM$M$pE, field),
                     error = function(e) paste0("P", q))
    if (length(q) == 1) Pstr <- list(Pstr)
  }

  # Pull per-subject densities
  pE <- vector("list", Ns)
  pC <- vector("list", Ns)
  qE <- vector("list", Ns)
  qC <- vector("list", Ns)
  iF <- numeric(Ns)

  for (i in seq_len(Ns)) {
    D <- P[[i]]
    pC[[i]] <- if (is.list(D$M$pC)) diag(spm_vec(D$M$pC)) else as.matrix(D$M$pC)
    pE[[i]] <- spm_vec(D$M$pE)
    qE[[i]] <- spm_vec(D$Ep)
    qC[[i]] <- as.matrix(D$Cp)
    iF[i]   <- as.numeric(D$F)
    if (i == 1) Ne <- length(pE[[i]])
    if (length(pE[[i]]) != Ne) stop("All DCMs must have same parameterisation.")
  }

  # Average priors over subjects (on selected indices q)
  PE <- numeric(length(q))
  PC <- matrix(0, length(q), length(q))
  for (i in seq_len(Ns)) {
    PE <- PE + pE[[i]][q]
    PC <- PC + pC[[i]][q, q, drop = FALSE]
  }
  PE <- PE / Ns
  PC <- PC / Ns

  # Rank-reduction (projector U)
  U <- if (Ns > 1) as.matrix(spm_svd(PC)$U) else diag(length(q))

  # Project densities and apply shrinkage
  for (i in seq_len(Ns)) {
    pE[[i]] <- as.vector(t(U) %*% pE[[i]][q])
    pC[[i]] <- t(U) %*% pC[[i]][q, q, drop = FALSE] %*% U
    qE[[i]] <- as.vector(t(U) %*% qE[[i]][q])
    qC[[i]] <- t(U) %*% qC[[i]][q, q, drop = FALSE] %*% U
    if (Ns > 1) {
      qC[[i]] <- spm_inv(spm_inv(qC[[i]]) + spm_inv(pC[[i]]) / 16)
    }
  }

  # Second-level design
  if (Ns > 1) {
    W <- t(U) %*% (M$W %||% diag(length(q))) %*% U
    X <- M$X
  } else {
    W <- M$W %||% diag(length(q))
    X <- matrix(1, 1, 1)
  }

  # Covariance component option
  OPTION <- "single"; Qcells <- list()
  if (!is.null(M$Q)) {
    if (is.list(M$Q) && is.numeric(M$Q[[1]])) { OPTION <- "manual"; Qcells <- M$Q }
    else { OPTION <- M$Q }
  }

  # Hyperpriors and base priors
  alpha <- M$alpha %||% 1
  beta  <- M$beta  %||% 16
  Nx <- ncol(X)
  Xnames <- M$Xnames %||% lapply(seq_len(Nx), function(i) sprintf("Covariate %d", i))

  # 3rd-level priors
  Np <- length(q)
  if (!is.null(M$bE)) {
    bE0 <- spm_vec(M$bE); if (Ns > 1 && length(bE0) > Np) bE0 <- bE0[q]
  } else bE0 <- PE
  if (!is.null(M$bC)) {
    bC0 <- if (is.list(M$bC)) diag(spm_vec(M$bC)) else as.matrix(M$bC)
    if (Ns > 1 && nrow(bC0) > Np) bC0 <- bC0[q, q, drop = FALSE]
  } else bC0 <- PC / alpha
  if (!is.null(M$pC)) {
    pC2 <- if (is.list(M$pC)) diag(spm_vec(M$pC)) else as.matrix(M$pC)
    if (Ns > 1 && nrow(pC2) > Np) pC2 <- pC2[q, q, drop = FALSE]
  } else if (isTRUE(beta == 0)) {
    qEmat <- do.call(cbind, qE)
    pC2 <- diag(apply(qEmat, 1, var))
  } else {
    pC2 <- PC / beta
  }

  # Base precision in projected space
  pQ <- spm_inv(t(U) %*% pC2 %*% U)

  # Build Q components in projected space
  Q <- list()
  if (OPTION == "single") {
    Q <- list(pQ)
  } else if (OPTION == "fields") {
    pq <- spm_inv(pC2)
    k <- 1; off <- 0L
    for (nm in names(DCM$M$pE)) {
      len <- length(spm_vec(DCM$M$pE[[nm]]))
      if (nm %in% field) {
        ids <- off + seq_len(len)
        ids <- which(ids %in% q)
        if (length(ids)) {
          Qi <- matrix(0, Np, Np)
          Qi[ids, ids] <- pq[ids, ids, drop = FALSE]
          Q[[k]] <- t(U) %*% Qi %*% U; k <- k + 1
        }
      }
      off <- off + len
    }
  } else if (OPTION == "all") {
    pq <- spm_inv(pC2)
    k <- 1
    for (i in seq_len(Np)) {
      Qi <- matrix(0, Np, Np); Qi[i, i] <- pq[i, i]
      qk <- t(U) %*% Qi %*% U
      if (any(qk != 0)) { Q[[k]] <- qk; k <- k + 1 }
    }
  } else if (OPTION == "manual") {
    pq <- spm_inv(pC2)
    k <- 1
    for (i in seq_along(Qcells)) {
      diag_idx <- which(diag(Qcells[[i]]) != 0)
      jj <- which(seq_len(length(spm_vec(DCM$M$pE))) %in% q[diag_idx])
      Qi <- matrix(0, Np, Np)
      if (length(jj)) Qi[jj, jj] <- pq[jj, jj, drop = FALSE]
      qk <- t(U) %*% Qi %*% U
      if (any(qk != 0)) { Q[[k]] <- qk; k <- k + 1 }
    }
  } # 'none' -> Q empty

  Ng <- length(Q); Nw <- ncol(W); Nb <- Nw * Nx

  # Hyperpriors
  gE <- rep(M$hE %||% 0, Ng)
  gC <- diag(M$hC %||% (1/16), nrow = Ng, ncol = Ng)

  # X scaling (Xc)
  Xc <- tryCatch(diag(nrow(X) / colSums(X^2)), error = function(e) diag(1, Nx))

  # Prior expectations and covariances for second-level parameters
  Xe <- matrix(1, nrow = Nx, ncol = 1)   # column of ones (spm_speye(Nx,1))
  bE <- as.vector(kronecker(Xe, t(U) %*% bE0))
  bC <- kronecker(Xc, t(U) %*% bC0 %*% U)
  bP <- spm_inv(bC)
  gP <- spm_inv(gC)

  # Initial values
  b <- bE
  g <- gE

  list(
    # Inputs
    P = P, M = M, field = field,

    # Core sets
    Ns = Ns, X = X, W = W, U = U, q = q, Pstr = Pstr,
    pE = pE, pC = pC, qE = qE, qC = qC, iF = iF,

    # Sizes
    Nx = ncol(X), Nw = ncol(W), Ng = Ng, Np = length(q), Nb = Nb,

    # Priors/components
    bE0 = bE0, bC0 = bC0, pC2 = pC2, pQ = pQ, Q = Q,
    gE = gE, gC = gC, Xnames = (M$Xnames %||% lapply(seq_len(ncol(X)), \(i) sprintf("Covariate %d", i))),
    Xc = Xc,

    # Second-level priors (projected)
    bE = bE, bC = bC, bP = bP, gP = gP,

    # Initial parameter vectors
    b = b, g = g,

    # Control
    maxit = maxit,
    is_peb_of_pebs = is_peb_of_pebs
  )
}
# 
# 
# 
# 
# # ======================================================================
# # Part 2: spm_dcm_peb_run
# # - Runs VL using spm_dx on the prepared structure
# # - Returns PEB and (optionally) reduced first-level DCMs
# # Assumes your helpers exist: spm_dx, spm_log_evidence_reduce, etc.
# # ======================================================================
spm_dcm_peb_run <- function(prep, verbose = TRUE) {
  # Unpack
  P   <- prep$P; M <- prep$M
  Ns  <- prep$Ns; X <- prep$X; W <- prep$W; U <- prep$U
  q   <- prep$q; Pstr <- prep$Pstr
  pE  <- prep$pE; pC <- prep$pC; qE <- prep$qE; qC <- prep$qC; iF <- prep$iF
  Nx  <- prep$Nx; Nw <- prep$Nw; Ng <- prep$Ng; Nb <- prep$Nb
  bE0 <- prep$bE0; bC0 <- prep$bC0; pC2 <- prep$pC2; pQ <- prep$pQ; Q <- prep$Q
  gE  <- prep$gE; gC <- prep$gC; Xc <- prep$Xc; Xnames <- prep$Xnames
  bE  <- prep$bE; bC <- prep$bC; bP <- prep$bP; gP <- prep$gP
  b   <- prep$b; g <- prep$g
  maxit <- prep$maxit
  is_peb_of_pebs <- prep$is_peb_of_pebs

  # Joint prior precision over (b,g)
  ipC <- spm_cat(matrix(list(bP,NULL, 0, spm_inv(gC)), 2, 2, byrow = TRUE))

  # VL loop
  t_fs <- -4
  F0 <- -Inf
  tmp <- NULL

  for (niter in seq_len(maxit)) {
    # rP = floor + sum exp(g_i) Q_i
    if (Ng > 0) {
      rP <- pQ * exp(-8)
      for (i in seq_len(Ng)) rP <- rP + exp(g[i]) * Q[[i]]
    } else {
      rP <- pQ * exp(-8)
    }
    rC <- spm_inv(rP)

    # Gradients and curvature
    F <- 0
    dFdb  <- -bP %*% (b - bE)
    dFdbb <- -bP
    dFdg  <- -spm_inv(gC) %*% (g - gE)
    dFdgg <- -spm_inv(gC)
    dFdbg <- matrix(0, nrow = Nb, ncol = Ng)

    for (i in seq_len(Ns)) {
      Xi <- kronecker(matrix(X[i, ], nrow = 1), W)
      rE_i <- as.vector(Xi %*% b)

      red <- spm_log_evidence_reduce(qE[[i]], qC[[i]], pE[[i]], pC[[i]], rE_i, rC)
      Fi <- red$F; sE <- red$sE; sC <- red$sC
      F <- F + Fi + iF[i]
      dE <- sE - rE_i

      dFdb  <- dFdb  + t(Xi) %*% (rP %*% dE)
      dFdbb <- dFdbb + t(Xi) %*% ((rP %*% sC %*% rP - rP)) %*% Xi

      if (Ng > 0) {
        for (j in seq_len(Ng)) {
          dFdgj <- exp(g[j]) * ( (sum(diag((rC - sC) %*% Q[[j]])) - t(dE) %*% Q[[j]] %*% dE) / 2 )
          dFdg[j] <- dFdg[j] + as.numeric(dFdgj)
          dFdgg[j, j] <- dFdgg[j, j] + as.numeric(dFdgj)

          dFdbgj <- exp(g[j]) * t( (Xi - sC %*% rP %*% Xi) ) %*% Q[[j]] %*% dE
          dFdbg[, j] <- dFdbg[, j] + as.vector(dFdbgj)

          for (k in seq_len(Ng)) {
            dFdggj <- exp(g[j] + g[k]) * (
              ( sum(diag((rC %*% Q[[k]] %*% rC - sC %*% Q[[k]] %*% sC) %*% Q[[j]])) / 2 )
              - t(dE) %*% Q[[k]] %*% sC %*% Q[[j]] %*% dE
            )
            dFdgg[j, k] <- dFdgg[j, k] - as.numeric(dFdggj)
          }
        }
      }
    }

    # Block matrices and covariances
    dFdp  <- rbind(dFdb, dFdg)
    dFdpp <- rbind(cbind(dFdbb, dFdbg), cbind(t(dFdbg), dFdgg))
    Cp    <- spm_inv(-dFdpp)
    Cb    <- spm_inv(-dFdbb)
    Cg    <- spm_inv(-dFdgg)

    # Complexity of second level
    Fb <- as.numeric(t(b) %*% bP %*% b)
    Fg <- as.numeric(t(g) %*% spm_inv(gC) %*% g)
    Fc <- Fb/2 + Fg/2 - spm_logdet(as.matrix(ipC) %*% as.matrix(Cp))/2
    F  <- F - Fc

    if (niter == 1) F0 <- F
    if (F >= F0) {
      dF <- F - F0; F0 <- F
      tmp <- list(b = b, g = g, F0 = F0, dFdb = dFdb, dFdbb = dFdbb, dFdg = dFdg, dFdgg = dFdgg)
      t_fs <- min(t_fs + 0.25, 2)
    } else {
      t_fs <- max(t_fs - 1, -4)
      b    <- tmp$b; g <- tmp$g; F0 <- tmp$F0
      dFdb <- tmp$dFdb; dFdbb <- tmp$dFdbb; dFdg <- tmp$dFdg; dFdgg <- tmp$dFdgg
    }

    # Fisher scoring step (Gauss-Newton update with spm_dx)
    dp <- spm_dx(as.matrix(dFdpp), as.matrix(dFdp), list(t_fs))
    if (sum(abs(dp)) >= 8) {
      dFdpp_blk <- rbind(cbind(dFdbb, matrix(0, nrow(dFdbb), ncol(dFdgg))),
                         cbind(matrix(0, nrow(dFdgg), ncol(dFdbb)), dFdgg))
      dp <- spm_dx(dFdpp_blk, dFdp, list(t_fs))
    }
    db <- dp[seq_len(length(b))]
    dg <- dp[(length(b) + 1):length(dp)]
    b  <- b + db
    g  <- g + tanh(dg)

    if (isTRUE(verbose)) {
      cat(sprintf("VL Iter %3d: F=%8.3f  dF=%8.4f  [t=%+.2f]\n", niter, F, ifelse(exists("dF"), dF, 0), t_fs))
    }
    if (niter > 4 && (t_fs <= -4 || (exists("dF") && dF < 1e-4))) break
  }

  # Assemble PEB
  Ub <- kronecker(diag(Nx), U)
  PEB <- list()
  PEB$Snames <- character(Ns)
  PEB$Pnames <- unlist(Pstr)
  PEB$Pind   <- q
  PEB$Xnames <- unlist(Xnames)
  if (is_peb_of_pebs) {
    PEB$Pind0 <- do.call(rbind, replicate(ncol(M$X), P[[1]]$Pind0, simplify = FALSE))
  } else {
    PEB$Pind0 <- PEB$Pind
  }

  PEB$M <- list()
  PEB$M$X  <- X
  PEB$M$W  <- W
  PEB$M$U  <- U
  PEB$M$pE <- kronecker(matrix(1, nrow = Nx, ncol = 1), M$bE %||% bE0)
  PEB$M$pC <- kronecker(Xc, M$bC %||% bC0)
  PEB$M$hE <- gE
  PEB$M$hC <- gC

  PEB$Ep <- U %*% matrix(b, nrow = Nw, ncol = Nx)
  PEB$Eh <- g
  PEB$Ch <- spm_inv(-dFdgg)
  PEB$Cp <- Ub %*% spm_inv(-dFdbb) %*% t(Ub)
  PEB$Ce <- U %*% spm_inv(pQ * exp(-8) + Reduce(`+`, lapply(seq_along(Q), function(i) exp(g[i]) * Q[[i]]))) %*% t(U)
  PEB$F  <- F
  PEB$M$Q <- lapply(Q, function(Qi) U %*% Qi %*% t(U))

  # Optional: update 1st-level DCMs with empirical priors
  P_out <- vector("list", Ns)
  for (i in seq_len(Ns)) {
    D <- P[[i]]
    RP <- spm_inv(if (is.list(D$M$pC)) diag(spm_vec(D$M$pC)) else as.matrix(D$M$pC))
    RP[q, q] <- U %*% (pQ * exp(-8) + Reduce(`+`, lapply(seq_along(Q), function(k) exp(g[k]) * Q[[k]]))) %*% t(U)
    RC <- spm_inv(RP)

    RE <- spm_vec(D$M$pE)
    # Use first-row covariate for reconstruction (matches MATLAB’s simple path)
    RE[q] <- U %*% (matrix(b, nrow = Nw, ncol = Nx) %*% matrix(X[1, ], ncol = 1))
    RE <- spm_unvec(RE, D$M$pE)

    red <- spm_log_evidence_reduce(spm_vec(D$Ep), D$Cp, spm_vec(D$M$pE), if (is.list(D$M$pC)) diag(spm_vec(D$M$pC)) else as.matrix(D$M$pC), spm_vec(RE), RC)

    D$M$pE <- RE
    D$M$pC <- RC
    D$Ep   <- spm_unvec(red$sE, D$Ep)
    D$Cp   <- red$sC
    # share of complexity (simple split as in MATLAB code)
    Fb <- as.numeric(t(b) %*% bP %*% b)
    Fg <- as.numeric(t(g) %*% spm_inv(gC) %*% g)
    Fc <- Fb/2 + Fg/2
    D$F    <- as.numeric(red$F + iF[i] - Fc/Ns)
    P_out[[i]] <- D
    #PEB$Snames[i] <- tryCatch(D$name, error = function(e) sprintf("Subject %d", i))
    
    nm <- NA_character_
    if (is.character(P[[i]]) && length(P[[i]]) == 1) {
      nm <- sub("\\.[Rr][Dd][Ss]$|\\.[Mm][Aa][Tt]$", "", basename(P[[i]]))
    } else if (!is.null(D$name) && length(D$name) >= 1 && nzchar(as.character(D$name[1]))) {
      nm <- as.character(D$name[1])
    } else if (!is.null(D$Snames) && length(D$Snames) >= 1 && nzchar(as.character(D$Snames[1]))) {
      nm <- paste0("PEB_", i, ":", as.character(D$Snames[1]))
    } else {
      nm <- sprintf("Subject %d", i)
    }
    PEB$Snames[i] <- nm
    
  }

  list(PEB = PEB, P = P_out)
}

# === Third-level design ===
build_third_level_design <- function(n_subj, covariates = NULL, Xnames = NULL) {
  if (is.null(covariates)) {
    X <- matrix(1, nrow = n_subj, ncol = 1)
    colnames(X) <- "Intercept"
  } else {
    stopifnot(is.data.frame(covariates), nrow(covariates) == n_subj)
    # model.matrix includes intercept by default
    X <- model.matrix(~ ., data = covariates)
  }
  if (!is.null(Xnames)) colnames(X) <- Xnames
  X
}

# === Third-level PEB-of-PEBs orchestrator ======================================
# Loads subject PEBs (.rds or .mat), builds the third-level design, runs the
# PEB-of-PEBs, saves the group PEB, and returns a tidy summary.
#   covariates : NULL -> intercept-only; else a data.frame aligned with subjects
#   Xnames     : optional custom design column names
#   field      : "all" is typical for PEB inputs
peb_of_pebs <- function(peb_dir,
                        subjects   = NULL,
                        covariates = NULL,
                        Xnames     = NULL,
                        save_path  = file.path(peb_dir, "PEB_GROUP.rds"),
                        M     = list(Q = "single", alpha = 1, beta = 16, maxit = 64),
                        field = "all",
                        verbose = TRUE) {

  # 1) Load subject PEBs (.rds or .mat, normalized to native shape)
  L      <- load_pebs_generic(peb_dir, subjects = subjects)
  pebs   <- unname(L$pebs)
  subj   <- L$subjects
  n_subj <- length(pebs)
  if (verbose) message(sprintf("Loaded %d subject PEBs.", n_subj))

  # 2) Third-level design
  X <- build_third_level_design(n_subj, covariates = covariates, Xnames = Xnames)
  if (verbose) {
    message("Third-level design X (rows=subjects, cols=covariates):")
    capture.output(print(X)) |> paste(collapse = "\n") |> message()
  }

  # 3) Third-level model options
  M2 <- modifyList(list(X = X), M)

  # 4) Prepare & run the third-level PEB-of-PEBs
  prep3 <- spm_dcm_peb_prepare(P = pebs, M = M2, field = field)
  run3  <- spm_dcm_peb_run(prep3, verbose = verbose)

  # 5) Save group PEB
  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(run3$PEB, save_path)
  if (verbose) message(sprintf("Saved group PEB: %s", save_path))

  # 6) Tidy summary
  Ep <- as.numeric(run3$PEB$Ep)
  Eh <- if (length(run3$PEB$Eh)) run3$PEB$Eh else NA_real_
  F  <- run3$PEB$F
  Pnames <- run3$PEB$Pnames
  if (!is.null(dim(Pnames))) Pnames <- as.vector(Pnames)

  summary_df <- data.frame(
    param = if (length(Pnames)) Pnames else paste0("P", seq_along(Ep)),
    Ep    = Ep,
    Eh    = ifelse(length(Eh), Eh[1], NA_real_),
    F     = F,
    stringsAsFactors = FALSE
  )

  list(
    group_PEB = run3$PEB,
    save_path = save_path,
    subjects  = subj,
    X         = X,
    summary   = summary_df
  )
}
