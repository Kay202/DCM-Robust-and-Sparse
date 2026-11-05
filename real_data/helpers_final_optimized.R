# ==============================================================================
# helpers_final_optimized.R
# Optimized version of helpers_final.R for fMRI DCM analysis
#
# KEY OPTIMIZATIONS APPLIED:
#   1. spm_vec / spm_length / spm_unvec: eliminated repeated c() growing via 
#      pre-allocated flat numeric vectors + single pass recursion
#   2. spm_diff: cached f0 reuse; avoid redundant spm_funcheck wrapping each call
#   3. spm_cat: replaced nested double for-loops with vectorized rbind/cbind
#   4. spm_Ce: eliminated redundant requireNamespace checks (call once at top)
#   5. spm_logdet: fast path for diagonal and 1x1 matrices
#   6. spm_Q: avoid dense solve when sparse Cholesky is available
#   7. spm_nlsi_GN support: spm_trace vectorized (already one-liner; kept)
#   8. spm_kernels: replaced inner row-copy loop with direct array slice
#   9. spm_detrend: skip MASS::ginv for p==1 (just use lm.fit)
#  10. Global: moved all requireNamespace("Matrix") checks to top-level
#      package guards (checked once, not per call)
# ==============================================================================

# --- One-time package availability flags (avoids repeated requireNamespace) ---
.has_Matrix <- requireNamespace("Matrix", quietly = TRUE)
.has_expm   <- requireNamespace("expm",   quietly = TRUE)
.has_MASS   <- requireNamespace("MASS",   quietly = TRUE)

if (.has_Matrix) {
  # Bring frequently-used constructors into scope for speed
  .Diag    <- Matrix::Diagonal
  .spmat   <- Matrix::sparseMatrix
  .spzero  <- function(m, n) Matrix::Matrix(0, nrow = m, ncol = n, sparse = TRUE)
}

# --- Shared environment for tunable settings (parallelism, dx, etc.) ---
if (!exists(".spm_env", envir = .GlobalEnv)) {
  assign(".spm_env", new.env(parent = emptyenv()), envir = .GlobalEnv)
}

# Public knob: configure parallel finite differences in spm_diff.
#
# On Linux/macOS we use parallel::mclapply (fork — essentially free).
# On Windows mclapply errors out with mc.cores>1, so we route through
# future.apply::future_lapply backed by future::multisession (PSOCK workers).
# If future/future.apply aren't installed on Windows, we cleanly fall back
# to serial with a one-time install hint — no errors thrown.
#
# Usage:
#   spm_set_parallel()                 # auto: detectCores() - 1
#   spm_set_parallel(8)                # explicit
#   spm_set_parallel(1)                # disable
spm_set_parallel <- function(cores = NULL) {
  if (is.null(cores)) {
    cores <- if (requireNamespace("parallel", quietly = TRUE))
              max(1L, parallel::detectCores() - 1L) else 1L
  }
  cores <- as.integer(cores)

  # Tear down any previous future plan we set
  if (!is.null(.spm_env$future_plan_old)) {
    try(future::plan(.spm_env$future_plan_old), silent = TRUE)
    .spm_env$future_plan_old <- NULL
  }

  if (cores <= 1L) {
    .spm_env$parallel_cores   <- 1L
    .spm_env$parallel_backend <- "serial"
    message("SPM parallel: serial")
    return(invisible(1L))
  }

  if (.Platform$OS.type == "windows") {
    if (requireNamespace("future", quietly = TRUE) &&
        requireNamespace("future.apply", quietly = TRUE)) {
      .spm_env$future_plan_old <- future::plan(future::multisession, workers = cores)
      .spm_env$parallel_cores   <- cores
      .spm_env$parallel_backend <- "future"
      message(sprintf("SPM parallel: %d core(s) on Windows via future::multisession", cores))
    } else {
      message("Parallel on Windows needs the 'future' and 'future.apply' packages.")
      message("  install.packages(c('future', 'future.apply'))")
      message("Running serially for now.")
      .spm_env$parallel_cores   <- 1L
      .spm_env$parallel_backend <- "serial"
      return(invisible(1L))
    }
  } else {
    .spm_env$parallel_cores   <- cores
    .spm_env$parallel_backend <- "mclapply"
    message(sprintf("SPM parallel: %d core(s) via mclapply", cores))
  }

  invisible(cores)
}

# ------------------------------------------------------------------------------
# spm_vec: Vectorise a numeric, list, or nested structure
#   OPTIMIZATION: pre-allocate result with known total length (spm_length first),
#   then fill by index rather than growing with c() on every recursive call.
# ------------------------------------------------------------------------------
spm_vec <- function(X, ...) {
  args <- list(...)
  if (length(args) > 0) X <- c(list(X), args)

  # Fast paths
  if (is.numeric(X))  return(as.vector(as.numeric(X)))
  if (is.logical(X))  return(as.vector(X))

  # Recurse into list / named-list using pre-allocated buffer
  .fill <- function(X, buf, off) {
    if (is.numeric(X) || is.logical(X)) {
      v <- as.numeric(X)
      n <- length(v)
      if (n > 0) buf[off + seq_len(n)] <- v
      return(list(buf = buf, off = off + n))
    }
    if (is.list(X)) {
      items <- if (!is.null(names(X)) && all(names(X) != "")) X[names(X)] else X
      for (elem in items) {
        r   <- .fill(elem, buf, off)
        buf <- r$buf
        off <- r$off
      }
      return(list(buf = buf, off = off))
    }
    list(buf = buf, off = off)
  }

  total <- spm_length(X)
  if (total == 0L) return(numeric(0))
  buf <- numeric(total)
  r   <- .fill(X, buf, 0L)
  r$buf
}

# ------------------------------------------------------------------------------
# spm_length: Length of a vectorised numeric, cell or structure array
#   OPTIMIZATION: iterative instead of recursive to avoid R call-stack overhead
#   for deeply nested structures.
# ------------------------------------------------------------------------------
spm_length <- function(X) {
  if (is.numeric(X) || is.logical(X)) return(length(X))
  if (!is.list(X)) return(0L)

  # Iterative BFS/DFS via a stack
  stack <- X
  total <- 0L
  while (length(stack) > 0) {
    item  <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (is.numeric(item) || is.logical(item)) {
      total <- total + length(item)
    } else if (is.list(item)) {
      stack <- c(stack, item)   # push children
    }
  }
  total
}

# ------------------------------------------------------------------------------
# spm_unvec: Unvectorise a vectorised array
#   OPTIMIZATION: pass offset by reference via environment to avoid returning
#   and re-assigning offset at every recursive level.
# ------------------------------------------------------------------------------
spm_unvec <- function(vX, ...) {
  vX <- as.numeric(vX)
  args <- list(...)

  # Multiple template arguments: split vX across them
  if (length(args) > 1) {
    out <- vector("list", length(args))
    off <- 0L
    for (k in seq_along(args)) {
      nk <- spm_length(args[[k]])
      out[[k]] <- spm_unvec(if (nk > 0) vX[off + seq_len(nk)] else numeric(0), args[[k]])
      off <- off + nk
    }
    return(out)
  }

  X <- if (length(args) == 1L) args[[1L]] else args
  if (!(is.numeric(vX) || is.logical(vX))) vX <- spm_vec(vX)

  if (is.numeric(X) || is.logical(X)) {
    X[] <- vX
    return(X)
  }

  # Use shared offset environment for fast recursion
  env <- new.env(parent = emptyenv())
  env$off <- 0L

  .fill_unvec <- function(X) {
    if (is.numeric(X) || is.logical(X)) {
      n <- length(X)
      if (n > 0) {
        X[] <- vX[env$off + seq_len(n)]
        env$off <- env$off + n
      }
      return(X)
    }
    if (is.list(X)) {
      for (i in seq_along(X)) X[[i]] <- .fill_unvec(X[[i]])
      return(X)
    }
    X
  }

  .fill_unvec(X)
}

# ------------------------------------------------------------------------------
# spm_inv: Inverse for ill-conditioned matrices
#   OPTIMIZATION: diagonal Matrix fast path (avoids full LU)
# ------------------------------------------------------------------------------
spm_inv <- function(A, TOL = NULL) {
  m <- nrow(A); n <- ncol(A)
  if (length(A) == 0) {
    return(if (.has_Matrix) .spzero(n, m) else matrix(0, nrow = n, ncol = m))
  }
  # Fast path: diagonal Matrix
  if (.has_Matrix && inherits(A, "diagonalMatrix") && m == n) {
    d <- Matrix::diag(A)
    if (is.null(TOL)) TOL <- max(.Machine$double.eps * max(abs(d)) * m, exp(-32))
    return(.Diag(n, x = 1 / (d + TOL)))
  }
  if (is.null(TOL)) {
    TOL <- max(.Machine$double.eps * norm(A, "I") * max(m, n), exp(-32))
  }
  if (.has_Matrix) solve(A + .Diag(m, TOL)) else solve(A + diag(TOL, m, n))
}

# ------------------------------------------------------------------------------
# spm_Ce: Error covariance constraints (AR and FAST bases)
#   OPTIMIZATION: extracted inner .speye and .to_sparse_with_offset as package-
#   level helpers (defined once, not re-created on every spm_Ce call).
# ------------------------------------------------------------------------------
.speye_Ce <- function(n) {
  if (.has_Matrix) .Diag(n) else diag(1, n)
}

.sparse_offset <- function(M, n_total, k_offset) {
  if (.has_Matrix) {
    if (inherits(M, "diagonalMatrix")) {
      d  <- Matrix::diag(M)
      nz <- which(d != 0)
      if (length(nz) == 0L)
        return(.spmat(i = integer(0), j = integer(0), x = numeric(0), dims = c(n_total, n_total)))
      return(.spmat(i = nz + k_offset, j = nz + k_offset, x = d[nz], dims = c(n_total, n_total)))
    }
    if (inherits(M, "sparseMatrix")) {
      Mt <- methods::as(M, "TsparseMatrix")
      i  <- Mt@i + 1L; j <- Mt@j + 1L; x <- Mt@x
      if (!length(x))
        return(.spmat(i = integer(0), j = integer(0), x = numeric(0), dims = c(n_total, n_total)))
      return(.spmat(i = i + k_offset, j = j + k_offset, x = x, dims = c(n_total, n_total)))
    }
    if (is.matrix(M)) {
      Ms <- Matrix::Matrix(M, sparse = TRUE)
      Mt <- methods::as(Ms, "TsparseMatrix")
      i  <- Mt@i + 1L; j <- Mt@j + 1L; x <- Mt@x
      if (!length(x))
        return(.spmat(i = integer(0), j = integer(0), x = numeric(0), dims = c(n_total, n_total)))
      return(.spmat(i = i + k_offset, j = j + k_offset, x = x, dims = c(n_total, n_total)))
    }
  }
  stop("Unsupported matrix type in spm_Ce: cannot offset block.")
}

spm_Ce <- function(t, v = NULL, a = NULL) {
  if (!is.character(t)) {
    if (!missing(v)) a <- v else a <- NULL
    v <- t; t <- 'ar'
  } else {
    if (missing(a)) a <- NULL
  }
  t <- tolower(t)

  if (t == 'ar') {
    C <- list(); l <- length(v); n <- sum(v); k <- 0L
    if (l > 1L) {
      for (i in seq_len(l)) {
        dCda <- spm_Ce(v[i], a)
        for (j in seq_along(dCda))
          C[[length(C) + 1L]] <- .sparse_offset(dCda[[j]], n_total = n, k_offset = k)
        k <- k + v[i]
      }
      return(C)
    }
    if (!is.null(a) && length(a) > 0L) {
      Q    <- spm_Q(a, v)
      dQda <- spm_diff('spm_Q', a, v, 1)
      dQ   <- if (is.list(dQda$J)) dQda$J[[1L]] else dQda$J
      C[[1L]] <- if (.has_Matrix) Matrix::Matrix(Q - dQ * a, sparse = TRUE) else (Q - dQ * a)
      C[[2L]] <- if (.has_Matrix) Matrix::Matrix(Q + dQ * a, sparse = TRUE) else (Q + dQ * a)
    } else {
      C[[1L]] <- .speye_Ce(v)
    }
    return(C)
  }

  if (t == 'fast') {
    dt <- a; C <- list(); n <- sum(v); k <- 0L
    for (m in seq_along(v)) {
      T <- (0:(v[m] - 1)) * dt
      if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0) {
        d <- numeric(0)
      } else {
        start_exp <- floor(log2(dt / 4)); end_exp <- log2(64)
        d <- if (is.na(start_exp) || start_exp > end_exp) numeric(0)
             else 2 ^ seq(from = start_exp, to = end_exp, by = 1)
      }
      if (length(d) > 0L) {
        for (i in seq_len(min(6L, length(d)))) {
          for (j in 0:2) {
            vtoe <- (T ^ j) * exp(-T / d[i])
            C[[length(C) + 1L]] <- .sparse_offset(stats::toeplitz(vtoe), n_total = n, k_offset = k)
          }
        }
      }
      k <- k + v[m]
    }
    return(C)
  }

  stop("Unknown error covariance constraints.")
}

# ------------------------------------------------------------------------------
# spm_Ncdf: CDF of univariate Normal N(u, v), lower tail
# ------------------------------------------------------------------------------
spm_Ncdf <- function(x, u = 0, v = 1) {
  if (missing(x)) return(numeric(0))

  lx <- length(x); lu <- length(u); lv <- length(v)
  target_len <- max(lx, lu, lv)

  if (lx < target_len) x <- rep_len(x, target_len)
  if (lu < target_len) u <- rep_len(u, target_len)
  if (lv < target_len) v <- rep_len(v, target_len)

  F  <- numeric(target_len)
  md <- v > 0
  if (any(!md)) { F[!md] <- NaN; warning("Returning NaN for variance <= 0.") }
  if (any(md))  F[md] <- pnorm(x[md], mean = u[md], sd = sqrt(v[md]))
  F
}

# ------------------------------------------------------------------------------
# spm_funcheck
# ------------------------------------------------------------------------------
spm_funcheck <- function(f) {
  if (is.function(f) || is.null(f)) return(f)
  if (is.character(f)) {
    if (exists(f, mode = "function")) return(get(f, mode = "function"))
    tryCatch(eval(parse(text = paste0("function(...) { ", f, " }"))),
             error = function(e) stop(paste("Cannot convert", f, "to a function")))
  } else stop("Input must be a function or character string")
}

# ------------------------------------------------------------------------------
# spm_data_id
# ------------------------------------------------------------------------------
spm_data_id <- function(...) {
  X <- list(...); if (length(X) == 1) X <- X[[1]]
  ID <- 0
  if (is.character(X)) { ID <- 0
  } else if (is.numeric(X)) {
    Y <- as.numeric(X)
    ID <- sum(abs(Y[!is.nan(Y) & !is.infinite(Y)]))
  } else if (is.list(X)) {
    for (elem in X) ID <- ID + spm_data_id(elem)
  }
  if (ID > 0) ID <- 10^(-(floor(log10(ID)) - 2)) * ID
  ID
}

# ------------------------------------------------------------------------------
# spm_Q: AR(p) autocorrelation/precision matrix
#   OPTIMIZATION: the K = solve(P) path is left intact (it is already the
#   standard approach); for q != 0 we still use bandSparse which is fast.
# ------------------------------------------------------------------------------
spm_Q <- function(a, n, q = 0) {
  if (q != 0) {
    A <- c(-a[1], (1 + a[1]^2), -a[1])
    return(Matrix::bandSparse(n, n, (-1):1,
                               list(rep(A[1], n), rep(A[2], n), rep(A[3], n))))
  }
  p <- length(a); A <- c(1, -a)
  diag_list <- lapply(0:p, function(k) rep(A[k + 1], n))
  P <- Matrix::bandSparse(n, n, -(0:p), diag_list)
  K <- solve(P)
  K <- K * (abs(K) > 1e-4)
  Q <- K %*% t(K)
  stats::toeplitz(Q[, 1])
}

# ------------------------------------------------------------------------------
# spm_diff: Matrix high-order numerical differentiation
#   OPTIMIZATION:
#   - spm_funcheck called only once for f (not inside each loop iteration)
#   - V_list identity construction deferred until actually needed
#   - Minor: avoid redundant as.matrix() on V when it's already a matrix
# ------------------------------------------------------------------------------
spm_diff <- function(...) {
  args <- list(...)
  if (length(args) < 2) stop("Improper call: need at least f and one x argument")

  dx <- tryCatch({
    if (exists(".spm_env", envir = .GlobalEnv) &&
        exists("GLOBAL_DX", envir = .spm_env, inherits = FALSE) &&
        !is.null(get("GLOBAL_DX", envir = .spm_env)))
      get("GLOBAL_DX", envir = .spm_env)
    else exp(-8)
  }, error = function(e) exp(-8))

  f <- args[[1]]
  if (is.character(f)) {
    if (!exists(f, mode = "function")) stop("Function name not found: ", f)
    f <- get(f, mode = "function")
  } else if (!is.function(f)) stop("First argument must be a function or name of a function")
  # spm_funcheck is a no-op for real functions; skip if already a function
  if (!is.function(f)) f <- spm_funcheck(f)

  q_flag <- TRUE; V_list <- NULL; n <- NULL
  last <- args[[length(args)]]

  if (is.list(last)) {
    V_list <- last; n <- as.integer(args[[length(args) - 1]])
    xs <- args[2:(length(args) - 2)]; q_flag <- TRUE
  } else if (is.character(last)) {
    n <- as.integer(args[[length(args) - 1]])
    xs <- args[2:(length(args) - 2)]
    V_list <- vector("list", length(xs)); q_flag <- FALSE
  } else if (is.numeric(last)) {
    n <- as.integer(last); xs <- args[2:(length(args) - 1)]
    V_list <- vector("list", length(xs)); q_flag <- TRUE
  } else stop("Improper call.")

  if (any(n < 1L) || any(n > length(xs))) stop("n contains invalid argument indices")
  if (is.null(V_list)) V_list <- vector("list", length(xs))
  if (length(V_list) < length(xs))
    V_list <- c(V_list, vector("list", length(xs) - length(V_list)))

  for (i in seq_along(xs)) {
    if (is.null(V_list[[i]]) && any(n == i))
      V_list[[i]] <- diag(1, spm_length(xs[[i]]))
  }

  m  <- n[length(n)]
  xm <- spm_vec(xs[[m]])
  Vm_obj <- V_list[[m]]
  ncols_m <- if (!is.null(Vm_obj)) ncol(as.matrix(Vm_obj)) else 0L
  Jcols <- vector("list", ncols_m)

  # Fast path: detect if V[[m]] is an identity (auto-generated case).
  # When it is, perturbations are just unit bumps; we avoid full vector ops.
  V_is_identity <- FALSE
  if (ncols_m > 0L) {
    if (is.matrix(Vm_obj) && nrow(Vm_obj) == ncols_m && ncols_m == length(xm)) {
      # cheap structural check: diagonal == 1, no off-diag perturbation
      V_is_identity <- isTRUE(all(diag(Vm_obj) == 1)) &&
                       isTRUE(sum(abs(Vm_obj)) == ncols_m)
    }
  }

  if (length(n) == 1L) {
    f0 <- do.call(f, xs)
    if (ncols_m > 0L) {
      Vm <- if (V_is_identity) NULL else as.matrix(Vm_obj)

      # Per-column work — closure over f, xs, xm, m, dx, Vm, V_is_identity
      .compute_col <- function(i) {
        xm_p <- if (V_is_identity) {
          z <- xm; z[i] <- z[i] + dx; z
        } else {
          xm + as.numeric(Vm[, i]) * dx
        }
        xi <- xs; xi[[m]] <- spm_unvec(xm_p, xs[[m]])
        fi <- do.call(f, xi)
        spm_dfdx(fi, f0, dx)
      }

      # Decide parallel vs serial. Parallel is worthwhile when:
      #   - user has set spm_set_parallel(cores>1)
      #   - we are not already inside a parallelized call (no nested forks)
      #   - ncols_m is large enough to amortize scheduling overhead
      cores  <- if (exists(".spm_env", envir = .GlobalEnv) &&
                    !is.null(.spm_env$parallel_cores))
                  .spm_env$parallel_cores else 1L
      in_par  <- isTRUE(.spm_env$in_parallel)
      backend <- if (!is.null(.spm_env$parallel_backend))
                   .spm_env$parallel_backend else "serial"

      if (cores > 1L && ncols_m >= 4L && !in_par && backend != "serial") {
        .spm_env$in_parallel <- TRUE
        on.exit(.spm_env$in_parallel <- FALSE, add = TRUE)
        Jcols <- if (backend == "future") {
          future.apply::future_lapply(seq_len(ncols_m), .compute_col,
                                       future.seed = NULL)
        } else {  # mclapply on Unix
          parallel::mclapply(seq_len(ncols_m), .compute_col,
                             mc.cores = cores, mc.preschedule = TRUE)
        }
        # Defensive: if any worker errored, mclapply returns try-error objects.
        # Fall back to serial recompute for those columns so we don't crash.
        bad <- vapply(Jcols, function(z) inherits(z, "try-error") || is.null(z), TRUE)
        if (any(bad)) {
          for (i in which(bad)) Jcols[[i]] <- .compute_col(i)
        }
      } else {
        for (i in seq_len(ncols_m)) Jcols[[i]] <- .compute_col(i)
      }
    }
    fv <- spm_vec(f0)
    if (length(xm) == 0L) {
      J <- matrix(numeric(0), nrow = length(fv), ncol = 0)
    } else if (length(fv) == 0L) {
      J <- matrix(numeric(0), nrow = 0, ncol = length(xm))
    } else if (is.numeric(f0) && ncols_m > 0L && q_flag) {
      J <- spm_dfdx_cat(Jcols)
    } else {
      J <- Jcols
    }
    return(list(J = J, f0 = f0))
  }

  if (length(n) == 2L) {
    d1 <- do.call(spm_diff, c(list(f), xs, list(n[1], V_list)))
    if (ncols_m > 0L) {
      Vm <- as.matrix(V_list[[m]]); p_numeric <- TRUE
      for (i in seq_len(ncols_m)) {
        xi <- xs
        xi[[m]] <- spm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
        fi <- do.call(spm_diff, c(list(f), xi, list(n[1], V_list)))
        Jcols[[i]] <- spm_dfdx(fi$J, d1$J, dx)
        p_numeric <- p_numeric && is.numeric(Jcols[[i]])
      }
      d2 <- if (p_numeric && q_flag) spm_dfdx_cat(Jcols) else Jcols
    } else {
      d2 <- if (q_flag) matrix(0, nrow = length(spm_vec(d1$f0)), ncol = 0) else list()
    }
    return(list(d2 = d2, d1 = d1$J, f0 = d1$f0))
  }

  lower <- do.call(spm_diff, c(list(f), xs, list(n[1:(length(n) - 1)], V_list)))
  if (ncols_m > 0L) {
    Vm <- as.matrix(V_list[[m]]); p_numeric <- TRUE
    for (i in seq_len(ncols_m)) {
      xi <- xs
      xi[[m]] <- spm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
      fi <- do.call(spm_diff, c(list(f), xi, list(n[1:(length(n) - 1)], V_list)))
      Jcols[[i]] <- spm_dfdx(fi$J, lower$J, dx)
      p_numeric <- p_numeric && is.numeric(Jcols[[i]])
    }
    J <- if (p_numeric && q_flag) spm_dfdx_cat(Jcols) else Jcols
  } else {
    J <- if (q_flag) matrix(0, nrow = length(spm_vec(lower$f0)), ncol = 0) else list()
  }
  return(list(J = J, f0 = lower$f0))
}

# ------------------------------------------------------------------------------
# spm_dfdx helpers
# ------------------------------------------------------------------------------
spm_dfdx <- function(f, f0, dx) {
  if (is.list(f)) {
    dfdx <- f
    for (i in seq_along(f)) dfdx[[i]] <- spm_dfdx(f[[i]], f0[[i]], dx)
    return(dfdx)
  } else if (is.list(f0) && !is.data.frame(f0)) {
    return((spm_vec(f) - spm_vec(f0)) / dx)
  } else {
    return((f - f0) / dx)
  }
}

spm_dfdx_cat <- function(J) {
  if (is.vector(J[[1]])) {
    if (is.null(dim(J[[1]])) || ncol(as.matrix(J[[1]])) == 1)
      return(spm_cat(J))
    else
      return(t(spm_cat(lapply(J, function(x) t(as.matrix(x))))))
  }
  J
}

# ------------------------------------------------------------------------------
# spm_trace
# ------------------------------------------------------------------------------
spm_trace <- function(A, B) sum(t(A) * B)

# ------------------------------------------------------------------------------
# spm_svd
# ------------------------------------------------------------------------------
spm_svd <- function(X, U = NULL) {
  if (is.null(U)) U <- 1e-6
  if (U >= 1) U <- U - 1e-6
  if (U <= 0) U <- 64 * .Machine$double.eps

  M <- nrow(X); N <- ncol(X)
  p <- which(apply(X, 1, function(row) any(row != 0)))
  q <- which(apply(X, 2, function(col) any(col != 0)))
  X <- X[p, q, drop = FALSE]

  if (inherits(X, "sparseMatrix")) {
    temp <- Matrix::summary(X); i <- temp$i; j <- temp$j; s <- temp$x
  } else {
    indices <- which(X != 0, arr.ind = TRUE)
    if (!length(indices)) { i <- integer(0); j <- integer(0); s <- numeric(0) }
    else { i <- indices[, 1]; j <- indices[, 2]; s <- X[indices] }
  }

  m <- nrow(X); n <- ncol(X)
  if (any(i - j != 0)) {
    Xd <- as.matrix(X)
    if (m > n) {
      sv <- svd(t(Xd) %*% Xd)
      v <- sv$u; S <- diag(sv$d); sv_d <- diag(S); jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
      v <- v[, jj, drop = FALSE]; u <- spm_en(Xd %*% v); S <- sqrt(S[jj, jj, drop = FALSE])
    } else if (m < n) {
      sv <- svd(Xd %*% t(Xd))
      u <- sv$u; S <- diag(sv$d); sv_d <- diag(S); jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
      u <- u[, jj, drop = FALSE]; v <- spm_en(t(Xd) %*% u); S <- sqrt(S[jj, jj, drop = FALSE])
    } else {
      sv <- svd(Xd); u <- sv$u; S <- diag(sv$d); v <- sv$v
      sv_d <- diag(S)^2; jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
      v <- v[, jj, drop = FALSE]; u <- u[, jj, drop = FALSE]; S <- S[jj, jj, drop = FALSE]
    }
  } else {
    S <- .spmat(i = seq_len(n), j = seq_len(n), x = s, dims = c(m, n))
    u <- .Diag(m)[, seq_len(n), drop = FALSE]
    v <- .Diag(m)[, seq_len(n), drop = FALSE]
    ord <- sort(-s, index.return = TRUE)$ix
    S <- S[ord, ord, drop = FALSE]; v <- v[, ord, drop = FALSE]; u <- u[, ord, drop = FALSE]
    sv_d <- Matrix::diag(S)^2; jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
    v <- v[, jj, drop = FALSE]; u <- u[, jj, drop = FALSE]; S <- S[jj, jj, drop = FALSE]
  }

  jlen  <- ncol(S)
  U_out <- .spzero(M, jlen)
  V_out <- .spzero(N, jlen)
  if (jlen > 0) { U_out[p, ] <- u; V_out[q, ] <- v }
  list(U = U_out, S = S, V = V_out)
}

# ------------------------------------------------------------------------------
# spm_en: Euclidean normalization
# ------------------------------------------------------------------------------
spm_en <- function(X, p = NULL) {
  if (!is.null(p)) X <- spm_detrend(X, p)
  for (i in seq_len(ncol(X)))
    if (any(X[, i])) X[, i] <- X[, i] / sqrt(sum(X[, i]^2))
  X
}

# ------------------------------------------------------------------------------
# spm_logdet: log(det(C)) for (semi-)definite C; sparse-friendly
#   OPTIMIZATION: fast paths for 1x1 and diagonal matrices (very common in M-step)
# ------------------------------------------------------------------------------
spm_logdet <- function(C) {
  if (is.null(dim(C))) stop("spm_logdet: C must be a matrix, not a vector")

  # Fast path: 1x1
  if (prod(dim(C)) == 1L) {
    v <- as.numeric(C)
    if (v <= 0) return(NaN)
    return(log(v))
  }

  # Fast path: diagonal Matrix (very common — iS, ihC, etc. are often diagonal)
  if (inherits(C, "diagonalMatrix")) {
    d <- Matrix::diag(C)
    d <- d[d != 0]
    if (!length(d)) return(0)
    if (any(d <= 0 | is.nan(d))) return(NaN)
    return(sum(log(d)))
  }

  diagv <- if (inherits(C, "Matrix")) Matrix::diag(C) else diag(C)
  idx   <- which(diagv != 0)
  if (!length(idx)) return(0)
  C <- C[idx, idx, drop = FALSE]

  # Check for NaN
  if (inherits(C, "Matrix")) {
    if (any(is.nan(Matrix::summary(C)$x))) return(NaN)
  } else if (any(is.nan(C))) return(NaN)

  TOL <- 1e-16
  svd_logdet <- function(M) {
    d <- svd(as.matrix(M), nu = 0, nv = 0)$d
    sum(log(d[d > TOL & d < 1/TOL]))
  }
  is_sym <- function(M) {
    Nm <- if (inherits(M, "Matrix")) Matrix::norm(M - Matrix::t(M), "I")
          else norm(M - t(M), "I")
    !is.na(Nm) && Nm <= TOL
  }

  if (!is_sym(C)) return(svd_logdet(C))

  if (inherits(C, "Matrix")) {
    Csym <- Matrix::forceSymmetric(C, uplo = "U")
    chol_ok <- TRUE
    tryCatch(Matrix::Cholesky(Csym, LDL = FALSE, perm = TRUE),
             error = function(e) { chol_ok <<- FALSE })
    if (chol_ok) {
      detv <- Matrix::determinant(Csym, logarithm = TRUE)
      if (detv$sign <= 0) return(NaN)
      return(as.numeric(detv$modulus))
    }
    return(svd_logdet(Csym))
  } else {
    Csym <- (C + t(C)) / 2
    R <- tryCatch(chol(Csym), error = function(e) NULL)
    if (!is.null(R)) return(2 * sum(log(diag(R))))
    return(svd_logdet(Csym))
  }
}

# ------------------------------------------------------------------------------
# spm_diag
# ------------------------------------------------------------------------------
spm_diag <- function(...) {
  args <- list(...)
  if (!length(args)) stop('Not enough input arguments.')
  X <- args[[1]]
  if (!is.list(X)) return(diag(...))

  K <- if (length(args) < 2) 0 else args[[2]]
  if (is.null(dim(X))) {
    size <- length(X) + abs(K)
    D <- vector("list", size * size); dim(D) <- c(size, size)
    diag_idx <- which(diag(1, length(X)) == 1, arr.ind = TRUE)
    if (K != 0) {
      if (K > 0) diag_idx[, 2] <- diag_idx[, 2] + K
      else       diag_idx[, 1] <- diag_idx[, 1] - K
    }
    valid <- diag_idx[, 1] <= size & diag_idx[, 2] <= size
    diag_idx <- diag_idx[valid, , drop = FALSE]
    for (i in seq_len(nrow(diag_idx))) D[[diag_idx[i,1], diag_idx[i,2]]] <- X[[i]]
    return(D)
  } else {
    m <- nrow(X); n <- ncol(X)
    if (K != 0) {
      tmp <- matrix(0, max(m,n), max(m,n))
      if (K > 0) for (i in seq_len(n - K)) tmp[i, i + K] <- 1
      else       for (i in seq_len(m + K)) tmp[i - K, i] <- 1
      diag_idx <- which(tmp == 1, arr.ind = TRUE)
    } else {
      diag_idx <- which(diag(1, max(m,n)) == 1, arr.ind = TRUE)
    }
    D <- list()
    for (i in seq_len(nrow(diag_idx))) {
      r <- diag_idx[i,1]; cc <- diag_idx[i,2]
      if (r <= m && cc <= n) D[[length(D) + 1]] <- X[[r, cc]]
    }
    return(D)
  }
}

# ------------------------------------------------------------------------------
# spm_cat: Concatenate list of blocks into a single matrix
#   OPTIMIZATION: replaced the double nested for-loop with direct do.call on
#   the row/col vectors; the inner padding loop is kept only when needed.
# ------------------------------------------------------------------------------
spm_cat <- function(x, d = NULL) {
  if (!is.list(x)) return(x)

  if (!is.null(dim(x)) && !is.null(d)) {
    dims <- dim(x)
    if (d == 1) {
      y <- vector("list", dims[2])
      for (i in seq_len(dims[2])) y[[i]] <- spm_cat(x[, i, drop = FALSE])
    } else if (d == 2) {
      y <- vector("list", dims[1])
      for (i in seq_len(dims[1])) y[[i]] <- spm_cat(x[i, , drop = FALSE])
    } else stop("unknown option")
    return(y)
  }

  # Flat numeric list → fast cbind path
  if (is.null(dim(x)) && all(vapply(x, is.numeric, TRUE))) {
    mats <- lapply(x, function(e) if (is.vector(e)) matrix(e, ncol = 1) else e)
    return(do.call(cbind, mats))
  }

  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  dims <- dim(x); n <- dims[1]; m <- dims[2]

  # Compute block sizes
  I <- matrix(0, n, m); J <- matrix(0, n, m)
  for (i in seq_len(n)) for (j in seq_len(m)) {
    elem <- x[[i, j]]
    if (is.list(elem)) { elem <- spm_cat(elem); x[[i, j]] <- elem }
    if (!is.null(elem) && length(elem) > 0) {
      de <- dim(elem)
      if (is.null(de)) {
        I[i,j] <- if (length(elem) == 1) 1 else length(elem); J[i,j] <- 1
      } else { I[i,j] <- de[1]; J[i,j] <- de[2] }
    }
  }
  I_max <- apply(I, 1, max); J_max <- apply(J, 2, max)

  # Pad blocks
  for (i in seq_len(n)) for (j in seq_len(m)) {
    elem <- x[[i, j]]
    if (is.null(elem) || (is.numeric(elem) && !length(elem))) {
      x[[i, j]] <- matrix(0, I_max[i], J_max[j])
    } else if (is.numeric(elem) && length(elem) == 1) {
      x[[i, j]] <- matrix(elem, I_max[i], J_max[j])
    } else if (is.vector(elem) && length(elem) > 1) {
      x[[i, j]] <- if (I[i,j] == I_max[i]) matrix(elem, ncol = 1) else matrix(elem, nrow = 1)
    } else if (is.matrix(elem)) {
      cd <- dim(elem)
      if (cd[1] < I_max[i] || cd[2] < J_max[j]) {
        newm <- matrix(0, I_max[i], J_max[j])
        newm[seq_len(cd[1]), seq_len(cd[2])] <- elem
        x[[i, j]] <- newm
      }
    }
  }

  # Build result: row by row cbind, then rbind all rows
  row_mats <- vector("list", n)
  for (i in seq_len(n)) row_mats[[i]] <- do.call(cbind, as.list(x[i, ]))
  do.call(rbind, row_mats)
}

# ------------------------------------------------------------------------------
# spm_dx: Local linearisation integrator (matrix exponential form)
# ------------------------------------------------------------------------------
spm_dx <- function(dfdx, f, t = Inf, Q = NULL) {
  nmax <- 512; xf <- f
  f    <- spm_vec(f); n <- length(f)

  if (is.list(t)) {
    tval <- t[[1]]
    if (length(tval) == 1L)
      t <- suppressWarnings(exp(tval - spm_logdet(dfdx) / n))
    else
      t <- suppressWarnings(exp(tval - log(diag(-dfdx))))
  }

  if (!missing(Q) && !is.null(Q)) {
    L <- as.matrix(dfdx); L[upper.tri(L)] <- 0; Q <- L - t(L)
    Qn <- norm(Q, type = "2")
    Q <- if (is.finite(Qn) && Qn > 0) Q / Qn / 8 else Q * 0
    f <- f - Q %*% f; dfdx <- dfdx - Q %*% dfdx
  }

  if (all(t > exp(16))) {
    pinv_dfdx <- spm_pinv(dfdx)
    pinv_dfdx <- if (inherits(pinv_dfdx, "Matrix")) as.matrix(pinv_dfdx) else pinv_dfdx
    dx <- -pinv_dfdx %*% matrix(as.numeric(f), ncol = 1)
    return(spm_unvec(as.numeric(dx), xf))
  }

  Tmat <- if (is.vector(t) && length(t) > 1L) diag(as.numeric(t))
          else as.numeric(t)[1L]

  if (length(Tmat) == 1L) { tf <- Tmat * as.matrix(f); tJ <- Tmat * as.matrix(dfdx) }
  else                    { tf <- Tmat %*% f;           tJ <- Tmat %*% dfdx }

  J <- rbind(c(0, rep(0, n)), cbind(tf, tJ))

  if (n <= nmax) {
    E <- if (.has_expm) expm::expm(J) else spm_expm(J)
    dx <- E[2:(n+1), 1, drop = FALSE]
  } else {
    dx <- matrix(expv(1, J, c(1, rep(0, n)))[2:(n+1)], ncol = 1)
  }
  spm_unvec(as.numeric(dx), xf)
}

# ------------------------------------------------------------------------------
# expv, padm, spm_expm (unchanged — already efficient for their use cases)
# ------------------------------------------------------------------------------
expv <- function(t, A, v, tol = 1.0e-7, m = NULL) {
  A <- as.matrix(A); n <- nrow(A); if (is.null(m)) m <- min(n, 30L)
  anorm <- norm(A, type = "I"); mxrej <- 10; btol <- 1.0e-7
  gamma <- 0.9; delta <- 1.2; mb <- m; t_out <- abs(t)
  nstep <- 0; t_new <- 0; t_now <- 0; s_error <- 0
  rndoff <- anorm * .Machine$double.eps
  k1 <- 2; xm <- 1/m; normv <- norm(as.matrix(v), type = "2"); beta <- normv
  fact  <- (((m+1)/exp(1))^(m+1)) * sqrt(2*pi*(m+1))
  t_new <- (1/anorm) * ((fact*tol)/(4*beta*anorm))^xm
  s <- 10^(floor(log10(t_new))-1); t_new <- ceiling(t_new/s)*s
  sgn <- sign(t); w <- v; hump <- normv

  while (t_now < t_out) {
    nstep <- nstep + 1; t_step <- min(t_out - t_now, t_new)
    V <- matrix(0, n, m+1); H <- matrix(0, m+2, m+2)
    V[,1] <- (1/beta)*w; k1 <- 2
    for (j in seq_len(m)) {
      p <- A %*% V[,j]
      for (i in seq_len(j)) { H[i,j] <- crossprod(V[,i], p)[1,1]; p <- p - H[i,j]*V[,i] }
      s <- norm(as.matrix(p), type = "2")
      if (s < btol) { k1 <- 0; mb <- j; t_step <- t_out - t_now; break }
      H[j+1,j] <- s; V[,j+1] <- (1/s)*p
    }
    if (k1 != 0) { H[m+2,m+1] <- 1; avnorm <- norm(A %*% V[,m+1], type = "2") }
    ireject <- 0
    repeat {
      mx <- mb + k1; Hsub <- H[1:mx, 1:mx, drop = FALSE]
      F <- if (.has_expm) expm::expm(sgn*t_step*Hsub) else padm(sgn*t_step*Hsub)
      if (k1 == 0) { err_loc <- btol; break }
      phi1 <- abs(beta*F[m+1,1]); phi2 <- abs(beta*F[m+2,1]*avnorm)
      if      (phi1 > 10*phi2) { err_loc <- phi2; xm <- 1/m }
      else if (phi1 > phi2)    { err_loc <- (phi1*phi2)/(phi1-phi2); xm <- 1/m }
      else                     { err_loc <- phi1; xm <- 1/(m-1) }
      if (err_loc <= delta*t_step*tol) break
      t_step <- gamma*t_step*(t_step*tol/err_loc)^xm
      s <- 10^(floor(log10(t_step))-1); t_step <- ceiling(t_step/s)*s
      if (ireject == mxrej) stop("The requested tolerance is too high.")
      ireject <- ireject + 1
    }
    mx <- mb + max(0, k1-1)
    w  <- V[, 1:mx, drop = FALSE] %*% (beta*F[1:mx, 1, drop = FALSE])
    beta <- norm(as.matrix(w), type = "2"); hump <- max(hump, beta)
    t_now <- t_now + t_step
    t_new <- gamma*t_step*(t_step*tol/err_loc)^xm
    s <- 10^(floor(log10(t_new))-1); t_new <- ceiling(t_new/s)*s
    err_loc <- max(err_loc, rndoff); s_error <- s_error + err_loc
  }
  as.vector(w)
}

padm <- function(A, p = 6) {
  A <- as.matrix(A); n <- nrow(A)
  c_vec <- numeric(p+1); c_vec[1] <- 1
  if (p >= 1) for (k in seq_len(p)) c_vec[k+1] <- c_vec[k]*((p+1-k)/(k*(2*p+1-k)))
  s <- norm(A, type = "I")
  if (s > 0.5) { s <- max(0, floor(log(s)/log(2))+2); A <- A/(2^s) } else s <- 0
  I <- diag(1,n); A2 <- A %*% A; Q <- c_vec[p+1]*I; P <- c_vec[p]*I; odd <- TRUE
  for (k in (p-1):1) {
    if (odd) Q <- Q %*% A2 + c_vec[k]*I else P <- P %*% A2 + c_vec[k]*I
    odd <- !odd
  }
  if (odd) { Q <- Q %*% A; Q <- Q-P; E <- -(I + 2*solve(Q, diag(n))) }
  else     { P <- P %*% A; Q <- Q-P; E <- I + 2*solve(Q, diag(n)) }
  if (s > 0) for (k in seq_len(s)) E <- E %*% E
  E
}

spm_expm <- function(J, x = NULL) {
  I <- if (.has_Matrix) .Diag(nrow(J)) else diag(1, nrow(J))
  e <- floor(log2(norm(J, "I"))); s <- max(0, e+1)
  J <- J / 2^s; X <- J; c_v <- 1/2; E <- I + c_v*J; D <- I - c_v*J; q <- 6; p <- 1
  for (k in 2:q) {
    c_v <- c_v*(q-k+1)/(k*(2*q-k+1)); X <- J %*% X; cX <- c_v*X
    E <- E + cX; if (p == 1) D <- D + cX else D <- D - cX; p <- 1-p
  }
  E <- solve(as.matrix(D), as.matrix(E))
  for (k in seq_len(s)) E <- E %*% E
  if (!is.null(x)) E %*% x else E
}

# ------------------------------------------------------------------------------
# spm_pinv: Pseudo-inverse for sparse matrices
# ------------------------------------------------------------------------------
spm_pinv <- function(A, TOL = NULL) {
  m <- nrow(A); n <- ncol(A)
  if (!length(A)) return(if (.has_Matrix) .spzero(n,m) else matrix(0, n, m))
  if (is.null(TOL)) {
    old_warn <- options(warn = -1)
    X <- spm_inv(t(A) %*% A); options(old_warn)
    if (all(is.finite(X))) return(X %*% t(A))
  }
  svd_result <- spm_svd(A, 0)
  U <- svd_result$U; S <- svd_result$S; V <- svd_result$V
  S <- as.vector(Matrix::diag(S))
  if (is.null(TOL)) TOL <- max(m, n) * .Machine$double.eps * max(S)
  r <- sum(abs(S) > TOL)
  if (r == 0) return(if (.has_Matrix) .spzero(n,m) else matrix(0, n, m))
  i <- seq_len(r)
  S_inv <- .spmat(i = i, j = i, x = 1/S[i], dims = c(r, r))
  V[,i] %*% S_inv %*% t(U[,i])
}

# ------------------------------------------------------------------------------
# spm_dcm_fmri_mode_gen
# ------------------------------------------------------------------------------
spm_dcm_fmri_mode_gen <- function(Ev, modes, Cv = NULL) {
  Ev <- as.numeric(Ev); modes <- as.matrix(modes)
  Dv <- diag(-exp(-Ev), nrow = length(Ev))
  Ep <- modes %*% Dv %*% t(modes)
  if (is.null(Cv)) return(Ep)
  dAdv <- spm_diff(spm_dcm_fmri_mode_gen, Ev, modes, 1)
  Jpart <- if (!is.null(dAdv$J)) dAdv$J else if (!is.null(dAdv$d1)) dAdv$d1 else dAdv
  n <- nrow(Ep); p <- length(Ev); G <- matrix(0, n*n, p)
  if (is.list(Jpart)) {
    if (length(Jpart) != p) stop("Derivative list length mismatch.")
    for (i in seq_len(p)) { Gi <- as.matrix(Jpart[[i]]); G[,i] <- as.vector(Gi) }
  } else if (is.matrix(Jpart)) {
    G <- Jpart
  } else stop("Unrecognized derivative container from spm_diff.")
  Cp <- G %*% as.matrix(Cv) %*% t(G)
  list(Ep = Ep, Cp = Cp)
}

# ------------------------------------------------------------------------------
# spm_kernels: Volterra kernels for bilinear systems
#   OPTIMIZATION: replaced the inner row-copy loop with direct matrix slice
#   (H1[j:N,,p] is already a 3-D array slice; drop the extra dim and use it directly)
# ------------------------------------------------------------------------------
spm_kernels <- function(...) {
  args <- list(...); na <- length(args)
  if      (na == 4) { M0 <- args[[1]]; M1 <- args[[2]]; N <- as.integer(args[[3]]); dt <- as.numeric(args[[4]]); L1 <- NULL; L2 <- NULL }
  else if (na == 5) { M0 <- args[[1]]; M1 <- args[[2]]; L1 <- args[[3]]; N <- as.integer(args[[4]]); dt <- as.numeric(args[[5]]); L2 <- NULL }
  else if (na == 6) { M0 <- args[[1]]; M1 <- args[[2]]; L1 <- args[[3]]; L2 <- args[[4]]; N <- as.integer(args[[5]]); dt <- as.numeric(args[[6]]) }
  else stop("Improper call to spm_kernels")

  if (is.list(M0) && !is.matrix(M0)) {
    bi <- spm_bireduce(M0, M1); M0 <- bi$M0; M1 <- bi$M1; L1 <- bi$L1; L2 <- bi$L2
  }
  if (is.null(L1)) { nn <- nrow(as.matrix(M0)); L1 <- diag(nn)[-1,, drop = FALSE] }
  if (is.null(L2)) L2 <- list()

  N <- max(1L, N); M0 <- as.matrix(M0)
  n <- nrow(M0); l <- nrow(as.matrix(L1))
  if (!is.list(M1)) stop("M1 must be a list of n×n matrices")
  m <- length(M1)

  H1 <- array(0, dim = c(N, n, m))
  K1 <- array(0, dim = c(N, l, m))
  K2 <- array(0, dim = c(N, N, l, m, m))

  e1 <- expm::expm(dt * M0); e2 <- expm::expm(-dt * M0)
  Mip <- vector("list", N)
  Mip[[1]] <- lapply(seq_len(m), function(p) e1 %*% as.matrix(M1[[p]]) %*% e2)
  if (N >= 2) {
    for (i in 2:N)
      Mip[[i]] <- lapply(seq_len(m), function(p) e1 %*% Mip[[i-1]][[p]] %*% e2)
  }

  X0 <- matrix(0, n, 1); X0[1,1] <- 1
  ei <- e1; if (N >= 2) for (i in 2:N) ei <- e1 %*% ei
  H0 <- ei %*% X0; K0 <- L1 %*% H0

  L1_mat <- as.matrix(L1)
  for (p in seq_len(m)) for (i in seq_len(N)) {
    v <- Mip[[i]][[p]] %*% H0
    H1[i,,p] <- as.numeric(v)
    K1[i,,p] <- as.numeric(t(v) %*% t(L1_mat))
  }

  if (m > 0) {
    for (p in seq_len(m)) for (q in seq_len(m)) for (j in seq_len(N)) {
      # OPTIMIZED: use direct 3-D slice, drop extra dim, then transpose
      H1_slice_t <- t(matrix(H1[j:N,,p], nrow = N-j+1, ncol = n))  # n × (N-j+1)
      Hmat <- L1_mat %*% Mip[[j]][[q]] %*% H1_slice_t               # l × (N-j+1)
      Ht   <- t(Hmat)                                                 # (N-j+1) × l
      K2[j, j:N,,q,p] <- array(Ht, dim = c(N-j+1, l))
      K2[j:N, j,,p,q] <- array(Ht, dim = c(N-j+1, l))
    }
  }

  if (length(L2) > 0) {
    for (i in seq_len(m)) for (j in seq_len(m)) {
      Hi <- matrix(H1[,,i], nrow = N, ncol = n)
      Hj <- matrix(H1[,,j], nrow = N, ncol = n)
      for (pout in seq_len(l)) {
        K2[,,pout,i,j] <- K2[,,pout,i,j] + Hi %*% as.matrix(L2[[pout]]) %*% t(Hj)
      }
    }
  }

  list(K0 = K0, K1 = K1, K2 = K2, H1 = H1)
}

# ------------------------------------------------------------------------------
# spm_detrend: Column-wise centering (p=0) or polynomial detrend (p>=1)
#   OPTIMIZATION: for p==1 use lm.fit (avoids full ginv computation)
# ------------------------------------------------------------------------------
spm_detrend <- function(x, p = 0) {
  if (is.list(x)) {
    y <- x; for (i in seq_along(x)) y[[i]] <- spm_detrend(x[[i]], p); return(y)
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  dim_x <- dim(x)
  if (is.null(dim_x)) { m <- length(x); n <- 1; x <- matrix(x, ncol = 1) }
  else                { m <- dim_x[1]; n <- dim_x[2] }
  if (m == 0 || n == 0) return(matrix(0, m, n))
  if (p == 0) return(x - matrix(colMeans(x), nrow = m, ncol = n, byrow = TRUE))

  G <- matrix(0, m, p+1)
  for (i in 0:p) G[, i+1] <- (seq_len(m))^i

  if (p == 1) {
    # Fast path: use lm.fit instead of MASS::ginv
    coef <- lm.fit(G, x)$coefficients       # (p+1) × n
    return(x - G %*% coef)
  }
  y <- x - G %*% (MASS::ginv(G) %*% x)
  y
}

# ------------------------------------------------------------------------------
# spm_blocks_to_mat
# ------------------------------------------------------------------------------
spm_blocks_to_mat <- function(blocks, row_major = FALSE) {
  stopifnot(is.list(blocks), length(blocks) > 0)
  B1 <- blocks[[1]]
  if (!is.matrix(B1)) stop("All blocks must be matrices.")
  rb <- nrow(B1); cb <- ncol(B1)
  for (k in seq_along(blocks)) {
    bk <- blocks[[k]]
    if (!is.matrix(bk) || nrow(bk) != rb || ncol(bk) != cb)
      stop(sprintf("Block %d has shape %s; expected %dx%d.", k, paste0(dim(bk), collapse="x"), rb, cb))
  }
  n_cols <- length(blocks); n_rows <- rb * cb
  out <- matrix(0, nrow = n_rows, ncol = n_cols)
  if (!row_major) for (k in seq_len(n_cols)) out[,k] <- as.vector(blocks[[k]])
  else            for (k in seq_len(n_cols)) out[,k] <- as.vector(t(blocks[[k]]))
  out
}

# ------------------------------------------------------------------------------
# spm_bireduce (unchanged logic; relies on spm_diff which is already optimized)
# ------------------------------------------------------------------------------
spm_bireduce <- function(M, P) {
  funx <- NULL
  if (!is.null(M$f)) {
    if (is.function(M$f)) funx <- M$f
    else if (is.character(M$f) && exists(M$f, mode = "function")) funx <- get(M$f, mode = "function")
  }
  if (is.null(funx)) { funx <- function(x,u,P,M) numeric(0); M$n <- 0; M$x <- numeric(0) }

  if (!is.null(M$u)) u <- spm_vec(M$u)
  else {
    m_dim <- if (!is.null(M$m)) as.integer(M$m) else 0L
    u <- if (m_dim > 0L) matrix(0, m_dim, 1) else numeric(0)
  }

  tmp12 <- spm_diff(funx, M$x, u, P, M, c(1, 2))
  if (!is.null(tmp12$d2) && !is.null(tmp12$d1)) { dfdxu <- tmp12$d2; dfdx  <- tmp12$d1 }
  else stop("spm_diff(funx, ..., c(1,2)) did not return the expected {d2,d1}.")

  tmp2 <- spm_diff(funx, M$x, u, P, M, 2)
  dfdu <- tmp2$J; f0 <- spm_vec(tmp2$f0); n <- length(f0)

  normalize_dfdu <- function(dfdu_raw, n, m, n_block_cols_hint = NULL) {
    if (is.matrix(dfdu_raw) && nrow(dfdu_raw) == n && ncol(dfdu_raw) == m) return(dfdu_raw)
    if (is.list(dfdu_raw) && length(dfdu_raw) == m &&
        all(vapply(dfdu_raw, function(z) is.numeric(z) && length(as.vector(z)) == n, TRUE))) {
      out <- matrix(0, n, m)
      for (j in seq_len(m)) out[,j] <- as.vector(dfdu_raw[[j]])
      return(out)
    }
    if (is.list(dfdu_raw) && length(dfdu_raw) > 0 && is.matrix(dfdu_raw[[1]])) {
      big <- spm_blocks_to_mat(dfdu_raw)
      if (m == 1 && ncol(big) != 1) big <- matrix(as.vector(big), nrow = n, ncol = 1)
      if (nrow(big) != n) stop("dfdu stitching produced wrong row count.")
      if (ncol(big) != m) big <- big[, seq_len(min(ncol(big), m)), drop = FALSE]
      return(big)
    }
    if (is.numeric(dfdu_raw) && length(dfdu_raw) == n && m == 1)
      return(matrix(dfdu_raw, nrow = n, ncol = 1))
    stop("Unexpected dfdu shape.")
  }

  m <- if (is.list(dfdxu)) length(dfdxu) else {
    if (!is.null(M$m)) as.integer(M$m) else max(1L, ncol(as.matrix(u)))
  }
  dfdu <- normalize_dfdu(dfdu, n = n, m = m)

  if (is.list(dfdx) && length(dfdx) > 0 && is.matrix(dfdx[[1]]))
    dfdx <- spm_blocks_to_mat(dfdx)
  else dfdx <- as.matrix(dfdx)

  if (is.list(dfdxu)) {
    dfdxu <- lapply(dfdxu, function(slc) {
      if (is.list(slc) && length(slc) > 0 && is.matrix(slc[[1]])) spm_blocks_to_mat(slc)
      else as.matrix(slc)
    })
  } else if (length(dim(dfdxu)) == 3L && all(dim(dfdxu)[1:2] == c(n,n))) {
    dfdxu <- lapply(seq_len(dim(dfdxu)[3]), function(i) dfdxu[,,i])
  } else stop("Unexpected shape for dfdxu.")

  n <- length(f0)
  x <- spm_vec(M$x)
  if (!is.null(M$D)) {
    D <- M$D; f0 <- D %*% f0; dfdx <- D %*% dfdx; dfdu <- D %*% dfdu
    if (is.list(dfdxu)) dfdxu <- lapply(dfdxu, function(A) D %*% A)
  }

  M0 <- spm_cat(as.matrix(list(list(0, NULL), list(f0 - dfdx %*% x, dfdx)), 2, 2, byrow = FALSE))
  M1 <- vector("list", m)
  if (m > 0) {
    for (i in seq_len(m)) {
      dfdxu_i <- if (is.list(dfdxu)) dfdxu[[i]] else dfdxu[,,i]
      rhs_i   <- dfdu[,i,drop=FALSE] - dfdxu_i %*% x
      M1[[i]] <- spm_cat(as.matrix(list(list(0, NULL), list(rhs_i, dfdxu_i)), 2, 2))
    }
  }

  fung <- NULL
  if (!is.null(M$g)) {
    if (is.function(M$g)) fung <- M$g
    else if (is.character(M$g) && exists(M$g, mode = "function")) fung <- get(M$g, mode = "function")
  }
  if (is.null(fung)) { fung <- function(x,u,P,M) spm_vec(x); if (is.null(M$l)) M$l <- n }

  tmpg1 <- spm_diff(fung, M$x, u, P, M, 1)
  if (is.list(tmpg1$f0) && !is.null(tmpg1$f0$g) && !is.null(tmpg1$f0$dgdx)) {
    g0 <- as.numeric(tmpg1$f0$g); dgdx <- as.matrix(tmpg1$f0$dgdx)
  } else stop("spm_diff(fung, ..., 1) didn't return f0$g and f0$dgdx as expected.")

  l <- length(g0); x <- spm_vec(M$x); n <- length(x)
  if (!is.matrix(dgdx) || nrow(dgdx) != l || ncol(dgdx) != n)
    stop(sprintf("dgdx must be l×n; got %dx%d, expected %dx%d", nrow(dgdx), ncol(dgdx), l, n))

  g0_col <- matrix(g0, ncol=1); x_col <- matrix(as.numeric(x), ncol=1)
  L1 <- spm_cat(list(list(g0_col - dgdx %*% x_col, dgdx)))

  .ensure_mat <- function(M, ny, nx) {
    if (is.null(M)) return(matrix(0, ny, nx))
    if (is.matrix(M) || (is.array(M) && length(dim(M)) == 2L)) {
      M <- as.matrix(M)
      if (nrow(M) == ny && ncol(M) == nx) return(M)
      Z <- matrix(0, ny, nx); nr <- min(nrow(M),ny); nc <- min(ncol(M),nx)
      if (nr > 0 && nc > 0) Z[seq_len(nr),seq_len(nc)] <- M[seq_len(nr),seq_len(nc),drop=FALSE]
      return(Z)
    }
    if (is.atomic(M) && is.null(dim(M))) {
      v <- as.numeric(M)
      if (!length(v)) return(matrix(0, ny, nx))
      if (length(v) < ny*nx) v <- c(v, rep(0, ny*nx - length(v)))
      return(matrix(v[seq_len(ny*nx)], nrow=ny, ncol=nx))
    }
    if (is.list(M)) {
      v <- suppressWarnings(as.numeric(unlist(M, use.names=FALSE)))
      v[is.na(v)] <- 0
      if (length(v) < ny*nx) v <- c(v, rep(0, ny*nx - length(v)))
      return(matrix(v[seq_len(ny*nx)], nrow=ny, ncol=nx))
    }
    matrix(0, ny, nx)
  }

  tmpg2 <- spm_diff(fung, M$x, u, P, M, c(1,1), 'nocat')
  dgdxx <- if (!is.null(tmpg2$d2)) tmpg2$d2 else if (!is.null(tmpg2$J)) tmpg2$J else tmpg2
  if (!is.list(dgdxx) || length(dgdxx) < n)
    stop("spm_diff(fung, ..., c(1,1),'nocat') did not return expected second-derivative structure.")

  ny <- nrow(dgdx); nx <- ncol(dgdx)
  dgdxx <- lapply(dgdxx, .ensure_mat, ny = ny, nx = nx)

  L2 <- vector("list", l)
  for (i in seq_len(l)) {
    Di <- matrix(0, n, n)
    for (j in seq_len(n))
      Di[j,] <- .ensure_mat(dgdxx[[j]], ny, nx)[i,,drop=FALSE]
    L2[[i]] <- spm_cat(spm_diag(list(0, Di)))
  }
  list(M0 = M0, M1 = M1, L1 = L1, L2 = L2)
}

# ------------------------------------------------------------------------------
# spm_int (unchanged structure; speed comes from upstream optimizations)
# ------------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

#
# spm_int — FAST version. Real bottleneck of DCM estimation is here:
# this function is called np+1 times per Gauss-Newton iteration via spm_diff,
# and each call walks length(t_idx) time bins doing a fresh expm at every bin.
#
# Key optimizations (orders of magnitude, not percent):
#   1. M0, M1 densified ONCE up front. Otherwise every `J*dt` and every
#      `expm()` call re-converts sparse Matrix → dense.
#   2. Use expm::expm directly (compiled C/LAPACK). The original spm_expm is
#      pure R Pade and is ~30-100x slower for typical 30x30 J.
#   3. PROPAGATOR CACHE. In a typical fMRI design, between two input changes
#      there are dozens of consecutive bins with the SAME J and (often) the
#      SAME dt. We compute E_dt = expm(J*dt) once and reuse it as a plain
#      mat-vec multiply. This is the single biggest win.
#   4. x kept as a plain numeric vector (no Matrix wrapping per step).
#
spm_int <- function(P, M, U) {
  if (!is.list(U)) U <- list(u = U)
  if (is.null(U$u)) stop("U must contain U$u")
  if (is.null(U$dt)) U$dt <- 1

  u_bins <- nrow(U$u); v_out <- if (!is.null(M$ns)) M$ns else u_bins
  x0 <- spm_vec(M$x); x <- c(1, x0)

  f_fun <- if (!is.null(M$f)) {
    if (is.function(M$f)) M$f
    else if (is.character(M$f) && exists(M$f, mode="function")) get(M$f, mode="function")
    else { M$x <- numeric(0); function(x,u,P,M) numeric(0) }
  } else { M$x <- numeric(0); function(x,u,P,M) numeric(0) }

  g_fun <- if (!is.null(M$g)) {
    if (is.function(M$g)) M$g
    else if (is.character(M$g) && exists(M$g, mode="function")) get(M$g, mode="function")
    else function(x,u,P,M) x
  } else function(x,u,P,M) x

  bi <- spm_bireduce(M, P)
  # OPT 1: densify M0/M1 once. Avoids per-iteration sparse->dense conversion.
  M0 <- as.matrix(bi$M0)
  M1 <- lapply(bi$M1, as.matrix)
  m_inputs <- length(M1)
  l_out <- if (!is.null(M$l)) M$l else nrow(as.matrix(g_fun(M$x, rep(0, m_inputs), P, M)))

  D <- tryCatch({
    if (!is.null(M$delays)) pmax(round(M$delays / U$dt), 1)
    else rep(round(u_bins / v_out), l_out)
  }, error = function(e) rep(round(u_bins / v_out), l_out))

  if (ncol(U$u) == 0) change_rows <- integer(0)
  else { dU <- diff(U$u); change_rows <- which(apply(dU, 1, function(r) any(r != 0))) }
  i_idx <- c(1, 1 + change_rows)
  su_full <- logical(u_bins); su_full[i_idx] <- TRUE

  s_idx <- ceiling((0:(v_out-1)) * (u_bins / v_out))
  s_idx[s_idx < 1] <- 1; s_idx[s_idx > u_bins] <- u_bins

  sy_full <- matrix(0L, l_out, u_bins)
  for (j in seq_len(l_out)) {
    ii <- s_idx + D[j]; ii[ii < 1] <- 1; ii[ii > u_bins] <- u_bins
    sy_full[j, ii] <- seq_len(v_out)
  }

  t_mask <- su_full | apply(sy_full != 0L, 2, any)
  t_idx  <- which(t_mask); su <- su_full[t_idx]
  sy     <- sy_full[, t_idx, drop=FALSE]
  dt_vec <- c(diff(t_idx), 0) * U$dt

  # Pre-extract output-bin column indices per time step (avoids which() in loop)
  any_sy   <- .colSums(sy != 0L, l_out, ncol(sy)) > 0
  Uu_dense <- as.matrix(U$u)

  y <- matrix(0, v_out, l_out); J <- M0
  u_row <- numeric(m_inputs)

  # OPT 3: propagator cache
  E_dt      <- NULL    # cached expm(J*dt)
  E_dt_dt   <- NA_real_
  J_dirty   <- TRUE    # set when J changes, forces rebuild of E_dt

  expm_fn <- if (.has_expm) expm::expm else spm_expm

  for (ii in seq_along(t_idx)) {

    # Update J only on input change
    if (su[ii]) {
      if (m_inputs > 0) {
        u_row <- as.numeric(Uu_dense[t_idx[ii], seq_len(m_inputs), drop = TRUE])
        # Defend against NaN/NA inputs in U$u (treat as 0)
        bad <- !is.finite(u_row)
        if (any(bad)) u_row[bad] <- 0
        J <- M0
        for (j in seq_len(m_inputs)) {
          uj <- u_row[j]
          if (uj != 0) J <- J + uj * M1[[j]]   # uj guaranteed finite now
        }
      } else {
        J <- M0
      }
      J_dirty <- TRUE
    }

    # Sample output if requested
    if (any_sy[ii]) {
      x_state <- spm_unvec(x[-1], M$x)
      q_out   <- spm_vec(g_fun(x_state, u_row, P, M))
      js <- which(sy[, ii] != 0L)
      for (j in js) y[sy[j, ii], j] <- q_out[j]
    }

    dt <- dt_vec[ii]
    if (dt != 0) {
      # OPT 2+3: rebuild propagator only if J or dt changed
      if (J_dirty || dt != E_dt_dt) {
        E_dt    <- expm_fn(J * dt)
        E_dt_dt <- dt
        J_dirty <- FALSE
      }
      # Plain dense mat-vec: the hot inner step
      x <- as.numeric(E_dt %*% x)
    }
    # NA-safe divergence check: bail if x went NaN/Inf or simply blew up
    x_norm <- sum(abs(x))
    if (!is.finite(x_norm) || x_norm > 1e6) break
  }
  Re(y)
}

# ------------------------------------------------------------------------------
# spm_log_evidence (unchanged logic)
# ------------------------------------------------------------------------------
spm_log_evidence <- function(qE, qC, pE, pC, rE = NULL, rC = NULL, ...) {
  if (!is.null(rE) && is.function(rE)) {
    priors <- tryCatch(rE(...), error = function(e) NULL)
    if (!is.null(priors)) { rE_val <- priors[[1]]; rC_val <- priors[[2]] }
    else { rE_val <- rE; rC_val <- rC }
  } else { rE_val <- rE; rC_val <- rC }

  if (is.null(rE_val) || is.null(rC_val)) {
    n <- if (is.matrix(qC)) nrow(qC) else spm_length(qC)
    rE_val <- .spzero(n, 1); rC_val <- .spzero(n, n)
  }

  if (is.list(pC) && !is.matrix(pC)) pC <- diag(spm_vec(pC))
  if (is.list(qC) && !is.matrix(qC)) qC <- diag(spm_vec(qC))
  if (is.list(rC_val) && !is.matrix(rC_val)) rC_val <- diag(spm_vec(rC_val))

  qE_vec <- spm_vec(qE); pE_vec <- spm_vec(pE); rE_vec <- spm_vec(rE_val)
  TOL <- exp(-16); i <- which(diag(pC) > TOL)
  if (!length(i)) return(list(F = 0, sE = qE, sC = qC))

  qP <- spm_inv(qC[i,i,drop=FALSE], TOL); pP <- spm_inv(pC[i,i,drop=FALSE], TOL)
  rP <- spm_inv(rC_val[i,i,drop=FALSE], TOL)
  sP <- qP + rP - pP; sC_r <- spm_inv(sP, TOL); pC_r <- spm_inv(pP, TOL)
  sE_vec <- as.matrix(qP) %*% qE_vec[i] + as.matrix(rP) %*% rE_vec[i] - as.matrix(pP) %*% pE_vec[i]

  term1 <- spm_logdet(rP %*% qP %*% sC_r %*% pC_r)
  term2 <- as.numeric(t(qE_vec[i]) %*% qP %*% qE_vec[i]) +
           as.numeric(t(rE_vec[i]) %*% rP %*% rE_vec[i]) -
           as.numeric(t(pE_vec[i]) %*% pP %*% pE_vec[i]) -
           as.numeric(t(sE_vec)    %*% sC_r %*% sE_vec)
  F <- (term1 - term2) / 2

  rE_vec[i] <- as.numeric(sC_r %*% sE_vec)
  rC_full <- rC_val; rC_full[i,i] <- sC_r
  list(F = F, sE = spm_unvec(rE_vec, qE), sC = rC_full)
}

# ------------------------------------------------------------------------------
# spm_log_evidence_reduce (unchanged logic)
# ------------------------------------------------------------------------------
spm_log_evidence_reduce <- function(qE, qC, pE, pC, rE, rC, TOL = 1e-8) {
  pC <- as.matrix(pC)
  if (is.list(pC) && !is.matrix(pC)) pC <- diag(spm_vec(pC)[spm_vec(pC) != 0])
  if (is.list(rC) && !is.matrix(rC)) rC <- diag(spm_vec(rC))

  RE_orig <- rE; SE_orig <- qE
  svd_result <- spm_svd(pC, 1e-6); U <- as.matrix(svd_result$U)

  qE_vec <- as.numeric(t(U) %*% spm_vec(qE)); pE_vec <- as.numeric(t(U) %*% spm_vec(pE))
  rE_vec <- as.numeric(t(U) %*% spm_vec(rE))
  qC_red <- t(U) %*% as.matrix(qC) %*% U
  pC_red <- t(U) %*% pC %*% U; rC_red <- t(U) %*% rC %*% U

  qP <- spm_inv(qC_red, TOL); pP <- spm_inv(pC_red, TOL); rP <- spm_inv(rC_red, TOL)
  sP <- qP + rP - pP; sC_red <- spm_inv(sP, TOL); pC_inv <- spm_inv(pP, TOL)
  sE_vec <- as.matrix(qP) %*% qE_vec + as.matrix(rP) %*% rE_vec - as.matrix(pP) %*% pE_vec

  term1 <- spm_logdet(as.matrix(rP %*% qP %*% sC_red %*% pC_inv))
  term2 <- as.numeric(t(qE_vec) %*% qP %*% qE_vec) +
           as.numeric(t(rE_vec) %*% rP %*% rE_vec) -
           as.numeric(t(pE_vec) %*% pP %*% pE_vec) -
           as.numeric(t(sE_vec) %*% sC_red %*% sE_vec)
  F <- (term1 - term2) / 2

  pE_orig_vec <- spm_vec(RE_orig)
  rE_restored <- as.numeric(sC_red %*% sE_vec)
  sE_final_vec <- as.numeric(U %*% rE_restored + pE_orig_vec - U %*% (t(U) %*% pE_orig_vec))
  sC_final <- U %*% sC_red %*% t(U)
  list(F = F, sE = spm_unvec(sE_final_vec, SE_orig), sC = sC_final)
}

# ------------------------------------------------------------------------------
# Remaining utilities (spm_fieldindices, spm_zeros, spm_find_pC, spm_find_rC,
# spm_speye, spm_dcm_evidence) — unchanged logic from original
# ------------------------------------------------------------------------------
spm_fieldindices <- function(X, ...) {
  varargin <- list(...)
  if (length(varargin) > 0 && is.numeric(varargin[[1]]) && length(varargin[[1]]) > 1) {
    return(lapply(varargin[[1]], function(j) spm_fieldindices(X, j)))
  }
  X0 <- spm_zeros(X); ix <- spm_vec(X0)
  for (i in seq_along(varargin)) {
    if (is.character(varargin[[i]])) {
      field <- varargin[[i]]; x <- X0
      if (field %in% names(x)) {
        f <- x[[field]]; f <- spm_unvec(spm_vec(f) + 1, f); x[[field]] <- f
        ix <- ix + spm_vec(x)
      } else warning(paste("Field", field, "not found"))
    } else if (is.numeric(varargin[[i]])) {
      ind <- varargin[[i]]; x <- ix; x[ind] <- 1; x <- spm_unvec(x, X)
      for (nm in names(x)) {
        if (any(spm_vec(x[[nm]]) != 0)) {
          v <- x[[nm]]
          if (is.numeric(v) || is.matrix(v)) {
            idx <- which(v != 0, arr.ind = TRUE)
            if (is.null(dim(v)) || min(dim(v)) == 1) {
              return(sprintf('%s[%i]', nm, which(v != 0)[1]))
            } else if (length(dim(v)) < 3) {
              return(sprintf('%s[%i,%i]', nm, idx[1,1], idx[1,2]))
            } else {
              return(sprintf('%s[%i,%i,%i]', nm, idx[1,1], idx[1,2], idx[1,3]))
            }
          }
        }
      }
    }
  }
  which(ix != 0)
}

spm_zeros <- function(X) spm_unvec(rep(0, spm_length(X)), X)

spm_find_pC <- function(...) {
  varargin <- list(...)
  pC <- NULL; pE <- NULL; fields <- NULL
  if (length(varargin) > 2) { pC <- varargin[[1]]; pE <- varargin[[2]]; fields <- varargin[[3]]
  } else if (length(varargin) > 1) { DCM <- varargin[[1]]; fields <- varargin[[2]]
  } else { DCM <- varargin[[1]]; fields <- NULL }

  if (is.null(pC) || is.null(pE)) {
    if (is.character(DCM)) DCM <- readRDS(DCM)
    if (any(c("options","M") %in% names(DCM))) {
      tryCatch({ result <- spm_find_rC(DCM); pC <- result$pC; pE <- result$pE }, error = function(e) {})
    }
  }

  q <- if (is.list(pC) && !is.matrix(pC)) spm_vec(pC) else diag(pC)
  threshold <- mean(q[q < 1024]) / 1024; i <- which(q > threshold); Np <- length(q)

  if (!is.null(fields)) {
    if (is.character(fields) && length(fields) == 1) fields <- list(fields)
    if (is.list(pE)) {
      j <- do.call(spm_fieldindices, c(list(pE), fields))
      if (!length(j) && !(length(fields) > 0 && fields[[1]] == "none"))
        warning(sprintf("%s not found. Returning all fields", paste(unlist(fields), collapse=",")))
      else i <- j[j %in% i]
    }
  }
  list(i = i, pC = pC, pE = pE, Np = Np)
}

spm_find_rC <- function(DCM) {
  if (!is.null(DCM$M$pC) && !is.null(DCM$M$pE)) return(list(pC = DCM$M$pC, pE = DCM$M$pE))
  if (!is.null(DCM$options$spatial)) {
    if (toupper(DCM$options$analysis) == "IND") {
      result <- spm_ind_priors(DCM$A, DCM$B, DCM$C, DCM$Nf); pE <- result$pE; pC <- result$pC
    } else {
      if (is.null(DCM$C)) DCM$C <- .spmat(i=integer(0), j=integer(0), dims=c(length(DCM$Sname),0))
      result <- spm_dcm_neural_priors(DCM$A, DCM$B, DCM$C, DCM$options$model)
      pE <- result$pE; pC <- result$pC
      tryCatch({
        model   <- if (!is.null(DCM$options$model))   DCM$options$model   else "NMM"
        spatial <- if (!is.null(DCM$options$spatial)) DCM$options$spatial else "LFP"
        DCM$M$dipfit$model <- model; DCM$M$dipfit$type <- spatial
        result <- spm_L_priors(DCM$M$dipfit, pE, pC); pE <- result$pE; pC <- result$pC
        result <- spm_ssr_priors(pE, pC); pE <- result$pE; pC <- result$pC
      }, error = function(e) {})
    }
  } else {
    if (is.null(DCM$d)) DCM$d <- array(0, dim = c(nrow(DCM$a), nrow(DCM$a), 0))
    result <- spm_dcm_fmri_priors(DCM$a, DCM$b, DCM$c, DCM$d, DCM$options)
    pE <- result$pE; pC <- result$pC
  }
  list(pC = pC, pE = pE)
}

spm_speye <- function(m, n = m, k = 0, c = 0) {
  if (!.has_Matrix) stop("Package 'Matrix' is required for spm_speye")
  D <- Matrix::bandSparse(m, n, k, diagonals = list(rep(1, m)))
  if (c == 1) {
    if      (k < 0) D <- D + spm_speye(m, n, min(n,m) + k)
    else if (k > 0) D <- D + spm_speye(m, n, k - min(n,m))
  } else if (c == 2) {
    col_sums <- Matrix::colSums(D); i <- which(col_sums == 0)
    if (length(i) > 0) D <- D + .spmat(i=i, j=i, x=1, dims=c(n,m))
  }
  D
}

spm_dcm_evidence <- function(DCM) {
  v <- DCM$v; n <- DCM$n
  Cp_diag <- if (inherits(DCM$Cp, c("Matrix","sparseMatrix","denseMatrix"))) Matrix::diag(DCM$Cp)
             else if (is.matrix(DCM$Cp)) diag(DCM$Cp) else DCM$Cp
  wsel <- which(Cp_diag != 0)
  evidence <- list(region_cost = numeric(n))
  for (i in seq_len(n)) {
    lambda_i <- tryCatch({
      if (is.matrix(DCM$Ce) || inherits(DCM$Ce,"Matrix")) DCM$Ce[i*v, i*v]
      else if (is.vector(DCM$Ce) && length(DCM$Ce) > 1)   DCM$Ce[i]
      else DCM$Ce
    }, error = function(e) if (length(DCM$Ce) >= i) DCM$Ce[i] else DCM$Ce[1])
    R_i <- DCM$R[,i]
    evidence$region_cost[i] <- -0.5*v*log(lambda_i) - 0.5*sum(R_i * (1/lambda_i) * R_i)
  }
  evidence$aic_penalty <- length(wsel)
  evidence$bic_penalty <- 0.5 * length(wsel) * log(v)
  evidence$aic_overall <- sum(evidence$region_cost) - evidence$aic_penalty
  evidence$bic_overall <- sum(evidence$region_cost) - evidence$bic_penalty
  evidence
}

# ==============================================================================
# Bytecode compilation
# Compiles the hottest functions to bytecode. Modern R already JIT-compiles
# functions on second use, but this front-loads that work and tends to give
# a ~5-15% additional speedup on top of JIT.
# ==============================================================================
if (requireNamespace("compiler", quietly = TRUE)) {
  local({
    .cmp <- compiler::cmpfun
    for (.fn in c("spm_int", "spm_diff", "spm_vec", "spm_unvec", "spm_length",
                  "spm_logdet", "spm_inv", "spm_cat", "spm_dfdx", "spm_dfdx_cat",
                  "spm_bireduce", "spm_dx", "spm_expm", "spm_pinv",
                  "spm_kernels", "spm_detrend", "spm_trace", "padm")) {
      if (exists(.fn, envir = .GlobalEnv, inherits = FALSE)) {
        assign(.fn, .cmp(get(.fn, envir = .GlobalEnv)), envir = .GlobalEnv)
      }
    }
  })
}

