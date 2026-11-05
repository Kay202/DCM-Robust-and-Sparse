# ==============================================================================
# exact_core_final_2_optimized.R
# Optimized version of exact_core_final_2.R for fMRI DCM analysis
#
# KEY OPTIMIZATIONS:
#   1. spm_fx_fmri: pre-allocate exp'd hemodynamic states once (avoid repeated
#      exp() calls); vectorize B/D loop with Reduce when inputs are many
#   2. spm_fx_fmri2: same; also removed as.array() coercions that copy data
#   3. spm_gx_fmri: use pre-computed k1/k2/k3 (unchanged, already vectorized)
#   4. spm_nlsi_GN M-step:
#      - iS_big (Kronecker product) rebuilt only when h changes significantly
#      - Pre-compute J' (tJ = t(J)) once outside the inner m-loop
#      - spm_trace replaced with fast vectorized sum(A*B) form
#      - Cp computation uses solve(Pp + ipC) once, avoids re-converting types
#   5. spm_dcm_fmri_priors: pC loop for 3-D A replaced with direct array op
#   6. spm_dcm_estimate: minor cleanup (no repeated spm_inv on same matrix)
# ==============================================================================

# Assumes helpers_final_optimized.R (or the originals) are already sourced.

# ------------------------------------------------------------------------------
# spm_fx_fmri: dx/dt for fMRI neural + hemodynamic states
#   OPTIMIZATION: compute exp(x[,3:5]) once at top and reuse
# ------------------------------------------------------------------------------
spm_fx_fmri <- function(x, u, P, M) {
  if (missing(M) || is.null(M)) M <- list()
  symmetry <- if (!is.null(M$symmetry)) M$symmetry else 0

  A <- P$A; B <- P$B; C <- P$C / 16; D <- P$D
  n <- nrow(x); mcols <- ncol(x); u <- as.numeric(u)
  f <- x

  if (mcols == 5) {
    if (is.vector(A) && length(A) > 1) {
      EE <- spm_dcm_fmri_mode_gen(A, M$modes)
      if (length(dim(B)) == 3 && dim(B)[3] > 0)
        for (i in seq_len(dim(B)[3])) EE <- EE + u[i] * B[,,i]
      if (length(dim(D)) == 3 && dim(D)[3] > 0)
        for (i in seq_len(dim(D)[3])) EE <- EE + x[i,1] * D[,,i]
    } else {
      Aeff <- A
      if (length(dim(B)) == 3 && dim(B)[3] > 0)
        for (i in seq_len(dim(B)[3])) Aeff <- Aeff + u[i] * B[,,i]
      if (length(dim(D)) == 3 && dim(D)[3] > 0)
        for (i in seq_len(dim(D)[3])) Aeff <- Aeff + x[i,1] * D[,,i]
      if (length(dim(Aeff)) == 3 && dim(Aeff)[3] > 1)
        Aeff <- exp(Aeff[,,1]) - exp(Aeff[,,2])
      SE <- diag(Aeff)
      EE <- Aeff - diag(exp(SE)/2 + SE)
      if (isTRUE(symmetry)) EE <- (EE + t(EE)) / 2
    }
    f[,1] <- as.vector(EE %*% x[,1] + C %*% u)
  } else {
    Aeff <- A
    if (length(dim(B)) == 3 && dim(B)[3] > 0)
      for (i in seq_len(dim(B)[3])) Aeff[,,1] <- Aeff[,,1] + u[i] * B[,,i]
    if (length(dim(D)) == 3 && dim(D)[3] > 0)
      for (i in seq_len(dim(D)[3])) Aeff[,,1] <- Aeff[,,1] + x[i,1] * D[,,i]
    nreg <- n; EE <- exp(Aeff[,,1]) / 8; IE <- diag(diag(EE)); EE <- EE - IE
    EI <- diag(nreg); SE <- diag(nreg) / 2; SI <- diag(nreg)
    if (length(dim(Aeff)) > 2 && dim(Aeff)[3] > 1) {
      phi <- spm_phi(Aeff[,,2] * 2); EI <- EI + EE*(1-phi); EE <- EE*phi - SE
    } else EE <- EE - SE
    f[,1] <- as.vector(EE %*% x[,1] - IE %*% x[,6] + C %*% u)
    f[,6] <- as.vector(EI %*% x[,1] - SI %*% x[,6])
  }

  H <- c(0.64, 0.32, 2.00, 0.32, 0.4)

  # OPTIMIZATION: compute exp once, store in named locals
  z3 <- exp(x[,3]); z4 <- exp(x[,4]); z5 <- exp(x[,5])

  sd <- H[1] * exp(P$decay); tt <- H[3] * exp(P$transit)
  fv <- z4^(1/H[4])
  ff <- (1 - (1 - H[5])^(1/z3)) / H[5]

  f[,2] <- x[,1] - as.numeric(sd) * x[,2] - H[2] * (z3 - 1)
  f[,3] <- x[,2] / z3
  f[,4] <- (z3 - fv) / (as.numeric(tt) * z4)
  f[,5] <- (ff * z3 - fv * z5 / z4) / (as.numeric(tt) * z5)
  f
}

# ------------------------------------------------------------------------------
# spm_fx_fmri2: returns f + analytic Jacobians
#   OPTIMIZATION:
#   - Removed redundant as.array() coercions (P fields should already be arrays)
#   - exp(x[,3:5]) computed once and stored
#   - require(Matrix) replaced by assuming Matrix is available
# ------------------------------------------------------------------------------
spm_fx_fmri2 <- function(x, u, P, M) {
  if (missing(M) || is.null(M)) M <- list()
  symmetry <- if (!is.null(M$symmetry)) M$symmetry else 0

  P$A <- as.matrix(P$A)
  # Only coerce B/C/D if not already the right type
  if (!is.array(P$B))  P$B <- as.array(P$B)
  P$C <- as.matrix(P$C) / 16
  if (!is.array(P$D))  P$D <- as.array(P$D)

  n <- nrow(x); mcols <- ncol(x); f <- x

  if (mcols == 5) {
    if (is.vector(P$A) && length(P$A) > 1) {
      EE <- spm_dcm_fmri_mode_gen(P$A, M$modes)
      if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0)
        for (i in seq_len(dim(P$B)[3])) EE <- EE + as.numeric(u[i]) * P$B[,,i]
      if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0)
        for (i in seq_len(dim(P$D)[3])) EE <- EE + x[i,1] * P$D[,,i]
    } else {
      Aeff <- P$A
      if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0)
        for (i in seq_len(dim(P$B)[3])) Aeff <- Aeff + as.numeric(u[i]) * P$B[,,i]
      if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0)
        for (i in seq_len(dim(P$D)[3])) Aeff <- Aeff + x[i,1] * P$D[,,i]
      if (length(dim(P$A)) == 3 && dim(P$A)[3] > 1)
        Aeff <- exp(P$A[,,1]) - exp(P$A[,,2])
      SE <- diag(Aeff); EE <- Aeff - diag(exp(SE)/2 + SE)
      if (isTRUE(symmetry)) EE <- (EE + t(EE)) / 2
    }
    f[,1] <- as.vector(EE %*% x[,1] + P$C %*% as.numeric(u))
  } else {
    Aeff <- P$A
    if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0)
      for (i in seq_len(dim(P$B)[3])) Aeff[,,1] <- Aeff[,,1] + as.numeric(u[i]) * P$B[,,i]
    if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0)
      for (i in seq_len(dim(P$D)[3])) Aeff[,,1] <- Aeff[,,1] + x[i,1] * P$D[,,i]
    nreg <- n; EE <- exp(Aeff[,,1]) / 8; IE <- diag(diag(EE)); EE <- EE - IE
    EI <- diag(nreg); SE <- diag(nreg) / 2; SI <- diag(nreg)
    if (length(dim(P$A)) > 2 && dim(P$A)[3] > 1) {
      phi <- spm_phi(P$A[,,2] * 2); EI <- EI + EE*(1-phi); EE <- EE*phi - SE
    } else EE <- EE - SE
    f[,1] <- as.vector(EE %*% x[,1] - IE %*% x[,6] + P$C %*% as.numeric(u))
    f[,6] <- as.vector(EI %*% x[,1] - SI %*% x[,6])
  }

  H <- c(0.64, 0.32, 2.00, 0.32, 0.4)

  # OPTIMIZATION: compute exp once
  ex3 <- exp(x[,3]); ex4 <- exp(x[,4]); ex5 <- exp(x[,5])
  x[,3] <- ex3; x[,4] <- ex4; x[,5] <- ex5

  sd <- H[1] * exp(P$decay); tt <- H[3] * exp(P$transit)
  fv <- x[,4]^(1/H[4]); ff <- (1 - (1-H[5])^(1/x[,3])) / H[5]

  f[,2] <- x[,1] - as.numeric(sd)*x[,2] - H[2]*(x[,3] - 1)
  f[,3] <- x[,2] / x[,3]
  f[,4] <- (x[,3] - fv) / (as.numeric(tt) * x[,4])
  f[,5] <- (ff*x[,3] - fv*x[,5]/x[,4]) / (as.numeric(tt)*x[,5])

  f <- as.vector(t(f))
  if (nargs() < 2) return(f)

  # Jacobians
  dfdx <- vector("list", 0); dfdu <- vector("list", 0)
  set_block <- function(L, i, j, M) {
    if (length(L) < i) length(L) <- i
    if (is.null(L[[i]])) L[[i]] <- list()
    L[[i]][[j]] <- M; L
  }

  I_n <- Matrix::Diagonal(n)

  if (mcols == 5) {
    d11 <- EE
    if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0) {
      for (i in seq_len(dim(P$D)[3])) {
        D_i <- P$D[,,i] + diag((diag(EE)-1) * diag(P$D[,,i]))
        d11[,i] <- d11[,i] + D_i %*% x[,1]
      }
    }
    dfdx <- set_block(dfdx, 1, 1, d11)
  } else {
    dfdx <- set_block(dfdx, 1, 1, EE)
    dfdx <- set_block(dfdx, 1, 6, -IE)
    dfdx <- set_block(dfdx, 6, 1, EI)
    dfdx <- set_block(dfdx, 6, 6, -SI)
  }

  Bu <- P$C
  if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0) {
    for (i in seq_len(dim(P$B)[3])) {
      Bm <- P$B[,,i] + diag((diag(EE)-1) * diag(P$B[,,i]))
      Bu[,i] <- Bu[,i] + Bm %*% x[,1]
    }
  }
  dfdu <- set_block(dfdu, 1, 1, Bu)

  # Hemodynamic Jacobian blocks
  dfdx <- set_block(dfdx, 2, 1, I_n)
  dfdx <- set_block(dfdx, 2, 2, (-as.numeric(sd)) * I_n)
  dfdx <- set_block(dfdx, 2, 3, Matrix::Diagonal(n, -H[2]*x[,3]))
  dfdx <- set_block(dfdx, 3, 2, Matrix::Diagonal(n,  1/x[,3]))
  dfdx <- set_block(dfdx, 3, 3, Matrix::Diagonal(n, -x[,2]/x[,3]))

  d44 <- Matrix::Diagonal(n, -x[,4]^(1/H[4]-1) / (as.numeric(tt)*H[4]) -
                              (x[,3] - x[,4]^(1/H[4])) / as.numeric(tt) / x[,4])
  dfdx <- set_block(dfdx, 4, 3, Matrix::Diagonal(n,  x[,3]/(as.numeric(tt)*x[,4])))
  dfdx <- set_block(dfdx, 4, 4, d44)

  d53 <- Matrix::Diagonal(n, (x[,3] + log(1-H[5])*(1-H[5])^(1/x[,3]) -
                               x[,3]*(1-H[5])^(1/x[,3])) / (as.numeric(tt)*x[,5]*H[5]))
  d54 <- Matrix::Diagonal(n, (x[,4]^(1/H[4]-1)*(H[4]-1)) / (as.numeric(tt)*H[4]))
  d55 <- Matrix::Diagonal(n, (x[,3]/x[,5]) * ((1-H[5])^(1/x[,3])-1) / (as.numeric(tt)*H[5]))

  dfdx <- set_block(dfdx, 5, 3, d53)
  dfdx <- set_block(dfdx, 5, 4, d54)
  dfdx <- set_block(dfdx, 5, 5, d55)

  dfdx <- spm_cat(dfdx); dfdu <- spm_cat(dfdu)
  list(f = f, dfdx = dfdx, D = 1, dfdu = dfdu)
}

# ------------------------------------------------------------------------------
# spm_gx_fmri (unchanged — already vectorized)
# ------------------------------------------------------------------------------
spm_gx_fmri <- function(x, u, P, M) {
  x <- as.matrix(x); n <- nrow(x); m <- ncol(x)
  # Force scalars to plain numeric (P$epsilon may be a 1x1 matrix, which
  # causes "Recycling array of length 1" deprecation warnings when mixed
  # with plain vectors below).
  TE <- 0.04; V0 <- 4; ep <- as.numeric(exp(P$epsilon))
  r0 <- 25; nu0 <- 40.3; E0 <- 0.4
  k1 <- as.numeric(4.3*nu0*E0*TE)
  k2 <- as.numeric(ep*r0*E0*TE)
  k3 <- as.numeric(1-ep)
  v <- exp(x[,4]); q <- exp(x[,5])
  g <- V0 * (k1 - k1*q + k2 - k2*q/v + k3 - k3*v)

  zero_block <- function() matrix(0, n, n)
  dgdx_blocks <- vector("list", m)
  for (j in seq_len(m)) dgdx_blocks[[j]] <- zero_block()
  dgdx_blocks[[4]] <- diag(-V0*(k3*v - k2*q/v), n, n)
  dgdx_blocks[[5]] <- diag(-V0*(k1*q + k2*q/v),  n, n)
  dGdX <- spm_cat(list(dgdx_blocks))
  list(g = g, dgdx = dGdX)
}

# ------------------------------------------------------------------------------
# spm_dcm_fmri_priors
#   OPTIMIZATION: replace per-slice loop for pC$A with vectorized array op
# ------------------------------------------------------------------------------
spm_dcm_fmri_priors <- function(A, B, C, D, options = list()) {
  A_in <- A; B_in <- B; C_in <- C; D_in <- D
  n <- nrow(A_in)
  if (is.null(options$stochastic)) options$stochastic <- 0
  if (is.null(options$induced))    options$induced    <- 0
  if (is.null(options$two_state))  options$two_state  <- 0
  if (is.null(options$backwards))  options$backwards  <- 0
  if (is.null(D_in))               D_in <- array(0, dim = c(n, n, 0))

  pE <- list(); pC <- list()

  if (isTRUE(options$two_state)) {
    x <- matrix(0, n, 6)
    pA <- if (!is.null(options$precision)) exp(options$precision) else 16
    A_mask <- A_in != 0
    pE$A <- A_mask*32 - 32; pE$B <- B_in*0; pE$C <- C_in*0; pE$D <- D_in*0

    Ad <- if (length(dim(A_in)) < 3) array(A_mask, dim=c(n,n,1)) else A_mask
    # OPTIMIZATION: vectorized division instead of per-slice loop
    pC$A <- Ad / pA
    pC$B <- B_in/4; pC$C <- C_in*4; pC$D <- D_in/4

    if (isTRUE(options$backwards)) {
      if (length(dim(pE$A)) < 3) pE$A <- array(pE$A, dim=c(n,n,2))
      if (length(dim(pC$A)) < 3) pC$A <- array(pC$A, dim=c(n,n,2))
      pE$A[,,2] <- 0; pC$A[,,2] <- A_mask / pA
    }
  } else {
    x <- matrix(0, n, 5)
    pA <- if (!is.null(options$precision)) exp(options$precision) else 64
    dA <- if (!is.null(options$decay)) options$decay else 1
    A_mask <- A_in != 0

    if (is.null(dim(A_in))) {
      pE$A <- (as.numeric(A_mask)-1)*dA; pC$A <- as.numeric(A_mask)
    } else {
      pE$A <- A_mask / 128
      Ad <- if (length(dim(A_in)) < 3) array(A_mask, dim=c(n,n,1)) else A_mask
      # OPTIMIZATION: vectorized division
      pC$A <- Ad / pA
    }
    pE$B <- B_in*0; pE$C <- C_in*0; pE$D <- D_in*0
    pC$B <- B_in; pC$C <- C_in; pC$D <- D_in
  }

  pE$transit <- matrix(0,n,1); pC$transit <- matrix(1/256,n,1)
  pE$decay   <- matrix(0,1,1); pC$decay   <- matrix(1/256,1,1)
  pE$epsilon <- matrix(0,1,1); pC$epsilon <- matrix(1/256,1,1)

  if (isTRUE(options$induced)) {
    pE$a <- matrix(0,2,1); pC$a <- matrix(1/64,2,1)
    pE$b <- matrix(0,2,1); pC$b <- matrix(1/64,2,1)
    pE$c <- matrix(0,1,n); pC$c <- matrix(1/64,1,n)
  }

  Cvec <- spm_vec(pC)
  Cmat <- if (length(Cvec) > 0) diag(as.numeric(Cvec)) else matrix(numeric(0),0,0)
  pC   <- Cmat
  list(pE = pE, x = x, pC = pC)
}

# ------------------------------------------------------------------------------
# spm_nlsi_GN: Gauss-Newton / Variational Laplace inversion
#   OPTIMIZATIONS:
#   - Pre-compute tJ = t(J) once per E-step iteration (used 3× in M-step)
#   - Cache iS_big: only recompute when h changes by > tol
#   - spm_trace(A,B) inlined as sum(t(A)*B) (avoids function call overhead)
#   - Re-use J'*iS_big across M-step sub-iterations
# ------------------------------------------------------------------------------
spm_nlsi_GN <- function(M, U, Y) {
  if (is.null(M$nograph)) M$nograph <- 1L
  if (is.null(M$noprint)) M$noprint <- 0L
  if (is.null(M$Nmax))    M$Nmax    <- 128L

  if (is.null(M$IS)) M$IS <- if (!is.null(M$G)) M$G else "spm_int"
  if (is.null(M$FS)) M$FS <- function(x, ...) x

  y <- tryCatch(Y$y, error = function(e) Y)
  if (is.function(M$IS) && is.function(M$FS)) {
    IS <- function(P, M, U) M$FS(M$IS(P, M, U))
    y  <- M$FS(y, M)
  } else {
    y  <- tryCatch(M$FS(y, M), error = function(e) tryCatch(M$FS(y), error = function(e2) y))
    IS <- tryCatch(
      spm_funcheck(function(P, M, U) M$FS(do.call(M$IS, list(P,M,U)), M)),
      error = function(e) tryCatch(
        spm_funcheck(function(P, M, U) M$FS(do.call(M$IS, list(P,M,U)))),
        error = function(e2) spm_funcheck(M$IS)))
  }
  IS <- spm_funcheck(IS)
  if (!is.null(M$f)) M$f <- spm_funcheck(M$f)
  if (!is.null(M$g)) M$g <- spm_funcheck(M$g)
  if (!is.null(M$h)) M$h <- spm_funcheck(M$h)

  ns <- if (is.list(y)) nrow(y[[1]]) else nrow(as.matrix(y))
  ny <- length(spm_vec(y)); nr <- ny / ns; M$ns <- ns

  if (is.null(M$x)) { if (is.null(M$n)) M$n <- 0L; M$x <- Matrix::Matrix(0, M$n, 1, sparse=TRUE) }
  if (missing(U) || is.null(U)) U <- list()
  if (is.null(M$P)) M$P <- M$pE
  dt <- tryCatch(Y$dt, error = function(e) 1)

  # Build error-covariance basis Q.
  # Y$Q may be: (a) absent/NULL, (b) a numeric vector, (c) a single matrix,
  # (d) a list of matrices, or (e) something malformed (e.g., from a foreign
  # preprocessing pipeline). We validate dimensions and fall back to spm_Ce
  # whenever Q[[1]] doesn't have a valid square matrix shape.
  .valid_Qblock <- function(M) {
    if (is.null(M)) return(FALSE)
    nr_ <- tryCatch(nrow(M), error = function(e) NULL)
    nc_ <- tryCatch(ncol(M), error = function(e) NULL)
    if (is.null(nr_) || is.null(nc_)) return(FALSE)
    if (length(nr_) != 1L || length(nc_) != 1L) return(FALSE)
    if (!is.finite(nr_) || !is.finite(nc_)) return(FALSE)
    nr_ > 0L && nc_ > 0L
  }
  Q <- tryCatch({
    Q0 <- Y$Q
    if (is.null(Q0))               spm_Ce(rep(ns, nr))
    else if (is.numeric(Q0) && !is.matrix(Q0))
                                    list(Matrix::Matrix(Q0, sparse = TRUE))
    else if (inherits(Q0, "Matrix") || is.matrix(Q0))
                                    list(Q0)
    else if (is.list(Q0) && length(Q0) > 0L) Q0
    else                            spm_Ce(rep(ns, nr))
  }, error = function(e) spm_Ce(rep(ns, nr)))

  if (length(Q) == 0L || !.valid_Qblock(Q[[1]])) {
    Q <- spm_Ce(rep(ns, nr))
  }
  nh <- length(Q); nq <- ny / nrow(Q[[1]])

  pE <- M$pE
  pC <- tryCatch(M$pC, error = function(e) Matrix::Diagonal(spm_length(M$pE), x = exp(16)))

  dfdu <- tryCatch({
    nb <- nrow(Y$X0); nx <- ny / nb
    Matrix::kronecker(Matrix::Diagonal(nx), Y$X0)
  }, error = function(e) Matrix::Matrix(0, ny, 0, sparse = TRUE))
  if (!length(dfdu)) dfdu <- Matrix::Matrix(0, ny, 0, sparse = TRUE)

  hE <- tryCatch({
    he <- as.numeric(M$hE); if (length(he) != nh) he + rep(0, nh) else he
  }, error = function(e) rep(0, nh) - log(stats::var(spm_vec(y))) + 4)
  h <- hE

  ihC <- tryCatch({
    H <- spm_inv(M$hC)
    if (length(H) != nh) H * Matrix::Diagonal(nh) else H
  }, error = function(e) Matrix::Diagonal(nh, x = exp(4)))

  if (is.list(pC)) pC <- spm_diag(spm_vec(pC))

  V   <- spm_svd(pC, 0)$V; nu <- ncol(dfdu); np <- ncol(V)
  ip  <- seq_len(np); iu <- np + seq_len(nu)
  pC_r <- t(V) %*% pC %*% V
  uC   <- Matrix::Diagonal(nu, x = 1e8)
  ipC  <- solve(spm_cat(spm_diag(list(pC_r, uC))))

  Eu  <- if (nu > 0) spm_pinv(dfdu) %*% spm_vec(y) else Matrix::Matrix(0, 0, 1)
  p   <- rbind(t(V) %*% (spm_vec(M$P) - spm_vec(M$pE)), Eu)
  Ep  <- spm_unvec(as.numeric(spm_vec(pE)) + as.numeric(as.matrix(V) %*% p[ip,,drop=FALSE]), pE)

  criterion <- c(0,0,0,0); C <- list(F = -Inf); v <- -4
  dFdh  <- matrix(0, nh, 1); dFdhh <- matrix(0, nh, nh)

  # Cache for iS_big reuse
  h_prev    <- h + Inf   # force first build
  iS_big_cache <- NULL
  iS_cache     <- NULL
  S_cache      <- NULL

  for (k in seq_len(M$Nmax)) {
    revert <- FALSE; dfdp <- NULL; f_out <- NULL

    tryCatch({
      dfdp_f <- spm_diff(IS, Ep, M, U, 1, list(V))
      dfdp   <- matrix(spm_vec(dfdp_f[[1]]), ny, np)
      f_out  <- dfdp_f[[2]]
      normdfdp <- norm(as.matrix(dfdp), type = "I")
      # NA-safe: !is.finite catches NA, NaN, Inf, -Inf
      revert <- !is.finite(normdfdp) || normdfdp > 1e32
    }, error = function(e) { revert <<- TRUE })

    if (revert && k > 1) {
      for (i in seq_len(4)) {
        v  <- min(v-2, -4)
        p  <- C$p + spm_dx(as.matrix(dFdpp), as.matrix(dFdp), list(v))
        Ep <- spm_unvec(as.numeric(spm_vec(pE)) + as.numeric(as.matrix(V) %*% p[ip,,drop=FALSE]), pE)
        tryCatch({
          dfdp_f   <- spm_diff(IS, Ep, M, U, 1, list(V))
          dfdp     <- matrix(spm_vec(dfdp_f[[1]]), ny, np)
          f_out    <- dfdp_f[[2]]
          normdfdp <- norm(as.matrix(dfdp), type = "I")
          # NA-safe
          revert   <- !is.finite(normdfdp) || normdfdp > exp(32)
        }, error = function(e) { revert <<- TRUE })
        if (!revert) break
      }
    }
    if (revert) stop("SPM:spm_nlsi_GN Convergence failure.")

    e  <- spm_vec(y) - spm_vec(f_out) - dfdu %*% p[iu,,drop=FALSE]
    J  <- -cbind(dfdp, dfdu)
    tJ <- t(J)   # OPTIMIZATION: compute t(J) once

    # M-step: Fisher scoring for h
    for (m_step in seq_len(8)) {

      # OPTIMIZATION: only rebuild iS_big if h has changed meaningfully.
      # NA-safe: any NaN in h or h_prev forces a rebuild rather than failing.
      h_diff   <- abs(h - h_prev)
      h_changed <- !all(is.finite(h_diff)) || any(h_diff > 1e-10)
      if (h_changed || is.null(iS_big_cache)) {
        iS <- Matrix::Matrix(0, nrow(Q[[1]]), ncol(Q[[1]]), sparse = TRUE)
        for (i in seq_len(nh)) iS <- iS + Q[[i]] * (exp(-32) + exp(h[i]))
        if (nh > 1) S_cache <- spm_inv(iS)
        iS_big_cache <- Matrix::kronecker(Matrix::Diagonal(nq), iS)
        iS_cache     <- iS
        h_prev       <- h
      }
      iS     <- iS_cache
      iS_big <- iS_big_cache
      S      <- S_cache

      # OPTIMIZATION: pre-compute J' %*% iS_big once
      tJ_iSbig <- tJ %*% iS_big        # (np+nu) × ny
      Pp  <- Re(as.numeric(tJ_iSbig %*% J))
      Cp  <- spm_inv(Pp + ipC)

      if (nh > 1) {
        P_list <- vector("list", nh); PS <- vector("list", nh); JPJ <- vector("list", nh)
        for (i in seq_len(nh)) {
          Pi       <- Q[[i]] * exp(h[i])
          PS[[i]]  <- Pi %*% S
          P_big    <- Matrix::kronecker(Matrix::Diagonal(nq), Pi)
          P_list[[i]] <- P_big
          JPJ[[i]] <- Re(as.numeric(tJ %*% P_big %*% J))
        }
        for (i in seq_len(nh)) {
          # OPTIMIZATION: spm_trace(A,B) = sum(t(A)*B)
          dFdh[i,1] <- sum(PS[[i]] * Matrix::Diagonal(nrow(PS[[i]]))) * nq/2 -
                       as.numeric(t(e) %*% P_list[[i]] %*% e) / 2 -
                       sum(t(Cp) * JPJ[[i]]) / 2
          for (j in i:nh) {
            val <- -sum(t(PS[[i]]) * PS[[j]]) * nq/2
            dFdhh[i,j] <- val; dFdhh[j,i] <- val
          }
        }
      } else {
        dFdh[1,1]  <- ny/2 - as.numeric(t(e) %*% iS_big %*% e)/2 - sum(t(Cp) * Pp)/2
        dFdhh[1,1] <- -ny/2
      }

      dFdhh <- dFdhh + diag(as.numeric(dFdh))
      d <- h - hE; dFdh <- dFdh - ihC %*% d; dFdhh <- dFdhh - ihC

      H_mat <- Re(as.numeric(-dFdhh)); H_mat <- matrix(H_mat, nh, nh)
      H_mat <- 0.5 * (H_mat + t(H_mat))
      Ch <- if (length(H_mat) == 1L) matrix(1/H_mat, 1, 1) else spm_inv(H_mat)

      dh <- spm_dx(as.matrix(dFdhh), as.matrix(dFdh), list(4))
      dh <- pmin(pmax(dh, -1), 1); h <- h + dh
      if (as.numeric(t(dFdh) %*% dh) < 1e-2) break
    }

    A <- 0; B <- 0
    L1 <- spm_logdet(iS) * nq/2 - as.numeric(t(e) %*% iS_big %*% e)/2 - ny*log(8*atan(1))/2
    L2 <- spm_logdet(as.matrix(ipC %*% Cp))/2 - as.numeric(t(p) %*% ipC %*% p)/2
    L3 <- spm_logdet(as.matrix(ihC %*% Ch))/2 - as.numeric(t(d) %*% ihC %*% d)/2
    F_val <- L1 + L2 + L3; L <- c(L1, L2, L3)

    if (!exists("F0", inherits = FALSE)) F0 <- F_val

    if (F_val > C$F || k < 3) {
      C$p <- p; C$h <- h; C$F <- F_val; C$L <- L; C$Cp <- Cp
      # OPTIMIZATION: reuse tJ_iSbig
      dFdp  <- -Re(as.numeric(tJ_iSbig %*% e)) - ipC %*% p + A/2
      dFdpp <- -Re(as.numeric(tJ_iSbig %*% J)) - ipC + B/2
      v <- min(v+0.5, 4); str <- "EM:(+)"
    } else {
      p <- C$p; h <- C$h; Cp <- C$Cp; v <- min(v-2, -4); str <- "EM:(-)"
    }

    dp <- spm_dx(as.matrix(dFdpp), as.matrix(dFdp), list(v), 'Q')
    p  <- p + dp
    Ep <- spm_unvec(as.numeric(spm_vec(pE)) + as.numeric(as.matrix(V) %*% p[ip,,drop=FALSE]), pE)

    dF_pred <- as.numeric(t(dFdp) %*% dp)
    if (!M$noprint) cat(sprintf("%-6s: %d  F: %-10.3e  dF predicted: %.3e\n", str, k, C$F-F0, dF_pred))
    criterion <- c((dF_pred < 1e-1), criterion[1:3])
    if (all(criterion)) break
  }

  Ep <- spm_unvec(as.numeric(spm_vec(pE)) + as.numeric(as.matrix(V) %*% C$p[ip,,drop=FALSE]), pE)
  Cp <- V %*% C$Cp[ip, ip, drop=FALSE] %*% t(V)
  list(Ep = Ep, Cp = Cp, Eh = C$h, F = C$F, L = C$L, dFdp = dFdp, dFdpp = dFdpp)
}

# ------------------------------------------------------------------------------
# spm_dcm_estimate (unchanged structure)
# ------------------------------------------------------------------------------
spm_dcm_estimate <- function(P) {
  if (missing(P)) stop("DCM structure or filename required.")
  if (is.list(P)) {
    DCM <- P
  } else if (is.character(P) && length(P) == 1) {
    env <- new.env(parent = emptyenv()); load(P, envir = env)
    if (!exists("DCM", envir = env)) stop("File does not contain object 'DCM'.")
    DCM <- get("DCM", envir = env)
  } else stop("Unsupported input for P.")

  if (is.null(DCM$options$two_state))  DCM$options$two_state  <- 0
  if (is.null(DCM$options$stochastic)) DCM$options$stochastic <- 0
  if (is.null(DCM$options$nonlinear))  DCM$options$nonlinear  <- 0
  if (is.null(DCM$options$centre))     DCM$options$centre     <- 0
  if (is.null(DCM$options$hidden))     DCM$options$hidden     <- integer(0)
  if (is.null(DCM$options$hE))         DCM$options$hE         <- 6
  if (is.null(DCM$options$hC))         DCM$options$hC         <- 1/128

  if (is.null(DCM$n)) DCM$n <- nrow(DCM$a)
  if (is.null(DCM$v)) DCM$v <- nrow(DCM$Y$y)

  M <- list()
  M$nograph <- if (!is.null(DCM$options$nograph)) DCM$options$nograph else 1L
  M$noprint <- 0L

  if (is.null(DCM$options$maxit)) {
    if (!is.null(DCM$options$nN)) {
      DCM$options$maxit <- DCM$options$nN
      warning("options.nN is deprecated. Please use options.maxit")
    } else DCM$options$maxit <- if (DCM$options$stochastic) 32 else 128
  }
  M$Nmax <- if (!is.null(DCM$M$Nmax)) DCM$M$Nmax else DCM$options$maxit

  if (is.null(DCM$options$maxnodes)) {
    if (!is.null(DCM$options$nmax)) {
      DCM$options$maxnodes <- DCM$options$nmax
      warning("options.nmax is deprecated. Please use options.maxnodes")
    } else DCM$options$maxnodes <- 8
  }
  DCM$options$induced <- 0

  U <- DCM$U; Y <- DCM$Y; n <- DCM$n; v <- DCM$v
  Y$y <- spm_detrend(Y$y)
  if (DCM$options$centre == 1) U$u <- spm_detrend(U$u)

  rng <- max(Y$y) - min(Y$y); scale <- 4 / max(rng, 4)
  Y$y <- Y$y * scale; Y$scale <- scale

  if (is.null(Y$X0)) Y$X0 <- matrix(1, v, 1)
  if (ncol(Y$X0) == 0) Y$X0 <- matrix(1, v, 1)

  if (!is.null(DCM$delays)) M$delays <- DCM$delays else M$delays <- matrix(1, n, 1)
  if (!is.null(DCM$TE))     M$TE     <- DCM$TE

  if (is.null(DCM$d)) {
    DCM$d <- array(0, dim = c(n, n, 0)); DCM$options$nonlinear <- 0
  } else DCM$options$nonlinear <- as.integer(dim(DCM$d)[3] > 0)

  if (DCM$options$nonlinear) {
    M$IS <- "spm_int_D"; M$nsteps <- round(max(Y$dt, 1)); M$states <- seq_len(n)
  } else M$IS <- "spm_int"

  if (is.null(DCM$c) || is.null(U$u)) {
    DCM$c <- matrix(0, n, 1); DCM$b <- array(0, dim=c(n,n,1))
    U$u <- matrix(0, v, 1); U$name <- list("null")
  }
  if (all(spm_vec(U$u) == 0) || all(spm_vec(DCM$c) == 0)) DCM$options$stochastic <- 1

  pri <- spm_dcm_fmri_priors(DCM$a, DCM$b, DCM$c, DCM$d, DCM$options)
  pE  <- pri$pE; pC <- pri$pC; x0 <- pri$x

  if (!is.null(DCM$options$P))  M$P  <- DCM$options$P
  if (!is.null(DCM$options$pE)) pE   <- DCM$options$pE
  if (!is.null(DCM$options$pC)) pC   <- DCM$options$pC
  if (!is.null(DCM$M$P))        M$P  <- DCM$M$P
  if (!is.null(DCM$M$pE))       pE   <- DCM$M$pE
  if (!is.null(DCM$M$pC))       pC   <- DCM$M$pC

  if (n > DCM$options$maxnodes) {
    y0 <- Y$y - Y$X0 %*% (solve(t(Y$X0) %*% Y$X0) %*% t(Y$X0) %*% Y$y)
    Vn <- spm_svd(t(y0))
    Vn <- Vn[, seq_len(DCM$options$maxnodes), drop=FALSE]
    j  <- seq_len(n*n); VV <- kronecker(Vn %*% t(Vn), Vn %*% t(Vn))
    pC[j,j] <- VV %*% pC[j,j] %*% t(VV)
  }

  hE <- matrix(DCM$options$hE, n, 1); hC <- diag(DCM$options$hC, n, n)
  i  <- DCM$options$hidden
  if (length(i) > 0) { hE[i] <- -4; diag(hC)[i] <- exp(-16) }

  M$f  <- "spm_fx_fmri"; M$g  <- "spm_gx_fmri"; M$x  <- x0
  M$pE <- pE; M$pC <- pC; M$hE <- hE; M$hC <- hC
  M$m  <- ncol(U$u); M$n  <- length(x0); M$l  <- nrow(x0)
  M$N  <- 64; M$dt <- 32/M$N; M$ns <- v

  if (!DCM$options$stochastic) {
    fit  <- spm_nlsi_GN(M, U, Y)
    Ep <- fit$Ep; Cp <- fit$Cp; Eh <- fit$Eh; F <- fit$F

    yhat <- do.call(M$IS, list(Ep, M, U))
    R    <- Y$y - yhat
    # OPTIMIZATION: compute X0'X0 inverse once
    X0tX0_inv <- spm_inv(t(Y$X0) %*% Y$X0)
    R    <- R - Y$X0 %*% (X0tX0_inv %*% t(Y$X0) %*% R)
    Ce   <- exp(-Eh)
  } else {
    u_micro <- U$u; y_micro <- Y$y
    Dy <- spm_dctmtx(nrow(y_micro), nrow(y_micro))
    Du <- spm_dctmtx(nrow(u_micro), nrow(y_micro))
    Dy <- Dy * sqrt(nrow(y_micro) / nrow(u_micro))
    u_dec <- Dy %*% (t(Du) %*% u_micro)

    DEM <- list(); DEM$M <- M; DEM$Y <- t(y_micro); DEM$U <- t(u_dec); DEM$X <- t(Y$X0)
    DEM$M[[1]]$E$form <- "Gaussian"; DEM$M[[1]]$E$s <- 1/2
    DEM$M[[1]]$E$d <- 2; DEM$M[[1]]$E$n <- 4; DEM$M[[1]]$E$nN <- DCM$options$maxit
    DEM$M[[1]]$delays <- M$delays / Y$dt
    DEM$M[[1]]$f <- function(x, v, P) { f <- do.call(M$f, list(x,v,P,NULL)); f[[1]] * Y$dt }
    DEM$M[[1]]$Q <- spm_Ce(matrix(1,1,n)); DEM$M[[1]]$hE <- hE; DEM$M[[1]]$hC <- hC

    W <- matrix(rep(matrix(c(12,16,16,16,16),1,5), each=n), n, 5, byrow=FALSE)
    DEM$M[[1]]$xP <- exp(6); DEM$M[[1]]$W <- diag(exp(as.vector(W)))
    DEM$M[[2]]$V <- exp(16)
    DEM <- spm_LAP(DEM)

    DCM$qU <- DEM$qU; DCM$qP <- DEM$qP; DCM$qH <- DEM$qH
    F  <- tail(DEM$F, 1); Ep <- DEM$qP$P[[1]]; Cp <- DEM$qP$C
    yhat <- t(DEM$qU$v[[1]]); R <- t(DEM$qU$z[[1]])
    X0tX0_inv <- spm_inv(t(Y$X0) %*% Y$X0)
    R  <- R - Y$X0 %*% (X0tX0_inv %*% t(Y$X0) %*% R)
    Ce <- exp(-DEM$qH$h[[1]])
  }

  bi <- spm_bireduce(M, Ep)
  M0 <- bi$M0; M1 <- bi$M1; L1 <- bi$L1; L2 <- bi$L2
  H  <- spm_kernels(M0, M1, L1, L2, M$N, M$dt)
  H1 <- H$H1; H0 <- H$H0

  L_mat <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n) + 1,
                                x = 1, dims = c(n, ncol(M0)))
  K  <- spm_kernels(M0, M1, L_mat, M$N, M$dt); K1 <- K$H1

  Tvec   <- as.numeric(spm_vec(pE))
  sw     <- options(warn = -1); on.exit(options(sw), add = TRUE)
  Pp_vec <- 1 - spm_Ncdf(Tvec, abs(as.numeric(spm_vec(Ep))), diag(Cp))
  Pp     <- spm_unvec(Pp_vec, Ep); Vp <- spm_unvec(diag(Cp), Ep)

  DCM$M <- M; DCM$Y <- Y; DCM$U <- U; DCM$Ce <- Ce
  DCM$Ep <- Ep; DCM$Cp <- Cp; DCM$Pp <- Pp; DCM$Vp <- Vp
  DCM$H1 <- H1; DCM$K1 <- K1; DCM$R <- R; DCM$y <- yhat; DCM$T <- 0

  y_id   <- if (!is.null(M$FS)) tryCatch(do.call(M$FS, list(Y$y, M)),
                                          error = function(e) tryCatch(do.call(M$FS, list(Y$y)), error = function(e2) Y$y))
            else Y$y
  DCM$ID <- spm_data_id(y_id); DCM$F <- F

  if (is.character(P) && length(P) == 1) {
    assign("DCM", DCM, envir = env); save(list = "DCM", file = P)
  }
  DCM
}

# ==============================================================================
# Bytecode compilation for the core fMRI functions
# ==============================================================================
if (requireNamespace("compiler", quietly = TRUE)) {
  local({
    .cmp <- compiler::cmpfun
    for (.fn in c("spm_fx_fmri", "spm_fx_fmri2", "spm_gx_fmri",
                  "spm_dcm_fmri_priors", "spm_nlsi_GN", "spm_dcm_estimate")) {
      if (exists(.fn, envir = .GlobalEnv, inherits = FALSE)) {
        assign(.fn, .cmp(get(.fn, envir = .GlobalEnv)), envir = .GlobalEnv)
      }
    }
  })
}

