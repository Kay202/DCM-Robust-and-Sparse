# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Godfred Arhin and Nilotpal Sanyal
#
# Export the corrected within-subject MATLAB PEB files to the R object used by
# the group-level NARPS analysis. The original res2.rds is never overwritten.
#
# Run from any directory with:
#   Rscript export_corrected_pebs.R

if (!requireNamespace("R.matlab", quietly = TRUE)) {
  stop("Package 'R.matlab' is required. Install it with install.packages('R.matlab').")
}

args <- commandArgs(trailingOnly = TRUE)
revision_root <- Sys.getenv(
  "NARPS_REVISION_ROOT",
  unset = if (length(args)) args[1] else ""
)
if (!nzchar(revision_root)) {
  stop(
    "Supply the first-level revision directory as the first argument or set ",
    "NARPS_REVISION_ROOT."
  )
}

peb_dir <- file.path(revision_root, "04_within_subject_PEB")
output_dir <- file.path(revision_root, "07_exported_R_objects")
output_file <- file.path(output_dir, "res2_corrected.rds")

subject_numbers <- c(
   1,  2,  3,  4,  5,  6,  8,  9, 10, 11, 13, 14, 15, 16, 17,
  18, 19, 22, 24, 25, 26, 27, 29, 30, 32, 33, 35, 36, 37, 38,
  39, 40, 41, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 54, 55,
  56, 57, 58
)

subject_ids <- sprintf("sub-%03d", subject_numbers)
expected_files <- file.path(
  peb_dir,
  sprintf("PEB_corrected_%s_runwise.mat", subject_ids)
)

if (!dir.exists(peb_dir)) {
  stop("Corrected PEB directory was not found: ", peb_dir)
}

missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files)) {
  stop(
    "Missing corrected PEB files:\n",
    paste(missing_files, collapse = "\n")
  )
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (file.exists(output_file)) {
  stop(
    "The corrected export already exists and will not be overwritten: ",
    output_file,
    "\nMove or rename it explicitly before rerunning this script."
  )
}

# R.matlab represents MATLAB structs as list arrays. This normalization is the
# same one used by spm_dcm_peb_port.R.
normalize_mat_peb <- function(x) {
  if (!is.null(dim(x)) && length(dim(x)) >= 3L) {
    x <- x[, , 1]
  }
  if (!is.null(x$M) && !is.null(dim(x$M)) && length(dim(x$M)) >= 3L) {
    x$M <- x$M[, , 1]
  }
  x
}

read_corrected_peb <- function(file) {
  mat <- R.matlab::readMat(file, sparseMatrix = "matrix")
  if (is.null(mat$PEB)) {
    stop("Object 'PEB' was not found in: ", file)
  }

  peb <- normalize_mat_peb(mat$PEB)
  required <- c("Ep", "Cp", "Pnames")
  missing_fields <- setdiff(required, names(peb))
  if (length(missing_fields)) {
    stop(
      "Missing PEB field(s) ", paste(missing_fields, collapse = ", "),
      " in: ", file
    )
  }
  peb
}

extract_parameter_names <- function(peb) {
  parameter_names <- unlist(peb$Pnames, recursive = TRUE, use.names = FALSE)
  parameter_names <- trimws(as.character(parameter_names))
  parameter_names[!is.na(parameter_names) & nzchar(parameter_names)]
}

peb_list <- lapply(expected_files, read_corrected_peb)
names(peb_list) <- subject_ids

parameter_names <- extract_parameter_names(peb_list[[1]])
p <- length(as.numeric(peb_list[[1]]$Ep))

if (p != 22L) {
  stop("Expected 22 DCM parameters but found ", p, " in the first PEB file.")
}
if (length(parameter_names) != p || anyDuplicated(parameter_names)) {
  stop("The first PEB file does not contain 22 unique parameter names.")
}

eta_rows <- vector("list", length(peb_list))
Cp_list <- vector("list", length(peb_list))
Ce_list <- vector("list", length(peb_list))

for (i in seq_along(peb_list)) {
  peb <- peb_list[[i]]
  ep <- as.numeric(peb$Ep)
  cp_values <- as.numeric(peb$Cp)
  this_parameter_names <- extract_parameter_names(peb)

  if (length(ep) != p) {
    stop("Unexpected number of posterior means in: ", expected_files[i])
  }
  if (length(cp_values) != p * p) {
    stop("Unexpected posterior covariance dimension in: ", expected_files[i])
  }
  if (!identical(this_parameter_names, parameter_names)) {
    stop("Parameter names or ordering differ in: ", expected_files[i])
  }
  if (any(!is.finite(ep)) || any(!is.finite(cp_values))) {
    stop("Nonfinite posterior values found in: ", expected_files[i])
  }

  cp <- matrix(cp_values, nrow = p, ncol = p)
  asymmetry <- max(abs(cp - t(cp)))
  if (asymmetry > 1e-6) {
    stop(
      "Posterior covariance is not symmetric in ", expected_files[i],
      "; maximum asymmetry = ", signif(asymmetry, 4)
    )
  }

  dimnames(cp) <- list(parameter_names, parameter_names)
  eta_rows[[i]] <- ep
  Cp_list[[i]] <- cp
  Ce_list[i] <- list(if (is.null(peb$Ce)) NULL else peb$Ce)
}

eta <- do.call(rbind, eta_rows)
rownames(eta) <- subject_ids
colnames(eta) <- parameter_names
names(Cp_list) <- subject_ids
names(Ce_list) <- subject_ids

res2_corrected <- list(
  eta = eta,
  subject_names = sprintf("PEB_%s_runwise", subject_ids),
  parameter_names = parameter_names,
  file_type = "mat",
  PEB_list = peb_list,
  Cp_list = Cp_list,
  Ce_list = Ce_list,
  source_files = normalizePath(expected_files, mustWork = TRUE),
  export_information = list(
    source = "corrected run-level DCM and within-subject PEB analysis",
    number_of_subjects = nrow(eta),
    number_of_parameters = ncol(eta),
    created = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
)

saveRDS(res2_corrected, output_file, compress = "xz")

# Read the saved object once and verify the dimensions before reporting success.
check <- readRDS(output_file)
stopifnot(
  identical(dim(check$eta), c(48L, 22L)),
  length(check$Cp_list) == 48L,
  length(check$PEB_list) == 48L,
  length(check$subject_names) == 48L,
  length(check$parameter_names) == 22L
)

message("Created: ", output_file)
message("Validated 48 subjects and 22 parameters in the corrected export.")
