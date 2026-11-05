# Run in each final analysis environment to record R and package versions.
packages <- c(
  "Rcpp", "Matrix", "MASS", "expm", "R.matlab",
  "snowfall", "snow", "Rmpi", "future", "future.apply"
)

versions <- vapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    "not installed"
  }
}, character(1))

output <- c(
  capture.output(sessionInfo()),
  "",
  "Selected package versions:",
  paste(names(versions), versions, sep = ": ")
)

writeLines(output, "sessionInfo.txt")
message("Wrote sessionInfo.txt")

