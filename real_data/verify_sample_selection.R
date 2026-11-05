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

full <- read.delim(
  file.path(script_dir, "participants_full.tsv"),
  stringsAsFactors = FALSE
)
included <- read.csv(
  file.path(script_dir, "included_48_participants.csv"),
  stringsAsFactors = FALSE
)

eligible <- full[!full$participant_id %in% c("sub-020", "sub-021"), ]
expected <- eligible[seq_len(48L), c("participant_id", "group", "gender", "age")]
names(expected)[names(expected) == "group"] <- "task_version"

rownames(expected) <- NULL
rownames(included) <- NULL
if (!identical(included, expected)) {
  stop("The included-participant file does not match the documented selection rule.")
}

message("Verified: the 48 included participants match the documented selection rule.")

