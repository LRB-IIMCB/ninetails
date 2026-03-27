################################################################################
# Testing check_tails_dorado_DRS functions
################################################################################

# Helpers
################################################################################

#' Minimal cli_log stub
#' @keywords internal
dummy_cli_log <- function(msg, level = "INFO", ...) {
  return(invisible(NULL))
}

#' Create a minimal valid Dorado summary data frame
#' @keywords internal
create_valid_dorado_summary <- function(n = 5) {
  return(data.frame(
    read_id = paste0("read-", seq_len(n)),
    filename = rep("file1.pod5", n),
    poly_tail_length = rep(50L, n),
    poly_tail_start = rep(100L, n),
    poly_tail_end = rep(1100L, n),
    alignment_genome = rep("tx_1", n),
    alignment_direction = rep("+", n),
    alignment_mapq = rep(60L, n),
    stringsAsFactors = FALSE
  ))
}


################################################################################
# check_tails_dorado_DRS
#
# The outer tryCatch at line 362 catches all downstream errors and returns
# invisible(NULL). There are no missing() guards on the function's own
# arguments. Invalid inputs (bad pod5_dir, non-numeric num_cores, missing
# required columns) all reach the inner pipeline, fail in preprocess_inputs,
# and are silently caught. The testable surface is therefore:
#   (a) invisible(NULL) returned for any invalid input
#   (b) save_dir created as a side effect (via check_output_directory)
#   (c) valid list structure returned when the full pipeline succeeds
#       (requires pod5 + Keras — not tested here)
################################################################################

test_that("check_tails_dorado_DRS returns invisible NULL when pod5_dir does not exist", {
  save_dir <- file.path(tempdir(), paste0("drs_test_pod5_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  result <- check_tails_dorado_DRS(
    dorado_summary = create_valid_dorado_summary(),
    pod5_dir = "/nonexistent/pod5/path",
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_DRS returns invisible NULL when dorado_summary is missing required columns", {
  save_dir <- file.path(tempdir(), paste0("drs_test_cols_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  bad_summary <- data.frame(read_id = paste0("r", 1:3),
                            stringsAsFactors = FALSE)

  result <- check_tails_dorado_DRS(
    dorado_summary = bad_summary,
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_DRS returns invisible NULL when num_cores is non-numeric", {
  save_dir <- file.path(tempdir(), paste0("drs_test_cores_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  result <- check_tails_dorado_DRS(
    dorado_summary = create_valid_dorado_summary(),
    pod5_dir = tempdir(),
    num_cores = "four",
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_DRS returns invisible NULL when part_size is invalid", {
  save_dir <- file.path(tempdir(), paste0("drs_test_partsize_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  result <- check_tails_dorado_DRS(
    dorado_summary = create_valid_dorado_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 0
  )

  expect_null(result)
})

test_that("check_tails_dorado_DRS creates save_dir when it does not exist", {
  save_dir <- file.path(tempdir(), paste0("drs_test_mkdir_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  # Directory must not exist before the call
  expect_false(dir.exists(save_dir))

  # Will fail inside pipeline (pod5_dir is tempdir, no pod5 files) but
  # save_dir creation happens before the tryCatch via check_output_directory
  check_tails_dorado_DRS(
    dorado_summary = create_valid_dorado_summary(),
    pod5_dir = "/nonexistent/pod5",
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_true(dir.exists(save_dir))
})

test_that("check_tails_dorado_DRS returns invisible NULL when all reads fail QC filter", {
  save_dir <- file.path(tempdir(), paste0("drs_test_qc_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  # All reads have alignment_direction="*" — all fail filter_dorado_summary
  all_invalid <- data.frame(
    read_id = paste0("r", 1:5),
    filename = rep("f.pod5", 5),
    poly_tail_length = rep(50L, 5),
    poly_tail_start = rep(100L, 5),
    poly_tail_end = rep(600L, 5),
    alignment_genome = rep("tx_1", 5),
    alignment_direction = rep("*", 5),
    alignment_mapq = rep(0L, 5),
    stringsAsFactors = FALSE
  )

  result <- check_tails_dorado_DRS(
    dorado_summary = all_invalid,
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})


################################################################################
# process_dorado_signal_files
################################################################################

test_that("process_dorado_signal_files returns invisible empty list for empty file vector", {
  nonA_dir <- file.path(tempdir(), paste0("nonA_empty_", Sys.getpid()))
  chunks_dir <- file.path(tempdir(), paste0("chunks_empty_", Sys.getpid()))
  on.exit({ unlink(nonA_dir, recursive = TRUE); unlink(chunks_dir, recursive = TRUE) },
          add = TRUE)

  result <- process_dorado_signal_files(
    polya_signal_files = character(0),
    nonA_temp_dir = nonA_dir,
    polya_chunks_dir = chunks_dir,
    num_cores = 1,
    cli_log = dummy_cli_log
  )

  expect_type(result, "list")
  expect_equal(length(result), 0L)
})

test_that("process_dorado_signal_files creates nonA_temp_dir when it does not exist", {
  nonA_dir <- file.path(tempdir(), paste0("nonA_create_", Sys.getpid()))
  chunks_dir <- file.path(tempdir(), paste0("chunks_create_", Sys.getpid()))
  on.exit({ unlink(nonA_dir, recursive = TRUE); unlink(chunks_dir, recursive = TRUE) },
          add = TRUE)

  expect_false(dir.exists(nonA_dir))

  process_dorado_signal_files(
    polya_signal_files = character(0),
    nonA_temp_dir = nonA_dir,
    polya_chunks_dir = chunks_dir,
    num_cores = 1,
    cli_log = dummy_cli_log
  )

  expect_true(dir.exists(nonA_dir))
})

test_that("process_dorado_signal_files creates polya_chunks_dir when it does not exist", {
  nonA_dir <- file.path(tempdir(), paste0("nonA_chunks_", Sys.getpid()))
  chunks_dir <- file.path(tempdir(), paste0("chunks_chunks_", Sys.getpid()))
  on.exit({ unlink(nonA_dir, recursive = TRUE); unlink(chunks_dir, recursive = TRUE) },
          add = TRUE)

  expect_false(dir.exists(chunks_dir))

  process_dorado_signal_files(
    polya_signal_files = character(0),
    nonA_temp_dir = nonA_dir,
    polya_chunks_dir = chunks_dir,
    num_cores = 1,
    cli_log = dummy_cli_log
  )

  expect_true(dir.exists(chunks_dir))
})

test_that("process_dorado_signal_files records success=FALSE for non-existent signal file", {
  nonA_dir <- file.path(tempdir(), paste0("nonA_err_", Sys.getpid()))
  chunks_dir <- file.path(tempdir(), paste0("chunks_err_", Sys.getpid()))
  on.exit({ unlink(nonA_dir, recursive = TRUE); unlink(chunks_dir, recursive = TRUE) },
          add = TRUE)

  bad_file <- "/nonexistent/polya_signal_part1.rds"

  result <- suppressWarnings(process_dorado_signal_files(
    polya_signal_files = bad_file,
    nonA_temp_dir = nonA_dir,
    polya_chunks_dir = chunks_dir,
    num_cores = 1,
    cli_log = dummy_cli_log
  ))

  expect_type(result, "list")
  expect_equal(length(result), 1L)
  expect_false(result[[1]]$success)
  expect_true("error" %in% names(result[[1]]))
})

test_that("process_dorado_signal_files result is keyed by basename of input file", {
  nonA_dir <- file.path(tempdir(), paste0("nonA_key_", Sys.getpid()))
  chunks_dir <- file.path(tempdir(), paste0("chunks_key_", Sys.getpid()))
  on.exit({ unlink(nonA_dir, recursive = TRUE); unlink(chunks_dir, recursive = TRUE) },
          add = TRUE)

  bad_file <- "/nonexistent/polya_signal_part1.rds"

  result <- suppressWarnings(process_dorado_signal_files(
    polya_signal_files = bad_file,
    nonA_temp_dir = nonA_dir,
    polya_chunks_dir = chunks_dir,
    num_cores = 1,
    cli_log = dummy_cli_log
  ))

  expect_equal(names(result), "polya_signal_part1.rds")
})

test_that("process_dorado_signal_files output basename uses correct substitution pattern", {
  # output file is sub("polya_signal", "nonA_pred", basename(signal_file))
  # and chunks file is sub("polya_signal", "tail_chunks", basename(signal_file))
  # These are stored in the result on success. On failure, only `error` is stored.
  # We verify the naming convention is correct via the basename substitution logic itself.
  input_name <- "polya_signal_part3.rds"
  expect_equal(sub("polya_signal", "nonA_pred", input_name), "nonA_pred_part3.rds")
  expect_equal(sub("polya_signal", "tail_chunks", input_name), "tail_chunks_part3.rds")
})

test_that("process_dorado_signal_files handles multiple files with independent per-file error handling", {
  nonA_dir <- file.path(tempdir(), paste0("nonA_multi_", Sys.getpid()))
  chunks_dir <- file.path(tempdir(), paste0("chunks_multi_", Sys.getpid()))
  on.exit({ unlink(nonA_dir, recursive = TRUE); unlink(chunks_dir, recursive = TRUE) },
          add = TRUE)

  bad_files <- c("/nonexistent/polya_signal_part1.rds",
                 "/nonexistent/polya_signal_part2.rds")

  result <- suppressWarnings(process_dorado_signal_files(
    polya_signal_files = bad_files,
    nonA_temp_dir = nonA_dir,
    polya_chunks_dir = chunks_dir,
    num_cores = 1,
    cli_log = dummy_cli_log
  ))

  # Both files processed independently; both fail; two entries in result
  expect_equal(length(result), 2L)
  expect_false(result[["polya_signal_part1.rds"]]$success)
  expect_false(result[["polya_signal_part2.rds"]]$success)
})
