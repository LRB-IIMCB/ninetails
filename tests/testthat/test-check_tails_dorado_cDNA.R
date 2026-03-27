################################################################################
# Testing check_tails_dorado_cDNA functions
################################################################################

# Helpers
################################################################################

#' Minimal cli_log stub
#' @keywords internal
dummy_cli_log_cdna_wrap <- function(msg, level = "INFO", ...) {
  return(invisible(NULL))
}

#' Create a minimal valid Dorado cDNA summary data frame
#' @keywords internal
create_valid_cdna_summary <- function(n = 5) {
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
# check_tails_dorado_cDNA
#
# The outer tryCatch at line 188 catches all downstream errors and returns
# invisible(NULL). There are no missing() guards on the function's own
# arguments. The save_dir creation block (lines 123-130) runs BEFORE the
# tryCatch, so an unwritable/invalid path there would throw directly.
# All other invalid inputs reach the inner pipeline and are silently caught.
#
# Testable surface:
#   (a) save_dir is created as a side effect (before tryCatch)
#   (b) invisible(NULL) returned for any invalid input that the pipeline catches
#
# Not tested:
#   - polyT processing (temporary polyA model used — excluded per user request)
#   - full pipeline (requires BAM + POD5 + Keras)
################################################################################

test_that("check_tails_dorado_cDNA creates save_dir when it does not exist", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_mkdir_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  expect_false(dir.exists(save_dir))

  # save_dir is created before the tryCatch; pipeline fails at bam_file
  # validation inside the tryCatch and returns NULL, but save_dir persists
  check_tails_dorado_cDNA(
    bam_file = "/nonexistent/file.bam",
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_true(dir.exists(save_dir))
})

test_that("check_tails_dorado_cDNA returns invisible NULL when bam_file does not exist", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_bam_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  result <- check_tails_dorado_cDNA(
    bam_file = "/nonexistent/file.bam",
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA returns invisible NULL when bam_file is not a string", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_bamtype_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  result <- check_tails_dorado_cDNA(
    bam_file = 12345,
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA returns invisible NULL when pod5_dir does not exist", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_pod5_", Sys.getpid()))
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit({ unlink(save_dir, recursive = TRUE); unlink(tmp_bam) }, add = TRUE)

  result <- check_tails_dorado_cDNA(
    bam_file = tmp_bam,
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = "/nonexistent/pod5/path",
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA returns invisible NULL when dorado_summary is missing required columns", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_cols_", Sys.getpid()))
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit({ unlink(save_dir, recursive = TRUE); unlink(tmp_bam) }, add = TRUE)

  bad_summary <- data.frame(read_id = paste0("r", 1:3),
                            stringsAsFactors = FALSE)

  result <- check_tails_dorado_cDNA(
    bam_file = tmp_bam,
    dorado_summary = bad_summary,
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA returns invisible NULL when dorado_summary is an empty data frame", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_emptydf_", Sys.getpid()))
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit({ unlink(save_dir, recursive = TRUE); unlink(tmp_bam) }, add = TRUE)

  result <- check_tails_dorado_cDNA(
    bam_file = tmp_bam,
    dorado_summary = data.frame(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA returns invisible NULL when num_cores is non-numeric", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_cores_", Sys.getpid()))
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit({ unlink(save_dir, recursive = TRUE); unlink(tmp_bam) }, add = TRUE)

  result <- check_tails_dorado_cDNA(
    bam_file = tmp_bam,
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = "four",
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA returns invisible NULL when qc is not logical", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_qc_", Sys.getpid()))
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit({ unlink(save_dir, recursive = TRUE); unlink(tmp_bam) }, add = TRUE)

  result <- check_tails_dorado_cDNA(
    bam_file = tmp_bam,
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = "yes",
    save_dir = save_dir,
    part_size = 100
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA returns invisible NULL when part_size is invalid", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_partsize_", Sys.getpid()))
  tmp_bam <- tempfile(fileext = ".bam")
  file.create(tmp_bam)
  on.exit({ unlink(save_dir, recursive = TRUE); unlink(tmp_bam) }, add = TRUE)

  result <- check_tails_dorado_cDNA(
    bam_file = tmp_bam,
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 0
  )

  expect_null(result)
})

test_that("check_tails_dorado_cDNA writes a log file to save_dir", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_log_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  # Pipeline fails inside tryCatch but the log file is created before
  # the first cli_log call inside the tryCatch at line 196
  check_tails_dorado_cDNA(
    bam_file = "/nonexistent/file.bam",
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  log_files <- list.files(save_dir, pattern = "_ninetails_cDNA\\.log$")
  expect_true(length(log_files) > 0)
})

test_that("check_tails_dorado_cDNA log file contains expected header content", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_logcontent_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  check_tails_dorado_cDNA(
    bam_file = "/nonexistent/file.bam",
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 100
  )

  log_files <- list.files(save_dir, pattern = "_ninetails_cDNA\\.log$",
                          full.names = TRUE)
  if (length(log_files) > 0) {
    log_content <- readLines(log_files[1])
    expect_true(any(grepl("Ninetails", log_content)))
  } else {
    skip("Log file not created — skipping content check")
  }
})

test_that("check_tails_dorado_cDNA prefix is reflected in log file name", {
  save_dir <- file.path(tempdir(), paste0("cdna_test_prefix_", Sys.getpid()))
  on.exit(unlink(save_dir, recursive = TRUE), add = TRUE)

  check_tails_dorado_cDNA(
    bam_file = "/nonexistent/file.bam",
    dorado_summary = create_valid_cdna_summary(),
    pod5_dir = tempdir(),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    prefix = "myexp",
    part_size = 100
  )

  log_files <- list.files(save_dir, pattern = "_ninetails_cDNA\\.log$")
  if (length(log_files) > 0) {
    expect_true(any(grepl("myexp", log_files)))
  } else {
    skip("Log file not created — skipping prefix check")
  }
})
