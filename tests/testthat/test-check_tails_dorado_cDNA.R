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
# The full pipeline (both orientations, Lachesis and dorado input) is tested
# at the end of this file on the shipped cDNA test data.
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


################################################################################
# check_tails_dorado_cDNA - full pipeline on the shipped test data
#
# Requires Python with the pod5 and tensorflow modules; skipped otherwise.
# Expected tables were produced with the models shipped in
# inst/extdata/cnn_model (see inst/extdata/test_data/test_cDNA_*.tsv).
################################################################################

#' Skip unless the full cDNA pipeline can run
#' @keywords internal
skip_if_no_cdna_runtime <- function() {
  skip_if_not_installed("vroom")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("keras")
  skip_if(nchar(system.file("extdata", "test_data", "pod5_cDNA", package = "ninetails")) == 0,
          "cDNA test data not found in installed package")
  skip_if_not(reticulate::py_module_available("pod5"), "Python module pod5 not available")
  skip_if_not(reticulate::py_module_available("tensorflow"), "Python module tensorflow not available")
}

#' Compare pipeline output with the expected tables
#' @keywords internal
expect_cdna_outputs <- function(result, expected_prefix) {
  test_data_dir <- system.file("extdata", "test_data", package = "ninetails")
  expected_classes <- as.data.frame(vroom::vroom(
    file.path(test_data_dir, paste0(expected_prefix, "_read_classes.tsv")), show_col_types = FALSE))
  expected_residues <- as.data.frame(vroom::vroom(
    file.path(test_data_dir, paste0(expected_prefix, "_nonadenosine_residues.tsv")), show_col_types = FALSE))

  expect_type(result, "list")
  expect_true(all(c("read_classes", "nonadenosine_residues") %in% names(result)))

  # read classes
  classes <- result$read_classes[match(expected_classes$readname, result$read_classes$readname), ]
  expect_equal(nrow(result$read_classes), nrow(expected_classes))
  expect_equal(classes$class, expected_classes$class)
  expect_equal(classes$comments, expected_classes$comments)
  expect_equal(classes$tail_type, expected_classes$tail_type)

  # non-adenosine residues: same reads and residues; positions may differ by
  # one window step (about 2.5 nt) between platforms because of float ties
  residues <- result$nonadenosine_residues
  residues <- residues[order(residues$readname, residues$est_nonA_pos), ]
  expected_residues <- expected_residues[order(expected_residues$readname, expected_residues$est_nonA_pos), ]
  expect_equal(residues$readname, expected_residues$readname)
  expect_equal(residues$prediction, expected_residues$prediction)
  expect_true(all(abs(residues$est_nonA_pos - expected_residues$est_nonA_pos) <= 3))
}

test_that("check_tails_dorado_cDNA reproduces expected outputs with Lachesis input", {
  skip_if_no_cdna_runtime()
  test_dir <- system.file("extdata", "test_data", "pod5_cDNA", package = "ninetails")
  save_dir <- tempfile("cdna_lachesis_")

  result <- suppressMessages(check_tails_dorado_cDNA(
    pod5_dir = test_dir,
    lachesis = file.path(test_dir, "lachesis_polya.tsv"),
    num_cores = 1,
    qc = TRUE,
    save_dir = save_dir,
    part_size = 10,
    cleanup = TRUE
  ))

  expect_cdna_outputs(result, "test_cDNA")
  # output tables are written and intermediate files are removed
  expect_equal(length(list.files(save_dir, pattern = "_read_classes\\.tsv$")), 1L)
  expect_equal(length(list.files(save_dir, pattern = "_nonadenosine_residues\\.tsv$")), 1L)
  expect_equal(length(list.dirs(save_dir, recursive = FALSE)), 0L)
})

test_that("check_tails_dorado_cDNA reproduces expected outputs with dorado input", {
  skip_if_no_cdna_runtime()
  skip_if_not_installed("Rsamtools")
  test_dir <- system.file("extdata", "test_data", "pod5_cDNA", package = "ninetails")

  result <- suppressMessages(check_tails_dorado_cDNA(
    bam_file = file.path(test_dir, "example_cDNA_dorado.bam"),
    dorado_summary = file.path(test_dir, "dorado_summary.txt"),
    pod5_dir = test_dir,
    num_cores = 1,
    qc = TRUE,
    save_dir = tempfile("cdna_dorado_")
  ))

  expect_cdna_outputs(result, "test_cDNA_dorado")
})
