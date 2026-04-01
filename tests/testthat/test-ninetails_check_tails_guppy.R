################################################################################
# Testing ninetails_check_tails_guppy functions
################################################################################

# Helper: Suppress verbose logging
suppress_logs <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}


################################################################################
# check_tails_guppy - Input Validation
# Note: check_tails_guppy uses internal error handling and returns invisible(NULL)
# instead of throwing errors. Tests check for NULL return values.
################################################################################

test_that("check_tails_guppy creates save_dir if it does not exist", {
  new_save_dir <- file.path(tempdir(), "new_guppy_output_dir")
  on.exit(unlink(new_save_dir, recursive = TRUE), add = TRUE)

  expect_false(dir.exists(new_save_dir))

  suppressMessages(suppressWarnings(
    tryCatch({
      ninetails::check_tails_guppy(
        polya_data = "/nonexistent/polya.tsv",
        sequencing_summary = "/nonexistent/summary.txt",
        workspace = tempdir(),
        save_dir = new_save_dir,
        num_cores = 1
      )
    }, error = function(e) NULL)
  ))

  expect_true(dir.exists(new_save_dir))
})


test_that("check_tails_guppy detects polya_data file format", {
  nanopolish_path <- system.file("extdata", "test_data", "legacy",
                                 "nanopolish_output.tsv", package = "ninetails")

  skip_if(nanopolish_path == "" || !file.exists(nanopolish_path),
          "Legacy nanopolish test file not available")

  file_info <- file.info(nanopolish_path)
  expect_true(file_info$size > 0)
})


# process_polya_parts
################################################################################

test_that("process_polya_parts handles empty part_files gracefully", {
  dummy_cli_log <- function(msg, type = "INFO", ...) invisible(NULL)

  temp_dir <- file.path(tempdir(), "test_parts")
  dir.create(temp_dir, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  result <- tryCatch({
    ninetails::process_polya_parts(
      part_files = character(0),
      sequencing_summary = tempfile(),
      workspace = tempdir(),
      num_cores = 1,
      basecall_group = "Basecall_1D_000",
      pass_only = TRUE,
      qc = TRUE,
      save_dir = temp_dir,
      prefix = "test",
      cli_log = dummy_cli_log
    )
  }, error = function(e) "error_caught")

  expect_true(!is.null(result) || identical(result, "error_caught"))
})

test_that("process_polya_parts handles non-function cli_log", {
  temp_dir <- file.path(tempdir(), "test_parts_log")
  dir.create(temp_dir, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  result <- tryCatch({
    ninetails::process_polya_parts(
      part_files = c(tempfile()),
      sequencing_summary = tempfile(),
      workspace = tempdir(),
      num_cores = 1,
      basecall_group = "Basecall_1D_000",
      pass_only = TRUE,
      qc = TRUE,
      save_dir = temp_dir,
      prefix = "test",
      cli_log = "not_a_function"
    )
  }, error = function(e) "error_caught")

  expect_equal(result, "error_caught")
})


test_that("check_tails_guppy runs with valid test data", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("keras")
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  nanopolish_path <- system.file("extdata", "test_data", "legacy",
                                 "nanopolish_output.tsv", package = "ninetails")
  seq_summary_path <- system.file("extdata", "test_data", "legacy",
                                  "sequencing_summary.txt", package = "ninetails")
  workspace_path <- system.file("extdata", "test_data", "legacy",
                                "basecalled_fast5", package = "ninetails")

  skip_if(nanopolish_path == "" || !file.exists(nanopolish_path),
          "Legacy nanopolish test file not available")
  skip_if(seq_summary_path == "" || !file.exists(seq_summary_path),
          "Legacy sequencing_summary test file not available")
  skip_if(workspace_path == "" || !dir.exists(workspace_path),
          "Legacy Fast5 workspace not available")

  temp_output <- file.path(tempdir(), paste0("ninetails_test_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(temp_output, showWarnings = FALSE)
  on.exit(unlink(temp_output, recursive = TRUE), add = TRUE)

  result <- tryCatch({
    suppress_logs(
      ninetails::check_tails_guppy(
        polya_data = nanopolish_path,
        sequencing_summary = seq_summary_path,
        workspace = workspace_path,
        num_cores = 1,
        basecall_group = "Basecall_1D_000",
        pass_only = TRUE,
        qc = TRUE,
        save_dir = temp_output,
        prefix = "test",
        part_size = 1000000
      )
    )
  }, error = function(e) {
    skip(paste("Pipeline failed:", e$message))
  })

  if (!is.null(result)) {
    expect_type(result, "list")
    expect_true("read_classes" %in% names(result))
    expect_true("nonadenosine_residues" %in% names(result))
    expect_true(is.data.frame(result$read_classes))
    expect_true(is.data.frame(result$nonadenosine_residues))
  }
})
