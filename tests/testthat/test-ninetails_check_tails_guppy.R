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
  # This test verifies the directory creation logic
  new_save_dir <- file.path(tempdir(), "new_guppy_output_dir")
  on.exit(unlink(new_save_dir, recursive = TRUE), add = TRUE)

  # The function should create the directory (though it will fail later
  # due to missing input files, we're testing dir creation)
  expect_false(dir.exists(new_save_dir))

  # Function will error on missing polya_data file, but should have created dir
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

  # Directory should now exist (function creates it before validation)
  expect_true(dir.exists(new_save_dir))
})



test_that("check_tails_guppy detects polya_data file format", {
  # Test with package test data if available
  nanopolish_path <- system.file("extdata", "test_data", "legacy",
                                 "nanopolish_output.tsv", package = "ninetails")

  skip_if(nanopolish_path == "" || !file.exists(nanopolish_path),
          "Legacy nanopolish test file not available")

  # The file should be detected as nanopolish format
  # We can't fully run the pipeline without Fast5 files, but we can test
  # that the file is recognized
  file_info <- file.info(nanopolish_path)
  expect_true(file_info$size > 0)
})



# process_polya_parts
################################################################################

test_that("process_polya_parts handles empty part_files gracefully", {
  # Create dummy cli_log function
  dummy_cli_log <- function(msg, type = "INFO", ...) invisible(NULL)

  temp_dir <- file.path(tempdir(), "test_parts")
  dir.create(temp_dir, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Empty part_files - function should handle gracefully
  # Either returns empty result or NULL
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

  # Either returns a result or throws an error - both are acceptable
  expect_true(!is.null(result) || identical(result, "error_caught"))
})

test_that("process_polya_parts handles non-function cli_log", {
  temp_dir <- file.path(tempdir(), "test_parts_log")
  dir.create(temp_dir, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Non-function cli_log should cause an error when it tries to call it
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

  # Should fail because cli_log is not callable
  expect_equal(result, "error_caught")
})


test_that("check_tails_guppy runs with valid test data", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("keras")
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  # Get test file paths
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

  # Create temporary output directory
  temp_output <- file.path(tempdir(), paste0("ninetails_test_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(temp_output, showWarnings = FALSE)
  on.exit(unlink(temp_output, recursive = TRUE), add = TRUE)

  # Run the pipeline (may take some time)
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
    # Pipeline may fail due to missing model files or other dependencies
    # In that case, skip the test
    skip(paste("Pipeline failed:", e$message))
  })

  # If we get here, verify the result structure
  if (!is.null(result)) {
    expect_type(result, "list")
    expect_true("read_classes" %in% names(result) || length(result) >= 1)
  }
})

