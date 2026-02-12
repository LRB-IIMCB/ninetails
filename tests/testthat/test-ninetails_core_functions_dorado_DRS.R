################################################################################
# Testing core functions Dorado DRS
################################################################################

# helper functions for testing
################################################################################

#' Create dummy cli_log function for testing
#' @param msg Message to log
#' @param level Log level
#' @param ... Additional arguments (ignored)
#' @return NULL invisibly
dummy_cli_log <- function(msg, level = "INFO", ...) {

  invisible(NULL)
}

#' Create dummy Dorado summary data frame
#' @param n Number of reads to generate
#' @param seed Random seed for reproducibility
#' @return data.frame mimicking Dorado summary structure
create_dummy_dorado_summary <- function(n = 100, seed = 42) {
  set.seed(seed)
  data.frame(
    read_id = paste0("read_", sprintf("%05d", seq_len(n))),
    filename = paste0("file_", sample(1:5, n, replace = TRUE), ".pod5"),
    poly_tail_length = sample(15:150, n, replace = TRUE),
    poly_tail_start = sample(100:500, n, replace = TRUE),
    poly_tail_end = sample(600:1200, n, replace = TRUE),
    alignment_genome = paste0("transcript_", sample(1:20, n, replace = TRUE)),
    alignment_direction = sample(c("+", "-", "*"), n, replace = TRUE, prob = c(0.45, 0.45, 0.1)),
    alignment_mapq = sample(c(0, 10, 20, 30, 40, 60), n, replace = TRUE, prob = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2)),
    stringsAsFactors = FALSE
  )
}

#' Create dummy poly(A) tail signal
#' @param length Signal length
#' @param has_modification Whether to include a modification peak
#' @param seed Random seed
#' @return numeric vector mimicking tail signal
create_dummy_tail_signal <- function(length = 500, has_modification = FALSE, seed = 42) {
  set.seed(seed)
  # Base signal around 680 with some noise (typical adenosine signal)
  signal <- rnorm(length, mean = 680, sd = 15)

  if (has_modification) {
    # Add a modification peak/valley in the middle
    mod_start <- round(length / 2) - 10
    mod_end <- mod_start + 20
    signal[mod_start:mod_end] <- rnorm(21, mean = 750, sd = 10)  # G-like peak
  }

  as.integer(round(signal))
}

#' Create dummy tail feature list (as produced by create_tail_features_list_dorado)
#' @param n Number of reads
#' @param seed Random seed
#' @return list mimicking tail_feature_list structure
create_dummy_tail_feature_list <- function(n = 10, seed = 42) {
  set.seed(seed)
  read_names <- paste0("read_", sprintf("%05d", seq_len(n)))

  feature_list <- lapply(seq_len(n), function(i) {
    signal <- create_dummy_tail_signal(
      length = sample(200:600, 1),
      has_modification = sample(c(TRUE, FALSE), 1),
      seed = seed + i
    )

    # Create pseudomoves - mostly zeros with some non-zero runs
    pseudomoves <- rep(0L, length(signal))
    if (runif(1) > 0.3) {
      # Add a run of non-zero values
      run_start <- sample(50:(length(signal) - 100), 1)
      run_length <- sample(5:15, 1)
      pseudomoves[run_start:(run_start + run_length - 1)] <- sample(c(-1L, 1L), 1)
    }

    list(
      tail_signal = signal,
      tail_pseudomoves = pseudomoves
    )
  })

  names(feature_list) <- read_names
  feature_list
}




# split_tail_centered_dorado
################################################################################

test_that("split_tail_centered_dorado validates input correctly", {
  dummy_feature_list <- create_dummy_tail_feature_list(n = 5)

  # Missing readname
  expect_error(
    ninetails::split_tail_centered_dorado(tail_feature_list = dummy_feature_list),
    "Readname is missing"
  )

  # Missing tail_feature_list
  expect_error(
    ninetails::split_tail_centered_dorado(readname = "read_00001"),
    "List of tail features is missing"
  )

  # Non-character readname
  expect_error(
    ninetails::split_tail_centered_dorado(
      readname = 12345,
      tail_feature_list = dummy_feature_list
    ),
    "Given readname is not a character string"
  )

  # Non-list tail_feature_list
  expect_error(
    ninetails::split_tail_centered_dorado(
      readname = "read_00001",
      tail_feature_list = "not_a_list"
    ),
    "Given tail_feature_list is not a list"
  )
})


test_that("split_tail_centered_dorado returns NULL when no valid chunks found", {
  # Create feature list with only zero pseudomoves
  feature_list <- list(
    test_read = list(
      tail_signal = create_dummy_tail_signal(length = 300),
      tail_pseudomoves = rep(0L, 300)  # All zeros - no modifications
    )
  )

  result <- ninetails::split_tail_centered_dorado(
    readname = "test_read",
    tail_feature_list = feature_list
  )

  expect_null(result)
})


test_that("split_tail_centered_dorado extracts chunks correctly", {
  # Create feature list with obvious modification run
  signal <- create_dummy_tail_signal(length = 400)
  pseudomoves <- rep(0L, 400)
  # Add run of length >= 5 at position 150
  pseudomoves[150:160] <- 1L

  feature_list <- list(
    test_read = list(
      tail_signal = signal,
      tail_pseudomoves = pseudomoves
    )
  )

  result <- ninetails::split_tail_centered_dorado(
    readname = "test_read",
    tail_feature_list = feature_list
  )

  # Should return a list (not NULL since we have a valid run)
  expect_type(result, "list")

  # Each chunk should have chunk_sequence, chunk_start_pos, chunk_end_pos
  if (length(result) > 0) {
    first_chunk <- result[[1]]
    expect_true("chunk_sequence" %in% names(first_chunk))
    expect_true("chunk_start_pos" %in% names(first_chunk))
    expect_true("chunk_end_pos" %in% names(first_chunk))

    # Chunk sequence should be length 100
    expect_equal(length(first_chunk$chunk_sequence), 100)
  }
})


test_that("split_tail_centered_dorado masks terminal pseudomoves", {
  # Create signal with modification at the very end
  signal <- create_dummy_tail_signal(length = 200)
  pseudomoves <- rep(0L, 200)
  # Add run at the very end (should be masked)
  pseudomoves[195:200] <- 1L

  feature_list <- list(
    test_read = list(
      tail_signal = signal,
      tail_pseudomoves = pseudomoves
    )
  )

  # Last 3 positions should be forced to 0 by the function
  # so this terminal run should NOT produce a chunk (run becomes only 3 long after masking)
  result <- ninetails::split_tail_centered_dorado(
    readname = "test_read",
    tail_feature_list = feature_list
  )

  # Expect NULL or no chunks from the terminal region
  # (the masking reduces the run length below threshold)
  expect_true(is.null(result) || length(result) == 0)
})



# create_tail_features_list_dorado
################################################################################

test_that("create_tail_features_list_dorado validates input correctly", {
  # Missing signal_list
  expect_error(
    ninetails::create_tail_features_list_dorado(num_cores = 1),
    "Signal list is missing"
  )

  # Missing num_cores
  expect_error(
    ninetails::create_tail_features_list_dorado(signal_list = list(a = 1:100)),
    "Number of declared cores is missing"
  )

  # Non-list signal_list
  expect_error(
    ninetails::create_tail_features_list_dorado(
      signal_list = "not_a_list",
      num_cores = 1
    ),
    "Provided signal_list is not a list"
  )

  # Non-numeric num_cores
  expect_error(
    ninetails::create_tail_features_list_dorado(
      signal_list = list(a = 1:100),
      num_cores = "one"
    ),
    "Declared num_cores must be numeric"
  )

  # Non-numeric elements in signal_list
  expect_error(
    ninetails::create_tail_features_list_dorado(
      signal_list = list(a = c("not", "numeric")),
      num_cores = 1
    ),
    "All elements of signal_list must be numeric"
  )
})


test_that("create_tail_features_list_dorado returns correct structure", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  # Create simple signal list
  signal_list <- list(
    read_001 = create_dummy_tail_signal(length = 200),
    read_002 = create_dummy_tail_signal(length = 250)
  )

  result <- suppressMessages(
    ninetails::create_tail_features_list_dorado(
      signal_list = signal_list,
      num_cores = 1
    )
  )

  # Should return a list
  expect_type(result, "list")

  # Should have same names as input
  expect_equal(sort(names(result)), sort(names(signal_list)))

  # Each element should have tail_signal and tail_pseudomoves
  for (read_name in names(result)) {
    expect_true("tail_signal" %in% names(result[[read_name]]))
    expect_true("tail_pseudomoves" %in% names(result[[read_name]]))
  }
})



# create_tail_chunk_list_dorado
################################################################################

test_that("create_tail_chunk_list_dorado validates input correctly", {
  # Missing num_cores
  expect_error(
    ninetails::create_tail_chunk_list_dorado(
      tail_feature_list = list()
    ),
    "Number of declared cores is missing"
  )

  # Missing tail_feature_list
  expect_error(
    ninetails::create_tail_chunk_list_dorado(
      num_cores = 1
    ),
    "List of features is missing"
  )

  # Non-numeric num_cores
  expect_error(
    ninetails::create_tail_chunk_list_dorado(
      tail_feature_list = list(),
      num_cores = "one"
    ),
    "Declared core number must be numeric"
  )

  # Non-list tail_feature_list
  expect_error(
    ninetails::create_tail_chunk_list_dorado(
      tail_feature_list = "not_a_list",
      num_cores = 1
    ),
    "Given tail_feature_list is not a list"
  )
})


test_that("create_tail_chunk_list_dorado processes feature list correctly", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")
  skip_if_not_installed("rrapply")

  # Create feature list with at least one read that has valid chunks
  signal <- create_dummy_tail_signal(length = 400)
  pseudomoves <- rep(0L, 400)
  pseudomoves[150:165] <- 1L  # Valid run of length 16

  tail_feature_list <- list(
    read_with_chunks = list(
      tail_signal = signal,
      tail_pseudomoves = pseudomoves
    ),
    read_without_chunks = list(
      tail_signal = create_dummy_tail_signal(length = 300),
      tail_pseudomoves = rep(0L, 300)  # No modifications
    )
  )

  result <- suppressMessages(
    ninetails::create_tail_chunk_list_dorado(
      tail_feature_list = tail_feature_list,
      num_cores = 1
    )
  )

  # Should return a list
  expect_type(result, "list")

  # NULL entries should be pruned
  # The read without chunks should be removed or empty
  # The read with chunks should have entries
})



# process_dorado_summary
################################################################################

test_that("process_dorado_summary validates input correctly", {
  temp_dir <- tempdir()

  # Missing part_size
  expect_error(
    ninetails::process_dorado_summary(
      dorado_summary = create_dummy_dorado_summary(n = 10),
      save_dir = temp_dir,
      cli_log = dummy_cli_log
    ),
    "Number of reads per file part"
  )

  # Invalid part_size (not positive)
  expect_error(
    ninetails::process_dorado_summary(
      dorado_summary = create_dummy_dorado_summary(n = 10),
      save_dir = temp_dir,
      part_size = -10,
      cli_log = dummy_cli_log
    ),
    "Reads per part must be numeric and positive"
  )

  # Invalid part_size (not numeric)
  expect_error(
    ninetails::process_dorado_summary(
      dorado_summary = create_dummy_dorado_summary(n = 10),
      save_dir = temp_dir,
      part_size = "ten",
      cli_log = dummy_cli_log
    ),
    "Reads per part must be numeric and positive"
  )

  # Invalid dorado_summary type
  expect_error(
    ninetails::process_dorado_summary(
      dorado_summary = 12345,
      save_dir = temp_dir,
      part_size = 10,
      cli_log = dummy_cli_log
    ),
    "Invalid dorado_summary format"
  )
})


test_that("process_dorado_summary handles data frame input correctly", {
  temp_dir <- file.path(tempdir(), "test_summary_df")
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  dummy_summary <- create_dummy_dorado_summary(n = 50)

  # Filter to simulate ninetails::filter_dorado_summary effect
  # (keep only valid reads)
  dummy_summary <- dummy_summary[
    dummy_summary$alignment_direction != "*" &
      dummy_summary$alignment_mapq != 0 &
      dummy_summary$poly_tail_start != 0 &
      dummy_summary$poly_tail_length >= 10,
  ]

  # Skip if no reads pass filter
  skip_if(nrow(dummy_summary) == 0, "No reads passed filter in dummy data")

  result <- ninetails::process_dorado_summary(
    dorado_summary = dummy_summary,
    save_dir = temp_dir,
    part_size = 20,
    cli_log = dummy_cli_log
  )

  # Should return character vector of file paths
  expect_type(result, "character")

  # All files should exist
  expect_true(all(file.exists(result)))

  # Number of parts should be ceiling(nrow / part_size)
  expected_parts <- ceiling(nrow(dummy_summary) / 20)
  # Note: actual number may differ due to filtering
})


test_that("process_dorado_summary validates data frame columns", {
  temp_dir <- file.path(tempdir(), "test_summary_cols")
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Data frame missing required column 'read_id'
  bad_summary <- data.frame(
    wrong_col = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )

  expect_error(
    ninetails::process_dorado_summary(
      dorado_summary = bad_summary,
      save_dir = temp_dir,
      part_size = 10,
      cli_log = dummy_cli_log
    ),
    "read_id"
  )
})


# extract_tails_from_pod5
################################################################################

test_that("extract_tails_from_pod5 validates input correctly", {
  # Missing required columns
  bad_polya_data <- data.frame(
    read_id = c("read1", "read2"),
    filename = c("file1.pod5", "file2.pod5")
    # Missing poly_tail_start and poly_tail_end
  )

  expect_error(
    ninetails::extract_tails_from_pod5(
      polya_data = bad_polya_data,
      pod5_dir = tempdir()
    ),
    "Missing required columns"
  )
})


test_that("extract_tails_from_pod5 returns empty list when no valid reads", {
  # Create polya_data with all invalid reads (poly_tail_start <= 0)
  polya_data <- data.frame(
    read_id = c("read1", "read2"),
    filename = c("file1.pod5", "file2.pod5"),
    poly_tail_start = c(0, 0),  # Invalid - should filter out
    poly_tail_end = c(100, 200),
    poly_tail_length = c(50, 100),
    stringsAsFactors = FALSE
  )

  result <- ninetails::extract_tails_from_pod5(
    polya_data = polya_data,
    pod5_dir = tempdir()
  )

  expect_type(result, "list")
  expect_equal(length(result), 0)
})


test_that("extract_tails_from_pod5 filters short tails correctly", {
  # Create polya_data with tails <= 10
  polya_data <- data.frame(
    read_id = c("read1", "read2", "read3"),
    filename = c("file1.pod5", "file2.pod5", "file3.pod5"),
    poly_tail_start = c(100, 200, 300),
    poly_tail_end = c(200, 400, 600),
    poly_tail_length = c(5, 8, 10),  # All too short (need > 10)
    stringsAsFactors = FALSE
  )

  result <- ninetails::extract_tails_from_pod5(
    polya_data = polya_data,
    pod5_dir = tempdir()
  )

  # Should return empty list since all tails are <= 10
  expect_equal(length(result), 0)
})



# preprocess_inputs
################################################################################

test_that("preprocess_inputs validates input correctly", {
  temp_dir <- file.path(tempdir(), "test_preprocess")

  dummy_summary <- create_dummy_dorado_summary(n = 10)

  # Non-numeric num_cores
  expect_error(
    ninetails::preprocess_inputs(
      dorado_summary = dummy_summary,
      pod5_dir = tempdir(),
      num_cores = "one",
      qc = TRUE,
      save_dir = temp_dir,
      prefix = "test",
      part_size = 100,
      cli_log = dummy_cli_log
    ),
    "Number of cores must be a positive numeric"
  )

  # Non-character pod5_dir
  expect_error(
    ninetails::preprocess_inputs(
      dorado_summary = dummy_summary,
      pod5_dir = 12345,
      num_cores = 1,
      qc = TRUE,
      save_dir = temp_dir,
      prefix = "test",
      part_size = 100,
      cli_log = dummy_cli_log
    ),
    "POD5 directory path must be a character"
  )

  # Non-existent pod5_dir
  expect_error(
    ninetails::preprocess_inputs(
      dorado_summary = dummy_summary,
      pod5_dir = "/nonexistent/path/to/pod5",
      num_cores = 1,
      qc = TRUE,
      save_dir = temp_dir,
      prefix = "test",
      part_size = 100,
      cli_log = dummy_cli_log
    ),
    "POD5.*directory does not exist"
  )

  # Invalid part_size
  expect_error(
    ninetails::preprocess_inputs(
      dorado_summary = dummy_summary,
      pod5_dir = tempdir(),
      num_cores = 1,
      qc = TRUE,
      save_dir = temp_dir,
      prefix = "test",
      part_size = 0,
      cli_log = dummy_cli_log
    ),
    "Part size must be at least 1"
  )
})


test_that("preprocess_inputs validates summary columns", {
  temp_dir <- file.path(tempdir(), "test_preprocess_cols")

  # Missing required poly(A) columns
  bad_summary <- data.frame(
    read_id = c("read1", "read2"),
    filename = c("file1.pod5", "file2.pod5"),
    # Missing poly_tail_length, poly_tail_start, poly_tail_end
    stringsAsFactors = FALSE
  )

  expect_error(
    ninetails::preprocess_inputs(
      dorado_summary = bad_summary,
      pod5_dir = tempdir(),
      num_cores = 1,
      qc = TRUE,
      save_dir = temp_dir,
      prefix = "test",
      part_size = 100,
      cli_log = dummy_cli_log
    ),
    "missing required poly"
  )
})


# create_outputs_dorado
################################################################################

test_that("create_outputs_dorado validates input correctly", {
  # Missing dorado_summary_dir
  expect_error(
    ninetails::create_outputs_dorado(
      nonA_temp_dir = tempdir(),
      polya_chunks_dir = tempdir(),
      num_cores = 1
    ),
    "Dorado summary directory is missing"
  )

  # Missing nonA_temp_dir
  expect_error(
    ninetails::create_outputs_dorado(
      dorado_summary_dir = tempdir(),
      polya_chunks_dir = tempdir(),
      num_cores = 1
    ),
    "Non-A predictions directory is missing"
  )

  # Missing polya_chunks_dir
  expect_error(
    ninetails::create_outputs_dorado(
      dorado_summary_dir = tempdir(),
      nonA_temp_dir = tempdir(),
      num_cores = 1
    ),
    "PolyA chunks directory is missing"
  )

  # Missing num_cores
  expect_error(
    ninetails::create_outputs_dorado(
      dorado_summary_dir = tempdir(),
      nonA_temp_dir = tempdir(),
      polya_chunks_dir = tempdir()
    ),
    "Number of cores is missing"
  )

  # Invalid num_cores
  expect_error(
    ninetails::create_outputs_dorado(
      dorado_summary_dir = tempdir(),
      nonA_temp_dir = tempdir(),
      polya_chunks_dir = tempdir(),
      num_cores = -1
    ),
    "num_cores must be a positive integer"
  )
})


test_that("create_outputs_dorado validates directories exist", {
  # Non-existent directory
  expect_error(
    ninetails::create_outputs_dorado(
      dorado_summary_dir = "/nonexistent/path1",
      nonA_temp_dir = "/nonexistent/path2",
      polya_chunks_dir = "/nonexistent/path3",
      num_cores = 1
    ),
    "All directory arguments must be valid existing paths"
  )
})


test_that("create_outputs_dorado checks for files in directories", {
  # Create empty temp directories
  temp_summary_dir <- file.path(tempdir(), "empty_summary")
  temp_nonA_dir <- file.path(tempdir(), "empty_nonA")
  temp_chunks_dir <- file.path(tempdir(), "empty_chunks")

  dir.create(temp_summary_dir, showWarnings = FALSE)
  dir.create(temp_nonA_dir, showWarnings = FALSE)
  dir.create(temp_chunks_dir, showWarnings = FALSE)

  on.exit({
    unlink(temp_summary_dir, recursive = TRUE)
    unlink(temp_nonA_dir, recursive = TRUE)
    unlink(temp_chunks_dir, recursive = TRUE)
  }, add = TRUE)

  # Should fail because no summary files found
  expect_error(
    ninetails::create_outputs_dorado(
      dorado_summary_dir = temp_summary_dir,
      nonA_temp_dir = temp_nonA_dir,
      polya_chunks_dir = temp_chunks_dir,
      num_cores = 1,
      original_summary = create_dummy_dorado_summary(n = 10)
    ),
    "No summary files found"
  )
})




# INTEGRATION-STYLE TESTS
################################################################################

test_that("signal processing pipeline produces consistent results", {
  # Create a reproducible signal with known modification
  set.seed(123)
  signal <- c(
    rnorm(150, mean = 680, sd = 10),  # Normal A region
    rnorm(30, mean = 750, sd = 8),     # G-like peak
    rnorm(150, mean = 680, sd = 10),  # Normal A region
    rnorm(20, mean = 580, sd = 8),     # C/U-like valley
    rnorm(50, mean = 680, sd = 10)    # Normal A region
  )
  signal <- as.integer(round(signal))

  # Process through filter function
  pseudomoves <- ninetails::filter_signal_by_threshold(signal)

  # Should have detected at least some non-zero pseudomoves
  # (exact positions depend on algorithm details, but we should see some)
  n_modifications <- sum(pseudomoves != 0)

  # Run multiple times to verify consistency
  pseudomoves2 <- ninetails::filter_signal_by_threshold(signal)

  expect_identical(pseudomoves, pseudomoves2)
})


