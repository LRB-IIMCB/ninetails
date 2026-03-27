################################################################################
# Testing core functions Dorado DRS
################################################################################

# Helpers
################################################################################

#' Create dummy cli_log function for testing
#' @param msg Message to log
#' @param level Log level
#' @param ... Additional arguments (ignored)
#' @return NULL invisibly
#' @keywords internal
dummy_cli_log <- function(msg, level = "INFO", ...) {
  return(invisible(NULL))
}


#' Create dummy Dorado summary data frame
#' @param n Number of reads to generate
#' @param seed Random seed for reproducibility
#' @param all_invalid Logical; if TRUE all reads have alignment_direction = "*"
#'   and alignment_mapq = 0, so all fail filter_dorado_summary
#' @return data.frame mimicking Dorado summary structure
#' @keywords internal
create_dummy_dorado_summary <- function(n = 100, seed = 42, all_invalid = FALSE) {
  set.seed(seed)
  if (all_invalid) {
    direction <- rep("*", n)
    mapq <- rep(0L, n)
  } else {
    direction <- sample(c("+", "-", "*"), n, replace = TRUE, prob = c(0.45, 0.45, 0.1))
    mapq <- sample(c(0L, 10L, 20L, 30L, 40L, 60L), n, replace = TRUE,
                   prob = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2))
  }
  return(data.frame(
    read_id = paste0("read_", sprintf("%05d", seq_len(n))),
    filename = paste0("file_", sample(1:5, n, replace = TRUE), ".pod5"),
    poly_tail_length = sample(15:150, n, replace = TRUE),
    poly_tail_start = sample(100:500, n, replace = TRUE),
    poly_tail_end = sample(600:1200, n, replace = TRUE),
    alignment_genome = paste0("transcript_", sample(1:20, n, replace = TRUE)),
    alignment_direction = direction,
    alignment_mapq = mapq,
    stringsAsFactors = FALSE
  ))
}


#' Create dummy poly(A) tail signal
#' @param length Signal length
#' @param has_modification Whether to include a modification peak
#' @param seed Random seed
#' @return numeric vector mimicking tail signal
#' @keywords internal
create_dummy_tail_signal <- function(length = 500, has_modification = FALSE, seed = 42) {
  set.seed(seed)
  signal <- rnorm(length, mean = 680, sd = 15)
  if (has_modification) {
    mod_start <- round(length / 2) - 10
    mod_end <- mod_start + 20
    signal[mod_start:mod_end] <- rnorm(21, mean = 750, sd = 10)
  }
  return(as.integer(round(signal)))
}


#' Create dummy tail feature list (as produced by create_tail_features_list_dorado)
#' @param n Number of reads
#' @param seed Random seed
#' @return list mimicking tail_feature_list structure
#' @keywords internal
create_dummy_tail_feature_list <- function(n = 10, seed = 42) {
  set.seed(seed)
  read_names <- paste0("read_", sprintf("%05d", seq_len(n)))
  feature_list <- lapply(seq_len(n), function(i) {
    signal <- create_dummy_tail_signal(
      length = sample(200:600, 1),
      has_modification = sample(c(TRUE, FALSE), 1),
      seed = seed + i
    )
    pseudomoves <- rep(0L, length(signal))
    if (runif(1) > 0.3) {
      run_start <- sample(50:(length(signal) - 100), 1)
      run_length <- sample(5:15, 1)
      pseudomoves[run_start:(run_start + run_length - 1)] <- sample(c(-1L, 1L), 1)
    }
    return(list(
      tail_signal = signal,
      tail_pseudomoves = pseudomoves
    ))
  })
  names(feature_list) <- read_names
  return(feature_list)
}


#' Create a complete set of temp-dir fixtures for create_outputs_dorado tests.
#'
#' Writes one summary TSV (filtered, QC-passed reads), three prediction RDS
#' files (one per supported format), one chunk RDS file, and returns an
#' original_summary data frame with both passing and failing reads.
#'
#' Layout:
#'   Decorated: read_00001 (prediction C, valid position)
#'   Blank: read_00002 (all-A predictions → MPU)
#'   UNM: read_00003 (alignment_direction = "*" & mapq = 0)
#'   IRL: read_00004 (poly_tail_length = 8 < 10)
#'   BAC: read_00005 (poly_tail_start = 0)
#'   MAU: read_00006 (none of the above)
#'
#' @param base_dir Parent temp directory (created by caller).
#' @return Named list: summary_dir, nonA_dir, chunks_dir, original_summary.
#' @keywords internal
create_outputs_dorado_fixtures <- function(base_dir) {
  summary_dir <- file.path(base_dir, "dorado_summary_dir")
  nonA_dir    <- file.path(base_dir, "nonA_temp_dir")
  chunks_dir  <- file.path(base_dir, "polya_chunks_dir")
  dir.create(summary_dir, recursive = TRUE)
  dir.create(nonA_dir,    recursive = TRUE)
  dir.create(chunks_dir,  recursive = TRUE)

  # Filtered summary (QC-passed reads only).
  # IMPORTANT: read IDs must NOT contain underscores because the pipeline
  # extracts read_id from chunkname via sub('_.*', '', chunkname), which
  # strips everything from the first underscore onward.
  # Hyphenated IDs (like real Dorado UUIDs) are safe.
  summary_df <- data.frame(filename = rep("file1.pod5", 2),
                           read_id = c("readA-00001", "readB-00002"),
                           poly_tail_length = c(100L, 100L),
                           poly_tail_start = c(100L, 100L),
                           poly_tail_end = c(1100L, 1100L),
                           alignment_genome = c("transcript_1", "transcript_1"),
                           alignment_direction = c("+", "+"),alignment_mapq = c(60L, 60L),
                           stringsAsFactors = FALSE)

  vroom::vroom_write(summary_df,
                     file.path(summary_dir, "summary_part1.txt"),
                     delim = "\t")

  # Prediction RDS — format 1: list(chunkname, prediction)
  # readA-00001_1: prediction 1 = C (non-A, will appear in nonadenosine_residues)
  # readB-00002_1: prediction 0 = A (blank read)
  saveRDS(list(chunkname = c("readA-00001_1", "readB-00002_1"),
               prediction = c(1L, 0L)),
          file.path(nonA_dir, "pred_format1.rds")
  )
  # Prediction RDS — format 2: named list
  # readB-00002_2: prediction 0 = A
  saveRDS(list("readB-00002_2" = 0L),
          file.path(nonA_dir, "pred_format2.rds")
  )
  # Prediction RDS — format 3: named vector
  # readB-00002_3: prediction 0 = A
  saveRDS(c("readB-00002_3" = 0L),
          file.path(nonA_dir, "pred_format3.rds")
  )

  # Chunk RDS — nested list keyed by read_id.
  # readA-00001: 1 chunk → chunkname readA-00001_1
  # readB-00002: 3 chunks → chunknames readB-00002_1, readB-00002_2, readB-00002_3
  # chunk_start_pos=50, chunk_end_pos=149 → centr_signal_pos = mean(50,149) = 99.5
  # signal_length = 0.2*(1100-100) = 200
  # est_nonA_pos = round(100 - (100*99.5/200), 0) = round(50.25, 0) = 50
  # 50 is within [2, 98] so qc=TRUE does not discard it
  chunk_list <- list(
    "readA-00001" = list(list(chunk_sequence = rep(680L, 100),
                              chunk_start_pos = 50L,
                              chunk_end_pos = 149L)
    ),
    "readB-00002" = list(list(chunk_sequence = rep(680L, 100),
                              chunk_start_pos = 50L,
                              chunk_end_pos = 149L),
                         list(chunk_sequence = rep(680L, 100),
                              chunk_start_pos = 200L,
                              chunk_end_pos = 299L),
                         list(chunk_sequence = rep(680L, 100),
                              chunk_start_pos = 350L,
                              chunk_end_pos = 449L)))

  saveRDS(chunk_list, file.path(chunks_dir, "chunks_part1.rds"))

  # Original (unfiltered) summary — all 6 reads representing every comment code:
  #   readA-00001: decorated (YAY)
  #   readB-00002: blank/MPU (all-A predictions)
  #   readC-00003: UNM (alignment_direction="*" & mapq=0)
  #   readD-00004: IRL (poly_tail_length=8 < 10)
  #   readE-00005: BAC (poly_tail_start=0)
  #   readF-00006: MAU (mapped, valid tail, no non-A detected)
  original_summary <- data.frame(
    filename = rep("file1.pod5", 6),
    read_id = c("readA-00001", "readB-00002", "readC-00003",
                "readD-00004", "readE-00005", "readF-00006"),
    poly_tail_length = c(100L, 100L, 5L, 8L, 50L, 60L),
    poly_tail_start = c(100L, 100L, 200L, 300L, 0L, 150L),
    poly_tail_end = c(1100L, 1100L, 700L, 800L, 550L, 1150L),
    alignment_genome = c("tx_1", "tx_1", "tx_2", "tx_2", "tx_3", "tx_3"),
    alignment_direction = c("+", "+", "*", "+", "+", "+"),
    alignment_mapq = c(60L, 60L, 0L, 60L, 60L, 60L),
    stringsAsFactors = FALSE
  )

  return(list(
    summary_dir      = summary_dir,
    nonA_dir         = nonA_dir,
    chunks_dir       = chunks_dir,
    original_summary = original_summary
  ))
}


################################################################################
# split_tail_centered_dorado
################################################################################

test_that("split_tail_centered_dorado validates input correctly", {
  dummy_feature_list <- create_dummy_tail_feature_list(n = 5)

  expect_error(
    split_tail_centered_dorado(tail_feature_list = dummy_feature_list),
    "Readname is missing")
  expect_error(
    split_tail_centered_dorado(readname = "read_00001"),
    "List of tail features is missing")
  expect_error(
    split_tail_centered_dorado(readname = 12345,
                               tail_feature_list = dummy_feature_list),
    "Given readname is not a character string")
  expect_error(
    split_tail_centered_dorado(readname = "read_00001",
                               tail_feature_list = "not_a_list"),
    "Given tail_feature_list is not a list")
})


test_that("split_tail_centered_dorado returns NULL when no valid chunks found", {
  feature_list <- list(test_read = list(tail_signal = create_dummy_tail_signal(length = 300),
                                        tail_pseudomoves = rep(0L, 300)))

  result <- split_tail_centered_dorado(readname = "test_read",
                                       tail_feature_list = feature_list)
  expect_null(result)
})


test_that("split_tail_centered_dorado extracts chunks correctly", {
  signal <- create_dummy_tail_signal(length = 400)
  pseudomoves <- rep(0L, 400)
  pseudomoves[150:160] <- 1L

  feature_list <- list(test_read = list(tail_signal = signal,
                                        tail_pseudomoves = pseudomoves))

  result <- split_tail_centered_dorado(readname = "test_read",
                                       tail_feature_list = feature_list)

  expect_type(result, "list")
  if (length(result) > 0) {
    first_chunk <- result[[1]]
    expect_true("chunk_sequence"  %in% names(first_chunk))
    expect_true("chunk_start_pos" %in% names(first_chunk))
    expect_true("chunk_end_pos"   %in% names(first_chunk))
    expect_equal(length(first_chunk$chunk_sequence), 100)
  }
})


test_that("split_tail_centered_dorado masks terminal pseudomoves", {
  signal <- create_dummy_tail_signal(length = 200)
  pseudomoves <- rep(0L, 200)
  pseudomoves[195:200] <- 1L  # terminal run — last 3 will be zeroed

  feature_list <- list(test_read = list(tail_signal = signal,
                                        tail_pseudomoves = pseudomoves))

  result <- split_tail_centered_dorado(readname = "test_read",
                                       tail_feature_list = feature_list)
  # After masking, run length drops below 5 → no valid chunks
  expect_true(is.null(result) || length(result) == 0)
})


test_that("split_tail_centered_dorado handles negative start_pos with NA padding", {
  # Run of 5 starting at position 1 forces start_pos = 1 + floor(5/2) - 50 = -47 ≤ 0.
  # Line 935 pads the chunk with NAs and then imputes from most-frequent values.
  signal <- create_dummy_tail_signal(length = 200)
  pseudomoves <- integer(200)
  pseudomoves[1:5] <- 1L  # run at very start, unaffected by terminal masking

  feature_list <- list(
    test_read = list(
      tail_signal      = signal,
      tail_pseudomoves = pseudomoves
    )
  )
  result <- split_tail_centered_dorado(
    readname = "test_read",
    tail_feature_list = feature_list)

  expect_type(result, "list")
  # If chunk was extracted the sequence must still be length 100 (imputed)
  if (length(result) > 0) {
    expect_equal(length(result[[1]]$chunk_sequence), 100)
    expect_false(any(is.na(result[[1]]$chunk_sequence)))
  }
})


################################################################################
# create_tail_features_list_dorado
################################################################################

test_that("create_tail_features_list_dorado validates input correctly", {
  expect_error(create_tail_features_list_dorado(num_cores = 1),"Signal list is missing")

  expect_error(create_tail_features_list_dorado(signal_list = list(a = 1:100)),"Number of declared cores is missing")

  expect_error(create_tail_features_list_dorado(signal_list = "not_a_list",
                                                num_cores = 1),"Provided signal_list is not a list")

  expect_error(create_tail_features_list_dorado(signal_list = list(a = 1:100),
                                                num_cores = "one"),
               "Declared num_cores must be numeric")

  expect_error(create_tail_features_list_dorado(
    signal_list = list(a = c("not", "numeric")),
    num_cores = 1),
    "All elements of signal_list must be numeric")
})


test_that("create_tail_features_list_dorado returns correct structure", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  signal_list <- list(read_001 = create_dummy_tail_signal(length = 200),
                      read_002 = create_dummy_tail_signal(length = 250))

  result <- suppressMessages(
    create_tail_features_list_dorado(signal_list = signal_list,num_cores = 1))

  expect_type(result, "list")
  expect_equal(sort(names(result)), sort(names(signal_list)))
  for (read_name in names(result)) {
    expect_true("tail_signal" %in% names(result[[read_name]]))
    expect_true("tail_pseudomoves" %in% names(result[[read_name]]))
  }
})


################################################################################
# create_tail_chunk_list_dorado
################################################################################

test_that("create_tail_chunk_list_dorado validates input correctly", {
  expect_error(create_tail_chunk_list_dorado(tail_feature_list = list()),
               "Number of declared cores is missing")

  expect_error(create_tail_chunk_list_dorado(num_cores = 1),"List of features is missing")

  expect_error(create_tail_chunk_list_dorado(tail_feature_list = list(),num_cores = "one"),
               "Declared core number must be numeric")

  expect_error(create_tail_chunk_list_dorado(tail_feature_list = "not_a_list",num_cores = 1),
               "Given tail_feature_list is not a list")
})


test_that("create_tail_chunk_list_dorado processes feature list correctly", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")
  skip_if_not_installed("rrapply")

  signal <- create_dummy_tail_signal(length = 400)
  pseudomoves <- rep(0L, 400)
  pseudomoves[150:165] <- 1L

  tail_feature_list <- list(
    read_with_chunks = list(
      tail_signal = signal,
      tail_pseudomoves = pseudomoves
    ),
    read_without_chunks = list(
      tail_signal = create_dummy_tail_signal(length = 300),
      tail_pseudomoves = rep(0L, 300)
    )
  )

  result <- suppressMessages(create_tail_chunk_list_dorado(tail_feature_list = tail_feature_list,num_cores = 1))

  expect_type(result, "list")
  # NULL entries (reads without chunks) must be pruned
  expect_false("read_without_chunks" %in% names(result))
})


################################################################################
# process_dorado_summary
################################################################################

test_that("process_dorado_summary validates input correctly", {
  temp_dir <- tempdir()

  expect_error(process_dorado_summary(dorado_summary = create_dummy_dorado_summary(n = 10),
                                      save_dir = temp_dir,
                                      cli_log = dummy_cli_log),
               "Number of reads per file part")

  expect_error(process_dorado_summary(dorado_summary = create_dummy_dorado_summary(n = 10),
                                      save_dir = temp_dir,
                                      part_size = -10,
                                      cli_log = dummy_cli_log),
               "Reads per part must be numeric and positive")

  expect_error(process_dorado_summary(dorado_summary = create_dummy_dorado_summary(n = 10),
                                      save_dir = temp_dir,
                                      part_size = "ten",
                                      cli_log = dummy_cli_log),
               "Reads per part must be numeric and positive")

  expect_error(process_dorado_summary(dorado_summary = 12345,
                                      save_dir = temp_dir,
                                      part_size = 10,
                                      cli_log = dummy_cli_log),
               "Invalid dorado_summary format")
})


test_that("process_dorado_summary handles data frame input correctly", {
  temp_dir <- file.path(tempdir(), paste0("test_summary_df_", Sys.getpid()))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Keep only reads that survive filter_dorado_summary
  dummy_summary <- create_dummy_dorado_summary(n = 50)
  dummy_summary <- dummy_summary[
    dummy_summary$alignment_direction != "*" &
      dummy_summary$alignment_mapq != 0 &
      dummy_summary$poly_tail_start != 0 &
      dummy_summary$poly_tail_length >= 10,
  ]
  skip_if(nrow(dummy_summary) == 0, "No reads passed filter in dummy data")

  result <- process_dorado_summary(dorado_summary = dummy_summary,
                                   save_dir = temp_dir,
                                   part_size = 20,
                                   cli_log = dummy_cli_log)

  expect_type(result, "character")
  expect_true(all(file.exists(result)))
  expect_equal(length(result), ceiling(nrow(dummy_summary) / 20))
})


test_that("process_dorado_summary validates data frame columns", {
  temp_dir <- file.path(tempdir(), paste0("test_summary_cols_", Sys.getpid()))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  bad_summary <- data.frame(wrong_col = c("a", "b", "c"),stringsAsFactors = FALSE)

  expect_error(process_dorado_summary(dorado_summary = bad_summary,
                                      save_dir = temp_dir,
                                      part_size = 10,
                                      cli_log = dummy_cli_log),"read_id")
})


test_that("process_dorado_summary reads valid TSV file path", {
  skip_if_not_installed("vroom")

  temp_dir <- file.path(tempdir(), paste0("test_summary_file_", Sys.getpid()))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Write a valid TSV that will survive filter_dorado_summary
  tsv_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(tsv_path), add = TRUE)

  valid_df <- data.frame(read_id = paste0("read_", 1:5),
                         filename = rep("file1.pod5", 5),
                         poly_tail_length = rep(50L, 5),
                         poly_tail_start = rep(100L, 5),
                         poly_tail_end = rep(600L, 5),
                         alignment_genome = rep("tx_1", 5),
                         alignment_direction = rep("+", 5),
                         alignment_mapq = rep(60L, 5),
                         stringsAsFactors = FALSE)

  vroom::vroom_write(valid_df, tsv_path, delim = "\t")

  result <- process_dorado_summary(dorado_summary = tsv_path,
                                   save_dir = temp_dir,
                                   part_size = 10,
                                   cli_log = dummy_cli_log)

  expect_type(result, "character")
  expect_true(all(file.exists(result)))
})


test_that("process_dorado_summary errors on TSV missing read_id column", {
  skip_if_not_installed("vroom")

  temp_dir <- tempdir()
  tsv_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(tsv_path), add = TRUE)

  # Two-column TSV so vroom can detect the tab delimiter, but without
  # read_id, which triggers the inner stop at lines 72-79.
  # A single-column file causes "Could not guess delimiter" before reaching
  # the read_id check.
  bad_df <- data.frame(wrong_col = c("a", "b"),
                       other_col = c("x", "y"),
                       stringsAsFactors = FALSE)

  vroom::vroom_write(bad_df, tsv_path, delim = "\t")

  expect_error(process_dorado_summary(dorado_summary = tsv_path,
                                      save_dir = temp_dir,
                                      part_size = 10,
                                      cli_log = dummy_cli_log), "read_id")
})


test_that("process_dorado_summary errors via tryCatch on non-existent file (lines 81-84)", {
  temp_dir <- tempdir()
  expect_error(process_dorado_summary(dorado_summary = "/nonexistent/path/summary.tsv",
                                      save_dir = temp_dir,
                                      part_size = 10,
                                      cli_log = dummy_cli_log))
})


################################################################################
# extract_tails_from_pod5
################################################################################

test_that("extract_tails_from_pod5 validates input correctly", {
  bad_polya_data <- data.frame(read_id  = c("read1", "read2"),
                               filename = c("file1.pod5", "file2.pod5"),
                               stringsAsFactors = FALSE)

  expect_error(extract_tails_from_pod5(polya_data = bad_polya_data, pod5_dir = tempdir()),
               "Missing required columns")
})


test_that("extract_tails_from_pod5 returns empty list when no valid reads", {
  polya_data <- data.frame(read_id= c("read1", "read2"),
                           filename= c("file1.pod5", "file2.pod5"),
                           poly_tail_start = c(0L, 0L),
                           poly_tail_end= c(100L, 200L),
                           poly_tail_length = c(50L, 100L),
                           stringsAsFactors = FALSE)

  result <- extract_tails_from_pod5(polya_data = polya_data, pod5_dir = tempdir())

  expect_type(result, "list")
  expect_equal(length(result), 0)
})


test_that("extract_tails_from_pod5 filters short tails correctly", {
  polya_data <- data.frame(
    read_id         = c("read1", "read2", "read3"),
    filename        = c("file1.pod5", "file2.pod5", "file3.pod5"),
    poly_tail_start = c(100L, 200L, 300L),
    poly_tail_end   = c(200L, 400L, 600L),
    poly_tail_length = c(5L, 8L, 10L),
    stringsAsFactors = FALSE
  )
  result <- extract_tails_from_pod5(polya_data = polya_data, pod5_dir = tempdir())
  expect_equal(length(result), 0)
})


test_that("extract_tails_from_pod5 filters by poly_tail_start when no poly_tail_length column", {
  # Without poly_tail_length column the else branch at line 275 runs:
  # polya_data <- polya_data[polya_data$poly_tail_start > 0, ]
  # All reads have poly_tail_start = 0 → filtered to empty → returns list()
  polya_data <- data.frame(read_id = c("read1", "read2"),
                           filename = c("file1.pod5", "file2.pod5"),
                           poly_tail_start = c(0L, 0L),
                           poly_tail_end = c(200L, 400L),
                           stringsAsFactors = FALSE)
  result <- extract_tails_from_pod5(polya_data = polya_data, pod5_dir = tempdir())
  expect_type(result, "list")
  expect_equal(length(result), 0)
})


# Lines 283-336 (Python script lookup and execution) require a working Python
# environment with pandas and pod5 installed and are not tested here.


################################################################################
# preprocess_inputs
################################################################################

test_that("preprocess_inputs validates input correctly", {
  temp_dir <- file.path(tempdir(), paste0("test_preprocess_", Sys.getpid()))
  dummy_summary <- create_dummy_dorado_summary(n = 10)

  expect_error(preprocess_inputs(dorado_summary = dummy_summary,
                                 pod5_dir = tempdir(),
                                 num_cores = "one",
                                 qc = TRUE,
                                 save_dir = temp_dir,
                                 prefix = "test",
                                 part_size = 100,
                                 cli_log = dummy_cli_log),
               "Number of cores must be a positive numeric")

  expect_error(preprocess_inputs(dorado_summary = dummy_summary,
                                 pod5_dir = 12345,
                                 num_cores = 1,
                                 qc = TRUE,
                                 save_dir = temp_dir,
                                 prefix = "test",
                                 part_size = 100,
                                 cli_log = dummy_cli_log),
               "POD5 directory path must be a character")

  expect_error(preprocess_inputs(dorado_summary = dummy_summary,
                                 pod5_dir = "/nonexistent/path/to/pod5",
                                 num_cores = 1,
                                 qc = TRUE,
                                 save_dir = temp_dir,
                                 prefix = "test",
                                 part_size = 100,
                                 cli_log = dummy_cli_log),
               "POD5.*directory does not exist")

  expect_error(preprocess_inputs(dorado_summary = dummy_summary,
                                 pod5_dir = tempdir(),
                                 num_cores = 1,
                                 qc = TRUE,
                                 save_dir = temp_dir,
                                 prefix = "test",
                                 part_size = 0,
                                 cli_log = dummy_cli_log),
               "Part size must be at least 1")
})


test_that("preprocess_inputs validates summary columns", {
  temp_dir <- file.path(tempdir(), paste0("test_preprocess_cols_", Sys.getpid()))

  bad_summary <- data.frame(read_id = c("read1", "read2"),
                            filename = c("file1.pod5", "file2.pod5"),
                            stringsAsFactors = FALSE)
  expect_error(
    preprocess_inputs(dorado_summary = bad_summary,
                      pod5_dir = tempdir(),
                      num_cores = 1,
                      qc = TRUE,
                      save_dir = temp_dir,
                      prefix = "test",
                      part_size = 100,
                      cli_log = dummy_cli_log),
    "missing required poly")
})


test_that("preprocess_inputs renames input_filename to filename (line 465)", {
  # Dorado >= 1.4.0 uses 'input_filename'. The function renames it to 'filename'
  # at line 464-465. If the rename fires, the subsequent required-columns error
  # will mention poly_tail_length/poly_tail_start but NOT filename.
  temp_dir <- file.path(tempdir(), paste0("test_rename_", Sys.getpid()))

  summary_new_format <- data.frame(
    read_id        = c("read1", "read2"),
    input_filename = c("file1.pod5", "file2.pod5"),  # Dorado >= 1.4.0 column name
    stringsAsFactors = FALSE
  )
  err <- tryCatch(preprocess_inputs(dorado_summary = summary_new_format,
                                    pod5_dir = tempdir(),
                                    num_cores = 1,
                                    qc = TRUE,
                                    save_dir = temp_dir,
                                    prefix = "test",
                                    part_size = 100,
                                    cli_log = dummy_cli_log),
                  error = function(e) conditionMessage(e)
  )
  # rename happened → filename not listed as missing
  expect_false(grepl("\\bfilename\\b", err))
  # missing poly_tail_length (and others) IS listed
  expect_true(grepl("poly_tail", err))
})


test_that("preprocess_inputs qc-filters and errors on empty output (lines 485-659 partial)", {
  # All reads fail filter_dorado_summary (alignment_direction = "*").
  # After filtering: 0 reads → 0 parts → pod5 loop is skipped →
  # polya_signal_files is empty → error at line 644.
  temp_dir <- file.path(tempdir(), paste0("test_qc_empty_", Sys.getpid()))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  all_invalid <- create_dummy_dorado_summary(n = 20, all_invalid = TRUE)

  expect_error(
    preprocess_inputs(dorado_summary = all_invalid,
                      pod5_dir = tempdir(),
                      num_cores = 1,
                      qc = TRUE,
                      save_dir = temp_dir,
                      prefix = "test",
                      part_size = 100,
                      cli_log = dummy_cli_log),
    "No poly\\(A\\) signals were successfully extracted")
})


test_that("preprocess_inputs with qc=FALSE splits summary and errors at empty pod5 output", {
  # qc=FALSE bypasses filter_dorado_summary (lines 494-519 skipped).
  # All reads have poly_tail_start=0 so extract_tails_from_pod5 returns list()
  # → empty signals → else branch at line 632 → polya_signal_files empty → error at 644.
  temp_dir <- file.path(tempdir(), paste0("test_no_qc_", Sys.getpid()))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  summary_zero_start <- data.frame(read_id = paste0("read_", 1:5),
                                   filename = rep("file1.pod5", 5),
                                   poly_tail_length = rep(50L, 5),
                                   poly_tail_start = rep(0L, 5),
                                   poly_tail_end = rep(600L, 5),
                                   alignment_genome = rep("tx_1", 5),
                                   alignment_direction = rep("+", 5),
                                   alignment_mapq = rep(60L, 5),
                                   stringsAsFactors = FALSE)

  expect_error(
    preprocess_inputs(dorado_summary = summary_zero_start,
                      pod5_dir = tempdir(),
                      num_cores = 1,
                      qc = FALSE,
                      save_dir = temp_dir,
                      prefix = "test",
                      part_size = 10,
                      cli_log = dummy_cli_log),
    "No poly\\(A\\) signals were successfully extracted"
  )
})


################################################################################
# create_outputs_dorado
################################################################################

test_that("create_outputs_dorado validates input correctly", {
  expect_error(create_outputs_dorado(nonA_temp_dir = tempdir(),
                                     polya_chunks_dir = tempdir(),
                                     num_cores = 1),
               "Dorado summary directory is missing")

  expect_error(create_outputs_dorado(dorado_summary_dir = tempdir(),
                                     polya_chunks_dir = tempdir(),
                                     num_cores = 1),
               "Non-A predictions directory is missing")

  expect_error(create_outputs_dorado(dorado_summary_dir = tempdir(),
                                     nonA_temp_dir = tempdir(),
                                     num_cores = 1),
               "PolyA chunks directory is missing")

  expect_error(create_outputs_dorado(dorado_summary_dir = tempdir(),
                                     nonA_temp_dir = tempdir(),
                                     polya_chunks_dir = tempdir()),
               "Number of cores is missing")

  expect_error(create_outputs_dorado(dorado_summary_dir = tempdir(),
                                     nonA_temp_dir = tempdir(),
                                     polya_chunks_dir = tempdir(),num_cores = -1),
               "num_cores must be a positive integer")
})


test_that("create_outputs_dorado validates directories exist", {
  expect_error(create_outputs_dorado(dorado_summary_dir = "/nonexistent/path1",
                                     nonA_temp_dir = "/nonexistent/path2",
                                     polya_chunks_dir = "/nonexistent/path3",
                                     num_cores = 1),
               "All directory arguments must be valid existing paths")
})


test_that("create_outputs_dorado checks for files in directories", {
  temp_summary_dir <- file.path(tempdir(), paste0("empty_summary_", Sys.getpid()))
  temp_nonA_dir <- file.path(tempdir(), paste0("empty_nonA_", Sys.getpid()))
  temp_chunks_dir <- file.path(tempdir(), paste0("empty_chunks_", Sys.getpid()))
  dir.create(temp_summary_dir, showWarnings = FALSE)
  dir.create(temp_nonA_dir, showWarnings = FALSE)
  dir.create(temp_chunks_dir, showWarnings = FALSE)
  on.exit({
    unlink(temp_summary_dir, recursive = TRUE)
    unlink(temp_nonA_dir, recursive = TRUE)
    unlink(temp_chunks_dir, recursive = TRUE)
  }, add = TRUE)

  expect_error(
    create_outputs_dorado(dorado_summary_dir = temp_summary_dir,
                          nonA_temp_dir = temp_nonA_dir,
                          polya_chunks_dir = temp_chunks_dir,
                          num_cores = 1,
                          original_summary = create_dummy_dorado_summary(n = 10)),
    "No summary files found")
})


test_that("create_outputs_dorado errors when nonA_temp_dir has no RDS files (lines 1331-1333)", {
  base_dir <- file.path(tempdir(), paste0("test_no_nonA_", Sys.getpid()))
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  fix <- create_outputs_dorado_fixtures(base_dir)
  # Remove prediction files to trigger the no-predictions error
  file.remove(list.files(fix$nonA_dir, full.names = TRUE))

  expect_error(
    create_outputs_dorado(dorado_summary_dir = fix$summary_dir,
                          nonA_temp_dir = fix$nonA_dir,
                          polya_chunks_dir = fix$chunks_dir,
                          num_cores = 1,
                          original_summary = fix$original_summary),
    "No prediction files found"
  )
})


test_that("create_outputs_dorado errors when polya_chunks_dir has no RDS files (lines 1328-1330)", {
  base_dir <- file.path(tempdir(), paste0("test_no_chunks_", Sys.getpid()))
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  fix <- create_outputs_dorado_fixtures(base_dir)
  file.remove(list.files(fix$chunks_dir, full.names = TRUE))

  expect_error(
    create_outputs_dorado(dorado_summary_dir = fix$summary_dir,
                          nonA_temp_dir = fix$nonA_dir,
                          polya_chunks_dir = fix$chunks_dir,
                          num_cores = 1,
                          original_summary = fix$original_summary),
    "No RDS files found in polya_chunks_dir"
  )
})


test_that("create_outputs_dorado errors on invalid original_summary type (line 1400-1401)", {
  base_dir <- file.path(tempdir(), paste0("test_bad_orig_", Sys.getpid()))
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  fix <- create_outputs_dorado_fixtures(base_dir)

  expect_error(
    create_outputs_dorado(dorado_summary_dir = fix$summary_dir,
                          nonA_temp_dir = fix$nonA_dir,
                          polya_chunks_dir = fix$chunks_dir,
                          num_cores = 1,
                          original_summary = 12345),
    "original_summary must be a file path or data frame"
  )
})


test_that("create_outputs_dorado returns correct output structure with qc=TRUE", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("vroom")

  base_dir <- file.path(tempdir(), paste0("test_outputs_qc_", Sys.getpid()))
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  fix <- create_outputs_dorado_fixtures(base_dir)

  result <- suppressMessages(suppressWarnings(
    create_outputs_dorado(dorado_summary_dir = fix$summary_dir,
                          nonA_temp_dir = fix$nonA_dir,
                          polya_chunks_dir = fix$chunks_dir,
                          num_cores = 1,
                          qc = TRUE,
                          original_summary = fix$original_summary)))

  # Output must be a named list with both tables
  expect_type(result, "list")
  expect_true("read_classes" %in% names(result))
  expect_true("nonadenosine_residues" %in% names(result))
  expect_true(is.data.frame(result[["read_classes"]]))
  expect_true(is.data.frame(result[["nonadenosine_residues"]]))

  # read_classes column structure
  rc_cols <- c("readname", "contig", "polya_length", "qc_tag", "class", "comments")
  expect_true(all(rc_cols %in% colnames(result[["read_classes"]])))

  # All 6 reads from original_summary must appear in read_classes
  expect_equal(nrow(result[["read_classes"]]), 6L)

  # readA-00001 must be decorated
  expect_equal(
    result[["read_classes"]]$class[result[["read_classes"]]$readname == "readA-00001"],
    "decorated"
  )

  # Comment code breakdown
  comments <- result[["read_classes"]]$comments
  expect_true("YAY" %in% comments)
  expect_true("MPU" %in% comments)
  expect_true("UNM" %in% comments)
  expect_true("IRL" %in% comments)
  expect_true("BAC" %in% comments)
  expect_true("MAU" %in% comments)

  # nonadenosine_residues column structure
  nr_cols <- c("readname", "contig", "prediction", "est_nonA_pos", "polya_length", "qc_tag")
  expect_true(all(nr_cols %in% colnames(result[["nonadenosine_residues"]])))

  # Only read_00001 has a non-A residue (C, pos 50)
  expect_equal(nrow(result[["nonadenosine_residues"]]), 1L)
  expect_equal(result[["nonadenosine_residues"]]$prediction, "C")
  expect_equal(result[["nonadenosine_residues"]]$est_nonA_pos, 50L)
})


test_that("create_outputs_dorado returns correct output structure with qc=FALSE", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("vroom")

  base_dir <- file.path(tempdir(), paste0("test_outputs_noqc_", Sys.getpid()))
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  fix <- create_outputs_dorado_fixtures(base_dir)

  result <- suppressMessages(suppressWarnings(
    create_outputs_dorado(dorado_summary_dir = fix$summary_dir,
                          nonA_temp_dir = fix$nonA_dir,
                          polya_chunks_dir = fix$chunks_dir,
                          num_cores = 1,
                          qc = FALSE,
                          original_summary = fix$original_summary)))

  expect_type(result, "list")
  expect_true(is.data.frame(result[["read_classes"]]))
  expect_true(is.data.frame(result[["nonadenosine_residues"]]))
  expect_equal(nrow(result[["read_classes"]]), 6L)
})


test_that("create_outputs_dorado accepts file path for original_summary", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("vroom")

  base_dir <- file.path(tempdir(), paste0("test_orig_path_", Sys.getpid()))
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  fix <- create_outputs_dorado_fixtures(base_dir)

  # Write original_summary to a TSV file and pass its path
  orig_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(orig_path), add = TRUE)
  vroom::vroom_write(fix$original_summary, orig_path, delim = "\t")

  result <- suppressMessages(suppressWarnings(
    create_outputs_dorado(dorado_summary_dir = fix$summary_dir,
                          nonA_temp_dir = fix$nonA_dir,
                          polya_chunks_dir = fix$chunks_dir,num_cores = 1,
                          qc = TRUE,
                          original_summary = orig_path)))

  expect_type(result, "list")
  expect_equal(nrow(result[["read_classes"]]), 6L)
})


################################################################################
# Integration: signal processing consistency
################################################################################

test_that("signal processing pipeline produces consistent results", {
  set.seed(123)
  signal <- c(rnorm(150, mean = 680, sd = 10),
              rnorm(30,  mean = 750, sd = 8),
              rnorm(150, mean = 680, sd = 10),
              rnorm(20,  mean = 580, sd = 8),
              rnorm(50,  mean = 680, sd = 10))

  signal <- as.integer(round(signal))

  pseudomoves  <- filter_signal_by_threshold(signal)
  pseudomoves2 <- filter_signal_by_threshold(signal)

  # Result must be deterministic
  expect_identical(pseudomoves, pseudomoves2)
  # Signal with embedded modification peaks should produce some non-zero pseudomoves
  expect_true(sum(pseudomoves != 0) > 0)
})

