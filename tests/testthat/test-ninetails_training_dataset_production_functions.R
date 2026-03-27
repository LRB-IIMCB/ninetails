################################################################################
# Testing training dataset production functions
################################################################################

# Helpers
################################################################################

#' Create a synthetic tail_feature_list mimicking the structure produced by
#' create_tail_feature_list_trainingset (for use with split_tail_centered_trainingset
#' and create_tail_chunk_list_trainingset).
#'
#' Structure: list[[1]][[readname]][[1..4]]
#'   [[1]] fast5_filename
#'   [[2]] tail_signal (numeric vector)
#'   [[3]] tail_moves (numeric vector)
#'   [[4]] tail_pseudomoves (integer vector)
#'
#' @keywords internal
make_training_feature_list <- function(readname = "test-read-001",
                                       signal_length = 300,
                                       add_run = TRUE,
                                       run_start = 100,
                                       run_length = 6,
                                       run_value = 1L) {
  signal <- as.integer(rnorm(signal_length, mean = 680, sd = 15))
  moves <- rep(1L, signal_length)
  pseudomoves <- integer(signal_length)
  if (add_run && run_start + run_length <= signal_length) {
    pseudomoves[run_start:(run_start + run_length - 1)] <- run_value
  }
  feature_list <- list(list(
    list("dummy.fast5", signal, moves, pseudomoves)
  ))
  names(feature_list[[1]]) <- readname
  return(feature_list)
}


#' Create a minimal tail_chunk_list for filter_nonA_chunks_trainingset tests.
#' Each chunk has chunk_sequence, chunk_start_pos, chunk_end_pos,
#' and pseudomoves with or without a qualifying run.
#' @keywords internal
make_training_chunk_list <- function(readname = "test-read-001",
                                     n_chunks = 2,
                                     pseudomove_value = 1L,
                                     run_length = 5L) {
  make_chunk <- function(pv, rl) {
    pm <- integer(100)
    pm[40:(40 + rl - 1)] <- pv
    return(list(chunk_sequence = as.integer(rnorm(100, 680, 15)),
                chunk_start_pos = 50L,
                chunk_end_pos = 149L,
                pseudomoves = pm))
  }
  chunks <- replicate(n_chunks, make_chunk(pseudomove_value, run_length),
                      simplify = FALSE)
  names(chunks) <- paste0(readname, "_", seq_len(n_chunks))
  chunk_list <- list(chunks)
  names(chunk_list) <- readname
  return(chunk_list)
}


################################################################################
# extract_tail_data_trainingset — guard tests only (fast5 I/O not tested)
################################################################################

test_that("extract_tail_data_trainingset errors on missing readname", {
  expect_error(
    extract_tail_data_trainingset(polya_summary = data.frame(readname = "r1"),
                                  workspace = tempdir(),
                                  basecall_group = "Basecall_1D_000"),
    "Readname is missing"
  )
})

test_that("extract_tail_data_trainingset errors on missing workspace", {
  expect_error(
    extract_tail_data_trainingset(readname = "read-001",
                                  polya_summary = data.frame(readname = "read-001"),
                                  basecall_group = "Basecall_1D_000"),
    "Directory with basecalled fast5s is missing"
  )
})

test_that("extract_tail_data_trainingset errors on missing basecall_group", {
  expect_error(
    extract_tail_data_trainingset(readname = "read-001",
                                  polya_summary = data.frame(readname = "read-001"),
                                  workspace = tempdir()),
    "Basecall group is missing"
  )
})

test_that("extract_tail_data_trainingset errors on missing polya_summary", {
  expect_error(
    extract_tail_data_trainingset(readname = "read-001",
                                  workspace = tempdir(),
                                  basecall_group = "Basecall_1D_000"),
    "Polya_summary is missing"
  )
})

test_that("extract_tail_data_trainingset errors on empty polya_summary", {
  expect_error(
    extract_tail_data_trainingset(readname = "read-001",
                                  polya_summary = data.frame(),
                                  workspace = tempdir(),
                                  basecall_group = "Basecall_1D_000"),
    "Empty data frame"
  )
})

test_that("extract_tail_data_trainingset errors on non-character readname", {
  expect_error(
    extract_tail_data_trainingset(readname = 12345,
                                  polya_summary = data.frame(readname = "r1"),
                                  workspace = tempdir(),
                                  basecall_group = "Basecall_1D_000"),
    "not a character string"
  )
})


################################################################################
# filter_signal_by_threshold_trainingset
#
# BUG NOTE: lines 515-519 reference `baseline` and `std_cutoff` without prior
# declaration inside the function body. In a clean R session both objects are
# absent, causing "object 'baseline' not found". The function appears to rely
# on globally pre-allocated vectors from interactive use. The happy-path tests
# below will pass only after the bug is fixed by adding, before the loop at
# line 523:
#   baseline <- numeric(length(adjusted_signal))
#   std_cutoff <- numeric(length(adjusted_signal))
################################################################################

test_that("filter_signal_by_threshold_trainingset errors when signal is missing", {
  expect_error(
    filter_signal_by_threshold_trainingset(),
    "Signal is missing"
  )
})

test_that("filter_signal_by_threshold_trainingset returns numeric vector of same length as input", {
  set.seed(42)
  signal <- as.integer(rnorm(200, mean = 680, sd = 15))
  result <- filter_signal_by_threshold_trainingset(signal)
  expect_type(result, "double")
  expect_equal(length(result), length(signal))
})

test_that("filter_signal_by_threshold_trainingset output values are in {-1, 0, 1}", {
  set.seed(42)
  signal <- as.integer(rnorm(200, mean = 680, sd = 15))
  result <- filter_signal_by_threshold_trainingset(signal)
  expect_true(all(result %in% c(-1, 0, 1)))
})

test_that("filter_signal_by_threshold_trainingset zeroes out terminal positions", {
  set.seed(42)
  signal <- as.integer(rnorm(200, mean = 680, sd = 15))
  result <- filter_signal_by_threshold_trainingset(signal)
  # First 5 and last 6 positions hardcoded to 0
  expect_true(all(result[1:5] == 0))
  expect_true(all(result[(length(result) - 5):length(result)] == 0))
})


################################################################################
# substitute_gaps
################################################################################

test_that("substitute_gaps fills isolated single-zero gap between same-sign runs", {
  pm <- c(1L, 1L, 1L, 0L, 1L, 1L, 1L)
  result <- substitute_gaps(pm)
  expect_true(all(result == 1L))
})

test_that("substitute_gaps does not fill multi-zero gaps", {
  pm <- c(1L, 1L, 0L, 0L, 1L, 1L)
  result <- substitute_gaps(pm)
  expect_equal(result[3:4], c(0L, 0L))
})

test_that("substitute_gaps returns vector of same length as input", {
  pm <- c(1L, 0L, 1L, -1L, 0L, -1L, 0L, 0L, 1L)
  expect_equal(length(substitute_gaps(pm)), length(pm))
})

test_that("substitute_gaps does not modify all-zero input", {
  pm <- integer(10)
  expect_equal(substitute_gaps(pm), pm)
})


################################################################################
# split_tail_centered_trainingset
################################################################################

test_that("split_tail_centered_trainingset errors when readname is missing", {
  tfl <- make_training_feature_list()
  expect_error(
    split_tail_centered_trainingset(tail_feature_list = tfl),
    "Readname is missing"
  )
})

test_that("split_tail_centered_trainingset errors when tail_feature_list is missing", {
  expect_error(
    split_tail_centered_trainingset(readname = "test-read-001"),
    "List of tail features is missing"
  )
})

test_that("split_tail_centered_trainingset errors on non-character readname", {
  tfl <- make_training_feature_list()
  expect_error(
    split_tail_centered_trainingset(readname = 12345,
                                    tail_feature_list = tfl),
    "not a character string"
  )
})

test_that("split_tail_centered_trainingset errors on non-list tail_feature_list", {
  expect_error(
    split_tail_centered_trainingset(readname = "test-read-001",
                                    tail_feature_list = "not_a_list"),
    "not a list"
  )
})

test_that("split_tail_centered_trainingset returns list with correct chunk structure", {
  tfl <- make_training_feature_list(add_run = TRUE, run_length = 6)

  result <- split_tail_centered_trainingset(readname = "test-read-001",
                                            tail_feature_list = tfl)

  expect_type(result, "list")
  if (length(result) > 0) {
    first_chunk <- result[[1]]
    expect_true("chunk_sequence" %in% names(first_chunk))
    expect_true("chunk_start_pos" %in% names(first_chunk))
    expect_true("chunk_end_pos" %in% names(first_chunk))
    expect_true("pseudomoves" %in% names(first_chunk))
    expect_equal(length(first_chunk$chunk_sequence), 100)
    expect_equal(length(first_chunk$pseudomoves), 100)
  }
})

test_that("split_tail_centered_trainingset returns empty list when no pseudomove run >= 4", {
  tfl <- make_training_feature_list(add_run = FALSE)

  result <- split_tail_centered_trainingset(readname = "test-read-001",
                                            tail_feature_list = tfl)

  expect_equal(length(result), 0L)
})

test_that("split_tail_centered_trainingset chunk names follow <readname>_<index> convention", {
  tfl <- make_training_feature_list(add_run = TRUE, run_length = 6)

  result <- split_tail_centered_trainingset(readname = "test-read-001",
                                            tail_feature_list = tfl)

  if (length(result) > 0) {
    expect_true(all(grepl("^test-read-001_", names(result))))
  }
})


################################################################################
# split_with_overlaps
################################################################################

test_that("split_with_overlaps returns list of numeric vectors", {
  tfl <- make_training_feature_list(signal_length = 400)

  result <- split_with_overlaps(readname = "test-read-001",
                                tail_feature_list = tfl,
                                segment = 100,
                                overlap = 50)

  expect_type(result, "list")
  expect_true(length(result) > 0)
  expect_true(all(sapply(result, is.numeric)))
})

test_that("split_with_overlaps all output vectors have length equal to segment", {
  tfl <- make_training_feature_list(signal_length = 400)

  result <- split_with_overlaps(readname = "test-read-001",
                                tail_feature_list = tfl,
                                segment = 100,
                                overlap = 50)

  expect_true(all(sapply(result, length) == 100))
})

test_that("split_with_overlaps output contains no NAs", {
  # Signal length 320 is not evenly divisible by step=50, so the last chunk
  # will have trailing NAs that must be imputed
  tfl <- make_training_feature_list(signal_length = 320)

  result <- split_with_overlaps(readname = "test-read-001",
                                tail_feature_list = tfl,
                                segment = 100,
                                overlap = 50)

  expect_false(any(sapply(result, function(x) any(is.na(x)))))
})

test_that("split_with_overlaps produces expected number of chunks", {
  # signal 400, step = 100-50 = 50 → starts: 1,51,101,...,351 → 8 chunks
  tfl <- make_training_feature_list(signal_length = 400)

  result <- split_with_overlaps(readname = "test-read-001",
                                tail_feature_list = tfl,
                                segment = 100,
                                overlap = 50)

  expected_n <- length(seq(1, 400, by = 100 - 50))
  expect_equal(length(result), expected_n)
})


################################################################################
# create_tail_chunk_list_trainingset — guard tests
################################################################################

test_that("create_tail_chunk_list_trainingset errors when num_cores is missing", {
  expect_error(
    create_tail_chunk_list_trainingset(tail_feature_list = list()),
    "Number of declared cores is missing"
  )
})

test_that("create_tail_chunk_list_trainingset errors when tail_feature_list is missing", {
  expect_error(
    create_tail_chunk_list_trainingset(num_cores = 1),
    "List of features is missing"
  )
})

test_that("create_tail_chunk_list_trainingset errors on non-numeric num_cores", {
  expect_error(
    create_tail_chunk_list_trainingset(tail_feature_list = list(),
                                       num_cores = "one"),
    "Declared core number must be numeric"
  )
})

test_that("create_tail_chunk_list_trainingset errors on non-list tail_feature_list", {
  expect_error(
    create_tail_chunk_list_trainingset(tail_feature_list = "not_a_list",
                                       num_cores = 1),
    "not a list"
  )
})


################################################################################
# create_tail_chunk_list_A — guard tests
################################################################################

test_that("create_tail_chunk_list_A errors when num_cores is missing", {
  expect_error(
    create_tail_chunk_list_A(tail_feature_list = list()),
    "Number of declared cores is missing"
  )
})

test_that("create_tail_chunk_list_A errors when tail_feature_list is missing", {
  expect_error(
    create_tail_chunk_list_A(num_cores = 1),
    "List of features is missing"
  )
})

test_that("create_tail_chunk_list_A errors on non-numeric num_cores", {
  expect_error(
    create_tail_chunk_list_A(tail_feature_list = list(), num_cores = "one"),
    "Declared core number must be numeric"
  )
})

test_that("create_tail_chunk_list_A errors on non-list tail_feature_list", {
  expect_error(
    create_tail_chunk_list_A(tail_feature_list = "not_a_list", num_cores = 1),
    "not a list"
  )
})


################################################################################
# create_gaf_list_A — guard tests
################################################################################

test_that("create_gaf_list_A errors when num_cores is missing", {
  expect_error(
    create_gaf_list_A(tail_chunk_list = list()),
    "Number of declared cores is missing"
  )
})

test_that("create_gaf_list_A errors when tail_chunk_list is missing", {
  expect_error(
    create_gaf_list_A(num_cores = 1),
    "List of tail chunks is missing"
  )
})

test_that("create_gaf_list_A errors on non-numeric num_cores", {
  expect_error(
    create_gaf_list_A(tail_chunk_list = list(), num_cores = "one"),
    "Declared core number must be numeric"
  )
})


################################################################################
# filter_nonA_chunks_trainingset — guard tests + computation
################################################################################

test_that("filter_nonA_chunks_trainingset errors when num_cores is missing", {
  expect_error(
    filter_nonA_chunks_trainingset(tail_chunk_list = list(), value = 1),
    "Number of declared cores is missing"
  )
})

test_that("filter_nonA_chunks_trainingset errors when tail_chunk_list is missing", {
  expect_error(
    filter_nonA_chunks_trainingset(value = 1, num_cores = 1),
    "List of tail chunks is missing"
  )
})

test_that("filter_nonA_chunks_trainingset errors on non-list tail_chunk_list", {
  expect_error(
    filter_nonA_chunks_trainingset(tail_chunk_list = "not_a_list",
                                   value = 1,
                                   num_cores = 1),
    "not a list"
  )
})

test_that("filter_nonA_chunks_trainingset errors on non-numeric num_cores", {
  expect_error(
    filter_nonA_chunks_trainingset(tail_chunk_list = list(),
                                   value = 1,
                                   num_cores = "one"),
    "Declared core number must be numeric"
  )
})

test_that("filter_nonA_chunks_trainingset errors on non-numeric value", {
  expect_error(
    filter_nonA_chunks_trainingset(tail_chunk_list = list(),
                                   value = "G",
                                   num_cores = 1),
    "value must be numeric"
  )
})

test_that("filter_nonA_chunks_trainingset errors when value is not 1 or -1", {
  expect_error(
    filter_nonA_chunks_trainingset(tail_chunk_list = list(),
                                   value = 0,
                                   num_cores = 1),
    "value must be either 1 or -1"
  )
  expect_error(
    filter_nonA_chunks_trainingset(tail_chunk_list = list(),
                                   value = 2,
                                   num_cores = 1),
    "value must be either 1 or -1"
  )
})

test_that("filter_nonA_chunks_trainingset keeps chunks with qualifying peak run (value=1)", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  tcl <- make_training_chunk_list(pseudomove_value = 1L, run_length = 5L)

  result <- suppressMessages(
    filter_nonA_chunks_trainingset(tail_chunk_list = tcl,
                                   value = 1,
                                   num_cores = 1)
  )

  expect_type(result, "list")
  expect_true(length(result) > 0)
})

test_that("filter_nonA_chunks_trainingset keeps chunks with qualifying valley run (value=-1)", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  tcl <- make_training_chunk_list(pseudomove_value = -1L, run_length = 5L)

  result <- suppressMessages(
    filter_nonA_chunks_trainingset(tail_chunk_list = tcl,
                                   value = -1,
                                   num_cores = 1)
  )

  expect_type(result, "list")
  expect_true(length(result) > 0)
})

test_that("filter_nonA_chunks_trainingset removes chunks when pseudomove sign does not match value", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  # Chunks have valley (-1); filtering for peaks (value=1) should yield nothing
  tcl <- make_training_chunk_list(pseudomove_value = -1L, run_length = 5L)

  result <- suppressMessages(
    filter_nonA_chunks_trainingset(tail_chunk_list = tcl,
                                   value = 1,
                                   num_cores = 1)
  )

  expect_equal(length(result), 0L)
})

test_that("filter_nonA_chunks_trainingset removes chunks with qualifying run shorter than 4", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")

  tcl <- make_training_chunk_list(pseudomove_value = 1L, run_length = 3L)

  result <- suppressMessages(
    filter_nonA_chunks_trainingset(tail_chunk_list = tcl,
                                   value = 1,
                                   num_cores = 1)
  )

  expect_equal(length(result), 0L)
})


################################################################################
# create_tail_feature_list_trainingset — guard tests only
################################################################################

test_that("create_tail_feature_list_trainingset errors when num_cores is missing", {
  expect_error(
    create_tail_feature_list_trainingset(nanopolish = "f.tsv",
                                         sequencing_summary = "s.txt",
                                         workspace = tempdir(),
                                         basecall_group = "Basecall_1D_000"),
    "Number of declared cores is missing"
  )
})

test_that("create_tail_feature_list_trainingset errors when basecall_group is missing", {
  expect_error(
    create_tail_feature_list_trainingset(nanopolish = "f.tsv",
                                         sequencing_summary = "s.txt",
                                         workspace = tempdir(),
                                         num_cores = 1),
    "Basecall group is missing"
  )
})

test_that("create_tail_feature_list_trainingset errors when workspace is missing", {
  expect_error(
    create_tail_feature_list_trainingset(nanopolish = "f.tsv",
                                         sequencing_summary = "s.txt",
                                         num_cores = 1,
                                         basecall_group = "Basecall_1D_000"),
    "Directory with basecalled fast5s.*is missing"
  )
})

test_that("create_tail_feature_list_trainingset errors on non-numeric num_cores", {
  expect_error(
    create_tail_feature_list_trainingset(nanopolish = "f.tsv",
                                         sequencing_summary = "s.txt",
                                         workspace = tempdir(),
                                         num_cores = "one",
                                         basecall_group = "Basecall_1D_000"),
    "Declared core number must be numeric"
  )
})


################################################################################
# create_tail_feature_list_A — guard tests only
################################################################################

test_that("create_tail_feature_list_A errors when num_cores is missing", {
  expect_error(
    create_tail_feature_list_A(nanopolish = "f.tsv",
                               sequencing_summary = "s.txt",
                               workspace = tempdir(),
                               basecall_group = "Basecall_1D_000"),
    "Number of declared cores is missing"
  )
})

test_that("create_tail_feature_list_A errors when basecall_group is missing", {
  expect_error(
    create_tail_feature_list_A(nanopolish = "f.tsv",
                               sequencing_summary = "s.txt",
                               workspace = tempdir(),
                               num_cores = 1),
    "Basecall group is missing"
  )
})

test_that("create_tail_feature_list_A errors when workspace is missing", {
  expect_error(
    create_tail_feature_list_A(nanopolish = "f.tsv",
                               sequencing_summary = "s.txt",
                               num_cores = 1,
                               basecall_group = "Basecall_1D_000"),
    "Directory with basecalled fast5s.*is missing"
  )
})

test_that("create_tail_feature_list_A errors on non-numeric num_cores", {
  expect_error(
    create_tail_feature_list_A(nanopolish = "f.tsv",
                               sequencing_summary = "s.txt",
                               workspace = tempdir(),
                               num_cores = "one",
                               basecall_group = "Basecall_1D_000"),
    "Declared core number must be numeric"
  )
})


################################################################################
# prepare_trainingset — invalid nucleotide guard only
# All valid nucleotide branches (A, C, G, U) require fast5 files and the
# full pipeline chain; they are not tested here.
################################################################################

test_that("prepare_trainingset errors on invalid nucleotide", {
  expect_error(
    prepare_trainingset(nucleotide = "X",
                        nanopolish = "f.tsv",
                        sequencing_summary = "s.txt",
                        workspace = tempdir(),
                        num_cores = 1),
    "Wrong nucleotide selected"
  )
})

test_that("prepare_trainingset errors on lowercase nucleotide", {
  expect_error(
    prepare_trainingset(nucleotide = "a",
                        nanopolish = "f.tsv",
                        sequencing_summary = "s.txt",
                        workspace = tempdir(),
                        num_cores = 1),
    "Wrong nucleotide selected"
  )
})
