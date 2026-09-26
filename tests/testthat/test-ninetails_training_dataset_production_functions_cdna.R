################################################################################
# Testing training dataset production functions (Dorado cDNA)
################################################################################

# Helpers
################################################################################

#' Create a synthetic signal list mimicking the output of
#' extract_tail_signals_trainingset_cdna()$polya_signals: a named list of
#' integer vectors (winsorized, downsampled tail signals).
#' Reads with add_run = TRUE carry a block of deviating values so that
#' filter_signal_by_threshold_trainingset() reports a pseudomove run.
#' @keywords internal
make_cdna_signal_list <- function(n_reads = 3,
                                  signal_length = 300,
                                  add_run = TRUE,
                                  run_start = 150,
                                  run_length = 12,
                                  run_shift = 200L) {
  set.seed(123)
  signal_list <- lapply(seq_len(n_reads), function(i) {
    signal <- as.integer(round(rnorm(signal_length, mean = 680, sd = 5)))
    if (add_run) {
      idx <- run_start:(run_start + run_length - 1)
      signal[idx] <- signal[idx] + run_shift
    }
    return(signal)
  })
  names(signal_list) <- paste0("test-read-", sprintf("%03d", seq_len(n_reads)))
  return(signal_list)
}


#' Create a synthetic tail_feature_list mimicking the structure produced by
#' create_tail_feature_list_trainingset_cdna().
#'
#' Structure: list[[1]][[readname]][[1..4]]
#'   [[1]] pod5_filename (NA)
#'   [[2]] tail_signal (numeric vector)
#'   [[3]] tail_moves (NA)
#'   [[4]] tail_pseudomoves (integer vector)
#' @keywords internal
make_cdna_feature_list <- function(readname = "test-read-001",
                                   signal_length = 300,
                                   peak_runs = 1,
                                   valley_runs = 1,
                                   run_length = 6) {
  signal <- as.integer(rnorm(signal_length, mean = 680, sd = 15))
  pseudomoves <- integer(signal_length)
  pos <- 50
  for (i in seq_len(peak_runs)) {
    pseudomoves[pos:(pos + run_length - 1)] <- 1L
    pos <- pos + 40
  }
  for (i in seq_len(valley_runs)) {
    pseudomoves[pos:(pos + run_length - 1)] <- -1L
    pos <- pos + 40
  }
  feature_list <- list(list(
    list(pod5_filename = NA_character_,
         tail_signal = signal,
         tail_moves = NA_real_,
         tail_pseudomoves = pseudomoves)
  ))
  names(feature_list[[1]]) <- readname
  return(list(tail_feature_list = feature_list[[1]],
              discarded_readnames = character(0)))
}


################################################################################
# extract_tail_signals_trainingset_cdna - guard tests only (BAM/POD5 I/O not tested)
################################################################################

test_that("extract_tail_signals_trainingset_cdna errors on missing dorado_summary", {
  expect_error(
    extract_tail_signals_trainingset_cdna(bam_file = "x.bam",
                                          pod5_dir = tempdir()),
    "Dorado summary is missing"
  )
})

test_that("extract_tail_signals_trainingset_cdna errors on missing bam_file", {
  expect_error(
    extract_tail_signals_trainingset_cdna(dorado_summary = data.frame(read_id = "r1"),
                                          pod5_dir = tempdir()),
    "BAM file is missing"
  )
})

test_that("extract_tail_signals_trainingset_cdna errors on missing pod5_dir", {
  expect_error(
    extract_tail_signals_trainingset_cdna(dorado_summary = data.frame(read_id = "r1"),
                                          bam_file = "x.bam"),
    "Directory with POD5 files is missing"
  )
})

test_that("extract_tail_signals_trainingset_cdna errors on non-numeric num_cores", {
  expect_error(
    extract_tail_signals_trainingset_cdna(dorado_summary = data.frame(read_id = "r1"),
                                          bam_file = "x.bam",
                                          pod5_dir = tempdir(),
                                          num_cores = "two"),
    "Declared core number must be numeric"
  )
})

test_that("extract_tail_signals_trainingset_cdna errors on non-existent BAM file", {
  expect_error(
    extract_tail_signals_trainingset_cdna(dorado_summary = data.frame(read_id = "r1"),
                                          bam_file = file.path(tempdir(), "does_not_exist.bam"),
                                          pod5_dir = tempdir(),
                                          num_cores = 1)
  )
})

test_that("extract_tail_signals_trainingset_cdna errors on missing summary columns", {
  bam_stub <- file.path(tempdir(), "stub_trainingset.bam")
  writeLines("stub", bam_stub)
  on.exit(unlink(bam_stub))
  expect_error(
    extract_tail_signals_trainingset_cdna(dorado_summary = data.frame(read_id = "r1"),
                                          bam_file = bam_stub,
                                          pod5_dir = tempdir(),
                                          num_cores = 1),
    "Required columns missing"
  )
})

test_that("extract_tail_signals_trainingset_cdna errors on empty summary data frame", {
  bam_stub <- file.path(tempdir(), "stub_trainingset.bam")
  writeLines("stub", bam_stub)
  on.exit(unlink(bam_stub))
  expect_error(
    extract_tail_signals_trainingset_cdna(dorado_summary = data.frame(),
                                          bam_file = bam_stub,
                                          pod5_dir = tempdir(),
                                          num_cores = 1),
    "Empty data frame"
  )
})


################################################################################
# create_tail_feature_list_trainingset_cdna
################################################################################

test_that("create_tail_feature_list_trainingset_cdna errors when signal_list is missing", {
  expect_error(
    create_tail_feature_list_trainingset_cdna(num_cores = 1, nucleotide = "C"),
    "Signal list is missing"
  )
})

test_that("create_tail_feature_list_trainingset_cdna errors when num_cores is missing", {
  expect_error(
    create_tail_feature_list_trainingset_cdna(signal_list = make_cdna_signal_list(),
                                              nucleotide = "C"),
    "Number of declared cores is missing"
  )
})

test_that("create_tail_feature_list_trainingset_cdna errors when nucleotide is missing", {
  expect_error(
    create_tail_feature_list_trainingset_cdna(signal_list = make_cdna_signal_list(),
                                              num_cores = 1),
    "Nucleotide is missing"
  )
})

test_that("create_tail_feature_list_trainingset_cdna errors on wrong nucleotide", {
  expect_error(
    create_tail_feature_list_trainingset_cdna(signal_list = make_cdna_signal_list(),
                                              num_cores = 1,
                                              nucleotide = "T"),
    "Wrong nucleotide selected"
  )
})

test_that("create_tail_feature_list_trainingset_cdna errors on unnamed signal_list", {
  expect_error(
    create_tail_feature_list_trainingset_cdna(signal_list = unname(make_cdna_signal_list()),
                                              num_cores = 1,
                                              nucleotide = "C"),
    "must be named by read IDs"
  )
})

test_that("create_tail_feature_list_trainingset_cdna errors on non-numeric signals", {
  expect_error(
    create_tail_feature_list_trainingset_cdna(signal_list = list(r1 = "a", r2 = "b"),
                                              num_cores = 1,
                                              nucleotide = "C"),
    "must be numeric vectors"
  )
})

test_that("create_tail_feature_list_trainingset_cdna keeps four-slot layout and retains runs for non-A", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")
  signal_list <- make_cdna_signal_list(n_reads = 3, add_run = TRUE)
  res <- suppressMessages(
    create_tail_feature_list_trainingset_cdna(signal_list = signal_list,
                                              num_cores = 1,
                                              nucleotide = "C")
  )
  expect_type(res, "list")
  expect_true(all(c("tail_feature_list", "discarded_readnames") %in% names(res)))
  expect_true(length(res$tail_feature_list) > 0)
  read <- res$tail_feature_list[[1]]
  expect_equal(names(read),
               c("pod5_filename", "tail_signal", "tail_moves", "tail_pseudomoves"))
  expect_equal(length(read[[2]]), length(read[[4]]))
  expect_true(all(read[[4]] %in% c(-1, 0, 1)))
  expect_true(any(with(rle(read[[4]]), lengths[values != 0] >= 4)))
})

test_that("create_tail_feature_list_trainingset_cdna discards run-carrying reads for A", {
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")
  signal_list <- make_cdna_signal_list(n_reads = 3, add_run = TRUE)
  res <- suppressMessages(
    create_tail_feature_list_trainingset_cdna(signal_list = signal_list,
                                              num_cores = 1,
                                              nucleotide = "A")
  )
  expect_equal(length(res$tail_feature_list) + length(res$discarded_readnames),
               length(signal_list))
  # every retained read has no qualifying run
  for (read in res$tail_feature_list) {
    expect_false(any(with(rle(read[[4]]), lengths[values != 0] >= 4)))
  }
})


################################################################################
# count_pseudomove_runs_trainingset_cdna
################################################################################

test_that("count_pseudomove_runs_trainingset_cdna errors when tail_feature_list is missing", {
  expect_error(
    count_pseudomove_runs_trainingset_cdna(),
    "List of tail features is missing"
  )
})

test_that("count_pseudomove_runs_trainingset_cdna errors on non-list input", {
  expect_error(
    count_pseudomove_runs_trainingset_cdna(tail_feature_list = "not a list"),
    "not a list"
  )
})

test_that("count_pseudomove_runs_trainingset_cdna errors on invalid min_run_length", {
  expect_error(
    count_pseudomove_runs_trainingset_cdna(tail_feature_list = make_cdna_feature_list(),
                                           min_run_length = 0),
    "positive number"
  )
})

test_that("count_pseudomove_runs_trainingset_cdna counts peaks and valleys per read", {
  tfl <- make_cdna_feature_list(peak_runs = 2, valley_runs = 1, run_length = 6)
  res <- count_pseudomove_runs_trainingset_cdna(tfl)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("readname", "signal_length", "peak_runs", "valley_runs") %in% colnames(res)))
  expect_equal(nrow(res), 1)
  expect_equal(res$peak_runs, 2)
  expect_equal(res$valley_runs, 1)
  expect_equal(res$signal_length, 300)
})

test_that("count_pseudomove_runs_trainingset_cdna ignores runs shorter than min_run_length", {
  tfl <- make_cdna_feature_list(peak_runs = 1, valley_runs = 1, run_length = 3)
  res <- count_pseudomove_runs_trainingset_cdna(tfl, min_run_length = 4)
  expect_equal(res$peak_runs, 0)
  expect_equal(res$valley_runs, 0)
})

test_that("count_pseudomove_runs_trainingset_cdna returns empty table for empty feature list", {
  res <- count_pseudomove_runs_trainingset_cdna(list(tail_feature_list = list(),
                                                     discarded_readnames = character(0)))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
})


################################################################################
# plot_tail_features_trainingset_cdna
################################################################################

test_that("plot_tail_features_trainingset_cdna errors when readname is missing", {
  expect_error(
    plot_tail_features_trainingset_cdna(tail_feature_list = make_cdna_feature_list()),
    "Readname is missing"
  )
})

test_that("plot_tail_features_trainingset_cdna errors when tail_feature_list is missing", {
  expect_error(
    plot_tail_features_trainingset_cdna(readname = "test-read-001"),
    "List of tail features is missing"
  )
})

test_that("plot_tail_features_trainingset_cdna errors on unknown readname", {
  expect_error(
    plot_tail_features_trainingset_cdna(readname = "no-such-read",
                                        tail_feature_list = make_cdna_feature_list()),
    "not present in the tail_feature_list"
  )
})

test_that("plot_tail_features_trainingset_cdna returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  p <- plot_tail_features_trainingset_cdna(readname = "test-read-001",
                                           tail_feature_list = make_cdna_feature_list())
  expect_s3_class(p, "gg")
})

test_that("plot_tail_features_trainingset_cdna works for reads without runs", {
  skip_if_not_installed("ggplot2")
  tfl <- make_cdna_feature_list(peak_runs = 0, valley_runs = 0)
  p <- plot_tail_features_trainingset_cdna(readname = "test-read-001",
                                           tail_feature_list = tfl)
  expect_s3_class(p, "gg")
})


################################################################################
# prepare_trainingset_cdna - guard tests only
################################################################################

test_that("prepare_trainingset_cdna errors when nucleotide is missing", {
  expect_error(
    prepare_trainingset_cdna(tail_type = "polyA",
                             dorado_summary = data.frame(read_id = "r1"),
                             bam_file = "x.bam",
                             pod5_dir = tempdir()),
    "Nucleotide is missing"
  )
})

test_that("prepare_trainingset_cdna errors when tail_type is missing", {
  expect_error(
    prepare_trainingset_cdna(nucleotide = "C",
                             dorado_summary = data.frame(read_id = "r1"),
                             bam_file = "x.bam",
                             pod5_dir = tempdir()),
    "Tail type is missing"
  )
})

test_that("prepare_trainingset_cdna errors on wrong nucleotide", {
  expect_error(
    prepare_trainingset_cdna(nucleotide = "T",
                             tail_type = "polyA",
                             dorado_summary = data.frame(read_id = "r1"),
                             bam_file = "x.bam",
                             pod5_dir = tempdir()),
    "Wrong nucleotide selected"
  )
})

test_that("prepare_trainingset_cdna errors on wrong tail_type", {
  expect_error(
    prepare_trainingset_cdna(nucleotide = "C",
                             tail_type = "polyU",
                             dorado_summary = data.frame(read_id = "r1"),
                             bam_file = "x.bam",
                             pod5_dir = tempdir()),
    "Wrong tail type selected"
  )
})

test_that("prepare_trainingset_cdna requires polarity value for non-A residues", {
  expect_error(
    prepare_trainingset_cdna(nucleotide = "G",
                             tail_type = "polyT",
                             dorado_summary = data.frame(read_id = "r1"),
                             bam_file = "x.bam",
                             pod5_dir = tempdir()),
    "must be either 1 or -1"
  )
})
