################################################################################
# Testing cDNA model functions (coordinates, signals, windows, GAFs, calling)
################################################################################

# Helpers
################################################################################

#' Path to the cDNA test data shipped with the package
#' @keywords internal
cdna_test_dir <- function() {
  return(system.file("extdata", "test_data", "pod5_cDNA", package = "ninetails"))
}

#' Create a minimal Lachesis-like table
#' @keywords internal
create_lachesis_table <- function() {
  return(data.frame(
    readname = c("r1", "r2", "r3"),
    contig = rep("beta_gal", 3),
    adapter_end = c(1000, 2000, 3000),
    polya_end = c(1700, 2600, 3900),
    polya_length = c(60, 55, 80),
    qc_tag = c("PASS", "PASS", "NOTFLAT"),
    strand = c("+", "-", "*"),
    stringsAsFactors = FALSE
  ))
}

#' Create a minimal dorado-like summary table
#' @keywords internal
create_dorado_table <- function() {
  return(data.frame(
    read_id = c("r1", "r2"),
    input_filename = c("a.pod5", "a.pod5"),
    poly_tail_start = c(1000, 2000),
    poly_tail_end = c(1700, 2600),
    poly_tail_length = c(60, 55),
    alignment_genome = c("beta_gal", "beta_gal"),
    alignment_mapping_quality = c(60, 0),
    stringsAsFactors = FALSE
  ))
}

#' Simulated raw signal with a homopolymer-like tail
#' @keywords internal
create_cdna_signal <- function(n = 1000) {
  set.seed(123) # reproducible noise
  return(round(stats::rnorm(n, mean = 600, sd = 8)))
}

#' Probability matrix (A, C, G, U) with nonA windows at given rows
#' @keywords internal
create_probabilities <- function(n, nonA_rows = integer(0), class_column = 3, p = 0.99) {
  probabilities <- matrix(0, nrow = n, ncol = 4)
  probabilities[, 1] <- 1
  probabilities[nonA_rows, 1] <- 1 - p
  probabilities[nonA_rows, class_column] <- p
  return(probabilities)
}


################################################################################
# .get_cdna_model_config
################################################################################

test_that(".get_cdna_model_config errors on wrong tail_type", {
  expect_error(ninetails:::.get_cdna_model_config("polyG"),
               "Wrong tail_type provided")
})

test_that(".get_cdna_model_config returns a complete configuration for both orientations", {
  for (tail_type in c("polyA", "polyT")) {
    config <- ninetails:::.get_cdna_model_config(tail_type)
    expect_type(config, "list")
    expect_true(all(c("model_file", "theta_mid", "theta_edge", "min_run") %in% names(config)))
    expect_true(config$theta_mid > 0 && config$theta_mid < 1)
    expect_true(config$theta_edge > 0 && config$theta_edge < 1)
  }
})

test_that(".get_cdna_model_config points to model files shipped with the package", {
  for (tail_type in c("polyA", "polyT")) {
    config <- ninetails:::.get_cdna_model_config(tail_type)
    model_path <- system.file("extdata", "cnn_model", config$model_file, package = "ninetails")
    expect_true(nchar(model_path) > 0)
  }
})


################################################################################
# read_tail_coordinates_cdna
################################################################################

test_that("read_tail_coordinates_cdna errors on missing tail_coordinates", {
  expect_error(read_tail_coordinates_cdna(),
               "Tail coordinates are missing")
})

test_that("read_tail_coordinates_cdna errors on empty data frame", {
  expect_error(read_tail_coordinates_cdna(data.frame()),
               "Empty data frame provided as an input")
})

test_that("read_tail_coordinates_cdna errors on wrong source", {
  expect_error(read_tail_coordinates_cdna(create_lachesis_table(), source = "nanopolish"),
               "Wrong source provided")
})

test_that("read_tail_coordinates_cdna errors when required columns are missing", {
  lachesis <- create_lachesis_table()
  lachesis$adapter_end <- NULL
  expect_error(read_tail_coordinates_cdna(lachesis, source = "lachesis"),
               "Required columns missing: adapter_end")
})

test_that("read_tail_coordinates_cdna converts Lachesis output", {
  result <- read_tail_coordinates_cdna(create_lachesis_table(), source = "lachesis")

  expect_true(is.data.frame(result))
  expect_true(all(c("read_id", "contig", "tail_start", "tail_end", "polya_length",
                    "qc_tag", "tail_type") %in% colnames(result)))
  expect_equal(result$read_id, c("r1", "r2", "r3"))
  expect_equal(result$tail_start, c(1000, 2000, 3000))
  expect_equal(result$tail_end, c(1700, 2600, 3900))
  # strand + is poly(A), strand - is poly(T), anything else is unidentified
  expect_equal(result$tail_type, c("polyA", "polyT", "unidentified"))
})

test_that("read_tail_coordinates_cdna converts dorado summary with new column names", {
  result <- read_tail_coordinates_cdna(create_dorado_table(), source = "dorado")

  expect_equal(result$read_id, c("r1", "r2"))
  expect_equal(result$filename, c("a.pod5", "a.pod5"))
  expect_equal(result$tail_start, c(1000, 2000))
  expect_equal(result$polya_length, c(60, 55))
  expect_equal(result$qc_tag, c("60", "0"))
  # orientation is not available in the dorado summary
  expect_true(all(is.na(result$tail_type)))
})

test_that("read_tail_coordinates_cdna keeps one row per read", {
  lachesis <- create_lachesis_table()
  result <- read_tail_coordinates_cdna(rbind(lachesis, lachesis[1, ]), source = "lachesis")
  expect_equal(nrow(result), 3L)
})

test_that("read_tail_coordinates_cdna reads the shipped Lachesis file", {
  skip_if_not_installed("vroom")
  skip_if(nchar(cdna_test_dir()) == 0, "cDNA test data not found in installed package")

  result <- read_tail_coordinates_cdna(file.path(cdna_test_dir(), "lachesis_polya.tsv"),
                                       source = "lachesis")
  expect_equal(nrow(result), 26L)
  expect_equal(as.integer(table(result$tail_type)[c("polyA", "polyT")]), c(13L, 13L))
})


################################################################################
# normalise_tail_cdna
################################################################################

test_that("normalise_tail_cdna errors on missing arguments", {
  signal <- create_cdna_signal()
  expect_error(normalise_tail_cdna(tail_start = 100, tail_end = 800),
               "Signal is missing")
  expect_error(normalise_tail_cdna(signal = signal, tail_end = 800),
               "Tail start is missing")
  expect_error(normalise_tail_cdna(signal = signal, tail_start = 100),
               "Tail end is missing")
})

test_that("normalise_tail_cdna errors on non-numeric signal", {
  expect_error(normalise_tail_cdna(signal = letters, tail_start = 1, tail_end = 10),
               "Signal must be numeric")
})

test_that("normalise_tail_cdna errors when coordinates lie outside the signal", {
  signal <- create_cdna_signal()
  expect_error(normalise_tail_cdna(signal, tail_start = 100, tail_end = 2000),
               "Tail coordinates must lie within the signal")
  expect_error(normalise_tail_cdna(signal, tail_start = 800, tail_end = 100),
               "Tail coordinates must lie within the signal")
})

test_that("normalise_tail_cdna returns values within (-1, 1) and keeps the length", {
  signal <- create_cdna_signal()
  result <- normalise_tail_cdna(signal, tail_start = 100, tail_end = 800)

  expect_equal(length(result), length(signal))
  expect_true(all(abs(result) < 1))
  # the tail median is mapped to zero
  expect_equal(stats::median(result[101:800]), 0, tolerance = 0.05)
})

test_that("normalise_tail_cdna does not depend on signal offset and scale", {
  signal <- create_cdna_signal()
  result <- normalise_tail_cdna(signal, tail_start = 100, tail_end = 800)
  rescaled <- normalise_tail_cdna(signal * 0.17 + 42, tail_start = 100, tail_end = 800)

  expect_equal(result, rescaled)
})

test_that("normalise_tail_cdna handles a flat tail", {
  result <- normalise_tail_cdna(rep(500, 300), tail_start = 50, tail_end = 250)
  expect_true(all(result == 0))
})


################################################################################
# create_tail_windows_cdna
################################################################################

test_that("create_tail_windows_cdna returns windows of 100 samples every 25 samples", {
  signal <- create_cdna_signal()
  result <- create_tail_windows_cdna(signal, tail_start = 100, tail_end = 800)

  expect_type(result, "list")
  expect_true(all(c("centres", "windows") %in% names(result)))
  expect_equal(ncol(result$windows), 100L)
  expect_equal(nrow(result$windows), length(result$centres))
  expect_true(all(diff(result$centres) == 25))
  # centres run from 50 samples before the tail to 50 samples behind it
  expect_equal(range(result$centres), c(50, 850))
})

test_that("create_tail_windows_cdna window content matches the normalised signal", {
  signal <- create_cdna_signal()
  result <- create_tail_windows_cdna(signal, tail_start = 100, tail_end = 800)
  normalised_signal <- normalise_tail_cdna(signal, tail_start = 100, tail_end = 800)

  centre <- result$centres[5]
  expect_equal(result$windows[5, ], normalised_signal[(centre - 50 + 1):(centre + 50)])
})

test_that("create_tail_windows_cdna skips windows leaving the signal", {
  signal <- create_cdna_signal(n = 500)
  # no flank on either side of the tail
  result <- create_tail_windows_cdna(signal, tail_start = 0, tail_end = 500)

  expect_true(all(result$centres - 50 >= 0))
  expect_true(all(result$centres + 50 <= 500))
})

test_that("create_tail_windows_cdna returns no windows for a signal shorter than a window", {
  result <- create_tail_windows_cdna(create_cdna_signal(n = 60), tail_start = 5, tail_end = 55)

  expect_equal(length(result$centres), 0L)
  expect_equal(dim(result$windows), c(0L, 100L))
})


################################################################################
# create_gaf_cdna
################################################################################

test_that("create_gaf_cdna errors on missing windows", {
  expect_error(create_gaf_cdna(),
               "Windows are missing")
})

test_that("create_gaf_cdna errors on non-matrix input", {
  expect_error(create_gaf_cdna(1:100),
               "Windows must be a numeric matrix")
})

test_that("create_gaf_cdna returns an array of windows x points x points x 2", {
  signal <- create_cdna_signal()
  windows <- create_tail_windows_cdna(signal, tail_start = 100, tail_end = 800)$windows
  result <- create_gaf_cdna(windows)

  expect_equal(dim(result), c(nrow(windows), 100L, 100L, 2L))
  expect_true(all(abs(result) <= 1 + 1e-8))
})

test_that("create_gaf_cdna computes summation and difference fields", {
  windows <- matrix(c(0, 0.5, -0.5, 1), nrow = 1)
  result <- create_gaf_cdna(windows)

  x <- windows[1, ]
  s <- sqrt(1 - x^2)
  expect_equal(result[1, , , 1], outer(x, x) - outer(s, s))
  expect_equal(result[1, , , 2], outer(x, s) + outer(s, x))
  # the summation field is symmetric
  expect_equal(result[1, , , 1], t(result[1, , , 1]))
})

test_that("create_gaf_cdna clips values outside the interval [-1, 1]", {
  result <- create_gaf_cdna(matrix(c(-3, 0, 3), nrow = 1))
  expect_false(any(is.nan(result)))
})


################################################################################
# call_nonA_cdna
################################################################################

test_that("call_nonA_cdna errors on missing arguments", {
  probabilities <- create_probabilities(10)
  expect_error(call_nonA_cdna(probabilities = probabilities, tail_start = 0, tail_end = 1000,
                              theta_mid = 0.9, theta_edge = 0.9),
               "Centres are missing")
  expect_error(call_nonA_cdna(centres = 1:10, tail_start = 0, tail_end = 1000,
                              theta_mid = 0.9, theta_edge = 0.9),
               "Probabilities are missing")
  expect_error(call_nonA_cdna(centres = 1:10, probabilities = probabilities,
                              tail_start = 0, tail_end = 1000),
               "Thresholds are missing")
})

test_that("call_nonA_cdna errors on malformed probabilities", {
  expect_error(call_nonA_cdna(centres = 1:10, probabilities = matrix(0, 10, 3),
                              tail_start = 0, tail_end = 1000, theta_mid = 0.9, theta_edge = 0.9),
               "Probabilities must be a matrix with 4 columns")
  expect_error(call_nonA_cdna(centres = 1:5, probabilities = create_probabilities(10),
                              tail_start = 0, tail_end = 1000, theta_mid = 0.9, theta_edge = 0.9),
               "Number of centres must match")
})

test_that("call_nonA_cdna returns an empty table when there are no windows", {
  result <- call_nonA_cdna(centres = numeric(0), probabilities = matrix(0, 0, 4),
                           tail_start = 0, tail_end = 1000, theta_mid = 0.9, theta_edge = 0.9)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0L)
  expect_true(all(c("prediction", "centr_signal_pos", "score") %in% colnames(result)))
})

test_that("call_nonA_cdna returns no calls for a pure poly(A) tail", {
  centres <- seq(50, 1050, by = 25)
  result <- call_nonA_cdna(centres, create_probabilities(length(centres)),
                           tail_start = 100, tail_end = 1000, theta_mid = 0.9, theta_edge = 0.9)
  expect_equal(nrow(result), 0L)
})

test_that("call_nonA_cdna calls a run of windows in the middle of the tail", {
  centres <- seq(50, 1050, by = 25)
  nonA_rows <- which(centres %in% c(500, 525, 550))
  probabilities <- create_probabilities(length(centres), nonA_rows, class_column = 3, p = 0.97)
  probabilities[nonA_rows[2], c(1, 3)] <- c(0.01, 0.99) # highest scoring window

  result <- call_nonA_cdna(centres, probabilities, tail_start = 100, tail_end = 1000,
                           theta_mid = 0.9, theta_edge = 0.98)

  expect_equal(nrow(result), 1L)
  expect_equal(result$prediction, "G")
  expect_equal(result$centr_signal_pos, 525)
  expect_equal(result$score, 0.99)
})

test_that("call_nonA_cdna decodes all three residues", {
  centres <- seq(50, 1050, by = 25)
  nonA_rows <- which(centres %in% c(500, 525))
  for (class_column in 2:4) {
    result <- call_nonA_cdna(centres, create_probabilities(length(centres), nonA_rows, class_column),
                             tail_start = 100, tail_end = 1000, theta_mid = 0.9, theta_edge = 0.9)
    expect_equal(result$prediction, c("C", "G", "U")[class_column - 1])
  }
})

test_that("call_nonA_cdna ignores runs shorter than min_run", {
  centres <- seq(50, 1050, by = 25)
  probabilities <- create_probabilities(length(centres), which(centres == 500))

  result <- call_nonA_cdna(centres, probabilities, tail_start = 100, tail_end = 1000,
                           theta_mid = 0.9, theta_edge = 0.9, min_run = 2)
  expect_equal(nrow(result), 0L)

  result <- call_nonA_cdna(centres, probabilities, tail_start = 100, tail_end = 1000,
                           theta_mid = 0.9, theta_edge = 0.9, min_run = 1)
  expect_equal(nrow(result), 1L)
})

test_that("call_nonA_cdna applies theta_edge near the tail boundaries", {
  centres <- seq(50, 1050, by = 25)
  # windows within 100 samples of the tail start
  nonA_rows <- which(centres %in% c(150, 175))
  probabilities <- create_probabilities(length(centres), nonA_rows, p = 0.95)

  strict <- call_nonA_cdna(centres, probabilities, tail_start = 100, tail_end = 1000,
                           theta_mid = 0.9, theta_edge = 0.98)
  lenient <- call_nonA_cdna(centres, probabilities, tail_start = 100, tail_end = 1000,
                            theta_mid = 0.9, theta_edge = 0.9)

  expect_equal(nrow(strict), 0L)
  expect_equal(nrow(lenient), 1L)
})

test_that("call_nonA_cdna drops calls positioned outside the tail", {
  centres <- seq(50, 1050, by = 25)
  nonA_rows <- which(centres %in% c(50, 75))
  result <- call_nonA_cdna(centres, create_probabilities(length(centres), nonA_rows),
                           tail_start = 100, tail_end = 1000, theta_mid = 0.9, theta_edge = 0.9)
  expect_equal(nrow(result), 0L)
})

test_that("call_nonA_cdna reports separate runs as separate calls", {
  centres <- seq(50, 1050, by = 25)
  nonA_rows <- which(centres %in% c(300, 325, 700, 725))
  result <- call_nonA_cdna(centres, create_probabilities(length(centres), nonA_rows),
                           tail_start = 100, tail_end = 1000, theta_mid = 0.9, theta_edge = 0.9)
  expect_equal(nrow(result), 2L)
})


################################################################################
# extract_tail_signals_cdna
################################################################################

test_that("extract_tail_signals_cdna errors on missing arguments", {
  expect_error(extract_tail_signals_cdna(pod5_dir = tempdir()),
               "Tail coordinates are missing")
  expect_error(extract_tail_signals_cdna(tail_coordinates = data.frame(read_id = "r1")),
               "Directory with POD5 files is missing")
})

test_that("extract_tail_signals_cdna errors on empty data frame", {
  expect_error(extract_tail_signals_cdna(tail_coordinates = data.frame(), pod5_dir = tempdir()),
               "Empty data frame provided as an input")
})

test_that("extract_tail_signals_cdna errors on non-numeric num_cores", {
  expect_error(extract_tail_signals_cdna(tail_coordinates = data.frame(read_id = "r1"),
                                         pod5_dir = tempdir(), num_cores = "two"))
})

test_that("extract_tail_signals_cdna errors when required columns are missing", {
  expect_error(extract_tail_signals_cdna(tail_coordinates = data.frame(read_id = "r1"),
                                         pod5_dir = tempdir()),
               "Required columns missing")
})

test_that("extract_tail_signals_cdna returns flanked raw signals for the shipped data", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("reticulate")
  skip_if(nchar(cdna_test_dir()) == 0, "cDNA test data not found in installed package")
  skip_if_not(reticulate::py_module_available("pod5"), "Python module pod5 not available")

  tail_coordinates <- read_tail_coordinates_cdna(file.path(cdna_test_dir(), "lachesis_polya.tsv"),
                                                 source = "lachesis")
  tail_coordinates <- tail_coordinates[tail_coordinates$qc_tag == "PASS", ][1:4, ]

  signal_list <- extract_tail_signals_cdna(tail_coordinates, pod5_dir = cdna_test_dir())

  expect_type(signal_list, "list")
  expect_true(all(tail_coordinates$read_id %in% names(signal_list)))
  # tail plus 100 samples of flank on each side
  expected_length <- tail_coordinates$tail_end - tail_coordinates$tail_start + 200
  expect_equal(as.numeric(lengths(signal_list[tail_coordinates$read_id])), expected_length)
})


################################################################################
# predict_nonA_cdna
################################################################################

test_that("predict_nonA_cdna errors on missing arguments", {
  tail_coordinates <- data.frame(read_id = "r1", tail_start = 100, tail_end = 800, polya_length = 60)
  expect_error(predict_nonA_cdna(tail_coordinates = tail_coordinates, tail_type = "polyA"),
               "Signal list is missing")
  expect_error(predict_nonA_cdna(signal_list = list(), tail_type = "polyA"),
               "Tail coordinates are missing")
  expect_error(predict_nonA_cdna(signal_list = list(), tail_coordinates = tail_coordinates),
               "Tail type is missing")
})

test_that("predict_nonA_cdna errors on malformed inputs", {
  tail_coordinates <- data.frame(read_id = "r1", tail_start = 100, tail_end = 800, polya_length = 60)
  expect_error(predict_nonA_cdna(signal_list = "not_a_list", tail_coordinates = tail_coordinates,
                                 tail_type = "polyA"),
               "Signal list must be a list")
  expect_error(predict_nonA_cdna(signal_list = list(r1 = 1:10), tail_coordinates = data.frame(read_id = "r1"),
                                 tail_type = "polyA"),
               "Required columns missing")
  expect_error(predict_nonA_cdna(signal_list = list(r1 = 1:10), tail_coordinates = tail_coordinates,
                                 tail_type = "polyG"),
               "Wrong tail_type provided")
})

test_that("predict_nonA_cdna classifies the shipped poly(A) and poly(T) reads", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("keras")
  skip_if(nchar(cdna_test_dir()) == 0, "cDNA test data not found in installed package")
  skip_if_not(reticulate::py_module_available("pod5"), "Python module pod5 not available")
  skip_if_not(reticulate::py_module_available("tensorflow"), "Python module tensorflow not available")

  tail_coordinates <- read_tail_coordinates_cdna(file.path(cdna_test_dir(), "lachesis_polya.tsv"),
                                                 source = "lachesis")
  tail_coordinates <- tail_coordinates[tail_coordinates$qc_tag == "PASS", ]
  signal_list <- extract_tail_signals_cdna(tail_coordinates, pod5_dir = cdna_test_dir())

  for (tail_type in c("polyA", "polyT")) {
    coordinates <- tail_coordinates[tail_coordinates$tail_type == tail_type, ]
    result <- predict_nonA_cdna(signal_list = signal_list[coordinates$read_id],
                                tail_coordinates = coordinates,
                                tail_type = tail_type)

    expect_true(is.data.frame(result))
    expect_true(all(c("read_id", "prediction", "est_nonA_pos", "score") %in% colnames(result)))
    expect_true(all(result$prediction %in% c("C", "G", "U")))
    expect_true(all(result$read_id %in% coordinates$read_id))
    # the spike-ins carry a non-A in every tail
    expect_gte(length(unique(result$read_id)), 10)
    # positions lie within the tail
    polya_length <- coordinates$polya_length[match(result$read_id, coordinates$read_id)]
    expect_true(all(result$est_nonA_pos >= 0 & result$est_nonA_pos <= polya_length))
  }
})
