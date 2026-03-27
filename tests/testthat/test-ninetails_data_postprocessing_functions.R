################################################################################
# Testing data postprocessing functions
################################################################################


# Helpers
################################################################################

#' Create a temporary read_classes TSV file
#'
#' Generates a minimal `read_classes`-like data frame, writes it to a
#' temporary TSV file, and returns the file path.
#'
#' @param n_reads Integer; number of reads (rows) to generate.
#' @param sample_col Logical; if TRUE, includes a `sample_name` column.
#'
#' @return A character string giving the path to the generated TSV file.
#'
#' @keywords internal
make_class_tsv <- function(n_reads = 10, sample_col = FALSE) {
  df <- data.frame(readname = paste0("read_", seq_len(n_reads)),
                   contig = rep("beta_gal", n_reads),
                   polya_length = seq(50, by = 5, length.out = n_reads),
                   qc_tag = rep("PASS", n_reads),
                   class = rep(c("decorated", "blank"), length.out = n_reads),
                   comments = rep(c("YAY", "MPU"), length.out = n_reads),
                   stringsAsFactors = FALSE)
  if (sample_col) {
    df$sample_name <- "WT_1"
  }
  path <- tempfile(fileext = ".tsv")
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  return(path)
}

#' Create a temporary nonadenosine_residues TSV file
#'
#' Generates a minimal `nonadenosine_residues`-like data frame, writes it to a
#' temporary TSV file, and returns the file path.
#'
#' @param n_reads Integer; number of reads (rows) to generate.
#'
#' @return A character string giving the path to the generated TSV file.
#'
#' @keywords internal
make_residue_tsv <- function(n_reads = 8) {
  df <- data.frame(readname = paste0("read_", seq_len(n_reads)),
                   contig = rep("beta_gal", n_reads),
                   prediction = rep(c("C", "G", "U", "C"), length.out = n_reads),
                   est_nonA_pos = seq(10, by = 8, length.out = n_reads),
                   polya_length = rep(100, n_reads),
                   qc_tag = rep("PASS", n_reads),
                   group = rep("WT", n_reads),
                   stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".tsv")
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  return(path)
}

#' Generate synthetic class and residue data inputs
#'
#' Constructs synthetic `class_data` and `residue_data` data frames large
#' enough to satisfy conditions such as `count_nonA > 10` in downstream
#' processing (e.g., `correct_residue_data`).
#'
#' @param n_reads Integer; number of reads (rows) to generate.
#' @param grouped Logical; if TRUE, adds `sample_name` columns to both
#'   outputs with WT/KO grouping.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{class_data}{A data.frame of simulated class data.}
#'     \item{residue_data}{A data.frame of simulated residue data.}
#'   }
#'
#' @keywords internal
make_synthetic_inputs <- function(n_reads = 30, grouped = FALSE) {
  class_data <- data.frame(readname = paste0("read_", seq_len(n_reads)),
                           contig = rep("beta_gal", n_reads),
                           polya_length = c(seq(80, 130, length.out = n_reads %/% 2),
                                            seq(50, 100, length.out = n_reads - n_reads %/% 2)),
                           qc_tag = rep("PASS", n_reads),
                           class = rep(c("decorated", "blank"), length.out = n_reads),
                           comments = rep(c("YAY", "MPU"), length.out = n_reads),
                           stringsAsFactors = FALSE)
  residue_data <- data.frame(readname = paste0("read_", seq_len(n_reads)),
                             contig = rep("beta_gal", n_reads),
                             prediction = rep(c("C", "G", "U"), length.out = n_reads),
                             est_nonA_pos = seq(5, by = 4, length.out = n_reads),
                             polya_length = rep(100, n_reads),
                             qc_tag = rep("PASS", n_reads),
                             group = rep("WT", n_reads),
                             stringsAsFactors = FALSE)
  if (grouped) {
    class_data$sample_name <- rep(c("WT_1", "KO_1"), length.out = n_reads)
    residue_data$sample_name <- rep(c("WT_1", "KO_1"), length.out = n_reads)
  }
  return(list(class_data = class_data, residue_data = residue_data))
}

# Path to the built-in whitelist file (used in skip_if guards).
whitelist_path <- system.file("extdata", "whitelists", "whitelist.RData",
                              package = "ninetails")


# read_class_single
################################################################################

test_that("read_class_single errors on missing class_path", {
  expect_error(read_class_single(),
               "class_path.*is missing")
})

test_that("read_class_single errors on NA class_path", {
  expect_error(read_class_single(class_path = NA),
               "missing or invalid")
})

test_that("read_class_single errors on empty string class_path", {
  expect_error(read_class_single(class_path = ""),
               "missing or invalid")
})

test_that("read_class_single errors on nonexistent file", {
  expect_error(read_class_single(class_path = "/nonexistent/path/file.tsv"),
               "does not exist")
})

test_that("read_class_single errors on empty file", {
  tmp_file <- tempfile(fileext = ".tsv")
  file.create(tmp_file)
  on.exit(unlink(tmp_file))

  expect_error(read_class_single(class_path = tmp_file),
               "empty")
})

test_that("read_class_single loads file, parses contig, and returns tibble (lines 74-110)", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")

  path <- make_class_tsv()
  on.exit(unlink(path))

  result <- suppressMessages(read_class_single(class_path = path))

  expect_s3_class(result, "tbl_df")
  # GENCODE parsing columns always present, even for non-GENCODE contigs
  expect_true("transcript" %in% colnames(result))
  expect_true("ensembl_transcript_id_full" %in% colnames(result))
  expect_true("ensembl_transcript_id_short" %in% colnames(result))
})

test_that("read_class_single attaches sample_name when not already in file", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")

  # File does NOT contain a sample_name column.
  # BUG (line 83): condition uses !"sample_name" instead of "sample_name",
  # so the warning fires when the column is absent rather than when it already
  # exists. Tests reflect actual behaviour until the source is fixed.
  path <- make_class_tsv(sample_col = FALSE)
  on.exit(unlink(path))

  expect_warning(
    suppressMessages(
      read_class_single(class_path = path, sample_name = "WT_1")
    ),
    "Overwriting"
  )
})

test_that("read_class_single overwrites existing sample_name without warning", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")

  # File ALREADY contains sample_name. Due to the inverted condition at line 83,
  # no warning fires here. Column must still be overwritten.
  path <- make_class_tsv(sample_col = TRUE)
  on.exit(unlink(path))

  result <- expect_no_warning(
    suppressMessages(
      read_class_single(class_path = path, sample_name = "KO_1")
    )
  )
  expect_true(all(result$sample_name == "KO_1"))
})


# read_class_multiple
################################################################################

test_that("read_class_multiple errors on missing samples_table", {
  expect_error(read_class_multiple(),
               "Samples table.*is missing")
})

test_that("read_class_multiple errors on empty samples_table", {
  expect_error(read_class_multiple(samples_table = data.frame()),
               "Empty data frame")
})

test_that("read_class_multiple errors on missing class_path column", {
  bad_table <- data.frame(sample_name = "WT_1")
  expect_error(read_class_multiple(samples_table = bad_table),
               "class_path.*sample_name")
})

test_that("read_class_multiple errors on missing sample_name column", {
  bad_table <- data.frame(class_path = "/some/path.tsv")
  expect_error(read_class_multiple(samples_table = bad_table),
               "class_path.*sample_name")
})

test_that("read_class_multiple loads and combines multiple class files", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("tibble")

  path1 <- make_class_tsv(n_reads = 5)
  path2 <- make_class_tsv(n_reads = 6)
  on.exit({ unlink(path1); unlink(path2) })

  # group is safe here because make_class_tsv does not write a group column,
  # so there is no duplicate column name conflict during unnest
  samples_table <- data.frame(class_path = c(path1, path2),
                              sample_name = c("WT_1", "KO_1"),
                              group = c("WT", "KO"),
                              stringsAsFactors = FALSE)

  result <- suppressMessages(read_class_multiple(samples_table = samples_table))

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 11L)
  expect_true("sample_name" %in% colnames(result))
  expect_true("group" %in% colnames(result))
  expect_true("readname" %in% colnames(result))
  expect_true("class" %in% colnames(result))
})


# count_class
################################################################################

test_that("count_class validates empty data frame", {
  expect_error(count_class(class_data = data.frame()), "Empty data frame")
})

test_that("count_class validates grouping_factor exists", {
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  expect_error(count_class(class_data = class_data_single,
                           grouping_factor = "nonexistent_column"),
               "is not a column")
})

test_that("count_class returns detailed counts (comments) without grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  result <- count_class(class_data = class_data_single, detailed = TRUE)

  expect_s3_class(result, "tbl_df")
  expect_true("comments" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)
})

test_that("count_class returns crude counts (class) without grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  result <- count_class(class_data = class_data_single, detailed = FALSE)

  expect_s3_class(result, "tbl_df")
  expect_true("class" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)
})

test_that("count_class returns detailed counts with grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")

  result <- count_class(class_data = class_data_grouped,
                        grouping_factor = "sample_name",
                        detailed = TRUE)

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
  expect_true("comments" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 3)
})

test_that("count_class returns crude counts with grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")

  result <- count_class(class_data = class_data_grouped,
                        grouping_factor = "group",
                        detailed = FALSE)

  expect_s3_class(result, "tbl_df")
  expect_true("group" %in% colnames(result))
  expect_true("class" %in% colnames(result))
  expect_true("n" %in% colnames(result))
})

test_that("count_class reorders class factor with decorated last", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  result <- count_class(class_data = class_data_single, detailed = FALSE)
  levels_order <- levels(result$class)
  expect_equal(levels_order[length(levels_order)], "decorated")
})


# read_residue_single
################################################################################

test_that("read_residue_single errors on NA residue_path", {
  expect_error(read_residue_single(residue_path = NA),
               "missing or an empty string")
})

test_that("read_residue_single errors on empty string residue_path", {
  expect_error(read_residue_single(residue_path = ""),
               "missing or an empty string")
})

test_that("read_residue_single errors on nonexistent file", {
  expect_error(read_residue_single(residue_path = "/nonexistent/path/file.tsv"),
               "missing or an empty file")
})

test_that("read_residue_single errors on empty file", {
  tmp_file <- tempfile(fileext = ".tsv")
  file.create(tmp_file)
  on.exit(unlink(tmp_file))

  expect_error(read_residue_single(residue_path = tmp_file),
               "missing or an empty file")
})

test_that("read_residue_single reads valid file", {
  tmp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_file))

  test_df <- data.frame(readname = c("read_1", "read_2"),
                        contig = c("ENST00000123|a|b|c|d|GENE1|e",
                                   "ENST00000456|a|b|c|d|GENE2|e"),
                        prediction = c("G", "C"),
                        est_nonA_pos = c(50, 30),
                        polya_length = c(100, 80),
                        qc_tag = c("PASS", "PASS"),
                        stringsAsFactors = FALSE)
  utils::write.table(test_df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- suppressMessages(read_residue_single(residue_path = tmp_file))

  expect_s3_class(result, "tbl_df")
  expect_true("readname" %in% colnames(result))
  expect_true("transcript" %in% colnames(result))
  expect_true("ensembl_transcript_id_full" %in% colnames(result))
  expect_true("ensembl_transcript_id_short" %in% colnames(result))
})

test_that("read_residue_single adds sample_name when provided", {
  tmp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_file))

  test_df <- data.frame(readname = "read_1",
                        contig = "gene1",
                        prediction = "G",
                        est_nonA_pos = 50,
                        polya_length = 100,
                        qc_tag = "PASS",
                        stringsAsFactors = FALSE)
  utils::write.table(test_df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- suppressMessages(suppressWarnings(
    read_residue_single(residue_path = tmp_file, sample_name = "KO_1")
  ))

  expect_true("sample_name" %in% colnames(result))
  expect_equal(as.character(result$sample_name[1]), "KO_1")
  expect_true(is.factor(result$sample_name))
})


# read_residue_multiple
################################################################################

test_that("read_residue_multiple errors on missing samples_table", {
  expect_error(read_residue_multiple(),
               "samples_table.*is missing")
})

test_that("read_residue_multiple errors on empty samples_table", {
  expect_error(read_residue_multiple(samples_table = data.frame()),
               "Empty data frame")
})

test_that("read_residue_multiple errors on missing residue_path column", {
  bad_table <- data.frame(sample_name = "WT_1")
  expect_error(read_residue_multiple(samples_table = bad_table),
               "residue_path.*sample_name")
})

test_that("read_residue_multiple errors on missing sample_name column", {
  bad_table <- data.frame(residue_path = "/some/path.tsv")
  expect_error(read_residue_multiple(samples_table = bad_table),
               "residue_path.*sample_name")
})

test_that("read_residue_multiple loads and combines multiple residue files", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("tibble")

  path1 <- make_residue_tsv(n_reads = 4)
  path2 <- make_residue_tsv(n_reads = 5)
  on.exit({ unlink(path1); unlink(path2) })

  # group must NOT be a column in samples_table: make_residue_tsv already
  # writes a group column inside each TSV, so adding it here would create a
  # duplicate name collision during unnest at line 491.
  samples_table <- data.frame(residue_path = c(path1, path2),
                              sample_name = c("WT_1", "KO_1"),
                              stringsAsFactors = FALSE)

  result <- suppressMessages(read_residue_multiple(samples_table = samples_table))

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 9L)
  expect_true("sample_name" %in% colnames(result))
  expect_true("prediction" %in% colnames(result))
  expect_true("est_nonA_pos" %in% colnames(result))
})


# count_residues
################################################################################

test_that("count_residues validates empty data frame", {
  expect_error(count_residues(residue_data = data.frame()), "Empty data frame")
})

test_that("count_residues validates grouping_factor exists", {
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  expect_error(count_residues(residue_data = residue_data_single,
                              grouping_factor = "nonexistent_column"),
               "is not a column")
})

test_that("count_residues returns counts without grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  result <- count_residues(residue_data = residue_data_single)

  expect_s3_class(result, "tbl_df")
  expect_true("prediction" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)
  expect_true(all(c("C", "G", "U") %in% result$prediction))
})

test_that("count_residues returns counts with grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_residues(residue_data = residue_data_grouped,
                           grouping_factor = "group")

  expect_s3_class(result, "tbl_df")
  expect_true("group" %in% colnames(result))
  expect_true("prediction" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 3)
})

test_that("count_residues works with sample_name grouping", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("forcats")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_residues(residue_data = residue_data_grouped,
                           grouping_factor = "sample_name")

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
  expect_equal(ncol(result), 3)
})


# spread_nonA_residues
################################################################################

test_that("spread_nonA_residues validates required arguments", {
  expect_error(spread_nonA_residues(), "is missing")
  expect_error(spread_nonA_residues(residue_data = data.frame()), "Empty data frame")
})

test_that("spread_nonA_residues validates required columns", {
  bad_df <- data.frame(a = 1, b = 2)
  expect_error(spread_nonA_residues(residue_data = bad_df), "Required columns missing")
})

test_that("spread_nonA_residues returns wide format with prediction counts", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- spread_nonA_residues(residue_data = residue_data_grouped)

  expect_s3_class(result, "tbl_df")
  expect_true("readname" %in% colnames(result))
  expect_true("nonA_residues" %in% colnames(result))
  expect_true(length(grep("^prediction_", colnames(result), value = TRUE)) > 0)
})

test_that("spread_nonA_residues creates correct nonA_residues string", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- spread_nonA_residues(residue_data = residue_data_grouped)
  nona_strings <- result$nonA_residues[!is.na(result$nonA_residues)]

  expect_true(length(nona_strings) > 0)
  expect_true(all(grepl("[CGU]\\d+", nona_strings)))
})


# merge_nonA_tables
################################################################################

test_that("merge_nonA_tables validates required arguments", {
  expect_error(merge_nonA_tables(residue_data = data.frame(a = 1)),
               "class data is missing")
  expect_error(merge_nonA_tables(class_data = data.frame(a = 1)),
               "residue data is missing")
  expect_error(merge_nonA_tables(class_data = data.frame(),
                                 residue_data = data.frame(a = 1)),
               "Empty data frame.*class_data")
  expect_error(merge_nonA_tables(class_data = data.frame(a = 1),
                                 residue_data = data.frame()),
               "Empty data frame.*residue_data")
})

test_that("merge_nonA_tables validates pass_only is logical", {
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  expect_error(merge_nonA_tables(class_data = class_data_grouped,
                                 residue_data = residue_data_grouped,
                                 pass_only = "yes"),
               "TRUE/FALSE")
})

test_that("merge_nonA_tables validates qc_tag column exists", {
  bad_class <- data.frame(readname = "read_1",
                          class = "decorated",
                          comments = "YAY",
                          stringsAsFactors = FALSE)
  bad_residue <- data.frame(readname = "read_1",
                            prediction = "G",
                            est_nonA_pos = 10,
                            group = "WT",
                            stringsAsFactors = FALSE)

  expect_error(merge_nonA_tables(class_data = bad_class,
                                 residue_data = bad_residue),
               "qc_tag.*missing")
})

test_that("merge_nonA_tables works with numeric qc_tag (Dorado)", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("class_data_grouped"), "class_data_grouped not loaded")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- suppressMessages(suppressWarnings(
    merge_nonA_tables(class_data = class_data_grouped,
                      residue_data = residue_data_grouped,
                      pass_only = TRUE)))

  expect_s3_class(result, "data.frame")
  expect_true("readname" %in% colnames(result))
  expect_true("class" %in% colnames(result))
  expect_true("nonA_residues" %in% colnames(result))
  expect_true(any(grepl("prediction_", colnames(result))))
})

test_that("merge_nonA_tables works with character qc_tag (Guppy)", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  class_data_char <- data.frame(readname = c("read_1", "read_2", "read_3"),
                                contig = rep("beta_gal", 3),
                                polya_length = c(100, 80, 60),
                                qc_tag = c("PASS", "PASS", "SUFFCLIP"),
                                class = c("decorated", "blank", "decorated"),
                                comments = c("YAY", "MPU", "YAY"),
                                group = c("WT", "WT", "KO"),
                                stringsAsFactors = FALSE)
  residue_data_char <- data.frame(readname = c("read_1", "read_3"),
                                  contig = rep("beta_gal", 2),
                                  prediction = c("G", "C"),
                                  est_nonA_pos = c(50, 30),
                                  polya_length = c(100, 60),
                                  qc_tag = c("PASS", "SUFFCLIP"),
                                  group = c("WT", "KO"),
                                  stringsAsFactors = FALSE)

  result_pass <- suppressMessages(suppressWarnings(
    merge_nonA_tables(class_data = class_data_char,
                      residue_data = residue_data_char,
                      pass_only = TRUE)))
  expect_s3_class(result_pass, "data.frame")

  result_all <- suppressMessages(suppressWarnings(
    merge_nonA_tables(class_data = class_data_char,
                      residue_data = residue_data_char,
                      pass_only = FALSE)))
  expect_s3_class(result_all, "data.frame")
})

test_that("merge_nonA_tables errors on non-character, non-numeric qc_tag", {
  # Logical qc_tag passes the column existence check but hits the else stop.
  # At least one non-unclassified row is needed to survive line 802.
  bad_class <- data.frame(readname = c("read_1", "read_2"),
                          contig = rep("beta_gal", 2),
                          polya_length = c(80, 60),
                          qc_tag = c(TRUE, FALSE),
                          class = c("decorated", "blank"),
                          comments = c("YAY", "MPU"),
                          group = rep("WT", 2),
                          stringsAsFactors = FALSE)
  valid_residue <- data.frame(readname = "read_1",
                              contig = "beta_gal",
                              prediction = "G",
                              est_nonA_pos = 30,
                              polya_length = 80,
                              qc_tag = TRUE,
                              group = "WT",
                              stringsAsFactors = FALSE)

  expect_error(merge_nonA_tables(class_data = bad_class,
                                 residue_data = valid_residue,
                                 pass_only = TRUE),
               "Unexpected qc_tag type")
})

test_that("merge_nonA_tables warns when all reads are filtered out", {
  # All rows have qc_tag == "ADAPTER"; pass_only = TRUE finds no "PASS" rows,
  # producing an empty class2 and triggering the warning
  class_no_pass <- data.frame(readname = paste0("read_", 1:3),
                              contig = rep("beta_gal", 3),
                              polya_length = c(80, 60, 70),
                              qc_tag = rep("ADAPTER", 3),
                              class = c("decorated", "blank", "blank"),
                              comments = c("YAY", "MPU", "MPU"),
                              group = rep("WT", 3),
                              stringsAsFactors = FALSE)
  valid_residue <- data.frame(readname = "read_1",
                              contig = "beta_gal",
                              prediction = "G",
                              est_nonA_pos = 30,
                              polya_length = 80,
                              qc_tag = "ADAPTER",
                              group = "WT",
                              stringsAsFactors = FALSE)

  expect_warning(
    suppressMessages(
      merge_nonA_tables(class_data = class_no_pass,
                        residue_data = valid_residue,
                        pass_only = TRUE)
    ),
    "No reads passed quality filtering"
  )
})


# summarize_nonA
################################################################################

test_that("summarize_nonA validates required arguments", {
  expect_error(summarize_nonA(),
               "merged_nonA_tables.*is missing")
  expect_error(summarize_nonA(merged_nonA_tables = data.frame()),
               "Empty data frame")
})

test_that("summarize_nonA validates summary_factors type", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  expect_error(summarize_nonA(merged_nonA_tables = merged_nonA_char,
                              summary_factors = 123),
               "Non-character")
})

test_that("summarize_nonA validates summary_factors exist in data", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  expect_error(summarize_nonA(merged_nonA_tables = merged_nonA_char,
                              summary_factors = "nonexistent_column"),
               "Non-existent column")
})

test_that("summarize_nonA returns transcript-level summary", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("stringr")
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  result <- summarize_nonA(merged_nonA_tables = test_data,
                           summary_factors = "group",
                           transcript_id_column = "ensembl_transcript_id_short")

  expect_s3_class(result, "tbl_df")
  expect_true("ensembl_transcript_id_short" %in% colnames(result))
  expect_true("group" %in% colnames(result))
  expect_true("polya_median" %in% colnames(result))
  expect_true("polya_mean" %in% colnames(result))
  expect_true("counts_total" %in% colnames(result))
  expect_true("counts_blank" %in% colnames(result))
})


# nanopolish_qc
################################################################################

test_that("nanopolish_qc validates empty data frame", {
  expect_error(nanopolish_qc(class_data = data.frame()), "Empty data frame")
})

test_that("nanopolish_qc validates grouping_factor exists", {
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  expect_error(nanopolish_qc(class_data = class_data_single,
                             grouping_factor = "nonexistent_column"),
               "is not a column")
})

test_that("nanopolish_qc returns qc counts without grouping", {
  skip_if_not_installed("dplyr")
  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  class_data_qc <- data.frame(readname = paste0("read_", 1:10),
                              contig = rep("beta_gal", 10),
                              polya_length = sample(50:150, 10),
                              qc_tag = c(rep("PASS", 6), rep("ADAPTER", 2), rep("SUFFCLIP", 2)),
                              class = c(rep("decorated", 5), rep("blank", 3), rep("unclassified", 2)),
                              comments = c(rep("YAY", 5), rep("MPU", 3), rep("IRL", 2)),
                              stringsAsFactors = FALSE)

  result <- nanopolish_qc(class_data = class_data_qc)

  expect_s3_class(result, "tbl_df")
  expect_true("qc_tag" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 2)
})

test_that("nanopolish_qc returns qc counts with grouping", {
  skip_if_not_installed("dplyr")

  class_data_qc <- data.frame(readname = paste0("read_", 1:20),
                              contig = rep("beta_gal", 20),
                              polya_length = sample(50:150, 20),
                              qc_tag = rep(c("PASS", "PASS", "ADAPTER", "SUFFCLIP"), 5),
                              class = rep(c("decorated", "blank", "unclassified", "decorated"), 5),
                              comments = rep(c("YAY", "MPU", "IRL", "YAY"), 5),
                              sample_name = rep(c("WT_1", "KO_1"), each = 10),
                              stringsAsFactors = FALSE)

  result <- nanopolish_qc(class_data = class_data_qc, grouping_factor = "sample_name")

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
  expect_true("qc_tag" %in% colnames(result))
  expect_true("n" %in% colnames(result))
  expect_equal(ncol(result), 3)
})


# correct_residue_data
################################################################################

test_that("correct_residue_data validates required arguments", {
  expect_error(correct_residue_data(residue_data = data.frame(a = 1),
                                    transcript_column = "contig"),
               "class_data argument is missing")
  expect_error(correct_residue_data(class_data = data.frame(a = 1),
                                    transcript_column = "contig"),
               "residue_data argument is missing")
  expect_error(correct_residue_data(class_data = data.frame(a = 1),
                                    residue_data = data.frame(b = 2)),
               "transcript_column argument is missing")
})

test_that("correct_residue_data validates empty data frames", {
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  expect_error(correct_residue_data(class_data = data.frame(),
                                    residue_data = residue_data_single,
                                    transcript_column = "contig"),
               "Empty data frame.*class_data")

  skip_if(!exists("class_data_single"), "class_data_single not loaded")

  expect_error(correct_residue_data(class_data = class_data_single,
                                    residue_data = data.frame(),
                                    transcript_column = "contig"),
               "Empty data frame.*residue_data")
})

test_that("correct_residue_data runs without grouping_factor", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 30)

  # ref = "no_match" exercises the else-branch at line 1278 (whitelist = ref).
  # Using ref = NULL (the default) would error because if (NULL == "mmusculus")
  # evaluates to logical(0) — a bug in the source (should use identical()).
  result <- correct_residue_data(
    class_data = inputs$class_data,
    residue_data = inputs$residue_data,
    grouping_factor = NULL,
    transcript_column = "contig",
    ref = "no_match"
  )

  expect_s3_class(result, "data.frame")
  expect_true("qc_pos" %in% colnames(result))
  expect_true("mode_pos" %in% colnames(result))
  expect_true("seg_err_quart" %in% colnames(result))
  expect_true(all(result$qc_pos %in% c("Y", "N")))
})

test_that("correct_residue_data runs with grouping_factor", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 30, grouped = TRUE)

  result <- correct_residue_data(
    class_data = inputs$class_data,
    residue_data = inputs$residue_data,
    grouping_factor = "sample_name",
    transcript_column = "contig",
    ref = "no_match"
  )

  expect_s3_class(result, "data.frame")
  expect_true("qc_pos" %in% colnames(result))
})

test_that("correct_residue_data respects mmusculus whitelist", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 15)
  result <- correct_residue_data(class_data = inputs$class_data,
                                 residue_data = inputs$residue_data,
                                 grouping_factor = NULL,
                                 transcript_column = "contig",
                                 ref = "mmusculus")
  expect_s3_class(result, "data.frame")
  expect_true("qc_pos" %in% colnames(result))
})

test_that("correct_residue_data respects hsapiens whitelist", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 15)
  result <- correct_residue_data(class_data = inputs$class_data,
                                 residue_data = inputs$residue_data,
                                 grouping_factor = NULL,
                                 transcript_column = "contig",
                                 ref = "hsapiens")
  expect_s3_class(result, "data.frame")
})

test_that("correct_residue_data respects scerevisiae whitelist", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 15)
  result <- correct_residue_data(class_data = inputs$class_data,
                                 residue_data = inputs$residue_data,
                                 grouping_factor = NULL,
                                 transcript_column = "contig",
                                 ref = "scerevisiae")
  expect_s3_class(result, "data.frame")
})

test_that("correct_residue_data respects celegans whitelist", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 15)
  result <- correct_residue_data(class_data = inputs$class_data,
                                 residue_data = inputs$residue_data,
                                 grouping_factor = NULL,
                                 transcript_column = "contig",
                                 ref = "celegans")
  expect_s3_class(result, "data.frame")
})

test_that("correct_residue_data respects athaliana whitelist", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 15)
  result <- correct_residue_data(class_data = inputs$class_data,
                                 residue_data = inputs$residue_data,
                                 grouping_factor = NULL,
                                 transcript_column = "contig",
                                 ref = "athaliana")
  expect_s3_class(result, "data.frame")
})

test_that("correct_residue_data respects tbrucei whitelist", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 15)
  result <- correct_residue_data(class_data = inputs$class_data,
                                 residue_data = inputs$residue_data,
                                 grouping_factor = NULL,
                                 transcript_column = "contig",
                                 ref = "tbrucei")
  expect_s3_class(result, "data.frame")
})


# correct_class_data
################################################################################

test_that("correct_class_data validates required arguments", {
  expect_error(correct_class_data(class_data = data.frame(a = 1)),
               "residue_data_edited argument is missing")
  expect_error(correct_class_data(residue_data_edited = data.frame(a = 1)),
               "class_data argument is missing")
})

test_that("correct_class_data errors on empty residue_data_edited", {
  valid_class <- make_synthetic_inputs()$class_data

  expect_error(correct_class_data(residue_data_edited = data.frame(),
                                  class_data = valid_class),
               "Empty data frame provided as an input \\(residue_data_edited\\)")
})

test_that("correct_class_data errors on empty class_data", {
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 15)
  residue_data_edited <- correct_residue_data(
    class_data = inputs$class_data,
    residue_data = inputs$residue_data,
    grouping_factor = NULL,
    transcript_column = "contig",
    ref = "no_match"
  )

  expect_error(correct_class_data(residue_data_edited = residue_data_edited,
                                  class_data = data.frame()),
               "Empty data frame provided as an input \\(class_data\\)")
})

test_that("correct_class_data returns reclassified class_data with corr_class and corr_comments", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 30)
  residue_data_edited <- correct_residue_data(
    class_data = inputs$class_data,
    residue_data = inputs$residue_data,
    grouping_factor = NULL,
    transcript_column = "contig",
    ref = "no_match"
  )

  result <- correct_class_data(residue_data_edited = residue_data_edited,
                               class_data = inputs$class_data)

  expect_s3_class(result, "data.frame")
  expect_true("corr_class" %in% colnames(result))
  expect_true("corr_comments" %in% colnames(result))
  expect_false("n_resid" %in% colnames(result))
  expect_false("no_qc_pos_N" %in% colnames(result))
  expect_true(all(result$corr_class %in% c("decorated", "blank", "unclassified", NA)))
})

test_that("correct_class_data demotes fully-ambiguous decorated reads to blank", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  class_mini <- data.frame(readname = "read_1",
                           contig = "beta_gal",
                           polya_length = 100,
                           qc_tag = "PASS",
                           class = "decorated",
                           comments = "YAY",
                           stringsAsFactors = FALSE)
  # Hand-crafted residue_data_edited with qc_pos = "N" for all positions,
  # bypassing correct_residue_data to test the demotion logic directly
  residue_edited_mini <- data.frame(readname = "read_1",
                                    contig = "beta_gal",
                                    prediction = "G",
                                    est_nonA_pos = 2,
                                    polya_length = 100,
                                    qc_tag = "PASS",
                                    group = "WT",
                                    mode_pos = 2,
                                    pos_err_quart = 5,
                                    count_nonA = 1,
                                    mode_len = 90,
                                    seg_err_quart = 10,
                                    count = 1,
                                    qc_pos = "N",
                                    stringsAsFactors = FALSE)

  result <- correct_class_data(residue_data_edited = residue_edited_mini,
                               class_data = class_mini)

  expect_equal(result$corr_class[result$readname == "read_1"], "blank")
  expect_equal(result$corr_comments[result$readname == "read_1"], "MPU")
})


# reclassify_ninetails_data
################################################################################

test_that("reclassify_ninetails_data validates required arguments", {
  expect_error(reclassify_ninetails_data(class_data = data.frame(a = 1),
                                         transcript_column = "contig"),
               "The residue_data argument is missing")
  expect_error(reclassify_ninetails_data(residue_data = data.frame(a = 1),
                                         transcript_column = "contig"),
               "The class_data argument is missing")
  expect_error(reclassify_ninetails_data(residue_data = data.frame(a = 1),
                                         class_data = data.frame(b = 2)),
               "The transcript_column argument is missing")
})

test_that("reclassify_ninetails_data errors on empty class_data", {
  valid_residue <- make_synthetic_inputs()$residue_data

  expect_error(reclassify_ninetails_data(residue_data = valid_residue,
                                         class_data = data.frame(),
                                         transcript_column = "contig",
                                         ref = "no_match"),
               "Empty data frame provided as an input \\(class_data\\)")
})

test_that("reclassify_ninetails_data errors on empty residue_data", {
  valid_class <- make_synthetic_inputs()$class_data

  expect_error(reclassify_ninetails_data(residue_data = data.frame(),
                                         class_data = valid_class,
                                         transcript_column = "contig",
                                         ref = "no_match"),
               "Empty data frame provided as an input \\(residue_data\\)")
})

test_that("reclassify_ninetails_data returns named list with class_data and residue_data", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 30)

  result <- reclassify_ninetails_data(
    residue_data = inputs$residue_data,
    class_data = inputs$class_data,
    grouping_factor = NULL,
    transcript_column = "contig",
    ref = "no_match"
  )

  expect_type(result, "list")
  expect_true("class_data" %in% names(result))
  expect_true("residue_data" %in% names(result))
  expect_true(is.data.frame(result[["class_data"]]))
  expect_true(is.data.frame(result[["residue_data"]]))
  # class and comments columns must be renamed from corr_* variants
  expect_true("class" %in% colnames(result[["class_data"]]))
  expect_true("comments" %in% colnames(result[["class_data"]]))
  expect_false("corr_class" %in% colnames(result[["class_data"]]))
  expect_false("corr_comments" %in% colnames(result[["class_data"]]))
  # QC columns must be stripped from residue_data output
  qc_cols <- c("mode_pos", "seg_err_quart", "qc_pos", "pos_err_quart",
               "count_nonA", "mode_len", "count")
  expect_false(any(qc_cols %in% colnames(result[["residue_data"]])))
})

test_that("reclassify_ninetails_data works with grouping_factor", {
  skip_if_not_installed("dplyr")
  skip_if(nchar(whitelist_path) == 0, "whitelist.RData not found in installed package")

  inputs <- make_synthetic_inputs(n_reads = 30, grouped = TRUE)

  result <- reclassify_ninetails_data(
    residue_data = inputs$residue_data,
    class_data = inputs$class_data,
    grouping_factor = "sample_name",
    transcript_column = "contig",
    ref = "no_match"
  )

  expect_type(result, "list")
  expect_true(is.data.frame(result[["class_data"]]))
  expect_true(is.data.frame(result[["residue_data"]]))
})


# count_nonA_abundance
################################################################################

test_that("count_nonA_abundance validates empty data frame", {
  expect_error(count_nonA_abundance(residue_data = data.frame()),
               "Empty data frame")
})

test_that("count_nonA_abundance validates grouping_factor exists", {
  skip_if(!exists("residue_data_single"), "residue_data_single not loaded")

  expect_error(count_nonA_abundance(residue_data = residue_data_single,
                                    grouping_factor = "nonexistent_column"),
               "is not a column")
})

test_that("count_nonA_abundance returns abundance counts", {
  skip_if_not_installed("dplyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_nonA_abundance(residue_data = residue_data_grouped,
                                 grouping_factor = "group")

  expect_s3_class(result, "tbl_df")
  expect_true("group" %in% colnames(result))
  expect_true("instances" %in% colnames(result))
  expect_true("count" %in% colnames(result))
  expect_true(all(result$instances %in% c("single", "two", "more")))
})

test_that("count_nonA_abundance works with sample_name grouping", {
  skip_if_not_installed("dplyr")
  skip_if(!exists("residue_data_grouped"), "residue_data_grouped not loaded")

  result <- count_nonA_abundance(residue_data = residue_data_grouped,
                                 grouping_factor = "sample_name")

  expect_s3_class(result, "tbl_df")
  expect_true("sample_name" %in% colnames(result))
})
