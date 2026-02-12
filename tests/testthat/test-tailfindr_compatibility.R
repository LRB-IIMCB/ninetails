################################################################################
# Testing tailfindr compatibility
################################################################################

# convert_tailfindr_output
################################################################################

test_that("convert_tailfindr_output converts data frame with correct column mapping", {
  input_df <- data.frame(
    read_id = c("read_001", "read_002", "read_003"),
    tail_start = c(100L, 200L, 300L),
    tail_end = c(130L, 225L, 340L),
    tail_length = c(30, 25, 40),
    stringsAsFactors = FALSE
  )

  result <- convert_tailfindr_output(input_df)

  # renamed columns present
  expect_true("readname" %in% colnames(result))
  expect_true("polya_start" %in% colnames(result))
  expect_true("polya_length" %in% colnames(result))

  # original column names should not persist (read_id, tail_start, tail_length)
  expect_false("read_id" %in% colnames(result))
  expect_false("tail_start" %in% colnames(result))
  expect_false("tail_length" %in% colnames(result))
})

test_that("convert_tailfindr_output preserves read identifiers correctly", {
  input_df <- data.frame(
    read_id = c("abc123", "def456"),
    tail_start = c(10L, 20L),
    tail_end = c(30L, 50L),
    tail_length = c(20, 30),
    stringsAsFactors = FALSE
  )

  result <- convert_tailfindr_output(input_df)
  expect_equal(result$readname, c("abc123", "def456"))
})



test_that("convert_tailfindr_output assigns qc_tag PASS for tails >= 10 nt", {
  input_df <- data.frame(
    read_id = c("r1", "r2", "r3"),
    tail_start = c(10L, 20L, 30L),
    tail_end = c(20L, 30L, 40L),
    tail_length = c(10, 15, 50),
    stringsAsFactors = FALSE
  )

  result <- convert_tailfindr_output(input_df)
  expect_true(all(result$qc_tag == "PASS"))
})


test_that("convert_tailfindr_output output contains all expected columns", {
  input_df <- data.frame(
    read_id = "r1",
    tail_start = 100L,
    tail_end = 130L,
    tail_length = 30,
    stringsAsFactors = FALSE
  )

  result <- convert_tailfindr_output(input_df)
  expected_cols <- c("readname", "polya_start", "tail_end",
                     "polya_length", "transcript_start", "contig", "qc_tag")
  expect_true(all(expected_cols %in% colnames(result)))
})



# check_polya_length_filetype
################################################################################

test_that("check_polya_length_filetype detects nanopolish format via qc_tag column", {
  nanopolish_df <- data.frame(
    readname = c("r1", "r2"),
    contig = c("tx1", "tx2"),
    polya_start = c(100L, 200L),
    polya_length = c(30, 40),
    qc_tag = c("PASS", "SUFFCLIP"),
    stringsAsFactors = FALSE
  )

  result <- check_polya_length_filetype(nanopolish_df)

  expect_type(result, "list")
  expect_equal(result$file_type, "nanopolish")
  expect_true(is.data.frame(result$data))
  expect_equal(nrow(result$data), 2)
})

test_that("check_polya_length_filetype returns list with 'data' and 'file_type' elements", {
  nanopolish_df <- data.frame(
    readname = "r1",
    qc_tag = "PASS",
    polya_length = 25,
    stringsAsFactors = FALSE
  )

  result <- check_polya_length_filetype(nanopolish_df)
  expect_named(result, c("data", "file_type"))
})



test_that("check_polya_length_filetype detects tailfindr DRS format and converts", {
  tailfindr_drs_df <- data.frame(
    read_id = c("r1", "r2"),
    tail_start = c(100L, 200L),
    tail_end = c(130L, 240L),
    tail_length = c(30, 40),
    samples_per_nt = c(10.5, 11.2),
    stringsAsFactors = FALSE
  )

  result <- check_polya_length_filetype(tailfindr_drs_df)

  expect_equal(result$file_type, "tailfindr_drs")
  expect_true(is.data.frame(result$data))
  # converted output should have nanopolish-like columns
  expect_true("readname" %in% colnames(result$data))
  expect_true("polya_length" %in% colnames(result$data))
  expect_true("qc_tag" %in% colnames(result$data))
})



test_that("check_polya_length_filetype rejects tailfindr cDNA format", {
  tailfindr_cdna_df <- data.frame(
    read_id = c("r1"),
    tail_is_valid = c(TRUE),
    tail_length = c(30),
    stringsAsFactors = FALSE
  )

  expect_error(
    check_polya_length_filetype(tailfindr_cdna_df),
    "cDNA output is not compatible"
  )
})

test_that("check_polya_length_filetype errors on unrecognized column layout", {
  unknown_df <- data.frame(
    some_column = c("a", "b"),
    another_col = c(1, 2),
    stringsAsFactors = FALSE
  )

  expect_error(
    check_polya_length_filetype(unknown_df),
    "Unrecognized input format"
  )
})

