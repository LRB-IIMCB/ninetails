################################################################################
# Testing helper functions
################################################################################

# assert_condition
#################################################################################

test_that("assert_condition passes on TRUE", {
  expect_invisible(assert_condition(TRUE))
  expect_true(assert_condition(TRUE, "should not fire"))
})

test_that("assert_condition stops on FALSE", {
  expect_error(assert_condition(FALSE), "Assertion failed")
  expect_error(assert_condition(FALSE, "custom msg"), "custom msg")
})


# is_string
################################################################################

test_that("is_string returns TRUE for valid single character strings", {
  expect_true(is_string("hello world"))
  expect_true(is_string("a"))
  expect_true(is_string("/path/to/file.tsv"))
})

test_that("is_string returns FALSE for non-strings", {
  expect_false(is_string(123))
  expect_false(is_string(TRUE))
  expect_false(is_string(list("a")))
  expect_false(is_string(NA))
  expect_false(is_string(NA_character_))
})

test_that("is_string handles NULL via null.ok parameter", {
  expect_false(is_string(NULL, null.ok = FALSE))
  expect_true(is_string(NULL, null.ok = TRUE))
})

# assert_file_exists
################################################################################

test_that("assert_file_exists passes for existing file", {
  tmp <- tempfile()
  writeLines("test", tmp)
  on.exit(unlink(tmp))

  expect_invisible(assert_file_exists(tmp))
  expect_true(assert_file_exists(tmp))
})

test_that("assert_file_exists stops for nonexistent file", {
  fake_path <- "/nonexistent/path/file.txt"
  expect_error(assert_file_exists(fake_path), "File does not exist")
})

# assert_dir_exists
################################################################################

test_that("assert_dir_exists passes for existing directory", {
  tmp_dir <- tempdir()

  expect_invisible(assert_dir_exists(tmp_dir))
  expect_true(assert_dir_exists(tmp_dir))
})

test_that("assert_dir_exists stops for nonexistent directory", {
  fake_dir <- "/nonexistent/directory/path"
  expect_error(assert_dir_exists(fake_dir), "Directory does not exist")
})

# no_na
################################################################################

test_that("no_na returns TRUE for vectors without NAs", {
  expect_true(no_na(1:10))
  expect_true(no_na(c("a", "b", "c")))
  expect_true(no_na(c(TRUE, FALSE)))
  expect_true(no_na(numeric(0)))
})

test_that("no_na returns FALSE for vectors with NAs", {
  expect_false(no_na(c(1, NA, 3)))
  expect_false(no_na(c(NA)))
  expect_false(no_na(c("a", NA)))
  expect_false(no_na(c(TRUE, NA)))
})


# is_multifast5
################################################################################
test_that("is_multifast5 returns TRUE for multi-read structure", {
  mock_structure <- data.frame(
    name = c("read_abc123", "read_def456", "read_ghi789"),
    stringsAsFactors = FALSE
  )
  expect_true(is_multifast5(mock_structure))
})

test_that("is_multifast5 returns FALSE for single-read structure", {
  mock_structure <- data.frame(
    name = c("Raw", "UniqueGlobalKey", "Analyses"),
    stringsAsFactors = FALSE
  )
  expect_false(is_multifast5(mock_structure))
})


# winsorize_signal
################################################################################
test_that("winsorize_signal works correctly", {
  empty_vector <- vector()
  test_vector <- c(-1:11)
  test_vector_with_NAs <- c(-1:11, NA)
  test_vector_with_character <- c(-1:11,"A")
  expect_vector(winsorize_signal(test_vector), size=13)
  expect_equal(winsorize_signal(test_vector), c(0,0,1,2,3,4,5,6,7,8,9,10,10))
  expect_error(winsorize_signal(empty_vector))
  expect_error(winsorize_signal(test_vector_with_NAs))
  expect_error(winsorize_signal(test_vector_with_character))

})

test_that("winsorize_signal errors on non-atomic input", {
  expect_error(winsorize_signal(list(1, 2, 3)), "must be numeric|must be atomic")
})

test_that("winsorize_signal errors on missing argument", {
  expect_error(winsorize_signal(), "Signal vector is missing")
})


test_that("winsorize_signal preserves length", {
  test_vector <- c(100:300)
  result <- winsorize_signal(test_vector)
  expect_equal(length(result), length(test_vector))
})

# substitute_gaps
################################################################################

test_that("substitute_gaps fills short zero-gaps between identical nonzero values", {
  # gap of length 1 between two 1's → should be filled
  pseudomoves <- c(1, 1, 1, 0, 1, 1, 1)
  result <- substitute_gaps(pseudomoves)
  expect_equal(result, c(1, 1, 1, 1, 1, 1, 1))
})

test_that("substitute_gaps preserves longer zero-gaps", {
  # gap of length 2 between two 1's → should NOT be filled
  pseudomoves <- c(1, 1, 0, 0, 1, 1)
  result <- substitute_gaps(pseudomoves)
  expect_equal(result, pseudomoves)
})

test_that("substitute_gaps does not fill gaps between different nonzero values", {
  # gap between 1 and -1 → different flanking values, should not fill
  pseudomoves <- c(1, 1, 0, -1, -1)
  result <- substitute_gaps(pseudomoves)
  expect_equal(result, pseudomoves)
})


# get_mode
################################################################################

test_that("get_mode with method='value' returns most frequent value", {
  x <- c(rep(2, 5), rep(3, 4), rep(1, 4), rep(8, 2))
  result <- get_mode(x, method = "value")
  expect_equal(result, 2)
})

test_that("get_mode with method='value' returns multiple modes for multimodal data", {
  x <- c(rep(2, 5), rep(3, 5), rep(1, 1))
  result <- get_mode(x, method = "value")
  expect_true(all(c(2, 3) %in% result))
  expect_equal(length(result), 2)
})

test_that("get_mode errors on non-numeric input", {
  expect_error(get_mode(c("a", "b", "a")), "must be numeric")
})

# correct_labels
################################################################################

test_that("correct_labels replaces legacy label 'modified' with 'decorated'", {
  df <- data.frame(
    readname = c("read1", "read2"),
    class = c("modified", "unmodified"),
    stringsAsFactors = FALSE
  )
  result <- correct_labels(df)
  expect_true("decorated" %in% result$class)
})

test_that("correct_labels replaces 'unmodified' with 'blank'", {
  df <- data.frame(
    readname = c("read1"),
    class = c("unmodified"),
    stringsAsFactors = FALSE
  )
  result <- correct_labels(df)
  expect_equal(result$class, "blank")
})

test_that("correct_labels replaces 'blank' with 'unclassified'", {
  df <- data.frame(
    readname = c("read1"),
    class = c("blank"),
    stringsAsFactors = FALSE
  )
  result <- correct_labels(df)
  expect_equal(result$class, "unclassified")
})


# filter_dorado_summary
################################################################################

test_that("filter_dorado_summary filters reads correctly from data frame", {
  df <- data.frame(
    read_id = c("r1", "r2", "r3", "r4", "r5"),
    alignment_direction = c("+", "-", "*", "+", "+"),
    alignment_mapq = c(60, 30, 10, 0, 60),
    poly_tail_start = c(100, 200, 150, 300, 0),
    poly_tail_length = c(25, 15, 50, 20, 30),
    stringsAsFactors = FALSE
  )

  result <- filter_dorado_summary(df)

  # r1: passes all criteria
  # r2: passes all criteria
  # r3: fails - direction is "*"
  # r4: fails - mapq is 0
  # r5: fails - poly_tail_start is 0
  expect_equal(nrow(result), 2)
  expect_true(all(result$read_id %in% c("r1", "r2")))
})


test_that("filter_dorado_summary errors on non-data-frame non-string input", {
  expect_error(filter_dorado_summary(123), "must be a data.frame")
  expect_error(filter_dorado_summary(list(read_id = "r1")), "must be a data.frame")
})

test_that("filter_dorado_summary errors when required columns are missing", {
  # missing read_id
  df <- data.frame(poly_tail_length = c(20))
  expect_error(filter_dorado_summary(df), "Required columns missing.*read_id")

  # missing poly_tail_length
  df <- data.frame(read_id = c("r1"))
  expect_error(filter_dorado_summary(df), "Required columns missing.*poly_tail_length")
})

# count_trailing_chars
################################################################################

test_that("count_trailing_chars counts trailing characters correctly", {
  expect_equal(count_trailing_chars("ACGTTTTT", "T"), 5)
  expect_equal(count_trailing_chars("AAAAACGT", "T"), 1)
  expect_equal(count_trailing_chars("AAAAACGT", "A"), 0)
  expect_equal(count_trailing_chars("AAAA", "A"), 4)
})

test_that("count_trailing_chars errors on non-string input", {
  expect_error(count_trailing_chars(12345, "A"), "must be a character string")
  expect_error(count_trailing_chars(c("ABC", "DEF"), "A"), "must be a character string")
})


# reverse_complement
################################################################################

test_that("reverse_complement produces correct reverse complement", {
  expect_equal(reverse_complement("ATCG"), "CGAT")
  expect_equal(reverse_complement("AAAA"), "TTTT")
  expect_equal(reverse_complement("TTTT"), "AAAA")
  expect_equal(reverse_complement("AAATTT"), "AAATTT")
})

test_that("reverse_complement is its own inverse", {
  seq <- "ACGTACGT"
  expect_equal(reverse_complement(reverse_complement(seq)), seq)
})


test_that("reverse_complement handles lowercase input", {
  # docstring says toupper is applied
  expect_equal(reverse_complement("atcg"), "CGAT")
  expect_equal(reverse_complement("aAnN"), "NNTT")
})

test_that("reverse_complement handles single nucleotide", {
  expect_equal(reverse_complement("A"), "T")
  expect_equal(reverse_complement("G"), "C")
})

test_that("reverse_complement errors on non-string input", {
  expect_error(reverse_complement(123), "must be a character string")
  expect_error(reverse_complement(c("ATCG", "GCTA")), "must be a character string")
})


# edit_distance_hw
################################################################################

test_that("edit_distance_hw returns 0 for exact match inside target", {
  expect_equal(edit_distance_hw("ATCG", "XXATCGXX"), 0)
  expect_equal(edit_distance_hw("ATCG", "ATCG"), 0)
})

test_that("edit_distance_hw returns correct distance for mismatches", {
  # single substitution
  expect_equal(edit_distance_hw("ATCG", "ATCC"), 1)
})

test_that("edit_distance_hw finds best alignment in sliding window", {
  # query is embedded in target with flanking noise
  expect_equal(edit_distance_hw("AAAA", "TTTTAAAATT"), 0)
  # one mismatch in best window
  expect_equal(edit_distance_hw("AAAA", "TTTTAAAGTT"), 1)
})

test_that("edit_distance_hw handles query longer than target", {
  # query > target: no sliding window positions, falls back to full alignment
  result <- edit_distance_hw("ATCGATCG", "ATCG")
  expect_true(is.numeric(result))
  expect_gte(result, 0)
})

test_that("edit_distance_hw handles identical query and target", {
  expect_equal(edit_distance_hw("ATCG", "ATCG"), 0)
})

test_that("edit_distance_hw handles completely different strings", {
  result <- edit_distance_hw("AAAA", "TTTT")
  expect_equal(result, 4)
})



# check_fast5_filetype — argument validation only
################################################################################

test_that("check_fast5_filetype errors on missing arguments", {
  expect_error(check_fast5_filetype(basecall_group = "Basecall_1D_000"),
               "workspace")
  expect_error(check_fast5_filetype(workspace = "/some/path"),
               "Basecall group is missing")
})


# is_RNA — argument validation only
################################################################################

test_that("is_RNA errors on nonexistent file path", {
  expect_error(is_RNA(fast5_file = "/nonexistent/file.fast5",
                      read_id = "read_abc"),
               "does not exist")
})


# check_output_directory
################################################################################
test_that("check_output_directory creates new directory", {
  tmp_dir <- file.path(tempdir(), paste0("ninetails_test_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE))

  # ensure it doesn't exist yet
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)

  mock_log <- function(msg, type = "INFO", section = NULL) invisible(NULL)
  result <- check_output_directory(tmp_dir, mock_log)

  expect_true(result)
  expect_true(dir.exists(tmp_dir))
})







