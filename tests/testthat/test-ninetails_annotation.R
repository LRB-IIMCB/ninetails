################################################################################
# Testing ninetails_annotation
################################################################################

# annotate_with_biomart
################################################################################

test_that("annotate_with_biomart validates missing input_data", {
  expect_error(annotate_with_biomart(organism = "mmusculus"),"provide a dataframe")
  expect_error(annotate_with_biomart(input_data = data.frame(),
                                     organism = "mmusculus"), "Empty data frame")
})


test_that("annotate_with_biomart rejects empty attributes_to_get", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")
  skip_if_not_installed("biomaRt")

  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  expect_error(annotate_with_biomart(input_data = test_data,
                                     attributes_to_get = character(0),
                                     organism = "mmusculus"),
               "please provide attributes to get from biomart")
})


test_that("annotate_with_biomart requires either organism or mart_to_use", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  expect_error(annotate_with_biomart(input_data = test_data,
                                     organism = NULL,
                                     mart_to_use = NULL),
               "Either 'organism' or 'mart_to_use' must be provided")
})


test_that("annotate_with_biomart rejects both organism and mart_to_use", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")
  skip_if_not_installed("biomaRt")

  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  mock_mart <- structure(list(), class = "Mart")

  expect_error(annotate_with_biomart(input_data = test_data,
                                     organism = "mmusculus",
                                     mart_to_use = mock_mart),
               "Only one of 'organism' or 'mart_to_use'")
})


test_that("annotate_with_biomart validates mart_to_use is Mart class", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  expect_error(annotate_with_biomart(input_data = test_data,
                                     mart_to_use = list(not = "a mart")),"must be a mart")
})


test_that("annotate_with_biomart rejects invalid organism name", {
  # The else stop fires before any biomaRt network call so this runs offline.
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")
  skip_if_not_installed("biomaRt")

  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  expect_error(annotate_with_biomart(input_data = test_data,
                                     organism = "invalid_organism"),
               "provide valid organism name")
})


test_that("annotate_with_biomart errors when ensembl_transcript_id_short column is absent", {
  # No explicit column guard exists in the source — the function currently
  # errors with a cryptic dplyr rename message when the column is missing.
  # This test documents the gap and ensures the function at least errors.
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")
  skip_if_not_installed("biomaRt")

  test_data <- merged_nonA_char[, setdiff(colnames(merged_nonA_char),
                                          "ensembl_transcript_id_short")]

  expect_error(annotate_with_biomart(input_data = test_data,organism = "mmusculus"))
})


test_that("annotate_with_biomart biomaRt version is at least 2.40", {
  skip_if_not_installed("biomaRt")

  biomart_version <- utils::packageVersion("biomaRt")
  expect_s3_class(biomart_version, "package_version")
  expect_true(biomart_version >= "2.40.0")
})


# Note: full integration tests for organism branches (mmusculus, hsapiens,
# athaliana, scerevisiae) and the getBM + join return path require a live
# Ensembl connection and are not included here.
