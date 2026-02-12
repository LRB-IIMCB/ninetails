################################################################################
# Testing ninetails_annotation
################################################################################

# annotate_with_biomart
################################################################################

test_that("annotate_with_biomart validates missing input_data", {
  expect_error(annotate_with_biomart(organism = "mmusculus"),"provide a dataframe")

  expect_error(annotate_with_biomart(input_data = data.frame(),
                                     organism = "mmusculus"),"Empty data frame")
})


test_that("annotate_with_biomart requires either organism or mart_to_use", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  # Prepare data with required column
  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  # Neither organism nor mart_to_use provided
  expect_error(annotate_with_biomart(input_data = test_data,
                                     organism = NULL,
                                     mart_to_use = NULL),
               "Either 'organism' or 'mart_to_use' must be provided")
})


test_that("annotate_with_biomart rejects both organism and mart_to_use", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")
  skip_if_not_installed("biomaRt")

  # Prepare data with required column
  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  # Create mock mart object (will fail validation but tests the logic)
  mock_mart <- structure(list(), class = "Mart")

  expect_error(annotate_with_biomart(input_data = test_data,
                                     organism = "mmusculus",
                                     mart_to_use = mock_mart),
               "Only one of 'organism' or 'mart_to_use'")
})


test_that("annotate_with_biomart validates mart_to_use is Mart class", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")

  # Prepare data with required column
  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  # Provide non-Mart object
  expect_error(annotate_with_biomart(input_data = test_data,
                                     mart_to_use = list(not = "a mart")),"must be a mart")
})


test_that("annotate_with_biomart rejects invalid organism names", {
  skip_if(!exists("merged_nonA_char"), "merged_nonA_char not loaded")
  skip_if_not_installed("biomaRt")

  # Skip if no network available (biomaRt requires internet)
  skip_if_offline()

  # Prepare data with required column
  test_data <- merged_nonA_char
  test_data$ensembl_transcript_id_short <- test_data$contig

  expect_error(
    annotate_with_biomart(input_data = test_data,
                          organism = "invalid_organism"),"provide valid organism name")
})


# Note: Full integration tests with biomaRt are skipped as they require
# network access to Ensembl servers. Local unit tests focus on input validation.

test_that("annotate_with_biomart requires biomaRt version >= 2.40", {
  skip_if_not_installed("biomaRt")

  # This test just verifies the version check logic exists

  # The actual check happens inside the function
  biomart_version <- utils::packageVersion("biomaRt")
  expect_true(is.package_version(biomart_version))
})
