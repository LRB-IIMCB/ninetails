################################################################################
# Testing plotting functions
################################################################################

# HELPER FUNCTIONS
################################################################################
# Suppress ggplot2/cowplot warnings during testing
suppress_plot_warnings <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

# Helper: Create minimal processing_info for plot_nanopolish_qc tests
# (not included in test_data.RData)
create_test_processing_info <- function(grouped = FALSE) {
  qc_tags <- c("PASS", "ADAPTER", "SUFFCLIP", "NOREGION", "READ_FAILED_LOAD")

  if (!grouped) {
    data.frame(
      qc_tag = qc_tags,
      n = c(500, 50, 30, 15, 5),
      stringsAsFactors = FALSE
    )
  } else {
    rbind(
      data.frame(
        sample_name = "WT_1",
        qc_tag = qc_tags,
        n = c(250, 25, 15, 8, 2),
        stringsAsFactors = FALSE
      ),
      data.frame(
        sample_name = "KO_1",
        qc_tag = qc_tags,
        n = c(300, 30, 20, 10, 4),
        stringsAsFactors = FALSE
      )
    )
  }
}


# plot_squiggle_fast5
################################################################################

test_that("plot_squiggle_fast5 correctly parses data & draws a signal plot", {
  empty_tempfile = tempfile()
  nanopolish = system.file('extdata', 'test_data', 'legacy', 'nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary = system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  test_workspace = system.file('extdata', 'test_data', 'legacy','basecalled_fast5', package = 'ninetails')
  uncalled_workspace= system.file('extdata', 'test_data', 'legacy', 'uncalled_fast5', package = 'ninetails')
  test_basecall_group = 'Basecall_1D_000'
  wrong_basecall_group = 'Basecall_1D_003'

  # wrong readname
  expect_error(plot_squiggle_fast5(readname = "this_is_a_wrong_readname",
                                   nanopolish = nanopolish,
                                   sequencing_summary = sequencing_summary,
                                   workspace = test_workspace,
                                   basecall_group = test_basecall_group,
                                   moves=FALSE,
                                   rescale=FALSE))
  # wrong workspace
  expect_error(plot_squiggle_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                   nanopolish = nanopolish,
                                   sequencing_summary = sequencing_summary,
                                   workspace = uncalled_workspace,
                                   basecall_group = test_basecall_group,
                                   moves=FALSE,
                                   rescale=FALSE))
  # wrong nanopolish
  expect_error(plot_squiggle_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                   nanopolish = empty_tempfile,
                                   sequencing_summary = sequencing_summary,
                                   workspace = test_workspace,
                                   basecall_group = test_basecall_group,
                                   moves=FALSE,
                                   rescale=FALSE))
  # wrong seq_sum
  expect_error(plot_squiggle_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                   nanopolish = nanopolish,
                                   sequencing_summary = empty_tempfile,
                                   workspace = test_workspace,
                                   basecall_group = test_basecall_group,
                                   moves=FALSE,
                                   rescale=FALSE))
  # wrong basecall group
  expect_error(plot_squiggle_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                   nanopolish = nanopolish,
                                   sequencing_summary = sequencing_summary,
                                   workspace = test_workspace,
                                   basecall_group = wrong_basecall_group,
                                   moves=FALSE,
                                   rescale=FALSE))

})



# plot_tail_range_fast5
################################################################################

test_that("plot_tail_range correctly parses data & draws a signal plot", {
  empty_tempfile = tempfile()
  nanopolish = system.file('extdata', 'test_data', 'legacy','nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary = system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  test_workspace = system.file('extdata', 'test_data', 'legacy','basecalled_fast5', package = 'ninetails')
  uncalled_workspace= system.file('extdata', 'test_data', 'legacy', 'uncalled_fast5', package = 'ninetails')
  test_basecall_group = 'Basecall_1D_000'
  wrong_basecall_group = 'Basecall_1D_003'

  # wrong readname
  expect_error(plot_tail_range_fast5(readname = "this_is_a_wrong_readname",
                                     nanopolish = nanopolish,
                                     sequencing_summary = sequencing_summary,
                                     workspace = test_workspace,
                                     basecall_group = test_basecall_group,
                                     moves=FALSE,
                                     rescale=FALSE))
  # wrong workspace
  expect_error(plot_tail_range_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                     nanopolish = nanopolish,
                                     sequencing_summary = sequencing_summary,
                                     workspace = uncalled_workspace,
                                     basecall_group = test_basecall_group,
                                     moves=FALSE,
                                     rescale=FALSE))
  # wrong nanopolish
  expect_error(plot_tail_range_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                     nanopolish = empty_tempfile,
                                     sequencing_summary = sequencing_summary,
                                     workspace = test_workspace,
                                     basecall_group = test_basecall_group,
                                     moves=FALSE,
                                     rescale=FALSE))
  # wrong seq_sum
  expect_error(plot_tail_range_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                     nanopolish = nanopolish,
                                     sequencing_summary = empty_tempfile,
                                     workspace = test_workspace,
                                     basecall_group = test_basecall_group,
                                     moves=FALSE,
                                     rescale=FALSE))
  # wrong basecall group
  expect_error(plot_tail_range_fast5(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                                     nanopolish = nanopolish,
                                     sequencing_summary = sequencing_summary,
                                     workspace = test_workspace,
                                     basecall_group = wrong_basecall_group,
                                     moves=FALSE,
                                     rescale=FALSE))


})


# plot_tail_chunk
################################################################################

test_that("plot_tail_chunk validates required arguments", {
  expect_error(
    ninetails::plot_tail_chunk(chunk_name = "test_chunk"),
    "List of tail chunks is missing"
  )

  expect_error(
    ninetails::plot_tail_chunk(tail_chunk_list = list()),
    "Chunk_name is missing"
  )

  expect_error(
    ninetails::plot_tail_chunk(
      chunk_name = "test",
      tail_chunk_list = "not_a_list"
    ),
    "not a list"
  )
})


test_that("plot_tail_chunk returns ggplot object for valid input", {
  skip_if_not_installed("ggplot2")
  skip_if(!exists("test_tail_chunk_list"), "test_tail_chunk_list not loaded")

  # Get first valid chunk name from the test data
  chunk_names <- unlist(lapply(test_tail_chunk_list, names))
  skip_if(length(chunk_names) == 0, "No chunks in test_tail_chunk_list")

  first_chunk <- chunk_names[1]

  result <- suppress_plot_warnings(
    ninetails::plot_tail_chunk(
      chunk_name = first_chunk,
      tail_chunk_list = test_tail_chunk_list
    )
  )

  expect_s3_class(result, "gg")
  expect_s3_class(result, "ggplot")
})


# plot_gaf
################################################################################



# plot_multiple_gaf
################################################################################



# plot_class_counts
################################################################################




# plot_residue_counts
################################################################################



# plot_nanopolish_qc
################################################################################




# plot_tail_distribution
################################################################################




# plot_panel_characteristics
################################################################################




# plot_rug_density
################################################################################



# plot_nonA_abundance
################################################################################







