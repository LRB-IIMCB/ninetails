# testing plot_squiggle

test_that("plot_squiggle correctly parses data & draws a signal plot", {
  empty_tempfile = tempfile()
  nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt', package = 'ninetails')
  test_workspace = system.file('extdata', 'test_data', 'basecalled_fast5', package = 'ninetails')
  uncalled_workspace= system.file('extdata', 'test_data', 'uncalled_fast5', package = 'ninetails')
  test_basecall_group = 'Basecall_1D_000'
  wrong_basecall_group = 'Basecall_1D_003'

  # wrong readname
  expect_error(plot_squiggle(readname = "this_is_a_wrong_readname",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong workspace
  expect_error(plot_squiggle(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong nanopolish
  expect_error(plot_squiggle(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = empty_tempfile,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong seq_sum
  expect_error(plot_squiggle(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = nanopolish,
                             sequencing_summary = empty_tempfile,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong basecall group
  expect_error(plot_squiggle(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = wrong_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))

})


test_that("plot_tail_range correctly parses data & draws a signal plot", {
  empty_tempfile = tempfile()
  nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt', package = 'ninetails')
  test_workspace = system.file('extdata', 'test_data', 'basecalled_fast5', package = 'ninetails')
  uncalled_workspace= system.file('extdata', 'test_data', 'uncalled_fast5', package = 'ninetails')
  test_basecall_group = 'Basecall_1D_000'
  wrong_basecall_group = 'Basecall_1D_003'

  # wrong readname
  expect_error(plot_tail_range(readname = "this_is_a_wrong_readname",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong workspace
  expect_error(plot_tail_range(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong nanopolish
  expect_error(plot_tail_range(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = empty_tempfile,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong seq_sum
  expect_error(plot_tail_range(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = nanopolish,
                             sequencing_summary = empty_tempfile,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong basecall group
  expect_error(plot_tail_range(readname = "9c11d71e-eaaa-413f-958e-4ca1254e0369",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = wrong_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))


})










