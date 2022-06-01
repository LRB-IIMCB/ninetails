# testing plot_squiggle

test_that("plot_squiggle correctly parses data & draws a signal plot", {
  empty_tempfile = tempfile()

  # wrong readname
  expect_error(plot_squiggle(readname = "this_is_a_wrong_readname",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong workspace
  expect_error(plot_squiggle(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong nanopolish
  expect_error(plot_squiggle(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = empty_tempfile,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong seq_sum
  expect_error(plot_squiggle(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = nanopolish,
                             sequencing_summary = empty_tempfile,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong basecall group
  expect_error(plot_squiggle(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = wrong_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))


})




test_that("plot_squiggle correctly parses data & draws a signal plot", {
  empty_tempfile = tempfile()

  # wrong readname
  expect_error(plot_tail_range(readname = "this_is_a_wrong_readname",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = test_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong workspace
  expect_error(plot_tail_range(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong nanopolish
  expect_error(plot_tail_range(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = empty_tempfile,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong seq_sum
  expect_error(plot_tail_range(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = nanopolish,
                             sequencing_summary = empty_tempfile,
                             workspace = uncalled_workspace,
                             basecall_group = test_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))
  # wrong basecall group
  expect_error(plot_tail_range(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                             nanopolish = nanopolish,
                             sequencing_summary = sequencing_summary,
                             workspace = uncalled_workspace,
                             basecall_group = wrong_basecall_group,
                             moves=FALSE,
                             rescale=FALSE))


})










