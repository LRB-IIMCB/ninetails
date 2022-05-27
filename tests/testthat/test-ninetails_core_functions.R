# Testing core functions

# testing extract_polya_data

test_that("extract_polya_data correctly reads from provided files", {
  empty_tempfile = tempfile()
  test_nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', package = 'ninetails')
  test_sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt',
                                        package = 'ninetails')
  # empty nanopolish
  expect_error(extract_polya_data(nanopolish = empty_tempfile, sequencing_summary = test_sequencing_summary, pass_only = TRUE))
  # empty sequencing summary
  expect_error(extract_polya_data(nanopolish = test_nanopolish, sequencing_summary = empty_tempfile, pass_only = TRUE))
  # all reads parsed correctly
  expect_length(rownames(extract_polya_data(nanopolish = test_nanopolish,
                                            sequencing_summary = test_sequencing_summary,
                                            pass_only = TRUE)), 20)
  # all reads parsed correctly
  expect_length(rownames(extract_polya_data(nanopolish = test_nanopolish,
                                            sequencing_summary = test_sequencing_summary,
                                            pass_only = FALSE)), 22)
})



test_that("extract_tail_data correctly reads from provided files", {

  empty_tempfile = tempfile()
  test_polya_summary <- extract_polya_data(nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv',
                                                                    package = 'ninetails'),
                                           sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt',
                                                                            package = 'ninetails'),
                                           pass_only = TRUE)

  test_workspace = system.file('extdata', 'test_data', 'basecalled_fast5',
                               package = 'ninetails')
  uncalled_workspace = system.file('extdata', 'test_data', 'uncalled_fast5',
                                   package = 'ninetails')
  test_basecall_group = 'Basecall_1D_000'
  wrong_basecall_group = 'Basecall_1D_003'

  # no required readname
  expect_error(extract_polya_data(readname = "",
                                  polya_summary = test_polya_summary,
                                  workspace = test_workspace,
                                  basecall_group = test_basecall_group))
  # wrong fast5 files (uncalled)
  expect_error(extract_polya_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                  polya_summary = test_polya_summary,
                                  workspace = uncalled_workspace,
                                  basecall_group = test_basecall_group))
  # nonexistent basecall group
  expect_error(extract_polya_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                  polya_summary = test_polya_summary,
                                  workspace = test_workspace,
                                  basecall_group = wrong_basecall_group))
  # empty polya_summary file provided
  expect_error(extract_polya_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                  polya_summary = empty_tempfile,
                                  workspace = test_workspace,
                                  basecall_group = test_basecall_group))

})




