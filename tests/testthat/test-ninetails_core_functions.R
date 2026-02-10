# Testing core functions

# testing extract_polya_data

test_that("extract_polya_data correctly reads from provided files", {
  empty_tempfile = tempfile()
  nanopolish = system.file('extdata', 'test_data', 'legacy', 'nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary = system.file('extdata', 'test_data', 'legacy', 'sequencing_summary.txt', package = 'ninetails')

  # empty nanopolish
  expect_error(extract_polya_data(nanopolish = empty_tempfile, sequencing_summary = sequencing_summary, pass_only = TRUE))
  # empty sequencing summary
  expect_error(extract_polya_data(nanopolish = nanopolish, sequencing_summary = empty_tempfile, pass_only = TRUE))
  # all reads parsed correctly
  expect_length(rownames(extract_polya_data(nanopolish = nanopolish,
                                            sequencing_summary = sequencing_summary,
                                            pass_only = TRUE)), 17)
  # all reads parsed correctly
  expect_length(rownames(extract_polya_data(nanopolish = nanopolish,
                                            sequencing_summary = sequencing_summary,
                                            pass_only = FALSE)), 20)
})




test_that("extract_tail_data correctly reads from provided files", {

  empty_tempfile = tempfile()
  test_workspace = system.file('extdata', 'test_data', 'legacy','basecalled_fast5', package = 'ninetails')
  uncalled_workspace= system.file('extdata', 'test_data','legacy', 'uncalled_fast5', package = 'ninetails')
  test_basecall_group = 'Basecall_1D_000'
  wrong_basecall_group = 'Basecall_1D_003'


  # no required readname
  expect_error(extract_tail_data(readname = "",
                                 polya_summary = test_polya_summary,
                                 workspace = test_workspace,
                                 basecall_group = test_basecall_group))
  # wrong fast5 files (uncalled)
  expect_error(extract_tail_data(readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
                                 polya_summary = test_polya_summary,
                                 workspace = uncalled_workspace,
                                 basecall_group = test_basecall_group))
  # nonexistent basecall group
  expect_error(extract_tail_data(readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
                                 polya_summary = test_polya_summary,
                                 workspace = test_workspace,
                                 basecall_group = wrong_basecall_group))
  # empty polya_summary file provided
  expect_error(extract_tail_data(readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
                                 polya_summary = empty_tempfile,
                                 workspace = test_workspace,
                                 basecall_group = test_basecall_group))

  #correct input provided, test output length
  expect_length(extract_tail_data(readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
                                  polya_summary = test_polya_summary,
                                  workspace = test_workspace,
                                  basecall_group = test_basecall_group), 4)
  #check if returned object is a list
  expect_type(extract_tail_data(readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
                                polya_summary = test_polya_summary,
                                workspace = test_workspace,
                                basecall_group = test_basecall_group), "list")

})


test_that("create_tail_feature_list correctly produces the output", {

  empty_tempfile = tempfile()
  num_cores=2
  nanopolish = system.file('extdata', 'test_data', 'legacy','nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary = system.file('extdata', 'test_data','legacy', 'sequencing_summary.txt', package = 'ninetails')
  test_workspace = system.file('extdata', 'test_data', 'legacy','basecalled_fast5', package = 'ninetails')
  uncalled_workspace= system.file('extdata', 'test_data', 'legacy','uncalled_fast5', package = 'ninetails')
  test_basecall_group = 'Basecall_1D_000'
  wrong_basecall_group = 'Basecall_1D_003'


  # wrong workspace (uncalled)
  expect_error(invisible(capture.output(create_tail_feature_list(nanopolish = nanopolish,
                                                                 sequencing_summary = sequencing_summary,
                                                                 workspace = uncalled_workspace,
                                                                 num_cores=num_cores,
                                                                 basecall_group = test_basecall_group,
                                                                 pass_only=TRUE))))
  # wrong basecall group
  expect_error(invisible(capture.output(create_tail_feature_list(nanopolish = nanopolish,
                                                                 sequencing_summary = sequencing_summary,
                                                                 workspace = test_workspace,
                                                                 num_cores=num_cores,
                                                                 basecall_group = wrong_basecall_group,
                                                                 pass_only=TRUE))))
  # nonexistent file
  expect_error(invisible(capture.output(create_tail_feature_list(nanopolish = empty_tempfile,
                                                                 sequencing_summary = sequencing_summary,
                                                                 workspace = test_workspace,
                                                                 num_cores=num_cores,
                                                                 basecall_group = wrong_basecall_group,
                                                                 pass_only=TRUE))))
  #check if it returns a nonempty file
  expect_type(invisible(capture.output(create_tail_feature_list(nanopolish = nanopolish,
                                                                sequencing_summary = sequencing_summary,
                                                                workspace = test_workspace,
                                                                num_cores=num_cores,
                                                                basecall_group = test_basecall_group,
                                                                pass_only=TRUE))), "character")

})




test_that("create_gaf correctly produces the output", {

  test_chunk <- test_tail_chunk_list[[1]][[1]][[1]]
  test_chunk_empty <- vector()

  test_chunk_noninteger <- c(708, 675, 672, 666, 667, 712, 692, 684, 670, 678, 647, 682, 665, 673, 691, 695, 687, 694, 655, 680, 672, 676, 656, 663, 686,
                             669, 702, 677, 666, 691, 682, 694, "ABCD", 651, 674, 658, 650, 678, 673, 672, 664, 663, 648, 668, 672, 647, 644, 654, 670, 687,
                             695, 687, 687, 687, 687, 687, 691, 691, 687, 695, 691, 687, 691, 695, 691, 695, 687, 691, 687, 687, 687, 691, 695, 695, 695,
                             695, 687, 687, 691, 695, 695, 687, 691, 687, 695, 695, 687, 691, 691, 691, 687, 687, 691, 695, 687, 687, 695, 691, 695, 691)

  test_chunk_short <- c(100:110)


  #check if returns numeric matrix
  expect_type(create_gaf(test_chunk, method="s"), "double")

  # noninteger vec as an input
  expect_error(create_gaf(test_chunk_noninteger, method="s"))

  # too short vec as an input
  expect_error(create_gaf(test_chunk_short, method="s"))

  # empty vec as an input
  expect_error(create_gaf(test_chunk_empty, method="s"))


})



test_that("filter_signal_by_threshold correctly produces the output", {

  test_signal <- test_tail_feature_list[["tail_feature_list"]][[1]][[2]]
  test_signal_empty <- vector()

  test_signal_noninteger <- c(708, 675, 672, 666, 667, 712, 692, 684, 670, 678, 647, 682, 665, 673, 691, 695, 687, 694, 655, 680, 672, 676, 656, 663, 686,
                             669, 702, 677, 666, 691, 682, 694, "ABCD", 651, 674, 658, 650, 678, 673, 672, 664, 663, 648, 668, 672, 647, 644, 654, 670, 687,
                             695, 687, 687, 687, 687, 687, 691, 691, 687, 695, 691, 687, 691, 695, 691, 695, 687, 691, 687, 687, 687, 691, 695, 695, 695,
                             695, 687, 687, 691, 695, 695, 687, 691, 687, 695, 695, 687, 691, 691, 691, 687, 687, 691, 695, 687, 687, 695, 691, 695, 691)



  #check if returns numeric vector
  expect_type(filter_signal_by_threshold(test_signal), "double")

  # noninteger vec as an input
  expect_error(filter_signal_by_threshold(test_signal_noninteger))

  # empty vec as an input
  expect_error(filter_signal_by_threshold(test_signal_empty))


})


test_that("create_tail_chunk_list correctly produces the output", {

  #num cores is character
  expect_error(create_tail_chunk_list(tail_feature_list =test_tail_feature_list, num_cores="Putin is an asshole"))

})


test_that("split_tail_centered correctly produces the output", {

  #readname is absent in tail features list
  expect_error(split_tail_centered(readname="Wladimir Wladimirowicz", tail_feature_list=test_tail_feature_list))

  #check if returned object is a list
  expect_type(split_tail_centered(readname = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b",
                                  tail_feature_list = test_tail_feature_list), "list")

})



