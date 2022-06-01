# Testing core functions

# testing extract_polya_data

test_that("extract_polya_data correctly reads from provided files", {
  empty_tempfile = tempfile()

  # empty nanopolish
  expect_error(extract_polya_data(nanopolish = empty_tempfile, sequencing_summary = sequencing_summary, pass_only = TRUE))
  # empty sequencing summary
  expect_error(extract_polya_data(nanopolish = nanopolish, sequencing_summary = empty_tempfile, pass_only = TRUE))
  # all reads parsed correctly
  expect_length(rownames(extract_polya_data(nanopolish = nanopolish,
                                            sequencing_summary = sequencing_summary,
                                            pass_only = TRUE)), 20)
  # all reads parsed correctly
  expect_length(rownames(extract_polya_data(nanopolish = nanopolish,
                                            sequencing_summary = sequencing_summary,
                                            pass_only = FALSE)), 22)
})



test_that("extract_tail_data correctly reads from provided files", {

  empty_tempfile = tempfile()


  # no required readname
  expect_error(extract_tail_data(readname = "",
                                  polya_summary = test_polya_summary,
                                  workspace = test_workspace,
                                  basecall_group = test_basecall_group))
  # wrong fast5 files (uncalled)
  expect_error(extract_tail_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                  polya_summary = test_polya_summary,
                                  workspace = uncalled_workspace,
                                  basecall_group = test_basecall_group))
  # nonexistent basecall group
  expect_error(extract_tail_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                  polya_summary = test_polya_summary,
                                  workspace = test_workspace,
                                  basecall_group = wrong_basecall_group))
  # empty polya_summary file provided
  expect_error(extract_tail_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                  polya_summary = empty_tempfile,
                                  workspace = test_workspace,
                                  basecall_group = test_basecall_group))

  #correct input provided, test output length
  expect_length(extract_tail_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                 polya_summary = test_polya_summary,
                                 workspace = test_workspace,
                                 basecall_group = test_basecall_group), 6)
  #check if returned object is a list
  expect_type(extract_tail_data(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                  polya_summary = test_polya_summary,
                                  workspace = test_workspace,
                                  basecall_group = test_basecall_group), "list")

})



test_that("create_tail_feature_list correctly produces the output", {

  empty_tempfile = tempfile()


  # wrong workspace (uncalled)
  expect_error(invisible(capture.output(create_tail_feature_list(nanopolish = nanopolish,
                                                                 sequencing_summary = sequencing_summary,
                                                                 workspace = uncalled_workspace,
                                                                 num_cores=2,
                                                                 basecall_group = test_basecall_group,
                                                                 pass_only=TRUE))))
  # wrong basecall group
  expect_error(invisible(capture.output(create_tail_feature_list(nanopolish = nanopolish,
                                                                 sequencing_summary = sequencing_summary,
                                                                 workspace = test_workspace,
                                                                 num_cores=2,
                                                                 basecall_group = wrong_basecall_group,
                                                                 pass_only=TRUE))))
  # nonexistent file
  expect_error(invisible(capture.output(create_tail_feature_list(nanopolish = empty_tempfile,
                                                                 sequencing_summary = sequencing_summary,
                                                                 workspace = test_workspace,
                                                                 num_cores=2,
                                                                 basecall_group = wrong_basecall_group,
                                                                 pass_only=TRUE))))
  #check if it returns a nonempty file
  expect_type(invisible(capture.output(create_tail_feature_list(nanopolish = nanopolish,
                                       sequencing_summary = sequencing_summary,
                                       workspace = test_workspace,
                                       num_cores=2,
                                       basecall_group = test_basecall_group,
                                       pass_only=TRUE))), "character")

})



test_that("split_with_overlaps_moved correctly produces the output", {



  # segment is character
  expect_error(split_with_overlaps_moved(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                        tail_feature_list = test_tail_feature_list,
                                        segment = "A",
                                        overlap = 50))
  # overlap larger than segment
  expect_error(split_with_overlaps_moved(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                         tail_feature_list = test_tail_feature_list,
                                         segment = 1,
                                         overlap = 50))

  # overlap not an integer
  expect_error(split_with_overlaps_moved(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                         tail_feature_list = test_tail_feature_list,
                                         segment = 1,
                                         overlap = 50.5))
  # segment not an integer
  expect_error(split_with_overlaps_moved(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                         tail_feature_list = test_tail_feature_list,
                                         segment = 100.5,
                                         overlap = 50))


  # expect list as output
  expect_type(split_with_overlaps_moved(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                         tail_feature_list = test_tail_feature_list,
                                         segment = 100,
                                         overlap = 50), "list")


  #expect length of 7
  expect_length(split_with_overlaps_moved(readname = "002c1032-528c-4d79-98fc-857f0f0bd566",
                                          tail_feature_list = test_tail_feature_list,
                                          segment = 100,
                                          overlap = 50), 7)

})


test_that("create_tail_chunk_list_moved correctly produces the output", {

  #num cores is character
  expect_error(create_tail_chunk_list_moved(tail_feature_list =test_tail_feature_list, num_cores="A"))

})



test_that("create_gasf correctly produces the output", {

  test_chunk <- test_tail_chunk_list_moved[[1]][[1]]
  test_chunk_empty <- vector()

  test_chunk_noninteger <- c(708, 675, 672, 666, 667, 712, 692, 684, 670, 678, 647, 682, 665, 673, 691, 695, 687, 694, 655, 680, 672, 676, 656, 663, 686,
                  669, 702, 677, 666, 691, 682, 694, "ABCD", 651, 674, 658, 650, 678, 673, 672, 664, 663, 648, 668, 672, 647, 644, 654, 670, 687,
                  695, 687, 687, 687, 687, 687, 691, 691, 687, 695, 691, 687, 691, 695, 691, 695, 687, 691, 687, 687, 687, 691, 695, 695, 695,
                  695, 687, 687, 691, 695, 695, 687, 691, 687, 695, 695, 687, 691, 691, 691, 687, 687, 691, 695, 687, 687, 695, 691, 695, 691)

  test_chunk_short <- c(100:110)


  #check if returns numeric matrix
  expect_type(create_gasf(test_chunk), "double")

  # noninteger vec as an input
  expect_error(create_gasf(test_chunk_noninteger))

  # too short vec as an input
  expect_error(create_gasf(test_chunk_short))

  # empty vec as an input
  expect_error(create_gasf(test_chunk_empty))


})

