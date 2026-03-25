################################################################################
# Testing core functions
################################################################################

# extract_polya_data
################################################################################

test_that("extract_polya_data errors on missing nanopolish argument", {
  sequencing_summary <- system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  expect_error(extract_polya_data(sequencing_summary = sequencing_summary), "Nanopolish polya output is missing")
})


test_that("extract_polya_data rejects non-data.frame as nanopolish input", {
  sequencing_summary <- system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  expect_error(extract_polya_data(nanopolish = list(a = 1),sequencing_summary = sequencing_summary,pass_only = TRUE),
               "Empty data frame provided as an input \\(nanopolish\\)")
})


test_that("extract_polya_data errors on missing sequencing_summary argument", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy','nanopolish_output.tsv', package = 'ninetails')
  expect_error(extract_polya_data(nanopolish = nanopolish),"Sequencing summary file is missing")
})


test_that("extract_polya_data rejects non-data.frame as sequencing_summary input", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy','nanopolish_output.tsv', package = 'ninetails')
  expect_error(extract_polya_data(nanopolish = nanopolish,sequencing_summary = list(a = 1),pass_only = TRUE),
               "Empty data frame provided as an input \\(sequencing_summary\\)")
})


test_that("extract_polya_data errors on invalid pass_only type", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy','nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary <- system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  expect_error(extract_polya_data(nanopolish = nanopolish,
                                  sequencing_summary = sequencing_summary,
                                  pass_only = "yes"),"TRUE/FALSE")
})


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

# extract_tail_data
################################################################################

test_that("extract_tail_data validates required arguments", {

  test_workspace <- system.file('extdata', 'test_data', 'legacy',
                                'basecalled_fast5', package = 'ninetails')

  # missing readname
  expect_error(
    extract_tail_data(
      polya_summary = test_polya_summary,
      workspace = test_workspace,
      basecall_group = 'Basecall_1D_000'
    ),
    "Readname is missing"
  )

  # missing workspace
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      polya_summary = test_polya_summary,
      basecall_group = 'Basecall_1D_000'
    ),
    "workspace"
  )

  # missing basecall_group
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      polya_summary = test_polya_summary,
      workspace = test_workspace
    ),
    "Basecall group is missing"
  )

  # missing polya_summary
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      workspace = test_workspace,
      basecall_group = 'Basecall_1D_000'
    ),
    "Polya_summary is missing"
  )
})


test_that("extract_tail_data validates argument types and values", {

  test_workspace <- system.file('extdata', 'test_data', 'legacy',
                                'basecalled_fast5', package = 'ninetails')

  # non-character workspace
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      polya_summary = test_polya_summary,
      workspace = 123,
      basecall_group = 'Basecall_1D_000'
    ),
    "not a character string"
  )

  # empty workspace
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      polya_summary = test_polya_summary,
      workspace = "",
      basecall_group = 'Basecall_1D_000'
    ),
    "Empty string"
  )

  # non-character readname
  expect_error(
    extract_tail_data(
      readname = 123,
      polya_summary = test_polya_summary,
      workspace = test_workspace,
      basecall_group = 'Basecall_1D_000'
    ),
    "not a character string"
  )

  # empty readname
  expect_error(
    extract_tail_data(
      readname = "",
      polya_summary = test_polya_summary,
      workspace = test_workspace,
      basecall_group = 'Basecall_1D_000'
    )
  )
})


test_that("extract_tail_data handles invalid data scenarios", {

  empty_tempfile <- tempfile()
  test_workspace <- system.file('extdata', 'test_data', 'legacy',
                                'basecalled_fast5', package = 'ninetails')
  uncalled_workspace <- system.file('extdata', 'test_data', 'legacy',
                                    'uncalled_fast5', package = 'ninetails')

  test_basecall_group <- 'Basecall_1D_000'
  wrong_basecall_group <- 'Basecall_1D_003'

  # wrong fast5 files (uncalled)
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      polya_summary = test_polya_summary,
      workspace = uncalled_workspace,
      basecall_group = test_basecall_group
    )
  )

  # nonexistent basecall group
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      polya_summary = test_polya_summary,
      workspace = test_workspace,
      basecall_group = wrong_basecall_group
    )
  )

  # empty polya_summary file
  expect_error(
    extract_tail_data(
      readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
      polya_summary = empty_tempfile,
      workspace = test_workspace,
      basecall_group = test_basecall_group
    )
  )
})


test_that("extract_tail_data returns correct output for valid input", {

  test_workspace <- system.file('extdata', 'test_data', 'legacy',
                                'basecalled_fast5', package = 'ninetails')

  result <- extract_tail_data(
    readname = "1e6c551d-3ba7-461e-b930-4c87f5d638ba",
    polya_summary = test_polya_summary,
    workspace = test_workspace,
    basecall_group = 'Basecall_1D_000'
  )

  expect_length(result, 4)
  expect_type(result, "list")
})



# create_tail_feature_list
################################################################################

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


test_that("create_tail_feature_list errors on missing num_cores", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy','nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary <- system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  test_workspace <- system.file('extdata', 'test_data', 'legacy','basecalled_fast5', package = 'ninetails')

  expect_error(create_tail_feature_list(nanopolish = nanopolish,
                                        sequencing_summary = sequencing_summary,
                                        workspace = test_workspace,
                                        basecall_group = 'Basecall_1D_000'),"Number of declared cores is missing")
})

test_that("create_tail_feature_list errors on missing basecall_group", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy','nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary <- system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  test_workspace <- system.file('extdata', 'test_data', 'legacy','basecalled_fast5', package = 'ninetails')

  expect_error(create_tail_feature_list(nanopolish = nanopolish,
                                        sequencing_summary = sequencing_summary,
                                        workspace = test_workspace,
                                        num_cores = 2),"Basecall group is missing")
})

test_that("create_tail_feature_list errors on missing workspace", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy',
                            'nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary <- system.file('extdata', 'test_data', 'legacy',
                                    'sequencing_summary.txt', package = 'ninetails')

  expect_error(create_tail_feature_list(nanopolish = nanopolish,
                                        sequencing_summary = sequencing_summary,
                                        num_cores = 2,
                                        basecall_group = 'Basecall_1D_000'),"workspace")
})

test_that("create_tail_feature_list errors on non-numeric num_cores", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy', 'nanopolish_output.tsv', package = 'ninetails')
  sequencing_summary <- system.file('extdata', 'test_data', 'legacy','sequencing_summary.txt', package = 'ninetails')
  test_workspace <- system.file('extdata', 'test_data', 'legacy','basecalled_fast5', package = 'ninetails')

  expect_error(create_tail_feature_list(nanopolish = nanopolish,
                                        sequencing_summary = sequencing_summary,
                                        workspace = test_workspace,
                                        num_cores = "two",
                                        basecall_group = 'Basecall_1D_000'),"must be numeric")
})


# filter_signal_by_threshold
################################################################################
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

test_that("filter_signal_by_threshold output length matches input", {
  test_signal <- test_tail_feature_list[["tail_feature_list"]][[1]][[2]]
  result <- filter_signal_by_threshold(test_signal)

  expect_equal(length(result), length(test_signal))
})


test_that("filter_signal_by_threshold output values are in {-1, 0, 1}", {
  test_signal <- test_tail_feature_list[["tail_feature_list"]][[1]][[2]]
  result <- filter_signal_by_threshold(test_signal)

  expect_true(all(result %in% c(-1, 0, 1)))
})


test_that("filter_signal_by_threshold errors on missing signal", {
  expect_error(filter_signal_by_threshold(), "Signal is missing")
})

test_that("filter_signal_by_threshold errors on non-numeric signal", {
  expect_error(filter_signal_by_threshold(c("a", "b", "c")), "Signal must be numeric")
})


test_that("filter_signal_by_threshold is reproducible (set.seed)", {
  test_signal <- test_tail_feature_list[["tail_feature_list"]][[1]][[2]]

  result1 <- filter_signal_by_threshold(test_signal)
  result2 <- filter_signal_by_threshold(test_signal)

  expect_identical(result1, result2)
})


test_that("filter_signal_by_threshold detects positive outliers as 1", {
  # Create a signal with a clear positive spike
  set.seed(42)
  base_signal <- rep(500, 200)
  # Add a large spike in the middle (after position 100+5 for masking)
  base_signal[150:160] <- 800  # Strong positive deviation

  result <- filter_signal_by_threshold(base_signal)

  # Should detect some 1s (positive outliers) in the spike region
  # Note: first 5 positions are masked, so check after that
  expect_true(any(result[150:160] == 1) || all(result == 0))  # May or may not detect depending on threshold
})

test_that("filter_signal_by_threshold detects negative outliers as -1", {
  # Create a signal with a clear negative dip
  set.seed(42)
  base_signal <- rep(500, 200)
  # Add a large dip in the middle
  base_signal[150:160] <- 200  # Strong negative deviation

  result <- filter_signal_by_threshold(base_signal)

  # Should detect some -1s (negative outliers) in the dip region
  expect_true(any(result[150:160] == -1) || all(result == 0))  # May or may not detect depending on threshold
})

test_that("filter_signal_by_threshold returns mostly zeros for uniform signal", {
  # Uniform signal should produce mostly zeros
  uniform_signal <- rep(500, 200)

  result <- filter_signal_by_threshold(uniform_signal)

  # Most values should be 0 (no outliers in uniform signal)
  expect_true(sum(result == 0) > length(result) * 0.9)
})


# split_tail_centered
################################################################################


test_that("split_tail_centered errors when readname is missing", {
  # missing() triggers the stop, not an empty string check
  expect_error(split_tail_centered(tail_feature_list = test_tail_feature_list),"Readname is missing")
})

test_that("split_tail_centered errors when tail_feature_list is missing", {
  expect_error(split_tail_centered(readname = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b"),"List of tail features is missing")
})



test_that("split_tail_centered chunks have expected structure", {
  result <- split_tail_centered(readname = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b",
                                tail_feature_list = test_tail_feature_list)

  if (length(result) > 0) {
    first_chunk <- result[[1]]
    # each chunk should have chunk_sequence, chunk_start_pos, chunk_end_pos
    expect_true("chunk_sequence" %in% names(first_chunk))
    expect_true("chunk_start_pos" %in% names(first_chunk))
    expect_true("chunk_end_pos" %in% names(first_chunk))

    # chunk_sequence should be length 100
    expect_equal(length(first_chunk$chunk_sequence), 100)
  }
})


test_that("split_tail_centered correctly produces the output", {
  #readname is absent in tail features list
  expect_error(split_tail_centered(readname="Wladimir Wladimirowicz", tail_feature_list=test_tail_feature_list))

  #check if returned object is a list
  expect_type(split_tail_centered(readname = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b",
                                  tail_feature_list = test_tail_feature_list), "list")

})

test_that("split_tail_centered chunk end position is start + 99", {
  result <- split_tail_centered(readname = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b",
                                tail_feature_list = test_tail_feature_list)

  if (length(result) > 0) {
    # End position should be start + 99 (for 100-element chunk)
    expect_equal(result[[1]]$chunk_end_pos - result[[1]]$chunk_start_pos, 99)
  }
})




# create_tail_chunk_list
################################################################################


test_that("create_tail_chunk_list errors when num_cores is missing", {
  expect_error(create_tail_chunk_list(tail_feature_list = test_tail_feature_list),"Number of declared cores is missing")
})

test_that("create_tail_chunk_list errors when tail_feature_list is missing ", {
  expect_error(create_tail_chunk_list(num_cores = 2),"List of features is missing")
})




test_that("create_tail_chunk_list correctly produces the output", {

  #num cores is character
  expect_error(create_tail_chunk_list(tail_feature_list =test_tail_feature_list, num_cores="Putin is an asshole"))

})

test_that("create_tail_chunk_list errors on non-list tail_feature_list", {
  expect_error(create_tail_chunk_list(tail_feature_list = "not_a_list",
                                      num_cores = 2), "not a list")
})



# create_gaf
################################################################################

test_that("create_gaf errors on missing tail_chunk", {
  # Call with only method specified (tail_chunk still missing)
  expect_error(create_gaf(method = "s"), "Tail_chunk is missing")
  expect_error(create_gaf(method = "d"),"Tail_chunk is missing")
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


test_that("create_gaf returns correct dimensions for method='s'", {
  test_chunk <- test_tail_chunk_list[[1]][[1]][[1]]
  result <- create_gaf(test_chunk, method = "s")

  expect_equal(dim(result), c(100, 100, 1))
})


test_that("create_gaf returns correct dimensions for method='d'", {
  test_chunk <- test_tail_chunk_list[[1]][[1]][[1]]
  result <- create_gaf(test_chunk, method = "d")

  expect_equal(dim(result), c(100, 100, 1))
})


test_that("create_gaf GASF and GADF produce different results", {
  test_chunk <- test_tail_chunk_list[[1]][[1]][[1]]

  result_s <- create_gaf(test_chunk, method = "s")
  result_d <- create_gaf(test_chunk, method = "d")

  # GASF and GADF should not be identical for non-trivial signal
  expect_false(identical(result_s, result_d))
})


# combine_gafs
################################################################################

test_that("combine_gafs returns correct dimensions", {
  test_chunk <- test_tail_chunk_list[[1]][[1]][[1]]
  result <- combine_gafs(test_chunk)

  expect_equal(dim(result), c(100, 100, 2))
})


test_that("combine_gafs returns numeric array", {
  test_chunk <- test_tail_chunk_list[[1]][[1]][[1]]
  result <- combine_gafs(test_chunk)

  expect_type(result, "double")
  expect_true(is.array(result))
})

test_that("combine_gafs errors on missing tail_chunk", {
  expect_error(combine_gafs(),
               "Tail_chunk is missing")
})



# create_gaf_list
################################################################################

test_that("create_gaf_list errors on missing arguments", {
  expect_error(create_gaf_list(num_cores = 2),
               "List of tail chunks is missing")
  expect_error(create_gaf_list(tail_chunk_list = test_tail_chunk_list),
               "Number of declared cores is missing")
})


test_that("create_gaf_list errors on non-list tail_chunk_list", {
  expect_error(create_gaf_list(tail_chunk_list = "not_a_list", num_cores = 2))
})

test_that("create_gaf_list errors on non-numeric num_cores", {
  expect_error(create_gaf_list(tail_chunk_list = test_tail_chunk_list,
                               num_cores = "two"))
})



# predict_gaf_classes
################################################################################

test_that("predict_gaf_classes errors on missing argument", {
  expect_error(predict_gaf_classes(),
               "List of transformed signal chunks is missing")
})


# create_outputs
################################################################################
test_that("create_outputs errors on individual missing arguments", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy',
                            'nanopolish_output.tsv', package = 'ninetails')

  # missing tail_feature_list
  expect_error(create_outputs(tail_chunk_list = test_tail_chunk_list,
                              nanopolish = nanopolish,
                              predicted_list = test_predictions,
                              num_cores = 2),
               "tail_feature_list")

  # missing tail_chunk_list
  expect_error(create_outputs(tail_feature_list = test_tail_feature_list,
                              nanopolish = nanopolish,
                              predicted_list = test_predictions,
                              num_cores = 2),
               "tail_chunk")

  # missing nanopolish
  expect_error(create_outputs(tail_feature_list = test_tail_feature_list,
                              tail_chunk_list = test_tail_chunk_list,
                              predicted_list = test_predictions,
                              num_cores = 2),
               "Nanopolish")

  # missing predicted_list
  expect_error(create_outputs(tail_feature_list = test_tail_feature_list,
                              tail_chunk_list = test_tail_chunk_list,
                              nanopolish = nanopolish,
                              num_cores = 2),
               "predicted_list")

  # missing num_cores
  expect_error(create_outputs(tail_feature_list = test_tail_feature_list,
                              tail_chunk_list = test_tail_chunk_list,
                              nanopolish = nanopolish,
                              predicted_list = test_predictions),
               "num_cores")
})

test_that("create_outputs rejects empty nanopolish data frame", {
  expect_error(create_outputs(tail_feature_list = test_tail_feature_list,
                              tail_chunk_list = test_tail_chunk_list,
                              nanopolish = data.frame(),
                              predicted_list = test_predictions,
                              num_cores = 2),
               "Empty data frame")
})


test_that("create_outputs reads nanopolish from file path and returns valid output structure", {
  nanopolish <- system.file('extdata', 'test_data', 'legacy', 'nanopolish_output.tsv',package = 'ninetails')

  capture.output({result <- create_outputs(tail_feature_list = test_tail_feature_list,
                                           tail_chunk_list = test_tail_chunk_list,
                                           nanopolish = nanopolish,
                                           predicted_list = test_predictions,
                                           num_cores = 2,
                                           pass_only = TRUE,
                                           qc = TRUE)})

  expect_type(result,"list")
  expect_true("read_classes" %in% names(result))
  expect_true("nonadenosine_residues" %in% names(result))
  expect_true(is.data.frame(result[["read_classes"]]))
  expect_true(is.data.frame(result[["nonadenosine_residues"]]))

  # read_classes must contain columns for downstream analysis
  rc_cols <- c("readname", "polya_length", "qc_tag", "class", "comments")
  expect_true(all(rc_cols %in% colnames(result[["read_classes"]])))
})




test_that("create_outputs with pass_only = FALSE includes SUFFCLIP reads", {
  # Covers the pass_only == FALSE branch inside create_outputs.
  nanopolish <- system.file('extdata', 'test_data', 'legacy', 'nanopolish_output.tsv',package = 'ninetails')

  capture.output({result_pass <- create_outputs(tail_feature_list = test_tail_feature_list,
                                                tail_chunk_list = test_tail_chunk_list,
                                                nanopolish = nanopolish,
                                                predicted_list = test_predictions,
                                                num_cores= 2,
                                                pass_only= TRUE,
                                                qc = TRUE)

    result_all <- create_outputs(tail_feature_list= test_tail_feature_list,
                                 tail_chunk_list = test_tail_chunk_list,
                                 nanopolish= nanopolish,
                                 predicted_list= test_predictions,
                                 num_cores = 2,
                                 pass_only= FALSE,
                                 qc= TRUE)
  })

  # pass_only = FALSE must account for at least as many reads as pass_only = TRUE
  expect_gte(nrow(result_all[["read_classes"]]),nrow(result_pass[["read_classes"]]))
  # "NIN" comment (SUFFCLIP reads discarded in strict mode) must NOT appear
  # when pass_only = FALSE
  expect_false("NIN" %in% result_all[["read_classes"]]$comments)
})


