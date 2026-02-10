# testing postprocessing fun

test_that("create_outputs works correctly", {

  empty_tempfile = tempfile()
  tail_feature_list = test_tail_feature_list
  tail_chunk_list = test_tail_chunk_list
  predicted_list = test_predictions
  num_cores=2
  bad_num_cores= "ABC"
  nanopolish = system.file('extdata', 'test_data','legacy', 'nanopolish_output.tsv', package = 'ninetails')

  #tail feature list wrong
  expect_error(create_outputs(tail_feature_list= empty_tempfile,
                              tail_chunk_list=tail_chunk_list,
                              nanopolish=nanopolish,
                              predicted_list = predicted_list,
                              num_cores=num_cores))

  # check if there is an expected output
  expect_type(create_outputs(tail_feature_list= tail_feature_list,
                             tail_chunk_list=tail_chunk_list,
                             nanopolish=nanopolish,
                             predicted_list = predicted_list,
                             num_cores=num_cores), "list")

})

