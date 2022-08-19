# Testing postprocessing functions


test_that("create_coordinate_dataframe works correctly", {

  empty_tempfile = tempfile()

  #tail feature list wrong
  expect_error(create_coordinate_dataframe(tail_feature_list= empty_tempfile,
                                           num_cores=2))
  # missing num cores arg
  expect_error(create_coordinate_dataframe(tail_feature_list = test_tail_feature_list))
  #num cores not numeric
  expect_error(create_coordinate_dataframe(tail_feature_list = test_tail_feature_list,
                                           num_cores="two"))
  # check if there is an expected output
  expect_type(create_coordinate_dataframe(tail_feature_list = test_tail_feature_list,
                                          num_cores=2), "list")

})



test_that("analyze_results works correctly", {
  empty_tempfile = tempfile()

  # nanopolish broken
  expect_error(analyze_results(nanopolish=empty_tempfile,
                               coordinate_df=test_coordinate_df,
                               predicted_list=test_predictions,
                               pass_only = TRUE))
  # predicted list broken
  expect_error(analyze_results(nanopolish=nanopolish,
                               coordinate_df=test_coordinate_df,
                               predicted_list=empty_tempfile,
                               pass_only = TRUE))
  # coordinate df broken
  expect_error(analyze_results(nanopolish=nanopolish,
                               coordinate_df=empty_tempfile,
                               predicted_list=test_predictions,
                               pass_only = TRUE))
  # check if output is a list
  expect_type(analyze_results(nanopolish=nanopolish,
                              coordinate_df=test_coordinate_df,
                              predicted_list=test_predictions,
                              pass_only = TRUE), "list")
  # check if output is a list of length 2 (list of 2 lists)
  expect_length(analyze_results(nanopolish=nanopolish,
                              coordinate_df=test_coordinate_df,
                              predicted_list=test_predictions,
                              pass_only = TRUE), 2)


})











