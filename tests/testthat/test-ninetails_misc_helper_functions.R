# Testing helper functions

#testing winsorize_signal
test_that("winsorize_signal works correctly", {
  empty_vector <- vector()
  test_vector <- c(-1:11)
  test_vector_with_NAs <- c(-1:11, NA)
  test_vector_with_character <- c(-1:11,"A")
  expect_vector(winsorize_signal(test_vector), size=13)
  expect_equal(winsorize_signal(test_vector), c(0,0,1,2,3,4,5,6,7,8,9,10,10))
  expect_error(winsorize_signal(empty_vector))
  expect_error(winsorize_signal(test_vector_with_NAs))
  expect_error(winsorize_signal(test_vector_with_character))

})



