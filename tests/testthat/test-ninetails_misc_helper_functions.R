# Testing helper functions

# testing count_overlaps
test_that("count_overlaps works correctly", {
  test_vector <- c(1:200)
  empty_vector <- vector()
  expect_equal(count_overlaps(test_vector, 10, 2), 25)
  expect_error(count_overlaps(test_vector, 2, 10))
  expect_error(count_overlaps(test_vector, 0, 2))
  expect_error(count_overlaps(test_vector, -1, 22))
  expect_error(count_overlaps(test_vector, 20, -2))
  expect_error(count_overlaps(test_vector, 20, c(1,2,3)))
  expect_error(count_overlaps(test_vector, c(10,20,30), 1))
  expect_error(count_overlaps(test_vector, "A", 1))
  expect_error(count_overlaps(empty_vector, 10, 2))

})


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



