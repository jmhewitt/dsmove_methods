context('RookDirections.cpp C++ methods')

test_that('Validating orientation gives expected results', {
  
  # vector dir. is in the first coordinate (0-indexed) points away from neg.
  expect_identical(
    c(0,1),
    TestRookOrientation(head = c(10, 10), tail = c(9, 10))
  )
  
  # vector dir. is in the first coordinate (0-indexed) points toward neg.
  expect_identical(
    c(0,0),
    TestRookOrientation(head = c(10, 10), tail = c(11, 10))
  )
  
  # vector dir. is in the second coordinate (0-indexed), points toward neg.
  expect_identical(
    c(1,0),
    TestRookOrientation(head = c(10, 10, 3), tail = c(10, 11, 3))
  )
  
})

test_that('Validating AR covariate computed as intended', {
  
  # vector dir. is in first coordinate (0-indexed), points away from neg.
  tail = c(9,10)
  head = c(10,10)
  
  # movement is orthogonal to vector direction
  expect_true(
    all(
      TestRookDot(head = head, tail = tail, nextHead = head + c(0,1)) == 0,
      TestRookDot(head = head, tail = tail, nextHead = head + c(0,-1)) == 0
    )
  )
  
  # movement is aligned with vector direction
  expect_equal(
    1,
    TestRookDot(head = head, tail = tail, nextHead = head + c(-1,0))
  )
  
  # movement is opposite with vector direction
  expect_equal(
    -1,
    TestRookDot(head = head, tail = tail, nextHead = head + c(1,0))
  )
  
})