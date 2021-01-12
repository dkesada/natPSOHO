# Also covers Causlist initialization
test_that("velocity initialization works", {
  ordering <- c("A_t_0", "B_t_0", "C_t_0")
  size <- 3

  vl <- natVelocity$new(ordering, size)
  res <- c(0,0,0,0,0,0,0,0,0)

  expect_equal(vl$get_cl(), res)
  expect_equal(vl$get_cl_neg(), res)
})

test_that("random velocity generation works", {
  ordering <- c("A_t_0", "B_t_0", "C_t_0")
  size <- 3

  vl <- natVelocity$new(ordering, size)
  set.seed(42)
  vl$randomize_velocity(c(15,60,25))
  res_cl <- c(0,2,0,0,3,0,0,0,0)
  res_cl_neg <- c(3,0,0,0,0,0,3,0,0)
  res_n_ops <- 7

  expect_equal(vl$get_cl(), res_cl)
  expect_equal(vl$get_cl_neg(), res_cl_neg)
  expect_equal(vl$get_abs_op(), res_n_ops)
})

test_that("velocity addition works", {
  ordering <- c("A_t_0", "B_t_0", "C_t_0")
  size <- 3

  vl1 <- natVelocity$new(ordering, size)
  set.seed(42)
  vl1$randomize_velocity(c(15,60,25))

  vl2 <- natVelocity$new(ordering, size)
  set.seed(43)
  vl2$randomize_velocity(c(15,60,25))

  vl1$add_velocity(vl2)

  res_cl <- c(0,2,0,0,2,0,0,0,0)
  res_cl_neg <- c(3,0,0,0,0,0,3,0,0)

  expect_equal(vl1$get_cl(), res_cl)
  expect_equal(vl1$get_cl_neg(), res_cl_neg)
  expect_equal(vl1$get_abs_op(), 6)
})

test_that("cte times velocity works", {
  ordering <- c("A_t_0", "B_t_0", "C_t_0")
  size <- 3

  vl <- natVelocity$new(ordering, size)
  set.seed(42)
  vl$randomize_velocity(c(15,60,25))
  
  vl$cte_times_velocity(1.3)

  res_cl <- c(0,2,0,0,3,0,0,0,2)
  res_cl_neg <- c(3,0,0,0,0,1,3,0,0)

  expect_equal(vl$get_cl(), res_cl)
  expect_equal(vl$get_cl_neg(), res_cl_neg)
  expect_equal(vl$get_abs_op(), 9)
  
  vl$cte_times_velocity(-1)

  res_cl <- c(3,0,0,0,0,1,3,0,0)
  res_cl_neg <- c(0,2,0,0,3,0,0,0,2)

  expect_equal(vl$get_cl(), res_cl)
  expect_equal(vl$get_cl_neg(), res_cl_neg)
  expect_equal(vl$get_abs_op(), 9)

  vl$cte_times_velocity(-5.3)

  res_cl <- c(0,3,3,0,3,2,0,0,2)
  res_cl_neg <- c(3,0,0,3,0,1,3,3,1)

  expect_equal(vl$get_cl(), res_cl)
  expect_equal(vl$get_cl_neg(), res_cl_neg)
  expect_equal(vl$get_abs_op(), 18)

  vl$cte_times_velocity(-0.3)

  res_cl <- c(1,0,0,1,0,0,1,0,1)
  res_cl_neg <- c(0,0,0,0,1,0,0,0,0)
  
  expect_equal(vl$get_cl(), res_cl)
  expect_equal(vl$get_cl_neg(), res_cl_neg)
  expect_equal(vl$get_abs_op(), 5)

  vl$cte_times_velocity(0.1)

  res_cl <- c(0,0,0,0,0,0,0,0,0)
  res_cl_neg <- c(0,0,0,0,0,0,0,0,0)
  
  expect_equal(vl$get_cl(), res_cl)
  expect_equal(vl$get_cl_neg(), res_cl_neg)
  expect_equal(vl$get_abs_op(), 0)
})

