test_that("velocity initialization works", {
  ordering <- c("A_t_0", "B_t_0", "C_t_0")
  size <- 3
  
  vl <- Velocity$new(ordering, size)
  res <- list(
    list(
      list(c("A_t_1", "B_t_1", "C_t_1"),
           c(0,0,0)),
      list(c("A_t_1", "B_t_1", "C_t_1"),
           c(0,0,0)),
      list(c("A_t_1", "B_t_1", "C_t_1"),
           c(0,0,0))),
    list(
      list(c("A_t_2", "B_t_2", "C_t_2"),
           c(0,0,0)),
      list(c("A_t_2", "B_t_2", "C_t_2"),
           c(0,0,0)),
      list(c("A_t_2", "B_t_2", "C_t_2"),
           c(0,0,0)))
  )
  
  expect_equal(vl$get_vl(), res)
})