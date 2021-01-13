test_that("position translation works", {
  net <- bnlearn::model2network("[A_t_2][B_t_2][C_t_2][A_t_1][B_t_1][C_t_1][A_t_0|A_t_1:B_t_2:C_t_1][B_t_0|A_t_1:B_t_1][C_t_0|B_t_2:C_t_2]")
  class(net) <- c("dbn", class(net))
  size <- 3

  ps <- natPosition$new(net, size)
  res <- c(1,2,1,1,1,0,0,2,2)
  expect_equal(ps$get_cl(), res)
})

test_that("random network generation works", {
  ordering <- c("A", "B", "C")
  size <- 3

  set.seed(51)
  ps <- natPosition$new(NULL, size, ordering)
  res <- c(2,0,0,2,0,2,2,1,2)

  expect_equal(ps$get_cl(), res)
})


test_that("translation from Position to bn works", {
  net <- bnlearn::model2network("[A_t_2][B_t_2][C_t_2][A_t_1][B_t_1][C_t_1][A_t_0|A_t_1:B_t_2:C_t_1][B_t_0|A_t_1:B_t_1][C_t_0|B_t_2:C_t_2]")
  class(net) <- c("dbn", class(net))
  size <- 3

  ps <- natPosition$new(net, size)
  res_net <- ps$bn_translate()

  res <- apply(net$arcs, 1, function(x){paste0(x[1], x[2])})
  res_net <- apply(res_net$arcs, 1, function(x){paste0(x[1], x[2])})


  expect_setequal(res_net, res)
})

test_that("position plus velocity works", {
  net <- bnlearn::model2network("[A_t_2][B_t_2][C_t_2][A_t_1][B_t_1][C_t_1][A_t_0|A_t_1:B_t_2:C_t_1][B_t_0|A_t_1:B_t_1][C_t_0|B_t_2:C_t_2]")
  class(net) <- c("dbn", class(net))
  ordering <- c("A_t_0", "B_t_0", "C_t_0")
  size <- 3

  ps <- natPosition$new(net, size)
  vl <- natVelocity$new(ordering, size)
  set.seed(42)
  vl$randomize_velocity(c(15,60,25))
  ps$add_velocity(vl)

  res <- c(0,2,1,1,3,0,0,2,2)
  expect_equal(ps$get_cl(), res)
  expect_equal(ps$get_n_arcs(), 7)
})

test_that("position minus position works", {
  net <- bnlearn::model2network("[A_t_2][B_t_2][C_t_2][A_t_1][B_t_1][C_t_1][A_t_0|A_t_1:B_t_2:C_t_1][B_t_0|A_t_1:B_t_1][C_t_0|B_t_2:C_t_2]")
  class(net) <- c("dbn", class(net))
  ordering <- c("A_t_0", "B_t_0", "C_t_0")
  size <- 3
  ps1 <- natPosition$new(net, size)

  ordering <- c("A", "B", "C")
  set.seed(51)
  ps2 <- natPosition$new(NULL, size, ordering)

  vl <- ps1$subtract_position(ps2)

  res_cl <- c(1,2,1,1,1,0,0,2,0)
  res_cl_neg <- c(2,0,0,2,0,2,2,1,0)

  expect_equal(vl$get_cl(), res_cl)
  expect_equal(vl$get_cl_neg(), res_cl_neg)
  expect_equal(vl$get_abs_op(), 11)
})
