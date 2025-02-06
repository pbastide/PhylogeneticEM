context("Function merge alpha grids")

###############################################################################
test_that("Merge alpha grids", {
  testthat::skip_on_cran()
  
  data(monkeys)
  
  ## First fit with coarse grid
  expect_message(expect_warning(
    res1 <- PhyloEM(Y_data = monkeys$dat,
                    phylo = monkeys$phy,
                    process = "scOU",
                    random.root = TRUE,
                    stationary.root = TRUE,
                    K_max = 10,
                    alpha = c(0.2, 0.3),
                    parallel_alpha = FALSE,
                    progress.bar = FALSE),
    "Models selected by DDSE and Djump are different"),
    "There are some equivalent solutions")
  
  ## Second fit with finner grid
  expect_message(
    res2 <- PhyloEM(Y_data = monkeys$dat,
                    phylo = monkeys$phy,
                    process = "scOU",
                    random.root = TRUE,
                    stationary.root = TRUE,
                    K_max = 10,
                    alpha = c(0.22, 0.24),
                    parallel_alpha = FALSE,
                    progress.bar = FALSE),
    "There are some equivalent solutions")
  
  ## Third fit with redundancies
  expect_message(
    res3 <- PhyloEM(Y_data = monkeys$dat,
                    phylo = monkeys$phy,
                    process = "scOU",
                    random.root = TRUE,
                    stationary.root = TRUE,
                    K_max = 10,
                    alpha = c(0.26, 0.3),
                    parallel_alpha = FALSE,
                    progress.bar = FALSE),
    "There are some equivalent solutions")
  
  ## merge results
  expect_message(expect_warning(
    res_merge <- merge_alpha_grids(res1, res2, res3),
    "Models selected by DDSE and Djump are different"),
    "There are some equivalent solutions")
  
  ## large fit all alpha
  expect_message(expect_warning(
    resall <- PhyloEM(Y_data = monkeys$dat,
                      phylo = monkeys$phy,
                      process = "scOU",
                      random.root = TRUE,
                      stationary.root = TRUE,
                      K_max = 10,
                      alpha = c(0.2, 0.22, 0.24, 0.26, 0.3),
                      parallel_alpha = FALSE,
                      progress.bar = FALSE),
    "Models selected by DDSE and Djump are different"),
    "There are some equivalent solutions")
  
  ## compare
  compare_ignore_time <- function(a, b) {
    a$results_summary$time <- 0
    a$DDSE_BM1$results_summary$time <- 0
    a$Djump_BM1$results_summary$time <- 0
    a$BGHml$results_summary$time <- 0
    a$BGHmlraw$results_summary$time <- 0
    a$BGHlsq$results_summary$time <- 0
    a$BGHlsqraw$results_summary$time <- 0
    
    b$results_summary$time <- 0
    b$DDSE_BM1$results_summary$time <- 0
    b$Djump_BM1$results_summary$time <- 0
    b$BGHml$results_summary$time <- 0
    b$BGHmlraw$results_summary$time <- 0
    b$BGHlsq$results_summary$time <- 0
    b$BGHlsqraw$results_summary$time <- 0
    
    expect_equal(a, b)
  }
  compare_ignore_time(resall$alpha_max, res_merge$alpha_max)
  compare_ignore_time(resall$alpha_min, res_merge$alpha_min)
  compare_ignore_time(resall$alpha_min_raw, res_merge$alpha_min_raw)
})

###############################################################################
test_that("Merge alpha grids - errors", {
  testthat::skip_on_cran()
  data(monkeys)
  ## First fit with coarse grid
  expect_message(expect_warning(
    res1 <- PhyloEM(Y_data = monkeys$dat,
                    phylo = monkeys$phy,
                    process = "scOU",
                    random.root = TRUE,
                    stationary.root = TRUE,
                    K_max = 10,
                    alpha = c(0.2, 0.3),
                    parallel_alpha = FALSE,
                    progress.bar = FALSE),
    "Models selected by DDSE and Djump are different"),
    "There are some equivalent solutions")
  
  ## Second fit with finner grid
  monkeys$dat[1, 10] <- 0.2
  expect_message(
    res2 <- PhyloEM(Y_data = monkeys$dat,
                    phylo = monkeys$phy,
                    process = "scOU",
                    random.root = TRUE,
                    stationary.root = TRUE,
                    K_max = 10,
                    alpha = c(0.22, 0.24),
                    parallel_alpha = FALSE,
                    progress.bar = FALSE),
    "There are some equivalent solutions")
  
  expect_error(merge_alpha_grids(), "There should be at least")
  expect_equal(merge_alpha_grids(res1), res1)
  expect_error(merge_alpha_grids(res1, res2), "Component Y_data of PhyloEM results number 1 and number 2 are different.")
})