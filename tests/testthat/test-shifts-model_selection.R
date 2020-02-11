context("Model Selection")

test_that("Model Selection", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920902)
  ntaxa = 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 2
  variance <- diag(0.2, p, p) +  matrix(0.8, p, p)
  selection.strength <- 3
  independent <- FALSE
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = rep(1, p),
                     var.root = compute_stationary_variance(variance, selection.strength))
  shifts = list(edges = c(25, 10, 31),
                values = cbind(c(5, 5),
                               c(-5, -5),
                               c(3, 3)),
                relativeTimes = 0)
  paramsSimu <- list(variance = variance,
                     shifts = shifts,
                     root.state = root.state,
                     selection.strength = selection.strength,
                     optimal.value = rep(1, p))
  attr(paramsSimu, "p_dim") <- p
  
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "scOU",
                          variance = variance,
                          shifts = shifts,
                          selection.strength = selection.strength,
                          optimal.value = paramsSimu$optimal.value)
  
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  
  expect_warning(res_slope <- PhyloEM(phylo = tree,
                                      Y_data = traits,
                                      process = "scOU",
                                      K_max = 10,
                                      random.root = TRUE,
                                      stationary.root = TRUE,
                                      alpha = selection.strength,
                                      save_step = FALSE,
                                      Nbr_It_Max = 2,
                                      method.variance = "upward_downward",
                                      method.init = "lasso",
                                      use_previous = FALSE,
                                      method.selection = "DDSE",
                                      progress.bar = FALSE,
                                      K_lag_init = 0,
                                      light_result = TRUE))
  
  expect_warning(res_LIN <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "scOU",
                                    K_max = 10,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    alpha = selection.strength,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "upward_downward",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "LINselect",
                                    progress.bar = FALSE,
                                    K_lag_init = 0,
                                    light_result = TRUE))
  
  res <- model_selection(res_LIN, method.selection = "DDSE")
  res_bis <- model_selection(res_slope, method.selection = "LINselect")
  
  ## Time is different
  res$alpha_3$results_summary$time <- 0
  res_bis$alpha_3$results_summary$time <- 0
  res$alpha_max$results_summary$time <- 0
  res_bis$alpha_max$results_summary$time <- 0
  res$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_bis$alpha_max$DDSE_BM1$results_summary$time <- 0
  res$alpha_max$Djump_BM1$results_summary$time <- 0
  res_bis$alpha_max$Djump_BM1$results_summary$time <- 0
  res$alpha_max$BGHml$results_summary$time <- 0
  res_bis$alpha_max$BGHml$results_summary$time <- 0
  res$alpha_max$BGHmlraw$results_summary$time <- 0
  res_bis$alpha_max$BGHmlraw$results_summary$time <- 0
  res$alpha_min$results_summary$time <- 0
  res_bis$alpha_min$results_summary$time <- 0
  res$alpha_min_raw$results_summary$time <- 0
  res_bis$alpha_min_raw$results_summary$time <- 0
  res$alpha_min$BGHlsq$results_summary$time <- 0
  res_bis$alpha_min$BGHlsq$results_summary$time <- 0
  res$alpha_min_raw$BGHlsqraw$results_summary$time <- 0
  res_bis$alpha_min_raw$BGHlsqraw$results_summary$time <- 0
  ## Order is different
  res_bis$alpha_max <- res_bis$alpha_max[names(res$alpha_max)]
  res_bis$alpha_max$results_summary <- res_bis$alpha_max$results_summary[names(res$alpha_max$results_summary)]
  res_bis$alpha_max$K_select <- res_bis$alpha_max$K_select[names(res$alpha_max$K_select)]
  res_bis$alpha_max$BGHml$results_summary <- res_bis$alpha_max$BGHml$results_summary[names(res$alpha_max$BGHml$results_summary)]
  res_bis$alpha_max$BGHmlraw$results_summary <- res_bis$alpha_max$BGHmlraw$results_summary[names(res$alpha_max$BGHmlraw$results_summary)]
  res$alpha_max$DDSE_BM1$results_summary <- res$alpha_max$DDSE_BM1$results_summary[names(res_bis$alpha_max$DDSE_BM1$results_summary)]
  res$alpha_max$Djump_BM1$results_summary <- res$alpha_max$Djump_BM1$results_summary[names(res_bis$alpha_max$Djump_BM1$results_summary)]
  
  expect_equal(res, res_bis)
})