context("Trait Imputations")

test_that("imputations- scOU - random root", {
  testthat::skip_on_cran()
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
  
  res_heavy <- PhyloEM(phylo = tree,
                       Y_data = traits,
                       process = "scOU",
                       K_max = 5,
                       random.root = TRUE,
                       stationary.root = TRUE,
                       alpha = selection.strength,
                       save_step = FALSE,
                       Nbr_It_Max = 2000,
                       method.variance = "upward_downward",
                       method.init = "lasso",
                       use_previous = FALSE,
                       method.selection = c("BGHlsq", "BGHml", "BGHmlraw", "BGHlsqraw"),
                       progress.bar = FALSE,
                       K_lag_init = 2,
                       light_result = FALSE,
                       check.tips.names = FALSE)
  
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength, 
                              trait = 1:2, where = "tips", what = "expectations"),
               res_heavy$alpha_max$m_Y_estim$'4')
  
  expect_equal(sum(residuals(res_heavy, K=4)^2),
               res_heavy$alpha_max$results_summary$least_squares_raw[5])
  
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength,
                              trait = 1:2, where = "tips", what = "variances"),
               res_heavy$alpha_max$Yvar$'4')
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength,
                              trait = 1:2, where = "tips", what = "imput"),
               res_heavy$alpha_max$Yhat$'4')
  
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength,
                              trait = 1:2, where = "nodes", what = "variances"),
               res_heavy$alpha_max$Zvar$'4')
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength,
                              trait = 1:2, where = "nodes", what = "imput"),
               res_heavy$alpha_max$Zhat$'4')
  
  ## Light
  res_light <- enlight(res_heavy)
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength,
                              trait = 1:2, where = "tips", what = "expectations"),
               res_heavy$alpha_max$m_Y_estim$'4')
  
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength,
                              trait = 1:2, where = "tips", what = "variances"),
               res_heavy$alpha_max$Yvar$'4')
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength,
                              trait = 1:2, where = "tips", what = "imput"),
               res_heavy$alpha_max$Yhat$'4')
  
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength,
                              trait = 1:2, where = "nodes", what = "variances"),
               res_heavy$alpha_max$Zvar$'4')
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength,
                              trait = 1:2, where = "nodes", what = "imput"),
               res_heavy$alpha_max$Zhat$'4')
  
  ## Select
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHlsq"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHlsq"))
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHml"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHml"))
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHlsqraw"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHlsqraw"))
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHmlraw"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips",
                              what = "expectations", model.selection = "BGHmlraw"))
  
  ## rBM
  pp <- params_process(res_light, K = 3, alpha = 3, rBM = TRUE)
  pp$selection.strength <- NULL
  pp$process <- NULL
  class(pp) <- NULL
  rr <- compute_raw_parameters(res_light$phylo, 
                               res_light$alpha_3$params_estim$`3`)
  expect_equal(pp, rr)
})

test_that("imputations- scOU - fixed root", {
  testthat::skip_on_cran()
  set.seed(17920902)
  ntaxa = 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 2
  variance <- diag(0.2, p, p) +  matrix(0.8, p, p)
  selection.strength <- 3
  independent <- FALSE
  root.state <- list(random = FALSE,
                     stationary.root = FALSE,
                     value.root = rep(1, p),
                     exp.root = NA,
                     var.root = NA)
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
  
  res_heavy <- PhyloEM(phylo = tree,
                       Y_data = traits,
                       process = "scOU",
                       K_max = 5,
                       random.root = FALSE,
                       stationary.root = FALSE,
                       alpha = selection.strength,
                       save_step = FALSE,
                       Nbr_It_Max = 2000,
                       method.variance = "upward_downward",
                       method.init = "lasso",
                       use_previous = FALSE,
                       method.selection = c("BGHlsq", "BGHml", "BGHmlraw", "BGHlsqraw"),
                       progress.bar = FALSE,
                       K_lag_init = 2,
                       light_result = FALSE,
                       check.tips.names = FALSE)
  
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength, trait = 1:2, where = "tips", what = "expectations"),
               res_heavy$alpha_max$m_Y_estim$'4')
  
  expect_equal(sum(residuals(res_heavy, K=4)^2),
               res_heavy$alpha_max$results_summary$least_squares_raw[5])
  
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength, trait = 1:2, where = "tips", what = "variances"),
               res_heavy$alpha_max$Yvar$'4')
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength, trait = 1:2, where = "tips", what = "imput"),
               res_heavy$alpha_max$Yhat$'4')
  
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength, trait = 1:2, where = "nodes", what = "variances"),
               res_heavy$alpha_max$Zvar$'4')
  expect_equal(imputed_traits(res_heavy, K=4, alpha = selection.strength, trait = 1:2, where = "nodes", what = "imput"),
               res_heavy$alpha_max$Zhat$'4')
  
  ## Light
  res_light <- enlight(res_heavy)
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength, trait = 1:2, where = "tips", what = "expectations"),
               res_heavy$alpha_max$m_Y_estim$'4')
  
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength, trait = 1:2, where = "tips", what = "variances"),
               res_heavy$alpha_max$Yvar$'4')
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength, trait = 1:2, where = "tips", what = "imput"),
               res_heavy$alpha_max$Yhat$'4')
  
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength, trait = 1:2, where = "nodes", what = "variances"),
               res_heavy$alpha_max$Zvar$'4')
  expect_equal(imputed_traits(res_light, K=4, alpha = selection.strength, trait = 1:2, where = "nodes", what = "imput"),
               res_heavy$alpha_max$Zhat$'4')
  
  ## Select
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHlsq"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHlsq"))
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHml"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHml"))
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHlsqraw"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHlsqraw"))
  expect_equal(imputed_traits(res_light, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHmlraw"),
               imputed_traits(res_heavy, trait = 1:2, where = "tips", what = "expectations", model.selection = "BGHmlraw"))
})