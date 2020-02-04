context("Log Likelihood")

test_that("log-likelihood - scOU - random root", {
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
  
  res_new <- PhyloEM(phylo = tree,
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
                     method.selection = "LINselect",
                     progress.bar = FALSE,
                     K_lag_init = 2,
                     check.tips.names = FALSE)
  
  expect_equal(log_likelihood(res_new, K = 1, alpha = "max"),
               res_new$alpha_max$results_summary$log_likelihood[2])
  
  ## Test output
  params <- params_process.PhyloEM(res_new, K = 0, alpha = "max")
  expect_output(print(params),
                "This process has no shift.")
  expect_output(print(params),
                "2 dimensional scOU process with a random stationary root.")
  
  ## tree quantities
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  T_tree <- incidence.matrix(tree)
  ## params and correlation tree matrix
  params <- params_process.PhyloEM(res_new, K = 1, alpha = "max")
  params$selection.strength <- unique(diag(params$selection.strength))
  C <- compute_tree_correlations_matrix.scOU(times_shared, distances_phylo, params)
  C <- 1/(2*selection.strength) * extract.variance_covariance(C, what="YY",
                                   masque_data = c(rep(TRUE, ntaxa),
                                                   rep(FALSE, dim(C)[1] - ntaxa)))
  C_inv <- solve(C)
  ## regression matrix
  ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                    selection.strength = as.vector(params$selection.strength),
                                                    times_shared = times_shared)
  Tr <- T_tree * ac_tree
  Tr <- Tr[, params$shifts$edges, drop = F]
  Tr <- cbind(Tr, rep(1, dim(Tr)[1]))
  
  P <- C_inv %*% Tr %*% solve(t(Tr) %*% C_inv %*% Tr) %*% t(Tr) %*% C_inv
  ## Max ll parameters
  R_max <- traits %*% (C_inv - P) %*% t(traits) / ntaxa
  Delta_max <- solve(t(Tr) %*% C_inv %*% Tr) %*% t(Tr) %*% C_inv %*% t(traits)
  
  expect_equal(Delta_max[2, ], params$root.state$exp.root)
  expect_equal(Delta_max[1, ], params$shifts$values[, 1])
  expect_equal(as.matrix(R_max),
               as.matrix(params$variance), tolerance = 10^(-1))
  
  ## Max likelihood
  LL_max <- -ntaxa*p/2*(log(2*pi)+1) - ntaxa/2*determinant(R_max, logarithm = TRUE)$modulus - p/2 * determinant(C, logarithm = T)$modulus
  expect_equal(log_likelihood(res_new, K = 1, alpha = "max"),
               as.vector(LL_max), tolerance = 10^(-2))
  
  ## Least squares
  lsq_1 <- sum(apply(t(traits) - Tr %*% Delta_max, 2, function(z) as.vector(t(z) %*% C_inv %*% z)))
  lsq_2 <- ntaxa * sum(diag(R_max))
  expect_equal(lsq_1, lsq_2)
})

test_that("log-likelihood - BM - p=1 - fixed root", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920902)
  ntaxa = 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 1
  variance <- diag(0.2, p, p) +  matrix(0.8, p, p)
  independent <- FALSE
  root.state <- list(random = FALSE,
                     value.root = rep(1, p),
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(25, 10, 31),
                values = cbind(c(5),
                               c(-5),
                               c(3)),
                relativeTimes = 0)
  paramsSimu <- list(variance = variance,
                     shifts = shifts,
                     root.state = root.state)
  attr(paramsSimu, "p_dim") <- p
  
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "BM",
                          variance = variance,
                          shifts = shifts)
  
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  
  res_new <- PhyloEM(phylo = tree,
                     Y_data = traits,
                     process = "BM",
                     K_max = 5,
                     random.root = FALSE,
                     save_step = FALSE,
                     Nbr_It_Max = 2000,
                     method.variance = "upward_downward",
                     method.init = "lasso",
                     use_previous = FALSE,
                     method.selection = "LINselect",
                     progress.bar = FALSE,
                     K_lag_init = 2,
                     check.tips.names = FALSE)
  
  expect_equal(log_likelihood(res_new, K = 1, alpha = "max"),
               res_new$alpha_max$results_summary$log_likelihood[2])
  
  ## Test output
  params <- params_process.PhyloEM(res_new, K = 0, alpha = "max")
  expect_output(print(params),
                "This process has no shift.")
  expect_output(print(params),
                "1 dimensional BM process with a fixed root.")
  
  ## tree quantities
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  T_tree <- incidence.matrix(tree)
  ## params and correlation tree matrix
  params <- params_process.PhyloEM(res_new, K = 1, alpha = "max")
  params$selection.strength <- unique(diag(params$selection.strength))
  C <- compute_tree_correlations_matrix.BM(times_shared, params)
  C <- extract.variance_covariance(C, what="YY",
                                   masque_data = c(rep(TRUE, ntaxa),
                                                   rep(FALSE, dim(C)[1] - ntaxa)))
  C_inv <- solve(C)
  ## regression matrix
  ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                    selection.strength = as.vector(params$selection.strength),
                                                    times_shared = times_shared)
  Tr <- T_tree * ac_tree
  Tr <- Tr[, params$shifts$edges, drop = F]
  Tr <- cbind(Tr, rep(1, dim(Tr)[1]))
  
  P <- C_inv %*% Tr %*% solve(t(Tr) %*% C_inv %*% Tr) %*% t(Tr) %*% C_inv
  ## Max ll parameters
  R_max <- traits %*% (C_inv - P) %*% t(traits) / ntaxa
  Delta_max <- solve(t(Tr) %*% C_inv %*% Tr) %*% t(Tr) %*% C_inv %*% t(traits)
  
  expect_equal(Delta_max[2, ], params$root.state$value.root)
  expect_equal(Delta_max[1, ], params$shifts$values[, 1])
  expect_equal(as.matrix(R_max),
               as.matrix(params$variance), tolerance = 10^(-1))
  
  ## Max likelihood
  LL_max <- -ntaxa*p/2*(log(2*pi)+1) - ntaxa/2*determinant(R_max, logarithm = TRUE)$modulus - p/2 * determinant(C, logarithm = T)$modulus
  expect_equal(log_likelihood(res_new, K = 1, alpha = "max"),
               as.vector(LL_max), tolerance = 10^(-2))
  
  ## Least squares
  lsq_1 <- sum(apply(t(traits) - Tr %*% Delta_max, 2, function(z) as.vector(t(z) %*% C_inv %*% z)))
  lsq_2 <- ntaxa * sum(diag(R_max))
  expect_equal(lsq_1, lsq_2)
})

# context("Log Likelihood Computation")
# 
# test_that("independent computations", {
#   set.seed(586)
#   ntaxa <- 50
#   tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
#                           lambda = 1, mu = 0,
#                           age = 1, mrca = TRUE)[[1]]
#   
#   p <- 3
#   variance <- diag(0.5, p, p)
#   optimal.value <- c(-3, 5, 0)
#   selection.strength <- diag(3, p, p)
#   exp.stationary <- optimal.value
#   var.stationary  <- compute_stationary_variance(variance, selection.strength)
#   root.state <- list(random = TRUE,
#                      stationary.root = TRUE,
#                      value.root = NA,
#                      exp.root = exp.stationary,
#                      var.root = var.stationary)
#   shifts = list(edges = c(18, 32),
#                 values=cbind(c(4, -10, 3),
#                              c(-5, 5, 0)),
#                 relativeTimes = 0)
#   paramsSimu <- list(variance = variance,
#                      optimal.value = optimal.value,
#                      selection.strength = selection.strength,
#                      shifts = shifts,
#                      root.state = root.state)
#   attr(paramsSimu, "p_dim") <- p
#   
#   X1 <- simulate_internal(tree,
#                  p = p,
#                  root.state = root.state,
#                  process = "OU",
#                  variance = variance,
#                  optimal.value = optimal.value,
#                  selection.strength = selection.strength,
#                  shifts = shifts)
#   Y_sim <- extract_simulate_internal(X1, "tips", "states")
#   
#   temps <- system.time(moments <- compute_mean_variance.simple(phylo = tree,
#                                           times_shared = compute_times_ca(tree),
#                                           distances_phylo = compute_dist_phy(tree),
#                                           process = "OU",
#                                           params_old = paramsSimu))
#   
#   temps <- temps + system.time(log_likelihood.true <- compute_log_likelihood.simple(phylo = tree,
#                                                        Y_data_vec = as.vector(Y_sim),
#                                                        sim = moments$sim,
#                                                        Sigma = moments$Sigma,
#                                                        Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv))
#   temps <- temps + system.time(maha_data_mean <- compute_mahalanobis_distance.simple(phylo = tree,
#                                                         Y_data_vec = as.vector(Y_sim),
#                                                         sim = moments$sim,
#                                                         Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv))
#   
#   temps2 <- system.time(tmp <- lik_maha_ind(phylo = tree,
#                times_shared = compute_times_ca(tree),
#                distances_phylo = compute_dist_phy(tree),
#                process = "OU",
#                params_old = paramsSimu,
#                independent = TRUE,
#                Y_sim))
#   
#   expect_that(log_likelihood.true, equals(sum(sapply(tmp, function(z) z$log_likelihood))))
#   expect_that(maha_data_mean, equals(sum(sapply(tmp, function(z) z$maha_data_mean))))
# })