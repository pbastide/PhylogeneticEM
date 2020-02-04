context("E step upward-downward")

test_that("Upward Downward - BM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920921)
  ntaxa = 100
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 4
  variance <- matrix(0.8, p, p) + diag(0.2, p, p)
  independent <- FALSE
  root.state <- list(random = FALSE,
                     value.root = rep(1, p),
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(12, 24, 30),
                values=matrix(2*c(1, -0.5), nrow = p, ncol = 3),
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
  nMiss <- floor(ntaxa * p * 0.5)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  # traits[2, 3] <- traits[2, 17] <- traits[2, 18] <- traits[, 6] <- NA
  
  ## Log lik old way
  miss <- as.vector(is.na(traits))
  masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
  masque_data[1:(p * ntaxa)] <- !miss
  
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  T_tree <- incidence.matrix(tree)
  U_tree <- incidence.matrix.full(tree)
  h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  root_edge_length <- 0
  if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
  F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)
  
  tmp_old <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "BM",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            F_moments = F_moments,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.simple)
  
  ll_old <- tmp_old$log_likelihood_old
  conditional_law_X_old <- tmp_old$conditional_law_X
  
  tmp_new <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "BM",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            F_moments = F_moments,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.upward_downward)
  
  ll_new <- tmp_new$log_likelihood_old
  conditional_law_X_new <- tmp_new$conditional_law_X
  
  cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
  cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))

  
  expect_that(ll_new, equals(as.vector(ll_old)))
  expect_that(conditional_law_X_new$expectations, equals(conditional_law_X_old$expectations))
  expect_that(conditional_law_X_new$variances, equals(conditional_law_X_old$variances))
  expect_that(cov_old, equals(cov_new))
})


test_that("Upward Downward - BM - no missing", {
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920902)
  ntaxa = 500
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 10
  variance <- matrix(0.8, p, p) + diag(0.2, p, p)
  independent <- FALSE
  root.state <- list(random = FALSE,
                     value.root = rep(1, p),
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(18, 32, 45, 109, 254, 398),
                values=cbind(c(4, -10, 3, 12, 32, -5),
                             c(-5, 5, 0, -1.1, 32.89, 16),
                             c(4, -10, 3, 12, 32, -5),
                             c(4, -10, 3, 12, 32, -5),
                             c(4, -10, 3, 12, 32, -5),
                             c(4, -10, 3, 12, 32, -5)),
                relativeTimes = 0)
  shifts$values <- apply(shifts$values, 1, function(z) rep(z, length.out = p))
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
  
  ## Log lik old way
  miss <- as.vector(is.na(traits))
  masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
  masque_data[1:(p * ntaxa)] <- !miss
  
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  T_tree <- incidence.matrix(tree)
  U_tree <- incidence.matrix.full(tree)
  h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  root_edge_length <- 0
  if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
  F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)
  
  tmp_old <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "BM",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            F_moments = F_moments,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.simple.nomissing.BM)
  
  ll_old <- tmp_old$log_likelihood_old
  conditional_law_X_old <- tmp_old$conditional_law_X
  
  tmp_new <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "BM",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            F_moments = F_moments,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.upward_downward)
  
  ll_new <- tmp_new$log_likelihood_old
  conditional_law_X_new <- tmp_new$conditional_law_X
  
  cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
  cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))
  
  
  expect_that(ll_new, equals(as.vector(ll_old)))
  expect_that(conditional_law_X_new$expectations, equals(conditional_law_X_old$expectations))
  expect_that(conditional_law_X_new$variances, equals(conditional_law_X_old$variances))
  expect_that(cov_old, equals(cov_new))
  
  # library(microbenchmark)
  # microbenchmark(wrapper_E_step(phylo = tree,
  #                               times_shared = times_shared,
  #                               distances_phylo = distances_phylo,
  #                               process = "BM",
  #                               params_old = paramsSimu,
  #                               masque_data = masque_data,
  #                               F_moments = F_moments,
  #                               independent = independent,
  #                               Y_data_vec_known = as.vector(traits[!miss]),
  #                               miss = miss,
  #                               Y_data = traits,
  #                               U_tree = U_tree,
  #                               compute_E = compute_E.simple.nomissing.BM),
  #                wrapper_E_step(phylo = tree,
  #                               process = "BM",
  #                               params_old = paramsSimu,
  #                               independent = independent,
  #                               Y_data = traits,
  #                               compute_E = compute_E.upward_downward),
  #                times = 100)
})

test_that("Upward Downward - estimateEM - BM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920902)
  ntaxa = 100
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 4
  variance <- matrix(0.8, p, p) + diag(0.2, p, p)
  independent <- FALSE
  root.state <- list(random = TRUE,
                     value.root = NA,
                     exp.root = rep(1, p),
                     var.root = variance / 2)
  shifts = list(edges = c(12, 24, 30),
                values=matrix(4*c(1, -0.5), nrow = p, ncol = 3),
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
  nMiss <- floor(ntaxa * p * 0.1)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  
  expect_warning(res_old <- estimateEM(phylo = tree, 
                                       Y_data= traits, 
                                       process = "BM", 
                                       method.variance = "simple", 
                                       method.init = "lasso",
                                       Nbr_It_Max = 2,
                                       nbr_of_shifts = 3,
                                       random.root = FALSE,
                                       convergence_mode = "relative",
                                       K_lag_init = 2),
                 "The maximum number of iterations")
  
  expect_warning(res_new <- estimateEM(phylo = tree, 
                                       Y_data= traits, 
                                       process = "BM", 
                                       method.variance = "upward_downward", 
                                       method.init = "lasso",
                                       Nbr_It_Max = 2,
                                       nbr_of_shifts = 3,
                                       random.root = FALSE,
                                       convergence_mode = "relative",
                                       K_lag_init = 2),
                 "The maximum number of iterations")
  
  expect_that(res_new, equals(as.vector(res_old)))
})

test_that("Upward Downward - PhyloEM - BM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920902)
  ntaxa = 50
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 2
  variance <- matrix(0.8, p, p) + diag(0.2, p, p)
  independent <- FALSE
  root.state <- list(random = FALSE,
                     value.root = rep(1, p),
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(25, 31, 59),
                values = cbind(c(5, -5),
                               c(-5, -5),
                               c(3, 3)),
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
  nMiss <- floor(ntaxa * p * 0.1)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  
  expect_warning(res_old <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "BM",
                                    independent = FALSE,
                                    K_max = 10,
                                    Nbr_It_Max = 2,
                                    use_previous = FALSE,
                                    order = TRUE,
                                    method.selection = "BirgeMassart1",
                                    method.variance = "simple",
                                    method.init = "lasso",
                                    method.init.alpha = "estimation",
                                    method.init.alpha.estimation = c("regression", "regression.MM", "median"), 
                                    methods.segmentation = c("lasso", "best_single_move"),
                                    random.root = FALSE,
                                    progress.bar = FALSE,
                                    sBM_variance = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")
  
  expect_warning(res_new <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "BM",
                                    independent = FALSE,
                                    K_max = 10,
                                    Nbr_It_Max = 2,
                                    use_previous = FALSE,
                                    order = TRUE,
                                    method.selection = "BirgeMassart1",
                                    method.variance = "upward_downward",
                                    method.init = "lasso",
                                    method.init.alpha = "estimation",
                                    method.init.alpha.estimation = c("regression", "regression.MM", "median"), 
                                    methods.segmentation = c("lasso", "best_single_move"),
                                    random.root = FALSE,
                                    progress.bar = FALSE,
                                    sBM_variance = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")
  
  ## Time is different
  res_new$alpha_0$results_summary$time <- 0
  res_old$alpha_0$results_summary$time <- 0
  res_new$alpha_max$results_summary$time <- 0
  res_old$alpha_max$results_summary$time <- 0
  res_new$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_old$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_new$alpha_max$Djump_BM1$results_summary$time <- 0
  res_old$alpha_max$Djump_BM1$results_summary$time <- 0
  res_new$alpha_min$results_summary$time <- 0
  res_old$alpha_min$results_summary$time <- 0
  res_new$alpha_min_raw$results_summary$time <- 0
  res_old$alpha_min_raw$results_summary$time <- 0
  
  ## Sames objects
  expect_equal(res_new, res_old, tolerance = .Machine$double.eps ^ 0.3)
  
  ## Log Likelihood
  expect_equal(log_likelihood(res_new),
               res_new$alpha_max$DDSE_BM1$results_summary$log_likelihood)
  ll <- log_likelihood(res_new, K = 5, alpha = 0)
  expect_equal(ll, res_new$alpha_0$results_summary$log_likelihood[6])
})


test_that("Upward Downward - scOU - fixed root", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920902)
  ntaxa = 100
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 4
  variance <- diag(0.2, p, p) +  matrix(0.8, p, p)
  selection.strength <- 3
  independent <- FALSE
  root.state <- list(random = FALSE,
                     stationary.root = FALSE,
                     value.root = rep(1, p),
                     exp.root = NA,
                     var.root = NA) # compute_stationary_variance(variance, selection.strength))
  shifts <- list(edges = c(12, 26, 165),
                 values = cbind(c(4, -10, 3, 1),
                                c(-5, 5, 0, 1),
                                c(4, -10, 3, -1)),
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
  nMiss <- floor(ntaxa * p * 0.5)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  # traits[2, 3] <- traits[2, 17] <- traits[2, 18] <- traits[, 6] <- NA

  ## Log lik old way
  miss <- as.vector(is.na(traits))
  masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
  masque_data[1:(p * ntaxa)] <- !miss

  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  T_tree <- incidence.matrix(tree)
  U_tree <- incidence.matrix.full(tree)
  h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])

  tmp_old <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "scOU",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.simple)

  ll_old <- tmp_old$log_likelihood_old
  conditional_law_X_old <- tmp_old$conditional_law_X

  tmp_new <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "scOU",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.upward_downward)

  ll_new <- tmp_new$log_likelihood_old
  conditional_law_X_new <- tmp_new$conditional_law_X

  cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
  cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))


  expect_that(ll_new, equals(as.vector(ll_old)))
  expect_that(conditional_law_X_new$expectations, equals(conditional_law_X_old$expectations))
  expect_that(conditional_law_X_new$optimal.values, equals(conditional_law_X_old$optimal.values))
  expect_that(conditional_law_X_new$variances, equals(conditional_law_X_old$variances))
  expect_that(cov_old, equals(cov_new))
})

test_that("Upward Downward - scOU - random root", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(17920902)
  ntaxa = 100
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 4
  variance <- diag(0.2, p, p) +  matrix(0.8, p, p)
  selection.strength <- 3
  independent <- FALSE
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = rep(1, p),
                     var.root = compute_stationary_variance(variance, selection.strength))
  shifts <- list(edges = c(12, 26, 165),
                 values = cbind(c(4, -10, 3, 1),
                                c(-5, 5, 0, 1),
                                c(4, -10, 3, -1)),
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
  nMiss <- floor(ntaxa * p * 0.5)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  # traits[2, 3] <- traits[2, 17] <- traits[2, 18] <- traits[, 6] <- NA
  
  ## Log lik old way
  miss <- as.vector(is.na(traits))
  masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
  masque_data[1:(p * ntaxa)] <- !miss
  
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  T_tree <- incidence.matrix(tree)
  U_tree <- incidence.matrix.full(tree)
  h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  
  tmp_old <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "scOU",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.simple)
  
  ll_old <- tmp_old$log_likelihood_old
  conditional_law_X_old <- tmp_old$conditional_law_X
  
  tmp_new <- wrapper_E_step(phylo = tree,
                            times_shared = times_shared,
                            distances_phylo = distances_phylo,
                            process = "scOU",
                            params_old = paramsSimu,
                            masque_data = masque_data,
                            independent = independent,
                            Y_data_vec_known = as.vector(traits[!miss]),
                            miss = miss,
                            Y_data = traits,
                            U_tree = U_tree,
                            compute_E = compute_E.upward_downward)
  
  ll_new <- tmp_new$log_likelihood_old
  conditional_law_X_new <- tmp_new$conditional_law_X
  
  cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
  cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))
  
  
  expect_that(ll_new, equals(as.vector(ll_old)))
  expect_that(conditional_law_X_new$expectations, equals(conditional_law_X_old$expectations))
  expect_that(conditional_law_X_new$optimal.values, equals(conditional_law_X_old$optimal.values))
  expect_that(conditional_law_X_new$variances, equals(conditional_law_X_old$variances))
  expect_that(cov_old, equals(cov_new))
})

test_that("Upward Downward - PhyloEM - scOU - fixed root", {
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
  nMiss <- floor(ntaxa * p * 0.1)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }

  expect_warning(res_old <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "scOU",
                                    K_max = 10,
                                    random.root = FALSE,
                                    stationary.root = FALSE,
                                    alpha = selection.strength,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "simple",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")
  
  expect_warning(res_new <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "scOU",
                                    K_max = 10,
                                    random.root = FALSE,
                                    stationary.root = FALSE,
                                    alpha = selection.strength,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "upward_downward",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")

  ## Time is different
  res_new$alpha_3$results_summary$time <- 0
  res_old$alpha_3$results_summary$time <- 0
  res_new$alpha_max$results_summary$time <- 0
  res_old$alpha_max$results_summary$time <- 0
  res_new$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_old$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_new$alpha_max$Djump_BM1$results_summary$time <- 0
  res_old$alpha_max$Djump_BM1$results_summary$time <- 0
  res_new$alpha_min$results_summary$time <- 0
  res_old$alpha_min$results_summary$time <- 0
  res_new$alpha_min_raw$results_summary$time <- 0
  res_old$alpha_min_raw$results_summary$time <- 0

  expect_that(res_new, equals(res_old))
  
  expect_equal(log_likelihood(res_new),
               res_new$alpha_max$DDSE_BM1$results_summary$log_likelihood)
  expect_equal(log_likelihood(res_new, K = 5, alpha = "max"),
               res_new$alpha_max$results_summary$log_likelihood[6])
})

test_that("Upward Downward - PhyloEM - scOU - random root", {
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
  nMiss <- floor(ntaxa * p * 0.1)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  
  expect_warning(res_old <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "scOU",
                                    K_max = 10,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    alpha = selection.strength,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "simple",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")
  
  expect_warning(res_new <- PhyloEM(phylo = tree,
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
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")
  
  ## Time is different
  res_new$alpha_3$results_summary$time <- 0
  res_old$alpha_3$results_summary$time <- 0
  res_new$alpha_max$results_summary$time <- 0
  res_old$alpha_max$results_summary$time <- 0
  res_new$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_old$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_new$alpha_max$Djump_BM1$results_summary$time <- 0
  res_old$alpha_max$Djump_BM1$results_summary$time <- 0
  res_new$alpha_min$results_summary$time <- 0
  res_old$alpha_min$results_summary$time <- 0
  res_new$alpha_min_raw$results_summary$time <- 0
  res_old$alpha_min_raw$results_summary$time <- 0
  
  expect_that(res_new, equals(res_old))
  
  expect_equal(log_likelihood(res_new),
               res_new$alpha_max$DDSE_BM1$results_summary$log_likelihood)
  expect_equal(log_likelihood(res_new, K = 5, alpha = "max"),
               res_new$alpha_max$results_summary$log_likelihood[6])
})

test_that("Upward Downward - PhyloEM - scOU - random root - un-ordered", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  ntaxa = 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                          age = 1, mrca = TRUE)[[1]]
  # tree <- reorder(tree, order = "postorder")
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
  # nMiss <- floor(ntaxa * p * 0.1)
  # miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  # chars <- (miss - 1) %% p + 1
  # tips <- (miss - 1) %/% p + 1
  # for (i in 1:nMiss){
  #   traits[chars[i], tips[i]] <- NA
  # }
  
  expect_warning(res_old <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "scOU",
                                    K_max = 10,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    alpha = c(1.33, selection.strength),
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "simple",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")
  
  expect_warning(res_new <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "scOU",
                                    K_max = 10,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    alpha = c(1.33, selection.strength),
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "upward_downward",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 2),
                 "The maximum number of iterations")
  
  ## Time is different
  res_new$alpha_1.33$results_summary$time <- 0
  res_old$alpha_1.33$results_summary$time <- 0
  res_new$alpha_3$results_summary$time <- 0
  res_old$alpha_3$results_summary$time <- 0
  res_new$alpha_max$results_summary$time <- 0
  res_old$alpha_max$results_summary$time <- 0
  res_new$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_old$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_new$alpha_max$Djump_BM1$results_summary$time <- 0
  res_old$alpha_max$Djump_BM1$results_summary$time <- 0
  res_new$alpha_min$results_summary$time <- 0
  res_old$alpha_min$results_summary$time <- 0
  res_new$alpha_min_raw$results_summary$time <- 0
  res_old$alpha_min_raw$results_summary$time <- 0
  ## Init is not in the same order
  res_new$alpha_1.33$params_init_estim <- NULL
  res_old$alpha_1.33$params_init_estim <- NULL
  res_new$alpha_3$params_init_estim <- NULL
  res_old$alpha_3$params_init_estim <- NULL
  res_new$alpha_max$params_init_estim <- NULL
  res_old$alpha_max$params_init_estim <- NULL
  res_new$alpha_min$params_init_estim <- NULL
  res_old$alpha_min$params_init_estim <- NULL
  res_new$alpha_min_raw$params_init_estim <- NULL
  res_old$alpha_min_raw$params_init_estim <- NULL
  res_new$alpha_max$DDSE_BM1$params_init_estim <- NULL
  res_old$alpha_max$DDSE_BM1$params_init_estim <- NULL
  res_new$alpha_max$Djump_BM1$params_init_estim <- NULL
  res_old$alpha_max$Djump_BM1$params_init_estim <- NULL
  
  expect_equal(res_new, res_old, check.attributes = FALSE)
  
  expect_equal(log_likelihood(res_new),
               res_new$alpha_max$DDSE_BM1$results_summary$log_likelihood)
  expect_equal(log_likelihood(res_new, K = 5, alpha = "max"),
               res_new$alpha_max$results_summary$log_likelihood[6])
})

test_that("Upward Downward - PhyloEM - OU - independent", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  ntaxa = 50
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 2
  variance <- diag(0.2, p, p)
  selection.strength <- diag(3, p, p)
  independent <- FALSE
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = rep(1, p),
                     var.root = compute_stationary_variance(variance, selection.strength))
  shifts = list(edges = c(25, 10, 59),
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
                          process = "OU",
                          variance = variance,
                          shifts = shifts,
                          selection.strength = selection.strength,
                          optimal.value = paramsSimu$optimal.value)
  
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  nMiss <- floor(ntaxa * p * 0.1)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  
  expect_warning(res_old <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "OU",
                                    independent = TRUE,
                                    alpha_grid = FALSE,
                                    K_max = 10,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "simple",
                                    method.init.alpha.estimation = c("regression.MM", "median"),
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 1))
  
  expect_warning(res_new <- PhyloEM(phylo = tree,
                                    Y_data = traits,
                                    process = "OU",
                                    independent = TRUE,
                                    alpha_grid = FALSE,
                                    K_max = 10,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "upward_downward",
                                    method.init.alpha.estimation = c("regression.MM", "median"),
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "BirgeMassart1",
                                    progress.bar = FALSE,
                                    K_lag_init = 1))
  
  ## Time is different
  res_new$alpha_max$results_summary$time <- 0
  res_old$alpha_max$results_summary$time <- 0
  res_new$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_old$alpha_max$DDSE_BM1$results_summary$time <- 0
  res_new$alpha_max$Djump_BM1$results_summary$time <- 0
  res_old$alpha_max$Djump_BM1$results_summary$time <- 0
  
  expect_that(res_new, equals(res_old))
  
  ll <- suppressWarnings(log_likelihood(res_new))
  
  expect_equal(ll,
               res_new$alpha_max$DDSE_BM1$results_summary$log_likelihood)
  expect_warning(ll <- log_likelihood(res_new, K = 5, alpha = "max"),
                 "There are several equivalent solutions for this shift position.")

  expect_equal(ll, res_new$alpha_max$results_summary$log_likelihood[6])
})

# test_that("Upward Downward - full OU - random root", {
#   set.seed(17920902)
#   ntaxa = 100
#   tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
#                           age = 1, mrca = TRUE)[[1]]
#   tree <- reorder(tree, order = "postorder")
#   p <- 4
#   variance <- diag(0.2, p, p) +  matrix(0.8, p, p)
#   selection.strength <- diag(0.2, p, p) +  matrix(2.8, p, p)
#   independent <- FALSE
#   root.state <- list(random = TRUE,
#                      stationary.root = TRUE,
#                      value.root = NA,
#                      exp.root = rep(1, p),
#                      var.root = compute_stationary_variance(variance, selection.strength))
#   shifts <- list(edges = c(12, 26, 165),
#                  values = cbind(c(4, -10, 3, 1),
#                                 c(-5, 5, 0, 1),
#                                 c(4, -10, 3, -1)),
#                  relativeTimes = 0)
#   paramsSimu <- list(variance = variance,
#                      shifts = shifts,
#                      root.state = root.state,
#                      selection.strength = selection.strength,
#                      optimal.value = rep(1, p))
#   attr(paramsSimu, "p_dim") <- p
#   
#   X1 <- simulate_internal(tree,
#                  p = p,
#                  root.state = root.state,
#                  process = "OU",
#                  variance = variance,
#                  shifts = shifts,
#                  selection.strength = selection.strength,
#                  optimal.value = paramsSimu$optimal.value)
#   
#   traits <- extract_simulate_internal(X1, where = "tips", what = "state")
#   # nMiss <- floor(ntaxa * p * 0.5)
#   # miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
#   # chars <- (miss - 1) %% p + 1
#   # tips <- (miss - 1) %/% p + 1
#   # for (i in 1:nMiss){
#   #   traits[chars[i], tips[i]] <- NA
#   # }
#   # traits[2, 3] <- traits[2, 17] <- traits[2, 18] <- traits[, 6] <- NA
#   
#   ## Log lik old way
#   miss <- as.vector(is.na(traits))
#   masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
#   masque_data[1:(p * ntaxa)] <- !miss
#   
#   times_shared <- compute_times_ca(tree)
#   distances_phylo <- compute_dist_phy(tree)
#   T_tree <- incidence.matrix(tree)
#   U_tree <- incidence.matrix.full(tree)
#   h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
#   root_edge_length <- 0
#   if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
#   F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)
#   
#   tmp_old <- wrapper_E_step(phylo = tree,
#                             times_shared = times_shared,
#                             distances_phylo = distances_phylo,
#                             process = "OU",
#                             params_old = paramsSimu,
#                             masque_data = masque_data,
#                             F_moments = F_moments,
#                             independent = independent,
#                             Y_data_vec_known = as.vector(traits[!miss]),
#                             miss = miss,
#                             Y_data = traits,
#                             U_tree = U_tree,
#                             compute_E = compute_E.simple)
#   
#   ll_old <- tmp_old$log_likelihood_old
#   conditional_law_X_old <- tmp_old$conditional_law_X
#   
#   tmp_new <- wrapper_E_step(phylo = tree,
#                             times_shared = times_shared,
#                             distances_phylo = distances_phylo,
#                             process = "OU",
#                             params_old = paramsSimu,
#                             masque_data = masque_data,
#                             F_moments = F_moments,
#                             independent = independent,
#                             Y_data_vec_known = as.vector(traits[!miss]),
#                             miss = miss,
#                             Y_data = traits,
#                             U_tree = U_tree,
#                             compute_E = compute_E.upward_downward)
#   
#   ll_new <- tmp_new$log_likelihood_old
#   conditional_law_X_new <- tmp_new$conditional_law_X
#   
#   cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
#   cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))
#   
#   
#   expect_that(ll_new, equals(as.vector(ll_old)))
#   expect_that(conditional_law_X_new$expectations, equals(conditional_law_X_old$expectations))
#   expect_that(conditional_law_X_new$variances, equals(conditional_law_X_old$variances))
#   expect_that(cov_old, equals(cov_new))
# })

# test_that("Upward Downward - estimateEM - OU re-scaled", {
#   set.seed(17920902)
#   ntaxa = 100
#   tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
#                           age = 1, mrca = TRUE)[[1]]
#   tree <- reorder(tree, order = "postorder")
#   p <- 4
#   variance <- matrix(0.8, p, p) + diag(0.2, p, p)
#   selection.strength <- 3
#   independent <- FALSE
#   root.state <- list(random = TRUE,
#                      stationary.root = TRUE,
#                      value.root = NA,
#                      exp.root = rep(1, p),
#                      var.root = compute_stationary_variance(variance, selection.strength))
#   shifts = list(edges = c(12, 24, 30),
#                 values=matrix(4*c(1, -0.5), nrow = p, ncol = 3),
#                 relativeTimes = 0)
#   paramsSimu <- list(variance = variance,
#                      shifts = shifts,
#                      root.state = root.state,
#                      selection.strength = selection.strength,
#                      optimal.value = root.state$exp.root)
#   attr(paramsSimu, "p_dim") <- p
#   
#   X1 <- simulate_internal(tree,
#                  p = p,
#                  root.state = root.state,
#                  process = "scOU",
#                  variance = variance,
#                  shifts = shifts,
#                  selection.strength = selection.strength,
#                  optimal.value = root.state$exp.root)
#   
#   traits <- extract_simulate_internal(X1, where = "tips", what = "state")
#   nMiss <- floor(ntaxa * p * 0.1)
#   miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
#   chars <- (miss - 1) %% p + 1
#   tips <- (miss - 1) %/% p + 1
#   for (i in 1:nMiss){
#     traits[chars[i], tips[i]] <- NA
#   }
#   
#   expect_warning(res_old <- estimateEM(phylo = tree, 
#                                        Y_data = traits, 
#                                        process = "scOU", 
#                                        Nbr_It_Max = 2, 
#                                        method.variance = "simple", 
#                                        method.init = "lasso",
#                                        method.init.alpha = "estimation",
#                                        method.init.alpha.estimation = c("regression", 
#                                                                         "regression.MM", 
#                                                                         "median"),
#                                        nbr_of_shifts = 3,
#                                        random.root = TRUE,
#                                        stationary.root = TRUE,
#                                        alpha_known = TRUE,
#                                        known.selection.strength = selection.strength,
#                                        methods.segmentation = c("max_costs_0", 
#                                                                 "lasso", 
#                                                                 "same_shifts", 
#                                                                 "same_shifts_same_values",
#                                                                 "best_single_move", 
#                                                                 "lasso_one_move"),
#                                        sBM_variance = FALSE,
#                                        method.OUsun = "rescale",
#                                        K_lag_init = 2),
#                  "The maximum number of iterations")
#   
#   expect_warning(res_new <- estimateEM(phylo = tree, 
#                                        Y_data = traits, 
#                                        process = "scOU", 
#                                        Nbr_It_Max = 2, 
#                                        method.variance = "upward_downward", 
#                                        method.init = "lasso",
#                                        method.init.alpha = "estimation",
#                                        method.init.alpha.estimation = c("regression", 
#                                                                         "regression.MM", 
#                                                                         "median"),
#                                        nbr_of_shifts = 3,
#                                        random.root = TRUE,
#                                        stationary.root = TRUE,
#                                        alpha_known = TRUE,
#                                        known.selection.strength = selection.strength,
#                                        methods.segmentation = c("max_costs_0", 
#                                                                 "lasso", 
#                                                                 "same_shifts", 
#                                                                 "same_shifts_same_values",
#                                                                 "best_single_move", 
#                                                                 "lasso_one_move"),
#                                        sBM_variance = FALSE,
#                                        method.OUsun = "rescale",
#                                        K_lag_init = 2),
#                  "The maximum number of iterations")
#   
#   expect_that(res_new, equals(as.vector(res_old)))
# })