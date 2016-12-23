context("E step no missing data")

test_that("Likelihood missing/no missing methods", {
  testthat::skip_on_cran()
  set.seed(586)
  ntaxa <- 200
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
                                   lambda = 0.1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  
  ## Parameters
  p <- 6
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  root.state <- list(random = FALSE,
                     value.root = c(1, -1, 2, -10, 2.58, -13.2),
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
  
  params = list(variance = variance,
                root.state = root.state,
                shifts = shifts)
  attr(params, "p_dim") <- p
  
  ## Fixed Quantities 
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  subtree.list <- enumerate_tips_under_edges(tree)
  T_tree <- incidence.matrix(tree)
  h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  F_moments <- compute_fixed_moments(times_shared, ntaxa)
  
  ## Data
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "BM",
                          variance = variance,
                          shifts = shifts)
  Y_data <- extract_simulate_internal(X1,"tips","states")
  
  ## Old Method
  moments_old <- compute_mean_variance.simple(phylo = tree,
                                              times_shared = times_shared,
                                              distances_phylo = distances_phylo,
                                              process = "BM",
                                              params_old = params)
  
  log_likelihood_old <- compute_log_likelihood.simple(phylo = tree,
                                                      Y_data_vec = as.vector(Y_data),
                                                      sim = moments_old$sim,
                                                      Sigma = moments_old$Sigma,
                                                      Sigma_YY_chol_inv = moments_old$Sigma_YY_chol_inv)
  
  maha_data_mean_old <- compute_mahalanobis_distance.simple(phylo = tree,
                                                            Y_data_vec = as.vector(Y_data),
                                                            sim = moments_old$sim,
                                                            Sigma_YY_chol_inv = moments_old$Sigma_YY_chol_inv)
  
  conditional_law_X_old <- compute_cond_law.simple(phylo = tree,
                                                   Y_data_vec = as.vector(Y_data),
                                                   sim = moments_old$sim,
                                                   Sigma = moments_old$Sigma,
                                                   Sigma_YY_chol_inv = moments_old$Sigma_YY_chol_inv)
  
  ## New Method
  moments_new <- compute_mean_variance.simple.nomissing.BM(phylo = tree,
                                                           params_old = params,
                                                           F_moments = F_moments)
  
  log_likelihood_new <- compute_log_likelihood.simple.nomissing.BM(
    phylo = tree,
    Y_data = Y_data,
    sim = moments_new$sim,
    C_YY = F_moments$C_YY,
    C_YY_chol_inv = F_moments$C_YY_chol_inv,
    R = params$variance)
  
  maha_data_mean_new <- compute_mahalanobis_distance.simple.nomissing.BM(
    phylo = tree,
    sim = moments_new$sim,
    Y_data = Y_data,
    C_YY_chol_inv = F_moments$C_YY_chol_inv,
    R = params$variance)
  
  conditional_law_X_new <- compute_cond_law.simple.nomissing.BM(phylo = tree,
                                                                sim = moments_new$sim,
                                                                F_means = F_moments$F_means,
                                                                F_vars = F_moments$F_vars,
                                                                R = params$variance,
                                                                Y_data = Y_data)
  
  expect_that(as.vector(log_likelihood_new), equals(as.vector(log_likelihood_old)))
  expect_that(as.vector(maha_data_mean_new), equals(as.vector(maha_data_mean_old)))
  expect_that(conditional_law_X_new$expectations, equals(conditional_law_X_old$expectations))
  expect_that(conditional_law_X_new$variances, equals(conditional_law_X_old$variances))
})