# context("Log Likelihood Computation")
# 
# test_that("independent computations", {
#   set.seed(586)
#   ntaxa <- 50
#   tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
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
#   X1 <- simulate(tree,
#                  p = p,
#                  root.state = root.state,
#                  process = "OU",
#                  variance = variance,
#                  optimal.value = optimal.value,
#                  selection.strength = selection.strength,
#                  shifts = shifts)
#   Y_sim <- extract.simulate(X1, "tips", "states")
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