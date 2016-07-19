context("Functions for multivariate independent traits")

test_that("split/merge independent parameters ", {
  # Dimension p - BM - non random root
  p <- 4
  ntaxa <- 236
  
  params <- init.EM.default.BM(p = p,
                               variance.init = diag(1:p, p, p),
                               random.init = FALSE,
                               value.root.init = 1:p,
                               exp.root.init = 1:p,
                               var.root.init = 1:p,
                               edges.init = c(10, 23, 27),
                               values.init = matrix(1:(p*3), p, 3))
  
  params_split <- split_params_independent(params)
  params_bis <- merge_params_independent(params_split)
  
  expect_that(params, equals(params_bis))
  
  # Dimension p - BM - random root
  p <- 4
  ntaxa <- 236
  
  params <- init.EM.default.BM(p = p,
                               variance.init = diag(1:p),
                               random.init = TRUE,
                               value.root.init = 1:p,
                               exp.root.init = 1:p,
                               var.root.init = diag(1:p),
                               edges.init = c(10, 23, 27),
                               values.init = matrix(1:(p*3), p, 3))
  
  params_split <- split_params_independent(params)
  params_bis <- merge_params_independent(params_split)
  
  expect_that(params, equals(params_bis))
  
  # Dimension p - OU- stationary root
  p <- 4
  ntaxa <- 236
  
  params <- init.EM.default.OU(p = p,
                               variance.init = diag(1:p),
                               random.init = TRUE,
                               stationary.root.init = TRUE,
                               value.root.init = 1:p,
                               exp.root.init = 1:p,
                               var.root.init = diag(1:p),
                               edges.init = c(11, 14, 23),
                               values.init = matrix(1:(p*3), p, 3),
                               selection.strength.init = diag(1:p),
                               optimal.value.init = 1:p)
  
  params_split <- split_params_independent(params)
  params_bis <- merge_params_independent(params_split)
  
  expect_that(params, equals(params_bis))
  
  # Dimension p - OU - random
  p <- 4
  ntaxa <- 236
  
  params <- init.EM.default.OU(p = p,
                               variance.init = diag(1:p),
                               random.init = TRUE,
                               stationary.root.init = FALSE,
                               value.root.init = 1:p,
                               exp.root.init = 1:p,
                               var.root.init = diag(1:p),
                               edges.init = c(11, 14, 23),
                               values.init = matrix(1:(p*3), p, 3),
                               selection.strength.init = diag(1:p),
                               optimal.value.init = 1:p)
  
  params_split <- split_params_independent(params)
  params_bis <- merge_params_independent(params_split)
  
  expect_that(params, equals(params_bis))
  
  # Dimension p - OU - non random
  p <- 4
  ntaxa <- 236
  
  expect_warning(params <- init.EM.default.OU(p = p,
                                              variance.init = diag(1:p),
                                              random.init = FALSE,
                                              stationary.root.init = FALSE,
                                              value.root.init = 1:p,
                                              exp.root.init = 1:p,
                                              var.root.init = diag(1:p),
                                              edges.init = c(11, 14, 23),
                                              values.init = matrix(1:(p*3), p, 3),
                                              selection.strength.init = diag(1:p),
                                              optimal.value.init = 1:p))
  
  params_split <- split_params_independent(params)
  params_bis <- merge_params_independent(params_split)
  
  expect_that(params, equals(params_bis))
  
  # Dimension p - BM - no shift
  p <- 4
  ntaxa <- 236
  
  params <- init.EM.default.BM(p = p,
                               variance.init = diag(1:p, p, p),
                               random.init = FALSE,
                               value.root.init = 1:p,
                               exp.root.init = 1:p,
                               var.root.init = 1:p,
                               edges.init = NULL,
                               values.init = NULL)
  
  params_split <- split_params_independent(params)
  params_bis <- merge_params_independent(params_split)
  
  expect_that(params, equals(params_bis))
  
  # Dimension p - OU - one shift
  p <- 4
  ntaxa <- 236
  
  params <- init.EM.default.OU(p = p,
                               variance.init = diag(1:p),
                               random.init = FALSE,
                               stationary.root.init = FALSE,
                               value.root.init = 1:p,
                               exp.root.init = NA,
                               var.root.init = NA,
                               edges.init = c(11),
                               values.init = matrix(1:(p*1), p, 1),
                               selection.strength.init = diag(1:p),
                               optimal.value.init = 1:p)
  
  params_split <- split_params_independent(params)
  params_bis <- merge_params_independent(params_split)
  
  expect_that(params, equals(params_bis))
})

test_that("compute_mean_variance.simple", {
  # Dimension p - OU- stationary root
  set.seed(586)
  ntaxa <- 20
  tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
                          lambda = 0.1, mu = 0,
                          age = 1, mrca = TRUE)[[1]]
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  p <- 4
  
  Y_data <- matrix(1:(ntaxa * p), ncol = ntaxa)
  Y_data[1, 1] <- NA; Y_data[2, 4] <- NA; Y_data[4, 15] <- NA;
  Y_data_vec <- as.vector(Y_data)
  miss <- as.vector(is.na(Y_data))
  Y_data_vec_known <- as.vector(Y_data[!miss])
  # Vectorized Data Mask
  masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
  masque_data[1:(p*ntaxa)] <- !miss
  
  params <- init.EM.default.OU(p = p,
                               variance.init = diag(1:p),
                               random.init = TRUE,
                               stationary.root.init = TRUE,
                               value.root.init = 1:p,
                               exp.root.init = 1:p,
                               var.root.init = diag(1:p),
                               edges.init = c(11, 14, 23),
                               values.init = matrix(1:(p*3), p, 3),
                               selection.strength.init = diag(1:p),
                               optimal.value.init = 1:p)
  attr(params, "p_dim") <- p
  
  params_split <- split_params_independent(params)
  
  res_1 <- wrapper_E_step(phylo = tree,
                          times_shared = times_shared,
                          distances_phylo = distances_phylo,
                          process = "OU",
                          params_old = params,
                          masque_data = masque_data, 
                          independent = FALSE,
                          F_moments = NULL,
                          Y_data_vec_known = Y_data_vec_known,
                          miss = miss,
                          Y_data = Y_data,
                          compute_mean_variance = compute_mean_variance.simple,
                          compute_log_likelihood = compute_log_likelihood.simple,
                          compute_mahalanobis_distance = compute_mahalanobis_distance.simple,
                          compute_E = compute_E.simple)
  
  res_2 <- wrapper_E_step(phylo = tree,
                          times_shared = times_shared,
                          distances_phylo = distances_phylo,
                          process = "OU",
                          params_old = params_split,
                          masque_data = masque_data, 
                          independent = TRUE,
                          F_moments = NULL,
                          Y_data_vec_known = Y_data_vec_known,
                          miss = miss,
                          Y_data = Y_data,
                          compute_mean_variance = compute_mean_variance.simple,
                          compute_log_likelihood = compute_log_likelihood.simple,
                          compute_mahalanobis_distance = compute_mahalanobis_distance.simple,
                          compute_E = compute_E.simple)
  
  expect_equal(as.vector(res_1$log_likelihood_old),
              sum(sapply(res_2, function(z) return(z$log_likelihood_old))))
  expect_equal(as.vector(res_1$maha_data_mean),
              sum(sapply(res_2, function(z) return(z$maha_data_mean))))
})