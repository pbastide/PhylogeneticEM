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
  
  params <- init.EM.default.OU(p = p,
                               variance.init = diag(1:p),
                               random.init = FALSE,
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
})

compute_mean_variance <- compute_mean_variance.simple
compute_log_likelihood <- compute_log_likelihood.simple
compute_mahalanobis_distance <- compute_mahalanobis_distance.simple
compute_E <- compute_E.simple

wrapper_E_step <- function(phylo,
                           times_shared,
                           distances_phylo,
                           process,
                           params_old,
                           masque_data,
                           F_moments,
                           independent,
                           Y_data_vec_known,
                           miss,
                           Y_data){
  if (independent){
    # if independent, params_old is a list of p params
    masque_data_matr <- matrix(masque_data,
                               ncol = length(phylo$tip.label) + phylo$Nnode)
    miss_matr <- matrix(miss,
                        ncol = length(phylo$tip.label))
    res <- vector(mode = "list", length = length(params_old))
    for (i in 1:length(res)){
      res[[i]] <- wrapper_E_step(phylo = phylo,
                                 times_shared = times_shared,
                                 distances_phylo = distances_phylo,
                                 process = process,
                                 params_old = params_old[[i]],
                                 masque_data = masque_data_matr[i, ],
                                 F_moments = F_moments,
                                 independent = FALSE,
                                 Y_data_vec_known = Y_data[i, !miss_matr[i, ]],
                                 miss = miss_matr[i, ],
                                 Y_data = Y_data[i, , drop = F])
    }
    return(res)
  }
  moments <- compute_mean_variance(phylo = phylo,
                                   times_shared = times_shared,
                                   distances_phylo = distances_phylo,
                                   process = process,
                                   params_old = params_old,
                                   masque_data = masque_data,
                                   F_moments = F_moments)
  log_likelihood_old <- compute_log_likelihood(phylo = phylo,
                                               Y_data_vec = Y_data_vec_known,
                                               sim = moments$sim,
                                               Sigma = moments$Sigma,
                                               Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                               miss = miss, 
                                               masque_data = masque_data,
                                               C_YY = F_moments$C_YY,
                                               Y_data = Y_data,
                                               C_YY_chol_inv = F_moments$C_YY_chol_inv,
                                               R = params_old$variance)
  ## Compute Mahalanobis norm between data and mean at tips
  maha_data_mean <- compute_mahalanobis_distance(phylo = phylo,
                                                 Y_data_vec = Y_data_vec_known,
                                                 sim = moments$sim,
                                                 Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                 miss = miss,
                                                 Y_data = Y_data,
                                                 C_YY_chol_inv = F_moments$C_YY_chol_inv,
                                                 R = params_old$variance)
  ########## E step #########################################################
  conditional_law_X <- compute_E.simple(phylo = phylo,
                                        Y_data_vec = Y_data_vec_known,
                                        sim = moments$sim,
                                        Sigma = moments$Sigma,
                                        Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                        miss = miss,
                                        masque_data = masque_data,
                                        F_means = F_moments$F_means,
                                        F_vars = F_moments$F_vars,
                                        R = params_old$variance,
                                        Y_data = Y_data)
  return(list(log_likelihood_old = log_likelihood_old,
              maha_data_mean = maha_data_mean,
              conditional_law_X = conditional_law_X))
}

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
                          Y_data = Y_data)
  
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
                          Y_data = Y_data)
  
  expect_that(res_1$log_likelihood_old,
              equals(sum(sapply(res_2, function(z) return(z$log_likelihood_old)))))
  expect_that(res_1$maha_data_mean,
              equals(sum(sapply(res_2, function(z) return(z$maha_data_mean)))))
})