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
                               var.root.init = diag((1:p)/(2 * 1:p)),
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
  testthat::skip_on_cran()
  # Dimension p - OU- stationary root
  set.seed(586)
  ntaxa <- 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
                                   lambda = 0.1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  U_tree <- incidence.matrix.full(tree)
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
                          compute_E = compute_E.simple,
                          U_tree = U_tree)
  
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
                          compute_E = compute_E.simple,
                          U_tree = U_tree)
  
  expect_equal(as.vector(res_1$log_likelihood_old),
              sum(sapply(res_2, function(z) return(z$log_likelihood_old))))
  # expect_equal(as.vector(res_1$maha_data_mean),
  #             sum(sapply(res_2, function(z) return(z$maha_data_mean))))
})

test_that("Independent OU - uni/multi", {
  testthat::skip_on_cran()
  set.seed(586)
  ntaxa <- 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
                                   lambda = 0.1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  U_tree <- incidence.matrix.full(tree)
  p <- 2
  
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
  
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = params$root.state,
                          process = "OU",
                          variance = params$variance,
                          shifts = params$shifts,
                          selection.strength = params$selection.strength,
                          optimal.value = params$optimal.value)
  
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  
  expect_warning(res <- PhyloEM(phylo = tree,
                                Y_data = traits,
                                process = "OU",
                                K_max = 0,
                                random.root = TRUE,
                                stationary.root = TRUE,
                                independent = TRUE,
                                save_step = FALSE,
                                Nbr_It_Max = 5,
                                method.variance = "upward_downward",
                                method.init = "lasso",
                                use_previous = FALSE,
                                method.selection = "LINselect",
                                impute_init_Rphylopars = FALSE,
                                progress.bar = FALSE,
                                K_lag_init = 0,
                                check.tips.names = FALSE,
                                methods.segmentation = c(#"max_costs_0", 
                                                         "lasso", 
                                                         "same_shifts", 
                                                         #"same_shifts_same_values",
                                                         "best_single_move"
                                                         #"lasso_one_move"
                                                         )))
  
  expect_warning(res_uni <- PhyloEM(phylo = tree,
                                    Y_data = traits[1, ],
                                    process = "OU",
                                    K_max = 0,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    alpha_grid = FALSE,
                                    save_step = FALSE,
                                    Nbr_It_Max = 5,
                                    method.variance = "upward_downward",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "LINselect",
                                    impute_init_Rphylopars = FALSE,
                                    progress.bar = FALSE,
                                    K_lag_init = 0,
                                    check.tips.names = FALSE))
  
  ## Test extractors
  par_multi <- res$alpha_max$params_estim$`0`
  par_multi_bis <- params_process(res, K = 0)
  class(par_multi) <- "params_process"
  par_multi$process <- "OU"
  
  expect_equal(par_multi, par_multi_bis)
  
  par_uni <- res_uni$alpha_max$params_estim$`0`
  par_uni_bis <- params_process(res_uni, K = 0)
  class(par_uni) <- "params_process"
  par_uni$process <- "OU"
  
  expect_equal(par_uni_bis, par_uni)
  
  ## Test parameters K=0
  par_multi <- split_params_independent(res$alpha_max$params_estim$`0`)[[1]]
  par_uni <- res_uni$alpha_max$params_estim$`0`
  attr(par_uni, "log_likelihood") <- NULL
  attr(par_uni, "Neq") <- NULL
  
  expect_equal(par_multi, par_uni)
  
  ## Test imputed states
  tmp <- imputed_traits(res, K = 0, what = "expectations", trait = 1)
  tmp_uni <- imputed_traits(res_uni, K = 0, what = "expectations")
  expect_equal(tmp, tmp_uni)
  
  tmp <- imputed_traits(res, K = 0, what = "variances", trait = 1)
  tmp_uni <- imputed_traits(res_uni, K = 0, what = "variances")
  expect_equal(tmp, tmp_uni)
  
  tmp <- imputed_traits(res, K = 0, what = "imputed", trait = 1)
  tmp_uni <- imputed_traits(res_uni, K = 0, what = "imputed")
  expect_equal(tmp, tmp_uni)
  
  ## Simple test of grid univariate case
  traits[1, c(10, 15, 16)] <- NA
  expect_warning(res_uni <- PhyloEM(phylo = tree,
                                    Y_data = traits[1, ],
                                    process = "OU",
                                    K_max = 1,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    alpha_grid = TRUE,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "upward_downward",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "LINselect",
                                    impute_init_Rphylopars = FALSE,
                                    progress.bar = FALSE,
                                    K_lag_init = 0,
                                    check.tips.names = FALSE))
  
  res_uni_bis <- model_selection(res_uni, "BGHuni")
  
  expect_equal(res_uni, res_uni_bis)
})

test_that("Independent BM - uni/multi", {
  testthat::skip_on_cran()
  set.seed(586)
  ntaxa <- 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
                                   lambda = 0.1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  U_tree <- incidence.matrix.full(tree)
  p <- 2
  
  params <- init.EM.default.BM(p = p,
                               variance.init = diag(1:p),
                               random.init = FALSE,
                               value.root.init = 1:p,
                               exp.root.init = 1:p,
                               var.root.init = diag(1:p),
                               edges.init = c(11, 14, 23),
                               values.init = matrix(1:(p*3), p, 3))
  attr(params, "p_dim") <- p
  
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = params$root.state,
                          process = "BM",
                          variance = params$variance,
                          shifts = params$shifts)
  
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  
  expect_warning(res <- PhyloEM(phylo = tree,
                                Y_data = traits,
                                process = "BM",
                                K_max = 0,
                                random.root = FALSE,
                                independent = TRUE,
                                save_step = FALSE,
                                Nbr_It_Max = 50,
                                method.variance = "upward_downward",
                                method.init = "lasso",
                                use_previous = FALSE,
                                method.selection = "LINselect",
                                impute_init_Rphylopars = FALSE,
                                progress.bar = FALSE,
                                K_lag_init = 0,
                                check.tips.names = FALSE))
  
  res_uni <- PhyloEM(phylo = tree,
                     Y_data = traits[1, ],
                     process = "BM",
                     K_max = 0,
                     random.root = FALSE,
                     alpha_grid = FALSE,
                     save_step = FALSE,
                     Nbr_It_Max = 50,
                     method.variance = "upward_downward",
                     method.init = "lasso",
                     use_previous = FALSE,
                     method.selection = "LINselect",
                     impute_init_Rphylopars = FALSE,
                     progress.bar = FALSE,
                     K_lag_init = 0,
                     check.tips.names = FALSE)
  
  ## Test extractors
  par_multi <- res$alpha_max$params_estim$`0`
  par_multi_bis <- params_process(res, K = 0)
  class(par_multi) <- "params_process"
  par_multi$process <- "BM"
  
  expect_equal(par_multi, par_multi_bis)
  
  par_uni <- res_uni$alpha_max$params_estim$`0`
  par_uni_bis <- params_process(res_uni, K = 0)
  class(par_uni) <- "params_process"
  par_uni$process <- "BM"
  
  expect_equal(par_uni_bis, par_uni)
  
  ## Test parameters K=0
  par_multi <- split_params_independent(res$alpha_max$params_estim$`0`)[[1]]
  expect_null(par_multi$optimal.value)
  par_multi$optimal.value <- NULL
  
  par_uni <- res_uni$alpha_max$params_estim$`0`
  attr(par_uni, "ntaxa") <- NULL
  attr(par_uni, "log_likelihood") <- NULL
  attr(par_uni, "Neq") <- NULL
  
  expect_equal(par_multi, par_uni, tolerance = 10^(-3))
  
  ## Test imputed states
  tmp <- imputed_traits(res, K = 0, what = "expectations", trait = 1)
  tmp_uni <- imputed_traits(res_uni, K = 0, what = "expectations")
  expect_equal(tmp, tmp_uni, tolerance = 10^(-3))
  
  tmp <- imputed_traits(res, K = 0, what = "variances", trait = 1)
  tmp_uni <- imputed_traits(res_uni, K = 0, what = "variances")
  expect_equal(tmp, tmp_uni, tolerance = 10^(-3))
  
  tmp <- imputed_traits(res, K = 0, what = "imputed", trait = 1)
  tmp_uni <- imputed_traits(res_uni, K = 0, what = "imputed")
  expect_equal(tmp, tmp_uni, tolerance = 10^(-3))
  
  ## Simple test of grid univariate case
  traits[1, c(10, 15, 16)] <- NA
  expect_warning(res_uni <- PhyloEM(phylo = tree,
                                    Y_data = traits[1, ],
                                    process = "BM",
                                    K_max = 1,
                                    random.root = TRUE,
                                    stationary.root = TRUE,
                                    # independent = TRUE,
                                    save_step = FALSE,
                                    Nbr_It_Max = 2,
                                    method.variance = "upward_downward",
                                    method.init = "lasso",
                                    use_previous = FALSE,
                                    method.selection = "LINselect",
                                    impute_init_Rphylopars = FALSE,
                                    progress.bar = FALSE,
                                    K_lag_init = 0,
                                    check.tips.names = FALSE))
  
  res_uni_bis <- model_selection(res_uni, "BGHuni")
  
  expect_equal(res_uni, res_uni_bis)
})