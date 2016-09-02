context("E step upward-downward")

test_that("Upward Downward", {
  set.seed(17920902)
  ntaxa = 100
  tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
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
  
  X1 <- simulate(tree,
                 p = p,
                 root.state = root.state,
                 process = "BM",
                 variance = variance,
                 shifts = shifts)
  
  traits <- extract.simulate(X1, where = "tips", what = "state")
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


test_that("Upward Downward - no missing", {
  set.seed(17920902)
  ntaxa = 500
  tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                          age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  p <- 6
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
  paramsSimu <- list(variance = variance,
                     shifts = shifts,
                     root.state = root.state)
  attr(paramsSimu, "p_dim") <- p
  
  X1 <- simulate(tree,
                 p = p,
                 root.state = root.state,
                 process = "BM",
                 variance = variance,
                 shifts = shifts)
  traits <- extract.simulate(X1, where = "tips", what = "state")
  
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
})