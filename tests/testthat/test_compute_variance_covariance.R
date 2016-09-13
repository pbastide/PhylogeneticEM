context("Covariances computations")

test_that("OU and scOU, fixed root", {
  set.seed(586)
  ntaxa <- 20
  tree <- ape::rtree(ntaxa)
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  
  ## Process parameters
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  root.state <- list(random = FALSE,
                     stationary.root = FALSE,
                     value.root = c(1, -1, 2),
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(18, 32),
                values=cbind(c(4, -10, 3),
                             c(-5, 5, 0)),
                relativeTimes = 0)
  alpha <- 3
  
  # OU style
  params_OU <- list(variance = variance,
                    root.state = root.state,
                    shifts = shifts,
                    selection.strength = diag(rep(alpha, p)))
  varr_OU <- compute_variance_covariance.OU(times_shared, distances_phylo,
                                            params_OU)
  
  # sOU style
  params_scOU <- list(variance = variance,
                    root.state = root.state,
                    shifts = shifts,
                    selection.strength = alpha)
  varr_scOU <- compute_variance_covariance.scOU(times_shared, distances_phylo,
                                                params_scOU)
  
  expect_that(as.vector(varr_OU), equals(as.vector(varr_scOU)))
})

test_that("OU and scOU, random root", {
  set.seed(586)
  ntaxa <- 20
  tree <- ape::rtree(ntaxa)
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  
  ## Process parameters
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  root.state <- list(random = TRUE,
                     stationary.root = FALSE,
                     value.root = c(1, -1, 2),
                     exp.root = NA,
                     var.root = matrix(0.1, p, p))
  shifts = list(edges = c(18, 32),
                values=cbind(c(4, -10, 3),
                             c(-5, 5, 0)),
                relativeTimes = 0)
  alpha <- 3
  
  # OU style
  params_OU <- list(variance = variance,
                    root.state = root.state,
                    shifts = shifts,
                    selection.strength = diag(rep(alpha, p)))
  varr_OU <- compute_variance_covariance.OU(times_shared, distances_phylo, 
                                            params_OU)
  
  # sOU style
  params_scOU <- list(variance = variance,
                      root.state = root.state,
                      shifts = shifts,
                      selection.strength = alpha)
  varr_scOU <- compute_variance_covariance.scOU(times_shared, distances_phylo, params_scOU)
  
  expect_that(as.vector(varr_OU), equals(as.vector(varr_scOU)))
})

test_that("OU and scOU, stationary root", {
  set.seed(586)
  ntaxa <- 20
  tree <- ape::rtree(ntaxa)
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  
  ## Process parameters
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  alpha <- 3
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = c(1, -1, 2),
                     exp.root = NA,
                     var.root = 1 / (2 * alpha) * variance)
  shifts = list(edges = c(18, 32),
                values=cbind(c(4, -10, 3),
                             c(-5, 5, 0)),
                relativeTimes = 0)
  
  # OU style
  params_OU <- list(variance = variance,
                    root.state = root.state,
                    shifts = shifts,
                    selection.strength = diag(rep(alpha, p)))
  varr_OU <- compute_variance_covariance.OU(times_shared, distances_phylo,
                                            params_OU)
  
  # sOU style
  params_scOU <- list(variance = variance,
                      root.state = root.state,
                      shifts = shifts,
                      selection.strength = alpha)
  varr_scOU <- compute_variance_covariance.scOU(times_shared, distances_phylo, params_scOU)
  
  expect_that(as.vector(varr_OU), equals(as.vector(varr_scOU)))
})