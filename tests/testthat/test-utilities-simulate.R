context("Function simulate with shifts")

###############################################################################
test_that("Mean of the BM", {
  set.seed(586)
  ntaxa <- 20
  tree <- rtree(ntaxa)
  
  ## Simulate Process
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  root.state <- list(random = FALSE,
                     value.root = c(1, -1, 2),
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(18, 32),
                values=cbind(c(4, -10, 3),
                             c(-5, 5, 0)),
                relativeTimes = 0)
  
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "BM",
                          variance = variance,
                          shifts = shifts)
  
  X1.tips.exp <- extract_simulate_internal(X1, where = "tips", what = "exp")
  
  ## Compute expectations with tree matrix
  T_tree <- incidence.matrix(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  
  X1.tips.exp.mat <- tcrossprod(Delta, T_tree) + root.state$value.root
  
  expect_that(X1.tips.exp, equals(X1.tips.exp.mat))
  
  ## Compute expectations of internal nodes
  U_tree <- incidence.matrix.full(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  X1.all.exp.mat <- tcrossprod(Delta, U_tree) + root.state$value.root
  
  X1.nodes.exp <- extract_simulate_internal(X1, where = "nodes", what = "exp")
  X1.all.exp <- cbind(X1.tips.exp, X1.nodes.exp)
  
  expect_that(X1.all.exp.mat, equals(X1.all.exp))
  
  ## Without simulating
  X2 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "BM",
                          variance = variance,
                          shifts = shifts,
                          simulate_random = FALSE)
  X2.tips.exp <- extract_simulate_internal(X2, where = "tips", what = "exp")
  expect_that(X1.tips.exp, equals(X2.tips.exp))
})

###############################################################################
test_that("Mean of the BM - random root", {
  set.seed(586)
  ntaxa <- 20
  tree <- rtree(ntaxa)
  
  ## Simulate Process
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  root.state <- list(random = TRUE,
                     value.root = NA,
                     exp.root = c(1, -1, 2),
                     var.root = matrix(0.1, p, p) + diag(0.1, p, p))
  shifts = list(edges = c(18, 32),
                values=cbind(c(4, -10, 3),
                             c(-5, 5, 0)),
                relativeTimes = 0)
  
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "BM",
                          variance = variance,
                          shifts = shifts)
  
  X1.tips.exp <- extract_simulate_internal(X1, where = "tips", what = "exp")
  
  ## Does adding a root edge change anything ?
  tree$root.edge <- 1
  X2 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "BM",
                          variance = variance,
                          shifts = shifts)
  
  X2.tips.exp <- extract_simulate_internal(X2, where = "tips", what = "exp")
  expect_that(X1.tips.exp, equals(X2.tips.exp))
  
  ## Compute expectations with tree matrix
  T_tree <- incidence.matrix(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  
  X1.tips.exp.mat <- tcrossprod(Delta, T_tree) + root.state$exp.root
  
  expect_that(X1.tips.exp, equals(X1.tips.exp.mat))
  
  ## Compute expectations of internal nodes
  U_tree <- incidence.matrix.full(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  X1.all.exp.mat <- tcrossprod(Delta, U_tree) + root.state$exp.root
  
  X1.nodes.exp <- extract_simulate_internal(X1, where = "nodes", what = "exp")
  X1.all.exp <- cbind(X1.tips.exp, X1.nodes.exp)
  
  expect_that(X1.all.exp.mat, equals(X1.all.exp))
  
  ## Without simulating
  X2 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "BM",
                          variance = variance,
                          shifts = shifts,
                          simulate_random = FALSE,
                          U_tree = U_tree)
  X2.tips.exp <- extract_simulate_internal(X2, where = "tips", what = "exp")
  expect_that(X1.tips.exp, equals(X2.tips.exp))
})

###############################################################################
test_that("Mean of the OU", {
  set.seed(1899)
  ntaxa <- 32
  tree <- rcoal(ntaxa)
  
  ## Simulate Process
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  optimal.value <- c(-3, 5, 0)
  selection.strength <- diag(3, p, p) + tcrossprod(c(0.1, 0.2, 0.3))
  exp.stationary <- optimal.value
  var.stationary  <- compute_stationary_variance(variance, selection.strength)
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = exp.stationary,
                     var.root = var.stationary)
  shifts = list(edges = c(11, 35, 44),
                values=cbind(c(4, -10, 3),
                             c(-5, 5.4, 0),
                             c(2, -8, 0.3)),
                relativeTimes = 0)
  
  expect_warning(X1 <- simulate_internal(tree,
                                         p = p,
                                         root.state = root.state,
                                         process = "OU",
                                         variance = variance,
                                         optimal.value = optimal.value,
                                         selection.strength = selection.strength,
                                         shifts = shifts),
                 "root variance")
  
  X1.tips.exp <- extract_simulate_internal(X1, where = "tips", what = "exp")
  
  ## Compute expectations with tree matrix
  T_tree <- incidence.matrix(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  W <- compute_actualization_matrix_ultrametric(tree, selection.strength)
  
  vec_Y <- kronecker(T_tree, diag(1, p, p)) %*% W %*% as.vector(Delta)
  
  X1.tips.exp.mat <- matrix(vec_Y, p, ntaxa) + optimal.value
  
  expect_that(X1.tips.exp, equals(X1.tips.exp.mat))
  
  ## Without simulate
  expect_warning(X2 <- simulate_internal(tree,
                                         p = p,
                                         root.state = root.state,
                                         process = "OU",
                                         variance = variance,
                                         optimal.value = optimal.value,
                                         selection.strength = selection.strength,
                                         shifts = shifts,
                                         simulate_random = FALSE),
                 "root variance")
  X2.tips.exp <- extract_simulate_internal(X2, where = "tips", what = "exp")
  expect_that(X2.tips.exp, equals(X2.tips.exp))
})

###############################################################################
test_that("OU - fixed root", {
  set.seed(1899)
  ntaxa <- 32
  tree <- rcoal(ntaxa)
  
  ## Simulate Process
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  optimal.value <- c(-3, 5, 0)
  selection.strength <- diag(3, p, p) + tcrossprod(c(0.1, 0.2, 0.3))
  exp.stationary <- optimal.value
  # var.stationary  <- compute_stationary_variance(variance, selection.strength)
  root.state <- list(random = FALSE,
                     stationary.root = FALSE,
                     value.root = exp.stationary,
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(11, 35, 44),
                values=cbind(c(4, -10, 3),
                             c(-5, 5.4, 0),
                             c(2, -8, 0.3)),
                relativeTimes = 0)
  
  X1 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "OU",
                          variance = variance,
                          optimal.value = optimal.value,
                          selection.strength = selection.strength,
                          shifts = shifts)
  
  X1.tips.exp <- extract_simulate_internal(X1, where = "tips", what = "exp")
  
  ## Compute expectations with tree matrix
  T_tree <- incidence.matrix(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  W <- compute_actualization_matrix_ultrametric(tree, selection.strength)
  
  vec_Y <- kronecker(T_tree, diag(1, p, p)) %*% W %*% as.vector(Delta)
  
  X1.tips.exp.mat <- matrix(vec_Y, p, ntaxa) + optimal.value
  
  expect_that(X1.tips.exp, equals(X1.tips.exp.mat))
  
  ## Without simulate
  X2 <- simulate_internal(tree,
                          p = p,
                          root.state = root.state,
                          process = "OU",
                          variance = variance,
                          optimal.value = optimal.value,
                          selection.strength = selection.strength,
                          shifts = shifts,
                          simulate_random = FALSE)
  X2.tips.exp <- extract_simulate_internal(X2, where = "tips", what = "exp")
  expect_that(X2.tips.exp, equals(X2.tips.exp))
})

###############################################################################
test_that("Multivariate Scalar (scOU)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(586)
  ntaxa <- 100
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
                                   lambda = 1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  
  p <- 6
  variance <- diag(0.5, p, p)
  optimal.value <- c(-3, 5, 0, 7, 3, 0)
  selection.strength <- diag(3, p, p)
  exp.stationary <- optimal.value
  var.stationary  <- compute_stationary_variance(variance, selection.strength)
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = exp.stationary,
                     var.root = var.stationary)
  shifts = list(edges = c(18, 32),
                values=cbind(c(4, -10, 3, 2, 2, 9),
                             c(-5, 5, 0, -7, 5, -4)),
                relativeTimes = 0)
  
  ## Forget that it is scalar
  Xnot <- simulate_internal(tree,
                            p = p,
                            root.state = root.state,
                            process = "OU",
                            variance = variance,
                            optimal.value = optimal.value,
                            selection.strength = selection.strength,
                            shifts = shifts)
  ## Use that it is scalar
  Xsc <- simulate_internal(tree,
                           p = p,
                           root.state = root.state,
                           process = "scOU",
                           variance = variance,
                           optimal.value = optimal.value,
                           selection.strength = selection.strength,
                           shifts = shifts)
  
  expect_that(Xnot[,,2:3], equals(Xsc[,,2:3]))
  ## Do not simulate
  Xsc2 <- simulate_internal(tree,
                            p = p,
                            root.state = root.state,
                            process = "scOU",
                            variance = variance,
                            optimal.value = optimal.value,
                            selection.strength = selection.strength,
                            shifts = shifts,
                            simulate_random = FALSE)
  
  expect_that(Xnot[,,2:3], equals(Xsc2[,,2:3]))
})

###############################################################################
test_that("Multivariate Scalar (scOU) - Fixed Root", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TreeSim")
  set.seed(586)
  ntaxa <- 100
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
                                   lambda = 1, mu = 0,
                                   age = 1, mrca = TRUE)[[1]]
  
  p <- 6
  variance <- diag(0.5, p, p)
  optimal.value <- c(-3, 5, 0, 7, 3, 0)
  selection.strength <- diag(3, p, p)
  exp.stationary <- optimal.value
  root.state <- list(random = FALSE,
                     stationary.root = FALSE,
                     value.root = exp.stationary,
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(18, 32),
                values=cbind(c(4, -10, 3, 2, 2, 9),
                             c(-5, 5, 0, -7, 5, -4)),
                relativeTimes = 0)
  
  ## Forget that it is scalar
  Xnot <- simulate_internal(tree,
                            p = p,
                            root.state = root.state,
                            process = "OU",
                            variance = variance,
                            optimal.value = optimal.value,
                            selection.strength = selection.strength,
                            shifts = shifts)
  ## Use that it is scalar
  Xsc <- simulate_internal(tree,
                           p = p,
                           root.state = root.state,
                           process = "scOU",
                           variance = variance,
                           optimal.value = optimal.value,
                           selection.strength = selection.strength,
                           shifts = shifts)
  
  expect_that(Xnot[,,2:3], equals(Xsc[,,2:3]))
  ## Do not simulate
  Xsc2 <- simulate_internal(tree,
                            p = p,
                            root.state = root.state,
                            process = "scOU",
                            variance = variance,
                            optimal.value = optimal.value,
                            selection.strength = selection.strength,
                            shifts = shifts,
                            simulate_random = FALSE)
  
  expect_that(Xnot[,,2:3], equals(Xsc2[,,2:3]))
})

###############################################################################
test_that("Interval vs simul", {
  ntaxa <- 32
  tree <- rcoal(ntaxa)
  
  ## Simulate Internal
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  optimal.value <- c(-3, 5, 0)
  selection.strength <- diag(3, p, p) + tcrossprod(c(0.1, 0.2, 0.3))
  exp.stationary <- optimal.value
  var.stationary  <- compute_stationary_variance(variance, selection.strength)
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = exp.stationary,
                     var.root = var.stationary)
  shifts = list(edges = c(11, 35, 44),
                values=cbind(c(4, -10, 3),
                             c(-5, 5.4, 0),
                             c(2, -8, 0.3)),
                relativeTimes = 0)
  
  set.seed(1899)
  expect_warning(X1 <- simulate_internal(tree,
                                         p = p,
                                         root.state = root.state,
                                         process = "OU",
                                         variance = variance,
                                         optimal.value = optimal.value,
                                         selection.strength = selection.strength,
                                         shifts = shifts),
                 "root variance")
  
  ## Simulate External
  expect_warning(para <- params_process("OU", p = p, variance = variance,
                                        selection.strength = selection.strength,
                                        optimal.value = optimal.value, random = TRUE,
                                        stationary.root = TRUE, exp.root = exp.stationary,
                                        var.root = var.stationary, edges = shifts$edges,
                                        values = shifts$values),
                 "root variance")
  
  set.seed(1899)
  X2 <- expect_warning(simul_process(para, phylo = tree), "root variance")
  
  ## Comparisons and extractors
  X1.tips.exp <- extract_simulate_internal(X1, where = "tips", what = "exp")
  X2.tips.exp <- extract(X2, where = "tips", what = "exp")
  expect_equal(X1.tips.exp, unname(X2.tips.exp))
  
  X1.nodes.exp <- extract_simulate_internal(X1, where = "nodes", what = "exp")
  X2.nodes.exp <- extract(X2, where = "nodes", what = "exp")
  expect_equal(X1.nodes.exp, unname(X2.nodes.exp))
  
  X1.tips.states <- extract_simulate_internal(X1, where = "tips", what = "states")
  X2.tips.states <- extract(X2, where = "tips", what = "states")
  expect_equal(X1.tips.states, unname(X2.tips.states))
  
  X1.nodes.states <- extract_simulate_internal(X1, where = "nodes", what = "states")
  X2.nodes.states <- extract(X2, where = "nodes", what = "states")
  expect_equal(X1.nodes.states, unname(X2.nodes.states))
  
  X1.tips.optim <- extract_simulate_internal(X1, where = "tips", what = "optim")
  X2.tips.optim <- extract(X2, where = "tips", what = "optim")
  expect_equal(X1.tips.optim, unname(X2.tips.optim))
  
  X1.nodes.optim <- extract_simulate_internal(X1, where = "nodes", what = "optim")
  X2.nodes.optim <- extract(X2, where = "nodes", what = "optim")
  expect_equal(X1.nodes.optim, unname(X2.nodes.optim))
})

# test_that("Multivariate Independent (BM)", {
#   set.seed(586)
#   ntaxa <- 200
#   tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
#                           lambda = 1, mu = 0,
#                           age = 1, mrca = TRUE)[[1]]
#   
#   p <- 3
#   variance <- diag(0.5, p, p)
#   root.state <- list(random = FALSE,
#                      stationary.root = FALSE,
#                      value.root = c(1, -6, 2),
#                      exp.root = NA,
#                      var.root = NA)
#   shifts = list(edges = c(18, 32),
#                 values=cbind(c(4, -10, 3),
#                              c(-5, 5, 0)),
#                 relativeTimes = 0)
#   paramsSimu <- list(variance = variance,
#                      selection.strength = selection.strength,
#                      shifts = shifts,
#                      root.state = root.state)
#   
#   Xnot <- simulate_internal(tree,
#                    p = p,
#                    independent = FALSE,
#                    root.state = root.state,
#                    process = "BM",
#                    variance = variance,
#                    optimal.value = optimal.value,
#                    selection.strength = selection.strength,
#                    shifts = shifts)
#   
#   Xind <- simulate_internal(tree,
#                    p = p,
#                    independent = TRUE,
#                    root.state = root.state,
#                    process = "OU",
#                    variance = variance,
#                    optimal.value = optimal.value,
#                    selection.strength = selection.strength,
#                    shifts = shifts)
#   
#   expect_that(Xnot[,,2], equals(Xind[,,2]))
# })

###############################################################################
test_that("OU/BM", {
  set.seed(1899)
  ntaxa <- 32
  tree <- rcoal(ntaxa)
  
  ## Simulate Process
  p <- 3
  variance <- matrix(0.2, p, p) + diag(0.3, p, p)
  optimal.value <- c(-3, 5, 0)
  selection.strength <- diag(0:2)
  exp.stationary <- optimal.value
  # var.stationary  <- compute_stationary_variance(variance, selection.strength)
  root.state <- list(random = FALSE,
                     stationary.root = FALSE,
                     value.root = exp.stationary,
                     exp.root = NA,
                     var.root = NA)
  shifts = NULL
  
  X1 <- expect_warning(simulate_internal(tree,
                                         p = p,
                                         root.state = root.state,
                                         process = "OUBM",
                                         variance = variance,
                                         optimal.value = optimal.value,
                                         selection.strength = selection.strength,
                                         shifts = shifts),
                       "All the eigen values of the selection strengh do not have a strictly positive real part.")
  
  X1.tips.exp <- extract_simulate_internal(X1, where = "tips", what = "exp")
  
  ## Compute expectations with tree matrix
  T_tree <- incidence.matrix(tree)
  expect_error(shifts.list_to_matrix(tree, shifts), "the dimension p must be specified when shift is NULL")
  Delta <- shifts.list_to_matrix(tree, shifts, p = p)
  W <- compute_actualization_matrix_ultrametric(tree, selection.strength)
  
  vec_Y <- kronecker(T_tree, diag(1, p, p)) %*% W %*% as.vector(Delta)
  
  X1.tips.exp.mat <- matrix(vec_Y, p, ntaxa) + optimal.value
  
  expect_that(X1.tips.exp, equals(X1.tips.exp.mat))
  
  ## Without simulate
  expect_warning(X2 <- simulate_internal(tree,
                                         p = p,
                                         root.state = root.state,
                                         process = "OU",
                                         variance = variance,
                                         optimal.value = optimal.value,
                                         selection.strength = selection.strength,
                                         shifts = shifts,
                                         simulate_random = FALSE),
                 "strictly positive real part")
  X2.tips.exp <- extract_simulate_internal(X2, where = "tips", what = "exp")
  expect_that(X2.tips.exp, equals(X2.tips.exp))
  
  ## Variances
  set.seed(1899)
  ntaxa <- 5
  p <- 2
  tree <- rcoal(ntaxa)
  variance <- diag(0.5, p)
  optimal.value <- c(-3, 5)
  selection.strength <- diag(0:1)
  exp.stationary <- optimal.value
  root.state <- list(random = FALSE,
                     stationary.root = FALSE,
                     value.root = exp.stationary,
                     exp.root = NA,
                     var.root = NA)
  tree_heigth <- max(node.depth.edgelength(tree))
  var_1 <- tree_heigth * variance[1, 1]
  var_2 <- variance[2, 2] / (2 * selection.strength[2, 2])
  X1.tips.states <- sapply(1:50, function(x) extract_simulate_internal(expect_warning(simulate_internal(tree,
                                                                                                        p = p,
                                                                                                        root.state = root.state,
                                                                                                        process = "OUBM",
                                                                                                        variance = variance,
                                                                                                        optimal.value = optimal.value,
                                                                                                        selection.strength = selection.strength,
                                                                                                        shifts = shifts),
                                                                                      "strictly positive real part"),
                                                                       where = "tips", what = "states"))
  X1.tips.states.var <- apply(X1.tips.states, 1, var)
  expect_equal(mean(X1.tips.states.var[1 + p * 0:(ntaxa - 1)]), var_1, 0.05)
  expect_equal(mean(X1.tips.states.var[2 + p * 0:(ntaxa - 1)]), var_2, 0.05)
})
