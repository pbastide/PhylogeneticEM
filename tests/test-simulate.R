context("Function simulate with shifts")

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
  
  X1 <- simulate(tree,
                 p = p,
                 root.state = root.state,
                 process = "BM",
                 variance = variance,
                 shifts = shifts)
  
  X1.tips.exp <- extract.simulate(X1, where = "tips", what = "exp")
  
  ## Compute expectations with tree matrix
  T_tree <- incidence.matrix(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  
  X1.tips.exp.mat <- tcrossprod(Delta, T_tree) + root.state$value.root
  
  expect_that(X1.tips.exp, equals(X1.tips.exp.mat))
})

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
  var.stationary  <- compute_stationnary_variance(variance, selection.strength)
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
  
  X1 <- simulate(tree,
                 p = p,
                 root.state = root.state,
                 process = "OU",
                 variance = variance,
                 optimal.value = optimal.value,
                 selection.strength = selection.strength,
                 shifts = shifts)
  
  X1.tips.exp <- extract.simulate(X1, where = "tips", what = "exp")
  
  ## Compute expectations with tree matrix
  T_tree <- incidence.matrix(tree)
  Delta <- shifts.list_to_matrix(tree, shifts)
  W <- compute_actualization_matrix_ultrametric(tree, selection.strength)
  
  vec_Y <- kronecker(T_tree, diag(1, p, p)) %*% W %*% as.vector(Delta)
  
  X1.tips.exp.mat <- matrix(vec_Y, p, ntaxa) + optimal.value
  
  expect_that(X1.tips.exp, equals(X1.tips.exp.mat))
})
# 
# test_that("Multivariate Independent (OU)", {
#   set.seed(586)
#   ntaxa <- 1000
#   tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
#                           lambda = 1, mu = 0,
#                           age = 1, mrca = TRUE)[[1]]
#   
#   p <- 6
#   variance <- diag(0.5, p, p)
#   optimal.value <- c(-3, 5, 0, 7, 3, 0)
#   selection.strength <- diag(3, p, p)
#   exp.stationary <- optimal.value
#   var.stationary  <- compute_stationnary_variance(variance, selection.strength)
#   root.state <- list(random = TRUE,
#                      stationary.root = TRUE,
#                      value.root = NA,
#                      exp.root = exp.stationary,
#                      var.root = var.stationary)
#   shifts = list(edges = c(18, 32),
#                 values=cbind(c(4, -10, 3, 2, 2, 9),
#                              c(-5, 5, 0, -7, 5, -4)),
#                 relativeTimes = 0)
#   paramsSimu <- list(variance = variance,
#                      optimal.value = optimal.value,
#                      selection.strength = selection.strength,
#                      shifts = shifts,
#                      root.state = root.state)
#   
#   tempsnot <- system.time(Xnot <- simulate(tree,
#                  p = p,
#                  independent = FALSE,
#                  root.state = root.state,
#                  process = "OU",
#                  variance = variance,
#                  optimal.value = optimal.value,
#                  selection.strength = selection.strength,
#                  shifts = shifts))
#   
#   tempsind <- system.time(Xind <- simulate(tree,
#                    p = p,
#                    independent = TRUE,
#                    root.state = root.state,
#                    process = "OU",
#                    variance = variance,
#                    optimal.value = optimal.value,
#                    selection.strength = selection.strength,
#                    shifts = shifts))
#   
#   expect_that(Xnot[,,2:3], equals(Xind[,,2:3]))
# })
# 
# test_that("Multivariate Independent (BM)", {
#   set.seed(586)
#   ntaxa <- 200
#   tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, 
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
#   Xnot <- simulate(tree,
#                    p = p,
#                    independent = FALSE,
#                    root.state = root.state,
#                    process = "BM",
#                    variance = variance,
#                    optimal.value = optimal.value,
#                    selection.strength = selection.strength,
#                    shifts = shifts)
#   
#   Xind <- simulate(tree,
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