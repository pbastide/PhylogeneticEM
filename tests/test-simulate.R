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
  ntaxa <- 64
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