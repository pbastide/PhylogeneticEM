context("Functions handling parameters")

test_that("check_dimensions.root.state ", {
  # Dimension 1 - structure
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = 3,
                          exp.root = 2,
                          var.root = 5)
  root.state_correct <- list(random = TRUE,
                             stationary.root = TRUE, 
                             value.root = 3,
                             exp.root = 2,
                             var.root = as(Matrix(5, 1, 1), "dpoMatrix"))
  
  root.state_test <- check_dimensions.root.state(1, root.state_test)
  expect_that(root.state_test, equals(root.state_correct))
  
  # Dimension p - var
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = rep(3, 6),
                          exp.root = rep(3, 6),
                          var.root = matrix(1, 5, 6))
  
  expect_that(check_dimensions.root.state(6, root.state_test),
              throws_error())
  
  # Dimension p - exp
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = rep(3, 6),
                          exp.root = rep(3, 7),
                          var.root = matrix(1, 6, 6))
  
  expect_that(check_dimensions.root.state(6, root.state_test),
              throws_error())
  
  # Dimension p - value
  root.state_test <- list(random = FALSE,
                          stationary.root = TRUE, 
                          value.root = rep(3, 5),
                          exp.root = rep(3, 6),
                          var.root = matrix(1, 6, 6))
  
  expect_that(check_dimensions.root.state(6, root.state_test),
              throws_error())
})

test_that("check_dimensions.shifts ", {
  # Dimension 1 - structure
  shifts_test = list(edges = c(18, 32, 45),
                     values = c(6, 4, -2),
                     relativeTimes = 0)
  
  shifts_correct = list(edges = c(18, 32, 45),
                        values = matrix(c(6, 4, -2), 1, 3),
                        relativeTimes = c(0, 0, 0))
  
  shifts_test <- check_dimensions.shifts(1, shifts_test)
  expect_that(shifts_test, equals(shifts_correct))
  
  # Dimension p - Everything fine
  p <- 5
  shifts_test = list(edges = c(18, 32, 45),
                     values = matrix(rep(1, 3*p), p, 3),
                     relativeTimes = c(0.5, 0, 0.2))
  
  expect_that(shifts_test, 
              equals(check_dimensions.shifts(p, shifts_test)))
  
  # Dimension p - number of shifts
  p <- 5
  shifts_test = list(edges = c(18, 32, 45),
                     values = matrix(rep(1, 2*p), p, 2),
                     relativeTimes = 0)
  
  expect_that(check_dimensions.shifts(p, shifts_test),
              throws_error())
  
  # Dimension p - dimension
  p <- 5
  shifts_test = list(edges = c(18, 32, 45),
                     values = matrix(rep(1, 3*(p+2)), p + 2, 3),
                     relativeTimes = 0)
  
  expect_that(check_dimensions.shifts(p, shifts_test),
              throws_error())
  
  # Dimension p - relativeTimes
  p <- 5
  shifts_test = list(edges = c(18, 32, 45),
                     values = matrix(rep(1, 3*p), p, 3),
                     relativeTimes = c(0.5, 0))
  
  expect_that(check_dimensions.shifts(p, shifts_test),
              throws_error())
})

test_that("check_dimensions.matrix", {
  # Dimension 1 - structure
  variance <- 23.4
  
  expect_that(matrix(variance, 1, 1), 
              equals(check_dimensions.matrix(1, 1, variance)))
  
  # Dimension 1 - problem
  variance <- c(23.4, 21.8)
  
  expect_that(check_dimensions.matrix(1, 1, variance), 
              throws_error())
  
  # Dimension p - Everything fine
  p <- 5
  q <- 3
  variance <- matrix(1.3, p, q)
  
  expect_that(variance, 
              equals(check_dimensions.matrix(p, q, variance)))
  
  # Dimension p - mismatch
  p <- 5
  q <- 2
  variance <- matrix(1.3, p, q-1)
  
  expect_that(check_dimensions.matrix(p, q, variance), 
              throws_error())
  
  # Dimension p - wrong dimension
  p <- 5
  q <- 2
  variance <- matrix(1.3, p-1, q)
  
  expect_that(check_dimensions.matrix(p, q, variance), 
              throws_error())
})

test_that("test.root.state", {
  ## BM
  # Dimension 1
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = 3,
                          exp.root = 2,
                          var.root = as.matrix(5, 1, 1))
  
  root.state_correct <- list(random = TRUE,
                             value.root = NA,
                             exp.root = 2,
                             var.root = as(as.matrix(5, 1, 1), "dpoMatrix"))
  
  expect_warning(root.state_test <- test.root.state(root.state_test, "BM"))
  expect_that(root.state_test, equals(root.state_correct))
  
  # Dimension p
  p <- 5
  root.state_test <- list(random = FALSE,
                          stationary.root = TRUE, 
                          value.root = rep(3, p),
                          exp.root = rep(2, p),
                          var.root = as.matrix(5, p, p))
  
  root.state_correct <- list(random = FALSE,
                             value.root = rep(3, p),
                             exp.root = NA,
                             var.root = NA)
  
  expect_warning(root.state_test <- test.root.state(root.state_test, "BM"))
  expect_that(root.state_test, equals(root.state_correct))
  
  ## OU
  # Dimension 1
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = 3,
                          exp.root = 2,
                          var.root = as.matrix(5, 1, 1))
  optimal.value <- 2
  selection.strength <- matrix(3, 1, 1)
  variance <- compute_variance_from_stationary(root.state_test$var.root, selection.strength)
  root.state_correct <- list(random = TRUE,
                             stationary.root = TRUE, 
                             value.root = NA,
                             exp.root = 2,
                             var.root = as(as.matrix(5, 1, 1), "dpoMatrix"))
  
  expect_warning(root.state_test <- test.root.state(root.state_test, "OU",
                                                    optimal.value = optimal.value,
                                                    variance = variance,
                                                    selection.strength = selection.strength),
                 "As root state is supposed random, its value is not defined and set to NA")
  expect_that(root.state_test, equals(root.state_correct))
  
  ## OU
  # Dimension p
  p <- 5
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = rep(3, p),
                          exp.root = rep(2, p),
                          var.root = as.matrix(5, p, p))
  optimal.value <- rep(2, p) 
  selection.strength <- matrix(3, 1, 1)
  variance <- compute_variance_from_stationary(root.state_test$var.root, selection.strength)
  root.state_correct <- list(random = TRUE,
                             stationary.root = TRUE, 
                             value.root = NA,
                             exp.root = rep(2, p),
                             var.root = as(as.matrix(5, p, p), "dpoMatrix"))
  
  expect_warning(root.state_test <- test.root.state(root.state_test, "OU",
                                                    optimal.value = optimal.value,
                                                    variance = variance,
                                                    selection.strength = selection.strength),
                 "As root state is supposed random, its value is not defined and set to NA")
  expect_that(root.state_test, equals(root.state_correct))
  
  # Dimension p
  p <- 5
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = rep(3, p),
                          exp.root = rep(2, p),
                          var.root = as.matrix(5, p, p))
  optimal.value <- rep(2, p) 
  selection.strength <- matrix(3, 1, 1)
  variance <- as.matrix(5, p, p)
  var.root <- compute_stationary_variance(variance, selection.strength)
  root.state_correct <- list(random = TRUE,
                             stationary.root = TRUE, 
                             value.root = NA,
                             exp.root = rep(2, p),
                             var.root = as(var.root, "dpoMatrix"))
  
  expect_warning(root.state_test <- test.root.state(root.state_test, "OU",
                                                    optimal.value = optimal.value,
                                                    variance = variance,
                                                    selection.strength = selection.strength))
  expect_that(root.state_test, equals(root.state_correct))
  
  # optimal value
  p <- 5
  root.state_test <- list(random = TRUE,
                          stationary.root = TRUE, 
                          value.root = rep(3, p),
                          exp.root = rep(2, p),
                          var.root = as.matrix(5, p, p))
  optimal.value <- rep(5, p) 
  selection.strength <- matrix(3, 1, 1)
  variance <- compute_variance_from_stationary(root.state_test$var.root, selection.strength)
  root.state_correct <- list(random = TRUE,
                             stationary.root = TRUE, 
                             value.root = NA,
                             exp.root = rep(5, p),
                             var.root = as(as.matrix(5, p, p), "dpoMatrix"))
  
  expect_warning(root.state_test <- test.root.state(root.state_test, "OU",
                                                    optimal.value = optimal.value,
                                                    variance = variance,
                                                    selection.strength = selection.strength))
  expect_that(root.state_test, equals(root.state_correct))
})

test_that("check data",{
  p <- 4
  ntaxa <- 236
  
  set.seed(1958)
  tree <- rtree(ntaxa)
  
  ## Missing names
  Y_data <- matrix(1, p, ntaxa)
  expect_that(check_data(tree, Y_data, TRUE), gives_warning())
  
  ## Wrong dimensions
  Y_data <- matrix(1, ntaxa, p)
  expect_that(check_data(tree, Y_data, TRUE), throws_error())
  
  ## Reordering
  Y_data <- matrix(1, p, ntaxa)
  colnames(Y_data) <- tree$tip.label
  expect_that(check_data(tree, Y_data, TRUE), equals(Y_data))
  
  colnames(Y_data) <- sample(tree$tip.label, ntaxa)
  expect_that(data_new <- check_data(tree, Y_data, TRUE), gives_warning())
  expect_that(data_new, equals(Y_data[ , tree$tip.label]))
  expect_that(check_data(tree, Y_data, FALSE), equals(Y_data))
})

test_that("check tree",{
  set.seed(1958)
  ntaxa <- 20
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                                   age = 1, mrca = TRUE)[[1]]
  Y_data <- matrix(1, 2, ntaxa)
  
  ## Not ultrametric
  tree$edge.length[2 * ntaxa - 2] <- tree$edge.length[2 * ntaxa - 2] * 3/2
  expect_that(PhyloEM(phylo = tree, Y_data = Y_data, process = "scOU"), throws_error("The tree must be ultrametric."))
  tree$edge.length[2 * ntaxa - 2] <- tree$edge.length[2 * ntaxa - 2] * 2/3
  
  ## Zero length branch
  ll <- tree$edge.length[2]
  tree$edge.length[2] <- 0
  tmp <- extract.clade(tree, tree$edge[2, 2])
  tips <- match(tmp$tip.label, tree$tip.label)
  edges <- match(tips, tree$edge[, 2])
  tree$edge.length[edges] <- tree$edge.length[edges] + ll

  expect_that(PhyloEM(phylo = tree, Y_data = Y_data), throws_error("The tree has zero-length branches."))
})