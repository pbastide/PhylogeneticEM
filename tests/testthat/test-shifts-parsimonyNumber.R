context("Number of Parsimonious Solutions")

test_that("parsimonyNumber and enumerate_parsimony on simple exemples", {
  ## First exemple (tricky one)
  tree <- read.tree(text="(((T,T),C),C);")
  clusters <- c(3.21, 26, 5.3, 5.3)
  Nbr <- extract(parsimonyNumber(tree, clusters))
  Allocs <- extract(enumerate_parsimony(tree,clusters))
  # Same number of solutions
  expect_that(Nbr, equals(dim(Allocs)[1]))
  # Right number of solutions
  expect_that(Nbr, equals(3))
  # Right solutions
  right_allocs <- rbind(c(3.21,26,5.3,5.3,5.3,5.3,3.21),
                        c(3.21,26,5.3,5.3,5.3,5.3,26),
                        c(3.21,26,5.3,5.3,5.3,5.3,5.3))
  expect_that(Allocs, equals(right_allocs))
  
  ## Second exemple
  tree <- read.tree(text="(((T,T),C),(A,A));")
  clusters=c(1,1,2,3,3)
  Nbr <- extract(parsimonyNumber(tree, clusters))
  Allocs <- extract(enumerate_parsimony(tree,clusters))
  expect_that(Nbr, equals(dim(Allocs)[1]))
  expect_that(Nbr, equals(5))
  
  ## Third exemple
  tree <- read.tree(text="(((T,T),C),(A,A));")
  clusters=c(1,1,1,1,1)
  Nbr <- extract(parsimonyNumber(tree, clusters))
  Allocs <- extract(enumerate_parsimony(tree,clusters))
  expect_that(Nbr, equals(dim(Allocs)[1]))
  expect_that(Nbr, equals(1))
  right_allocs <- matrix(rep(1, 9), nrow = 1)
  expect_that(Allocs, equals(right_allocs))
})

test_that("Equivalent Shift - BM", {
  set.seed(17920902)
  ntaxa = 20
  phylo <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                                   age = 1, mrca = TRUE)[[1]]
  ## Fixed Root
  params <- params_BM(p = 4, 
                      edges = c(4, 17, 22),
                      values = cbind(c(1, 2, 3, 4),
                                     c(-1, -2, -3, -4),
                                     c(1, 1, 1, 1)))
  ed_shifts <- equivalent_shifts(phylo, params)
  # plot(ed_shifts)
  ## Random Root
  params_rand <- params_BM(p = 4, 
                           edges = c(4, 17, 22), 
                           values = cbind(c(1, 2, 3, 4),
                                          c(-1, -2, -3, -4),
                                          c(1, 1, 1, 1)),
                           random = TRUE,
                           exp.root = rep(0, 4))
  ed_shifts_rand <- equivalent_shifts(phylo, params_rand)
  
  expect_equal(ed_shifts, ed_shifts_rand)
  expect_equal(extract(ed_shifts), extract(ed_shifts_rand))
  expect_equal(extract(ed_shifts, what = "root"), extract(ed_shifts_rand, what = "root"))
  
  ## Check with enumerate parsimony
  clusters <- allocate_regimes_from_shifts(phylo, params$shifts$edges)
  Allocs <- extract(enumerate_parsimony(phylo, clusters))
  eq_shifts_bis <- apply(Allocs, 1, function(z) allocate_shifts_from_regimes(phylo, z))
  expect_equal(ed_shifts$eq_shifts_edges, eq_shifts_bis)
})

test_that("Equivalent Shift - OU", {
  set.seed(17920902)
  ntaxa = 20
  phylo <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                                    age = 1, mrca = TRUE)[[1]]
  ## Fixed Root
  params <- params_OU(p = 4, 
                      edges = c(4, 17, 22),
                      values = cbind(c(1, 2, 3, 4), c(-1, -2, -3, -4), c(1, 1, 1, 1)))
  ed_shifts <- equivalent_shifts(phylo, params)
  # plot(ed_shifts)
  ## Random Root
  params_rand <- params_OU(p = 4, 
                           edges = c(4, 17, 22), 
                           values = cbind(c(1, 2, 3, 4), c(-1, -2, -3, -4), c(1, 1, 1, 1)),
                           random = TRUE, exp.root = rep(0, 4))
  ed_shifts_rand <- equivalent_shifts(phylo, params_rand)
  
  expect_equal(ed_shifts, ed_shifts_rand)
  
  ## Check with enumerate parsimony
  clusters <- allocate_regimes_from_shifts(phylo, params$shifts$edges)
  Allocs <- extract(enumerate_parsimony(phylo, clusters))
  eq_shifts_bis <- apply(Allocs, 1, function(z) allocate_shifts_from_regimes(phylo, z))
  expect_equal(ed_shifts$eq_shifts_edges, eq_shifts_bis)
})

# test_that("parsimonyNumber and enumerate_parsimony on random exemples", {
#   require(ape)
#   require(plyr)
#   ## First exemple (one cluster)
#   set.seed(141114)
#   tree <- rtree(100)
#   clusters <- sample(1, 100, replace = TRUE)
#   Nbr <- extract(parsimonyNumber(tree, clusters))
#   Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
#   expect_that(Nbr, equals(dim(Allocs)[1]))
#   expect_that(Nbr, equals(1))
#   
#   ## Second exemple (3 clusters)
#   set.seed(131114)
#   tree <- rtree(50)
#   clusters <- sample(3, 50, replace = TRUE)
#   Nbr <- extract(parsimonyNumber(tree, clusters))
#   Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
#   expect_that(Nbr, equals(dim(Allocs)[1]))
# })