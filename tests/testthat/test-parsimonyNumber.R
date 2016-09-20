context("Number of Parsimonious Solutions")

test_that("parsimonyNumber and enumerate_parsimony on simple exemples", {
  ## First exemple (tricky one)
  tree <- read.tree(text="(((T,T),C),C);")
  clusters <- c(3.21, 26, 5.3, 5.3)
  Nbr <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
  Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
  # Same number of solutions
  expect_that(Nbr, equals(dim(Allocs)[1]))
  # Right number of solutions
  expect_that(Nbr, equals(3))
  # Right solutions
  right_allocs <- rbind(c(3.21,26,5.3,5.3,5.3,5.3,3.21), c(3.21,26,5.3,5.3,5.3,5.3,26), c(3.21,26,5.3,5.3,5.3,5.3,5.3))
  expect_that(Allocs, equals(right_allocs))
  
  ## Second exemple
  tree <- read.tree(text="(((T,T),C),(A,A));")
  clusters=c(1,1,2,3,3)
  Nbr <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
  Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
  expect_that(Nbr, equals(dim(Allocs)[1]))
  expect_that(Nbr, equals(5))
  
  ## Third exemple
  tree <- read.tree(text="(((T,T),C),(A,A));")
  clusters=c(1,1,1,1,1)
  Nbr <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
  Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
  expect_that(Nbr, equals(dim(Allocs)[1]))
  expect_that(Nbr, equals(1))
  right_allocs <- matrix(rep(1, 9), nrow = 1)
  expect_that(Allocs, equals(right_allocs))
})

# test_that("parsimonyNumber and enumerate_parsimony on random exemples", {
#   require(ape)
#   require(plyr)
#   ## First exemple (one cluster)
#   set.seed(141114)
#   tree <- rtree(100)
#   clusters <- sample(1, 100, replace = TRUE)
#   Nbr <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
#   Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
#   expect_that(Nbr, equals(dim(Allocs)[1]))
#   expect_that(Nbr, equals(1))
#   
#   ## Second exemple (3 clusters)
#   set.seed(131114)
#   tree <- rtree(50)
#   clusters <- sample(3, 50, replace = TRUE)
#   Nbr <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
#   Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
#   expect_that(Nbr, equals(dim(Allocs)[1]))
# })