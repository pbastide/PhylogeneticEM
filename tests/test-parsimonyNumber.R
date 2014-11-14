context("Number of Parsimonious Solutions")

test_that("parsimonyNumber and enumerate_parsimony work on simple exemples", {
  require(ape)
  require(plyr)
  ## First exemple (tricky one)
  tree <- read.tree(text="(((T,T),C),C);")
  clusters=c(1,2,3,3)
  Nbr <- extract.parcimonyNumber(parcimonyNumber(tree, clusters))
  Allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
  # Same number of solutions
  expect_that(Nbr, equals(dim(Allocs)[1]))
  expect_that(Nbr, equals(3))
})