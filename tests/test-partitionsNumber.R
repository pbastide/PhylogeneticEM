context("Number Tree Compatible Paritions")

test_that("partitionsNumber in Binary Case", {
  require(ape)
  p <- 6
  n <- 2^p
  k <- 30
  N <- 2:n
  
  ## Exact values
  val <- sapply(choose(2*N-1-k, k-1), function(z) max(0, z))
  val_mark <- sapply(choose(2*N-k, k-1), function(z) max(0, z))
  
  ## Comb Tree
  CombTree <- rtree.comb(n)
  val_comb <- extract.partitionsNumber(partitionsNumber(CombTree, k), node = rev((n+1):(2*n-1)))
  val_comb_mark <- extract.partitionsNumber(partitionsNumber(CombTree, k), node = rev((n+1):(2*n-1)), marqued = TRUE)
  # Not Marqued
  expect_that(val, equals(val_comb))
  # Marked
  expect_that(val_mark, equals(val_comb_mark))
  
  ## Symetric tree
  SymTree <- rtree.sym(p)
  val_sym <- extract.partitionsNumber(partitionsNumber(SymTree, k))
  val_sym_mark <- extract.partitionsNumber(partitionsNumber(SymTree, k), marqued = TRUE)
  # Not Marqued
  expect_that(val[n-1], equals(val_sym))
  # Marked
  expect_that(val_mark[n-1], equals(val_sym_mark))
  
  ## Random Tree
  randomTree <- rtree(n)
  val_rand <- extract.partitionsNumber(partitionsNumber(randomTree, k))
  val_rand_mark <- extract.partitionsNumber(partitionsNumber(randomTree, k), marqued = TRUE)
  # Not Marqued
  expect_that(val[n-1], equals(val_rand))
  # Marked
  expect_that(val_mark[n-1], equals(val_rand_mark))
})