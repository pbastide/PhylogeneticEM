context("Number Tree Compatible Paritions")

test_that("partitionsNumber in Binary Case", {
  p <- 6
  n <- 2^p
  k <- 30
  N <- 2:n
  
  ## Exact values
  val <- sapply(choose(2*N-1-k, k-1), function(z) max(0, z))
  val_mark <- sapply(choose(2*N-k, k-1), function(z) max(0, z))
  
  ## Comb Tree
  CombTree <- rtree.comb(n)
  val_comb <- extract(partitionsNumber(CombTree, k), node = rev((n+1):(2*n-1)))
  val_comb_mark <- extract(partitionsNumber(CombTree, k), node = rev((n+1):(2*n-1)), marked = TRUE)
  # Not marked
  expect_that(val, equals(val_comb))
  # Marked
  expect_that(val_mark, equals(val_comb_mark))
  
  ## Symetric tree
  SymTree <- rtree.sym(p)
  val_sym <- extract(partitionsNumber(SymTree, k))
  val_sym_mark <- extract(partitionsNumber(SymTree, k), marked = TRUE)
  # Not marked
  expect_that(val[n-1], equals(val_sym))
  # Marked
  expect_that(val_mark[n-1], equals(val_sym_mark))
  
  ## Random Tree
  randomTree <- rtree(n)
  val_rand <- extract(partitionsNumber(randomTree, k))
  val_rand_mark <- extract(partitionsNumber(randomTree, k), marked = TRUE)
  # Not marked
  expect_that(val[n-1], equals(val_rand))
  # Marked
  expect_that(val_mark[n-1], equals(val_rand_mark))
})

# test_that("break point in binary case", {
#   N <- 2:500
#   
#   ## Exact values
#   fun <- function(n){
#     kk <- 1:n
#     cc <- choose(2*n - 2 - kk, kk)
#     max <- which(cc == max(cc))
#     return(max[length(max)])
#   }
#   val <- sapply(N, fun)
#   
#   ## Break Points
#   bb <- sapply(N, complexity_break_point)
#   
#   expect_that(val, equals(bb))
# })

test_that("partitionsNumber in General Case", {
  testthat::skip_if_not_installed("combinat")
  K <- 5
  tree <- read.tree(text = "(A,(A,A,A));")
  xx <- partitionsNumber(tree, K)
  # Not marked
  NN <- extract.partitionsNumber(xx, npart = 1:K)
  expect_equal(NN, c(1, 4, 6, 1, 0))
  # Marked
  MM <- extract.partitionsNumber(xx, npart = 1:K, marked = TRUE)
  expect_equal(MM, c(1, 5, 9, 4, 0))
})