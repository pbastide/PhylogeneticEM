if (requireNamespace("testthat", quietly = TRUE)) {
  
  library(testthat)
  library(PhylogeneticEM)
  
  test_check("PhylogeneticEM", filter = "upward")
  
}