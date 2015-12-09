rm(list=ls())

library(ape)
library(LINselect)
library(Rphylopars)
library(glmnet)
# library(plyr)
# library(quadrupen)
# library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
# library(TreeSim)
library(doParallel)
library(foreach)

## Load Ericaceae data
load("../data/ericaceae_data.RData")

source("R/simulate.R")
source("R/estimateEM.R")
source("R/init_EM.R")
source("R/E_step.R")
source("R/M_step.R")
source("R/shutoff.R")
source("R/generic_functions.R")
source("R/shifts_manipulations.R")
source("R/plot_functions.R")
source("R/parsimonyNumber.R")
source("R/partitionsNumber.R")
source("R/model_selection.R")

datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

###############################################################################
## EM Univariate ##############################################################
###############################################################################
## Select data
phylo <- subtree_traits[[1]]
trait_matrix <- trait_matrix_all[[1]]
# 0 values
set.seed(17910402)
temp_cl <- replace_zeros(trait_matrix[1, ])
trait_matrix[1, ] <- temp_cl$x
# Log transform sizes
trait_matrix_transform <- trait_matrix
trait_matrix_transform[1:2, ] <- log(trait_matrix_transform[1:2, ])

## Estimate one trait
estimate_univariate <- function(i){
  # Drop NAs
  data_uni <- trait_matrix_transform[i, !is.na(trait_matrix_transform[i, ])]
  Overlap_traits <- name.check(phylo, data_uni)
  tree_uni <- drop.tip(phylo, Overlap_traits$Tree.not.data)
  tree_char <- tree_uni
  
  ## Re-scale tree to one
  height_tree <- node.depth.edgelength(tree_uni)[1]
  tree_uni$edge.length <- tree_uni$edge.length / height_tree
  
  # Alpha Values
  alpha_grid <- find_grid_alpha(tree_uni,
                                nbr_alpha = 10,
                                factor_up_alpha = 2,
                                factor_down_alpha = 4,
                                quantile_low_distance = 0.001,
                                log_transform = TRUE)
  
  res <- PhyloEM(phylo = tree_uni,
                 Y_data = data_uni,
                 process = "scOU",
                 random.root = FALSE,
                 K_max = 35,
                 alpha_known = TRUE,
                 alpha = alpha_grid,
                 tol = list(variance = 10^(-2), 
                            value.root = 10^(-2),
                            log_likelihood = 10^(-2)),
                 use_previous = FALSE,
                 method.init = "lasso",
                 method.selection = "BGH")
  
  res$tree <- tree_char
  res$alpha_grid <- alpha_grid
  
  save.image(file = paste0("../Results/Test_Cases/trait_scOU_fixed_root_trait_", i, ".RData"))
  
  return(res)
}

## Setting for parallel computations
reqpckg <- c("ape", "glmnet", "robustbase", "Rphylopars", "LINselect")
Ncores <- 2
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
res_uni_all <- foreach(i = 1:6, .packages = reqpckg) %dopar%
{
  estimate_univariate(i)
}
# Stop the cluster (parallel)
stopCluster(cl)

save.image(file = paste0("../Results/Test_Cases/trait_scOU_fixed_root_trait_all.RData"))