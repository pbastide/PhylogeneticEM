rm(list=ls())

library(doParallel)
library(foreach)
library(PhylogeneticEM)
reqpckg <- c("PhylogeneticEM")

## Load Ericaceae data
load("../data/ericaceae_data_2016-10-13.RData")

datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

## Prepare Data
phylo <- comptree_all
trait_matrix <- t(sortedData_all[,c(2:7)])
colnames(trait_matrix) <- sortedData_all[,1]

## Re-scale tree to unit length
height_tree <- node.depth.edgelength(phylo)[1]
phylo$edge.length <- phylo$edge.length / height_tree

## Find grid on alpha
alpha_grid <- find_grid_alpha(phylo,
                              nbr_alpha = 10,
                              factor_up_alpha = 2,
                              factor_down_alpha = 4,
                              quantile_low_distance = 0.0005,
                              log_transform = TRUE)

## Run algorithm (about 10 hours on 5 cores)
res <- PhyloEM(phylo = phylo,
               Y_data = trait_matrix,
               process = "scOU",            ## scalar OU
               random.root = TRUE,          ## Root is stationary
               stationary.root = TRUE,
               alpha = alpha_grid[-1],      ## On a grid of alpha (non 0)
               K_max = 40,                  ## Maximal number of shifts
               K_lag_init = 5,              ## Tuning for lasso init
               Nbr_It_Max = 2000,           ## Max number of EM steps
               parallel_alpha = TRUE,       ## Parallelize on alpha values
               Ncores = 5)

save.image(file = paste0("../Results/Test_Cases/ericaceae_all_", datestamp, ".RData"))