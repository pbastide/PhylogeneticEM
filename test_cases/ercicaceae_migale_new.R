rm(list=ls())

library(doParallel)
library(foreach)
library(PhylogeneticEM)
reqpckg <- c("PhylogeneticEM")

## Load Ericaceae data
load("../data/ericaceae_data_2016-09-20.RData")

datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

phylo <- comptree_all
trait_matrix <- t(sortedData_all[,c(2:7)])
colnames(trait_matrix) <- sortedData_all[,1]

## Alpha values
# Re-scale tree
height_tree <- node.depth.edgelength(phylo)[1]
phylo$edge.length <- phylo$edge.length / height_tree

alpha_grid <- find_grid_alpha(phylo,
                              nbr_alpha = 10,
                              factor_up_alpha = 2,
                              factor_down_alpha = 4,
                              quantile_low_distance = 0.0005,
                              log_transform = TRUE)

## Inference (no model selection)
# Root Stationary
# Lasso init
# K_max = 35
res <- PhyloEM(phylo = phylo,
               Y_data = trait_matrix,
               process = "scOU",
               K_max = 40,
               random.root = TRUE,
               stationary.root = TRUE,
               alpha = alpha_grid[2],
               save_step = FALSE,
               Nbr_It_Max = 2000,
               method.variance = "upward_downward",
               method.init = "lasso",
               use_previous = FALSE,
               method.selection = c("BirgeMassart1", "BirgeMassart2"),
               K_lag_init = 5)

save.image(file = paste0("../Results/Test_Cases/ericaceae_all_", datestamp, ".RData"))