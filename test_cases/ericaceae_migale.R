rm(list=ls())

library(ape)
library(plyr)
# library(quadrupen)
library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
# library(TreeSim)
library(Matrix)


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

exportFunctions <- ls() # All the functions for parallel computations.

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


# ## Get a single data frame ####################################################
# 
# ## Select Columns: Corrola lenght (2), anther length (7)
# linear_measures <- data_reduced[, c(1, 2, 7)]
# ## Compute species mean (one value per species)
# mean_linear_measures_length <- aggregate(cbind(corolla_tube_length, total_large_anther_length) ~ species,
#                                          data = linear_measures,
#                                          mean,
#                                          na.action = na.pass)
# ## Add PC1 and PC2 from Corolla and Anther Shape
# # anther
# # pca_anthers <- data.frame(pca_efou_anthers$x[, 1:2])
# pca_anthers <- data.frame(all_data_pcs_anthers_unique[, 1:3])
# colnames(pca_anthers) <- c("species", "anther.PC1", "anther.PC2")
# # pca_anthers$species <- rownames(pca_anthers)
# # Corolla
# # pca_corolla <- data.frame(pca_efou_flowers$x[, 1:2])
# pca_corolla <- data.frame(all_data_pcs_flowers_unique[, 1:3])
# colnames(pca_corolla) <- c("species", "corolla.PC1", "corolla.PC2")
# # pca_corolla$species <- rownames(pca_corolla)
# # Merge
# dd <- merge(mean_linear_measures_length, pca_anthers, by = "species", all = TRUE) # Match: 4 !!
# traits_data <- merge(dd, pca_corolla, by = "species", all = TRUE)
# ## Drop NA (entire ligne is NA)
# traits_data <- traits_data[rowSums(is.na(traits_data)) < (dim(traits_data)[2]-1),]
# rownames(traits_data) <- traits_data[,1]
# percentage_missing <- sum(is.na(traits_data))/prod(dim(traits_data)) * 100
# ## Check compatibility with the tree
# Overlap_traits <- name.check(tree, traits_data)
# subtree_traits <- drop.tip(tree, Overlap_traits$Tree.not.data)
# # Sort data in same order as tree
# match_traits <- match(subtree_traits$tip.label, rownames(traits_data))
# sortedData_traits <- traits_data[match_traits, -1]
# trait_matrix <- t(as.matrix(sortedData_traits))


###############################################################################
## EM inferences - without model selection ####################################
###############################################################################

## Select data
phylo <- subtree_traits[[6]]
trait_matrix <- trait_matrix_all[[6]]
# 0 values
set.seed(17910402)
temp_cl <- replace_zeros(trait_matrix[1, ])
trait_matrix[1, ] <- temp_cl$x
# Log transform sizes
trait_matrix_transform <- trait_matrix
trait_matrix_transform[1:2, ] <- log(trait_matrix_transform[1:2, ])

## Alpha values
# Re-scale tree
height_tree <- node.depth.edgelength(phylo)[1]
phylo$edge.length <- phylo$edge.length / height_tree

alpha_grid <- find_grid_alpha(phylo,
                              nbr_alpha = 10,
                              factor_up_alpha = 2,
                              factor_down_alpha = 4,
                              quantile_low_distance = 0.0006,
                              log_transform = TRUE)

## Inference (no model selection)
# Root fixed
# Lasso init
# K_max = 35
res <- PhyloEM(phylo = phylo,
               Y_data = trait_matrix_transform,
               process = "scOU",
               random.root = FALSE,
               K_max = 30,
               alpha_known = TRUE,
               alpha = alpha_grid,
               tol = list(variance = 10^(-2), 
                          value.root = 10^(-2),
                          log_likelihood = 10^(-2)),
               use_previous = FALSE,
               method.init = "lasso",
               method.selection = c("BirgeMassart1", "BirgeMassart2"))
               # parallel_alpha = TRUE, Ncores = 5,
               # exportFunctions = exportFunctions)

save.image(file = paste0("../Results/Test_Cases/ericaceae_migale_", datestamp, ".RData"))