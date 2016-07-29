rm(list=ls())

library(doParallel)
library(foreach)
library(ape)
library(glmnet) # For Lasso initialization
library(robustbase) # For robust fitting of alpha
library(gglasso)
library(capushe)
library(Matrix)
reqpckg <- c("ape", "glmnet", "robustbase", "gglasso", "Matrix", "capushe")

## Load Ericaceae data
load("../data/ericaceae_data_2016-07-27.RData")

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


## Get a single data frame ####################################################
## Select Columns: Corrola lenght, anther length
linear_measures <- data_reduced[, c("species",
                                    "corolla_tube_length",
                                    "total_large_anther_length")]
## Compute species mean (one value per species)
mean_linear_measures_length <- aggregate(cbind(corolla_tube_length, total_large_anther_length) ~ species,
                                         data = linear_measures,
                                         mean,
                                         na.action = na.pass)
## Add PC1 and PC2 from Corolla and Anther Shape
# anther
pca_anthers <- data.frame(all_data_pcs_anthers_unique[, 1:3])
colnames(pca_anthers) <- c("species", "anther.PC1", "anther.PC2")
# Corolla
pca_corolla <- data.frame(all_data_pcs_flowers_unique[, 1:3])
colnames(pca_corolla) <- c("species", "corolla.PC1", "corolla.PC2")
# Merge
dd <- merge(mean_linear_measures_length, pca_anthers,
            by = "species", all = TRUE)
traits_data <- merge(dd, pca_corolla, by = "species", all = TRUE)
## Drop NA (entire ligne is NA)
traits_data_all <- vector("list", 6)
trait_matrix_all <- vector("list", 6)
subtree_traits <- vector("list", 6)
percentage_missing <- vector(length = 6)
for (i in 1:6){ # i is the minimal number of non-NA trait
  traits_data_all[[i]] <- traits_data[rowSums(is.na(traits_data)) < (dim(traits_data)[2] - i),]
  rownames(traits_data_all[[i]]) <- traits_data_all[[i]][,1]
  ## Check compatibility with the tree
  Overlap_traits <- name.check(tree, traits_data_all[[i]])
  subtree_traits[[i]] <- drop.tip(tree, Overlap_traits$Tree.not.data)
  # Sort data in same order as tree
  match_traits <- match(subtree_traits[[i]]$tip.label,
                        rownames(traits_data_all[[i]]))
  sortedData_traits <- traits_data_all[[i]][match_traits, -1]
  trait_matrix_all[[i]] <- t(as.matrix(sortedData_traits))
  percentage_missing[i] <- sum(is.na(trait_matrix_all[[i]]))/prod(dim(trait_matrix_all[[i]])) * 100
}

replace_zeros <- function(x){
  zeros_x <- (!is.na(x) & x == 0)
  m <- min(x[x > 0], na.rm = TRUE) / 2
  sd <- m / 2.33
  values <- rnorm(sum(zeros_x), m, sd)
  counter <- 1
  while(counter < 100 && any(values <= 0)){
    values <- rnorm(sum(zeros_x), m, sd)
    counter <- counter + 1
  }
  if (counter < 100){
    x[zeros_x] <- rnorm(sum(zeros_x), m, sd)
    return(list(x = x,
                zeros_x = zeros_x,
                m = m,
                sd = sd))
  } else {
    stop("Please consider a smaller variance, as I could not find only positive values.")
  }
}

###############################################################################
## EM inferences ######################### ####################################
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
                              quantile_low_distance = 0.0005,
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
               K_lag_init = 5,
               alpha = alpha_grid[11],
               tol = list(variance = 10^(-2), 
                          value.root = 10^(-2),
                          log_likelihood = 10^(-2)),
               save_step = FALSE,
               Nbr_It_Max = 2000,
               use_previous = FALSE,
               method.init = "lasso",
               method.selection = c("BirgeMassart1", "BirgeMassart2"),
               impute_init_Rphylopars = FALSE)
               # parallel_alpha = TRUE, Ncores = 5,
               # exportFunctions = exportFunctions)

save.image(file = paste0("../Results/Test_Cases/ericaceae_1_NA_", datestamp, ".RData"))
