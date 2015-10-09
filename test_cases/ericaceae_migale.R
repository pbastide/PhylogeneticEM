rm(list=ls())

library(ape)
library(plyr)
library(quadrupen)
library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(TreeSim)

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

## Load Ericaceae data
load("../data/ericaceae_data.RData")

datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

## Get a single data frame ################################################

## Select Columns: Corrola lenght (2), anther length (7)
linear_measures <- data_reduced[, c(1, 2, 7)]
## Compute species mean (one value per species)
mean_linear_measures_length <- aggregate(cbind(corolla_tube_length, total_large_anther_length) ~ species,
                                         data = linear_measures,
                                         mean,
                                         na.action = na.pass)
## Add PC1 and PC2 from Corolla and Anther Shape
# anther
# pca_anthers <- data.frame(pca_efou_anthers$x[, 1:2])
pca_anthers <- data.frame(all_data_pcs_anthers_unique[, 1:3])
colnames(pca_anthers) <- c("species", "anther.PC1", "anther.PC2")
# pca_anthers$species <- rownames(pca_anthers)
# Corolla
# pca_corolla <- data.frame(pca_efou_flowers$x[, 1:2])
pca_corolla <- data.frame(all_data_pcs_flowers_unique[, 1:3])
colnames(pca_corolla) <- c("species", "corolla.PC1", "corolla.PC2")
# pca_corolla$species <- rownames(pca_corolla)
# Merge
dd <- merge(mean_linear_measures_length, pca_anthers, by = "species", all = TRUE) # Match: 4 !!
traits_data <- merge(dd, pca_corolla, by = "species", all = TRUE)
## Drop NA (entire ligne is NA)
traits_data <- traits_data[rowSums(is.na(traits_data)) < (dim(traits_data)[2]-1),]
rownames(traits_data) <- traits_data[,1]
percentage_missing <- sum(is.na(traits_data))/prod(dim(traits_data)) * 100
## Check compatibility with the tree
Overlap_traits <- name.check(tree, traits_data)
subtree_traits <- drop.tip(tree, Overlap_traits$Tree.not.data)
# Sort data in same order as tree
match_traits <- match(subtree_traits$tip.label, rownames(traits_data))
sortedData_traits <- traits_data[match_traits, -1]
trait_matrix <- t(as.matrix(sortedData_traits))

## EM inferences - without model selection #######################################
##################################################################################
## Fixed quantities
times_shared <- compute_times_ca(phylo)
distances_phylo <- compute_dist_phy(phylo)
subtree.list <- enumerate_tips_under_edges(phylo)
T_tree <- incidence.matrix(phylo)
h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])

## Inference
set.seed(17920920)
X <- list(Y_data = Y_data,
          K_try = 0:K_max,
          ntaxa = ntaxa)
## Estimations
X <- Phylo_EM_sequencial(phylo = subtree_traits,
                         Y_data = trait_matrix,
                         process = "BM",
                         K_max = 10,
                         curent = X,
                         random.root = FALSE,
                         alp = 0,
                         progress.bar = FALSE,
                         times_shared = times_shared,
                         distances_phylo = distances_phylo,
                         subtree.list = subtree.list,
                         T_tree = T_tree,
                         h_tree = h_tree,
                         save_step = FALSE,
                         tol = list(variance = 10^(-2), 
                                    value.root = 10^(-2)),
                         Nbr_It_Max = 100)

save.image(file = paste0("../Results/Test_Cases/ericaceae_migale_", datestamp, ".RData"))