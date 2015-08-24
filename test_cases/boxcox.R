##################################################
## BoxCox transform
##################################################

##########
## Chelonia EM

rm(list=ls())

PATH <- "../Results/Chelonia/"

load(paste0(PATH, "chelonia_OUwie", ".RData"))

# Dependencies of EM
library(ape)
library(quadrupen)
library(robustbase)
library(plyr)
# Tree Sim
library(TreeSim)
# Model selection
library(LINselect)
# Plot
library(reshape2)
library(ggplot2)

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

## Import data set form package bayou
library(bayou)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

## Fixed quantities
times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)
subtree.list <- enumerate_tips_under_edges(tree)
T_tree <- incidence.matrix(tree)

#############################################
## Select some regression variables
#############################################
## Shifts to ensure habitats
shifts_edges <- allocate_shifts_from_regimes(tree, ParSols[1,])

# Correction for independence
Sigma <- times_shared + 0.01 ## BM
#Sigma <- exp(-0.01*distances_phylo) ## OU
Sigma_YY <- extract.variance_covariance(Sigma, what = "YY")

L <- chol(Sigma_YY)
Linv <- solve(L)

############################################
## Box Cox
############################################
library(MASS)
# Corrected
X <- Linv %*% T_tree[, shifts_edges] + 0
Y <- Linv %*% exp(data)
mu <- - min(Y) + 1
fit <- lm(Y + mu ~ X)
bc <- boxcox(fit)

# Not corrected
X <- T_tree[, shifts_edges] + 0
Y <- exp(data)
fit <- lm(Y ~ X)
bc <- boxcox(fit)
