##########
## Test Case : Chelonian carapace evolution.
## Used in :
## - Uyeda 2014
## - Eastman 2011 (link to data : https://github.com/eastman/auteur/tree/master/auteur/data)
## -> Jaffe 2011 (original data)

rm(list=ls())
WD_mac <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
WD_unb <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
WD <- ifelse(file.exists(WD_mac), WD_mac, WD_unb)
setwd(WD)

PATH <- paste(WD, "/Results/Chelonia/", sep="")
library(ape)
library(glmnet) # For Lasso initialization
library(quadrupen)
#library(nlme) # For second derivative computation
#library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
library(TreeSim)
library(plyr)
source("Phylogenetic-EM/simulate.R")
source("Phylogenetic-EM/estimateEM.R")
source("Phylogenetic-EM/init_EM.R")
source("Phylogenetic-EM/E_step.R")
source("Phylogenetic-EM/M_step.R")
source("Phylogenetic-EM/shutoff.R")
source("Phylogenetic-EM/generic_functions.R")
source("Phylogenetic-EM/shifts_manipulations.R")
source("Phylogenetic-EM/plot_functions.R")
source("Phylogenetic-EM/parsimonyNumber.R")

load("data/chelonia.rda")

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

plot(tree, show.tip.label = FALSE)

# EM
K <- 16
time <- system.time(
  results_estim_EM <- estimateEM(phylo = tree, 
                                 Y_data = data, 
                                 tol = list(variance=10^(-4), 
                                            value.root=10^(-4), 
                                            exp.root=10^(-4), 
                                            var.root=10^(-4), 
                                            selection.strength=10^(-3)), 
                                 process = "OU", 
                                 method.variance = "simple", 
                                 method.init = "lasso",
                                 method.init.alpha = "estimation",
                                 Nbr_It_Max = 1000, 
                                 nbr_of_shifts = K, 
                                 alpha_known = FALSE,
                                 min_params=list(variance = 10^(-4), 
                                                 value.root = -10^(4), 
                                                 exp.root = -10^(4), 
                                                 var.root = 10^(-4),
                                                 selection.strength = 10^(-4)),
                                 max_params=list(variance = 10^(4), 
                                                 value.root = 10^(4), 
                                                 exp.root = 10^(4), 
                                                 var.root = 10^(4),
                                                 selection.strength = 10^(4)),
                                 methods.segmentation = c("same_shifts",
                                                          "best_single_move",
                                                          "lasso"))
)

name <- paste0("_K=", K)
history <- list_to_table.history(results_estim_EM$params_history)
write.csv2(history, paste0(PATH, "boite_noire", name, ".csv"))
plot.history.OU.stationnary(params_history = results_estim_EM$params_history,
                            tree = tree,
                            Y_data_ref = data,
                            PATH = PATH,
                            name = paste0("history_plot", name))
plot.process.actual(Y.state = data,
                    Z.state = results_estim_EM$ReconstructedNodesStates,
                    phylo = tree, 
                    paramsEstimate = results_estim_EM$params)
save.image(paste0(PATH, "estimation", name, ".RData"))
