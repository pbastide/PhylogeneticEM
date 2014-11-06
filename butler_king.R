#######################################
## Test of EM - Butler King
#######################################
rm(list=ls())
WD <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
# WD <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
setwd(WD)
PATH <- paste(WD, "/Results/Butler_King/", sep="")
library(ape)
#library(glmnet) # For Lasso initialization
library(quadrupen)
#library(nlme) # For second derivative computation
#library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
source("Phylogenetic-EM/simulate.R")
source("Phylogenetic-EM/estimateEM.R")
source("Phylogenetic-EM/init_EM.R")
source("Phylogenetic-EM/E_step.R")
source("Phylogenetic-EM/M_step.R")
source("Phylogenetic-EM/shutoff.R")
source("Phylogenetic-EM/generic_functions.R")
source("Phylogenetic-EM/shifts_manipulations.R")
source("Phylogenetic-EM/plot_functions.R")

library(ouch)
source("R/convert.R")
# Anolis bimaculatus lizard size data
data(bimac)
bimac_ouch_tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
bimac_ape_tree <- convert.pmc(bimac_ouch_tree)
plot(bimac_ape_tree); edgelabels()
Y_data_bimac <- bimac$size[23:45]
names(Y_data_bimac) <- bimac$species[23:45]
Y_data_bimac <- log(unname(Y_data_bimac[bimac_ape_tree$tip.label]))

# EM
K <- 4
time <- system.time(
  results_estim_EM <- estimateEM(phylo = bimac_ape_tree, 
                                       Y_data = Y_data_bimac, 
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

# Infered parameters by Butler and King
root.stateBK <- list(random = TRUE, 
                     stationary.root = TRUE, 
                     value.root = NA,
                     exp.root = 0.86, 
                     var.root = (0.22)^2/(2*2.49))
shiftsBK <- list(edges = c(1, 12, 14, 28),
                 values = c(2.75-0.86, 3.24-0.86, 3.56-3.24, 3.56-3.24),
                 relativeTimes = c(0, 0, 0, 0))
processBK <- "OU"
varianceBK <- (0.22)^2
selection.strengthBK <- 2.49
optimal.valueBK <- 0.86
# Test root.state
root.stateBK <- test.root.state(process = processBK, 
                                root.state = root.stateBK, 
                                variance = varianceBK, 
                                selection.strength = selection.strengthBK, 
                                optimal.value = optimal.valueBK)
# Display parameters
paramsBK = list(variance = varianceBK,
              root.state = root.stateBK,
              shifts = shiftsBK,
              selection.strength = selection.strengthBK, 
              optimal.value = optimal.valueBK)

## Save table history
name <- paste0("_K=", K)

results_estim_EM$params_history[["B-K"]] <- paramsBK
history <- list_to_table.history(results_estim_EM$params_history)
history[,"B-K"]["log_likelihood"] <-log_likelihood.OU(Y_data_bimac, bimac_ape_tree, paramsBK)
write.csv2(history, paste0(PATH, "boite_noire_BK", name, ".csv"))
plot.history.OU.stationnary(results_estim_EM$params_history, paramsBK, bimac_ape_tree, Y_data_bimac, paramsBK, PATH=PATH, paste0("history_plot", name), "B-K")

## Save plots
plot.process.actual(Y.state = Y_data_bimac,
                    Z.state = results_estim_EM$ReconstructedNodesStates,
                    phylo = bimac_ape_tree, 
                    paramsEstimate = results_estim_EM$params)
plot.process.actual(Y.state = Y_data_bimac,
                    Z.state = 0,
                    phylo = bimac_ape_tree, 
                    paramsEstimate = paramsBK)
