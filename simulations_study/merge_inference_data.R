#############################################
## Bind files together several K
#############################################

##############
## Parameters
##############

require(doParallel)
require(foreach)
require(ape)
#require(glmnet) # For Lasso initialization
require(quadrupen) # For Lasso initialization
require(robustbase) # For robust fitting of alpha
library(TreeSim) # For simulation of the tree
require(ggplot2) # Plot
library(scales) # plot
library(reshape2) # Plot
library(grid) # Plot
library(plyr)

## Source functions
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

saveresultfile = "../Results/Simulations_Several_K/several_K_estimations"
datestamp_day <- "2015-03-06"

simestimations_alpha_known <- NULL

for (inference.index in 1:40){
  load(paste0(saveresultfile, "_alpha_known-", datestamp_day, "_", inference.index, ".RData"))
  siminf <- as.name(paste0("simestimations_alpha_known_", inference.index))
  simestimations_alpha_known <- c(simestimations_alpha_known, 
                                  eval(siminf))
  rm(list = paste0("simestimations_alpha_known_", inference.index))
}

save.image(paste0(saveresultfile, "_alpha_known-", datestamp_day, "_all", ".RData"))