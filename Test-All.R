rm(list=ls())

WD <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
#WD <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
setwd(WD)

library(ape)
library(plyr)
library(quadrupen)
library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(TreeSim)

library(testthat)

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
source("Phylogenetic-EM/partitionsNumber.R")

test_dir(paste0(WD, "/Phylogenetic-EM/tests"))
