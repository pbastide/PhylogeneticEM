rm(list=ls())

WD <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
#WD <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
setwd(WD)

library(ape)
library(plyr)

library(testthat)

source("Phylogenetic-EM/generic_functions.R")
source("Phylogenetic-EM/parcimonyNumber.R")

test_dir(paste0(WD, "/Phylogenetic-EM/tests"))