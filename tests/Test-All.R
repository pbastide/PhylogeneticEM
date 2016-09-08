rm(list=ls())

library(ape)
library(plyr)
# library(quadrupen)
# library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(TreeSim)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)

library(testthat)

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
sourceCpp("src/upward_downward.cpp")

test_dir("tests")
