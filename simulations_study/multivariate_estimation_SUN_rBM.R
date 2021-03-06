####################
## Parameters
####################
library(doParallel)
library(foreach)
library(PhylogeneticEM)
# library(ape)
# library(glmnet) # For Lasso initialization
# library(robustbase) # For robust fitting of alpha
# library(gglasso)
# library(capushe)
# library(Matrix)
# library(Rcpp)
# library(RcppArmadillo)
# reqpckg <- c("ape", "glmnet", "robustbase", "gglasso", "Matrix", "capushe", "PhylogeneticEM")
reqpckg <- c("PhylogeneticEM")

## Set number of parallel cores
Ncores <- 5

## Define date-stamp for file names
datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
datestamp_day <- format(Sys.time(), "%Y-%m-%d")

## Load simulated data
datestamp_data <- "2016-09-14" # 
savedatafile = "../Results/Simulations_Multivariate/multivariate_simlist"
saveresultfile <- "../Results/Simulations_Multivariate/multivariate_estimations_SUN_rBM"
load(paste0(savedatafile, "_", datestamp_data, "_light.RData"))

# source("R/simulate.R")
# source("R/estimateEM.R")
# source("R/init_EM.R")
# source("R/E_step.R")
# source("R/M_step.R")
# source("R/shutoff.R")
# source("R/generic_functions.R")
# source("R/shifts_manipulations.R")
# source("R/plot_functions.R")
# source("R/parsimonyNumber.R")
# source("R/partitionsNumber.R")
# source("R/model_selection.R")
# sourceCpp("src/upward_downward.cpp")

## These values should be erased by further allocations (generate_inference_files)
n.range <- nrep
inference.index <- 0

## Select data (according to the value of nrep)
nrep <- nrep

## Here n.range should be defined by generate_inference_files.R
simulations2keep <- sapply(simlist, function(x) { x$nrep %in% n.range }, simplify = TRUE)
simlist <- simlist[simulations2keep]
nbrSim <- length(simlist)

# ## Log file
# logfile <- paste0(savedatafile, "_alpha_known-", datestamp_day, "_", inference.index,"_log.txt")
# 
# log <- function(it){
#   txt <- paste0(Sys.time(), " : on batch ", inference.index, ", iteration ", it, " on ", nbrSim, " completed.")
#   writeLines(txt, logfile)
# }

######################
## Estimation Function
######################
estimations_several_K <- function(X, pPCA = FALSE){
  if (pPCA){
    Y_data <- t(X$Y_data)
    if (anyNA(Y_data))   return(list(sim = X, res = NULL))
    rownames(Y_data) <- trees[[paste0(X$ntaxa)]]$tip.label
    ## Do a pPCA
    Y_data_pPCA <- phytools::phyl.pca(trees[[paste0(X$ntaxa)]], Y_data)
    Y_data <- t(Y_data_pPCA$S)
  } else {
    Y_data <- X$Y_data
  }
  alpha_grid <- find_grid_alpha(trees[[paste0(X$ntaxa)]],
                                nbr_alpha = 10,
                                factor_up_alpha = 2,
                                factor_down_alpha = 3,
                                quantile_low_distance = 0.0001,
                                log_transform = TRUE)
  time_SUN <- system.time(
  res <- PhyloEM(phylo = trees[[paste0(X$ntaxa)]],
                 Y_data = Y_data,
                 process = "scOU",
                 K_max = max(K_try[[paste0(X$ntaxa)]]) + 5,
                 random.root = TRUE,
                 stationary.root = TRUE,
                 alpha = alpha_grid[-1],
                 save_step = FALSE,
                 Nbr_It_Max = 2000,
                 tol = list(variance = 10^(-2), 
                            value.root = 10^(-2),
                            log_likelihood = 10^(-2)),
                 min_params = list(variance = 0, 
                                   value.root = -10^(5), 
                                   exp.root = -10^(5), 
                                   var.root = 0,
                                   selection.strength = 0),
                 method.variance = "upward_downward",
                 method.init = "lasso",
                 use_previous = FALSE,
                 method.selection = c("BirgeMassart1", "BirgeMassart2", "BGHml",
                                      "BGHlsq", "BGHlsqraw", "BGHmlraw"),
                 K_lag_init = 5,
                 light = FALSE)
  )
  res <- add_total_time(res, time_SUN["elapsed"])
  ret <- list(sim = X,
              res = res)
  return(ret)
}

# enlight_res <- function(res){
#   lres <- vector("list", 4)
#   lres[1:3] <- res[1:3]
#   lmax <- res$alpha_max$BGH[c("params_select", "params_raw", "params_init_estim",
#                               "results_summary",
#                               # "Yhat", "Zhat", "Yvar", "Zvar",
#                               "m_Y_estim", "edge.quality")]
#   lres$alpha_max <- res$alpha_max
#   lres$alpha_max$BGH <- lmax
#   return(lres)
# }

add_total_time <- function(res, tot_time){
  sum_times <- sum(sapply(res[grep("alpha_[[:digit:]]", names(res))],
                          function(z) z$results_summary$time))
  res$alpha_max$results_summary$total_time <- tot_time
  res$alpha_max$results_summary$sum_times <- sum_times
  sels <- names(res$alpha_max)[c(grep("DDSE", names(res$alpha_max)),
                                 grep("Djump", names(res$alpha_max)),
                                 grep("BGH", names(res$alpha_max)))]
  for (i in sels){
    res$alpha_max[[i]]$results_summary$total_time <- tot_time
    res$alpha_max[[i]]$results_summary$sum_times <- sum_times
    # res$alpha_max[[i]]$results_summary <- as.data.frame(res$alpha_max[i]$results_summary)
  }
  sels <- names(res$alpha_min)[c(grep("BGH", names(res$alpha_min)))]
  for (i in sels){
    res$alpha_min[[i]]$results_summary$total_time <- tot_time
    res$alpha_min[[i]]$results_summary$sum_times <- sum_times
    # res$alpha_max[[i]]$results_summary <- as.data.frame(res$alpha_max[i]$results_summary)
  }
  sels <- names(res$alpha_min_raw)[c(grep("BGH", names(res$alpha_min_raw)))]
  for (i in sels){
    res$alpha_min_raw[[i]]$results_summary$total_time <- tot_time
    res$alpha_min_raw[[i]]$results_summary$sum_times <- sum_times
    # res$alpha_max[[i]]$results_summary <- as.data.frame(res$alpha_max[i]$results_summary)
  }
  return(res)
}


estimations_several_K_ak <- function(X){
  alpha_grid <- X$alpha
  res <- PhyloEM(phylo = trees[[paste0(X$ntaxa)]],
                 Y_data = X$Y_data,
                 process = "scOU",
                 K_max = max(K_try[[paste0(X$ntaxa)]]) + 5,
                 random.root = TRUE,
                 stationary.root = TRUE,
                 alpha = alpha_grid,
                 save_step = FALSE,
                 Nbr_It_Max = 2000,
                 tol = list(variance = 10^(-2), 
                            value.root = 10^(-2),
                            log_likelihood = 10^(-2)),
                 min_params = list(variance = 0, 
                                   value.root = -10^(5), 
                                   exp.root = -10^(5), 
                                   var.root = 0,
                                   selection.strength = 0),
                 method.init = "lasso",
                 use_previous = FALSE,
                 method.selection = c("BirgeMassart1", "BirgeMassart2"),
                 K_lag_init = 5)
  res <- add_total_time(res)
  ret <- list(sim = X,
              res = res)
  return(ret)
}

############
## Estimations (alpha on a grid)
############

## Separate "favorable" values from others
simparams_keep <- subset(simparams, nrep %in% n.range)
favorables <- simparams_keep$alpha >= 3 & simparams_keep$K <= 7 & simparams_keep$NA_per == 0

## FAVORABLES ##
## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time_alpha_gird_fav <- system.time(
  simestimations_fav <- foreach(i = simlist[favorables], .packages = reqpckg) %dopar%
  {
    estimations_several_K(i, pPCA=TRUE)
  }
)
# Stop the cluster (parallel)
stopCluster(cl)

## rename object and save
assign(paste0("simestimations_fav_", inference.index),
       simestimations_fav)
rm(simestimations_fav)

save.image(paste0(saveresultfile, "favorables-", datestamp_day, "_", inference.index, ".RData"))

## NOT FAVORABLES ##
## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time_alpha_gird_unfav <- system.time(
  simestimations_unfav <- foreach(i = simlist[!favorables], .packages = reqpckg) %dopar%
  {
    estimations_several_K(i, pPCA=TRUE)
  }
)
# Stop the cluster (parallel)
stopCluster(cl)

## group favorables and unfavorables
simestimations <- vector(mode = "list", length = length(favorables))
simestimations[favorables] <- eval(as.name(paste0("simestimations_fav_", inference.index)))
simestimations[!favorables] <- simestimations_unfav

rm(simestimations_unfav)
rm(list = paste0("simestimations_fav_", inference.index))

## rename object and save
assign(paste0("simestimations_", inference.index),
       simestimations)
rm(simestimations)

save.image(paste0(saveresultfile, "-", datestamp_day, "_", inference.index, ".RData"))

# ############
# ## Estimations (alpha known)
# ############
# 
# ## Register parallel backend for computing
# cl <- makeCluster(Ncores)
# registerDoParallel(cl)
# 
# ## Parallelized estimations
# time_alpha_known <- system.time(
#   simestimations <- foreach(i = simlist, .packages = reqpckg) %dopar%
#   {
#     estimations_several_K_ak(i)
#   }
# )
# # Stop the cluster (parallel)
# stopCluster(cl)
# 
# ## rename object and save
# assign(paste0("simestimations_alpha_known_", inference.index), 
#        simestimations)
# rm(simestimations)
# 
# save.image(paste0(saveresultfile, "alpha_known_", datestamp_day, "_", inference.index, ".RData"))
