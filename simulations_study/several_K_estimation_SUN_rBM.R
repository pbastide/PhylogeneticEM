####################
## Parameters
####################
library(doParallel)
library(foreach)
library(ape)
library(glmnet) # For Lasso initialization
library(robustbase) # For robust fitting of alpha
reqpckg <- c("ape", "glmnet", "robustbase")

## Set number of parallel cores
Ncores <- 3

## Define date-stamp for file names
datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
datestamp_day <- format(Sys.time(), "%Y-%m-%d")

## Load simulated data
datestamp_data <- "2015-03-17" #format(Sys.time(), "%Y-%m-%d")
savedatafile = "../Results/Simulations_Several_K/several_K_simlist"
saveresultfile <- "../Results/Simulations_Several_K/several_K_estimations_SUN_rBM_"
load(paste0(savedatafile, "_", datestamp_data, ".RData"))

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

## These values should be erased by further allocations (generate_inference_files)
n.range <- n
inference.index <- 0

## Select data (according to the value of n)
n <- n

## Here n.range should be defined by generate_inference_files.R
simulations2keep <- sapply(simlist, function(x) { x$n %in% n.range }, simplify = TRUE)
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
estimations_several_K <- function(X){
  alpha_grid <- find_grid_alpha(trees[[paste0(X$ntaxa)]],
                                nbr_alpha = 10,
                                factor_up_alpha = 2,
                                factor_down_alpha = 3,
                                quantile_low_distance = 0.0001,
                                log_transform = TRUE)
  res <- PhyloEM(phylo = trees[[paste0(X$ntaxa)]],
                 Y_data = X$Y_data,
                 process = "scOU",
                 K_max = max(K_try[[paste0(X$ntaxa)]]),
                 random.root = TRUE,
                 stationnary.root = TRUE,
                 alpha = alpha_grid[-1],
                 save_step = FALSE,
                 Nbr_It_Max = 2000,
                 tol = list(variance = 10^(-2), 
                            value.root = 10^(-2),
                            log_likelihood = 10^(-2)),
                 method.init = "lasso",
                 use_previous = FALSE,
                 method.selection = "BGH")
  res <- add_total_time(res)
  res <- enlight_res(res)
  ret <- list(sim = X,
              res = res)
  return(ret)
}

enlight_res <- function(res){
  lres <- vector("list", 4)
  lres[1:3] <- res[1:3]
  lmax <- res$alpha_max$BGH[c("params_select", "params_raw", "params_init_estim",
                              "results_summary",
                              # "Yhat", "Zhat", "Yvar", "Zvar",
                              "m_Y_estim", "edge.quality")]
  lres$alpha_max <- res$alpha_max
  lres$alpha_max$BGH <- lmax
  return(lres)
}

add_total_time <- function(res){
  tot_time <- sum(sapply(res[grep("alpha_[[:digit:]]", names(res))],
                         function(z) z$results_summary$time))
  res$alpha_max$results_summary$total_time <- tot_time
  res$alpha_max$BGH$results_summary$total_time <- tot_time
  res$alpha_max$BGH$results_summary <- as.data.frame(res$alpha_max$BGH$results_summary)
  return(res)
}


estimations_several_K_ak <- function(X){
  alpha_grid <- X$alpha
  res <- PhyloEM_core(phylo = trees[[paste0(X$ntaxa)]],
                      Y_data = X$Y_data,
                      process = "scOU",
                      K_max = max(K_try[[paste0(X$ntaxa)]]),
                      random.root = TRUE,
                      stationnary.root = TRUE,
                      alpha = alpha_grid,
                      save_step = FALSE,
                      Nbr_It_Max = 2000,
                      tol = list(variance = 10^(-2), 
                                 value.root = 10^(-2),
                                 log_likelihood = 10^(-2)),
                      method.init = "lasso",
                      use_previous = FALSE)
  ret <- list(sim = X,
              res = res)
  return(ret)
}

############
## Estimations (alpha on a grid)
############

## Separate "favorable" values from others
simparams_keep <- subset(simparams, n %in% n.range)
favorables <- simparams_keep$gamma <= 1 & simparams_keep$alpha >= 3 & simparams_keep$K <= 5

## FAVORABLES ##
## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time_alpha_gird_fav <- system.time(
  simestimations_fav <- foreach(i = simlist[favorables][1:3], .packages = reqpckg) %dopar%
  {
    estimations_several_K(i)
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
  simestimations_unfav <- foreach(i = simlist[!favorables][1:3], .packages = reqpckg) %dopar%
  {
    estimations_several_K(i)
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
# ## Separate "favorable" values from others
# simparams_keep <- subset(simparams, n %in% n.range)
# favorables <- simparams_keep$gamma <= 1 & simparams_keep$alpha >= 3 & simparams_keep$K <= 5
# 
# ## FAVORABLES ##
# ## Register parallel backend for computing
# cl <- makeCluster(Ncores)
# registerDoParallel(cl)
# 
# ## Parallelized estimations
# time_alpha_known <- system.time(
#   simestimations_fav <- foreach(i = simlist[favorables], .packages = reqpckg) %dopar%
#   {
#     estimations_several_K_ak(i)
#   }
# )
# # Stop the cluster (parallel)
# stopCluster(cl)
# 
# ## rename object and save
# assign(paste0("simestimations_fav_alpha_known_", inference.index), 
#        simestimations_fav)
# rm(simestimations_fav)
# 
# save.image(paste0(saveresultfile, "favorables_alpha_known_", datestamp_day, "_", inference.index, ".RData"))
# 
# ## NOT FAVORABLES ##
# ## Register parallel backend for computing
# cl <- makeCluster(Ncores)
# registerDoParallel(cl)
# 
# ## Parallelized estimations
# time_alpha_known_unfav <- system.time(
#   simestimations_unfav <- foreach(i = simlist[!favorables], .packages = reqpckg) %dopar%
#   {
#     estimations_several_K_ak(i)
#   }
# )
# # Stop the cluster (parallel)
# stopCluster(cl)
# 
# ## group favorables and unfavorables
# simestimations <- vector(mode = "list", length = length(favorables))
# simestimations[favorables] <- eval(as.name(paste0("simestimations_fav_alpha_known_", inference.index)))
# simestimations[!favorables] <- simestimations_unfav
# 
# rm(simestimations_unfav)
# rm(list = paste0("simestimations_fav_alpha_known_", inference.index))
# 
# ## rename object and save
# assign(paste0("simestimations_alpha_known_", inference.index), 
#        simestimations)
# rm(simestimations)
# 
# save.image(paste0(saveresultfile, "alpha_known_", datestamp_day, "_", inference.index, ".RData"))

