#############################################
## Model selection
#############################################

##############
## Parameters
##############

library(doParallel)
library(foreach)
library(PhylogeneticEM)
library(plyr)
# library(capushe)

saveresultdir <- "../Results/Simulations_Multivariate"
saveresultname <- "multivariate_estimations_SUN_rBM"
datestamp_day_cut <- "2016-09-2"
datestamp_day <- "2016-09-28"
ak <- "" #"_alpha_known" # ""

files_list <- list.files(path = saveresultdir,
                         pattern = paste0(saveresultname, ak, "-", datestamp_day_cut),
                         full.names = TRUE)

fct_to_keep <- c("saveresultdir", "saveresultname",
                 "datestamp_day_cut", "datestamp_day", "files_list")
##############
## Add least squares
##############
U_trees <- lapply(trees, incidence.matrix.full)
names(U_trees) <- names(trees)

add_lsq <- function(X){
  tree_one <- trees[[paste0(X$sim$ntaxa)]]
  U_tree_one <- U_trees[[paste0(X$sim$ntaxa)]]
  times_shared_one <- times_shared[[paste0(X$sim$ntaxa)]]
  fun <- function(params){
    tmpsim <- PhylogeneticEM:::simulate_internal(phylo = tree_one,
                                                 process = "scOU",
                                                 p = attr(params, "p_dim"),
                                                 root.state = params$root.state,
                                                 shifts = params$shifts,
                                                 variance = params$variance,
                                                 optimal.value = params$optimal.value,
                                                 selection.strength = params$selection.strength,
                                                 simulate_random = FALSE,
                                                 U_tree = U_tree_one,
                                                 times_shared = times_shared_one)
    ## Mean at tips with estimated parameters
    m_Y_estim <- PhylogeneticEM:::extract_simulate_internal(tmpsim, 
                                                            where="tips",
                                                            what="expectations") 
  }
  
  nums <- grep("alpha_", names(X$res))
  for (i in nums){
    m_Y_estim <- lapply(X$res[[i]]$params_estim, fun)
    X$res[[i]]$m_Y_estim_new <- m_Y_estim
    X$res[[i]]$results_summary$least_squares_raw <- sapply(m_Y_estim,
    function(z) sum((X$sim$Y_data - z)^2, na.rm = TRUE))
    # X$res[[i]]$results_summary$least_squares_bis <- sapply(X$res[[i]]$m_Y_estim,
                                                       # function(z) sum((X$sim$Y_data - z)^2, na.rm = TRUE))
    X$res[[i]]$results_summary$least_squares <- sapply(X$res[[i]]$params_estim, function(z) sum(diag(z$variance)/(2*z$selection.strength)))
  }
  return(X)
}

fct_to_keep <- c(fct_to_keep, "U_trees", "add_lsq")
##############
## Model Selection
##############
model_selection_tmp <- function(simres){
  X <- simres$sim
  res <- PhyloEM(phylo = trees[[paste0(X$ntaxa)]],
                 Y_data = X$Y_data,
                 process = "scOU",
                 K_max = max(K_try[[paste0(X$ntaxa)]]) + 5,
                 random.root = TRUE,
                 stationary.root = TRUE,
                 estimates = simres$res,
                 method.selection = c("BGHlsq", "BGHml", "BGHlsqraw", "BGHmlraw"))
  simres$res <- res
  return(simres)
}

fct_to_keep <- c(fct_to_keep, "model_selection_tmp")
######################
## Loop on inference files - PhyloEM
######################
for (file in files_list){
  load(file)
  gc()
  ## simestimations
  simests <- as.name(paste0("simestimations_", ak, inference.index))
  
  ## Add least squares
  assign(paste(simests), lapply(eval(simests), add_lsq))
  
  ## Model Selection
  assign(paste(simests), lapply(eval(simests), model_selection_tmp))
  
  save.image(file)
  rm(list = ls()[-which(ls() %in% fct_to_keep)])
}