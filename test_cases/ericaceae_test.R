rm(list=ls())

library(ape)
library(plyr)
library(quadrupen)
library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(TreeSim)

## Load Ericaceae data
load("../data/ericaceae_data.RData")

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

datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

## Descriptive Statistics ########################################
# Number of species
plot(1:6, sapply(trait_matrix_all, function(z) dim(z)[2]),
     xlab = "Min number of non-NA traits",
     ylab = "Number of tips")

# Percentage missing data
plot(1:6, sapply(trait_matrix_all, function(z) sum(is.na(z))/prod(dim(z))),
     xlab = "Min number of non-NA traits",
     ylab = "Percentage of missing data")

# Summary
trait_matrix <- trait_matrix_all[[1]]
summary(t(trait_matrix))
  
# Histograms non transform
par(mfrow = c(2, 3))
for (i in 1:6){
  hist(trait_matrix[i, ], 100, main = rownames(trait_matrix)[i])
}
# Histograms log transform
par(mfrow = c(2, 3))
for (i in 1:6){
  hist(log(trait_matrix[i, ]), 100, main = paste0("log(", rownames(trait_matrix)[i], ")"))
}
# Histograms log transform only sizes
par(mfrow = c(2, 3))
for (i in 1:2){
  hist(log(trait_matrix[i, ]), 100, main = paste0("log(", rownames(trait_matrix)[i], ")"))
}
for (i in 3:6){
  hist(trait_matrix[i, ], 100, main = rownames(trait_matrix)[i])
}
par(mfrow = c(1,1))

###############################################################################
## EM Univariate ##############################################################
###############################################################################
## Select data
phylo <- subtree_traits[[1]]
trait_matrix <- trait_matrix_all[[1]]

# 0 values
set.seed(17910402)
temp_cl <- replace_zeros(trait_matrix[1, ])
trait_matrix[1, ] <- temp_cl$x
# Log transform sizes
trait_matrix_transform <- trait_matrix
trait_matrix_transform[1:2, ] <- log(trait_matrix_transform[1:2, ])

## Trait i
i <- 1

# Drop NAs
data_uni <- trait_matrix_transform[i, !is.na(trait_matrix_transform[i, ])]
Overlap_traits <- name.check(phylo, data_uni)
tree_uni <- drop.tip(phylo, Overlap_traits$Tree.not.data)

## Re-scale tree to one
height_tree <- node.depth.edgelength(tree_uni)[1]
tree_uni$edge.length <- tree_uni$edge.length / height_tree

# Alpha Values
alpha_grid <- find_grid_alpha(tree_uni,
                              nbr_alpha = 10,
                              factor_up_alpha = 2,
                              factor_down_alpha = 4,
                              quantile_low_distance = 0.001,
                              log_transform = TRUE)

res <- PhyloEM(phylo = tree_uni,
               Y_data = data_uni,
               process = "scOU",
               random.root = FALSE,
               K_max = 35,
               alpha_known = TRUE,
               alpha = alpha_grid,
               tol = list(variance = 10^(-2), 
                          value.root = 10^(-2),
                          log_likelihood = 10^(-2)),
               use_previous = FALSE,
               method.init = "lasso",
               method.selection = "BGH")

save.image(file = paste0("../Results/Test_Cases/trait_scOU_fixed_root_trait_", i, ".RData"))

library(lineprof)
l1 <- lineprof(results_estim_EM <- estimateEM(phylo = subtree_traits,
                                              Y_data = trait_matrix,
                                              process = "BM",
                                              method.init = "default",
                                              Nbr_It_Max = 2,
                                              nbr_of_shifts = 0,
                                              random.root = FALSE,
                                              tol = list(variance = 10^(-2), 
                                                         value.root = 10^(-2))))
shine(l1)

sapply(results_estim_EM$params_history, function(z) attr(z, "log_likelihood")[1])
sapply(results_estim_EM$params_history, function(z) z$root.state$value.root)
sapply(results_estim_EM$params_history, function(z) z$variance[1, 1])

results_estim_EM$params_history$`11`$variance
results_estim_EM$params_history$`12`$variance
results_estim_EM$params_history$`13`$variance

plot(res$alpha_max$capushe_outputBM2@DDSE, newwindow = FALSE)

params_estim_EM <- res$alpha_max$params_select_DDSE_BM2
par(mfrow = c(1, dim(trait_matrix)[1]), mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
for (l in 1:dim(trait_matrix)[1]){
  params <- params_estim_EM
  params$shifts$values <- round(params_estim_EM$shifts$values[l, ], 2)
  params$root.state$value.root <- round(params_estim_EM$root.state$value.root[l], 2)
  plot.data.process.actual(Y.state = trait_matrix[l, ],
                           phylo = subtree_traits, 
                           params = params,
                           adj.root = 0,
                           automatic_colors = TRUE,
                           margin_plot = NULL,
                           cex = 2,
                           bg_shifts = "azure2",
                           bg_beta_0 = "azure2")
  title(rownames(trait_matrix)[l])
}

plot.data.process.actual(Y.state = trait_matrix[2, ],
                         phylo = subtree_traits, 
                         params = NULL)

results_estim_EM_total_large_anther_length <- estimateEM(phylo = subtree_traits,
                                                         Y_data = trait_matrix[2, ],
                                                         process = "BM",
                                                         method.init = "default",
                                                         Nbr_It_Max = 500,
                                                         nbr_of_shifts = 0,
                                                         random.root = FALSE)
sapply(results_estim_EM_total_large_anther_length$params_history, function(z) attr(z, "log_likelihood")[1])

results_estim_EM_total_corolla_tube_length <- estimateEM(phylo = subtree_traits,
                                                         Y_data = trait_matrix[1, ],
                                                         process = "BM",
                                                         method.init = "default",
                                                         Nbr_It_Max = 500,
                                                         nbr_of_shifts = 0,
                                                         random.root = FALSE)
sapply(results_estim_EM_total_corolla_tube_length$params_history, function(z) attr(z, "log_likelihood")[1])

results_estim_EM_measurment_traits <- estimateEM(phylo = subtree_traits,
                                                 Y_data = trait_matrix[1:2, ],
                                                 process = "BM",
                                                 method.init = "default",
                                                 Nbr_It_Max = 500,
                                                 nbr_of_shifts = 0,
                                                 random.root = FALSE)

sapply(results_estim_EM_measurment_traits$params_history, function(z) attr(z, "log_likelihood")[1])
sapply(results_estim_EM_measurment_traits$params_history, function(z) z$root.state$value.root)
sapply(results_estim_EM_measurment_traits$params_history, function(z) z$variance[2, 2])

## Characteristics of the dataset ################################
summary(traits_data)
sum(is.na(traits_data)) / prod(dim(traits_data))

toy_ntaxa <- 128
set.seed(1296)
toy_tree <- rcoal(toy_ntaxa)

toy_p <- 6

toy_variance <- matrix(1, toy_p, toy_p) + diag(9, toy_p) #var(t(trait_matrix), na.rm = TRUE)

toy_root.state <- list(random = FALSE,
                   value.root = c(6, 2, 0.03, 0.007, 0.06, -0.01),
                   # apply(trait_matrix, 1, function(z) median(z, na.rm = TRUE))
                   exp.root = NA,
                   var.root = NA)

toy_shifts = list(edges = c(48, 167),
              values=cbind(c(1, 1, 0.01, 0.001, 0.01, -0.01),
                           c(10, 10, 0.1, 0.01, 0.1, -0.1)),
              relativeTimes = 0)

toy_paramsSimu <- list(variance = toy_variance,
                   shifts = toy_shifts,
                   root.state = toy_root.state)

## Simulate Process
X1 <- simulate_internal(toy_tree,
               p = toy_p,
               root.state = toy_root.state,
               process = "BM",
               variance = toy_variance,
               shifts = toy_shifts)

toy_Y_data <- extract_simulate_internal(X1,"tips","states")

toy_Y_data_miss <- toy_Y_data
set.seed(1122)
nMiss <- floor(toy_ntaxa * toy_p / 100) * 30
miss <- sample(1:(toy_p * toy_ntaxa), nMiss, replace = FALSE)
chars <- (miss - 1) %% toy_p + 1
tips <- (miss - 1) %/% toy_p + 1
# chars <- sample(1:p, nMiss, replace = TRUE)
# tips <- sample(1:ntaxa, nMiss, replace = TRUE)
for (i in 1:nMiss){
  toy_Y_data_miss[chars[i], tips[i]] <- NA
}

par(mfrow = c(1,toy_p), mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
for (l in 1:toy_p){
  toy_params <- toy_paramsSimu
  toy_params$shifts$values <- toy_paramsSimu$shifts$values[l, ]
  toy_params$root.state$value.root <- toy_paramsSimu$root.state$value.root[l]
  plot.data.process.actual(Y.state = toy_Y_data_miss[l, ],
                           phylo = toy_tree, 
                           params = toy_params,
                           adj.root = 0,
                           automatic_colors = TRUE,
                           margin_plot = NULL,
                           cex = 2)
}

set.seed(17920920)
res <- PhyloEM(phylo = subtree_traits,
               Y_data = trait_matrix,
               process = "BM",
               K_max = 10, random.root = FALSE)
save.image(file = paste0("../Results/Test_Cases/trait_BM_toy.RData"))

set.seed(17920920)
results_estim_EM <- estimateEM(phylo = toy_tree,
                               Y_data = toy_Y_data_miss,
                               process = "BM",
                               method.init = "default",
                               Nbr_It_Max = 500,
                               nbr_of_shifts = 2,
                               random.root = FALSE,
                               edges.init = toy_shifts$edges,
                               values.init = matrix(1, toy_p, 2))

par(mfrow = c(1,p), mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
for (l in 1:toy_p){
  params <- results_estim_EM$params
  params$shifts$values <- round(results_estim_EM$params$shifts$values[l, ], 2)
  params$root.state$value.root <- round(results_estim_EM$params$root.state$value.root[l], 2)
  plot.data.process.actual(Y.state = toy_Y_data_miss[l, ],
                           phylo = toy_tree, 
                           params = params,
                           adj.root = 0,
                           automatic_colors = TRUE,
                           margin_plot = NULL,
                           cex = 2)
}

######################################################################
## Pre-Processed EM ##################################################
######################################################################
load("../../../ericaceae_flowers/R_functions/ericaceae_migale_2015-10-10_00-09-12.RData")
set.seed(17920920)
res <- PhyloEM(phylo = subtree_traits,
               Y_data = trait_matrix,
               process = "BM",
               K_max = 10,
               random.root = FALSE,
               tol = list(variance = 10^(-2), 
                          value.root = 10^(-2)),
               Nbr_It_Max = 100,
               estimates = X)
save.image(file = paste0("../Results/Test_Cases/trait_BM_fith_try.RData"))
