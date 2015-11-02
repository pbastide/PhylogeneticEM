rm(list=ls())

library(ape)
library(plyr)
library(quadrupen)
library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(TreeSim)

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

## Load Ericaceae data
load("../data/ericaceae_data.RData")

## Get a single data frame ################################################

## Select Columns: Corrola lenght (2), anther length (7)
linear_measures <- data_reduced[, c(1, 2, 7)]
## Compute species mean (one value per species)
mean_linear_measures_length <- aggregate(cbind(corolla_tube_length, total_large_anther_length) ~ species,
                                         data = linear_measures,
                                         mean,
                                         na.action = na.pass)
## Add PC1 and PC2 from Corolla and Anther Shape
# anther
# pca_anthers <- data.frame(pca_efou_anthers$x[, 1:2])
pca_anthers <- data.frame(all_data_pcs_anthers_unique[, 1:3])
colnames(pca_anthers) <- c("species", "anther.PC1", "anther.PC2")
# pca_anthers$species <- rownames(pca_anthers)
# Corolla
# pca_corolla <- data.frame(pca_efou_flowers$x[, 1:2])
pca_corolla <- data.frame(all_data_pcs_flowers_unique[, 1:3])
colnames(pca_corolla) <- c("species", "corolla.PC1", "corolla.PC2")
# pca_corolla$species <- rownames(pca_corolla)
# Merge
dd <- merge(mean_linear_measures_length, pca_anthers, by = "species", all = TRUE) # Match: 4 !!
traits_data <- merge(dd, pca_corolla, by = "species", all = TRUE)
## Drop NA (entire ligne is NA)
traits_data_all <- vector("list", 6)
trait_matrix_all <- vector("list", 6)
percentage_missing <- vector(length = 6)
for (i in 1:6){ # i is the minimal number of non-NA trait
  traits_data_all[[i]] <- traits_data[rowSums(is.na(traits_data)) < (dim(traits_data)[2] - i),]
  rownames(traits_data_all[[i]]) <- traits_data_all[[i]][,1]
  ## Check compatibility with the tree
  Overlap_traits <- name.check(tree, traits_data_all[[i]])
  subtree_traits <- drop.tip(tree, Overlap_traits$Tree.not.data)
  # Sort data in same order as tree
  match_traits <- match(subtree_traits$tip.label, rownames(traits_data_all[[i]]))
  sortedData_traits <- traits_data_all[[i]][match_traits, -1]
  trait_matrix_all[[i]] <- t(as.matrix(sortedData_traits))
  percentage_missing[i] <- sum(is.na(trait_matrix_all[[i]]))/prod(dim(trait_matrix_all[[i]])) * 100
}

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

# Zero values for lengths ?
sum(trait_matrix[1, ] == 0, na.rm = TRUE)
sum(trait_matrix[2, ] == 0, na.rm = TRUE)
# Replace 0 with small value
corrola_lg <- trait_matrix[1, ]
nbrzeros <- sum(trait_matrix[1, ] == 0, na.rm = TRUE)
histogram <- hist(corrola_lg, 1000, plot = FALSE)
mean <- min(corrola_lg[corrola_lg > 0], na.rm = TRUE) / 10
sd <- mean / 10
corrola_lg[corrola_lg == 0] <- rnorm(nbrzeros, mean, sd)
  
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
par(mfrow = 1)

## EM ############################################################
set.seed(17920920)
res <- PhyloEM(phylo = subtree_traits,
               Y_data = trait_matrix,
               process = "BM",
               K_max = 10,
               random.root = FALSE,
               tol = list(variance = 10^(-2), 
                        value.root = 10^(-2)),
               Nbr_It_Max = 100)
save.image(file = paste0("../Results/Test_Cases/trait_BM_fith_try.RData"))

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
X1 <- simulate(toy_tree,
               p = toy_p,
               root.state = toy_root.state,
               process = "BM",
               variance = toy_variance,
               shifts = toy_shifts)

toy_Y_data <- extract.simulate(X1,"tips","states")

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
