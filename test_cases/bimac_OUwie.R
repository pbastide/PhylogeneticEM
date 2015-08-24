rm(list=ls())

PATH <- "../Results/Test_Cases/"

library(ape)
library(glmnet) # For Lasso initialization
library(quadrupen)
#library(nlme) # For second derivative computation
#library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
library(TreeSim)
library(plyr)
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

###############################################################################
## Import DataSet with habitat clustering
###############################################################################
## Import data set form package geiger
library(ouch)
source("R/convert.R")
# Anolis bimaculatus lizard size data
data(bimac)
data_type <- "bimac"
bimac_ouch_tree <- with(bimac, ouchtree(node, ancestor, time, species))
tree <- convert.pmc(bimac_ouch_tree)
data <- log(bimac$size[23:45])
names(data) <- bimac$species[23:45]

body_pattern <- bimac$OU.LP[23:45]

plot(tree, show.tip.label = FALSE)

## Clustering Habitat
clusters <- as.vector(body_pattern)
Nbr_ParSols <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
ParSols <- extract.enumerate_parsimony(enumerate_parsimony(tree, clusters))
Nbr_shifts_OUwie <- extract.parsimonyCost(parsimonyCost(tree, clusters), ntaxa + 1)

## Plot clustering
colors_habitat <- body_pattern
levels(colors_habitat) <- c("#400080FF", "#800000FF", "#008080FF")

colors_regimes <- as.factor(ParSols[1, tree$edge[,2]])
levels(colors_regimes) <- c("#400080FF", "#800000FF", "#008080FF")

plot.data.process.actual(Y.state = data,
                         normalize = FALSE,
                         phylo = tree, 
                         params = NULL,
                         color_characters = colors_habitat,
                         color_edges = colors_regimes)

###############################################################################
## OUwie
###############################################################################
library(OUwie)
PATH <- "../Results/Test_Cases/"

ntaxa <- length(tree$tip.label)
## format data for OUwie
treeWIE <- tree
treeWIE$node.label <- ParSols[1,(length(tree$tip.label)+1):dim(ParSols)[2]]
dataWIE <- bimac[23:45 ,c("species", "OU.LP", "size")]
dataWIE$size <- log(dataWIE$size)

time_OUwie <- system.time(res_OUwie <- OUwie(phy = treeWIE,
                                             data = dataWIE,
                                             model = "OUM"))

time_BMwie <- system.time(res_OUwie_BM <- OUwie(phy = treeWIE,
                                                data = dataWIE,
                                                model = "BM1"))

OU_wie <- list(Nbr_shifts = Nbr_shifts_OUwie,
               Nbr_regimes = 3,
               lnL = res_OUwie$loglik,
               MlnL = NaN,
               alpha = res_OUwie$solution["alpha", 1],
               half_life = log(2)/res_OUwie$solution["alpha", 1],
               sigma = res_OUwie$solution["sigma.sq", 1],
               gamma = res_OUwie$solution["sigma.sq", 1]/(2 * res_OUwie$solution["alpha", 1]),
               time = unname(time_OUwie[3]))

save.image(paste0(PATH, "bimac_OUwie", ".RData"))
save(bimac, clusters, Nbr_ParSols, ParSols, Nbr_shifts_OUwie,
     colors_habitat, colors_regimes,
     OU_wie, file = paste0(PATH, "bimac_OUwie_summary", ".RData"))
#rm(vv, treeWIE, dataWIE, res_OUwie)

## Equivalent solutions
times_shared <- compute_times_ca(tree)
t_tree <-  min(node.depth.edgelength(tree)[1:ntaxa])
T_tree <- incidence.matrix(tree)
ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                  selection.strength = res_OUwie$solution[1,1],
                                                  times_shared = times_shared)
T_tree_ac <- T_tree * ac_tree

betas_wie <- res_OUwie$theta[,1]
names(betas_wie) <- colnames(res_OUwie$solution)
betas_wie <- unname(betas_wie[ParSols[1,]])
shifts_wie <- compute_shifts_from_betas(phylo = tree, betas = betas_wie)
#shifts_wie_all_edges <- apply(ParSols, 1, function(z) allocate_shifts_from_regimes(phylo = tree, regimes = z))

eq_shifts_edges_wie <- equivalent_shifts_edges(tree, 
                                               shifts_wie$edges)
eq_shifts_values_wie <- equivalent_shifts_values(tree,
                                                 shifts = shifts_wie,
                                                 beta_0 = betas_wie[length(tree$tip.label)+1],
                                                 eq_shifts_edges = eq_shifts_edges_wie,
                                                 T_tree_ac = T_tree_ac)

fun <- function(edges, values){
  return(list(edges = edges,
              values = values,
              relativeTimes = rep(0, length(values))))
}

shifts_eq <- mapply(fun, alply(eq_shifts_edges_wie, 2), alply(eq_shifts_values_wie$shifts_values, 2), SIMPLIFY = FALSE)

betas_eq <- mapply( function(beta_0, shifts) compute_betas(tree, beta_0, shifts), as.list(eq_shifts_values_wie$betas_0), shifts_eq)

plot_equivalent_shifts.actual(tree, eq_shifts_edges_wie, eq_shifts_values_wie, use.edge.length = FALSE, adj = 0)

save.image(paste0(PATH, "chelonia_OUwie", ".RData"))
# save.image(paste0(PATH, "estimation_OUwie", ".RData"))