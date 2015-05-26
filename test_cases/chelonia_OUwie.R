##########
## Test Case : Chelonian carapace evolution.
## Used in :
## - Uyeda 2014
## - Eastman 2011 (link to data : https://github.com/eastman/auteur/tree/master/auteur/data)
## -> Jaffe 2011 (original data)

rm(list=ls())

PATH <- "../Results/Chelonia/"

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
library(bayou)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

plot(tree, show.tip.label = FALSE)

## Import Data from Jaffe
vv <- read.csv("../data/Chelonia_habitats_vector.txt", header = FALSE)
chel_df <- NULL
chel_df$species <- as.character(vv[c(2:44, 178:222, 362:406, 546:590, 730:774, 914:916), ])
chel_df$accession_number <-  as.character(vv[c(46:88, 224:268, 408:452, 592:636, 776:820, 918:920), ])
tmp <- vv[c(90:132, 270:314, 454:498, 638:682, 822:866, 922:924), ]
chel_df$length <- as.numeric(levels(tmp))[tmp]
chel_df$habitat <-  factor(vv[c(134:176, 316:360, 500:544, 684:728, 868:912, 926:928), ])
chel_df <- as.data.frame(chel_df)
# Correspondance
chel_df <- chel_df[match(names(data), chel_df$species), ]

## Clustering Habitat
clusters <- as.vector(chel_df$habitat)
Nbr_ParSols <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
ParSols <- extract.enumerate_parsimony(enumerate_parsimony(tree, clusters))
Nbr_shifts_OUwie <- extract.parsimonyCost(parsimonyCost(tree, clusters), 227)

## Plot clustering
colors_habitat <- chel_df$habitat
levels(colors_habitat) <- 1:length(levels(colors_habitat))

colors_regimes <- as.factor(ParSols[1, tree$edge[,2]])
levels(colors_regimes) <- 1:length(levels(colors_regimes))

plot(tree, show.tip.label = FALSE, edge.color = colors_regimes, x.lim = 300)
tiplabels(frame = "circle", col = colors_habitat, pch = 16, adj = 4)
legend("topright", legend = levels(chel_df$habitat), col = levels(colors_habitat), pch = 16)
axisPhylo()

###############################################################################
## OUwie
###############################################################################
library(OUwie)
PATH <- "../Results/Chelonia/"

ntaxa <- length(tree$tip.label)
## format data for OUwie
treeWIE <- tree
treeWIE$node.label <- ParSols[1,(length(tree$tip.label)+1):dim(ParSols)[2]]
dataWIE <- chel_df[ ,c("species", "habitat", "length")]
dataWIE$length <- log(dataWIE$length)

time_OUwie <- system.time(res_OUwie <- OUwie(phy = treeWIE,
                                             data = dataWIE,
                                             model = "OUM"))

time_BMwie <- system.time(res_OUwie_BM <- OUwie(phy = treeWIE,
                                                data = dataWIE,
                                                model = "BM1"))

OU_wie <- list(Nbr_shifts = Nbr_shifts_OUwie,
               Nbr_regimes = 4,
               lnL = res_OUwie$loglik,
               MlnL = NaN,
               alpha = res_OUwie$solution["alpha", 1],
               half_life = log(2)/res_OUwie$solution["alpha", 1],
               sigma = res_OUwie$solution["sigma.sq", 1],
               gamma = res_OUwie$solution["sigma.sq", 1]/(2 * res_OUwie$solution["alpha", 1]),
               time = unname(time_OUwie[3]))

save.image(paste0(PATH, "chelonia_OUwie", ".RData"))
save(chel_df, clusters, Nbr_ParSols, ParSols, Nbr_shifts_OUwie,
     colors_habitat, colors_regimes,
     OU_wie, file = paste0(PATH, "chelonia_OUwie_summary", ".RData"))
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