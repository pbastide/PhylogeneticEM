##########
## Chelonia EM
rm(list=ls())

PATH <- "../Results/Chelonia/"

# Dependencies of EM
library(phylolm)

###################################################
## Import data
###################################################

## Import data set form package bayou
library(bayou)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

K_max <- 20
data_type <- "chelonia"
ntaxa <- length(data)
K_true <- 16

## Random root ##############################################################################

time_OUshifts <- system.time(res_OUshifts <- OUshifts(y = data,
                                                      phy = tree,
                                                      method = "mbic",
                                                      nmax = K_max,
                                                      check.pruningwise = TRUE))

save.image(paste0(PATH, "chelonia_phylolm", ".RData"))

plot.OUshifts(res_OUshifts, show.tip.label=FALSE)

## Re-fit the model found for likelihood #####################################################
# Tree matrix
tree_prun <- reorder(tree, "pruningwise")
ntaxa <- length(tree_prun$tip.label)
Nedges <- dim(tree_prun$edge)[1]
ROOT <- ntaxa + 1
anc <- tree_prun$edge[, 1]
des <- tree_prun$edge[, 2]
v <- matrix(0, ntaxa + tree_prun$Nnode, ntaxa)
for (i in 1:Nedges) {
  if (des[i] <= ntaxa) v[des[i], des[i]] = 1
  v[anc[i], ] = v[anc[i], ] + v[des[i], ]
}

# Subset corresponding to the shifts found by OUshifts
correspondanceEdges <- function(edges, from, to){
  mm <- match(from$edge[, 2], to$edge[, 2])
  newEdges <- mm[edges]
  return(newEdges)
}
X <- t(v[c(ROOT, des[correspondanceEdges(res_OUshifts$pshift, from = tree, to = tree_prun)]), ])

# Fit
fit <- phylolm(data[match(tree_prun$tip.label, names(data))] ~ X - 1,
               phy = tree_prun,
               model = "OUfixedRoot")

OU_phylolm <- list(Nbr_shifts = res_OUshifts$nshift,
                   Nbr_regimes = res_OUshifts$nshift + 1,
                   lnL = fit$logLik,
                   MlnL = NaN,
                   alpha = res_OUshifts$alpha,
                   half_life = log(2) / res_OUshifts$alpha,
                   sigma = res_OUshifts$sigma2,
                   gamma = res_OUshifts$sigma2 / (2 * res_OUshifts$alpha),
                   time = time_OUshifts[3])

params_OUshifts <- list(variance = res_OUshifts$sigma2,
                        shifts = list(edges = res_OUshifts$pshift,
                                      values = res_OUshifts$shift,
                                      relativeTimes = 0),
                        root.state = list(random = FALSE,
                                          stationary.root = FALSE,
                                          value.root = res_OUshifts$mean,
                                          exp.root = NA,
                                          var.root = NA),
                        selection.strength = res_OUshifts$alpha,
                        optimal.value = res_OUshifts$mean)

save(OU_phylolm, params_OUshifts, file = paste0(PATH, "chelonia_phylolm_summary.RData"))

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = params_OUshifts,
                         adj.root = 1.3,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 1.3,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2",
                         no.margin = TRUE)