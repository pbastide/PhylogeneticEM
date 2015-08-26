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

###########################################################################
###########################################################################

#############################
## Plot tree for parsimony figures
tree <- read.tree(text="(((C,(T,T)),C),(A,A));")
ntaxa <- 6
clusters_tips <- c(1, 2, 2, 1, 0, 0)
colors_tip <- as.factor(clusters_tips)
levels(colors_tip) <- c("black", "gray", "white")

shifts_pars_1 <- c(1, 4)
colors_shifts_pars_1 <- c("gray", "white")
clusters_nodes_pars_1 <- allocate_regimes_from_shifts(tree, shifts_pars_1)
clusters_nodes_pars_1 <- as.factor(clusters_nodes_pars_1[(ntaxa + 1):(2*ntaxa - 1)])
levels(clusters_nodes_pars_1) <- c("black", "gray", "white")
plot(tree, 
     type = "cladogram",
     show.tip.label = FALSE,
     edge.width = 10,
     no.margin = TRUE,
     direction = "downwards")
tiplabels(text = rep("", ntaxa),
          frame = "circle",
          cex = 1.5,
          bg = as.vector(colors_tip))
nodelabels(text = rep("", ntaxa-1),
                frame = "circle",
                cex = 1.5,
                bg = as.vector(clusters_nodes_pars_1))
edgelabels(text = rep("         ", 2),
           edge = shifts_pars_1,
           cex = 1,
           bg = colors_shifts_pars_1)

shifts_pars_2 <- c(4, 8)
colors_shifts_pars_2 <- c("white", "black")
clusters_pars_2 <- allocate_regimes_from_shifts(tree, shifts_pars_2)
clusters_nodes_pars_2 <- as.factor(clusters_pars_2[(ntaxa + 1):(2*ntaxa - 1)])
levels(clusters_nodes_pars_2) <- c("gray", "white", "black")
plot(tree, 
     type = "cladogram",
     show.tip.label = FALSE,
     edge.width = 10,
     no.margin = TRUE,
     direction = "downwards")
tiplabels(text = rep("", ntaxa),
          frame = "circle",
          cex = 1.5,
          bg = as.vector(colors_tip))
nodelabels(text = rep("", ntaxa-1),
           frame = "circle",
           cex = 1.5,
           bg = as.vector(clusters_nodes_pars_2))
edgelabels(text = rep("         ", 2),
           edge = shifts_pars_2,
           cex = 1,
           bg = colors_shifts_pars_2)

shifts_non_pars <- c(1, 3, 7)
colors_shifts_non_pars <- c("white", "gray", "gray")
clusters_non_pars <- allocate_regimes_from_shifts(tree, shifts_non_pars)
clusters_nodes_non_pars <- as.factor(clusters_non_pars[(ntaxa + 1):(2*ntaxa - 1)])
levels(clusters_nodes_non_pars) <- c("black", "white")
plot(tree, 
     type = "cladogram",
     show.tip.label = FALSE,
     edge.width = 10,
     no.margin = TRUE,
     direction = "downwards")
tiplabels(text = rep("", ntaxa),
          frame = "circle",
          cex = 1.5,
          bg = as.vector(colors_tip))
nodelabels(text = rep("", ntaxa-1),
           frame = "circle",
           cex = 1.5,
           bg = as.vector(clusters_nodes_non_pars))
edgelabels(text = rep("         ", 3),
           edge = shifts_non_pars,
           cex = 1,
           bg = colors_shifts_non_pars)


#############################
## Test of several miscelaneous functions
#############################

SimuBM <- rTrait(Cetacea_Autocorrelated,model="BM",parameters=list(ancestral.state=0,sigma2=1,optimal.value=0,alpha=0))

plot(Cetacea_Autocorrelated,show.tip.label = FALSE)

# tr <- reorder(Cetacea_Autocorrelated, order = "cladewise", index.only = FALSE) #Put in cladwise order
# nod <- 206 # We want the subtree of nod
# i <- which(tr$edge[, 2] == nod) # number of the branch of anc
# anc <- getAncestor(tr,nod) # ancestor of nod
# tmp <- which(tr$edge[, 1] == anc) # number of the branch for wich anc is an ancestor
# j <- tmp[which(tmp == i) + 1] # Next time anc appears
# tr$edge[(i+1):(j-1), ] # Clade
# tr[tr$edge[(i+1):(j-1), ],]

clade <- extract.clade(Cetacea_Autocorrelated, 210, root.edge = 0, interactive = FALSE)
rest <- drop.tip(Cetacea_Autocorrelated, clade$tip, trim.internal = TRUE, subtree = FALSE,root.edge = 0, rooted = is.rooted(Cetacea_Autocorrelated), interactive = FALSE)
plot(clade)
plot(rest,show.tip.label = FALSE)

tree <- rtree(10)
plot(tree)
edgelabels()
tiplabels()
tr <- reorder(tree,order="postorder")
plot(tr)
edgelabels()
tiplabels()
tr <- reorder(tree,order="pruningwise")
plot(tr)
edgelabels()

tree <- rtree(7)
par(mar=c(0,0,0,0))
plot(tree, type="cladogram", use.edge.length=FALSE, show.tip.label=FALSE)

########################
## Plot BM on branches
########################
rm(list=ls())
# setwd("/Users/paulb/Dropbox/These/Code") # Dossier de travail (Mac)
# setwd("/home/bastide/Dropbox/These/Code/Phylogenetic-EM") # Dossier de travail (Ubuntu)
library(ape)
library(phytools)
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
load("../data/Several_Trees.RData")

tree <- read.tree(text="(((E:3,D:3):1,C:4):4,(B:2,A:2):6):0.5;");
plot(tree);
tree$edge.length <- tree$edge.length*100
sig2=0.01;
## simulate evolution along each edge
set.seed(465)
X <- lapply(tree$edge.length, function(x) c(0, cumsum(rnorm(n = x, sd = sqrt(sig2)))))
## reorder the edges of the tree for pre-order traversal
cw <- reorder(tree)
## now simulate on the tree
ll <- tree$edge.length + 1
for (i in 1:nrow(cw$edge)) {
  pp <- which(cw$edge[, 2] == cw$edge[i, 1])
  if (length(pp) > 0) 
    X[[i]] <- X[[i]] + X[[pp]][ll[pp]] else X[[i]] <- X[[i]] + X[[1]][1]
}
## get the starting and ending points of each edge for plotting
H <- nodeHeights(tree)
## plot the simulation
col <- rainbow(length(X), start = 0, s = 0.5, v = 0.8)
col <- col[c(1, 3, 5, 7, 2, 4, 6, 8)]
plot(H[1, 1], X[[1]][1], ylim = range(X), xlim = range(H), xlab = "time", ylab = "phenotype")
for (i in 1:length(X)){
  lines(H[i, 1]:H[i, 2], X[[i]], col=col[i])
  points(H[i,1], X[[i]][1], col="black", pch=19)
}
## add tip labels if desired
yy <- sapply(1:length(tree$tip.label), function(x, y) which(x == y), y = cw$edge[,2])
yy <- sapply(yy, function(x, y) y[[x]][length(y[[x]])], y = X)
text(x = max(H)+20, y = yy, cw$tip.label)
text(H[1, 1] - 20, X[[1]][1], "R")
## Plot tree with colors
plot(tree, edge.color = col, edge.width = 3, no.margin = TRUE,
     y.lim = c(0.7, 5), x.lim = c(0, 850.4207))
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
arrows(lastPP$xx[6], lastPP$yy[9] + 0.1, lastPP$xx[9], lastPP$yy[9] + 0.1, lwd = 2, code =3)
text((lastPP$xx[6] + lastPP$xx[9])/2, lastPP$yy[9] + 0.2, expression(t[A][B]))
arrows(lastPP$xx[6], lastPP$yy[1] - 0.1, lastPP$xx[1], lastPP$yy[1] - 0.1, lwd = 2, code =3)
text((lastPP$xx[6] + lastPP$xx[1])/2, lastPP$yy[1] - 0.2, expression(t))
text(lastPP$xx[6] - 20, lastPP$yy[6], "R")

####################################
## Test of function simulate with no shift
####################################
set.seed(586)
tree <- rtree(20)
TreeType <- "_rtree(20)"
plot(tree); edgelabels()

variance=1
optimal.value=-3
selection.strength=2
exp.stationary <- optimal.value
var.stationary  <- variance/(2*selection.strength)
root.state <- list(random=TRUE,stationary.root=TRUE, value.root=3,exp.root=exp.stationary,var.root=var.stationary)
shifts = list(edges=c(18),values=c(6),relativeTimes=c(0))
paramsSimu <- list(variance=variance, optimal.value=optimal.value, selection.strength=selection.strength, shifts=shifts, root.state=root.state)

X1 <- simulate(tree, root.state = root.state, process = "BM", variance = variance, shifts = shifts)

X1 <- simulate(tree, root.state = root.state, process = "OU", variance=variance, optimal.value=optimal.value, selection.strength=selection.strength, shifts=shifts)

plot(tree)
X1.tips <- extract.simulate(X1,"tips","states")
X1.nodes <- extract.simulate(X1,"nodes","states")
tiplabels(pch = 19, cex = abs(X1.tips), col = ifelse(X1.tips >= 0, "orangered", "lightblue"))
nodelabels(pch = 19, cex = abs(X1.nodes), col = ifelse(X1.nodes >= 0, "orangered", "lightblue"))
plot.process("Plot_sim_OU_shift", TreeType, X1.tips, X1.nodes, tree, process="OU", paramsSimu=paramsSimu)

## Is the distribution correct ?
X.tips <- NULL
alpha <- 1
sigma2 <- 1
beta <- -3
for (i in 1:1000) {
  XX <- simulate(tree, root.state, process = "OU", variance=sigma2, optimal.value=beta, selection.strenght=alpha)
  X.tips <- rbind(X.tips,extract.simulate(XX,"tips","states"))
}
MeanEmp <- colMeans(X.tips)
MeanTh <- extract.simulate(XX,"tips","exp")  # Theoretical Mean vector
MeanTh - MeanEmp
qqnorm(X.tips[1,]); qqline(X.tips[1,])
pairs(X.tips)
CovEmp <- cov(X.tips) # Empirical variance covariance matrix
ntaxa <- attr(XX,"ntaxa")
CovTh <- matrix(NA,ncol=ntaxa,nrow=ntaxa) # Theoretical variance covariance matrix
DistNodes <- dist.nodes(tree)
Ancestors <- mrca(tree)
for (i in 1:ntaxa) {
  for (j in i:ntaxa) {
    CovTh[i,j] <- sigma2/(2*alpha)*(1-exp(-alpha*(DistNodes[ntaxa+1,Ancestors[i,j]])))*exp(-alpha*DistNodes[i,j])
    CovTh[j,i] <- CovTh[i,j]
  }
}
M <- solve(t(chol(CovTh)))
qqnorm(M%*%X.tips[1,]); qqline(M%*%X.tips[1,])

library(mvnormtest)
mshapiro.test(X.tips)

############################
## Multivariate Simulate
############################

set.seed(586)
ntaxa <- 20
tree <- rtree(ntaxa)
TreeType <- "_rtree(20)"
plot(tree); edgelabels()

p <- 3
variance <- matrix(0.5, p, p)
optimal.value <- c(-3, 5, 0)
selection.strength <- diag(3, p, p)
exp.stationary <- optimal.value
var.stationary  <- compute_stationnary_variance(variance, selection.strength)
root.state <- list(random = TRUE,
                   stationary.root = TRUE,
                   value.root = 3,
                   exp.root = exp.stationary,
                   var.root = var.stationary)
shifts = list(edges = c(18, 32),
              values=cbind(c(4, -10, 3),
                           c(-5, 5, 0)),
              relativeTimes = 0)
paramsSimu <- list(variance = variance,
                   optimal.value = optimal.value,
                   selection.strength = selection.strength,
                   shifts = shifts,
                   root.state = root.state)

X1 <- simulate(tree,
               p = p,
               root.state = root.state,
               process = "BM",
               variance = variance,
               shifts = shifts)

X1 <- simulate(tree,
               p = p,
               root.state = root.state,
               process = "OU",
               variance = variance,
               optimal.value = optimal.value,
               selection.strength = selection.strength,
               shifts = shifts)

plot(tree)
X1.tips <- extract.simulate(X1,"tips","expe")
X1.nodes <- extract.simulate(X1,"nodes","states")

par(mfrow = c(1,p), mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
for (l in 1:p){
  params <- paramsSimu
  params$shifts$values <- paramsSimu$shifts$values[l, ]
  params$optimal.value <- paramsSimu$optimal.value[l]
  plot.data.process.actual(Y.state = X1.tips[l, ],
                           phylo = tree, 
                           params = params,
                           adj.root = 2,
                           automatic_colors = TRUE,
                           margin_plot = NULL,
                           cex = 2)
}

par(mfrow = c(1,p), mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
for (l in 1:p){
  params <- paramsSimu
  params$shifts$values <- paramsSimu$shifts$values[l, ]
  params$optimal.value <- paramsSimu$optimal.value[l]
  plot(tree)
  tiplabels(text = round(X1.tips[l,], 2))
  edgelabels_home(text = round(params$shifts$values, 1), edge = params$shifts$edges)
  nodelabels(text = round(params$optimal.value, 1), node = ntaxa + 1)
}

tiplabels(pch = 19, cex = abs(X1.tips), col = ifelse(X1.tips >= 0, "orangered", "lightblue"))
nodelabels(pch = 19, cex = abs(X1.nodes), col = ifelse(X1.nodes >= 0, "orangered", "lightblue"))
plot.process("Plot_sim_OU_shift", TreeType, X1.tips, X1.nodes, tree, process="OU", paramsSimu=paramsSimu)


###########################
## Test of function simulate with shifts
###########################
tree <- rtree(10)
plot(tree)
edgelabels()

root.state <- list(random=FALSE,value.root=-1,exp.root=3,var.root=2)
shifts = list(edges=c(2),values=c(2),relativeTimes=c(0))
X1 <- simulate(tree, root.state, process = "BM", shifts, variance=1)
plot(tree)
X1.tips <- extract.simulate(X1,"tips","expectations")
X1.nodes <- extract.simulate(X1,"nodes","exp")
tiplabels(pch = 19, cex = abs(X1.tips), col = ifelse(X1.tips >= 0, "orangered", "lightblue"))
nodelabels(pch = 19, cex = abs(X1.nodes), col = ifelse(X1.nodes >= 0, "orangered", "lightblue"))

#######################################
## Test of function parsimonyNumber
#######################################
tree <- read.tree(text="(((T,T),C),(A,A));")
ntaxa <- 5
plot(tree); tiplabels(); nodelabels()
clusters=c(1,1,2,3,3)
extract.parsimonyNumber(parsimonyNumber(tree,clusters))
extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
clusters=c(1,1,1,1,1)
extract.parsimonyNumber(parsimonyNumber(tree,clusters))
extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
clusters=c(1,2,3,4,4)
extract.parsimonyNumber(parsimonyNumber(tree,clusters))
pos <- extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
for (k in 1:dim(pos)[1]){
  plot(tree, show.tip.label=FALSE)
  tiplabels(text = pos[k, 1:ntaxa])
  nodelabels(text = pos[k,-(1:ntaxa)])
}
clusters=c(1,2,2,1,1)
extract.parsimonyNumber(parsimonyNumber(tree,clusters))
extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))

tree <- rtree(10)
plot(tree); tiplabels()
clusters <- c(1,1,1,1,1,2,2,2,2,2)
extract.parsimonyNumber(parsimonyNumber(tree,clusters))
extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))

tree <- read.tree(text="((((T,T),C),(A,A)),(((T,T),C),(A,A)));")
plot(tree); tiplabels(); nodelabels()
clusters=c(1,1,2,3,3,1,1,2,3,3)
extract.parsimonyNumber(parsimonyNumber(tree,clusters))
extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))

tree <- read.tree(text="(((T,T),C),C);")
ntaxa <- 4
plot(tree); tiplabels(); nodelabels()
clusters=c(1,2,3,3)
extract.parsimonyNumber(parsimonyNumber(tree,clusters))
extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))

#######################################
## Test of function partitionsNumber - Binary
#######################################
set.seed(23); phy <- rtree(25); plot(phy); nodelabels(); tiplabels()
extract.partitionsNumber(partitionsNumber(phy,2))

n <- 3

set.seed(456); phy <- rtree(2^n); plot(phy, show.tip.label = FALSE); nodelabels()
extract.partitionsNumber(partitionsNumber(phy,10), npart=1:10)

SymTree <- rtree.sym(n); plot(SymTree,show.tip.label = FALSE); nodelabels()
extract.partitionsNumber(partitionsNumber(SymTree,10), npart=1:10)

CombTree <- rtree.comb(20); plot(CombTree,show.tip.label = FALSE); nodelabels()
extract.partitionsNumber(partitionsNumber(CombTree,3))

n <- 50
k <- 30
N <- 2:n
CombTree <- rtree.comb(n)
choose(2*N-1-k,k-1)
extract.partitionsNumber(partitionsNumber(CombTree,k),node=rev((n+1):(2*n-1)))
choose(2*N-k,k-1)
extract.partitionsNumber(partitionsNumber(CombTree,k),node=rev((n+1):(2*n-1)),marqued=TRUE)

randomTree <- rtree(n); plot(randomTree)
extract.partitionsNumber(partitionsNumber(randomTree,k))
extract.partitionsNumber(partitionsNumber(randomTree,k), marqued=TRUE)

#######################################
## Test of function partitionsNumber - General
#######################################

K <- 5
tree <- read.tree(text="(A,(A,A,A));"); plot(tree, show.tip.label=FALSE)
extract.partitionsNumber(partitionsNumber(tree,K), npart=1:K)
CombTree <- rtree.comb(4); plot(CombTree,show.tip.label = FALSE); nodelabels()
extract.partitionsNumber(partitionsNumber(CombTree,K), npart=1:K)

K <- 8
tree <- read.tree(text="(A,(A,(A,A,A),A,A));"); plot(tree, show.tip.label=FALSE)
extract.partitionsNumber(partitionsNumber(tree,K), npart=1:K)
CombTree <- rtree.comb(7); plot(CombTree,show.tip.label = FALSE); nodelabels()
extract.partitionsNumber(partitionsNumber(CombTree,K), npart=1:K)


K <- 3
set.seed(654)
tree <- rtree(10); plot(tree, show.tip.label=FALSE, type="cladogram");
tree <- di2multi(tree,0.3); plot(tree, show.tip.label=FALSE, type="cladogram");
extract.partitionsNumber(partitionsNumber(tree,K))
choose(length(tree$edge.length)-2,K-1)

K <- 4
n <- 13
tree <- read.tree(text="(A,(A,A,A),(A,(A,(A,A,A),A),(A,A,A)));"); tree$root.edge <- 1; plot(tree, show.tip.label=FALSE); nodelabels()
extract.partitionsNumber(partitionsNumber(tree,K), npart=K, node=c(14,16,17,18))
choose((3*c(n,n-4,n-8,n-10)-3)/2,K-1)
CombTree <- rtree.comb(9); plot(CombTree,show.tip.label = FALSE); nodelabels()
extract.partitionsNumber(partitionsNumber(CombTree,K), npart=K, node=c(10,14,16))

K <- 9
tr1 <- read.tree(text="((A,A,A),(A,A,A),(A,A,A));"); tr1$root.edge <- 1; plot(tr1, show.tip.label=FALSE);
tr2 <- read.tree(text="(A,A,(A,A,(A,A,(A,A,A))));"); tr2$root.edge <- 1; plot(tr2, show.tip.label=FALSE);
extract.partitionsNumber(partitionsNumber(tr1,K), npart=1:K)
extract.partitionsNumber(partitionsNumber(tr2,K), npart=1:K)

#######################################
## Test of function incidence.matrix
#######################################
set.seed(249); phy <- rtree(4); plot(phy); edgelabels(); tiplabels(); nodelabels()
T <- incidence.matrix(phy); T
U <- incidence.matrix.full(phy); U
shifts <- list(edges=c(2,3),values=c(5,5),relativeTimes=c(0,0))
delta <- shifts.list_to_vector(phy, shifts); delta
T%*%delta
U%*%delta

set.seed(249); phy <- rtree(20); plot(phy); edgelabels(); tiplabels(); nodelabels()
T <- incidence.matrix(phy); T
U <- incidence.matrix.full(phy); U
shifts <- list(edges=c(12,15,17,21),values=c(5,5,5,5),relativeTimes=c(0,0,0,0))
delta <- shifts.list_to_vector(phy, shifts); delta
T%*%delta
U%*%delta

phy <- rtree(20);
shifts <- list(edges=c(12,15,17,21),values=c(1,2,3,4),relativeTimes=c(0,0,0,0))
delta <- shifts.list_to_vector(phy, shifts);
all.equal(shifts, shifts.vector_to_list(delta))


#######################################
## Test of allocate shifts
#######################################
set.seed(249);
ntaxa <- 40
phy <- rtree(ntaxa);
shifts <- list(edges=c(12,15,17,21),values=c(5,5,5,5),relativeTimes=c(0,0,0,0))
regimes <- allocate_regimes_from_shifts(phy, shifts$edges)
x11()
plot(phy, show.tip.label=FALSE); edgelabels(ed  = shifts$edges, bg = 1:4); tiplabels(bg = regimes[1:ntaxa]); nodelabels(bg = regimes[-(1:ntaxa)])

shift_edges <- allocate_shifts_from_regimes(phy, regimes)
all.equal(shift_edges, shifts$edges)
plot(phy, show.tip.label=FALSE); edgelabels(ed  = shift_edges, bg = 1:4); tiplabels(bg = regimes[1:ntaxa]); nodelabels(bg = regimes[-(1:ntaxa)])

#######################################
## Test of compute betas
#######################################
set.seed(249);
ntaxa <- 40
phy <- rtree(ntaxa);
shifts <- list(edges=c(12,15,17,21), values=c(1,-1,2,-2), relativeTimes=c(0,0,0,0))
betas <- compute_betas(phy, 0, shifts)
x11()
plot(phy, show.tip.label=FALSE); edgelabels(ed  = shifts$edges, text = shifts$values); tiplabels(text = betas[1:ntaxa]); nodelabels(text = betas[-(1:ntaxa)])

shifts_re <- compute_shifts_from_betas(phy, betas)
all.equal(shifts, shifts_re)
plot(phy, show.tip.label=FALSE); edgelabels(ed  = shifts_re$edges, text = shifts_re$values); tiplabels(text = betas[1:ntaxa]); nodelabels(text = betas[-(1:ntaxa)])

#######################################
## Initialization - lasso
#######################################

## Tree
#set.seed(152);
#phy <- rtree(20);
phy <- reorder(Cetacea_Autocorrelated, order="cladewise")
plot(phy); edgelabels(); nodelabels(); tiplabels()

## Matrix T
T <- incidence.matrix(phy);

## Simulate a process on the tree
set.seed(456)
# Parameters
root.state <- list(random=TRUE, stationary.root=TRUE, value.root=0,exp.root=0,var.root=1)
shifts <- list(edges=c(90,77),values=c(5,5),relativeTimes=c(0,0))
process <- "BM"
variance <- 1
selection.strength <- 2
optimal.value <- 0
# Test root.state
root.state <- test.root.state(process=process, root.state=root.state, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
# Display parameters
paramsSimu <- list(root.state=root.state, shifts=shifts, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
paramsSimu
# Simulate
XX <- simulate(phy, root.state, process = process, variance=variance, shifts=shifts, selection.strength=selection.strength, optimal.value=optimal.value)
Y_data <- extract.simulate(XX, what="states", where="tips")
Z.nodes <- extract.simulate(XX,"nodes","states"); Z.nodes

## Nombre de ruptures pour l'initialisation
K <- 2

## Lasso simple
# Fit
fit.simple <- glmnet(T, Y_data, alpha=1)
plot(fit.simple); print(fit.simple)
#coef(fit); fit$df
# Trouver un lambda qui donne le bon nombre de ruptures
delta.simple <- coef(fit.simple, s=fit.simple$lambda[min(which(fit.simple$df==K))])
E0.simple <- delta.simple[1]; # Intercept
delta.simple <- delta.simple[-1];
shifts.init.simple <- shifts.vector_to_list(delta.simple)
# Gauss Lasso
projection.simple <- which(delta.simple != 0)
fit.simple.gauss <- lm(Y_data ~ T[,projection.simple])
delta.simple.gauss <- rep(0, dim(T)[2])
EO.simple.gauss <- coef(fit.simple.gauss)[1]
delta.simple.gauss[projection.simple] <- coef(fit.simple.gauss)[-1]
shifts.init.simple.gauss <- shifts.vector_to_list(delta.simple.gauss)

## Lasso avec Sigma_YY
# Initialiser Sigma avec les paramètres par défauts
params.default <- init.EM.default.OU()
times_shared <- compute_times_ca(phy)
distances_phylo <- compute_dist_phy(phy)
Sigma <- compute_variance_covariance.OU(times_shared=times_shared, distances_phylo=distances_phylo, params_old=params.default)
Sigma_YY <- extract.variance_covariance(Sigma, what="YY")
# Cholesky
Sig_chol <- chol(Sigma_YY)
Sig_chol_inv <- t(solve(Sig_chol)) # Sigma_YY_inv = t(Sig_chol_inv)%*%Sig_chol_inv
# Transform Y_data and T
Tp <- Sig_chol_inv%*%T
Yp <- Sig_chol_inv%*%Y_data
# Fit
fit.sig <- glmnet(Tp, Yp, alpha=1)
plot(fit.sig);
# Trouver un lambda qui donne le bon nombre de ruptures
delta.sig <- coef(fit.sig, s=fit.sig$lambda[min(which(fit.sig$df==K))])
E0.sig <- delta.sig[1]; # Intercept
delta.sig <- delta.sig[-1];
shifts.init.sig <- shifts.vector_to_list(delta.sig);
# Gauss lasso
projection.sig <- which(delta.sig != 0)
fit.sig.gauss <- lm(Yp ~ Tp[,projection.sig])
delta.sig.gauss <- rep(0, dim(T)[2])
EO.sig.gauss <- coef(fit.sig.gauss)[1]
delta.sig.gauss[projection.sig] <- coef(fit.sig.gauss)[-1]
shifts.init.sig.gauss <- shifts.vector_to_list(delta.sig.gauss);

## Compute gamma^2 (hyp : OU + root is stationnary)
Lineages <- rowSums(T[,projection.sig]) > 0 # All lineages that have at least one shift in their history
gamma2 <- var(Y_data[!Lineages]) # Variance on all the other tips, that should be gamma^2 if stationnary.
# Estimation of sigma2 if alpha is known
sigma2 <- gamma2 * 2 * selection.strength

## Compute sigma^2 (Hyp :BM)
Lineages <- rowSums(T[,projection.sig]) > 0
depths <- node.depth.edgelength(phy)[1:length(phy$tip.label)]
sigma2 <- var(Y_data[!Lineages] / sqrt(depths[!Lineages]))

cbind(Simple.Lasso = unlist(shifts.init.simple), Simple.Lasso.Gauss = unlist(shifts.init.simple.gauss), Sigma.Lasso = unlist(shifts.init.sig), Sigma.Lasso.Gauss = unlist(shifts.init.sig.gauss))

#######################################
## Test of EM - BM
#######################################
set.seed(123)
ntips <- 50
tree <- rtree(ntips)
TreeType <- paste("_tree=rtree(",ntips,")",sep="")
# tree <- reorder(Cetacea_Autocorrelated, order="cladewise")
# TreeType <- "_Cetacea_Meredith-Autocorrelated"
pdf(paste("Results/Miscellaneous_Evals/Plot_tree",TreeType,".pdf",sep=""), height=20,width=10)
plot(tree, show.tip.label = FALSE); edgelabels();axisPhylo()
dev.off()

# Simulate a process on the tree
set.seed(456)
root.state <- list(random=TRUE,value.root=0,exp.root=0,var.root=1)
shifts <- list(edges=c(65),values=c(5),relativeTimes=c(0))
process <- "BM"
variance <- 1
root.state <- test.root.state(process=process, root.state=root.state, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
paramsSimu <- list(root.state=root.state, shifts=shifts, variance=variance)
XX <- simulate(tree, root.state, process = process, variance=variance, shifts=shifts)
Y_data <- extract.simulate(XX, what="states", where="tips")
Z.nodes <- extract.simulate(XX,"nodes","states"); Z.nodes
plot.process("Plot_simulation", TreeType, Y_data, Z.nodes, tree, process = process, paramsSimu=paramsSimu, directory="Results/Miscellaneous_Evals/")
save.process("Data_simulated",TreeType, XX, process, paramsSimu, directory="Results/Miscellaneous_Evals/")

# Estimate parameters from the data
tol=list(variance=10^(-4), value.root=10^(-4), exp.root=10^(-4), var.root=10^(-4))
params_algo_EM <- list(process=process, tol=tol, method.variance="simple", method.init="default", nbr_of_shifts=1)
results_estim_EM <- estimateEM(phylo=tree, Y_data=Y_data, tol=params_algo_EM$tol, process="BM", method.variance=params_algo_EM$method.variance, method.init=params_algo_EM$method.init, Nbr_It_Max=1000, nbr_of_shifts=params_algo_EM$nbr_of_shifts, random=TRUE)
results_estim_EM
params_estim_EM <- results_estim_EM$params
Z_reconstructed <- results_estim_EM$ReconstructedNodesStates

# Plot the reconstructed states
plot.process("Plot_reconstructed", TreeType, Y_data, Z_reconstructed, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM, estimate=TRUE, params_algo_EM=params_algo_EM, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstructed",TreeType, XX, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM, estimate=TRUE, directory="Results/Miscellaneous_Evals/")

# Simulate a process with estimated parameters
XX_sim <- simulate(tree, root.state=params_estim_EM$root.state, process = "BM", variance=params_estim_EM$variance)
plot(tree)
Y_sim <- extract.simulate(XX_sim, what="states", where="tips")
Z_sim <- extract.simulate(XX_sim,"nodes","states")
tiplabels(pch = 19, cex = abs(Y_sim)/mean(abs(Y_sim)), col = ifelse(Y_sim >= 0, "orangered", "lightblue"))
nodelabels(pch = 19, cex = abs(Z_sim)/mean(abs(Z_sim)), col = ifelse(Z_sim >= 0, "orangered", "lightblue"))

#######################################
## Test of EM - OU
#######################################
# set.seed(123)
# ntips <- 10
# tree <- rtree(ntips)
# TreeType <- paste("_tree=rtree(",ntips,")",sep="")
tree <- reorder(Cetacea_Autocorrelated, order="cladewise")
TreeType <- "_Cetacea_Meredith-Autocor"
pdf(paste("Results/Miscellaneous_Evals/Plot_tree",TreeType,".pdf",sep=""), height=20,width=10)
plot(tree, show.tip.label = FALSE); edgelabels();axisPhylo()
dev.off()
plot(tree); edgelabels(); nodelabels(); tiplabels()

## Simulate a process on the tree
set.seed(456)
  # Parameters
root.state <- list(random=TRUE, stationary.root=TRUE, value.root=0,exp.root=0,var.root=1)
#shifts <- list(edges=c(90),values=c(5),relativeTimes=c(0))
shifts <- list(edges=c(5),values=c(3),relativeTimes=c(0))
process <- "OU"
variance <- 2
selection.strength <- 2
optimal.value <- 0
  # Test root.state
root.state <- test.root.state(process=process, root.state=root.state, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
  # Display parameters
paramsSimu <- list(root.state=root.state, shifts=shifts, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
paramsSimu
  # Simulate
XX <- simulate(tree, root.state, process = process, variance=variance, shifts=shifts, selection.strength=selection.strength, optimal.value=optimal.value)
  # Plot and save
Y_data <- extract.simulate(XX, what="states", where="tips")
Z.nodes <- extract.simulate(XX,"nodes","states"); Z.nodes
plot.process("Plot_sim", TreeType, Y_data, Z.nodes, tree, process = process, paramsSimu=paramsSimu, directory="Results/Miscellaneous_Evals/")
save.process("Data_sim",TreeType, XX, process, paramsSimu, directory="Results/Miscellaneous_Evals/")

## Estimate parameters from the data
  # Parameters for the algorithm
tol=list(variance=10^(-4), value.root=10^(-4), exp.root=10^(-4), var.root=10^(-4))
params_algo_EM <- list(process="OU", tol=tol, method.variance="simple", method.init="default", nbr_of_shifts=1, selection.strength=selection.strength)
  # EM
system.time(results_estim_EM <- estimateEM(phylo=tree, Y_data=Y_data, tol=params_algo_EM$tol, process=params_algo_EM$process, method.variance=params_algo_EM$method.variance, method.init=params_algo_EM$method.init, Nbr_It_Max=500, nbr_of_shifts=params_algo_EM$nbr_of_shifts, specialCase=TRUE, selection.strength=params_algo_EM$selection.strength))
  # Display results
results_estim_EM
params_init <- results_estim_EM$params_init
params_init <- replaceInList(params_init, function(x) if(is.null(x))NA else x)
res <- cbind(params_simu = unlist(paramsSimu)[names(unlist(results_estim_EM$params))], params_init = unlist(params_init), params_old = unlist(results_estim_EM$params_old), params = unlist(results_estim_EM$params))
res

params_estim_EM <- results_estim_EM$params
Z_reconstructed <- results_estim_EM$ReconstructedNodesStates

## Plot the reconstructed states and save
plot.process("Plot_reconstruct", TreeType, Y_data, Z_reconstructed, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM, estimate=TRUE, params_algo_EM=params_algo_EM, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstruct",TreeType, XX, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
write.table.results("Table_Results",TreeType, res, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM, estimate=TRUE, directory="Results/Miscellaneous_Evals/")


#######################################
## Test of EM - OU with lasso initialization
#######################################
rm(list=ls())
#setwd("/Users/paulb/Dropbox/These/Code") # Dossier de travail (Mac)
setwd("/home/bastide/Dropbox/These/Code") # Dossier de travail (Ubuntu)
library(ape)
library(glmnet) # For Lasso initialization
source("R/functions.R")

load("data/Several_Trees.RData")
# set.seed(123)
# ntips <- 10
# tree <- rtree(ntips)
# TreeType <- paste("_tree=rtree(",ntips,")",sep="")
tree <- reorder(Cetacea_Autocorrelated, order="cladewise")
TreeType <- "_Cetacea_Meredith-Autocor"
pdf(paste("Results/Miscellaneous_Evals/Plot_tree",TreeType,".pdf",sep=""), height=20,width=10)
plot(tree, show.tip.label = FALSE); edgelabels();axisPhylo()
dev.off()
plot(tree, show.tip.label=FALSE); #edgelabels(); nodelabels(); tiplabels()

## Simulate a process on the tree
set.seed(456)
# Parameters
root.state <- list(random=TRUE, stationary.root=TRUE, value.root=0,exp.root=0,var.root=1)
shifts <- list(edges=c(5),values=c(5),relativeTimes=c(0))
#shifts <- list(edges=c(137, 138),values=c(2, -4),relativeTimes=c(0, 0))
process <- "OU"
variance <- 2
selection.strength <- 2
optimal.value <- 0
# Test root.state
root.state <- test.root.state(process=process, root.state=root.state, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
# Display parameters
paramsSimu <- list(root.state=root.state, shifts=shifts, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
paramsSimu
# Simulate
XX <- simulate(tree, root.state, process = process, variance=variance, shifts=shifts, selection.strength=selection.strength, optimal.value=optimal.value)
# Plot and save
Y_data <- extract.simulate(XX, what="states", where="tips")
Z.nodes <- extract.simulate(XX,"nodes","states"); Z.nodes
plot.process("Plot_sim", TreeType, Y_data, Z.nodes, tree, process = process, paramsSimu=paramsSimu, directory="Results/Miscellaneous_Evals/")
save.process("Data_sim",TreeType, XX, process, paramsSimu, directory="Results/Miscellaneous_Evals/")
rm(root.state, shifts, process, variance, selection.strength, optimal.value, XX, Z.nodes)

## Estimate parameters from the data
# Parameters for the algorithm
tol=list(variance=10^(-4), value.root=10^(-4), exp.root=10^(-4), var.root=10^(-4))
# Default init
params_algo_EM.default <- list(process="OU", tol=tol, method.variance="simple", method.init="default", nbr_of_shifts=1, selection.strength=paramsSimu$selection.strength)
# Lasso init
params_algo_EM.lasso <- params_algo_EM.default
params_algo_EM.lasso$method.init <- "lasso" 
# EM default init
time.default <- system.time(
  results_estim_EM.default <- estimateEM(
    phylo=tree, 
    Y_data=Y_data, 
    tol=params_algo_EM.default$tol, 
    process=params_algo_EM.default$process, 
    method.variance=params_algo_EM.default$method.variance, 
    method.init=params_algo_EM.default$method.init, 
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.default$nbr_of_shifts, 
    alpha_known=TRUE, 
    known.selection.strength=params_algo_EM.default$selection.strength
    )
  )
# EM lasso init
time.lasso <- system.time(
  results_estim_EM.lasso <- estimateEM(
    phylo=tree, 
    Y_data=Y_data, 
    tol=params_algo_EM.lasso$tol, 
    process=params_algo_EM.lasso$process,
    method.variance=params_algo_EM.lasso$method.variance, 
    method.init=params_algo_EM.lasso$method.init,
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.lasso$nbr_of_shifts,
    alpha_known=TRUE, 
    known.selection.strength=params_algo_EM.lasso$selection.strength
  )
)
# Display results
#results_estim_EM.default; results_estim_EM.lasso
params_init.default <- results_estim_EM.default$params_init
params_init.default <- replaceInList(params_init.default, function(x) if(is.null(x))rep(NA,params_algo_EM.default$nbr_of_shifts) else x)
res <- cbind(params_simu = unlist(paramsSimu)[names(unlist(results_estim_EM.default$params))],
             params_init.default = unlist(params_init.default),
             params_old.default = unlist(results_estim_EM.default$params_old),
             params.default = unlist(results_estim_EM.default$params),
             params_init.lasso = unlist(results_estim_EM.lasso$params_init),
             params_old.lasso = unlist(results_estim_EM.lasso$params_old),
             params.lasso = unlist(results_estim_EM.lasso$params))
res

params_estim_EM.default <- results_estim_EM.default$params
Z_reconstructed.default <- results_estim_EM.default$ReconstructedNodesStates
params_estim_EM.lasso <- results_estim_EM.lasso$params
Z_reconstructed.lasso <- results_estim_EM.lasso$ReconstructedNodesStates

## Plot the reconstructed states and save
process <- "OU"
# Default init
plot.process("Plot_reconstruct", TreeType, Y_data, Z_reconstructed.default, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, params_algo_EM=params_algo_EM.default, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstruct",TreeType, results_estim_EM.default, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
write.table.results("Table_Results",TreeType, res, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
# Lasso init
plot.process("Plot_reconstruct", TreeType, Y_data, Z_reconstructed.lasso, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, params_algo_EM=params_algo_EM.lasso, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstruct",TreeType, results_estim_EM.lasso, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
write.table.results("Table_Results",TreeType, res, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, directory="Results/Miscellaneous_Evals/")


#######################################
## Test of EM - BM with lasso initialization
#######################################
rm(list=ls())
#setwd("/Users/paulb/Dropbox/These/Code") # Dossier de travail (Mac)
setwd("/home/bastide/Dropbox/These/Code") # Dossier de travail (Ubuntu)
library(ape)
library(glmnet) # For Lasso initialization
source("R/functions.R")

load("data/Several_Trees.RData")
# set.seed(123)
# ntips <- 10
# tree <- rtree(ntips)
# TreeType <- paste("_tree=rtree(",ntips,")",sep="")
tree <- reorder(Cetacea_Autocorrelated, order="cladewise")
TreeType <- "_Cetacea_Meredith-Autocor"
pdf(paste("Results/Miscellaneous_Evals/Plot_tree",TreeType,".pdf",sep=""), height=20,width=10)
plot(tree, show.tip.label = FALSE); edgelabels();axisPhylo()
dev.off()
plot(tree, show.tip.label=FALSE); #edgelabels(); nodelabels(); tiplabels()

## Simulate a process on the tree
set.seed(456)
# Parameters
root.state <- list(random=TRUE, stationary.root=TRUE, value.root=0,exp.root=0,var.root=1)
#shifts <- list(edges=c(90),values=c(5),relativeTimes=c(0))
shifts <- list(edges=c(137, 138),values=c(2, -4),relativeTimes=c(0, 0))
process <- "BM"
variance <- 1
# Test root.state
root.state <- test.root.state(process=process, root.state=root.state, variance=variance)
# Display parameters
paramsSimu <- list(root.state=root.state, shifts=shifts, variance=variance)
paramsSimu
# Simulate
XX <- simulate(tree, root.state, process = process, variance=variance, shifts=shifts)
# Plot and save
Y_data <- extract.simulate(XX, what="states", where="tips")
Z.nodes <- extract.simulate(XX,"nodes","states"); Z.nodes
plot.process("Plot_sim", TreeType, Y_data, Z.nodes, tree, process = process, paramsSimu=paramsSimu, directory="Results/Miscellaneous_Evals/")
save.process("Data_sim",TreeType, XX, process, paramsSimu, directory="Results/Miscellaneous_Evals/")
rm(root.state, shifts, process, variance, XX, Z.nodes)

## Estimate parameters from the data
# Parameters for the algorithm
tol=list(variance=10^(-4), value.root=10^(-4), exp.root=10^(-4), var.root=10^(-4))
# Default init
params_algo_EM.default <- list(process="BM", tol=tol, method.variance="simple", method.init="default", nbr_of_shifts=2)
# Lasso init
params_algo_EM.lasso <- params_algo_EM.default
params_algo_EM.lasso$method.init <- "lasso" 
# EM default init
time.default <- system.time(
  results_estim_EM.default <- estimateEM(
    phylo=tree, 
    Y_data=Y_data, 
    tol=params_algo_EM.default$tol, 
    process=params_algo_EM.default$process, 
    method.variance=params_algo_EM.default$method.variance, 
    method.init=params_algo_EM.default$method.init, 
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.default$nbr_of_shifts, 
    specialCase=TRUE, 
    selection.strength=params_algo_EM.default$selection.strength
  )
)
# EM lasso init
time.lasso <- system.time(
  results_estim_EM.lasso <- estimateEM(
    phylo=tree, 
    Y_data=Y_data, 
    tol=params_algo_EM.lasso$tol, 
    process=params_algo_EM.lasso$process,
    method.variance=params_algo_EM.lasso$method.variance, 
    method.init=params_algo_EM.lasso$method.init,
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.lasso$nbr_of_shifts,
    specialCase=TRUE, 
    selection.strength=params_algo_EM.lasso$selection.strength
  )
)
# Display results
#results_estim_EM.default; results_estim_EM.lasso
params_init.default <- results_estim_EM.default$params_init
params_init.default <- replaceInList(params_init.default, function(x) if(is.null(x))rep(NA,params_algo_EM.default$nbr_of_shifts) else x)
res <- cbind(params_simu = unlist(paramsSimu)[names(unlist(results_estim_EM.default$params))],
             params_init.default = unlist(params_init.default),
             params_old.default = unlist(results_estim_EM.default$params_old),
             params.default = unlist(results_estim_EM.default$params),
             params_init.lasso = unlist(results_estim_EM.lasso$params_init),
             params_old.lasso = unlist(results_estim_EM.lasso$params_old),
             params.lasso = unlist(results_estim_EM.lasso$params))
res

params_estim_EM.default <- results_estim_EM.default$params
Z_reconstructed.default <- results_estim_EM.default$ReconstructedNodesStates
params_estim_EM.lasso <- results_estim_EM.lasso$params
Z_reconstructed.lasso <- results_estim_EM.lasso$ReconstructedNodesStates

## Plot the reconstructed states and save
process <- "BM"
# Default init
plot.process("Plot_reconstruct", TreeType, Y_data, Z_reconstructed.default, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, params_algo_EM=params_algo_EM.default, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstruct",TreeType, results_estim_EM.default, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
write.table.results("Table_Results",TreeType, res, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
# Lasso init
plot.process("Plot_reconstruct", TreeType, Y_data, Z_reconstructed.lasso, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, params_algo_EM=params_algo_EM.lasso, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstruct",TreeType, results_estim_EM.lasso, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
write.table.results("Table_Results",TreeType, res, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, directory="Results/Miscellaneous_Evals/")


#######################################
## Test of EM - OU with lasso initialization - alpha unknown
#######################################
rm(list=ls())
#setwd("/Users/paulb/Dropbox/These/Code") # Dossier de travail (Mac)
WD <- "/home/bastide/Dropbox/These/Code"
setwd(WD) # Dossier de travail (Ubuntu)
library(ape)
library(glmnet) # For Lasso initialization
library(nlme) # For second derivative computation
library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
source("R/simulate.R")
source("R/estimateEM.R")
source("R/init_EM.R")
source("R/E_step.R")
source("R/M_step.R")
source("R/shutoff.R")
source("R/generic_functions.R")
source("R/shifts_manipulations.R")
source("R/plot_functions.R")

load("data/Several_Trees.RData")
# set.seed(123)
# ntips <- 10
# tree <- rtree(ntips)
# TreeType <- paste("_tree=rtree(",ntips,")",sep="")
tree <- reorder(Cetacea_Autocorrelated, order="cladewise")
TreeType <- "_Cetacea_Meredith-Autocor"
pdf(paste("Results/Miscellaneous_Evals/Plot_tree",TreeType,".pdf",sep=""), height=10,width=20)
plot(tree, show.tip.label = FALSE, font = 3, cex=1);axisPhylo(cex.axis=2) #edgelabels();
dev.off()
plot(tree, show.tip.label=FALSE); #edgelabels(); nodelabels(); tiplabels()

## Simulate a process on the tree
set.seed(456)
# Parameters
root.state <- list(random=TRUE, stationary.root=TRUE, value.root=0,exp.root=0,var.root=1)
shifts <- list(edges=c(5, 90),values=c(2, -2),relativeTimes=c(0,0))
#shifts <- list(edges=c(137, 138),values=c(2, -4),relativeTimes=c(0, 0))
process <- "OU"
variance <- 1
selection.strength <- 5
optimal.value <- 0
# Test root.state
root.state <- test.root.state(process=process, root.state=root.state, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
# Display parameters
paramsSimu <- list(root.state=root.state, shifts=shifts, variance=variance, selection.strength=selection.strength, optimal.value=optimal.value)
paramsSimu
# Simulate
XX <- simulate(tree, root.state, process = process, variance=variance, shifts=shifts, selection.strength=selection.strength, optimal.value=optimal.value)
# Plot and save
Y_data <- extract.simulate(XX, what="states", where="tips")
Z.nodes <- extract.simulate(XX,"nodes","states"); Z.nodes
plot.process("Plot_sim", TreeType, Y_data, Z.nodes, tree, process = process, paramsSimu=paramsSimu, directory="Results/Miscellaneous_Evals/")
save.process("Data_sim",TreeType, XX, process, paramsSimu, directory="Results/Miscellaneous_Evals/")
rm(root.state, shifts, process, variance, selection.strength, optimal.value, XX, Z.nodes)

## Estimate parameters from the data
# Parameters for the algorithm
tol=list(variance=10^(-4), value.root=10^(-4), exp.root=10^(-4), var.root=10^(-4), selection.strength=10^(-3))
# Default init
params_algo_EM.default <- list(process="OU", tol=tol, method.variance="simple", method.init="default", method.init.alpha="estimation", nbr_of_shifts=2, alpha_known=FALSE)
# Lasso init
params_algo_EM.lasso <- params_algo_EM.default
params_algo_EM.lasso$method.init <- "lasso" 
# EM default init
time.default <- system.time(
  results_estim_EM.default <- estimateEM(
    phylo=tree, 
    Y_data=Y_data, 
    tol=params_algo_EM.default$tol, 
    process=params_algo_EM.default$process, 
    method.variance=params_algo_EM.default$method.variance, 
    method.init=params_algo_EM.default$method.init,
    method.init.alpha=params_algo_EM.default$method.init.alpha,
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.default$nbr_of_shifts, 
    alpha_known=params_algo_EM.default$alpha_known
  )
)
# EM lasso init
time.lasso <- system.time(
  results_estim_EM.lasso <- estimateEM(
    phylo=tree, 
    Y_data=Y_data, 
    tol=params_algo_EM.lasso$tol, 
    process=params_algo_EM.lasso$process,
    method.variance=params_algo_EM.lasso$method.variance, 
    method.init=params_algo_EM.lasso$method.init,
    method.init.alpha=params_algo_EM.lasso$method.init.alpha,
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.lasso$nbr_of_shifts,
    alpha_known=params_algo_EM.lasso$alpha_known
  )
)
plot.history.OU.stationnary(results_estim_EM.default$params_history, paramsSimu, PATH=paste(WD, "/Results/Miscellaneous_Evals/", sep=""), params_algo_EM.default)
plot.history.OU.stationnary(results_estim_EM.lasso$params_history, paramsSimu, PATH=paste(WD, "/Results/Miscellaneous_Evals/", sep=""), params_algo_EM.lasso)
# Display results
#results_estim_EM.default; results_estim_EM.lasso
params_init.default <- results_estim_EM.default$params_init
params_init.default <- replaceInList(params_init.default, function(x) if(is.null(x))rep(NA,params_algo_EM.default$nbr_of_shifts) else x)
res <- cbind(params_simu = unlist(paramsSimu)[names(unlist(results_estim_EM.default$params))],
             params_init.default = unlist(params_init.default),
             params_old.default = unlist(results_estim_EM.default$params_old),
             params.default = unlist(results_estim_EM.default$params),
             params_init.lasso = unlist(results_estim_EM.lasso$params_init),
             params_old.lasso = unlist(results_estim_EM.lasso$params_old),
             params.lasso = unlist(results_estim_EM.lasso$params))
res

params_estim_EM.default <- results_estim_EM.default$params
Z_reconstructed.default <- results_estim_EM.default$ReconstructedNodesStates
params_estim_EM.lasso <- results_estim_EM.lasso$params
Z_reconstructed.lasso <- results_estim_EM.lasso$ReconstructedNodesStates

## Plot the reconstructed states and save
process <- "OU"
# Default init
plot.process("Plot_reconstruct", TreeType, Y_data, Z_reconstructed.default, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, params_algo_EM=params_algo_EM.default, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstruct",TreeType, results_estim_EM.default, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
write.table.results("Table_Results",TreeType, res, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.default, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
# Lasso init
plot.process("Plot_reconstruct", TreeType, Y_data, Z_reconstructed.lasso, tree, process = process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, params_algo_EM=params_algo_EM.lasso, directory="Results/Miscellaneous_Evals/")
save.process("Data_reconstruct",TreeType, results_estim_EM.lasso, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, directory="Results/Miscellaneous_Evals/")
write.table.results("Table_Results",TreeType, res, process, paramsSimu=paramsSimu, paramsEstimate=params_estim_EM.lasso, estimate=TRUE, directory="Results/Miscellaneous_Evals/")


#######################################
## Test of EM - Butler King
#######################################
rm(list=ls())
WD <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
# WD <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
setwd(WD)
PATH <- paste(WD, "/Results/Butler_King/", sep="")
library(ape)
library(glmnet) # For Lasso initialization
#library(nlme) # For second derivative computation
#library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
source("R/functions.R")

library(ouch)
source("R/convert.R")
# Anolis bimaculatus lizard size data
data(bimac)
bimac_ouch_tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
bimac_ape_tree <- convert.pmc(bimac_ouch_tree)
plot(bimac_ape_tree); edgelabels()
Y_data_bimac <- bimac$size[23:45]
names(Y_data_bimac) <- bimac$species[23:45]
Y_data_bimac <- log(unname(Y_data_bimac[bimac_ape_tree$tip.label]))

tol=list(variance=10^(-4), value.root=10^(-4), exp.root=10^(-4), var.root=10^(-4), selection.strength=10^(-3))
# Default init
params_algo_EM.default <- list(process="OU", tol=tol, method.variance="simple", method.init="default", method.init.alpha="estimation", nbr_of_shifts=2, alpha_known=FALSE, known.selection.strength=2.5)
# Lasso init
params_algo_EM.lasso <- params_algo_EM.default
params_algo_EM.lasso$method.init <- "lasso" 
# EM default init
time.default <- system.time(
  results_estim_EM.default <- estimateEM(
    phylo=bimac_ape_tree, 
    Y_data=Y_data_bimac, 
    tol=params_algo_EM.default$tol, 
    process=params_algo_EM.default$process, 
    method.variance=params_algo_EM.default$method.variance, 
    method.init=params_algo_EM.default$method.init,
    method.init.alpha=params_algo_EM.default$method.init.alpha,
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.default$nbr_of_shifts, 
    alpha_known=params_algo_EM.default$alpha_known,
    known.selection.strength=params_algo_EM.default$known.selection.strength
  )
)
# EM lasso init
time.lasso <- system.time(
  results_estim_EM.lasso <- estimateEM(
    phylo=bimac_ape_tree, 
    Y_data=Y_data_bimac, 
    tol=params_algo_EM.lasso$tol, 
    process=params_algo_EM.lasso$process,
    method.variance=params_algo_EM.lasso$method.variance, 
    method.init=params_algo_EM.lasso$method.init,
    method.init.alpha=params_algo_EM.lasso$method.init.alpha,
    Nbr_It_Max=500, 
    nbr_of_shifts=params_algo_EM.lasso$nbr_of_shifts,
    alpha_known=params_algo_EM.lasso$alpha_known,
    known.selection.strength=params_algo_EM.default$known.selection.strength
  )
)

# Infered parameters by Butler and King
root.stateBK <- list(random=TRUE, stationary.root=TRUE, value.root=NA,exp.root=0.86, var.root=(0.22)^2/(2*2.49))
shiftsBK <- list(edges=c(1, 12, 14, 28),values=c(2.75-0.86, 3.24-0.86, 3.56-3.24, 3.56-3.24),relativeTimes=c(0, 0, 0, 0))
processBK <- "OU"
varianceBK <- (0.22)^2
selection.strengthBK <- 2.49
optimal.valueBK <- 0.86
# Test root.state
root.stateBK <- test.root.state(process=processBK, root.state=root.stateBK, variance=varianceBK, selection.strength=selection.strengthBK, optimal.value=optimal.valueBK)
# Display parameters
paramsBK <- list(root.state=root.stateBK, shifts=shiftsBK, variance=varianceBK, selection.strength=selection.strengthBK, optimal.value=optimal.valueBK)

## Save table history
history.default <- list_to_table.history(results_estim_EM.default$params_history)
write.table.history(history.default, params_algo_EM.default, PATH)
history.lasso <- list_to_table.history(results_estim_EM.lasso$params_history)
write.table.history(history.lasso, params_algo_EM.lasso, PATH)

## Save plots
plot.history.OU.stationnary(results_estim_EM.default$params_history, paramsBK, PATH=PATH, params_algo_EM.default)
plot.history.OU.stationnary(results_estim_EM.lasso$params_history, paramsBK, PATH=PATH, params_algo_EM.lasso)

#######################################
## Test of Divergence - 1 - random
#######################################
rm(list=ls())
#WD <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
WD <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
setwd(WD)
PATH <- paste(WD, "/Results/Miscellaneous_Evals/", sep="")
library(ape)
library(glmnet) # For Lasso initialization
#library(nlme) # For second derivative computation
#library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
library(TreeSim)
source("R/functions.R")

set.seed(25865)

sigma_delta_base <- 18
beta_0 <- 0
process <- "OU"

### Tree

ntaxa <- 64
lambda <- 0.1
tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
plot(tree)

### Functions

datasetsim <- function(alpha, gamma, K, n, grp) {
  shifts <- sample_shifts(tree, sigma_delta_base, K)
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = beta_0,
                     var.root = gamma)
  XX <- simulate(phylo = tree,
                 process = process,
                 root.state = root.state, 
                 variance = 2*alpha*gamma,
                 shifts = shifts, 
                 selection.strength = alpha, 
                 optimal.value = beta_0)
  return(list(alpha = alpha,
              gamma = gamma, 
              K = K,
              n = n,
              grp = grp,
              shifts = shifts,
              Y_data = extract.simulate(XX, what="states", where="tips"),
              Z_data = extract.simulate(XX, what = "states", where = "nodes"),
              params = list(variance = 2*alpha*gamma,
                            root.state = root.state,
                            shifts = shifts,
                            selection.strength = alpha, 
                            optimal.value = beta_0)))
}

sample_shifts <- function(tree, sigma_delta, K){
  shifts_edges <- sample_shifts_edges(tree, K)
  shifts_values <- sample_shifts_values(sigma_delta, K)
  shifts <- list(edges = shifts_edges, 
                 values = shifts_values, 
                 relativeTimes = rep(0, K))
  return(shifts)
}

sample_shifts_edges <- function(tree, K){
  Nbranch <- nrow(tree$edge)
  ntaxa <- length(tree$tip.label)
  if (K > ntaxa) stop("A parcimonious repartition of the shifts cannot have more shifts than tips !")
  ## Generate shifts so that they are parcimonious
  Tr <- incidence.matrix(tree)
  Nbr_grps <- 0
  It <- 0
  while ((Nbr_grps < K+1 && It < 500)) {
    ## Generate K branches
    edges <- sample_edges_intervals(tree, K)
    ## Generate Groups
    groups <- rep(0, ntaxa); names(groups) <- tree$tip.label
    for (i in order(edges)) { # Do it in order
      ed <- edges[i]
      groups[tree$tip.label[Tr[,ed]]] <- i
    }
    Nbr_grps <- length(unique(groups))
    It <- It + 1
  }
  if (It==500) stop("I could not find a parcimonious repartition of the shifts after 500 iterations. Please consider taking a smaller number of shifts.")
  return(edges)
}

sample_edges_intervals <- function(tree, K){
  pas <- 1/K * 0:K
  node_heights <- node.depth.edgelength(tree)
  groups <- split(1:length(node_heights), findInterval(node_heights, pas))
  sh <- NULL; p <- 1
  for (k in 1:K) {
    grp <- groups[[paste(k)]]
    if (is.null(grp)){
      p <- p + 1
    } else {
      if (p <= length(grp)){
        if (length(grp) == 1) {
          sh <- c(sh, grp)
        } else {
          sh <- c(sh, sample(grp, p))
        }
        p <- 1
      } else {
        sh <- c(sh, grp)
        p <- p - length(grp) + 1
      }
    }
  }
  sh <- sapply(sh, function(x) which(tree$edge[,1]==x)[1])
  if (length(sh) < K) {
    p <- K-length(sh)
    sh <- c(sh, sample(tree$edge[-sh], p))
  }
  return(sh)
}

sample_shifts_values <- function(sigma_delta, K){
  return(rnorm(K, mean=0, sd=sqrt(sigma_delta)))
}

estimationfunction <- function(X) {
  ## If an estimation fails, catch error with "try" and try again
  results_estim_EM <- estimateEM(phylo=tree, 
                                 Y_data=X$Y_data, 
                                 tol = list(variance=10^(-4), 
                                            value.root=10^(-4), 
                                            exp.root=10^(-4), 
                                            var.root=10^(-4), 
                                            selection.strength=10^(-3)), 
                                 process = process, 
                                 method.variance = "simple", 
                                 method.init = "lasso",
                                 method.init.alpha = "estimation",
                                 Nbr_It_Max = 1000, 
                                 nbr_of_shifts = X$K, 
                                 alpha_known = FALSE,
                                 min_params=list(variance = 10^(-3), 
                                                 value.root = -10^(3), 
                                                 exp.root = -10^(3), 
                                                 var.root = 10^(-3),
                                                 selection.strength = 10^(-3)),
                                 max_params=list(variance = 10^(3), 
                                                 value.root = 10^(3), 
                                                 exp.root = 10^(3), 
                                                 var.root = 10^(3),
                                                 selection.strength = 10^(3)))
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    z <- matrix(z, nrow = X$K)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X$alpha_estim = params$selection.strength
  X$gamma_estim = params$root.state$var.root
  X$shifts_estim = params$shifts
  X$EM_steps <- attr(results_estim_EM, "Nbr_It")
  X$DV_estim <- attr(results_estim_EM, "Divergence")
  X$CV_estim <- (attr(results_estim_EM, "Nbr_It") != 500) && !X$DV_estim
  X$Zhat <- results_estim_EM$ReconstructedNodesStates
  compute.quality <- function(i) {
    res <- 0 + (selected.edges == i)
    return(mean(colSums(res)))
  }
  edge.quality <- unlist(lapply(params$shifts$edges, compute.quality))
  names(edge.quality) <- params$shifts$edges
  X$edge.quality <- edge.quality
  X$history <- results_estim_EM$params_history
  X$beta_0_estim <- params$root.state$exp.root
  return(X)
}

estimationfunction_alpha_known <- function(X, alphaKN) {
  ## If an estimation fails, catch error with "try" and try again
  results_estim_EM <- estimateEM(phylo=tree, 
                                 Y_data=X$Y_data, 
                                 tol = list(variance=10^(-4), 
                                            value.root=10^(-4), 
                                            exp.root=10^(-4), 
                                            var.root=10^(-4), 
                                            selection.strength=10^(-3)), 
                                 process = process, 
                                 method.variance = "simple", 
                                 method.init = "lasso",
                                 method.init.alpha = "estimation",
                                 Nbr_It_Max = 1000, 
                                 nbr_of_shifts = X$K, 
                                 alpha_known = TRUE, ##
                                 known.selection.strength = alphaKN,
                                 min_params=list(variance = 10^(-3), 
                                                 value.root = -10^(3), 
                                                 exp.root = -10^(3), 
                                                 var.root = 10^(-3),
                                                 selection.strength = 10^(-3)),
                                 max_params=list(variance = 10^(3), 
                                                 value.root = 10^(3), 
                                                 exp.root = 10^(3), 
                                                 var.root = 10^(3),
                                                 selection.strength = 10^(3)))
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    z <- matrix(z, nrow = X$K)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X$alpha_estim = params$selection.strength
  X$gamma_estim = params$root.state$var.root
  X$shifts_estim = params$shifts
  X$EM_steps <- attr(results_estim_EM, "Nbr_It")
  X$DV_estim <- attr(results_estim_EM, "Divergence")
  X$CV_estim <- (attr(results_estim_EM, "Nbr_It") != 500) && !X$DV_estim
  X$Zhat <- results_estim_EM$ReconstructedNodesStates
  compute.quality <- function(i) {
    res <- 0 + (selected.edges == i)
    return(mean(colSums(res)))
  }
  edge.quality <- unlist(lapply(params$shifts$edges, compute.quality))
  names(edge.quality) <- params$shifts$edges
  X$edge.quality <- edge.quality
  X$history <- results_estim_EM$params_history
  X$beta_0_estim <- params$root.state$exp.root
  return(X)
}

## Simulation - Estimation

set.seed(145)
datasim <- datasetsim(alpha = 3, gamma = 0.2, K = 5, 1, 1)
simest <- estimationfunction(datasim)
simest$history[["true"]] <- datasim$params
history <- list_to_table.history(simest$history)
write.csv2(history, paste0(PATH, "boite_noire_7.csv"))

simest_true_alpha <- estimationfunction_alpha_known(datasim, alphaKN = 3)
simest_true_alpha$history[["true"]] <- datasim$params
history_true_alpha <- list_to_table.history(simest_true_alpha$history)
write.csv2(history_true_alpha, paste0(PATH, "boite_noire_true_alpha_3.csv"))

plot(tree); edgelabels(text = round(datasim$shifts$values, 2), edge = datasim$shifts$edges)
x11();
plot(tree); edgelabels(); edgelabels(edge = datasim$shifts$edges, col="red")

#######################################
## Test of Divergence - 2 - simple
#######################################
rm(list=ls())
WD_mac <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
WD_unb <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
WD <- ifelse(file.exists(WD_mac), WD_mac, WD_unb)
# setwd(WD)
PATH <- paste(WD, "/Results/Miscellaneous_Evals/", sep="")
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

### Functions

datasetsim <- function(alpha, gamma, K, shifts, n, grp) {
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = beta_0,
                     var.root = gamma)
  XX <- simulate(phylo = tree,
                 process = process,
                 root.state = root.state, 
                 variance = 2*alpha*gamma,
                 shifts = shifts, 
                 selection.strength = alpha, 
                 optimal.value = beta_0)
  return(list(alpha = alpha,
              gamma = gamma, 
              K = K,
              n = n,
              grp = grp,
              shifts = shifts,
              Y_data = extract.simulate(XX, what="states", where="tips"),
              Z_data = extract.simulate(XX, what = "states", where = "nodes"),
              params = list(variance = 2*alpha*gamma,
                            root.state = root.state,
                            shifts = shifts,
                            selection.strength = alpha, 
                            optimal.value = beta_0)))
}

estimationfunction <- function(X, seg) {
  ## If an estimation fails, catch error with "try" and try again
  results_estim_EM <- estimateEM(phylo=tree, 
                                 Y_data=X$Y_data, 
                                 tol = list(variance=10^(-4), 
                                            value.root=10^(-4), 
                                            exp.root=10^(-4), 
                                            var.root=10^(-4), 
                                            selection.strength=10^(-3)), 
                                 process = process, 
                                 method.variance = "simple", 
                                 method.init = "lasso",
                                 method.init.alpha = "estimation",
                                 Nbr_It_Max = 1000, 
                                 nbr_of_shifts = X$K, 
                                 alpha_known = FALSE,
                                 min_params=list(variance = 10^(-3), 
                                                 value.root = -10^(3), 
                                                 exp.root = -10^(3), 
                                                 var.root = 10^(-3),
                                                 selection.strength = 10^(-3)),
                                 max_params=list(variance = 10^(3), 
                                                 value.root = 10^(3), 
                                                 exp.root = 10^(3), 
                                                 var.root = 10^(3),
                                                 selection.strength = 10^(3)),
                                 methods.segmentation = seg)
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    if (is.null(z)){
      return(NULL)
    } else {     
      return(matrix(z, nrow = X$K))
    }
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X$alpha_estim = params$selection.strength
  X$gamma_estim = params$root.state$var.root
  X$shifts_estim = params$shifts
  X$EM_steps <- attr(results_estim_EM, "Nbr_It")
  X$DV_estim <- attr(results_estim_EM, "Divergence")
  X$CV_estim <- (attr(results_estim_EM, "Nbr_It") != 500) && !X$DV_estim
  X$Zhat <- results_estim_EM$ReconstructedNodesStates
  compute.quality <- function(i) {
    res <- 0 + (selected.edges == i)
    return(mean(colSums(res)))
  }
  edge.quality <- unlist(lapply(params$shifts$edges, compute.quality))
  names(edge.quality) <- params$shifts$edges
  X$edge.quality <- edge.quality
  X$history <- results_estim_EM$params_history
  X$beta_0_estim <- params$root.state$exp.root
  #X$CLL_history <- results_estim_EM$CLL_history
  X$params_estim <- params
  X$number_equivalent_solutions <- results_estim_EM$number_equivalent_solutions
  return(X)
}

estimationfunction_alpha_known <- function(X, alphaKN, seg) {
  ## If an estimation fails, catch error with "try" and try again
  results_estim_EM <- estimateEM(phylo=tree, 
                                 Y_data=X$Y_data, 
                                 tol = list(variance=10^(-4), 
                                            value.root=10^(-4), 
                                            exp.root=10^(-4), 
                                            var.root=10^(-4), 
                                            selection.strength=10^(-3)), 
                                 process = process, 
                                 method.variance = "simple", 
                                 method.init = "lasso",
                                 method.init.alpha = "estimation",
                                 Nbr_It_Max = 1000, 
                                 nbr_of_shifts = X$K, 
                                 alpha_known = TRUE, ##
                                 known.selection.strength = alphaKN,
                                 min_params=list(variance = 10^(-3), 
                                                 value.root = -10^(3), 
                                                 exp.root = -10^(3), 
                                                 var.root = 10^(-3),
                                                 selection.strength = 10^(-3)),
                                 max_params=list(variance = 10^(3), 
                                                 value.root = 10^(3), 
                                                 exp.root = 10^(3), 
                                                 var.root = 10^(3),
                                                 selection.strength = 10^(3)),
                                 methods.segmentation = seg)
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    if (is.null(z)){
      return(NULL)
    } else {     
      return(matrix(z, nrow = X$K))
    }
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X$alpha_estim = params$selection.strength
  X$gamma_estim = params$root.state$var.root
  X$shifts_estim = params$shifts
  X$EM_steps <- attr(results_estim_EM, "Nbr_It")
  X$DV_estim <- attr(results_estim_EM, "Divergence")
  X$CV_estim <- (attr(results_estim_EM, "Nbr_It") != 500) && !X$DV_estim
  X$Zhat <- results_estim_EM$ReconstructedNodesStates
  compute.quality <- function(i) {
    res <- 0 + (selected.edges == i)
    return(mean(colSums(res)))
  }
  edge.quality <- unlist(lapply(params$shifts$edges, compute.quality))
  names(edge.quality) <- params$shifts$edges
  X$edge.quality <- edge.quality
  X$history <- results_estim_EM$params_history
  X$beta_0_estim <- params$root.state$exp.root
  #X$CLL_history <- results_estim_EM$CLL_history
  X$number_equivalent_solutions <- results_estim_EM$number_equivalent_solutions
  return(X)
}

### Params 

process <- "OU"
beta_0 <- 0
alpha <- 3
gamma <- 0.1
K <- 9
#shifts <- list(edges = NULL, values = NULL, relativeTimes = NULL)
# haut placées
#shifts <- list(edges=c(53, 110), values=c(2, -2), relativeTimes=c(0,0))
# dans les feuilles
#shifts <- list(edges=c(17, 118),values=c(10, -10),relativeTimes=c(0,0))
#shifts <- list(edges=c(17, 118, 23, 85, 53, 110, 56, 96, 7),values=c(0.5,1,1.5,-0.5,-1,-1.5,2,-2,5),relativeTimes=c(0,0,0,0,0,0,0,0,0))
shifts <- list(edges=c(7, 17, 23, 53, 56, 85, 96, 110, 118),values=c(2,2,2,2,-2,-2,-2,-2,5),relativeTimes=c(0,0,0,0,0,0,0,0,0))
#shifts <- list(edges=c(7, 17, 23),values=c(2,2,-2),relativeTimes=c(0,0,0))
#shifts <- NULL
#shifts <- list(edges=c(118, 28, 85, 53, 110, 56, 96, 7, 8, 59, 50, 5),values=c(1,1,2,2,-1,-1,-2,-2,5,-5,5,-5),relativeTimes=c(0,0,0,0,0,0,0,0,0,0,0,0))

#seg <- "max_costs_0"
#seg <- "lasso"
#seg <- "best_single_move"
seg <- c("same_shifts", "lasso")
#seg <- c("lasso", "same_shifts", "best_single_move")

name <- paste0("_", paste0(seg, collapse="_"), "_alpha=", alpha, "_gamma=", gamma, "_K=", K, "_edges=", paste0(shifts$edges, collapse="-"), "_values=", paste0(shifts$values, collapse="-"))

params_algo_EM <- list(process = process)

### Tree
set.seed(25865)
ntaxa <- 64
lambda <- 0.1
tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
plot(tree); edgelabels()

## Simulation - Estimation
set.seed(145)
datasim <- datasetsim(alpha = alpha, gamma = gamma, K = K, shifts = shifts, 1, 1)

system.time(simest <- estimationfunction(datasim, seg = seg))
simest$history[["true"]] <- datasim$params
history <- list_to_table.history(simest$history)
history[,"true"]["log_likelihood"] <-log_likelihood.OU(datasim$Y_data, tree, datasim$params)
#CLL_history <- cbind(simest$CLL_history, c(NA, NA))
#history <- rbind(history, CLL_history)
write.csv2(history, paste0(PATH, "boite_noire_alpha_unknown", name, ".csv"))
plot.history.OU.stationnary(simest$history, tree, datasim$params, datasim$Y_data, PATH=PATH, paste0("history_plot", name))

simest_true_alpha <- estimationfunction_alpha_known(datasim, alphaKN = alpha, seg = seg)
simest_true_alpha$history[["true"]] <- datasim$params
history_true_alpha <- list_to_table.history(simest_true_alpha$history)
history_true_alpha[,"true"]["log_likelihood"] <-log_likelihood.OU(datasim$Y_data, tree, datasim$params)
#CLL_history_true_alpha <- cbind(simest_true_alpha$CLL_history, c(NA, NA))
#history_true_alpha <- rbind(history_true_alpha, CLL_history_true_alpha)
write.csv2(history_true_alpha, paste0(PATH, "boite_noire_alpha_known", name, ".csv"))

x11()
plot.process.actual(Y.state = datasim$Y_data,
                    Z.state = datasim$Z_data,
                    phylo = tree, 
                    paramsEstimate = datasim$params)

plot.process.actual(Y.state = simest$Y_data,
                    Z.state = simest$Zhat,
                    phylo = tree, 
                    paramsEstimate = simest$params_estim)

plot(tree); edgelabels(text = round(datasim$shifts$values, 2), edge = datasim$shifts$edges)
x11();
plot(tree); edgelabels(); edgelabels(edge = datasim$shifts$edges, col="red")

## Equivalent solutions
times_shared <- compute_times_ca(tree)
t_tree <-  min(node.depth.edgelength(tree)[1:ntaxa])
Tr <- incidence.matrix(tree)

# Simulated data
eq_shifts_edges_sim <- equivalent_shifts_edges(tree, datasim$shifts$edges)
eq_shifts_values_sim <- equivalent_shifts_values(tree,
                                                 shifts = datasim$shifts,
                                                 beta_0 = beta_0,
                                                 eq_shifts_edges_sim,
                                                 selection.strength = datasim$alpha,
                                                 t_tree, times_shared, Tr)

plot_equivalent_shifts(tree, eq_shifts_edges_sim, eq_shifts_values_sim, paste0(PATH, "sim_"), name)

# Estimated shifts
eq_shifts_edges_estim <- equivalent_shifts_edges(tree, simest$shifts_estim$edges)
eq_shifts_values_estim <- equivalent_shifts_values(tree,
                                             shifts = simest$shifts_estim,
                                             beta_0 = simest$beta_0_estim,
                                             eq_shifts_edges_estim,
                                             selection.strength = simest$alpha_estim,
                                             t_tree, times_shared, Tr)

plot_equivalent_shifts(tree, eq_shifts_edges_estim, eq_shifts_values_estim, paste0(PATH, "estim_"), name)
