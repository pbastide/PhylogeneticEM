# {Parcimony Number}
# Copyright (C) {2014} {SR, MM, PB}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

###############################################################################
## Here is a function to compute the number of partimonious allocation of shifts
## on the tree, given a clustering of the tips.
## Dependencies : generic_functions.R
###############################################################################
##
# parcimonyNumber (phylo,clusters=rep(1,length(phy$tip.label)))
# PARAMETERS:
# @phylo (tree) imput tree
# @clusters (vector) : vector indicating the clusters of each tip
# RETURNS:
# (matrix) matrix with ncluster columns and Nnodes+ntaxa rows. Each row i contains the numbers of most parcimonious reconstruction of the subtree under the node i, if starting with state j (column)
# DEPENDENCIES:
# init.parcimonyNumber, update.parcimonyNumber (, extrac.parcimonyNumber)
# PURPOSE:
# Find the number of parcimonious repartition of the shift that can lead to a given clustering
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
# 21/05/14 - Gestion of case only one cluster. Extraction of the recursion.
##
parcimonyNumber <- function(phylo,clusters=rep(1,length(phylo$tip.label))){
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Initialization (number of parcimonious reconstructions)
  ## Note: those are NOT globally most parcimonious reconstructions, but locally most
  ## parcimonious reconstructions, i.e. given that node i (row) is in state j (column)
  nbrReconstructions <- init.parcimonyNumber(phy,clusters)
  ## Initialization (cost of parcimonious reconstructions)
  costReconstructions <- init.parcimonyCost(phy,clusters)
  ## Tree recursion
  ## For the parcimony cost
  costReconstructions <- recursionUp(phy, costReconstructions, update.parcimonyCost)
  ## For the number of (locally) most parsimonious allocations
  nbrReconstructions <- recursionUp(phy, nbrReconstructions,
                                    update.parcimonyNumber, costReconstructions)
  ## Return globally most parsimonious allocations from locally most parcimonious ones
  ## Apply local filter (at each node) to keep only states that achieves minimum cost
  state.filter <- as.integer(t(apply(costReconstructions, 1, function(x) { x == min(x)} )))
  nbrReconstructions <- nbrReconstructions * state.filter
  attr(nbrReconstructions, "ntaxa") <- ntaxa
  return(nbrReconstructions)
}

##
# init.parcimonyNumber (phy,clusters)
# PARAMETERS:
# @(phy,clusters) see note above
# RETURNS:
# (matrix) matrix with ncluster columns and Nnodes+ntaxa rows. For each tip i, row i is the vector rep(1, length(unique(clusters)))
# DEPENDENCIES:
# none
# PURPOSE:
# Initialise the matrix for function parcimonyNumber
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
init.parcimonyNumber <- function(phy,clusters){
  ntaxa <- length(phy$tip.label)
  clus <- unique(clusters)
  nclus <- length(clus)
  nbrReconstructions <- matrix(NA,nrow=1 + nrow(phy$edge),ncol=nclus)
  for (i in 1:ntaxa){
    nbrReconstructions[i,] <- as.integer(clusters[i] == clus)
  }
  return(nbrReconstructions)
}

##
# init.parcimonyCost (phy,clusters)
# PARAMETERS:
# @(phy,clusters) see note above
# RETURNS:
# (matrix) matrix with ncluster columns and Nnodes+ntaxa rows. For each tip i, row i is the vector 1 - Ind(i=j). All other rows are set to NAs
# DEPENDENCIES:
# none
# PURPOSE:
# Initialise the parcimony cost matrix for function parcimonyNumber
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
init.parcimonyCost <- function(phy,clusters){
  ntaxa <- length(phy$tip.label)
  clus <- unique(clusters)
  nclus <- length(clus)
  costReconstructions <- matrix(NA,nrow=1 + nrow(phy$edge),ncol=nclus)
  for (i in 1:ntaxa){
    costReconstructions[i, ] <- as.integer(clusters[i] != clus)
  }
  return(costReconstructions)
}


##
# update.parcimonyNumber (daughtersNbr)
# PARAMETERS:
# @daughtersNbr (matrix) matrix with ncluster columns. Each row contains the result of parcimonyNumber for the children of a given node
# @costReconstructions (matrix) matrix with ncluster columns. Each row the results of parcimonyCost for the daughters of a given node
# RETURNS:
# (vector) vector containing the number of parcimonious reconstructions from the current node if starting with cluster j
# DEPENDENCIES:
# none
# PURPOSE:
# Update
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
update.parcimonyNumber <- function(daughters, daughtersParams, cost, ...){
    nclus <- dim(daughtersParams)[2]
    nbrAdm <- rep(0,nclus)
    ## Cost of daughter nodes
    cost <- cost[daughters, , drop = F]
    ## Local function to compute the number of allocations for  subtree(parent) when
    ## parent is in state k
    allocation <- function(k) {
        ## List of potential daughter states (i.e that realize the minimum cost for the tree
        ## parent -> daughter -> subtree(daughter) ) when parent is in state k
        candidate.states <- function(x) {
            daughter.cost <- x + as.numeric(1:nclus != k)
            return(as.integer(daughter.cost == min(daughter.cost)))
        }
        ## binary matrix of candidate states for each daughter
        state.filter <- t(apply(cost, 1, candidate.states))
        ## Compute the number of parsimonious allocations in each
        ## parent -> daughter -> subtree(daughter) tree by summing the number
        ## of parsimonious allocations over all candidate daughter states
        result <- rowSums(daughtersParams * state.filter)
        ## The number of allocations in subtree(parent) is the product of allocations in
        ## all parent -> daughter -> subtree(daughter) tree
        return(prod(result))
    }
    for (k in 1:nclus) {
        nbrAdm[k] <- allocation(k)
    }
    return(nbrAdm)
}

##
# update.parcimonyCost (daughtersNbr)
# PARAMETERS:
# @daughtersNbr (matrix) matrix with ncluster columns. Each row contains the result of parcimonyCost for the children of a given node
# RETURNS:
# (vector) vector containing the cost of parcimonious reconstructions from the current node if starting with cluster j
# DEPENDENCIES:
# none
# PURPOSE:
# Update
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
update.parcimonyCost <- function(daughtersParams, ...){
  nclus <- dim(daughtersParams)[2]
  nbrAdm <- rep(0,nclus)
  ## Minimum cost parent -> daughter -> subtree(daughter), over all possible states of daughter
  cost.subtree <- function(x) {
      y <- (1 - diag(1, length(x))) + x
      return(apply(y, 2, min))
  }
  daughtersCost <- t(apply(daughtersParams, 1, cost.subtree))
  ## Sum costs over all daughters
  nbrAdm <- colSums(daughtersCost)
  return(nbrAdm)
}

##
# extract.parcimonyNumber (nbrReconstructions,node=attr(nbrReconstructions, "ntaxa")+1)
# PARAMETERS:
# @nbrReconstructions (matrix) return of the function parcimonyNumber
# @node (int) : node from which to start (default value : root)
# RETURNS:
# (int) Number of parcimonious reconstructions for the subtree starting at node node.
# DEPENDENCIES:
# parcimonyNumber
# PURPOSE:
# Extract the values wanted from the raw result of function parcimonyNumber
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
extract.parcimonyNumber <- function(nbrReconstructions,node=attr(nbrReconstructions, "ntaxa")+1){
  if (is.vector(nbrReconstructions)){
    return(1)
  } else{
    return(sum(nbrReconstructions[node,]))
  }
}

###############################################################################
## Find the clustering of the tips, given the shifts
###############################################################################

##
#' @title Tips descendants of nodes.
#'
#' @description
#' \code{enumerate_tips_under_edges} gives, for each edge of the tree, the labels
#' of the tips that have this edge as an ancestor.
#'
#' @details
#' This function uses function \code{prop.part} from package \code{ape}.
#'
#' @param tree phylogenetic tree
#' 
#' @return list of size n+m-1, entry i is the vector of tips bellow edge i.
#'
##
enumerate_tips_under_edges <- function (tree) {
  ntaxa <- length(tree$tip.label)
  temp <- prop.part(tree)
  subtree.list <- vector("list", nrow(tree$edge))
  for (i in 1:nrow(tree$edge)) {
    node <- tree$edge[i, 2]
    if (node > ntaxa) {
      subtree.list[[i]] <- temp[[node -ntaxa]]
    } else {
      subtree.list[[i]] <- node
    }
  }
  ## Compact alternative using the special structure of prop.part
  ## subtree.list <- c(lapply(1:ntaxa, function(x) x),
  ##                   temp)
  return(subtree.list)
}

##
#' @title Clusters of the tips corresponding to a list of shifts.
#'
#' @description
#' \code{clusters_from_shifts} take a vector of shifts edges, and gives the
#' clustering of the tips induced by them.
#'
#' @details
#' By default, this function uses \code{enumerate_tips_under_edges} to compute 
#' the list of tips under each edge.
#'
#' @param tree phylogenetic tree
#' @param edges a vector of edges of the tree, where the shifts are
#' @param part.list a list giving the descendant tips of each edge
#' 
#' @return list of size n+m-1, entry i is the vector of tips bellow edge i.
#'
##
clusters_from_shifts <- function (tree, edges, part.list = enumerate_tips_under_edges(tree)) {
  ntaxa <- length(tree$tip.label)
  part <- rep(0, ntaxa)
  if (length(edges) > 0 ){
    edges <- sort(edges)
    for (i in 1:length(edges)) {
      part[part.list[[edges[i]]]] <- i
    }
  }
  return(part)
}

###############################################################################
## Enumerate all the equivalent parsimonious solutions
###############################################################################

##
#' @title Enumerate all the possible regime allocations, given a culstering 
#' of the tips.
#'
#' @description
#' \code{enumerate_parsimony} enumerate all the equivalent allocation of the 
#' regimes in the tree, a clustering of the tips being given. The number of such
#' equivalent regimes is given by \code{parcimonyNumber} (which is faster).
#' 
#' @details
#' This functin does a recursion up the tree, using functions 
#' \code{init.enumerate_parsimony} for the initialization at the tips, 
#' \code{updateUp_list} for the effective recursion on the tree,
#' and \code{update.enumerate_parsimony} for the actualisation of the parameters.
#' The function \code{extract_enumerate_parsimony} furnishes the result in a nice
#' form of a matrix (for any subtree).
#' 
#' @param phylo Input tree.
#' @param clusters a vector representing the group of each tip.
#'
#' @return A list of size nNode + ntaxa. Each entry i of the list represents the
#'  solutions for the subtree starting at node i. It is a list with nclus entries,
#'  each entry being a matrix. A line of the kth matrix for the ith node is one
#'  possible allocation of the shifts, starting with regime k for node i.
##

enumerate_parsimony <- function(phylo, clusters = rep(1,length(phylo$tip.label))){
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Initialization
  allocations <- init.enumerate_parsimony(phy, clusters)
  ## Tree recursion
  allocations <- recursionUp_list(phy, allocations, update.enumerate_parsimony)
  attr(allocations, "ntaxa") <- ntaxa
  return(allocations)
}

##
#' @title Initialization for the enumeration of parsimonious solutions.
#'
#' @description
#' \code{init.enumerate_parsimony} is used in function \code{enumerate_parsimony}, 
#' and initialize the correct data structure.
#' 
#' @details
#' This function returns a list with nNodes + ntaxa entries. Entries correponding to
#' the tips are initialized with a list of nclus matrices. For tip i of group k, all
#' matrices are set to NULL, exept for the kth, set to a vector of size
#' nNodes + ntaxa, with entry i set to k, and all the others to NA.
#' 
#' @param phy Input tree.
#' @param clusters a vector representing the group of each tip.
#'
#' @return A list of size nNode + ntaxa, as described above.
##

init.enumerate_parsimony <- function(phy, clusters){
  ntaxa <- length(phy$tip.label)
  nclus <- length(unique(clusters))
  allocations <- vector("list", 1 + nrow(phy$edge))
  allocations <- lapply(allocations, function(z) return(vector("list", nclus)))
  temp <- rep(NA, 1 + nrow(phy$edge))
  for (i in 1:ntaxa){
    allocations[[i]][[clusters[i]]] <- matrix(temp, nrow = 1)
    allocations[[i]][[clusters[i]]][i] <- clusters[i]
  }
  return(allocations)
}

##
#' @title Actualization of the enumeration.
#'
#' @description
#' \code{update.enumerate_parsimony} is used in function \code{enumerate_parsimony}, 
#' and compute the solution for the parent node, given its children.
#' 
#' @details
#' This function takes a list with L entries corresponding to the children of a node,
#' and compute, for all the regimes, the possible allocations starting with parent
#' node in that regime. It uses functions \code{select.matrices} to select the
#' possible states of the children, and \code{matrix_of_possibles} to find the 
#' possible states.
#' 
#' @param daughters vector of daughters nodes.
#' @param daughtersParams list with length(daughters) entries, each entry being a 
#' list of k matrices representing the possible allocations starting from daughter.
#' @param parent the parent node.
#'
#' @return A list of size nclus, each entry being a matrix representing the possible 
#' allocations starting with node parent in state k.
##
update.enumerate_parsimony <- function(daughters, daughtersParams, parent, ...){
  # Number of nodes and clusters
  nedges <- max(sapply(daughtersParams[[1]], 
                       function(z) max(dim(z)[2], length(dim(z)[2]))))
  nclus <- length(daughtersParams[[1]])
  # Computation of costs : sum_l I[e_k^l = 0], for all groups k.
  costs <- sapply(1:nclus, 
                  function(k) sum(sapply(1:length(daughters), 
                                         function(z) is.null(daughtersParams[[z]][[k]]))))
  # Initialization of the list, with all the entries to NULL
  Krond <- which(costs == min(costs))
  possibles <- vector("list", nclus)
  # Update matrices for adequate regimes
  for (k in Krond){
    # Select the possible regimes for each child
    matlist <- lapply(1:length(daughters), 
                      function(z) select.matrices(daughtersParams[[z]], k))
    # From the list of possible regimes, compute the possible regimes staring with 
    # parent node i regime k.
    possibles[[k]] <- matrix_of_possibles(matlist)
    possibles[[k]][ , parent] <- k
  }
  return(possibles)
}

##
#' @title Select possible regimes.
#'
#' @description
#' \code{select.matrices} is used in function \code{update.enumerate_parsimony}, 
#' and find the adequate matrices to keep.
#' 
#' @details
#' If there is a solution for daughter node starting with group k, then return the 
#' corresponding matrix. Otherwise, return all the other solutions, starting with
#' different groups at daughter node.
#' 
#' @param z a list of nclust matrices.
#' @param k the group being considered.
#'
#' @return Matrix of all possible regimes for the node bellow node daughter, when
#' parent is in state k.
##
select.matrices <- function(z, k){
  if (is.null(z[[k]])) return(do.call(rbind,z))
  return(z[[k]])
}

##
#' @title Compute parent matrix from possibles daughter matrices.
#'
#' @description
#' \code{matrix_of_possibles} is used in function \code{update.enumerate_parsimony}
#' to compute, from the list of possible matrices for the daughters, the matrix for
#' the node (a group for the parent being fixed).
#' 
#' @details
#' This function select all possible combinations of rows from all the daughters, and
#' merge then into one using function \code{merge_complementary_vectors}.
#' 
#' @param matrices a list of matrices with ndaughters entries.
#'
#' @return Matrix of all possible regimes fot the subtree bellow node parent.
##
matrix_of_possibles <- function(matrices){
  # Bind all the matrices together
  XY <- do.call(rbind, matrices)
  # Number of possible solutions for each daughter, plus position in the 
  # binded matrix XY.
  nbrs <- sapply(matrices, function(z) dim(z)[1])
  nbrs <- rbind(nbrs, cumsum(nbrs))
  # All possible combinations, taking one row from each matrix representing
  # a daughter.
  comb <- alply(nbrs, 2, function(z) return(z[2]:(z[2] - z[1] + 1)))
  comb <- expand.grid(comb)
  # For each combination, merge all vectors in one.
  XY <- apply(comb, 1, merge_complementary_vectors, XY)
  return(t(XY))
}

##
#' @title Merge several complementary vectors into one.
#'
#' @description
#' \code{merge_complementary_vectors} take several vectors with complementary entries
#' (i.e, vector of same length, that are such that only one vector has a non NA value
#' for each entry), and merge them into one.
#' 
#' @details
#' The vectors are selected from the entry matrix by comb. At each entry, the vectors
#' are added using function \code{add_complementary}.
#' 
#' @param comb vector giving the rows to be kept.
#' @param mat matrice containing the vectors as rows.
#'
#' @return row vector of the same size as entry matrix.
##
merge_complementary_vectors <- function(comb, mat){
  z <- mat[comb,]
  return(apply(z, 2, add_complementary))
}

##
#' @title Add several entries, when only one is not NA.
#'
#' @description
#' \code{add_complementary} return the only non NA value of a given vector.
#' 
#' @details
#' This function is used in function \code{merge_complementary_vectors}
#' 
#' @param z a vector containing at most only one non-NA value.
#'
#' @return row vector of the same size as entry matrix.
##
add_complementary <- function(z){
  z <- na.omit(z)
  if (length(z) == 0) return(NA)
  return(z)
}

##
#' @title Extract the result of \code{enumerate_parsimony} at a node.
#'
#' @description
#' \code{extract.enumerate_parsimony} returns a matrix contianing all the possible
#' regime allocations for the nodes of a given subtree.
#' 
#' @param allocations a list of list of matrices, result of function 
#' \code{enumerate_parsimony}
#'
#' @return A matrix with ntaxa + nNode columns, and as many rows as the number of
#' possible parsimonious reconstructions.
##
extract.enumerate_parsimony <- function(allocations,
                                        node = attr(allocations, "ntaxa") + 1){
  return(do.call(rbind, allocations[[node]]))
}

###############################################################################
## Compute equivalent shifts given one solution
###############################################################################

##
#' @title Find all the equivalent shift edges allocations.
#'
#' @description
#' \code{equivalent_shifts_edges} uses function \code{enumerate_parsimony} to find all
#' the shifts positions that are equivalent to a given one.
#' 
#' @details
#' This function is uses functions \code{enumerate_parsimony} for the actual
#' computation of equivalent regimes, \code{clusters_from_shifts} for the clustering
#' of the tips induced by the original set of shifts given, and 
#' \code{allocate_shifts_from_regimes} to convert back a parametrization in term of
#' regimes to a parametrization in term of shifts.
#' 
#' @param phylo a phylogenetic tree.
#' @param shifts_edges a vector of shifts positions on the edges of the tree.
#'
#' @return a matrix with as many columns as equivalent allocation, each column
#'  representing a possible parsimonious allocation of shifts on the tree.
##
equivalent_shifts_edges <- function(phylo, shifts_edges){
  clusters <- clusters_from_shifts(phylo, shifts_edges) + 1
  regime_allocs <- extract.enumerate_parsimony(enumerate_parsimony(phylo, clusters))
  eq_shifts_edges <- apply(regime_allocs, 1, allocate_shifts_from_regimes, phylo = phylo)
  return(eq_shifts_edges)
}

##
#' @title Find values given edges. OU stationnary case. Ultrametric tree.
#'
#' @description
#' \code{equivalent_shifts_values} computes the values of the shifts given all the 
#' possible allocations computed by function \code{equivalent_shifts_edges}.
#' 
#' @details
#' This function uses the linear representation of the problem. It fist compute the 
#' mean at the tips given by the orgininal shifts positions and values, and then uses
#' function \code{qr.solve} (through function \code{find_actualized_shift_values}) to 
#' find back the values of the shifts, givent their various positions, and the means
#' at the tips. Function \code{compute_actualization_factors} is used to compute the actualization factor that multipies the shifts values at the tips. Carefull, only 
#' work for ultrametric trees.
#' 
#' @param phylo a phylogenetic tree.
#' @param shifts a list of positions and values of original shifts.
#' @param beta_0 value of the original optimal value at the root.
#' @param eq_shifts_edges matrix, result of function \code{equivalent_shifts_edges}.
#' @param selection.strength the selection strength of the OU process.
#' @param t_tree the depth of the ultrametric tree.
#' @param times_shared matrix of times of shared evolution, result of function
#' \code{compute_times_ca}.
#' @param Tr matrix of incidence of the tree, result of function 
#' \code{incidence.matrix}.
#'
#' @return shifts_values a matrix of shifts values corresponding to the shifts
#' positions in eq_shifts_edges.
#' @return betas_0 a vector of corresponding root optimal values.
##
equivalent_shifts_values <- function(phylo, 
                                     shifts, beta_0, 
                                     eq_shifts_edges = equivalent_shifts_edges(phylo, shifts$edges), 
                                     selection.strength, t_tree, times_shared, Tr) {
  ## Compute actualized shifts values
  shifts_ac <- shifts
  fun <- function(shift_pos){
    parents <- phylo$edge[shift_pos, 1]
    factors <- compute_actualization_factors(selection.strength = selection.strength, 
                                             t_tree = t_tree, 
                                             times_shared = times_shared, 
                                             parents = parents)
    return(factors)
  }
  factors_or <- fun(shifts_ac$edges)
  shifts_ac$values <- shifts_ac$values*factors_or
  ## corresponding values at tips
  delta_ac <- shifts.list_to_vector(phylo, shifts_ac)
  m_Y <- Tr%*%delta_ac + beta_0
  ## find the right coefficients for each combination of edges
  shifts_and_beta <- apply(eq_shifts_edges, 2, find_actualized_shift_values, Tr = Tr, m_Y = m_Y)
  ## Apply the correct factor back
  Factors <- apply(eq_shifts_edges, 2, fun)
  return(list(shifts_values = shifts_and_beta[-1,, drop = FALSE]/Factors,
              betas_0 = shifts_and_beta[1,]))
}

##
#' @title Find values given edges. OU stationnary case. Ultrametric tree.
#'
#' @description
#' \code{find_actualized_shift_values} computes the values of the shifts their
#' positions, and the mean values at the tips.
#' Warning : this function does not check for consistency. Please make sure that the
#' shifts postions and the mean values are compatible.
#' 
#' @details
#' This function uses \code{qr.solve} for rectangular linear system solving.
#' 
#' @param shifts_edges a vector of positions of shifts on the tree.
#' @param Tr matrix of incidence of the tree, result of function 
#' \code{incidence.matrix}.
#' @param m_Y the vector of values of the means at the tips.
#'
#' @return vector, with first entry the values at the root, and other entries the 
#' values of the shifts.
##
find_actualized_shift_values <- function(shifts_edges, Tr, m_Y){
  coefs <- qr.solve(cbind(rep(1, dim(Tr)[1]),Tr[,shifts_edges]), m_Y)
  return(coefs)
}

##
#' @title Plot all the equivalent solutions.
#'
#' @description
#' \code{plot_equivalent_shifts} plots and save a representation of all the equivalent
#' shifts allocations, with a representation of the shifts and their values, and a 
#' coloration of the branches in term of regimes. Black is the ancestral regime.
#' 
#' @details
#' This function uses \code{qr.solve} for rectangular linear system solving.
#' 
#' @param phylo a phylogenetic tree.
#' @param eq_shifts_edges matrix, result of function \code{equivalent_shifts_edges}.
#' @param eq_shifts_values, result of function \code{equivalent_shifts_values}.
#' @param PATH a string linking to the path were the plot is to be stored.
#' @param name a name to be given to the plot.
#' 
##
plot_equivalent_shifts <- function(phylo, eq_shifts_edges, eq_shifts_values, 
                                   PATH, name){
  pdf(paste(PATH, "equivalent_shifts_solutions", name, ".pdf", sep=""),
      width = 4*3, height = nbrLignes * 3)
  plot_equivalent_shifts.actual(phylo, eq_shifts_edges, eq_shifts_values)
  dev.off()
  
}

plot_equivalent_shifts.actual <- function(phylo, eq_shifts_edges, eq_shifts_values){
  nbrSol <- dim(eq_shifts_edges)[2]
  nbrLignes <- (nbrSol %/% 3) + 1
  nbrShifts <- dim(eq_shifts_edges)[1]
  colors <- c("black", palette(rainbow(nbrShifts)))
  scr <- split.screen(c(nbrLignes, 3))
  for (sol in 1:nbrSol) {
    ## Shifts and beta_0
    params <- list(optimal.value = eq_shifts_values$betas_0[sol],
                   shifts = list(edges = eq_shifts_edges[, sol],
                                 values = eq_shifts_values$shifts_values[, sol],
                                 relativeTimes = rep(0, nbrShifts)))
    ## Regimes
    regimes <- allocate_regimes_from_shifts(phylo, eq_shifts_edges[, sol])
    edges_regimes <- regimes[phylo$edge[,2]]
    ## Plot
    screen(scr[sol])
    plot.process.actual(0, 0, phylo, params, bg_shifts = colors[1 + 1:nbrShifts], edge.color = colors[1 + edges_regimes], bg_beta_0 = "white", edge.width = 2)
  }
  close.screen(all.screens = TRUE)
}
