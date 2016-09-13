# {parsimony Number}
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
## Function implementing the Sankoff alogrithm
###############################################################################
##
#' @title Minimal number of shifts needed to get a clustering.
#'
#' @description
#' \code{parsimonyCost} is an implementation of the Sankoff algorithm, when the cost of
#' transition between two state is always one.
#'
#' @details
#' This functin does a recursion up the tree, using functions 
#' \code{init.parsimonyCost} for the initialization at the tips, 
#' \code{updateUp} for the actual recursion on the tree,
#' and \code{update.parsimonyCost} for the actualisation of the parameters.
#'
#' @param phylo phylogenetic tree.
#' @param clusters the vector of the clusters of the tips.
#' 
#' @return A (ntaxa + nNodes) x (nclus) matrix of the total number of shifts needed 
#' to get the clustering, if starting from a node in state k.
##
parsimonyCost <- function(phylo, 
                          clusters = rep(1, length(phylo$tip.label))){
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Initialization (cost of parcimonious reconstructions)
  costReconstructions <- init.parsimonyCost(phy,clusters)
  ## Tree recursion
  costReconstructions <- recursionUp(phy, costReconstructions, update.parsimonyCost)
  attr(costReconstructions, "ntaxa") <- ntaxa
  return(costReconstructions)
}

##
#' @title Initialization for parsimonyCost.
#'
#' @description
#' \code{init.parsimonyCost} initialize a (ntaxa + nNodes) x (nclus) matrix with
#' NAs everywhere, except for the tips.
#' 
#' @details
#' At a tip i in state k, the line-vector is initialized as follow : 
#' (1 - Ind(k=p)_{1<=p<=nclus})*Inf (where Inf * 0 = 0)
#'
#' @param phylo phylogenetic tree.
#' @param clusters the vector of the clusters of the tips.
#' 
#' @return A (ntaxa + nNodes)x(nclus) matrix, with ntaxa first lines initialized as
#' described.
##
init.parsimonyCost <- function(phy,clusters){
  ntaxa <- length(phy$tip.label)
  clus <- unique(clusters)
  nclus <- length(clus)
  costReconstructions <- matrix(NA, nrow = 1 + nrow(phy$edge), ncol = nclus)
  for (i in 1:ntaxa){
    costReconstructions[i, ] <- as.integer(clusters[i] != clus)
  }
  return(costReconstructions)
}

##
#' @title Actualization for parsimonyCost.
#'
#' @description
#' \code{update.parsimonyCost} compute the line vector of a parent node, given
#' the vectors of its daughters.
#' 
#' @details
#' This function computes the cost of putting the parent in a state k, as the 
#' minimum number of shifts needed to get the given clustering of the trees bellow
#' parental node.
#'
#' @param daughtersParams a (ndaughters) x (nclus) matrix with the line vectors of the cost
#' for the daughters node.
#' 
#' @return A line vector corresponding to the parent node.
##
update.parsimonyCost <- function(daughtersParams, ...){
  # require(plyr)
  nclus <- dim(daughtersParams)[2]
  parCosts <- rep(0, nclus)
  ## Minimum cost parent -> daughter -> subtree(daughter), over all possible states of daughter
  cost.subtree <- function(x) {
    y <- (1 - diag(1, length(x))) + x
    return(apply(y, 2, min))
#    return(unname(aaply(y, 2, min, .drop = FALSE)))
  }
  daughtersCost <- plyr::aaply(daughtersParams, 1, cost.subtree, .drop = FALSE)
#  daughtersCost <- t(unname(daughtersCost))
  ## Sum costs over all daughters
  parCosts <- colSums(daughtersCost)
  return(parCosts)
}

extract.parsimonyCost <- function(costReconstructions, 
                                  node = attr(costReconstructions$nbrReconstructions, "ntaxa") + 1){
  return(min(costReconstructions[node, ]))
}

###############################################################################
## Here is a function to compute the number of partimonious allocation of shifts
## on the tree, given a clustering of the tips.
## Dependencies : generic_functions.R
###############################################################################
##
#' @title Number of equivalent parsimonious allocations.
#'
#' @description
#' \code{parsimonyNumber} aims at finding the number of equivalent allocations of
#' the shifts on the tree, i.e allocations that are parsimonious and compatible
#' with a given clustering of the tips.
#'
#' @details
#' This functin does a recursion up the tree, using functions 
#' \code{init.parsimonyNumber} for the initialization at the tips, 
#' \code{updateUp} for the actual recursion on the tree,
#' and \code{update.parsimonyNumber} for the actualisation of the parameters.
#' The function \code{extract_parsimonyNumber} furnishes the result seeked for any
#' subtree.
#' The matrix of costs of the states (number of shifts) is also required, it is computed
#' by function \code{parsimonyCost}.
#'
#' @param phylo phylogenetic tree.
#' @param clusters the vector of the clusters of the tips.
#' 
#' @return nbrReconstructions a (ntaxa + nNodes) x (nclus) matrix of loccaly parsimonious
#'  solutions starting from a cluster k at a given node.
#' @return costReconstructions a (ntaxa + nNodes) x (nclus) matrix of the total number of
#' shifts needed to get the clustering, if starting from a node in state k (result of function
#' \code{parsimonyCost}.
##
parsimonyNumber <- function(phylo, 
                            clusters = rep(1, length(phylo$tip.label))){
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Computation of costs
  costReconstructions <- parsimonyCost(phylo, clusters)
  ## Initialization (number of parcimonious reconstructions)
  ## Note: those are NOT globally most parcimonious reconstructions, but locally most
  ## parcimonious reconstructions, i.e. given that node i (row) is in state j (column)
  nbrReconstructions <- init.parsimonyNumber(phy,clusters)
  ## Tree recursion for the number of (locally) most parsimonious allocations
  nbrReconstructions <- recursionUp(phy, nbrReconstructions,
                                    update.parsimonyNumber, costReconstructions)
  ## Return reconstructions, with their cost
  attr(nbrReconstructions, "ntaxa") <- ntaxa
  return(list(nbrReconstructions = nbrReconstructions,
              costReconstructions = costReconstructions))
}

##
#' @title Initialization for parsimonyNumber.
#'
#' @description
#' \code{init.parsimonyNumber} initialize a (ntaxa + nNodes)x(nclus) matrix with
#' NAs everywhere, except for the tips.
#' 
#' @details
#' At a tip i in state k, the line-vector is initialized as follow : Ind(k=p)_{1<=p<=nclus}
#'
#' @param phy phylogenetic tree.
#' @param clusters the vector of the clusters of the tips.
#' 
#' @return A (ntaxa + nNodes)x(nclus) matrix, with ntaxa first lines initialized as
#' described.
##
init.parsimonyNumber <- function(phy, clusters){
  ntaxa <- length(phy$tip.label)
  clus <- unique(clusters)
  nclus <- length(clus)
  nbrReconstructions <- matrix(NA, nrow=1 + nrow(phy$edge), ncol = nclus)
  for (i in 1:ntaxa){
    nbrReconstructions[i, ] <- as.integer(clusters[i] == clus)
  }
  return(nbrReconstructions)
}

##
#' @title Actualization for parsimonyNumber.
#'
#' @description
#' \code{update.parsimonyNumber} compute the line vector of a parent node, given
#' the vectors of its daughters.
#' 
#' @details
#' This function uses function \code{compute_state_filter} to find all the admissible
#' states of the daughters, given a starting state for the parent.
#'
#' @param daughters the identifiers of the daughters nodes.
#' @param daughtersParams a ndaughters x (nclus) matrix with the line vectors of the number of
#' solutions for the daughters nodes.
#' @param cost the (ntaxa + nNode) x nclus matrix of costs (computed by \code{parsimonyCost}).
#' 
#' @return A line vector corresponding to the parent node.
##
update.parsimonyNumber <- function(daughters, daughtersParams, cost, ...){
    nclus <- dim(daughtersParams)[2]
    nbrAdm <- rep(0, nclus)
    ## Cost of daughter nodes
    cost <- cost[daughters, , drop = F]
    ## Local function to compute the number of allocations for  subtree(parent) when
    ## parent is in state k
    allocation <- function(k) {
        ## List of potential daughter states (i.e that realize the minimum cost for the tree
        ## parent -> daughter -> subtree(daughter) ) when parent is in state k
        state.filter <- as.integer(compute_state_filter(cost, k))
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
#' @title List of potential daughter states when parent is in state k.
#'
#' @description
#' \code{compute_state_filter} compute the admissible daughters states, i.e. states that 
#' realize the minimum cost for the tree parent -> daughter -> subtree(daughter), when 
#' the parent node is in state k.
#' 
#' @details
#' This function is used in functions \code{parsimonyNumber} and \code{enumerate_parsimony}.
#'
#' @param cost a (ndaughters) x (nclus) matrix of the cost of each state for the 
#' daughters nodes.
#' @param k the parental state considered.
#' 
#' @return A (ndaughters) x (nclus) binary matrix indicating the admissible states for
#' the daughters node when parent node is in state k.
##
compute_state_filter <- function (cost, k) {
  nclus <- dim(cost)[2]
  candidate.states <- function(x) {
    daughter.cost <- x + as.numeric(1:nclus != k)
    return(daughter.cost == min(daughter.cost))
  }
  ## binary matrix of candidate states for each daughter
  state.filter <- t(apply(cost, 1, candidate.states))
}

##
#' @title Extraction of the actual number of solutions.
#'
#' @description
#' \code{extract.parsimonyNumber} takes the two matrices computed by 
#' \code{parsimonyNumber}, and compute the atual number of parsimonious solution for
#' any subtree starting from a given node.
#' 
#' @details
#' The parsimonious solutions are the one with the minimum number of shifts (that are
#' given by matrice costReconstructions). This function sums the number of 
#' solutions (given in matrice nbrReconstructions) that have the minimum number of 
#' shifts.
#'
#' @param Reconstructions two (ntaxa + nNodes)x(nclus) matrices, result of function
#' \code{pasimonyNumber}, with :
#'    -- nbrReconstructions the matrix of number of local parsimonious solutions.
#'    -- costReconstructions the matrix of costs (i.e number of shifts).
#' @param node the root node of the subtree. By default, the root of the tree.
#' 
#' @return An integer giving the number of equivalent parsimonious solutions.
##
extract.parsimonyNumber <- function(Reconstructions, 
                                    node = attr(Reconstructions$nbrReconstructions, "ntaxa") + 1){
  cost <- Reconstructions$costReconstructions[node, ]
  nbr <- Reconstructions$nbrReconstructions[node, ]
  return(sum(nbr[which(cost == min(cost))]))
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
#   for (i in 1:nrow(tree$edge)) {
#     node <- tree$edge[i, 2]
#     if (node > ntaxa) {
#       subtree.list[[i]] <- temp[[node - ntaxa]]
#     } else {
#       subtree.list[[i]] <- node
#     }
#   }
  # Compact alternative using the special structure of prop.part
  # List of tips under nodes, ordered by nodes
  subtree.list <- c(as.list(1:ntaxa),temp)
  # Subset subtree.list by daughter nodes of edges, i.e. sort subtree.list in
  # edge order
  subtree.list <- subtree.list[tree$edge[ , 2]]
  return(subtree.list)
}

##
#' @title Clusters of the tips corresponding to a list of shifts, in an "infinite site model".
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
clusters_from_shifts_ism <- function (tree, edges, part.list = enumerate_tips_under_edges(tree)) {
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

clusters_from_shifts_primary_opt <- function (tree, shifts, T_tree = incidence.matrix(tree)) {
  delta <- shifts.list_to_vector(tree, shifts)
  O_Y <- T_tree %*% delta
#   clus <- as.factor(O_Y)
#   levels(O_Y) <- 1:length(levels(O_Y))
  return(O_Y)
}

clusters_from_shifts_expectations <- function (tree, 
                                               shifts, 
                                               T_tree = incidence.matrix(tree),
                                               ac = TRUE, ...) {
  if (ac){
    ac_tree <- incidence_matrix_actualization_factors(tree = tree, ...)
    T_tree <- T_tree * ac_tree
  }
  delta <- shifts.list_to_vector(tree, shifts)
  O_Y <- T_tree %*% delta
  #   clus <- as.factor(O_Y)
  #   levels(O_Y) <- 1:length(levels(O_Y))
  return(O_Y)
}

##
#' @title Check wether an allocation of the shifts is parsimonious, 
#' in the "infinite site model".
#'
#' @description
#' \code{check_parsimony} take a vector of shifts edges, and check wether the number of 
#' groups of the tips induced by this allocation is exactly the number of shifts plus one.
#'
#' @details
#' This function computes explicitely the clustering of the tips, using 
#' function \code{check_parsimony}.
#' By default, this function uses \code{enumerate_tips_under_edges} to compute 
#' the list of tips under each edge, but a list can be provided (if many tests are done).
#'
#' @param tree phylogenetic tree
#' @param edges a vector of edges of the tree, where the shifts are
#' @param ... possibly, a list giving the descendant tips of each edge
#' 
#' @return boolean : TRUE if the allocation is parsimonious.
#'
##
check_parsimony_ism <- function(tree, edges, ...){
  clusters <- clusters_from_shifts_ism(tree, edges, ...)
  nclus <- length(unique(clusters))
  nshifts <- length(edges)
  return((nclus - 1) == nshifts)
}

# check_parsimony <- function(tree, clusters){
#   npars <- extract.parsimonyCost(parsimonyCost(tree, clusters))
#   nshifts <- length(shifts$edges)
#   if (nshifts < npars){
#     stop("There are less shifts than the parsimonious minimum ! There is a problem.")
#   }
#   return(nshifts == npars)
# }

# check_infinite_site_model <- function(tree, shifts, ...){
#   clusters <- clusters_from_shifts_primary_opt(tree, shifts, ...)
#   if (!check_parsimony(tree, clusters)){
#     warning("This allocation is not parsimonious. Could not check for the infinite site model assumption.")
#     return(NA)
#   } else {
#     nclus <- length(unique(clusters))
#     nshifts <- length(shifts$edges)
#     return(nclus == nshifts + 1)
#   }
# }

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
#' equivalent regimes is given by \code{parsimonyNumber} (which is faster).
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
#' @return nbrReconstructions a (ntaxa + nNodes) x (nclus) matrix of loccaly parsimonious
#' solutions starting from a cluster k at a given node.
#' @return allocations a list of size nNode + ntaxa. Each entry i of the list represents the
#'  solutions for the subtree starting at node i. It is a list with nclus entries,
#'  each entry being a matrix. A line of the kth matrix for the ith node is one
#'  possible allocation of the shifts, starting with regime k for node i.
##
enumerate_parsimony <- function(phylo, clusters = rep(1, length(phylo$tip.label))){
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Re-order clusters if necessary
  clus <- unique(clusters)
  pos <- function(z) which(clus == z)
  ## Computation of costs
  costReconstructions <- parsimonyCost(phylo, clusters)
  ## Initialization
  allocations <- init.enumerate_parsimony(phy, clusters, pos)
  ## Tree recursion
  allocations <- recursionUp_list(phy, allocations, update.enumerate_parsimony, costReconstructions, clus, pos)
  attr(allocations, "ntaxa") <- ntaxa
  return(list(costReconstructions = costReconstructions,
              allocations = allocations))
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

init.enumerate_parsimony <- function(phy, clusters, pos){
  ntaxa <- length(phy$tip.label)
  nclus <- length(unique(clusters))
  allocations <- vector("list", 1 + nrow(phy$edge))
  allocations <- lapply(allocations, function(z) return(vector("list", nclus)))
  temp <- rep(NA, 1 + nrow(phy$edge))
  for (i in 1:ntaxa){
    allocations[[i]][[pos(clusters[i])]] <- matrix(temp, nrow = 1)
    allocations[[i]][[pos(clusters[i])]][i] <- clusters[i]
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
update.enumerate_parsimony <- function(daughters, daughtersParams, parent, cost, clus, pos, ...){
  # Number of nodes and clusters
  nedges <- max(sapply(daughtersParams[[1]], 
                       function(z) max(dim(z)[2], length(dim(z)[2]))))
  nclus <- length(daughtersParams[[1]])
  ## Cost of daughter nodes
  cost <- cost[daughters, , drop = F]
  # Initialization of the list, with all the entries to NULL
  possibles <- vector("list", nclus)
  # Update matrices for adequate regimes
  for (k in clus){
    # List of potential daughter states (i.e that realize the minimum cost for the tree
    # parent -> daughter -> subtree(daughter) ) when parent is in state k
    state.filter <- compute_state_filter(cost, pos(k))
    state.filter <- plyr::alply(state.filter, 1, function(z) z)
    # Select the possible regimes for each child
    matlist <- mapply(function(dpar, sfil) do.call(rbind, dpar[sfil]), daughtersParams, state.filter, SIMPLIFY = FALSE)
    # From the list of possible regimes, compute the possible regimes staring with 
    # parent node i regime k.
    possibles[[pos(k)]] <- matrix_of_possibles(matlist)
    possibles[[pos(k)]][ , parent] <- k
  }
  return(possibles)
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
  comb <- plyr::alply(nbrs, 2, function(z) return(z[2]:(z[2] - z[1] + 1)))
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
# extract.enumerate_parsimony <- function(allocations,
#                                         node = attr(allocations, "ntaxa") + 1){
#   return(do.call(rbind, allocations[[node]]))
# }

extract.enumerate_parsimony <- function(Reconstructions, 
                                    node = attr(Reconstructions$allocations, "ntaxa") + 1){
  cost <- Reconstructions$costReconstructions[node, ]
  allocations <- Reconstructions$allocations[[node]]
  return(do.call(rbind, allocations[which(cost == min(cost))]))
}

###############################################################################
## Compute equivalent shifts given one solution
###############################################################################

equivalent_shifts.BM <- function(tree, shifts, beta_0, ...){
  ntaxa <- length(tree$tip.label)
  clusters <- clusters_from_shifts_primary_opt(tree, shifts, ...)
  betas_allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree, clusters))
  eq_shifts <- apply(betas_allocs, 1, compute_shifts_from_betas, phylo = tree)
  betas_0 <- betas_allocs[, ntaxa + 1]
}

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
equivalent_shifts_edges <- function(phylo, shifts_edges, ...){
  if(!check_parsimony_ism(phylo, shifts_edges, ...)) stop("The edges entered are not parsimonious.")
  clusters <- clusters_from_shifts_ism(phylo, shifts_edges, ...) + 1
  regime_allocs <- extract.enumerate_parsimony(enumerate_parsimony(phylo, clusters))
  eq_shifts_edges <- apply(regime_allocs, 1, allocate_shifts_from_regimes, phylo = phylo)
  if (length(shifts_edges) == 1) eq_shifts_edges <- matrix(eq_shifts_edges, nrow = 1) # Deal with dimentions
  return(eq_shifts_edges)
}

##
#' @title Find values given edges. OU stationary case. Ultrametric tree.
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
#' at the tips. Function \code{compute_actualization_factors} is used to compute the
#' actualization factor that multipies the shifts values at the tips. Carefull, only 
#' work for ultrametric trees.
#' 
#' @param phylo a phylogenetic tree.
#' @param shifts a list of positions and values of original shifts.
#' @param beta_0 value of the original optimal value at the root.
#' @param eq_shifts_edges matrix, result of function \code{equivalent_shifts_edges}.
#' @param T_tree_ac matrix of incidence of the tree, result of function 
#' \code{incidence.matrix}, actualized with coeficients computed by function 
#' \code{incidence_matrix_actualization_factors}.
#'
#' @return shifts_values a matrix of shifts values corresponding to the shifts
#' positions in eq_shifts_edges.
#' @return betas_0 a vector of corresponding root optimal values.
##
equivalent_shifts_values <- function(phylo, 
                                     shifts, beta_0, 
                                     eq_shifts_edges = equivalent_shifts_edges(phylo,
                                                                               shifts$edges), 
                                     T_tree_ac) {
  ## corresponding values at tips
  delta <- shifts.list_to_vector(phylo, shifts)
  m_Y <- T_tree_ac %*% delta + beta_0
  ## find the right coefficients for each combination of edges (! stationary case)
  shifts_and_beta <- apply(eq_shifts_edges, 2,
                           find_shift_values, T_tree_ac = T_tree_ac, m_Y = m_Y)
  ## exclude NAs column (when qr.solve failed)
  shifts_and_beta <- t(na.omit(t(shifts_and_beta)))
  return(list(shifts_values = shifts_and_beta[-1,, drop = FALSE],
              betas_0 = shifts_and_beta[1,]))
}

##
#' @title Find values given edges. OU stationary case. Ultrametric tree.
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
find_shift_values <- function(shifts_edges, T_tree_ac, m_Y){
  mat <- cbind(rep(1, dim(T_tree_ac)[1]), T_tree_ac[, shifts_edges]) # stationary case assumption used here
  coefs <- try(qr.solve_exact(mat, m_Y), silent = TRUE)
  if (inherits(coefs, "try-error")){
    warning("Had a problem solving exactly the linear system.")
    return(rep(NA, length(shifts_edges) + 1))
  } else {
    return(coefs)
  }
}

##
#' @title exact qr.solve
#'
#' @description
#' This is the same function as \code{qr.solve}, but it throws an error if an exact fit cannot
#' be found (instead of returning a least square fitted value, which is the default behavior
#' of \code{qr.solve}). 
#' 
#' @param a a QR decomposition or a rectangular matrix.
#' @param b a vector or matrix of right-hand sides of equations.
#' \code{incidence.matrix}.
#' @param tol the tolerance for detecting linear dependencies in the columns of x.
#' Only used if LAPACK is false and x is real.
#'
##
qr.solve_exact <- function (a, b, tol = 1e-07) {
  if (!is.qr(a)) 
    a <- qr(a, tol = tol)
  nc <- ncol(a$qr)
  nr <- nrow(a$qr)
  if (a$rank != min(nc, nr)) 
    stop("singular matrix 'a' in solve")
  if (missing(b)) {
    if (nc != nr) 
      stop("only square matrices can be inverted")
    b <- diag(1, nc)
  }
  test <- qr.resid(a, b)
  if (!isTRUE(all.equal(as.vector(test), rep(0, nr)))) {
    stop("This linear system cannot be solve exactly")
  }
  res <- qr.coef(a, b)
  res[is.na(res)] <- 0
  return(res)
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

plot_equivalent_shifts.actual <- function(phylo,
                                          eq_shifts_edges,
                                          eq_shifts_values,
                                          numbering = FALSE,
                                          colors_tips = NULL,
                                          nbr_col = 3, 
                                          gray_scale = FALSE, ...){
  ntaxa <- length(phylo$tip.label)
  nbrSol <- dim(eq_shifts_edges)[2]
  nbrLignes <- (nbrSol %/% nbr_col) + 1
  if (nbrSol %% nbr_col == 0) nbrLignes <- nbrLignes - 1
  nbrShifts <- dim(eq_shifts_edges)[1]
  ## Colors
  if (is.null(colors_tips)){
    if (!gray_scale){
      colors <- c("black", rainbow(nbrShifts, start = 0, v = 0.5))
    } else {
      colors <- gray.colors(nbrShifts + 1, start = 0, end = 0.8)
    }
    cor_col_reg <- as.factor(colors)
    levels(cor_col_reg) <- 0:nbrShifts
  } else {
    colors <- unique(colors_tips)
  }
  scr <- split.screen(c(nbrLignes, nbr_col))
  for (sol in 1:nbrSol) {
    ## Shifts and beta_0
    params <- list(optimal.value = eq_shifts_values$betas_0[sol],
                   shifts = list(edges = eq_shifts_edges[, sol],
                                 values = eq_shifts_values$shifts_values[, sol],
                                 relativeTimes = rep(0, nbrShifts)))
    ## Regimes
    regimes <- allocate_regimes_from_shifts(phylo, eq_shifts_edges[, sol])
    regimes <- as.factor(regimes)
    cor_col_reg <- cbind(unique(regimes[1:ntaxa]), colors)
    levels(regimes)[as.numeric(cor_col_reg[,1])] <- colors
    edges_regimes <- regimes[phylo$edge[,2]]
    ## Shifts Colors
    makeLighter = function(..., alpha = 0.5, saut=100) {
      alpha = floor(255*alpha)  
      newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
      .makeTransparent = function(col, alpha) {
        rgb(red=col[1] + saut, green=col[2] + saut, blue=col[3] + saut, maxColorValue=255)
      }
      newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
      return(newColor)
    }
    if (!gray_scale){
      box_col <- as.vector(edges_regimes)
      box_col <- makeLighter(box_col)
      box_col_shifts <- box_col[params$shifts$edges]
      beta_0_col <- box_col[which(!(box_col %in% box_col_shifts))[1]]
    } else {
      box_col_shifts <- rep("white", length(params$shifts$edges))
      beta_0_col <- "white"
    }
    ## Plot
    screen(scr[sol])
    plot.process.actual(0, 0, phylo, params,
                        bg_shifts = box_col_shifts,
                        edge.color = as.vector(edges_regimes),
                        bg_beta_0 = beta_0_col,
                        edge.width = 2, quant.root = 0.7, ...)
    if(numbering){
      legend("topleft",
             legend = sol,
             cex = 1,
             bty = "n",
             text.font = 4)
#              x.intersp = 0,
#              y.intersp = 0)
    }
  }
  close.screen(all.screens = TRUE)
}


##
#' @title Transform the shift values
#'
#' @description
#' \code{transform_shifts_values} takes the shifts generating a given expectation structure
#' given an OU with alpha = from, and gives back the equivelent shifts values that produce the
#' same structure with an OU with alpha = to. If from or to is 0, then the process is supposed
#' to be a BM.
#' 
#' 
#' @param shifts the shifts on the original process
#' @param from alpha value of the original process. If equals 0, then the original process is 
#' taken to be a BM.
#' @param to alpha value of the destination process
#' @param phylo the phlogenetic tree (un-scaled)
#' 
##

transform_shifts_values <- function(shifts, from = 0, to, phylo){
  if (!is.ultrametric(phylo)) stop("The processes are not equivelent on a non-ultrametric tree.")
  depths <- ape::node.depth.edgelength(phylo)
  h_tree <- depths[1]
  parents <- phylo$edge[shifts$edges, 1]
  times_parents <- depths[parents]
  if (from == 0 && to == 0){
    actualizations <- rep(1, length(times_parents))
  } else if (from == 0){
    actualizations <- 1 / (1 - exp(- to * (h_tree - times_parents)))
  } else if (to == 0){
    actualizations <- (1 - exp(- from * (h_tree - times_parents)))
  } else {
    actualizations <- (1 - exp(- from * (h_tree - times_parents))) / (1 - exp(- to * (h_tree - times_parents)))
  }
  if (is.vector(shifts$values)){
    shifts$values <- shifts$values * actualizations
  } else {
    shifts$values <- sweep(shifts$values, 2, actualizations, '*')
  }
  return(shifts)
}