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

##
#' @title Extraction function
#'
#' @description
#' \code{extract} the needed quantities out of an S3 object.
#'
#' @param x an S3 object.
# @param node the node where to extract. Default to the root of the tree.
#' @param ... further arguments to be passed to the specific method.
#' 
#' @return An integer giving the number of equivalent parsimonious solutions.
#' 
#' @seealso \code{\link{extract.parsimonyNumber}},
#' \code{\link{extract.parsimonyCost}},
#' \code{\link{extract.enumerate_parsimony}},
#' \code{\link{extract.partitionsNumber}}
#' 
#' @export
##
extract <- function(x, ...) UseMethod("extract")

###############################################################################
## Function implementing the Sankoff alogrithm
###############################################################################
##
#' @title Minimal number of shifts needed to get a clustering.
#'
#' @description
#' \code{parsimonyCost} is an implementation of the Sankoff algorithm,
#' when the cost of transition between two state is always one. It is used
#' in functions \code{\link{parsimonyNumber}} and \code{\link{enumerate_parsimony}}
#' to count or enumerate all the parsimonious solutions given one clustering of the
#' tips.
#'
# @details
# This functin does a recursion up the tree, using functions 
# \code{init.parsimonyCost} for the initialization at the tips, 
# \code{updateUp} for the actual recursion on the tree,
# and \code{update.parsimonyCost} for the actualisation of the parameters.
#
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param clusters the vector of the clusters of the tips. (Default to all the tips
#' in a single group).
#' 
#' @return An S3 class "\code{parsimonyCost}" containing a 
#' (ntaxa + nNodes) x (nclus) matrix of the total number of shifts needed to
#' get the clustering, if starting from a node in state k. The cost can be
#' extract from any subtree with function \code{\link{extract.parsimonyCost}}.
#' 
#' @seealso \code{\link{extract.parsimonyCost}}, \code{\link{parsimonyNumber}}, 
#' \code{\link{enumerate_parsimony}}, \code{\link{partitionsNumber}},
#' \code{\link{equivalent_shifts}}
#' 
#' @examples
#' tree <- read.tree(text="(((1,1),2),2);")
#' plot(tree); nodelabels()
#' clusters <- c(1, 1, 2, 2)
#' costs <- parsimonyCost(tree, clusters)
#' costs
#' 
#' ## Extract the parsimony cost at the root
#' extract(costs)
#' 
#' ## Extract the cost for the sub-tree below node 7
#' extract(costs, 7)
#' 
#' @export
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
  class(costReconstructions) <- "parsimonyCost"
  return(costReconstructions)
}

##
# @title Display a parsimony cost
#
# @description
# \code{print.parsimonyCost} prints the parsimony cost at the root of the tree,
# using function \code{\link{extract.parsimonyCost}}.
#
# @param x an object of class \code{\link{parsimonyCost}}.
# @param ... unused
# 
# @return NULL
# 
# @seealso \code{\link{parsimonyCost}}, \code{\link{extract.parsimonyCost}}
# 
#' @export
#' @method print parsimonyCost
##
print.parsimonyCost <- function(x, ...){
  cat(paste0("\nParsimony cost: ", extract(x), ".\n\n"))
}

##
#' @title Extraction of the actual number of solutions.
#'
#' @description
#' \code{extract.parsimonyCost} takes an object of class "\code{parsimonyCost}",
#' result of function \code{\link{parsimonyCost}}, and computes the minimum cost
#' at the given node.
#'
#' @param x an object of class "\code{parsimonyCost}", result of function
#' \code{\link{parsimonyCost}}.
#' @param node the root node of the subtree. By default, the root of the tree.
#' @param ... unused
#' 
#' @return An integer giving the minimum cost of the subtree.
#' 
#' @seealso \code{\link{parsimonyCost}}
#' 
#' @export
##
extract.parsimonyCost <- function(x, 
                                  node = attr(x, "ntaxa") + 1, ...){
  return(min(x[node, ]))
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
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param clusters the vector of the clusters of the tips.
#' 
#' @return A (ntaxa + nNodes)x(nclus) matrix, with ntaxa first lines initialized as
#' described.
#' 
#' @keywords internal
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
#' 
#' @keywords internal
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
#' This function does a recursion up the tree.
# using functions 
# \code{init.parsimonyNumber} for the initialization at the tips, 
# \code{updateUp} for the actual recursion on the tree,
# and \code{update.parsimonyNumber} for the actualisation of the parameters.
#' The function \code{\link{extract.parsimonyNumber}} gives the result seeked for
#' any subtree.
#' The matrix of costs of the states (number of shifts) is also required, it is
#' computed by function \code{\link{parsimonyCost}}.
#'
#' @param phylo phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param clusters the vector of the clusters of the tips. Default to all the tips
#' in one single cluster.
#' 
#' @return an object of S3 class "\code{parsimonyNumber}" with:
#' \describe{
#'  \item{nbrReconstructions}{a (ntaxa + nNodes) x (nclus)
#'  matrix of loccaly parsimonious solutions starting from a cluster k at a
#'  given node}
#'  \item{costReconstructions}{an object of class "\code{parsimonyCost}",
#'  result of function \code{\link{parsimonyCost}}.}
#' }
#' 
#' @examples
#' tree <- read.tree(text="(((0,1),2),2);")
#' plot(tree); nodelabels()
#' clusters <- c(0, 1, 2, 2)
#' n_sols <- parsimonyNumber(tree, clusters)
#' n_sols
#' 
#' ## Extract the number of parsimonious solutions at the root
#' extract(n_sols)
#' 
#' ## Extract the cost of the solutions from the root
#' extract(n_sols, what = "cost")
#' extract(parsimonyCost(tree, clusters)) # same, more efficient
#' 
#' ## Extract for the sub-tree below node 7
#' extract(n_sols, 7) # Result: 2 (the ancestral state is either "0" or "1"). 
#' 
#' @seealso \code{\link{extract.parsimonyNumber}}, \code{\link{parsimonyCost}}, 
#' \code{\link{enumerate_parsimony}}, \code{\link{partitionsNumber}},
#' \code{\link{equivalent_shifts}}
#' 
#' @export
#' 
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
  res <- list(nbrReconstructions = nbrReconstructions,
              costReconstructions = costReconstructions)
  class(res) <- "parsimonyNumber"
  return(res)
}

##
# @title Display the number of parsimonious solutions
#
# @description
# \code{print.parsimonyNumber} prints the number of equivalent parsimonious
# allocations of shifts, from the root of the tree, using function
# \code{\link{extract.parsimonyNumber}}.
#
# @param x an object of class \code{\link{parsimonyNumber}}.
# @param ... unused
# 
# @return NULL
# 
# @seealso \code{\link{parsimonyNumber}}, \code{\link{extract.parsimonyNumber}}
# 
#' @export
#' @method print parsimonyNumber
##
print.parsimonyNumber <- function(x, ...){
  cat(paste0("\nNumber of parsimonious solutions: ", extract(x), ".\n\n"))
}

##
#' @title Extraction of the actual number of solutions.
#'
#' @description
#' \code{extract.parsimonyNumber} takes the two matrices computed by 
#' \code{\link{parsimonyNumber}}, and compute the actual number of parsimonious
#' solution for any subtree starting from a given node.
#' 
#' @details
#' The parsimonious solutions are the one with the minimum number of shifts (that
#' are given by matrice costReconstructions). This function sums the number of 
#' solutions (given in matrice nbrReconstructions) that have the minimum number of 
#' shifts.
#'
#' @param x an object of class "\code{parsimonyNumber}", result of function
#' \code{\link{parsimonyNumber}}.
#' @param node the root node of the subtree. By default, the root of the tree.
#' @param what the quantity to retrieve. Either "number" for the number of
#' solutions, or "cost" for the minimal cost of a solution. Default to "number".
#' @param ... unused
#' 
#' @return An integer giving the number of equivalent parsimonious solutions.
#' 
#' @seealso \code{\link{parsimonyNumber}}
#' 
#' @export
##
extract.parsimonyNumber <- function(x, 
                                    node = attr(x$nbrReconstructions, "ntaxa") + 1,
                                    what = c("number", "cost"), ...){
  what <- match.arg(what)
  cost <- x$costReconstructions[node, ]
  if (what == "cost") return(min(cost))
  nbr <- x$nbrReconstructions[node, ]
  return(sum(nbr[which(cost == min(cost))]))
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
#' @return A (ntaxa + nNodes)x(nclus) matrix, with ntaxa first lines initialized
#' as described.
#' 
#' @keywords internal
#' 
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
#' 
#' @keywords internal
#' 
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
#' @return A (ndaughters) x (nclus) binary matrix indicating the admissible
#' states for the daughters node when parent node is in state k.
#' 
#' @keywords internal
#' 
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
#' This function uses function \code{\link[ape]{prop.part}} from package \code{ape}.
#'
#' @param tree phylogenetic tree, class \code{\link[ape]{phylo}}.
#' 
#' @return list of size nEdges, entry i is the vector of tips bellow edge i.
#' 
#' @export
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
#' @title Clustering associated to a shift allocation, assuming no homoplasy.
#'
#' @description
#' \code{clusters_from_shifts} take a vector of shifts edges, and gives the
#' clustering of the tips induced by them, in a "no homoplasy" model (i.e. no
#' convergence is allowed).
#'
#' @details
#' By default, this function uses \code{\link{enumerate_tips_under_edges}} to compute 
#' the list of tips under each edge.
#'
#' @param tree phylogenetic tree
#' @param edges a vector of edges of the tree, where the shifts are
#' @param part.list a list giving the descendant tips of each edge
#' 
#' @return list of size n+m-1, entry i is the vector of tips bellow edge i.
#' 
#' @export
#'
##
clusters_from_shifts <- function (tree, edges,
                                  part.list = enumerate_tips_under_edges(tree)) {
  ntaxa <- length(tree$tip.label)
  part <- rep(0, ntaxa)
  if (length(edges) > 0 ){
    ## Re-order tree
    tree_post <- reorder(tree, order = "postorder")
    edges_post <- correspondanceEdges(edges = edges, 
                                      from = tree, to = tree_post)
    ## Visit higher edges before edges closer to the tips
    ed_order <- order(edges_post, decreasing = TRUE)
    for (i in ed_order) {
      part[part.list[[edges[i]]]] <- i
    }
  }
  return(part)
}

# clusters_from_shifts_primary_opt <- function (tree, shifts, T_tree = incidence.matrix(tree)) {
#   delta <- shifts.list_to_vector(tree, shifts)
#   O_Y <- T_tree %*% delta
# #   clus <- as.factor(O_Y)
# #   levels(O_Y) <- 1:length(levels(O_Y))
#   return(O_Y)
# }

# clusters_from_shifts_expectations <- function (tree, 
#                                                shifts, 
#                                                T_tree = incidence.matrix(tree),
#                                                ac = TRUE, ...) {
#   if (ac){
#     ac_tree <- incidence_matrix_actualization_factors(tree = tree, ...)
#     T_tree <- T_tree * ac_tree
#   }
#   delta <- shifts.list_to_vector(tree, shifts)
#   O_Y <- T_tree %*% delta
#   #   clus <- as.factor(O_Y)
#   #   levels(O_Y) <- 1:length(levels(O_Y))
#   return(O_Y)
# }

##
#' @title Check Parsimony, assuming no homoplasy
#'
#' @description
#' \code{check_parsimony} take a vector of shifts edges, and check whether the
#' number of groups of the tips induced by this allocation is exactly the number of
#' shifts plus one. This is equivalent to parsimony when there is no homplasy (i.e. no 
#' convergent regimes).
#'
#' @details
#' This function computes explicitely the clustering of the tips, using 
#' function \code{\link{clusters_from_shifts}}.
#' By default, this function uses \code{\link{enumerate_tips_under_edges}} to compute 
#' the list of tips under each edge, but a list can be provided (to avoid extra
#' computation, if many tests on the same tree are done).
#'
#' @param tree phylogenetic tree
#' @param edges a vector of edges of the tree, where the shifts are
#' @param ... possibly, a list giving the descendant tips of each edge
#' 
#' @return boolean : TRUE if the allocation is parsimonious.
#' 
#' @export
#'
##
check_parsimony <- function(tree, edges, ...){
  clusters <- clusters_from_shifts(tree, edges, ...)
  nclus <- length(unique(clusters))
  nshifts <- length(edges)
  return((nclus - 1) == nshifts)
}

##
#' @title Check wether an allocation of the shifts is parsimonious, 
#' in the "infinite site model".
#'
#' @description
#' \code{check_parsimony_clusters} takes a vector clusters of the tips, and
#' checks wether the number of groups of the tips induced by this allocation is
#' exactly the number of shifts plus one.
#'
#' @details
#' This function computes explicitely the clustering of the tips, using 
#' function \code{check_parsimony}.
#' By default, this function uses \code{enumerate_tips_under_edges} to compute 
#' the list of tips under each edge, but a list can be provided (if many tests are done).
#'
#' @param tree phylogenetic tree.
#' @param clusters a vector of clusters of the tips of the tree (result of
#' function \code{\link{clusters_from_shifts}}).
#' 
#' @return boolean : TRUE if the allocation is parsimonious.
#' 
#' @keywords internal
#'
##
check_parsimony_clusters <- function(tree, edges, clusters){
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
#' equivalent regimes is given by \code{\link{parsimonyNumber}} (which is faster).
#' 
#' @details
# This functin does a recursion up the tree, using functions 
# \code{init.enumerate_parsimony} for the initialization at the tips, 
# \code{updateUp_list} for the effective recursion on the tree,
# and \code{update.enumerate_parsimony} for the actualisation of the parameters.
#' Function \code{\link{extract.enumerate_parsimony}} furnishes the result in a
#' human readable form (for any subtree).
#' Function \code{\link{plot.enumerate_parsimony}} plots all the solutions found on
#' the tree.
#' 
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param clusters a vector representing the group of each tip. (Default to only one
#' group with all the tips.)
#' 
#' @return
#' an S3 object of class "\code{enumerate_parsimony}", with:
#' \describe{
#' \item{nbrReconstructions}{an object of class "\code{parsimonyCost}", result
#' of function \code{\link{parsimonyCost}}.}
#' \item{allocations}{a list of size nNode + ntaxa. Each entry i of the list
#' represents the solutions for the subtree starting at node i. It is a list with
#' nclus entries, each entry being a matrix. A line of the kth matrix for the
#' ith node is one possible allocation of the shifts, starting with regime k
#' for node i.}
#' \item{phylo}{the entry phylogenetic tree}
#' }
#' 
#' @seealso \code{\link{extract.enumerate_parsimony}},
#' \code{\link{plot.enumerate_parsimony}}, \code{\link{parsimonyCost}},
#' \code{\link{parsimonyNumber}}, \code{\link{partitionsNumber}},
#' \code{\link{equivalent_shifts}}
#'
#' @examples
#' tree <- read.tree(text="(((A,B),C),D);")
#' plot(tree)
#' clusters <- c(0, 1, 2, 2)
#' sols <- enumerate_parsimony(tree, clusters)
#' plot(sols)
#' 
#' ## Extract the parsimonious solutions from the root
#' extract(sols) # each line is a solution, with states of each node
#' 
#' ## Extract the number of solutions from the root
#' extract(sols, what = "number")
#' extract(parsimonyNumber(tree, clusters)) # same result, more efficient
#' 
#' ## Extract the cost of the solutions from the root
#' extract(sols, what = "cost")
#' extract(parsimonyCost(tree, clusters)) # same result, more efficient:
#' 
#' ## Extract for the sub-tree below node 7
#' extract(sols, 7) # NAs: non-existing nodes in the sub-tree
#' 
#' @export
##
enumerate_parsimony <- function(phylo,
                                clusters = rep(1, length(phylo$tip.label))){
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
  allocations <- recursionUp_list(phy, allocations,
                                  update.enumerate_parsimony,
                                  costReconstructions, clus, pos)
  attr(allocations, "ntaxa") <- ntaxa
  res <- list(costReconstructions = costReconstructions,
              allocations = allocations,
              phylo = phylo)
  class(res) <- "enumerate_parsimony"
  return(res)
}

##
# @title Display the number of parsimonious solutions
#
# @description
# \code{print.enumerate_parsimony} prints the number of equivalent parsimonious
# allocations of shifts, from the root of the tree, using function
# \code{\link{extract.enumerate_parsimony}}.
#
# @param x an object of class \code{\link{enumerate_parsimony}}.
# @param ... unused
# 
# @return NULL
# 
# @seealso \code{\link{enumerate_parsimony}},
# \code{\link{extract.enumerate_parsimony}},
# \code{\link{plot.enumerate_parsimony}}
# 
#' @export
#' @method print enumerate_parsimony
##
print.enumerate_parsimony <- function(x, ...){
  cat(paste0("\nNumber of parsimonious solutions: ", extract(x, what = "number"), ".\n\n"))
  cat(paste0("Use function 'plot' to see them all.\n\n"))
}

##
#' @title Extract the result of \code{enumerate_parsimony} at a node.
#'
#' @description
#' \code{extract.enumerate_parsimony} returns a matrix containing all the
#' possible regime allocations for the nodes of a given subtree.
#' 
#' @param x an object of class "\code{enumerate_parsimony}",
#' result of function \code{\link{enumerate_parsimony}}.
#' @param node the node where to retrive the parsimony number. Default to the
#' root of the tree.
#' @param what the quantity to retrieve. Either "solutions" for the full
#' solutions, "number" for the number of solutions, or "cost" for the minimal
#' cost of a solution. Default to "solutions"
#' @param ... unused
#'
#' @return A matrix with ntaxa + nNode columns, and as many rows as the number of
#' possible parsimonious reconstructions.
#' 
#' @seealso \code{\link{enumerate_parsimony}}, \code{\link{plot.enumerate_parsimony}}
#' 
#' @export
##
extract.enumerate_parsimony <- function(x, 
                                        node = attr(x$allocations,"ntaxa") + 1,
                                        what = c("solutions", "number", "cost"),
                                        ...){
  what <- match.arg(what)
  cost <- x$costReconstructions[node, ]
  if (what == "cost") return(min(cost))
  allocations <- x$allocations[[node]]
  res <- do.call(rbind, allocations[which(cost == min(cost))])
  if (what == "number") return(nrow(res))
  return(res)
}

##
#' @title Plot all the equivalent solutions.
#'
#' @description
#' \code{plot.enumerate_parsimony} plots a representation of all the equivalent
#' solutions.
#' 
#' @details
#' This function uses functions \code{\link[ape]{plot.phylo}} and
#' \code{\link[ape]{nodelabels}} for the actual plotting of the trees.
#' 
#' @param x an object of class \code{enumerate_parsimony}, result of
#' function \code{\link{enumerate_parsimony}}
#' @param numbering whether to number the solutions. Default to FALSE.
#' @param nbr_col the number of columns on which to display the plot.
#' Default to 3.
#' @param ... further arguments to be passed to \code{\link[ape]{plot.phylo}} or 
#' \code{\link[ape]{nodelabels}}.
#' 
#' @return A plot of the equivalent shifts allocations.
#' 
#' @seealso \code{\link[ape]{plot.phylo}}, \code{\link{enumerate_parsimony}},
#' \code{\link{plot.equivalent_shifts}}, \code{\link[ape]{nodelabels}}
#' 
#' @export
#' 
##
plot.enumerate_parsimony <- function(x,
                                     numbering = FALSE,
                                     nbr_col = 3, ...){
  phylo <- x$phylo
  ntaxa <- length(phylo$tip.label)
  nbrSol <- extract(x, what = "number")
  solutions <- extract(x)
  nbrLignes <- (nbrSol %/% nbr_col) + 1
  if (nbrSol %% nbr_col == 0) nbrLignes <- nbrLignes - 1
  scr <- split.screen(c(nbrLignes, nbr_col))
  for (sol in 1:nbrSol) {
    screen(scr[sol])
    col <- solutions[sol, ]
    col <- as.factor(col)
    levels(col) <- c("black", rainbow(length(levels(col)) - 1,
                                      start = 0, v = 0.5))
    plot.phylo(phylo, edge.color = col[phylo$edge[, 2]], ...)
    nodelabels(text = solutions[sol, (ntaxa+1):ncol(solutions)], ...)
    tiplabels(text = solutions[sol, 1:ntaxa], ...)
    if(numbering){
      legend("topleft",
             legend = sol,
             cex = 1,
             bty = "n",
             text.font = 4)
    }
  }
  close.screen(all.screens = TRUE)
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
#' 
#' @keywords internal
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
#' 
#' @keywords internal
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
#' 
#' @keywords internal
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
#' 
#' @keywords internal
#' 
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
#' 
#' @keywords internal
##
add_complementary <- function(z){
  z <- na.omit(z)
  if (length(z) == 0) return(NA)
  return(z)
}

###############################################################################
## Compute equivalent shifts given one solution
###############################################################################

# equivalent_shifts.BM <- function(tree, shifts, beta_0, ...){
#   ntaxa <- length(tree$tip.label)
#   clusters <- clusters_from_shifts_primary_opt(tree, shifts, ...)
#   betas_allocs <- extract.enumerate_parsimony(enumerate_parsimony(tree, clusters))
#   eq_shifts <- apply(betas_allocs, 1, compute_shifts_from_betas, phylo = tree)
#   betas_0 <- betas_allocs[, ntaxa + 1]
# }

##
#' @title Find all equivalent shifts allocations and values.
#'
#' @description
#' \code{equivalent_shifts} computes the equivalent shifts positions and their
#' corresponding values, assuming an ultrametric tree.
#' 
#' @details
#' This function is only valid for ultrametric trees, and for models: BM, OU with
#' fixed root or stationary root. It assumes that there are no homoplasies.
#' 
#' @param phylo a phylogenetic tree, of class \code{\link[ape]{phylo}}.
#' @param params an object of class \code{params_process}, result inference by
#' function \code{\link{PhyloEM}}, or constructed throught function
#' \code{\link{params_process}}
#' @param T_tree (optional) matrix of incidence of the tree, result of function 
#' \code{\link{incidence.matrix}}
#' @param part.list (optional) list of partition of the tree, result of function
#' \code{\link{enumerate_tips_under_edges}}.
#' @param times_shared (optional) a matrix, result of function
#' \code{\link{compute_times_ca}}.
#' @param ... further arguments to be passed to \code{\link[ape]{plot.phylo}}.
#'
#' @return object of class \code{equivalent_shifts}, whith entries:
#' \describe{
#' \item{eq_shifts_edges}{matrix of equivalent shifts}
#result of function \code{\link{equivalent_shifts_edges}}}
#' \item{shifts_and_betas}{matrix of corresponding shifts values}
#' \item{phylo}{the entry phylogenetic tree}
#' \item{p}{the dimention}
#' }
#' 
#' @seealso \code{\link{plot.equivalent_shifts}},
#' \code{\link{extract.equivalent_shifts}}, \code{\link{params_BM}}, 
#' \code{\link{params_OU}}, \code{\link{enumerate_parsimony}}
#' 
#' @examples
#' ## Simualte a tree
#' set.seed(17920902)
#' ntaxa = 20
#' phylo <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1,
#'                                   mu = 0, age = 1, mrca = TRUE)[[1]]
#' 
#' ## Define parameters (BM, fixed root)
#' params <- params_BM(p = 4, edges = c(4, 17, 22),
#'                     values = cbind(1:4, -(1:4), rep(1, 4)))
#' ## Find equivalent solutions and plot them
#' eq_shifts <- equivalent_shifts(phylo, params)
#' eq_shifts
#' plot(eq_shifts)
#' ## Extract the values
#' # Shifts values for trait 2, for the three shifts (rows), and three solutions (columns)
#' extract(eq_shifts, trait = 2, what = "shifts_values")
#' # Root values for trait 4, for the tree solutions (columns)
#' extract(eq_shifts, trait = 4, what = "root_values")
#' 
#' ## Define parameters (OU, stationary root)
#' params <- params_OU(p = 4, edges = c(4, 17, 22),
#'                     values = cbind(1:4, -(1:4), rep(1, 4)),
#'                     random = TRUE)
#' ## Find equivalent solutions and plot them
#' eq_shifts <- equivalent_shifts(phylo, params)
#' eq_shifts
#' plot(eq_shifts)
#' ## Extract the values
#' # Shifts values for trait 2, for the three shifts (rows), and three solutions (columns)
#' extract(eq_shifts, trait = 2, what = "shifts_values")
#' # Root values for trait 4, for the three solutions (columns)
#' extract(eq_shifts, trait = 4, what = "root_values")
#' 
#' @export
##
equivalent_shifts <- function(phylo, params,
                              T_tree = incidence.matrix(phylo),
                              part.list = enumerate_tips_under_edges(phylo),
                              times_shared = NULL){
  if (length(params$shifts$edges) == 0){
    stop("There are no shifts in the parameters !")
  }
  p <- nrow(params$shifts$values)
  ntaxa <- length(phylo$tip.label)
  nedges <- dim(phylo$edge)[1]
  ## Process
  process <- "BM"
  if (!is.null(params$selection.strength) && sum(abs(params$selection.strength)) > 0) process <- "OU"
  ## Check feasible
  if (process == "OU" && params$root.state$random && !params$root.state$stationary.root){
    stop("The random OU with non-stationary root is not implemented.")
  }
  ## Equivalent shifts positions
  eq_shifts_edges <- equivalent_shifts_edges(phylo, params$shifts$edges,
                                             part.list)
  ## Regression Matrix
  Delta <- shifts.list_to_matrix(phylo, params$shifts)
  W <- diag(1, p*nedges, p*nedges)
  if (process == "OU"){
    if (is.null(times_shared)) times_shared <- compute_times_ca(phylo)
    W <- compute_actualization_matrix_ultrametric(phylo,
                                                  params$selection.strength,
                                                  times_shared = times_shared)
  }
  T_tree_ac <- kronecker(T_tree, diag(1, p, p)) %*% W
  ## Value of the means
  if (params$root.state$random){
    reste <- params$root.state$exp.root
  } else {
    reste <- params$root.state$value.root
  }
  vect_Y <- T_tree_ac %*% as.vector(Delta) + reste
  shifts_and_betas <- equivalent_shifts_values(phylo, 
                                               eq_shifts_edges,
                                               T_tree_ac, vect_Y, p)
  res <- list(eq_shifts_edges = eq_shifts_edges,
              shifts_and_betas = shifts_and_betas,
              phylo = phylo,
              p = p)
  class(res) <- "equivalent_shifts"
  return(res)
}

##
#' @export
#' @method print equivalent_shifts
##
print.equivalent_shifts <- function(x, ...){
  cat(paste0("\nThere are ", ncol(extract(x)), " equivalent solutions.\n\n"))
  cat(paste0("Use function plot to see them all."))
}

##
#' @title Extract the shifts values for one trait.
#'
#' @description
#' \code{extract.equivalent_shifts} takes an object of class
#' \code{equivalent_shifts}, result of function \code{\link{equivalent_shifts}},
#' and returns the shifts of root values for a given trait.
#' 
#' @param x an object of class \code{equivalent_shifts}, result of
#' function \code{\link{equivalent_shifts}}
#' @param trait the number of the trait to be extracted. Default to 1.
#' @param what one of "shifts_values" or "root_values".
#' @param ... unused.
#'
#' @return A matrix with the values of the shifts (\code{what = "shifts_values"}) or
#' the root (\code{what = "root_values"}) for the trait for each equivalent
#' configuration. Each column is one configuration.
#' 
#' @seealso \code{\link{equivalent_shifts}}, \code{\link{plot.equivalent_shifts}},
#' \code{\link{equivalent_shifts_edges}}
#' 
#' @export
##
extract.equivalent_shifts <- function(x, trait = 1,
                                      what = c("shifts_values", "root_values"),
                                      ...){
  what <- match.arg(what)
  if (what == "shifts_values") return(extract_shifts_values(x, trait))
  if (what == "root_values") return(extract_root_values(x, trait))
  return(NULL)
}

extract_shifts_values <- function(eq_shifts, trait){
  nbrShifts <- dim(eq_shifts$eq_shifts_edges)[1]
  return(eq_shifts$shifts_and_betas[1:nbrShifts * eq_shifts$p + trait, , drop = F])
}

extract_root_values <- function(eq_shifts, trait){
  return(eq_shifts$shifts_and_betas[trait, , drop = F])
}

##
#' @title Plot all the equivalent solutions.
#'
#' @description
#' \code{plot.equivalent_shifts} plots a representation of all the equivalent
#' shifts allocations, with a representation of the shifts and their values,
#' and a coloration of the branches in term of regimes.
#' 
#' @details
#' This function uses function \code{\link[ape]{plot.phylo}} for the actual
#' plotting of the trees.
#' 
#' @param x an object of class \code{equivalent_shifts}, result of
#' function \code{\link{equivalent_shifts}}
#' @param trait (integer) the trait to be plotted, if multivariate. Default to 1.
#' @param numbering wheter to number the solutions. Default to FALSE.
#' @param colors_tips user-provided colors for the tips of the tree. A vector
#' vector with as many colors as there are tips. Will be automatically computed
#' if not provided.
#' @param nbr_col the number of columns on which to display the plot.
#' Default to 3.
#' @param gray_scale if TRUE, a gray scale is used instead of colors. Default to
#' FALSE.
#' @param ... further arguments to be passed to \code{\link[ape]{plot.phylo}}.
#' 
#' @return A plot of the equivalent shifts allocations.
#' 
#' @seealso \code{\link{equivalent_shifts}}, \code{\link[ape]{plot.phylo}}
#' 
#' @export
#' 
##
plot.equivalent_shifts <- function(x,
                                   trait = 1,
                                   numbering = FALSE,
                                   colors_tips = NULL,
                                   nbr_col = 3, 
                                   gray_scale = FALSE, ...){
  phylo <- x$phylo
  ntaxa <- length(phylo$tip.label)
  nbrSol <- dim(x$eq_shifts_edges)[2]
  nbrLignes <- (nbrSol %/% nbr_col) + 1
  if (nbrSol %% nbr_col == 0) nbrLignes <- nbrLignes - 1
  nbrShifts <- dim(x$eq_shifts_edges)[1]
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
  shifts_values <- extract_shifts_values(x, trait)
  root_values <- extract_root_values(x, trait)
  for (sol in 1:nbrSol) {
    ## Shifts and beta_0
    params <- list(optimal.value = root_values[, sol],
                   shifts = list(edges = x$eq_shifts_edges[, sol],
                                 values = shifts_values[, sol],
                                 relativeTimes = rep(0, nbrShifts)))
    ## Regimes
    regimes <- allocate_regimes_from_shifts(phylo,
                                            x$eq_shifts_edges[, sol])
    regimes <- as.factor(regimes)
    cor_col_reg <- cbind(sort(unique(regimes[1:ntaxa])), colors)
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
                        shifts_bg = box_col_shifts,
                        edge.color = as.vector(edges_regimes),
                        root_bg = beta_0_col,
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

# plot_equivalent_shifts <- function(phylo, eq_shifts_edges, eq_shifts_values, 
#                                    PATH, name){
#   pdf(paste(PATH, "equivalent_shifts_solutions", name, ".pdf", sep=""),
#       width = 4*3, height = nbrLignes * 3)
#   plot_equivalent_shifts.actual(phylo, eq_shifts_edges, eq_shifts_values)
#   dev.off()
#   
# }

##
#' @title Find all the equivalent shift edges allocations.
#'
#' @description
#' \code{equivalent_shifts_edges} uses function \code{\link{enumerate_parsimony}}
#' to find all the shifts positions that are equivalent to a given one.
#' 
#' @details
#' This function is uses functions \code{\link{enumerate_parsimony}} for the
#' actual computation of equivalent regimes,
#' \code{\link{clusters_from_shifts}} for the clustering of the tips induced
#' by the original set of shifts given, and
#' \code{\link{allocate_shifts_from_regimes}} to convert back a parametrization
#' in term of regimes to a parametrization in term of shifts.
#' 
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param shifts_edges a vector of shifts positions on the edges of the tree.
#' @param part.list (optional) list of partition of the tree, result of function
#' \code{\link{enumerate_tips_under_edges}}.
#'
#' @return a matrix with as many columns as equivalent allocation, each column
#'  representing a possible parsimonious allocation of shifts on the tree.
#'  
#' @seealso \code{\link{equivalent_shifts}}, \code{\link{enumerate_parsimony}}
#'  
#' @keywords internal
##
equivalent_shifts_edges <- function(phylo,
                                    shifts_edges,
                                    part.list = enumerate_tips_under_edges(phylo)){
  clusters <- clusters_from_shifts(phylo, shifts_edges,
                                   part.list = part.list) + 1
  if(!check_parsimony_clusters(phylo, shifts_edges, clusters))
    stop("The edges entered are not parsimonious.")
  regime_allocs <- extract.enumerate_parsimony(enumerate_parsimony(phylo, clusters))
  eq_shifts_edges <- apply(regime_allocs, 1,
                           allocate_shifts_from_regimes, phylo = phylo)
  if (length(shifts_edges) == 1)
    eq_shifts_edges <- matrix(eq_shifts_edges, nrow = 1) # Deal with dimentions
  return(eq_shifts_edges)
}

##
#' @title Find values given edges. OU stationary case. Ultrametric tree.
#'
#' @description
#' \code{equivalent_shifts_values} computes the values of the shifts given all
#' thepossible allocations computed by function
#' \code{\link{equivalent_shifts_edges}}.
#' 
#' @details
#' This function uses the linear representation of the problem. It fist compute
#' the mean at the tips given by the orgininal shifts positions and values, and
#' then uses function \code{\link{qr.solve}}
# (through function \code{find_actualized_shift_values})
#' to find back the values of the shifts, givent their various positions,
#' and the means at the tips. Function \code{compute_actualization_factors} is
#' used to compute the actualization factor that multipies the shifts values at
#' the tips. Carefull, only work for ULTRAMETRIC trees.
#' 
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param shifts a list of positions and values of original shifts.
#' @param beta_0 value of the original optimal value at the root.
#' @param eq_shifts_edges matrix (optional) result of function
#' \code{\link{equivalent_shifts_edges}}.
#' @param T_tree_ac matrix of incidence of the tree, result of function 
#' \code{\link{incidence.matrix}}, actualized with coeficients computed by
#' function \code{\link{incidence_matrix_actualization_factors}}.
#'
#' @return Named list, with "shifts_values" a matrix of shifts values
#' corresponding to the shifts positions in eq_shifts_edges; and "betas_0" a
#' vector of corresponding root optimal values.
#' 
#' @keywords internal
##
equivalent_shifts_values <- function(phylo, 
                                     eq_shifts_edges,
                                     T_tree_ac, vect_Y, p) {
  ## corresponding values at tips
  # delta <- shifts.list_to_vector(phylo, shifts)
  # m_Y <- T_tree_ac %*% delta + beta_0
  ## find the right coefficients for each combination of edges (! stationary case)
  shifts_and_beta <- apply(eq_shifts_edges, 2,
                           find_shift_values,
                           T_tree_ac = T_tree_ac, vect_Y = vect_Y,
                           p = p, ntaxa = length(phylo$tip.label))
  ## exclude NAs column (when qr.solve failed)
  shifts_and_beta <- t(na.omit(t(shifts_and_beta)))
  ## Create Object
  return(shifts_and_beta)
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
#' @param T_tree_ac matrix of incidence of the tree, result of function 
#' \code{incidence.matrix}.
#' @param m_Y the vector of values of the means at the tips.
#'
#' @return vector, with first entry the values at the root, and other entries the
#' values of the shifts.
#' 
#' @keywords internal
##
find_shift_values <- function(shifts_edges, T_tree_ac, vect_Y, p, ntaxa){
  pos_shifts <- as.vector(sapply(shifts_edges, function(z) p * (z-1) + 1:p))
  mat <- cbind(t(matrix(rep(diag(1, p, p), ntaxa), nrow = p)), T_tree_ac[, pos_shifts]) # stationary case assumption used here
  coefs <- try(qr.solve_exact(mat, vect_Y), silent = TRUE)
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
#' @keywords internal
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
#' @param phylo the phylogenetic tree (un-scaled)
#' 
#' @keywords internal
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