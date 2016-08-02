# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

# include "upward_downward.h"

//---------------------------------------------------------------------------//
// Class Moments ------------------------------------------------------------//
//---------------------------------------------------------------------------//

// Default constructor ------------------------------------------------------//
Moments::Moments(int nE, int p_d) {
  // nEdges = nE;
  // p_dim = p_d;
  edge.set_size(nE, 2);
  edge.fill(NA_REAL);
  exps.set_size(p_d, nE + 1);
  exps.fill(NA_REAL);
  vars.set_size(p_d, p_d, nE + 1);
  vars.fill(NA_REAL);
  covars.set_size(p_d, p_d, nE + 1);
  covars.fill(NA_REAL);
}

// Constructor from edge matrix and dimension -------------------------------//
Moments::Moments(arma::mat const & ed, int p_d){
  int nE = ed.n_rows; // number of edges
  edge = ed;
  exps.set_size(p_d, nE + 1);
  exps.fill(NA_REAL);
  vars.set_size(p_d, p_d, nE + 1);
  vars.fill(NA_REAL);
  covars.set_size(p_d, p_d, nE + 1);
  covars.fill(NA_REAL);
}

// Constructor from data and edge matrix ------------------------------------//
Moments::Moments(arma::mat const & data, arma::mat const & ed){
  // init all NAs
  int p_d = data.n_rows;
  *this = Moments(ed, p_d);
  // Fill exp with data
  int ntaxa = data.n_cols;
  exps.head_cols(ntaxa) = data;
  // Zeros variances / covariance for known values
  for (arma::uword i = 0; i < ntaxa; i++){
    for (arma::uword l = 0; l < p_d; l++){
      if (! Rcpp::NumericVector::is_na(exps(l, i))){
        vars.subcube(0, l, i, p_d - 1, l, i).zeros(); // column var
        vars.subcube(l, 0, i, l, p_d - 1, i).zeros(); // line var
        covars.subcube(0, l, i, p_d - 1, l, i).zeros(); // column covar
      }
    }
  }
}

// Access to fields ---------------------------------------------------------//
arma::mat Moments::Exps() const {
  return exps;
}
arma::cube Moments::Vars() const {
  return vars;
}
arma::cube Moments::Covars() const {
  return covars;
}

Rcpp::List Moments::exportMoments2R() const {
  return Rcpp::List::create(Rcpp::Named("expectation") = exps,
                            Rcpp::Named("variances") = vars,
                            Rcpp::Named("covariances") = covars);
}

//---------------------------------------------------------------------------//
// Class Model --------------------------------------------------------------//
//---------------------------------------------------------------------------//

// Constructor for BM -------------------------------------------------------//
// !! edge_length will have been previously re-numbered in the node orders !!//
Model::Model(arma::mat const & Delta, arma::mat const & Variance,
               arma::vec const & edge_length){
  int p_d = Delta.n_rows;
  int ntaxa = Delta.n_cols; // number of taxa
  int nNodes = edge_length.n_cols; // number of nodes (edge_length(ntaxa+1) = NA)
  r = Delta;
  q.set_size(p_d, p_d, nNodes);
  sigma.set_size(p_d, p_d, nNodes);
  for (arma::uword i = 0; i < nNodes; i++){
    q.slice(i).eye();
    sigma.slice(i) = edge_length(i) * Variance;
  }
  q.slice(ntaxa + 1).fill(NA_REAL); // Root
}

// Access to fields ---------------------------------------------------------//
arma::mat Model::R() const {
  return r;
}
arma::cube Model::Q() const {
  return q;
}
arma::cube Model::Sigma() const {
  return sigma;
}

//---------------------------------------------------------------------------//
// Class Upward -------------------------------------------------------------//
//---------------------------------------------------------------------------//
// Constructor from edge matrix and dimension -------------------------------//
Upward::Upward(arma::mat const & ed, int p_d){
  int nE = ed.n_rows; // number of edges
  edge = ed;
  cst.set_size(nE + 1);
  cst.fill(NA_REAL);
  condexp.set_size(p_d, nE + 1);
  condexp.fill(NA_REAL);
  condvar.set_size(p_d, p_d, nE + 1);
  condvar.fill(NA_REAL);
}

// Constructor for initialization -------------------------------------------//
Upward::Upward(arma::mat const & data, arma::mat const & ed){
  // init all NAs
  int p_d = data.n_rows;
  *this = Upward(ed, p_d);
  // Fill constants with ones
  cst.ones();
  // Fill exp with data
  int ntaxa = data.n_cols;
  condexp.head_cols(ntaxa) = data;
  // Fill variances with zeros
  for (arma::uword i = 0; i < ntaxa; i++){
    condvar.slice(i).zeros();
  }
}

// Access to fields ---------------------------------------------------------//
arma::vec Upward::Cst() const {
  return cst;
}
arma::mat Upward::Condexp() const {
  return condexp;
}
arma::cube Upward::Condvar() const {
  return condvar;
}


//---------------------------------------------------------------------------//
// Function per se ----------------------------------------------------------//
//---------------------------------------------------------------------------//

// [[Rcpp::export]]
Rcpp::List upward_downward(arma::mat const & data, arma::mat const & ed) {
  // Construct object and initialize
  Moments mom(data, ed);
  
  // Do stuff
  
  // Return to R in list format
  return mom.exportMoments2R();
}

/*** R
library(ape)
tree <- sim.bd.taxa.age(n = 50, numbsim = 1, lambda = 0.1, mu = 0, 
                        age = 1, mrca = TRUE)[[1]]
tree <- reorder(tree, order = "cladewise")
p <- 2
variance <- matrix(-0.9, p, p) + diag(1.9, p, p)
root.state <- list(random = FALSE,
                   value.root = c(0, 0),
                   exp.root = NA,
                   var.root = NA)
shifts = list(edges = c(12),
              values=matrix(2*c(1, 0.5), nrow = 2),
              relativeTimes = 0)
paramsSimu <- list(variance = variance,
                   shifts = shifts,
                   root.state = root.state)

set.seed(17920902)
X1 <- simulate(tree,
               p = p,
               root.state = root.state,
               process = "BM",
               variance = variance,
               shifts = shifts)
traits <- extract.simulate(X1, where = "tips", what = "state")
traits[2, 3] <- traits[1, 10] <- traits[, 6] <- NA

conditional_law_X <- upward_downward(traits, tree$edge)
*/
