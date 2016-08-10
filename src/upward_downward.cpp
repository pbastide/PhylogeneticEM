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
Moments::Moments(arma::umat const & ed, int p_d){
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
Moments::Moments(arma::mat const & data, arma::umat const & ed){
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

// Default Constructor ------------------------------------------------------//
Model_Node::Model_Node(){
  r.set_size(0);
  q.set_size(0, 0);
  sigma.set_size(0, 0);
}

// Constructor for BM Model _Node -------------------------------------------//
Model_Node::Model_Node(arma::vec const & delta, arma::mat const & Variance,
                       double const & edge_length, arma::uvec no_miss){
  int p_d = Variance.n_rows;
  // R
  r = delta( no_miss );
  // Q
  q.eye(p_d, p_d);
  q = q.rows(no_miss);
  sigma = edge_length * Variance;
  sigma = sigma(no_miss, no_miss);
}

// Constructor for BM Model -------------------------------------------------//
Model::Model(arma::mat const & Delta, arma::mat const & Variance,
             arma::vec const & edge_length, arma::mat data, arma::umat ed){
  int p_d = Delta.n_rows;
  int ntaxa = data.n_cols; // number of taxa
  int nEdges = edge_length.n_rows; // number of edges
  // Allocate Object
  size = nEdges;
  mod = new Model_Node [size];
  // Fill Model
  for (int i = 0; i < size; i++){
    arma::vec delta = Delta.col(i);
    int child = ed(i, 1) - 1;
    arma::uvec no_miss = arma::linspace<arma::uvec>(0, p_d - 1, p_d);
    if (child < ntaxa){
      no_miss = arma::find_finite(data.col(child));
    }
    mod[i] = Model_Node(delta, Variance, edge_length(i), no_miss);
  }
  // int p_d = Delta.n_rows;
  // int ntaxa = Delta.n_cols; // number of taxa
  // int nEdges = edge_length.n_rows; // number of edges
  // r = Delta;
  // q.set_size(p_d, p_d, nEdges);
  // sigma.set_size(p_d, p_d, nEdges);
  // for (arma::uword i = 0; i < nEdges; i++){
  //   q.slice(i).eye();
  //   sigma.slice(i) = edge_length(i) * Variance;
  // }
}

// Destructor (array in C++) ------------------------------------------------//
Model::~Model(){
  if (mod != NULL){
    delete[] mod;
  }
  size = 0;
}

// // Remove the missing values ------------------------------------------------//
// Model Model::removeMissing(arma::mat const & miss, arma::mat const & ed){
//   int nEdges = e.n_rows;
//   int ntaxa = miss.n_cols;
//   int p_d = miss.nrows;
//   
// }


// [[Rcpp::export]]
Rcpp::List create_model(arma::mat const & Delta, arma::mat const & Variance,
                        arma::vec const & edge_length,
                        arma::mat data, arma::umat ed, int i) {
  // Construct object and initialize
  Model mod(Delta, Variance, edge_length, data, ed);
  
  // Return to R in list format
  return mod.exportModel2R(i - 1);
}

// Access to fields Model_Node ----------------------------------------------//
arma::vec Model_Node::R() const {
  return r;
}
arma::mat Model_Node::Q() const {
  return q;
}
arma::mat Model_Node::Sigma() const {
  return sigma;
}

// Access to fields Model ---------------------------------------------------//
Model_Node Model::Mod(int i) const {
  if ((i >= size) || (i < 0)){
    Rcpp::stop("In Mod function, index out of bounds.");
  }
  return mod[i];
}
unsigned int Model::Size() const {
  return size;
}

// Export to R (test) -------------------------------------------------------//
Rcpp::List Model::exportModel2R(int i) const {
  if ((i >= size) || (i < 0)){
    Rcpp::stop("In exportModel, index out of bounds.");
  }
  return Rcpp::List::create(Rcpp::Named("r") = mod[i].R(),
                            Rcpp::Named("q") = mod[i].Q(),
                            Rcpp::Named("sigma") = mod[i].Sigma());
}

// Rcpp::List Model::exportModel2R() const {
//   return Rcpp::List::create(Rcpp::Named("r") = r,
//                             Rcpp::Named("q") = q,
//                             Rcpp::Named("sigma") = sigma);
// }

//---------------------------------------------------------------------------//
// Class Upward_Node --------------------------------------------------------//
//---------------------------------------------------------------------------//

// Default Constructor ------------------------------------------------------//
Upward_Node::Upward_Node(){
  cst.set_size(0);
  condexp.set_size(0);
  condvar.set_size(0, 0);
}
// Constructor from dimension -----------------------------------------------//
Upward_Node::Upward_Node(int p_d){
  cst.set_size(1);
  cst.fill(NA_REAL);
  condexp.set_size(p_d);
  condexp.fill(NA_REAL);
  condvar.set_size(p_d, p_d);
  condvar.fill(NA_REAL);
}

// Access to fields ---------------------------------------------------------//
arma::vec Upward_Node::Cst() const {
  return cst;
}
arma::vec Upward_Node::Condexp() const {
  return condexp;
}
arma::mat Upward_Node::Condvar() const {
  return condvar;
}

// Allocate fields ----------------------------------------------------------//
void Upward_Node::allocate_cst(arma::vec c) {
  cst = c;
}
void Upward_Node::cst_ones() {
  cst.ones();
}
void Upward_Node::allocate_condexp(arma::vec exp) {
  condexp = exp;
}
void Upward_Node::Condvar(arma::mat var) {
  condvar = var;
}
void Upward_Node::condvar_zeros() {
  condvar.zeros();
}



//---------------------------------------------------------------------------//
// Class Upward -------------------------------------------------------------//
//---------------------------------------------------------------------------//
// Constructor from edge matrix and dimension -------------------------------//
Upward::Upward(int siz, int p_d){
  size = siz;
  up = new Upward_Node [size];
  for (int i = 0; i < size; i++){
    up[i] = Upward_Node(p_d);
  }
}

// Destructor (array in C++) ------------------------------------------------//
Upward::~Upward(){
  if (up != NULL){
    delete[] up;
  }
  size = 0;
}

// Constructor for initialization -------------------------------------------//
arma::vec remove_na(arma::vec const & data_col){
  arma::vec res = data_col( arma::find_finite(data_col) );
  return res;
}

int findDimensionsNAs(arma::vec const & data_col){
  arma::uvec non_fin = arma::find_finite(data_col);
  return non_fin.n_rows;
}

Upward::Upward(arma::mat const & data, int nE){
  // Allocate Object
  size = nE + 1;
  up = new Upward_Node [size];
  // Fill data
  int ntaxa = data.n_cols;
  for (int i = 0; i < ntaxa; i++){
    int p_d_tip = findDimensionsNAs(data.col(i));
    up[i] = Upward_Node(p_d_tip);
    // Fill constants with ones
    up[i].cst_ones();
    // Fill exp with data
    up[i].allocate_condexp(remove_na(data.col(i)));
    // Fill variances with zeros
    up[i].condvar_zeros();
  }
  // Rest with NAs
  int p_d = data.n_rows;
  for (int  i = ntaxa; i < size; i++){
    up[i] = Upward_Node(p_d);
  }
}

// Access to fields ---------------------------------------------------------//
Upward_Node Upward::Up(int i) const {
  if ((i >= size) || (i < 0)){
    Rcpp::stop("In Up function, index out of bounds.");
  }
  return up[i];
}
unsigned int Upward::Size() const {
  return size;
}
arma::vec Upward::Likelihood() const {
  return up[size - 1].Cst();
}

// Recursion ----------------------------------------------------------------//
arma::uvec findChildren(int father, arma::umat const & ed){
  arma::uvec edges_children = arma::find(ed.col(0) == father);
  arma::umat rows_children = ed.rows( edges_children );
  return rows_children.col(1);
}

// Upward & actualize_upward_children(arma::uvec const & children){
//   int nChild = childern.n_rows;
//   
//   Upward *res = new Upward(nChild, p_d);
//   
//   
//   return *res;
// }

void Upward::recursion(Model const & mod, arma::umat const & ed) {
  int nEdges = ed.n_rows;
  for (arma::uword i = 0; i < nEdges; i++){ // Loop on the edges (rows of ed)
    int father = ed(i, 0) - 1;
    if (arma::is_finite(up[father].Cst())){
      break; // This node has already been visited.
    }
    // Find children of the node
    arma::uvec children = findChildren(father, ed) - 1;
    // Create an array with actualized quantities at each child
    // Upward up_child = actualize_upward_children(children);
    // Compute the quantity for the node
    // Up[father] = merge_upward(up_child);
  }
}

// Export to R (test) -------------------------------------------------------//
Rcpp::List Upward::exportUpward2R(int i) const {
  if ((i >= size) || (i < 0)){
    Rcpp::stop("In exportUpward, index out of bounds.");
  }
  return Rcpp::List::create(Rcpp::Named("cst") = up[i].Cst(),
                            Rcpp::Named("condexp") = up[i].Condexp(),
                            Rcpp::Named("condvar") = up[i].Condvar());
}

//---------------------------------------------------------------------------//
// Function per se ----------------------------------------------------------//
//---------------------------------------------------------------------------//

// [[Rcpp::export]]
Rcpp::List upward_downward(arma::mat const & data, arma::umat const & ed,
                           arma::mat const & Delta, arma::mat const & Variance,
                           arma::vec const & edge_length,
                           int i) {
  // Construct objects and initialize
  //Moments mom(data, ed);
  int nE = ed.n_rows;
  Upward upw(data, nE);
  // BM
  Model mod(Delta, Variance, edge_length, data, ed);
  // Upward Recursion
  //upw.recursion(mod, ed);
  
  
  // Return to R in list format
  // return upw.Likelihood();
  return upw.exportUpward2R(i-1);
  //return mom.exportMoments2R();
}

/*** R
library(ape)
library(TreeSim)
library(Matrix)
tree <- sim.bd.taxa.age(n = 50, numbsim = 1, lambda = 0.1, mu = 0, 
                        age = 1, mrca = TRUE)[[1]]
tree <- reorder(tree, order = "postorder")
p <- 2
variance <- matrix(-0.9, p, p) + diag(1.9, p, p)
root.state <- list(random = FALSE,
                   value.root = c(0, 0),
                   exp.root = NA,
                   var.root = NA)
shifts = list(edges = c(12),
              values=matrix(2*c(1, 0.5), nrow = p),
              relativeTimes = 0)
paramsSimu <- list(variance = variance,
                   shifts = shifts,
                   root.state = root.state)

set.seed(17920902)
source("~/Dropbox/These/Code/Phylogenetic-EM/R/simulate.R")
source("~/Dropbox/These/Code/Phylogenetic-EM/R/generic_functions.R")
source("~/Dropbox/These/Code/Phylogenetic-EM/R/shifts_manipulations.R")
X1 <- simulate(tree,
               p = p,
               root.state = root.state,
               process = "BM",
               variance = variance,
               shifts = shifts)
traits <- extract.simulate(X1, where = "tips", what = "state")
traits[2, 3] <- traits[1, 10] <- traits[, 6] <- NA

# conditional_law_X <- upward_downward(traits, tree$edge)

Delta <- shifts.list_to_matrix(tree, shifts)
edge_length <- tree$edge.length
create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 1))
create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 3))
create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 6))

upward_downward(traits, tree$edge, Delta, variance, edge_length, 1)
upward_downward(traits, tree$edge, Delta, variance, edge_length, 3)
upward_downward(traits, tree$edge, Delta, variance, edge_length, 6)
*/
