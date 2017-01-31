# include <RcppArmadillo.h>
// # include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

# include "upward_downward.h"

//---------------------------------------------------------------------------//
// Class Model --------------------------------------------------------------//
//---------------------------------------------------------------------------//

// Default Constructor ------------------------------------------------------//
Model::Model(int siz){
  rs.set_size(0, siz);
  qs.set_size(0, 0, siz);
  sigmas.set_size(0, 0, siz);
}

// Constructor for BM Model -------------------------------------------------//
Model::Model(arma::mat const & Delta, arma::mat const & Variance,
             arma::vec const & edge_length){
  int p_d = Variance.n_rows;
  int nEdges = edge_length.n_rows; // number of edges
  // R
  rs = Delta;
  // Q and Sigma
  qs.set_size(p_d, p_d, nEdges);
  sigmas.set_size(p_d, p_d, nEdges);
  for (int i = 0; i < nEdges; i++){
    qs.slice(i).eye();
    sigmas.slice(i) = edge_length(i) * Variance;
  }
}

// Constructor for OU Model -------------------------------------------------//
Model::Model(arma::mat const & Beta, arma::mat const & Stationary_Var,
             arma::vec const & edge_length, arma::mat const & Alpha){
  int p_d = Beta.n_rows;
  int nEdges = edge_length.n_rows; // number of edges
  arma::mat Id;
  Id.eye(p_d, p_d);
  // Fill Model
  rs.set_size(p_d, nEdges);
  qs.set_size(p_d, p_d, nEdges);
  sigmas.set_size(p_d, p_d, nEdges);
  for (int i = 0; i < nEdges; i++){
    qs.slice(i) = arma::expmat( - edge_length(i) * Alpha);
    rs.col(i) = (Id - qs.slice(i)) *  Beta.col(i);
    sigmas.slice(i) = Stationary_Var - qs.slice(i) * Stationary_Var * qs.slice(i).t();
  }
}

// // [[Rcpp::export]]
// Rcpp::List create_model(arma::mat const & Delta, arma::mat const & Variance,
//                         arma::vec const & edge_length,
//                         arma::umat ed) {
//   // Construct object and initialize
//   Model mod(Delta, Variance, edge_length);
//   
//   // Return to R in list format
//   return mod.exportModel2R();
// }

// Access to fields Model ---------------------------------------------------//
arma::vec Model::Rs(int edge) const {
  return rs.col(edge);
}
arma::mat Model::Qs(int edge) const {
  return qs.slice(edge);
}
arma::mat Model::Sigmas(int edge) const {
  return sigmas.slice(edge);
}

// // Export to R (test) -------------------------------------------------------//
// Rcpp::List Model::exportModel2R() const {
//   return Rcpp::List::create(Rcpp::Named("r") = rs,
//                             Rcpp::Named("q") = qs,
//                             Rcpp::Named("sigma") = sigmas);
// }

//---------------------------------------------------------------------------//
// Class Upward -------------------------------------------------------------//
//---------------------------------------------------------------------------//

// Constructor from edge matrix and dimension -------------------------------//
Upward::Upward(int siz, int p_d){
  csts.set_size(siz);
  csts.fill(NA_REAL);
  condexps.set_size(p_d, siz);
  condexps.fill(NA_REAL);
  condvars.set_size(p_d, p_d, siz);
  condvars.fill(NA_REAL);
  missing_datas.set_size(p_d, siz);
  missing_datas.zeros(); // none missing
}

// Constructor for initialization -------------------------------------------//
Upward::Upward(arma::mat const & data, int nE){
  // Allocate Object right size with NAs
  int size = nE + 1;
  int p_d = data.n_rows;
  *this = Upward(size, p_d);
  // Fill data
  int ntaxa = data.n_cols;
  for (int i = 0; i < ntaxa; i++){
    // Nas positions
    arma::vec data_col(data.col(i));
    arma::uvec na_position = arma::find_nonfinite(data_col);
    arma::uvec miss_data(p_d);
    miss_data.zeros();
    miss_data(na_position).ones(); // put one when missing
    missing_datas.col(i) = miss_data;
    // Fill exp with data
    data_col(na_position).zeros(); // replace na with zero
    condexps.col(i) = data_col;
    // Fill variances with zeros
    arma::uvec non_na_position = arma::find_finite(data_col);
    arma::mat var(p_d, p_d);
    var.fill(arma::datum::inf);
    var(non_na_position, non_na_position).zeros();
    condvars.slice(i) = var;
    // Fill constants with right coef
    int nMiss = arma::sum(miss_data);
    csts(i) = nMiss * std::log(2 * arma::datum::pi) / 2;
  }
}

// Access to fields
double Upward::Csts(int node) const{
  return csts(node);
}

arma::vec Upward::Condexps(int node) const{
  return condexps.col(node);
}

arma::mat Upward::Condvars(int node) const{
  return condvars.slice(node);
}

arma::uvec Upward::Missing_Datas(int node) const{
  return missing_datas.col(node);
}

int Upward::Size() const{
  return condexps.n_cols;
}

// Allocate Fields
void Upward::allocate_condvars(int i, arma::mat M){
  condvars.slice(i) = M;
}

void Upward::allocate_condexps(int i, arma::vec E){
  condexps.col(i) = E;
}

void Upward::allocate_csts(int i, double c){
  csts(i) = c;
}

void Upward::allocate_missing_datas(int i, arma::uvec M){
  missing_datas.col(i) = M;
}

// Log Likelihood---------------------------------------------------------//

double Log_Likelihood_Gauss(arma::vec mean, arma::mat var, arma::vec point){
  int p_d = var.n_rows;
  double sign;
  double res;
  arma::log_det(res, sign, var);
  res += p_d * std::log(2 * arma::datum::pi);
  arma::vec diff = point - mean;
  res += as_scalar(diff.t() * arma::inv_sympd(var) * diff);
  return -res/2;
}

double Upward::Log_Likelihood(Root_State root_state, int ntaxa) const {
  double res = csts(ntaxa + 1 - 1);
  arma::vec mean = condexps.col(ntaxa + 1 - 1);
  arma::mat var;
  arma::vec point = root_state.Exp();
  if (root_state.Random()){
    var = condvars.slice(ntaxa + 1 - 1) + root_state.Var();
  } else {
    var = condvars.slice(ntaxa + 1 - 1);
  }
  res += Log_Likelihood_Gauss(mean, var, point);
  return res;
}


// Recursion ----------------------------------------------------------------//
arma::uvec findEdges(int father, arma::umat const & ed){
  arma::uvec edges_children = arma::find(ed.col(0) == father + 1);
  return edges_children;
}

arma::uvec findChildren(int father, arma::umat const & ed){
  arma::uvec edges_children = findEdges(father, ed);
  arma::umat rows_children = ed.rows( edges_children );
  return rows_children.col(1);
}

arma::mat inv_na(arma::mat const & S, arma::uvec const & miss, double ff){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::mat S_sub = S(non_na_pos, non_na_pos);
  S_sub = arma::inv_sympd(S_sub);
  arma::mat res(size(S));
  res.fill(ff);
  res(non_na_pos, non_na_pos) = S_sub;
  return res;
}

// double det_na(arma::mat const & S, arma::uvec const & miss){
//   arma::uvec non_na_pos = find(miss == 0); // position of non missing
//   arma::mat S_sub = S(non_na_pos, non_na_pos);
//   double res = arma::det(S_sub);
//   return res;
// }

double log_det_na(arma::mat const & S, arma::uvec const & miss){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::mat S_sub = S(non_na_pos, non_na_pos);
  double val;
  double sign;
  arma::log_det(val, sign, S_sub);
  return val;
}

arma::vec prod_na(arma::mat const & S, arma::vec const & m, arma::uvec const & miss, double ff){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::vec res = ff * m;
  res(non_na_pos) = S(non_na_pos, non_na_pos) * m(non_na_pos);
  return res;
}

void Upward::recursion(Model const & mod, arma::umat const & ed,
                       int p_d, int ntaxa) {
  int nEdges = ed.n_rows;
  for (int i = 0; i < nEdges; i++){ // Loop on the edges (rows of ed)
    int father = ed(i, 0) - 1;
    if (! arma::is_finite(csts(father))){// This node has not already been visited.
      // Find children of the node
      arma::uvec child_nodes = findChildren(father, ed) - 1;
      int nChild = child_nodes.n_rows;
      // Find Models associated to edges
      arma::uvec child_edges = findEdges(father, ed);
      
      // Initialize accumulated quantities
      arma::mat S_sum;
      S_sum.zeros(p_d, p_d);
      arma::vec S_m_sum;
      S_m_sum.zeros(p_d);
      double cst_sum = 0;
      arma::uvec merge_missing(p_d);
      merge_missing.ones();
      // Compute the quantity for the node
      for (int j = 0; j < nChild; j++){
        int node = child_nodes(j);
        int edge = child_edges(j);
        
        // Compute Condvar, Condexp, constant and missing_data
        arma::uvec missing_data = missing_datas.col(node);
        arma::mat check_S = condvars.slice(node) + mod.Sigmas(edge);
        arma::mat check_S_inv = inv_na(check_S, missing_data, 0);
        arma::mat Q_edge = mod.Qs(edge);
        arma::mat tQ_checkS = Q_edge.t() * check_S_inv;
        arma::mat tQ_checkS_Q = tQ_checkS *  Q_edge;
        arma::vec diff_exp = condexps.col(node) - mod.Rs(edge);
        arma::vec checkSinv_diff = prod_na(check_S_inv, diff_exp, missing_data, 0);
        arma::vec tQ_checkSinv_diff = Q_edge.t() * checkSinv_diff;
        double constant = csts(node);
        constant -= log_det_na(check_S, missing_data) / 2;
        constant -= arma::as_scalar(diff_exp.t() * checkSinv_diff)/2;
      
        // Sums and products
        S_sum += tQ_checkS_Q;
        S_m_sum += tQ_checkSinv_diff;
        cst_sum += constant;
        merge_missing = merge_missing % missing_data; // intersection of missing
      }
      
      // Put everything together
      arma::mat S_bar = inv_na(S_sum, merge_missing, arma::datum::inf);
      arma::vec m_bar = prod_na(S_bar, S_m_sum, merge_missing, 0);
      double constant = - (nChild - 1) * p_d * std::log(2 * arma::datum::pi) / 2;
      constant += log_det_na(S_bar, merge_missing) / 2; // + log_det_sum;
      constant += arma::as_scalar(m_bar.t() * S_sum * m_bar) / 2;
      constant += cst_sum;
      
      // Allocate updated quantities
      condvars.slice(father) = S_bar;
      condexps.col(father) = m_bar;
      csts(father) = constant;
      missing_datas.col(father) = merge_missing;
    }
  }
}

// // Export to R (test) -------------------------------------------------------//
// Rcpp::List Upward::exportUpward2R() const {
//   return Rcpp::List::create(Rcpp::Named("cst") = csts,
//                             Rcpp::Named("condexp") = condexps,
//                             Rcpp::Named("condvar") = condvars,
//                             Rcpp::Named("missing_data") = missing_datas);
// }
// 
// // [[Rcpp::export]]
// Rcpp::List upward_test(arma::mat const & data, arma::umat const & ed) {
//   int nE = ed.n_rows;
//   Upward upw(data, nE);
//   return upw.exportUpward2R();
// }

//---------------------------------------------------------------------------//
// Root_State ---------------------------------------------------------------//
//---------------------------------------------------------------------------//

// Constructors -------------------------------------------------------------//
Root_State::Root_State(arma::vec state){
  random = false;
  exp = state;
  var.set_size(1, 1);
  var.fill(NA_REAL);
}

Root_State::Root_State(arma::vec expe, arma::mat vari){
  random = true;
  exp = expe;
  var = vari;
}

Root_State::Root_State(bool rootRand, arma::vec rootVal,
                       arma::vec rootExp, arma::mat rootVar){
  if (rootRand){
    *this = Root_State(rootExp, rootVar);
  } else {
    *this = Root_State(rootVal);
  }
}

Root_State::Root_State(Rcpp::List root_state){
  int rand_root_int = root_state["random"];
  bool rand_root = rand_root_int;
  if (rand_root){
    arma::vec rootExp = root_state["exp.root"];
    arma::mat rootVar = root_state["var.root"];
    *this = Root_State(rootExp, rootVar);
  } else {
    arma::vec rootVal = root_state["value.root"];
    *this = Root_State(rootVal);
  }
}

// Acces to fields ----------------------------------------------------------//
bool Root_State::Random() const{
  return random;
}
arma::vec Root_State::Exp() const{
  return exp;
}
arma::mat Root_State::Var() const{
  return var;
}

//---------------------------------------------------------------------------//
// Class Moments ------------------------------------------------------------//
//---------------------------------------------------------------------------//

// Default constructor ------------------------------------------------------//
Moments::Moments(int nE, int p_d) {
  exps.set_size(p_d, nE + 1);
  exps.fill(NA_REAL);
  vars.set_size(p_d, p_d, nE + 1);
  vars.fill(NA_REAL);
  covars.set_size(p_d, p_d, nE + 1);
  covars.fill(NA_REAL);
}

// Constructor for initialisation -------------------------------------------//
Moments::Moments(Upward const & up, Root_State const & root_state, int ntaxa){
  // init all NAs
  int n_E = up.Size() - 1;
  int p_d = root_state.Exp().n_rows;
  *this = Moments(n_E, p_d);
  // Root
  // (No covariance at the root)
  if (root_state.Random()){
    arma::uvec missing_data = up.Missing_Datas(ntaxa);
    arma::mat gamma_inv = inv_sympd(root_state.Var());
    arma::mat S_inv = inv_na(up.Condvars(ntaxa), missing_data, 0);
    vars.slice(ntaxa + 1 - 1) = inv_na(gamma_inv + S_inv, missing_data, arma::datum::inf);
    exps.col(ntaxa + 1 - 1) = vars.slice(ntaxa) * (gamma_inv * root_state.Exp() + S_inv * up.Condexps(ntaxa));
  } else {
    vars.slice(ntaxa + 1 - 1).zeros();
    exps.col(ntaxa + 1 - 1) = root_state.Exp();
  }
}

// Constructor from edge matrix and dimension (USELESS) ---------------------//
Moments::Moments(arma::umat const & ed, int p_d){
  int nE = ed.n_rows; // number of edges
  exps.set_size(p_d, nE + 1);
  exps.fill(NA_REAL);
  vars.set_size(p_d, p_d, nE + 1);
  vars.fill(NA_REAL);
  covars.set_size(p_d, p_d, nE + 1);
  covars.fill(NA_REAL);
}

// Constructor from data and edge matrix (USELESS) --------------------------//
Moments::Moments(arma::mat const & data, arma::umat const & ed){
  // init all NAs
  int p_d = data.n_rows;
  *this = Moments(ed, p_d);
  // Fill exp with data
  int ntaxa = data.n_cols;
  exps.head_cols(ntaxa) = data;
  // Zeros variances / covariance for known values
  for (int i = 0; i < ntaxa; i++){
    for (int l = 0; l < p_d; l++){
      if (! Rcpp::NumericVector::is_na(exps(l, i))){
        vars.subcube(0, l, i, p_d - 1, l, i).zeros(); // column var
        vars.subcube(l, 0, i, l, p_d - 1, i).zeros(); // line var
        covars.subcube(0, l, i, p_d - 1, l, i).zeros(); // column covar
      }
    }
  }
}

// Downward -----------------------------------------------------------------//

void Moments::actualize_downward(Upward const & up, 
                                 Model const & mod,
                                 int edge, int child, int father){
  
  arma::mat Sig_edge = mod.Sigmas(edge);
  arma::mat condvar_child = up.Condvars(child);
  
  arma::mat checkS_inv = inv_sympd(condvar_child + Sig_edge);
  arma::mat S_checkS_inv = condvar_child * checkS_inv;
  
  arma::mat Q_edge = mod.Qs(edge);
  covars.slice(child) = S_checkS_inv * Q_edge * vars.slice(father);
  covars.slice(child) = covars.slice(child).t();
  
  exps.col(child) = S_checkS_inv * (Q_edge * exps.col(father) + mod.Rs(edge));
  exps.col(child) += Sig_edge * checkS_inv * up.Condexps(child);
  
  vars.slice(child) = S_checkS_inv * Sig_edge;
  vars.slice(child) += S_checkS_inv * Q_edge * covars.slice(child);
}

arma::mat compute_Sigma_bar(arma::mat S, arma::mat Sigma_inv, arma::uvec miss, bool isTip){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::uvec na_pos = find(miss == 1); // position of missing
  // Sigma_bar_inv
  arma::mat res_inv = Sigma_inv;
  // Inv
  arma::mat res(size(res_inv), arma::fill::zeros);
  if (isTip){
    arma::uvec no_miss(size(miss), arma::fill::ones);
    no_miss -= miss;
    res = inv_na(res_inv, no_miss, 0); // Tips values have Inf precision
  } else {
    res_inv += inv_na(S, miss, 0);
    res = inv_sympd(res_inv);
  }
  return res;
}

void Moments::actualize_downward_miss(Upward const & up,
                                      Model const & mod,
                                      int edge, int child, int father, int ntaxa){
  arma::uvec missing_data = up.Missing_Datas(child);
  arma::uvec na_pos = find(missing_data == 1);
  
  arma::mat Sig_edge = mod.Sigmas(edge);
  arma::mat Sigma_inv = arma::inv_sympd(Sig_edge);
  arma::mat condvar_child = up.Condvars(child);
  arma::mat Sigma_bar = compute_Sigma_bar(condvar_child, Sigma_inv,
                                          missing_data, (child < ntaxa));

  arma::mat checkS_inv = inv_na(condvar_child + Sig_edge,
                                missing_data, 0);
  arma::mat S_checkS_inv = Sigma_bar * Sigma_inv;
  
  arma::mat Q_edge = mod.Qs(edge);
  arma::mat cova = Q_edge * vars.slice(father);
  cova = S_checkS_inv * cova;
  covars.slice(child) = cova.t();

  exps.col(child) = Sigma_bar * Sigma_inv * (Q_edge * exps.col(father) + mod.Rs(edge));
  exps.col(child) += Sig_edge * checkS_inv * up.Condexps(child);

  vars.slice(child) = Sigma_bar;
  vars.slice(child) += S_checkS_inv * Q_edge * covars.slice(child);
}

void Moments::downward(Upward const & up, Model const & mod,
                       arma::umat const & ed, int ntaxa) {
  int nEdges = ed.n_rows;
  for (int i = nEdges - 1; i >= 0; i--){ // Loop on the edges (rows of ed)
    int father = ed(i, 0) - 1; // Father node of the edge
    int child = ed(i, 1) - 1; // Child edge of the edge
    if (arma::sum(up.Missing_Datas(child)) > 0){ // Some missing data
      (*this).actualize_downward_miss(up, mod, i, child, father, ntaxa);
    } else {
      (*this).actualize_downward(up, mod, i, child, father);
    }
  }
}

// Access to fields ---------------------------------------------------------//
// arma::mat Moments::Exps() const {
//   return exps;
// }
// arma::cube Moments::Vars() const {
//   return vars;
// }
// arma::cube Moments::Covars() const {
//   return covars;
// }

Rcpp::List Moments::exportMoments2R() const {
  return Rcpp::List::create(Rcpp::Named("expectations") = exps,
                            Rcpp::Named("variances") = vars,
                            Rcpp::Named("covariances") = covars);
}

//---------------------------------------------------------------------------//
// Function per se ----------------------------------------------------------//
//---------------------------------------------------------------------------//

Rcpp::List upward_downward_mod(arma::mat const & data, arma::umat const & ed,
                               Model const & mod,
                               Rcpp::List root_state_list) {
  // Numbers
  int nE = ed.n_rows;
  int ntaxa = data.n_cols;
  int p_d = data.n_rows;
  // Construct objects and initialize
  Upward upw(data, nE);
  // Rcpp::Rcout << "J'ai construit l'Upward !" << std::endl;
  // Upward Recursion
  upw.recursion(mod, ed, p_d, ntaxa);
  // Rcpp::Rcout << "J'ai fait l'Upward !" << std::endl;
  // Likelihood Computation
  Root_State root_state = Root_State(root_state_list);
  double logLik = upw.Log_Likelihood(root_state, ntaxa);
  // Rcpp::Rcout << "J'ai calculé la vraisemblance !" << std::endl;
  // Downward Init
  Moments mom(upw, root_state, ntaxa);
  // Rcpp::Rcout << "J'ai construits les moments !" << std::endl;
  mom.downward(upw, mod, ed, ntaxa);
  // Rcpp::Rcout << "J'ai fait le downward !" << std::endl;
  Rcpp::List condLaw = mom.exportMoments2R();
  
  // Return to R in list format
  return Rcpp::List::create(Rcpp::Named("log_likelihood_old") = logLik,
                            Rcpp::Named("conditional_law_X") = condLaw);
}

// [[Rcpp::export]]
Rcpp::List upward_downward_BM(arma::mat const & data, arma::umat const & ed,
                              arma::mat const & Delta, arma::mat const & Variance,
                              arma::vec const & edge_length,
                              Rcpp::List root_state_list) {
  // BM
  Model mod(Delta, Variance, edge_length);
  Rcpp::List res = upward_downward_mod(data, ed, mod, root_state_list);
  return res;
}

// [[Rcpp::export]]
Rcpp::List upward_downward_OU(arma::mat const & data, arma::umat const & ed,
                              arma::mat const & Beta, arma::mat const & Stationary_Var,
                              arma::vec const & edge_length, arma::mat const & Alpha,
                              Rcpp::List root_state_list) {
  // OU
  Model mod(Beta, Stationary_Var, edge_length, Alpha);
  Rcpp::List res = upward_downward_mod(data, ed, mod, root_state_list);
  return res;
}

double log_likelihood_mod(arma::mat const & data, arma::umat const & ed,
                          Model const & mod,
                          Rcpp::List root_state_list) {
  // Numbers
  int nE = ed.n_rows;
  int ntaxa = data.n_cols;
  int p_d = data.n_rows;
  // Construct objects and initialize
  Upward upw(data, nE);
  // Upward Recursion
  upw.recursion(mod, ed, p_d, ntaxa);
  // Likelihood Computation
  Root_State root_state = Root_State(root_state_list);
  double logLik = upw.Log_Likelihood(root_state, ntaxa);
  // Result
  return logLik;
}

// [[Rcpp::export]]
double log_likelihood_BM(arma::mat const & data, arma::umat const & ed,
                         arma::mat const & Delta, arma::mat const & Variance,
                         arma::vec const & edge_length,
                         Rcpp::List root_state_list) {
  // BM
  Model mod(Delta, Variance, edge_length);
  double res = log_likelihood_mod(data, ed, mod, root_state_list);
  return res;
}

// [[Rcpp::export]]
double log_likelihood_OU(arma::mat const & data, arma::umat const & ed,
                         arma::mat const & Beta, arma::mat const & Stationary_Var,
                         arma::vec const & edge_length, arma::mat const & Alpha,
                         Rcpp::List root_state_list) {
  // OU
  Model mod(Beta, Stationary_Var, edge_length, Alpha);
  double res = log_likelihood_mod(data, ed, mod, root_state_list);
  return res;
}

/*** R
library(ape)
library(TreeSim)
library(Matrix)
set.seed(17920902)
ntaxa = 500
tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0,
                        age = 1, mrca = TRUE)[[1]]
tree <- reorder(tree, order = "postorder")
p <- 6
variance <- matrix(0.8, p, p) + diag(0.2, p, p)
independent <- FALSE
root.state <- list(random = FALSE,
                   value.root = rep(1, p),
                   exp.root = NA,
                   var.root = NA)
shifts = list(edges = c(12, 102, 15, 146),
              values = matrix(2*c(1, 0.5), nrow = p, ncol = 4),
              relativeTimes = 0)
paramsSimu <- list(variance = variance,
                   shifts = shifts,
                   root.state = root.state)
attr(paramsSimu, "p_dim") <- p

source("~/Dropbox/These/Code/Phylogenetic-EM/R/simulate.R")
source("~/Dropbox/These/Code/Phylogenetic-EM/R/generic_functions.R")
source("~/Dropbox/These/Code/Phylogenetic-EM/R/shifts_manipulations.R")
source("~/Dropbox/These/Code/Phylogenetic-EM/R/E_step.R")
X1 <- simulate_internal(tree,
               p = p,
               root.state = root.state,
               process = "BM",
               variance = variance,
               shifts = shifts)

traits <- extract_simulate_internal(X1, where = "tips", what = "state")
nMiss <- floor(ntaxa * p * 0.5)
miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
chars <- (miss - 1) %% p + 1
tips <- (miss - 1) %/% p + 1
for (i in 1:nMiss){
  traits[chars[i], tips[i]] <- NA
}
# traits[2, 3] <- traits[2, 17] <- traits[2, 18] <- traits[, 6] <- NA

## Log lik old way
miss <- as.vector(is.na(traits))
masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
masque_data[1:(p * ntaxa)] <- !miss

times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)
T_tree <- incidence.matrix(tree)
U_tree <- incidence.matrix.full(tree)
h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
root_edge_length <- 0
if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)

tmp_old <- wrapper_E_step(phylo = tree,
                          times_shared = times_shared,
                          distances_phylo = distances_phylo,
                          process = "BM",
                          params_old = paramsSimu,
                          masque_data = masque_data,
                          F_moments = F_moments,
                          independent = independent,
                          Y_data_vec_known = as.vector(traits[!miss]),
                          miss = miss,
                          Y_data = traits,
                          U_tree = U_tree,
                          compute_E = compute_E.simple)

ll_old <- tmp_old$log_likelihood_old
conditional_law_X_old <- tmp_old$conditional_law_X

tmp_new <- wrapper_E_step(phylo = tree,
                          times_shared = NULL,
                          distances_phylo = NULL,
                          process = "BM",
                          params_old = paramsSimu,
                          masque_data = masque_data,
                          F_moments = F_moments,
                          independent = independent,
                          Y_data_vec_known = as.vector(traits[!miss]),
                          miss = miss,
                          Y_data = traits,
                          U_tree = U_tree,
                          compute_E = compute_E.upward_downward)

ll_new <- tmp_new$log_likelihood_old
conditional_law_X_new <- tmp_new$conditional_law_X

all.equal(ll_new, as.vector(ll_old))
all.equal(conditional_law_X_old$expectations, conditional_law_X_new$expectations)
all.equal(conditional_law_X_old$variances, conditional_law_X_new$variances)
cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))
all.equal(cov_old, cov_new)


# Delta <- shifts.list_to_matrix(tree, shifts)
# edge_length <- tree$edge.length
# 
# up <- upward_test(traits, tree$edge)
# 
# model <- create_model(Delta, variance, edge_length, tree$edge)
# # # create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 3))
# # # create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 6))
# #
# # # upward_test(traits, tree$edge, 30)
# #
# ll_new <- log_likelihood(traits, tree$edge, Delta, variance, edge_length,
#                          root.state)
# all.equal(ll_new, as.vector(ll_old))
# #
# conditional_law_X_new <- upward_downward_BM(traits, tree$edge, Delta, variance, edge_length, root.state)
# all.equal(conditional_law_X_old, conditional_law_X_new)
# cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
# cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))
# all.equal(cov_old, cov_new)
# 
library(microbenchmark)
microbenchmark(
  # wrapper_E_step(phylo = tree,
  #                             times_shared = times_shared,
  #                             distances_phylo = distances_phylo,
  #                             process = "BM",
  #                             params_old = paramsSimu,
  #                             masque_data = masque_data,
  #                             F_moments = F_moments,
  #                             independent = independent,
  #                             Y_data_vec_known = as.vector(traits[!miss]),
  #                             miss = miss,
  #                             Y_data = traits,
  #                             U_tree = U_tree,
  #                             compute_E = compute_E.simple),
               wrapper_E_step(phylo = tree,
                              process = "BM",
                              params_old = paramsSimu,
                              masque_data = masque_data,
                              independent = independent,
                              Y_data = traits,
                              compute_E = compute_E.upward_downward),
               times = 10)
*/
