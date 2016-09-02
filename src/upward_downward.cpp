# include <RcppArmadillo.h>
// # include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

# include "upward_downward.h"

//---------------------------------------------------------------------------//
// Class Model --------------------------------------------------------------//
//---------------------------------------------------------------------------//

// Default Constructor ------------------------------------------------------//
Model_Node::Model_Node(){
  r.set_size(0);
  q.set_size(0, 0);
  sigma.set_size(0, 0);
}

Model::Model(int siz){
  size = siz;
  mod = new Model_Node [size];
  for (int i = 0; i < size; i++){
    mod[i] = Model_Node();
  }
}

// Constructor for BM Model _Node -------------------------------------------//
Model_Node::Model_Node(arma::vec const & delta, arma::mat const & Variance,
                       double const & edge_length){
  int p_d = Variance.n_rows;
  // R
  // r = delta( no_miss );
  r = delta;
  // Q
  q.eye(p_d, p_d);
  // q = q.rows(no_miss);
  sigma = edge_length * Variance;
  // sigma = sigma(no_miss, no_miss);
}

// Constructor for BM Model -------------------------------------------------//
Model::Model(arma::mat const & Delta, arma::mat const & Variance,
             arma::vec const & edge_length, arma::umat ed){
  int p_d = Delta.n_rows;
  // int ntaxa = data.n_cols; // number of taxa
  int nEdges = edge_length.n_rows; // number of edges
  // Allocate Object
  size = nEdges;
  mod = new Model_Node [size];
  // Fill Model
  for (int i = 0; i < size; i++){
    arma::vec delta = Delta.col(i);
    // int child = ed(i, 1) - 1;
    // arma::uvec no_miss = arma::linspace<arma::uvec>(0, p_d - 1, p_d);
    // if (child < ntaxa){
    //   no_miss = arma::find_finite(data.col(child));
    // }
    mod[i] = Model_Node(delta, Variance, edge_length(i));
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

// [[Rcpp::export]]
Rcpp::List create_model(arma::mat const & Delta, arma::mat const & Variance,
                        arma::vec const & edge_length,
                        arma::umat ed, int i) {
  // Construct object and initialize
  Model mod(Delta, Variance, edge_length, ed);
  
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

// Allocate to fields Model -------------------------------------------------//
void Model::allocate_edge(int i, Model_Node mod_node) {
  mod[i] = mod_node;
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
  cst = NA_REAL;
  condexp.set_size(0);
  condvar.set_size(0, 0);
  missing_data.set_size(0);
}
// Constructor from dimension -----------------------------------------------//
Upward_Node::Upward_Node(int p_d){
  cst = NA_REAL;
  condexp.set_size(p_d);
  condexp.fill(NA_REAL);
  condvar.set_size(p_d, p_d);
  condvar.fill(NA_REAL);
  missing_data.set_size(p_d);
  missing_data.zeros(); // all missing
}

// Access to fields ---------------------------------------------------------//
double Upward_Node::Cst() const {
  return cst;
}
arma::vec Upward_Node::Condexp() const {
  return condexp;
}
arma::mat Upward_Node::Condvar() const {
  return condvar;
}
arma::uvec Upward_Node::Missing_Data() const {
  return missing_data;
}

// Allocate fields ----------------------------------------------------------//
void Upward_Node::allocate_cst(double c) {
  cst = c;
}
void Upward_Node::cst_ones() {
  cst = 1;
}
void Upward_Node::allocate_condexp(arma::vec exp) {
  condexp = exp;
}
void Upward_Node::allocate_condvar(arma::mat var) {
  condvar = var;
}
void Upward_Node::condvar_zeros() {
  condvar.zeros();
}
void Upward_Node::allocate_missing_data(arma::uvec miss_data) {
  missing_data = miss_data;
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
  int p_d = data.n_rows;
  // Fill data
  int ntaxa = data.n_cols;
  for (int i = 0; i < ntaxa; i++){
    // int p_d_tip = findDimensionsNAs(data.col(i));
    up[i] = Upward_Node(p_d);
    // Nas positions
    arma::uvec na_position = arma::find_nonfinite(data.col(i));
    arma::uvec miss_data(p_d);
    miss_data.zeros();
    miss_data(na_position).ones(); // put one when missing
    up[i].allocate_missing_data(miss_data);
    // Fill exp with data
    arma::vec data_col(data.col(i));
    data_col(na_position).zeros(); // replace na with zero
    up[i].allocate_condexp(data_col);
    // Fill variances with zeros
    arma::mat var;
    var.zeros(p_d, p_d);
    var.rows(na_position).fill(arma::datum::inf);
    var.cols(na_position).fill(arma::datum::inf);
    up[i].allocate_condvar(var);
    // Fill constants with right coef
    int nMiss = arma::sum(miss_data);
    up[i].allocate_cst(nMiss * std::log(2 * arma::datum::pi) / 2);
  }
  // Rest with NAs
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

// Log Likelihood---------------------------------------------------------//

double Log_Likelihood_Gauss(arma::vec mean, arma::mat var, arma::vec point){
  int p_d = var.n_rows;
  double res = p_d * std::log(2 * arma::datum::pi) + std::log(arma::det(var));
  arma::vec diff = point - mean;
  res += as_scalar(diff.t() * arma::inv_sympd(var) * diff);
  // Rcpp::Rcout << "ll " << -res/2 << std::endl;
  return -res/2;
}

double Upward::Log_Likelihood(Root_State root_state, int ntaxa) const {
  double res = up[ntaxa + 1 - 1].Cst();
  // Rcpp::Rcout << "log(cst) " << res << std::endl;
  arma::vec mean = up[ntaxa + 1 - 1].Condexp();
  arma::mat var;
  arma::vec point = root_state.Exp();
  if (root_state.Random()){
    var = up[ntaxa + 1 - 1].Condvar() + root_state.Var();
  } else {
    var = up[ntaxa + 1 - 1].Condvar();
  }
   // Rcpp::Rcout << "var: " << var << std::endl;
  res += Log_Likelihood_Gauss(mean, var, point);
  return res;
}

// Allocate to fields Upward ------------------------------------------------//
void Upward::allocate_node(int i, Upward_Node up_node) {
  up[i] = up_node;
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

// Model & findModelEdges(int father, arma::umat const & ed, Model const & mod){
//   arma::uvec edges_children = arma::find(ed.col(0) == father + 1);
//   int nChild = edges_children.n_rows;
//   Model *res = new Model(nChild);
//   for (int i = 0; i < nChild; i++){
//     (*res).allocate_edge(i, mod.Mod(edges_children(i)));
//   }
//   return *res;
// }

Upward_Node & actualize_upward_missing(Upward_Node const &  up_child,
                                       Model_Node const & mod_edge,
                                       int p_d){
  Upward_Node *res = new Upward_Node(p_d);
  arma::mat check_S_inv = arma::inv_sympd(up_child.Condvar() + mod_edge.Sigma());
  arma::mat tQ_checkS = mod_edge.Q().t() * check_S_inv;
  arma::mat tQ_checkS_Q_inv = arma::inv_sympd(tQ_checkS * mod_edge.Q());
  arma::mat Q_pseudo_inv = tQ_checkS_Q_inv * tQ_checkS;
  arma::vec diff_exp = (up_child.Condexp() - mod_edge.R());
  arma::mat Qpseudo_diff_exp = Q_pseudo_inv * diff_exp;
  double diff_diff_diff = arma::as_scalar(diff_exp.t() * (tQ_checkS.t() * Q_pseudo_inv - check_S_inv) * diff_exp);
  double dets = arma::det(check_S_inv) * arma::det(tQ_checkS_Q_inv);
  double constant = std::exp(- diff_diff_diff / 2) * std::sqrt(dets);
  constant *= up_child.Cst() * std::sqrt(std::pow(2 * arma::datum::pi, (int)(p_d - diff_exp.n_rows)));
  
  (*res).allocate_condvar(tQ_checkS_Q_inv);
  (*res).allocate_condexp(Qpseudo_diff_exp);
  (*res).allocate_cst(constant);
  
  return *res;
}

arma::mat inv_na(arma::mat const & S, arma::uvec const & miss, double ff){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::mat S_sub = S(non_na_pos, non_na_pos);
  // Rcpp::Rcout << "miss " << miss << std::endl;
  // Rcpp::Rcout << "S_sub " << S_sub << std::endl;
  S_sub = arma::inv_sympd(S_sub);
  arma::mat res(size(S));
  res.fill(ff);
  res(non_na_pos, non_na_pos) = S_sub;
  return res;
}

double det_na(arma::mat const & S, arma::uvec const & miss){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::mat S_sub = S(non_na_pos, non_na_pos);
  double res = arma::det(S_sub);
  return res;
}

arma::mat crossprod_na(arma::mat const & Q, arma::mat const & cS, arma::uvec const & miss){
  arma::mat S_inv = inv_na(cS, miss, 0);
  arma::mat res = Q.t() * S_inv * Q;
  res = inv_na(res, miss, arma::datum::inf);
  return res;
}

arma::vec prod_na(arma::mat const & S, arma::vec const & m, arma::uvec const & miss, double ff){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::vec res = ff * m;
  res(non_na_pos) = S(non_na_pos, non_na_pos) * m(non_na_pos);
  return res;
}

arma::vec prod_na_non(arma::mat const & S, arma::vec const & m, arma::uvec const & miss){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::uvec na_pos = find(miss == 1); // position of missing
  arma::vec res(size(m), arma::fill::zeros);
  res(non_na_pos) = S(non_na_pos, non_na_pos) * m(non_na_pos);
  res(na_pos) = S(na_pos, na_pos) * m(na_pos);
  return res;
}

// Two matrices of same dimension
arma::mat prod_na_mat(arma::mat const & S, arma::mat const & M, arma::uvec const & miss, double ff){
  arma::uvec non_na_pos = find(miss == 0); // position of non missing
  arma::mat res = ff * M; 
  res(non_na_pos, non_na_pos) = S(non_na_pos, non_na_pos) * M(non_na_pos, non_na_pos);
  return res;
}

Upward_Node & actualize_upward_simple(Upward_Node const &  up_child,
                                      Model_Node const & mod_edge,
                                      int p_d){
  Upward_Node *res = new Upward_Node(p_d);
  arma::uvec missing_data = up_child.Missing_Data();
  arma::mat check_S = up_child.Condvar() + mod_edge.Sigma();
  // Rcpp::Rcout << "check S " << check_S << std::endl;
  arma::mat Q_inv = arma::inv_sympd(mod_edge.Q());
  arma::mat tQ_checkS_Q_inv = crossprod_na(mod_edge.Q(), check_S, missing_data);
  // Rcpp::Rcout << "isisng up " << up_child.Missing_Data() << std::endl;
  arma::vec Q_inv_diff = prod_na(Q_inv, (up_child.Condexp() - mod_edge.R()), missing_data, 0);
    // Q_inv * (up_child.Condexp() - mod_edge.R());
  double constant = - std::log(det_na(mod_edge.Q(), missing_data)) / 2 + up_child.Cst();
  
  (*res).allocate_condvar(tQ_checkS_Q_inv);
  (*res).allocate_condexp(Q_inv_diff);
  (*res).allocate_cst(constant);
  (*res).allocate_missing_data(up_child.Missing_Data());
  
  return *res;
}

Upward & actualize_upward_children(arma::uvec const & child_nodes,
                                   arma::uvec const & child_edges,
                                   Upward & upw, Model const & mod,
                                   int p_d, int ntaxa){
  int nChild = child_nodes.n_rows;
  Upward *res = new Upward(nChild, p_d);
  for (int i = 0; i < nChild; i++){
    // if (upw.Up(child_nodes(i)).Condvar().n_rows < p_d){
    //   // Rcpp::Rcout << "missing; child: " << child_nodes(i)+1 << " edge: " << child_edges(i)+1 << std::endl;
    //   (*res).allocate_node(i, actualize_upward_missing(upw.Up(child_nodes(i)),
    //                                                    mod.Mod(child_edges(i)),
    //                                                    p_d));
    // } else {
      // Rcpp::Rcout << "no missing; child: " << child_nodes(i)+1 << " edge: " << child_edges(i)+1 << std::endl;
      (*res).allocate_node(i, actualize_upward_simple(upw.Up(child_nodes(i)),
                                                      mod.Mod(child_edges(i)),
                                                      p_d));
    // }
  }
  return *res;
}

Upward_Node & merge_upward(Upward const & up_child){
  int nChild = up_child.Size();
  int p_d = up_child.Up(0).Condexp().n_rows;
  Upward_Node *res = new Upward_Node(p_d);
  
  arma::mat S_sum;
  S_sum.zeros(p_d, p_d);
  arma::vec S_m_sum;
  S_m_sum.zeros(p_d);
  double m_S_m_sum = 0;
  double log_det_sum = 0;
  double cst_sum = 0;
  arma::uvec merge_missing(p_d);
  merge_missing.ones();
  
  for (int i = 0; i < nChild; i++){
    arma::mat S_inv = inv_na(up_child.Up(i).Condvar(), up_child.Up(i).Missing_Data(), 0);
    arma::vec S_inv_m = S_inv * up_child.Up(i).Condexp();
    S_sum += S_inv;
    S_m_sum += S_inv_m;
    m_S_m_sum += arma::as_scalar(up_child.Up(i).Condexp().t() * S_inv_m);
    log_det_sum += - std::log(det_na(up_child.Up(i).Condvar(), up_child.Up(i).Missing_Data())) / 2;
    cst_sum += up_child.Up(i).Cst();
    merge_missing = merge_missing % up_child.Up(i).Missing_Data(); // intersection of missing
  }
  
  arma::mat S_bar = inv_na(S_sum, merge_missing, arma::datum::inf);
  arma::vec m_bar = prod_na(S_bar, S_m_sum, merge_missing, 0);
  double constant = - (nChild - 1) * p_d * std::log(2 * arma::datum::pi) / 2;
  constant += std::log(det_na(S_bar, merge_missing)) / 2 + log_det_sum;
  constant += - arma::as_scalar(m_S_m_sum -  m_bar.t() * S_sum * m_bar) / 2;
  constant += cst_sum;
  // // number of missing values
  // int nMiss = arma::sum(merge_missing);
  // constant *= std::sqrt(pow(2 * arma::datum::pi, -nMiss));
  
  (*res).allocate_condvar(S_bar);
  (*res).allocate_condexp(m_bar);
  (*res).allocate_cst(constant);
  (*res).allocate_missing_data(merge_missing);
  
  return *res;
}

void Upward::recursion(Model const & mod, arma::umat const & ed,
                       int p_d, int ntaxa) {
  int nEdges = ed.n_rows;
  for (arma::uword i = 0; i < nEdges; i++){ // Loop on the edges (rows of ed)
    int father = ed(i, 0) - 1;
    if (! arma::is_finite(up[father].Cst())){// This node has not already been visited.
      // Find children of the node
      arma::uvec child_nodes = findChildren(father, ed) - 1;
      // Find Models associated to edges
      arma::uvec child_edges = findEdges(father, ed);
      // Create an array with actualized quantities at each child
      // Rcpp::Rcout << "children are" << child_nodes + 1 << std::endl;
      Upward up_child = actualize_upward_children(child_nodes, child_edges,
                                                  *this, mod,
                                                  p_d, ntaxa);
      // Compute the quantity for the node
      up[father] = merge_upward(up_child);
    }
  }
}

// Export to R (test) -------------------------------------------------------//
Rcpp::List Upward::exportUpward2R(int i) const {
  if ((i >= size) || (i < 0)){
    Rcpp::stop("In exportUpward, index out of bounds.");
  }
  return Rcpp::List::create(Rcpp::Named("cst") = up[i].Cst(),
                            Rcpp::Named("condexp") = up[i].Condexp(),
                            Rcpp::Named("condvar") = up[i].Condvar(),
                            Rcpp::Named("missing_data") = up[i].Missing_Data());
}

// [[Rcpp::export]]
Rcpp::List upward_test(arma::mat const & data, arma::umat const & ed, int i) {
  int nE = ed.n_rows;
  Upward upw(data, nE);
  return upw.exportUpward2R(i-1);
}

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

// Acces t fields -----------------------------------------------------------//
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
  // nEdges = nE;
  // p_dim = p_d;
  // edge.set_size(nE, 2);
  // edge.fill(NA_REAL);
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
    arma::uvec missing_data = up.Up(ntaxa).Missing_Data();
    arma::mat gamma_inv = inv_sympd(root_state.Var());
    arma::mat S_inv = inv_na(up.Up(ntaxa).Condvar(), missing_data, 0);
    vars.slice(ntaxa + 1 - 1) = inv_na(gamma_inv + S_inv, missing_data, arma::datum::inf);
    exps.col(ntaxa + 1 - 1) = vars.slice(ntaxa) * (gamma_inv * root_state.Exp() + S_inv * up.Up(ntaxa).Condexp());
  } else {
    vars.slice(ntaxa + 1 - 1).zeros();
    exps.col(ntaxa + 1 - 1) = root_state.Exp();
  }
}

// Constructor from edge matrix and dimension (USELESS) ---------------------//
Moments::Moments(arma::umat const & ed, int p_d){
  int nE = ed.n_rows; // number of edges
  // edge = ed;
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

// Downward -----------------------------------------------------------------//

void Moments::actualize_downward(Upward_Node const & up_child, 
                                 Model_Node const & mod_edge,
                                 int child, int father){
  arma::mat checkS_inv = inv_sympd(up_child.Condvar() + mod_edge.Sigma());
  arma::mat S_checkS_inv = up_child.Condvar() * checkS_inv;
  
  covars.slice(child) = S_checkS_inv * mod_edge.Q() * vars.slice(father);
  covars.slice(child) = covars.slice(child).t();
  
  exps.col(child) = S_checkS_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R());
  exps.col(child) += mod_edge.Sigma() * checkS_inv * up_child.Condexp();
  
  vars.slice(child) = S_checkS_inv * mod_edge.Sigma();
  vars.slice(child) += S_checkS_inv * mod_edge.Q() * covars.slice(child);
}

// void Moments::actualize_downward(Upward_Node const & up_child, 
//                                  Model_Node const & mod_edge,
//                                  int child, int father){
//   arma::uvec missing_data = up_child.Missing_Data();
//   arma::mat checkS_inv = inv_na(up_child.Condvar() + mod_edge.Sigma(), missing_data, 0);
//   arma::mat S_checkS_inv = prod_na_mat(up_child.Condvar(), checkS_inv, missing_data, arma::datum::inf);
//   covars.slice(child) = S_checkS_inv * mod_edge.Q() * vars.slice(father);
//   exps.col(child) = S_checkS_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R());
//   exps.col(child) += prod_na(mod_edge.Sigma(), prod_na(checkS_inv, up_child.Condexp(), missing_data, 0), missing_data, 0);
//   vars.slice(child) = S_checkS_inv * mod_edge.Sigma();
//   vars.slice(child) += S_checkS_inv * mod_edge.Q() * covars.slice(child).t();
// }

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
  // // data
  // arma::mat S_d = S(non_na_pos, non_na_pos);
  // arma::mat Sigma_d = Sigma(non_na_pos, non_na_pos);
  // // mixte
  // // arma::mat S_d_na = S(non_na_pos, na_pos);
  // arma::mat Sigma_d_na = Sigma(non_na_pos, na_pos);
  // // nas
  // // arma::mat S_na = S(na_pos, na_pos);
  // arma::mat Sigma_na = Sigma(na_pos, na_pos);
  // // Inversions Sigma
  // arma::mat temp = Sigma_d_na.t() * arma::inv_sympd(Sigma_d) * Sigma_d_na;
  // arma::mat Sigma_inv_na_inv = Sigma_na - Sigma_d_na.t() * temp;
  // arma::mat Sigma_ing_d_na = - temp * Sigma_inv_na_inv;
  // // Block computations
  // // // Data
  // // arma::mat check_S_inv_d = arma::inv_sympd(S_d + Sigma_d);
  // // arma::mat S_d_check_S_inv_d = S_d * check_S_inv_d;
  // // res(non_na_pos, non_na_pos) = S_d_check_S_inv_d * Sigma_d;
  // // // Mixte
  // // arma::mat cross = S_d_check_S_inv_d * Sigma_d_na;
  // // res(non_na_pos, na_pos) = cross;
  // // res(na_pos, non_na_pos) = cross.t();
  // // NAs
  // res(na_pos, na_pos) = Sigma_na - Sigma_d_na.t() * arma::inv_sympd(Sigma_d) * Sigma_d_na;
  // // res(na_pos, na_pos) = S_d_na.t() * check_S_inv_d * Sigma_d_na;
  // // Rcpp::Rcout << "S_d_na " << S_d_na << std::endl;
  // // Rcpp::Rcout << "check_S_inv_d " << check_S_inv_d << std::endl;
  // // Rcpp::Rcout << "Sigma_d_na " << Sigma_d_na << std::endl;
  // return res;
}

void Moments::actualize_downward_miss(Upward_Node const & up_child,
                                      Model_Node const & mod_edge,
                                      int child, int father, int ntaxa){
  arma::uvec missing_data = up_child.Missing_Data();
  arma::uvec na_pos = find(missing_data == 1);
  
  arma::mat Sigma_inv = arma::inv_sympd(mod_edge.Sigma());
  arma::mat Sigma_bar = compute_Sigma_bar(up_child.Condvar(), Sigma_inv,
                                          missing_data, (child < ntaxa));

  arma::mat checkS_inv = inv_na(up_child.Condvar() + mod_edge.Sigma(),
                                missing_data, 0);
  arma::mat S_checkS_inv = Sigma_bar * Sigma_inv;

  arma::mat cova = mod_edge.Q() * vars.slice(father);
  cova = S_checkS_inv * cova;
  covars.slice(child) = cova.t();

  exps.col(child) = Sigma_bar * Sigma_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R());
  exps.col(child) += mod_edge.Sigma() * checkS_inv * up_child.Condexp();

  vars.slice(child) = Sigma_bar;
  vars.slice(child) += S_checkS_inv * mod_edge.Q() * covars.slice(child);
  // Rcpp::Rcout << "child " << child << std::endl;
  // Rcpp::Rcout << "missing_data " << missing_data << std::endl;
  // Rcpp::Rcout << "checkS_inv " << checkS_inv << std::endl;
  // Rcpp::Rcout << "up_child.Condvar() " << up_child.Condvar() << std::endl;
  // Rcpp::Rcout << "S_checkS_inv " << S_checkS_inv << std::endl;
  // Rcpp::Rcout << "Sigma_bar " << Sigma_bar << std::endl;
  // Rcpp::Rcout << "Sigma_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R()) " << Sigma_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R()) << std::endl;
  // Rcpp::Rcout << "mod_edge.Sigma() " << mod_edge.Sigma() << std::endl;
  // Rcpp::Rcout << "exp 1 " << Sigma_bar * Sigma_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R()) << std::endl;
  // Rcpp::Rcout << "exp 2 " <<  mod_edge.Sigma() *checkS_inv * up_child.Condexp() << std::endl;
  // Rcpp::Rcout << "exp" << exps.col(child) << std::endl;
}

// void Moments::actualize_downward_miss(Upward_Node const & up_child, 
//                                       Model_Node const & mod_edge,
//                                       int child, int father){
//   arma::uvec missing_data = up_child.Missing_Data();
//   
//   if (child == 10-1){
//     Rcpp::Rcout << "mod_edge.Sigma() " << mod_edge.Sigma() << std::endl; 
//   }
//   
//   arma::mat Sigma_inv = arma::inv_sympd(mod_edge.Sigma());
//   
//   arma::uvec na_pos = find(missing_data == 1);
//   arma::mat Sigma_bar(size(mod_edge.Sigma()), arma::fill::zeros);
//   Sigma_bar(na_pos, na_pos) = mod_edge.Sigma()(na_pos, na_pos);
//   Rcpp::Rcout << "child " << child + 1 << std::endl; 
//   arma::mat condvar_inv = inv_na(up_child.Condvar(), missing_data, 0);
//   Rcpp::Rcout << "child " << child + 1 << std::endl;
//   Sigma_bar += inv_na(condvar_inv + Sigma_inv, missing_data, 0);
//   Rcpp::Rcout << "child " << child + 1 << std::endl; 
//   
//   arma::mat S_checkS_inv = Sigma_bar * Sigma_inv;
//   covars.slice(child) = S_checkS_inv * mod_edge.Q() * vars.slice(father);
//   
//   exps.col(child) = S_checkS_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R());
//   arma::vec values = up_child.Condexp();
//   values.rows(na_pos).zeros();
//   exps.col(child) += values;
//   
//   vars.slice(child) = Sigma_bar;
//   vars.slice(child) += S_checkS_inv * mod_edge.Q() * covars.slice(child).t();
//   
//   if (child == 3-1){
//     Rcpp::Rcout << "missing_data " << missing_data << std::endl; 
//     Rcpp::Rcout << "up_child.Condvar() " << up_child.Condvar() << std::endl; 
//     Rcpp::Rcout << "Sigma_bar " << Sigma_bar << std::endl;  
//     Rcpp::Rcout << "S_checkS_inv " << S_checkS_inv << std::endl;  
//     Rcpp::Rcout << "mod_edge.Sigma() " << mod_edge.Sigma() << std::endl;  
//     // Rcpp::Rcout << "exp 1 " << S_checkS_inv * (mod_edge.Q() * exps.col(father) + mod_edge.R()) << std::endl;
//     // Rcpp::Rcout << "exp 2 " << mod_edge.Sigma() * checkS_inv * up_child.Condexp() << std::endl;
//     Rcpp::Rcout << "exp" << exps.col(child) << std::endl;
//   }
// }


void Moments::downward(Upward const & up, Model const & mod,
                       arma::umat const & ed, int ntaxa) {
  int nEdges = ed.n_rows;
  for (int i = nEdges - 1; i >= 0; i--){ // Loop on the edges (rows of ed)
    // Rcpp::Rcout << "i " << i << std::endl;
    int father = ed(i, 0) - 1; // Father node of the edge
    int child = ed(i, 1) - 1; // Child edge of the edge
    // Rcpp::Rcout << "child " << child << std::endl;
    if (arma::sum(up.Up(child).Missing_Data()) > 0){ // Some missing data
      (*this).actualize_downward_miss(up.Up(child), mod.Mod(i), child, father, ntaxa);
    } else {
      (*this).actualize_downward(up.Up(child), mod.Mod(i), child, father);
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
  return Rcpp::List::create(Rcpp::Named("expectations") = exps,
                            Rcpp::Named("variances") = vars,
                            Rcpp::Named("covariances") = covars);
}

//---------------------------------------------------------------------------//
// Function per se ----------------------------------------------------------//
//---------------------------------------------------------------------------//

// [[Rcpp::export]]
Rcpp::List upward_downward(arma::mat const & data, arma::umat const & ed,
                           arma::mat const & Delta, arma::mat const & Variance,
                           arma::vec const & edge_length,
                           Rcpp::List root_state_list) {
  // Numbers
  int nE = ed.n_rows;
  int ntaxa = data.n_cols;
  int p_d = data.n_rows;
  // Construct objects and initialize
  //Moments mom(data, ed);
  Upward upw(data, nE);
  // BM
  Model mod(Delta, Variance, edge_length, ed);
  // Upward Recursion
  upw.recursion(mod, ed, p_d, ntaxa);
  // Likelihood Computation
  Root_State root_state = Root_State(root_state_list);
  double logLik = upw.Log_Likelihood(root_state, ntaxa);
  // Downward Init
  Moments mom(upw, root_state, ntaxa);
  mom.downward(upw, mod, ed, ntaxa);
  Rcpp::List condLaw = mom.exportMoments2R();
  
  // Return to R in list format
  return Rcpp::List::create(Rcpp::Named("log_likelihood_old") = logLik,
                            Rcpp::Named("conditional_law_X") = condLaw);
}

// [[Rcpp::export]]
double log_likelihood(arma::mat const & data, arma::umat const & ed,
                      arma::mat const & Delta, arma::mat const & Variance,
                      arma::vec const & edge_length,
                      Rcpp::List const & root_state_list) {
  // Numbers
  int nE = ed.n_rows;
  int ntaxa = data.n_cols;
  int p_d = data.n_rows;
  // Construct objects and initialize
  //Moments mom(data, ed);
  Upward upw(data, nE);
  // BM
  Model mod(Delta, Variance, edge_length, ed);
  // Upward Recursion
  upw.recursion(mod, ed, p_d, ntaxa);
  
  
  // Return to R in list format
  Root_State root_state = Root_State(root_state_list);
  double res = upw.Log_Likelihood(root_state, ntaxa);
  return res;
  // return upw.exportUpward2R(i-1);
  //return mom.exportMoments2R();
}

// /*** R
// library(ape)
// library(TreeSim)
// library(Matrix)
// set.seed(17920902)
// ntaxa = 50
// tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
//                         age = 1, mrca = TRUE)[[1]]
// tree <- reorder(tree, order = "postorder")
// p <- 4
// variance <- matrix(0.8, p, p) + diag(0.2, p, p)
// independent <- FALSE
// root.state <- list(random = FALSE,
//                    value.root = rep(1, p),
//                    exp.root = NA,
//                    var.root = NA)
// shifts = list(edges = c(12),
//               values=matrix(2*c(1, 0.5), nrow = p),
//               relativeTimes = 0)
// paramsSimu <- list(variance = variance,
//                    shifts = shifts,
//                    root.state = root.state)
// attr(paramsSimu, "p_dim") <- p
// 
// source("~/Dropbox/These/Code/Phylogenetic-EM/R/simulate.R")
// source("~/Dropbox/These/Code/Phylogenetic-EM/R/generic_functions.R")
// source("~/Dropbox/These/Code/Phylogenetic-EM/R/shifts_manipulations.R")
// source("~/Dropbox/These/Code/Phylogenetic-EM/R/E_step.R")
// X1 <- simulate(tree,
//                p = p,
//                root.state = root.state,
//                process = "BM",
//                variance = variance,
//                shifts = shifts)
// 
// traits <- extract.simulate(X1, where = "tips", what = "state")
// nMiss <- floor(ntaxa * p * 0.5)
// miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
// chars <- (miss - 1) %% p + 1
// tips <- (miss - 1) %/% p + 1
// for (i in 1:nMiss){
//   traits[chars[i], tips[i]] <- NA
// }
// # traits[2, 3] <- traits[2, 17] <- traits[2, 18] <- traits[, 6] <- NA
// 
// ## Log lik old way
// miss <- as.vector(is.na(traits))
// masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
// masque_data[1:(p * ntaxa)] <- !miss
// 
// times_shared <- compute_times_ca(tree)
// distances_phylo <- compute_dist_phy(tree)
// T_tree <- incidence.matrix(tree)
// U_tree <- incidence.matrix.full(tree)
// h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
// root_edge_length <- 0
// if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
// F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)
// 
// tmp_old <- wrapper_E_step(phylo = tree,
//                           times_shared = times_shared,
//                           distances_phylo = distances_phylo,
//                           process = "BM",
//                           params_old = paramsSimu,
//                           masque_data = masque_data,
//                           F_moments = F_moments,
//                           independent = independent,
//                           Y_data_vec_known = as.vector(traits[!miss]),
//                           miss = miss,
//                           Y_data = traits,
//                           U_tree = U_tree,
//                           compute_E = compute_E.simple)
// 
// ll_old <- tmp_old$log_likelihood_old
// conditional_law_X_old <- tmp_old$conditional_law_X
// 
// tmp_new <- wrapper_E_step(phylo = tree,
//                           times_shared = times_shared,
//                           distances_phylo = distances_phylo,
//                           process = "BM",
//                           params_old = paramsSimu,
//                           masque_data = masque_data,
//                           F_moments = F_moments,
//                           independent = independent,
//                           Y_data_vec_known = as.vector(traits[!miss]),
//                           miss = miss,
//                           Y_data = traits,
//                           U_tree = U_tree,
//                           compute_E = compute_E.upward_downward)
// 
// ll_new <- tmp_new$log_likelihood_old
// conditional_law_X_new <- tmp_new$conditional_law_X
// 
// all.equal(ll_new, as.vector(ll_old))
// all.equal(conditional_law_X_old$expectations, conditional_law_X_new$expectations)
// all.equal(conditional_law_X_old$variances, conditional_law_X_new$variances)
// cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
// cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))
// all.equal(cov_old, cov_new)
// 
// 
// # # conditional_law_X <- upward_downward(traits, tree$edge)
// # 
// # Delta <- shifts.list_to_matrix(tree, shifts)
// # edge_length <- tree$edge.length
// # # create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 1))
// # # create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 3))
// # # create_model(Delta, variance, edge_length, traits, tree$edge, which(tree$edge[, 2] == 6))
// # 
// # # upward_test(traits, tree$edge, 30)
// # 
// # ll_new <- log_likelihood(traits, tree$edge, Delta, variance, edge_length,
// #                          root.state)
// # all.equal(ll_new, as.vector(ll_old))
// # 
// # conditional_law_X_new <- upward_downward(traits, tree$edge, Delta, variance, edge_length, root.state)
// # all.equal(conditional_law_X_old, conditional_law_X_new)
// # cov_new <- apply(conditional_law_X_new$covariances, 3, function(z) z + t(z))
// # cov_old <- apply(conditional_law_X_old$covariances, 3, function(z) z + t(z))
// # all.equal(cov_old, cov_new)
// 
// library(microbenchmark)
// microbenchmark(wrapper_E_step(phylo = tree,
//                               times_shared = times_shared,
//                               distances_phylo = distances_phylo,
//                               process = "BM",
//                               params_old = paramsSimu,
//                               masque_data = masque_data,
//                               F_moments = F_moments,
//                               independent = independent,
//                               Y_data_vec_known = as.vector(traits[!miss]),
//                               miss = miss,
//                               Y_data = traits,
//                               U_tree = U_tree,
//                               compute_E = compute_E.simple),
//                wrapper_E_step(phylo = tree,
//                               process = "BM",
//                               params_old = paramsSimu,
//                               masque_data = masque_data,
//                               independent = independent,
//                               Y_data = traits,
//                               compute_E = compute_E.upward_downward),
//                times = 10)
// */
