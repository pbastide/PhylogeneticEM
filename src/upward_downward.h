#ifndef updown_h
#define updown_h

//---------------------------------------------------------------------------//
// Class Root_State ---------------------------------------------------------//
//---------------------------------------------------------------------------//
/*
* The root state
*/

class Root_State
{
private:
  
  bool random;
  arma::vec exp;
  arma::mat var;
  
public:
  // Constructors
  Root_State(arma::vec state); // Default non random
  Root_State(arma::vec expe, arma::mat vari); // Default random
  Root_State(bool rootRand, arma::vec rootVal,
             arma::vec rootExp, arma::mat rootVar);
  Root_State(Rcpp::List root_state);
  
  
  // Access to fields
  bool Random() const;
  arma::vec Exp() const;
  arma::mat Var() const;
  
};

//---------------------------------------------------------------------------//
// Class Model --------------------------------------------------------------//
//---------------------------------------------------------------------------//
/*
 * The model is assumed to be written in the following form: 
 *
 * E[X_j | X_pa(j)] = Q_j * X_pa(j) + r_j
 * Var[X_j | X_pa(j)] = Sigma_j
 * 
 * The nEdges entries (in edges order).
 */
class Model
{
private:
  arma::mat rs; // matix of expectations at all the nodes
  arma::cube qs; // 3D array of variances at all the nodes
  arma::cube sigmas; // 3D array of covariances at all nodes/parents. NA at the root.
  
public:
  // Constructors
  Model(int siz);
  // For a BM:
  Model(arma::mat const & Delta, arma::mat const & Variance,
        arma::vec const & edge_length);
  // For an OU:
  Model(arma::mat const & Beta, arma::mat const & Stationary_Var,
        arma::vec const & edge_length, arma::mat const & Alpha);
  
  // Access to fields
  arma::mat Rs() const;
  arma::cube Qs() const;
  arma::cube Sigmas() const;
  
  // Export to R (test only)
  // Rcpp::List exportModel2R() const;
};

//---------------------------------------------------------------------------//
// Class Upward -------------------------------------------------------------//
//---------------------------------------------------------------------------//
/*
 * The quantities needed at the upward phase.
 * This is one instance at one node.
 */
class Upward_Node
{
private:
  
  double cst; // constant. This is the log.
  arma::vec condexp; // vector of expectations
  arma::mat condvar; // variance matrix
  arma::uvec missing_data; // Position of the missing data (bool vector)
  
public:
  // Constructors
  Upward_Node();
  Upward_Node(int p_d);
  
  // Access to fields
  double Cst() const;
  arma::vec Condexp() const;
  arma::mat Condvar() const;
  arma::uvec Missing_Data() const;
  
  // Allocate Fields
  void allocate_cst(double c);
  void cst_ones();
  void allocate_condexp(arma::vec exp);
  void allocate_condvar(arma::mat var);
  void condvar_zeros();
  void allocate_missing_data(arma::uvec miss_data);
  
  
  // Export to R (test only)
  //Rcpp::List exportUpward2R() const;
};



/*
 * These are all the instances, numbered by the nodes.
 */

class Upward
{
private:
  
  arma::vec csts; // constant. This is the log.
  arma::mat condexps; // vector of expectations
  arma::cube condvars; // variance matrix
  arma::umat missing_datas; // Position of the missing data (bool vector)
  
public:
  // Constructors
  Upward(int siz, int p_d);
  Upward(arma::mat const & data, int nE);
  
  // Access to fields
  arma::vec Csts() const;
  arma::mat Condexps() const;
  arma::cube Condvars() const;
  arma::umat Missing_Datas() const;
  int Size() const;
  double Log_Likelihood(Root_State root, int ntaxa) const;
  
  // Allocate Fields
  void allocate_condvars(int i, arma::mat M);
  void allocate_condexps(int i, arma::vec E);
  void allocate_csts(int i, double c);
  void allocate_missing_datas(int i, arma::uvec M);
  
  // Recursion
  void recursion(Model const & mod, arma::umat const & ed,
                 int p_d, int ntaxa);
  
  // Export to R (test only)
  // Rcpp::List exportUpward2R() const;
};

//---------------------------------------------------------------------------//
// Class Moments ------------------------------------------------------------//
//---------------------------------------------------------------------------//
/*
* The moments we are looking for. 
* See doc of E_step function for further details.
*/
class Moments
{
private:
  
  //int nEdges; // number of edges in the tree
  //int p_dim; // dimension of the trait vector
  //arma::umat miss; // a bool matrix size p_dim x ntaxa with TRUE if the data is missing.
  // arma::umat edge; // matrix of edges in post-order
  arma::mat exps; // matix of expectations at all the nodes
  arma::cube vars; // 3D array of variances at all the nodes
  arma::cube covars; // 3D array of covariances at all nodes/parents. NA at the root.
  
public:
  // Constructors
  Moments(int = 0, int = 0);
  Moments(Upward const & up, Root_State const & root_state, int ntaxa);
  Moments(arma::umat const & ed, int p_d);
  Moments(arma::mat const & data, arma::umat const & ed);
  
  // Export to R
  Rcpp::List exportMoments2R() const;
  
  // Access to fields
  arma::mat Exps() const;
  arma::cube Vars() const;
  arma::cube Covars() const;
  
  // Downward
  void actualize_downward(Upward const & up, 
                          Model const & mod,
                          int edge, int child, int father);
  void actualize_downward_miss(Upward const & up, 
                               Model const & mod,
                               int edge, int child, int father, int ntaxa);
  void downward(Upward const & up, Model const & mod,
                arma::umat const & ed, int ntaxa);
};


#endif
