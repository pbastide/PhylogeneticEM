#ifndef updown_h
#define updown_h

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
  arma::umat edge; // matrix of edges in post-order
  arma::mat exps; // matix of expectations at all the nodes
  arma::cube vars; // 3D array of variances at all the nodes
  arma::cube covars; // 3D array of covariances at all nodes/parents. NA at the root.
  
public:
  // Constructors
  Moments(int = 0, int = 0);
  Moments(arma::umat const & ed, int p_d);
  Moments(arma::mat const & data, arma::umat const & ed);
  
  // Export to R
  Rcpp::List exportMoments2R() const;
  
  // Access to fields
  arma::mat Exps() const;
  arma::cube Vars() const;
  arma::cube Covars() const;
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
 * One instance at one edge
 */
class Model_Node
{
private:
  
  arma::vec r; // matix of vectors offsets
  arma::mat q; // 3D array of transmission matrices
  arma::mat sigma; // 3D array of variances matrics
  
public:
  // Constructors
  Model_Node();
  // For a BM:
  Model_Node(arma::vec const & delta, arma::mat const & Variance,
             double const & edge_length, arma::uvec no_miss);
  // For an OU:
  Model_Node(arma::vec const & beta, arma::mat const & Stationary_Var,
             double const & edge_length, arma::mat const & Alpha);
  
  // Access to fields
  arma::vec R() const;
  arma::mat Q() const;
  arma::mat Sigma() const;
  
  // Export to R (test only)
  //Rcpp::List exportModel2R() const;
};

/*
 * The nEdges entries (in edges order).
 */
class Model
{
private:
  
  Model_Node * mod; // Array if Upward_node elements
  unsigned int size; // Size
  
public:
  // Constructors
  // For a BM:
  Model(arma::mat const & Delta, arma::mat const & Variance,
        arma::vec const & edge_length,
        arma::mat data, arma::umat ed);
  // For an OU:
  Model(arma::mat const & Beta, arma::mat const & Stationary_Var,
        arma::vec const & edge_length, arma::mat const & Alpha,
        arma::mat data, arma::umat ed);
  
  // Destructor
  ~Model();
  
  // Deal with missing data
  //Model removeMissing(arma::mat const & miss);
  
  // Access to fields
  Model_Node Mod(int i) const;
  unsigned int Size() const;
  
  
  // Export to R (test only)
  Rcpp::List exportModel2R(int i) const;
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
  
  arma::vec cst; // constant
  arma::vec condexp; // vector of expectations
  arma::mat condvar; // variance matrix
  
public:
  // Constructors
  Upward_Node();
  Upward_Node(int p_d);
  
  // Access to fields
  arma::vec Cst() const;
  arma::vec Condexp() const;
  arma::mat Condvar() const;
  
  // Allocate Fields
  void allocate_cst(arma::vec c);
  void cst_ones();
  void allocate_condexp(arma::vec exp);
  void Condvar(arma::mat var);
  void condvar_zeros();
  
  // Export to R (test only)
  //Rcpp::List exportUpward2R() const;
};



/*
 * These are all the instances, numbered by the nodes.
 */

class Upward
{
private:
  
  Upward_Node * up; // Array if Upward_node elements
  unsigned int size; // Size
  
public:
  // Constructors
  Upward(int siz, int p_d);
  Upward(arma::mat const & data, int nE);
  
  // Destructor
  ~Upward();
  
  // Access to fields
  Upward_Node Up(int i) const;
  unsigned int Size() const;
  arma::vec Likelihood() const;
  
  // Recursion
  void recursion(Model const & mod, arma::umat const & ed);
  
  // Export to R (test only)
  Rcpp::List exportUpward2R(int i) const;
};


#endif
