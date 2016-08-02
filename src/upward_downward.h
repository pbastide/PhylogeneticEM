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
  arma::mat edge; // matrix of edges in post-order
  arma::mat exps; // matix of expectations at all the nodes
  arma::cube vars; // 3D array of variances at all the nodes
  arma::cube covars; // 3D array of covariances at all nodes/parents. NA at the root.
  
public:
  // Constructors
  Moments(int = 0, int = 0);
  Moments(arma::mat const & ed, int p_d);
  Moments(arma::mat const & data, arma::mat const & ed);
  
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
 * These have (ntaxa + nNodes) entries, and NA values at the root.
 */
class Model
{
private:
  
  arma::mat r; // matix of vectors offsets
  arma::cube q; // 3D array of transmission matrices
  arma::cube sigma; // 3D array of variances matrics
  
public:
  // Constructors
  // For a BM:
  Model(arma::mat const & Delta, arma::mat const & Variance,
        arma::vec const & edge_length);
  // For an OU:
  Model(arma::mat const & Beta, arma::mat const & Stationary_Var,
        arma::vec const & edge_length, arma::mat const & Alpha);
  
  // Access to fields
  arma::mat R() const;
  arma::cube Q() const;
  arma::cube Sigma() const;
};

//---------------------------------------------------------------------------//
// Class Upward -------------------------------------------------------------//
//---------------------------------------------------------------------------//
/*
 * The quantities needed at the upward phase.
 * Numbered by the nodes.
 */
class Upward
{
private:
  
  arma::mat edge; // matrix of edges in post-order (moving from last to first row for preorder)
  arma::vec cst; // vector of constants
  arma::mat condexp; // matrix of expectations
  arma::cube condvar; // 3D array of variances matrices
  
public:
  // Constructors
  Upward(arma::mat const & ed, int p_d);
  Upward(arma::mat const & data, arma::mat const & ed);
  
  // Access to fields
  arma::vec Cst() const;
  arma::mat Condexp() const;
  arma::cube Condvar() const;
};



#endif
