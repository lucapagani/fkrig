#ifndef H_COV_FKRIGG__
#define H_COV_FKRIGG__

#include <vector>
#include <eigen3/Eigen/Dense>
#include <boost/math/special_functions.hpp>
#include <ceres/ceres.h>
// #include "util_fkrig.hpp"

using std::vector;
using std::pow;
using std::sqrt;
using std::exp;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_ptr;

namespace fkrig {

//! Model to buil di model matrix starting from the design matrix
enum model_type { CONSTANT, //!< Only the constant is included in the model
                  LINEAR, //!< Constant + linear terms
                  QUADRATIC, //!< Constant + linear terms + cross terms + quadratic terms
                  PURE_QUDRATIC //!< Constant + linear terms + quadratic terms
                };

//! Type of covariance matrix available              
enum cov_type { MATERN_5_2, //!< Matern correlation function with nu=2.5
                MATERN_3_2, //!< Matern correlation function with nu=1.5
                POW_EXP, //!< Powered exponential correlation function
                EXP, //!< Exponential correlation function
                GAUSS //!< Gaussian correlation function
              };

class Covariance {

public:

  //! Default constructor
  Covariance () = default;

  /*! @brief Constructor
   * 
   *  Constructor with all the parameters as input
   * 
   *  @param phi set manually the value of the parameter phi
   *  @param sigma_sq set manually the value of the parameter sigma square (the variance of the gaussian process)
   *  @param n dimension of the correlation matrix, the number of the locations in the design space
   *  @param type set the type of covariance function
   *  @param fix_nugget true if the nugget must be fixed
   *  @param nugget set the value of the nugget if fix_nugget is true  
   */
  Covariance ( double phi,
               double sigma_sq,
               int n,
               cov_type type = fkrig::MATERN_5_2,
               bool fix_nugget = false,
               double nugget = 0. );

  /*! @brief Constructor
   * 
   *  Constructor with no parameters as input
   * 
   *  @param n dimension of the correlation matrix, the number of the locations in the design space
   *  @param type set the type of covariance function
   *  @param fix_nugget true if the nugget must be fixed
   *  @param nugget set the value of the nugget if fix_nugget is true  
   */
  Covariance ( int n,
               cov_type type = fkrig::MATERN_5_2,
               bool fix_nugget = false,
               double nugget = 0. );

  // Destructor
  ~Covariance () {};

  /*! Compute the correlation vector
   * 
   *  @param u vector of distances between locations
   */
  void
  ComputeCorr ( VectorXd& u );

  /*! Compute the correalation at distance u
   * 
   *  @param u vector of distances between locations
   */
  double
  ComputeCorr ( const double& u ) const;

  /*! Compute the correlation matrix ad distance u and rearrange the value in a row-based matrix n_row x n_col
   * 
   *  @param u vector of distances between locations
   *  @param n_row number of rows of the correlation matrix
   *  @param n_col number of coloumns of the correlation matrix
   */
  MatrixXd
  ComputeCorr ( VectorXd& u,
                int n_row,
                int n_col ) const;

  /*! Compute the covariance matrix ad distance u and rearrange the value in a row-based matrix n_row x n_col
   * 
   *  @param u vector of distances between locations
   *  @param n_row number of rows of the correlation matrix
   *  @param n_col number of coloumns of the correlation matrix
   */
  MatrixXd
  ComputeCov ( VectorXd& u,
               int n_row,
               int n_col ) const;

  /*! Compute the covariance matrix
   * 
   *  @param u vector of distances between locations
   */
  void
  Eval ( VectorXd& u );

  /*! Estimate the parameter of the covariance function
   * 
   *  @param u vector of distances between locations
   *  @param v vector of distances between functions
   */
  void
  Estimate ( VectorXd& u,
             VectorXd& v );

  //! Set the phi parameter
  void
  set_phi ( double phi ) {
    phi_ = phi;
  };

  //! Return the correlation matrix
  MatrixXd
  get_corr () {
    return corr_matr_;
  };

  //! Return the covariance matrix
  MatrixXd
  get_cov () {
    return cov_matr_;
  };

  //! Return the phi parameter
  double
  get_phi () {
    return phi_;
  };

  //! Return the kappa parameter
  double
  get_kappa () {
    return kappa_;
  };

  //! Return the nugget value
  double
  get_nugget () {
    return nugget_;
  };

  //! Return the value of sigma_sq
  double
  get_sigma_sq () {
    return sigma_sq_;
  };

  //! Return the summary of the non linear optimization
  ceres::Solver::Summary
  get_summary () {
    return summary_;
  };

private:

  //! Value of the parameter phi of the Correlation function
  double phi_;
  //! Value of the covariance
  double sigma_sq_;
  //! Dimension of the square matrix
  int n_;
  //! Type of covariance
  cov_type type_;
  //! TRUE if the nugget must be fixed
  bool fix_nugget_;
  //! Value of the nugget
  double nugget_;
  //! Value of the parameter kappa
  double kappa_;
  //! Correlation vector
  VectorXd corr_vec_;
  //! Correlation matrix
  MatrixXd corr_matr_;
  //! Covariance matrix
  MatrixXd cov_matr_;
  //! Summary of the optimization
  ceres::Solver::Summary summary_; 

};

// template<typename C_type>
// void
// fkrig::covariance::compute_corr ( VectorXd u )
// {
//   std::cerr << "Non supported covariance function" << std::endl;
//   exit ( 1 );
// }

} // End of namespace

#endif
