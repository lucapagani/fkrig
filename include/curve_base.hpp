#ifndef H_FKRIG__
#define H_FKRIG__

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <GoTools/geometry/SplineCurve.h>
#include <GoTools/geometry/SplineInterpolator.h>
#include <GoTools/geometry/SplineApproximator.h>
#include <GoTools/geometry/GeometryTools.h>
#include <algorithm>
#include "covariance.hpp"

using std::vector;
using std::unique_ptr;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RVectorXd;

namespace fkrig {

class CurveBase {

public:

  //! Default constructor
  CurveBase () = default;

  //! Desctructor
  virtual ~CurveBase() {};

  /*! @brief Constructor with parameterization as input
   *
   *  Constructor of the curve with points and parameterization of each point as input
   *  
   *  @param points a vector with the vectors of the observed points
   *  @param param a vector with the vectors of the parameterization of the points
   *  @param dim the dimension of the points (only dim=1 supported)
   *  @param coord matrix of the coordinates of the design space
   *  @param model the type of the model matrix
   *  @param n_max number of maximum iterazions of the algorithm the estimate the functional kriging parameters
   *  @param tol tollerance of the diffenrence of the estimated parameters between two different steps of the algorithm
   */
  CurveBase ( vector< vector<double> >& points,
              vector< vector<double> >& param,
              int dim,
              MatrixXd coord,
              model_type model,
              double n_max,
              double tol );

  //! Interpolate the input data
  void
  Interpolate ();

  /*! @brief Smoothing the input data
   * 
   *  Smothing the input data with regression spline with equally space control points
   * 
   *  @param num_coefs number of coefficients of the regression spline (number of control points)
   */
  void
  Smoothing ( int num_coefs );

  //! Compute the functional krigin parameters
  virtual void
  Eval () = 0;

  /*! Set covariance type
   *  
   *  @param phi set manually the value of the parameter phi
   *  @param sigma_sq set manually the value of the parameter sigma square (the variance of the gaussian process)
   *  @param type set the type of covariance function
   *  @param fix_nugget true if the nugget must be fixed
   *  @param nugget set the value of the nugget if fix_nugget is true  
   */
  void
  set_covariance ( double phi,
                   double sigma_sq,
                   cov_type type = fkrig::MATERN_5_2,
                   bool fix_nugget = false,
                   double nugget = 0. );

  /*! Set covariance type
   *  
   *  @param type set the type of covariance function
   *  @param fix_nugget true if the nugget must be fixed
   *  @param nugget set the value of the nugget if fix_nugget is true  
   */
  void
  set_covariance ( cov_type type = fkrig::MATERN_5_2,
                   bool fix_nugget = false,
                   double nugget = 0. );

  /*! Set covariance matrix
   *  
   *  @param cov std::unique pointer to a covariance object
   */
  void
  set_covaraince ( std::unique_ptr<fkrig::Covariance> cov );

  /*! Predict the value at the design coordinates coord in the geometric value param
   *  @param coord matrix of coordinates of the new design locations
   *  @param param vector of parameter values for each new design location
   */
  virtual vector<vector<double> >
  Predict ( MatrixXd& coord,
            vector<double>& param ) const = 0;

  /*! Predict the value at the design coordinates coord in n_u equally spaced geomtric points
   *
   * @param coord matrix of coordinates of the new design locations
   * @param n_u number of points to predict in the geometric space
   */
  vector<vector<double> >
  Predict ( MatrixXd& coord,
            int n_u ) const;

  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord row vector of coordinates of the new desing locations
   */
  virtual Go::SplineCurve
  Predict ( RVectorXd& coord ) const = 0;            
            
  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord matrix of coordinates of the new desing locations
   */
  virtual vector<Go::SplineCurve>
  Predict ( MatrixXd& coord ) const = 0;

  /*! Predicted the total covarince matrix at the design coordinates coord
   * 
   *  @param coord matrix of coordinates of the new design locations
   */
  virtual MatrixXd 
  PredictCovariance ( MatrixXd& coord ) const = 0;  
  
  //! Return the points in the design space
  MatrixXd
  get_coord() const {
    return this->coord_;
  };
 
  //! Return the points in the geometric space
  vector< vector<double> >
  get_points() const {
    return this->points_;
  };

  //! Return the parameterization in the geometric space
  vector< vector<double> >
  get_param() const {
    return this->param_;
  };

  //! Return the functional dataset
  vector< shared_ptr < Go::SplineCurve> >
  get_curve() const {
    return this->curve_ptr_;
  };
  
  //! Return the model matrix
  MatrixXd
  get_model_matrix() const {
    return this->F_;
  };

  //! Return the number of iterations of the algorithm
  int
  get_n_iter () const {
    return this->n_iter_;
  };

  //! Return if the algortithm has converged
  bool
  get_convergence () const {
    return this->conv_;
  };

  //! Return the covariance object
  fkrig::Covariance
  get_covariance () const {
    return *cov_;
  };
  
  //! Retrun the range of the curves
  std::pair<double, double>
  get_range () const {
    return this->range_points_;
  };
 
  //! Return the range of the parametric space (max - min) 
  double
  get_domain_range () const {
    
    double value = 0.;
    
    value = range_points_.second - range_points_.first;
    
    return value;
  };  

  //! Return the distance between the points in the design space
  VectorXd
  get_dist_design () const {
    return this->dist_; 
  };
  
  //! Return the distance between the points in the geometric space
  VectorXd
  get_dist_geometric () const {
    return this->par_sq_dist_; 
  };  
  
protected:

  //! Compute the model matrix
  void
  MakeModelMatrix ();

  /*! Compute the model matrix based on the matrix coord
   * 
   *  @param coord matrix of coordinate in the design space
   */
  MatrixXd
  MakeModelMatrix ( MatrixXd coord ) const;

  //! Compute the matrix of coefficients of the term a with the Least Square method
  virtual void
  ComputeALs () = 0;

  //! Compute the matrix of coefficients of the term m with the Least Square method
  virtual void
  ComputeMLs () = 0;

  //! Compute the matrix of coefficients of the term a with the Generalized Least Square method
  virtual void
  ComputeAGls () = 0;

  //! Compute the matrix of coefficients of the term m with the Generalized Least Square method
  virtual void
  ComputeMGls () = 0;

  //! Compute the coefficients of the spline of the term m
  virtual void
  ComputeDiffCoefs () = 0;

  //! Make the same knots vector for the curves
  void
  MakeUniqueRange ();

  //! Copmuted the vector of the square of the parwise distances
  virtual void
  ComputeSqPairwiseDistances () = 0;

  //! Compute the parameters of the covariance matrix
  void
  ComputeCovPar ();
  
  //! Points in the geometric space
  vector< vector<double> > points_;
  //! Parameterization in the geometric space
  vector< vector<double> > param_;
  //! Dimension of the geometric space
  int dim_; 
  //! Coordinates of the design matrix
  MatrixXd coord_;
  //! Model for construct the model matrix 
  model_type model_;
  //! Shared pointer to the input curves (functional dataset)
  vector< shared_ptr< Go::SplineCurve> > curve_ptr_; 
  //! Object of class covariance
  unique_ptr<fkrig::Covariance> cov_; 
  //! Maximum number of iteration for algorithm of the parameters estimation
  double n_max_; 
  //! Tollerance value for the convergence of the algorithm of the parameters estimation
  double tol_;
  //! Range of the curves
  std::pair<double, double> range_points_; 
  //! Model matrix
  MatrixXd F_;
  //! Union of the knots of the curves
  vector<double> knots_union; 
  //! Vector of the parwise distances in the coordinate space [d(x_0, x_1), d(x_0, x_2), d(x_0, x_n-1), ..., d(x_n-2, x_n-1)]
  VectorXd dist_;
  //! Vector of the square of the parwise distances in the functional space [d(x_0, x_1), d(x_0, x_2), d(x_0, x_n-1), ..., d(x_n-2, x_n-1)]
  VectorXd par_sq_dist_; 
  //! Covariance matrix
  MatrixXd sigma_; 
   //! Number of iteration of the algorithm
  int n_iter_;
  //! true if the algorithm has converged
  bool conv_; 

};

} //! end of namespace

#endif
