#ifndef H_FKRIG_SURF_BASE__
#define H_FKRIG_SURF_BASE__

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/geometry/SplineInterpolator.h>
#include <GoTools/geometry/SplineApproximator.h>
#include <GoTools/geometry/GeometryTools.h>
#include <algorithm>
#include "covariance.hpp"
#include <mba/MBA.h>
#include <mba/UCButils.h>
#include <mba/PointAccessUtils.h>
#include <cubature/cubature.h>

using std::vector;
using std::unique_ptr;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

namespace fkrig {

struct Point {
  double u;
  double v;
};
  
class SurfBase {

public:

  //! Default constructor
  SurfBase () = default;

  //! Desctructor
  virtual ~SurfBase() {};

  /*! @brief Constructor with parameterization as input
  *
  *  Constructor of the surface with points and parameterization of each point as input
  *
  *  @param points a vector with the vectors of the observed points
  *  @param param_u a vector with the vectors of the parameterization of the points in the u direction
  *  @param param_v a vector with the vectors of the parameterization of the points in the v direction
  *  @param dim the dimension of the points (only dim=1 supported)
  *  @param coord matrix of the coordinates of the design space
  *  @param model the type of the model matrix
  *  @param n_max number of maximum iterazions of the algorithm the estimate the functional kriging parameters
  *  @param tol tollerance of the diffenrence of the estimated parameters between two different steps of the algorithm
  */
  SurfBase ( vector< vector<double> >& points,
             vector< vector<double> >& param_u,
             vector< vector<double> >& param_v,
             int dim,
             MatrixXd coord,
             model_type model,
             double n_max,
             double tol );

  /*! Approximate the input data with MBA algorithm
   *
   *  @param h the number of levels of the MBA algorithm
   *  @param m0 number of paths in u direction
   *  @param n0 number of paths in v direction
   *  @param smoothing_iterations number of smoothing iterations for each levels
   */
  void
  Approximate ( int h = 7,
                int m0 = 1,
                int n0 = 1,
                int smoothing_iterations = 0 );

  /*! @brief Approximate the input data with MBA algorithm
   *
   *  Approximate the input data with the MBA algorithm choosing the number of levels between the first value of h_range and the second value of h_range.
   *  The level choose is the level with minum Mean Square Error based on a 10% cross validation procedure
   *
   *  @param h_range std::pair with the minimum and the maximum level to try
   *  @param m0 number of paths in u direction
   *  @param n0 number of paths in v direction
   *  @param smoothing_iterations number of smoothing iterations for each levels
   */
  void
  Approximate ( std::pair<int, int> h_range,
                int m0 = 1,
                int n0 = 1,
                int smoothing_iterations = 0 );

  //! Smoothing the input data
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
   *  @param param_u vector of parameter values in the u direction for each new design location
   *  @param param_v vector of parameter values in the v direction for each new design location
   */
  virtual vector<vector<double> >
  Predict ( MatrixXd& coord,
            vector<double>& param_u,
            vector<double>& param_v ) const = 0;

  /*! Predict the value at the design coordinates coord in n_u equally spaced geomtric points
   *
   * @param coord matrix of coordinates of the new design locations
   * @param n_u number of points to predict in the geometric space in the u direction
   * @param n_v number of points to predict in the geometric space in the v direction
   */
  vector<vector<double> >
  Predict ( MatrixXd& coord,
            int n_u,
            int n_v ) const;

  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord matrix of coordinates of the new desing locations
   */
  virtual vector<Go::SplineSurface>
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

  //! Return the parameterization in the geometric space in the u direction
  vector< vector<double> >
  get_param_u() const {
    return this->param_u_;
  };

  //! Return the parameterization in the geometric space in the v direction
  vector< vector<double> >
  get_param_v() const {
    return this->param_v_;
  };

  //! Return the functional dataset
  vector< shared_ptr < Go::SplineSurface> >
  get_surface() const {
    return this->surf_ptr_;
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

  //! Return the area of the parametric space 
  double
  get_domain_range () const {
    
    double value = 0.;
    
    value = ( range_points_u_.second - range_points_u_.first ) * ( range_points_v_.second - range_points_v_.first );
    
    return value;
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

  /*! Compute the MBA approximation and return a shared pointer to a SplineSurface
   *
   *  @param points a vector with the vectors of the observed points
   *  @param param_u a vector with the vectors of the parameterization of the points in the u direction
   *  @param param_v a vector with the vectors of the parameterization of the points in the v direction
   *  @param h the number of levels of the MBA algorithm
   *  @param m0 number of paths in u direction
   *  @param n0 number of paths in v direction
   *  @param smoothing_iterations number of smoothing iterations for each levels  *
   */
  shared_ptr< Go::SplineSurface >
  ComputeMba ( vector<double> points,
               vector<double> param_u,
               vector<double> param_v,
               int h = 7,
               int m0 = 1,
               int n0 = 1,
               int smoothing_iterations = 0 );

  /*! Check if a point is on the polygon (the point is inside the polygon if it is on the border)
   *
   *  @param u value of the u parameter
   *  @param v value of the v parameter
   */
  bool
  Pnpoly ( double u,
           double v );

  //! Points in the geometric space
  vector< vector<double> > points_;
  //! Parameterization in the geometric space in u direction
  vector< vector<double> > param_u_;
  //! Parameterization in the geometric space in v direction
  vector< vector<double> > param_v_;
  //! Dimension of the geometric space
  int dim_;
  //! Coordinates of the design matrix
  MatrixXd coord_;
  //! Model for construct the model matrix
  model_type model_;
  //! Shared pointer to the input surfaces (functional dataset)
  vector< shared_ptr< Go::SplineSurface> > surf_ptr_;
  //! Object of class covariance
  unique_ptr<fkrig::Covariance> cov_;
  //! Maximum number of iteration for algorithm of the parameters estimation
  double n_max_;
  //! Tollerance value for the convergence of the algorithm of the parameters estimation
  double tol_;
  //! Range of the surfaces in the u direction
  std::pair<double, double> range_points_u_;
  //! Range of the surfaces in the v direction
  std::pair<double, double> range_points_v_;
  //! Model matrix
  MatrixXd F_;
  //! Union of the knots of the surfaces
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
  //! Vertex of poligons if the region is not a square, the point are in anticlockwise order
  vector< Point > polygon_;

};

} //! end of namespace

#endif
