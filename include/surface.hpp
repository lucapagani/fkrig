#ifndef H_FKRIG_SURFACE__
#define H_FKRIG_SURFACE__

#include "surface_base.hpp"

namespace fkrig {

class Surf: public SurfBase {

public:

  //! Default constructor
  Surf () = default;

  //! Desctructor
  ~Surf() {};

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
  Surf ( vector< vector<double> >& points,
         vector< vector<double> >& param_u,
         vector< vector<double> >& param_v,
         const int dim,
         MatrixXd coord,
         model_type model = CONSTANT,
         double n_max = 1e3,
         double tol = 1e-6 );

  //! Compute the functional krigin parameters
  void
  Eval ();

  /*! Predict the value at the design coordinates coord in the geometric value param
   *  @param coord matrix of coordinates of the new design locations
   *  @param param_u vector of parameter values in the u direction for each new design location
   *  @param param_v vector of parameter values in the v direction for each new design location 
   */
  vector<vector<double> >
  Predict ( MatrixXd& coord,
            vector<double>& param_u,
            vector<double>& param_v ) const;

  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord matrix of coordinates of the new desing locations
   */
  vector<Go::SplineSurface>
  Predict ( MatrixXd& coord ) const;

  /*! Predicted the total covarince matrix at the design coordinates coord
   * 
   *  @param coord matrix of coordinates of the new design locations
   */
  MatrixXd 
  PredictCovariance ( MatrixXd& coord ) const;  
  
protected:

  //! Compute the matrix of coefficients of the term a with the Least Square method
  void
  ComputeALs ();

  //! Compute the matrix of coefficients of the term m with the Least Square method
  void
  ComputeMLs ();

  //! Compute the matrix of coefficients of the term a with the Generalized Least Square method
  void
  ComputeAGls ();

  //! Compute the matrix of coefficients of the term m with the Generalized Least Square method
  void
  ComputeMGls ();

  //! Compute the coefficients of the spline of the term m
  void
  ComputeDiffCoefs ();

  //! Copmute the vector of the square of the parwise distances
  void
  ComputeSqPairwiseDistances ();

private:

  //! Matrix of coefficients of the term a
  MatrixXd a_coefs_;
  //! Matrix of coefficients of the term m
  MatrixXd m_coefs_;
  //! Union of the knots of the curves
  vector<double> knots_union;
  //! Vector of coefficients of the splines of the term a
  MatrixXd a_curve_coefs_;
  //! Vector of coefficients of the difference between the curves and the splines of the term m
  MatrixXd m_diff_coefs_;
  //! Spline coefficients of inv(S) (z - F a)
  MatrixXd dx_coefs_;

};

} //! end of namespace

#endif
