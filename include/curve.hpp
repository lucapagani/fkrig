#ifndef H_FKRIG_CURVE__
#define H_FKRIG_CURVE__

#include "curve_base.hpp"
// #include "ego_curve.hpp"

namespace fkrig {

class Curve: public CurveBase {
//   friend class EgoCurve;

public:

  //! Default constructor
  Curve () = default;

  //! Desctructor
  ~Curve() {};

  /*! Constructor with implicit parameterization (only if the dimension of the geometric space is 2)
   * 
   * @todo all
   */
  Curve ( vector< vector<double> >& points,
          const int dim,
          MatrixXd coord,
          model_type model = CONSTANT,
          double n_max = 1e3,
          double tol = 1e-6 );

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
  Curve ( vector< vector<double> >& points,
          vector< vector<double> >& param,
          const int dim,
          MatrixXd coord,
          model_type model = CONSTANT,
          double n_max = 1e3,
          double tol = 1e-6 );

//   //! Interpolate the input data
//   void
//   Interpolate ();
// 
//   //! Smoothing the input data
//   void
//   Smoothing ( int num_coefs );

  //! Compute the functional krigin parameters
  void
  Eval ();

  /*! Predict the value at the design coordinates coord in the geometric value param
   *  @param coord matrix of coordinates of the new design locations
   *  @param param vector of parameter values for each new design location
   */
  vector<vector<double> >
  Predict ( MatrixXd& coord,
            vector<double>& param ) const;

  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord row vector of coordinates of the new desing locations
   */
  Go::SplineCurve
  Predict ( RVectorXd& coord ) const;            
            
  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord matrix of coordinates of the new desing locations
   */
  vector<Go::SplineCurve>
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
