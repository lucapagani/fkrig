#ifndef H_UTIL_FKRIG__
#define H_UTIL_FKRIG__

#include "curve.hpp"
#include "surface.hpp"
#include "ego_base.hpp"
#include <boost/math/distributions/normal.hpp> // for normal_distribution

using boost::math::normal; // typedef provides default type is double.

namespace fkrig {

//! Recursive auxiliary function for adaptiveSimpsons() function below
double
adaptiveSimpsonsAux ( double ( *f ) ( Go::SplineCurve, double ),
                      Go::SplineCurve curve,
                      double a,
                      double b,
                      double epsilon,
                      double S,
                      double fa,
                      double fb,
                      double fc,
                      int bottom );

//! Adaptive Simpson's Rule
double
adaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double ),
                   Go::SplineCurve curve,
                   double a, double b,
                   double epsilon,
                   int maxRecursionDepth );


//! Recursive auxiliary function for adaptiveSimpsons() function below
double
adaptiveSimpsonsAux ( double ( *f ) ( Go::SplineCurve, double, double ),
                      Go::SplineCurve curve,
                      double sd,
                      double a,
                      double b,
                      double epsilon,
                      double S,
                      double fa,
                      double fb,
                      double fc,
                      int bottom );

//! Adaptive Simpson's Rule
double
adaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double, double ),
                   Go::SplineCurve curve,
                   double sd,
                   double a, double b,
                   double epsilon,
                   int maxRecursionDepth );

//! Recursive auxiliary function for adaptiveSimpsons() function below
double
adaptiveSimpsonsAux ( double ( *f ) ( vector< shared_ptr< Go::SplineCurve > >, double, vector<MatrixXd> ),
                      vector< shared_ptr< Go::SplineCurve > > curve,
                      vector<MatrixXd> llt_sigma_folded,
                      double a,
                      double b,
                      double epsilon,
                      double S,
                      double fa,
                      double fb,
                      double fc,
                      int bottom );

//! Adaptive Simpson's Rule
double
adaptiveSimpsons ( double ( *f ) ( vector< shared_ptr< Go::SplineCurve > >, double, vector<MatrixXd> ),
                   vector< shared_ptr< Go::SplineCurve > > curve,
                   vector<MatrixXd> llt_sigma_folded,
                   double a, double b,
                   double epsilon,
                   int maxRecursionDepth );

//! Compute the squre of the curve in the parametric point param
double
square_curve_point ( Go::SplineCurve curve,
                     double param );

//! Compute the absolute value of the curve in the parametric point param
double
abs_curve_point ( Go::SplineCurve curve,
                  double param );

//! Compute the expected value of the absolute value of the curve in the parametric point param
double
e_abs_curve_point ( Go::SplineCurve curve,
                    double param,
                    double sd );

//! Compute the expected value for the expected improvment
double
MeanEiCurve ( vector< shared_ptr <Go::SplineCurve> > curve_ptr,
              double param,
              vector<MatrixXd> llt_sigma_folded );

//! Compute the square of the surface in the parametric point param
int
square_surface_point ( unsigned ndim,
                       const double *x,
                       void *surface,
                       unsigned fdim,
                       double *fval );

//! Compute the absolute value of the surface in the parametric point param
int
abs_surface_point ( unsigned ndim,
                    const double *x,
                    void *surface,
                    unsigned fdim,
                    double *fval );

/*! Objective function for the maximization of the expected improvment
 *
 *  @param x value of the coordinates in the disegn space
 *  @param grad value of the gradient
 *  @param param vector with value of the mean and standard deviazion of the random variable
 */
double
ObjectiveFunction ( const vector<double> &x,
                    vector<double> &grad,
                    void* ego_obj );

// // struct covariance_residual_mater_5_2_fix_nugget { double x_; double y_; };
// struct covariance_residual_mater_5_2_fix_nugget {
//   covariance_residual_mater_5_2_fix_nugget(double x, double y)
//       : x_(x), y_(y) {}
//
//   template <typename T>
//   bool operator()(const T* const sigma_sq, const T* const phi, T* residual) const {
//     const double x_phi = x_ / phi;
//     const T* const value = 0.;
//     x_phi > 0. ? value = ( 1. + sqrt ( 5. ) * x_phi + 5 / 3 * pow ( x_phi, 2 ) ) * exp ( - sqrt ( 5 ) * x_phi ) : value = 1.;
//     residual[0] = T(y_) - sigma_sq * ( 1 - value );
//     return true;
//   }
//
//  private:
//   // Observations for a sample.
//   const double x_;
//   const double y_;
// };

} // End of namespace
#endif
