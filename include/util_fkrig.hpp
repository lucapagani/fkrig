#ifndef H_UTIL_FKRIG__
#define H_UTIL_FKRIG__

#include "curve.hpp"
#include "surface.hpp"
#include "ego_base.hpp"
#include "ego_curve.hpp"
#include "ego_surface.hpp"
#include <boost/math/distributions/normal.hpp> // for normal_distribution

using boost::math::normal; // typedef provides default type is double.

namespace fkrig {

//! Recursive auxiliary function for adaptiveSimpsons() function below
double
AdaptiveSimpsonsAux ( double ( *f ) ( Go::SplineCurve, double ),
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
AdaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double ),
                   Go::SplineCurve curve,
                   double a, double b,
                   double epsilon,
                   int maxRecursionDepth );


//! Recursive auxiliary function for adaptiveSimpsons() function below
double
AdaptiveSimpsonsAux ( double ( *f ) ( Go::SplineCurve, double, double ),
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
AdaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double, double ),
                   Go::SplineCurve curve,
                   double sd,
                   double a, double b,
                   double epsilon,
                   int maxRecursionDepth );

//! Recursive auxiliary function for adaptiveSimpsons() function below
double
AdaptiveSimpsonsAux ( double ( *f ) ( vector< shared_ptr< Go::SplineCurve > >, double, vector<MatrixXd> ),
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
AdaptiveSimpsons ( double ( *f ) ( vector< shared_ptr< Go::SplineCurve > >, double, vector<MatrixXd> ),
                   vector< shared_ptr< Go::SplineCurve > > curve,
                   vector<MatrixXd> llt_sigma_folded,
                   double a, double b,
                   double epsilon,
                   int maxRecursionDepth );

//! Compute the squre of the curve in the parametric point param
double
SquareCurvePoint ( Go::SplineCurve curve,
                   double param );

//! Compute the absolute value of the curve in the parametric point param
double
AbsCurvePoint ( Go::SplineCurve curve,
                double param );

//! Compute the expected value of the absolute value of the curve in the parametric point param
double
EAbsCurvePoint ( Go::SplineCurve curve,
                 double param,
                 double sd );

//! Compute the variance of the absolute value of the curve in the parametric point param
double
VarAbsCurvePoint ( Go::SplineCurve curve,
                   double param,
                   double sd );

//! Compute the expected value for the expected improvment
double
MeanEiCurve ( vector< shared_ptr <Go::SplineCurve> > curve_ptr,
              double param,
              vector<MatrixXd> llt_sigma_folded );

//! Compute the expected value of the 2d folded normal distribution
Eigen::Vector2d
MeanCurve ( vector< shared_ptr <Go::SplineCurve> > curve_ptr,
            double param,
            vector<MatrixXd> llt_sigma_folded );

//! Compute the variance for the expected improvment
double
VarianceEiCurve ( vector< shared_ptr <Go::SplineCurve> > curve_ptr,
                  double param,
                  vector<MatrixXd> llt_sigma_folded );

//! Compute the square of the surface in the parametric point param
int
SquareSurfacePoint ( unsigned ndim,
                     const double *x,
                     void *surface,
                     unsigned fdim,
                     double *fval );

//! Compute the absolute value of the surface in the parametric point param
int
AbsSurfPoint ( unsigned ndim,
               const double *x,
               void *surface,
               unsigned fdim,
               double *fval );

//! Compute the expected value of the absolute value of the surface in the parametric point param
int
EAbsSurfPoint ( unsigned ndim,
                const double *x,
                void *util_surf,
                unsigned fdim,
                double *fval );

//! Compute the variance of the absolute value of the curve in the parametric point param
int
VarAbsSurfPoint ( unsigned ndim,
                  const double *x,
                  void *util_surf,
                  unsigned fdim,
                  double *fval );

//! Compute the expected value for the expected improvment
int
MeanEiSurf ( unsigned ndim,
             const double *x,
             void *util_surf,
             unsigned fdim,
             double *fval );

//! Compute the expected value of the 2d folded normal distribution
Eigen::Vector2d
MeanSurf ( vector< shared_ptr <Go::SplineSurface> > surf_ptr,
           const double *x,
           vector<MatrixXd> llt_sigma_folded );

//! Compute the variance for the expected improvment
int
VarianceEiSurf ( unsigned ndim,
                 const double *x,
                 void *util_surf,
                 unsigned fdim,
                 double *fval );

//! Compute the square of the surface in the parametric point param if there is a polygonal boundary
int
SquareSurfacePointPoly ( unsigned ndim,
                         const double *x,
                         void *util,
                         unsigned fdim,
                         double *fval );

//! Compute the absolute value of the surface in the parametric point param if there is a polygonal boundary
int
AbsSurfPointPoly ( unsigned ndim,
                   const double *x,
                   void *util,
                   unsigned fdim,
                   double *fval );

//! Compute the expected value of the absolute value of the surface in the parametric point param if there is a polygonal boundary
int
EAbsSurfPointPoly ( unsigned ndim,
                    const double *x,
                    void *util_surf,
                    unsigned fdim,
                    double *fval );

//! Compute the variance of the absolute value of the curve in the parametric point param
int
VarAbsSurfPointPoly ( unsigned ndim,
                      const double *x,
                      void *util_surf,
                      unsigned fdim,
                      double *fval );

//! Compute the expected value for the expected improvment if there is a polygonal boundary
int
MeanEiSurfPoly ( unsigned ndim,
                 const double *x,
                 void *util_surf,
                 unsigned fdim,
                 double *fval );

//! Compute the variance for the expected improvment if there is a polygonal boundary
int
VarianceEiSurfPoly ( unsigned ndim,
                     const double *x,
                     void *util_surf,
                     unsigned fdim,
                     double *fval );

/*! Objective function for the maximization of the expected improvment
 *
 *  @param x value of the coordinates in the disegn space
 *  @param grad value of the gradient
 *  @param param pointer to ego object
 */
double
ObjectiveFunction ( const vector<double> &x,
                    vector<double> &grad,
                    void* ego_obj );

/*! Objective function for the minization of the distance between the predicted function and the nominal function
 *
 *  @param x value of the coordinates in the disegn space
 *  @param grad value of the gradient
 *  @param param pointer to ego object
 */
double
ObjectiveFunctionMinCurve ( const vector<double> &x,
                            vector<double> &grad,
                            void* ego_obj );

/*! Objective function for the minization of the distance between the predicted function and the nominal function
 *
 *  @param x value of the coordinates in the disegn space
 *  @param grad value of the gradient
 *  @param param pointer to ego object
 */
double
ObjectiveFunctionMinSurf ( const vector<double> &x,
                           vector<double> &grad,
                           void* ego_obj );

/*! Check if a point is on the polygon (the point is inside the polygon if it is on the border)
 *
 *  @param u value of the u parameter
 *  @param v value of the v parameter
 *  @param polygon boundary polygon, the point are in anticlockwise order
 */
bool
PnPoly ( double u,
         double v,
         vector<Point> polygon );

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
