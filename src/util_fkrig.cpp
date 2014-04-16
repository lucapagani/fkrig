#include "util_fkrig.hpp"

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
                      int bottom )
{
  double c = ( a + b ) /2, h = b - a;
  double d = ( a + c ) /2, e = ( c + b ) /2;
  double fd = f ( curve, d ), fe = f ( curve, e );
  double Sleft = ( h/12 ) * ( fa + 4*fd + fc );
  double Sright = ( h/12 ) * ( fc + 4*fe + fb );
  double S2 = Sleft + Sright;
  if ( bottom <= 0 || fabs ( S2 - S ) <= 15*epsilon )
    return S2 + ( S2 - S ) /15;
  return adaptiveSimpsonsAux ( f, curve, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1 ) +
         adaptiveSimpsonsAux ( f, curve, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1 );
}


//! Adaptive Simpson's Rule
double
adaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double ), // ptr to function
                   Go::SplineCurve curve, // curve
                   double a, double b,  // interval [a,b]
                   double epsilon,  // error tolerance
                   int maxRecursionDepth )  // recursion cap
{
  double c = ( a + b ) /2, h = b - a;
  double fa = f ( curve, a ), fb = f ( curve, b ), fc = f ( curve, c );
  double S = ( h/6 ) * ( fa + 4*fc + fb );
  return adaptiveSimpsonsAux ( f, curve, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth );
}

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
                      int bottom )
{
  double c = ( a + b ) /2, h = b - a;
  double d = ( a + c ) /2, e = ( c + b ) /2;
  double fd = f ( curve, d, sd ), fe = f ( curve, e, sd );
  double Sleft = ( h/12 ) * ( fa + 4*fd + fc );
  double Sright = ( h/12 ) * ( fc + 4*fe + fb );
  double S2 = Sleft + Sright;
  if ( bottom <= 0 || fabs ( S2 - S ) <= 15*epsilon )
    return S2 + ( S2 - S ) /15;
  return adaptiveSimpsonsAux ( f, curve, sd, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1 ) +
         adaptiveSimpsonsAux ( f, curve, sd, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1 );

}

//! Adaptive Simpson's Rule
double
adaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double, double ), // ptr to function
                   Go::SplineCurve curve, // curve
                   double sd, // standard deviation
                   double a, double b,  // interval [a,b]
                   double epsilon,  // error tolerance
                   int maxRecursionDepth )  // recursion cap
{
  double c = ( a + b ) /2, h = b - a;
  double fa = f ( curve, a, sd ), fb = f ( curve, b, sd ), fc = f ( curve, c, sd );
  double S = ( h/6 ) * ( fa + 4*fc + fb );
  return adaptiveSimpsonsAux ( f, curve, sd, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth );
}

//! Compute the squre of the curve in the parametric point param
double
square_curve_point ( Go::SplineCurve curve,
                     double param )
{

  Go::Point point;

  curve.point ( point, param );

  return pow ( point[0], 2 );

}

//! Compute the squre of the curve in the parametric point param
double
abs_curve_point ( Go::SplineCurve curve,
                  double param )
{

  Go::Point point;

  curve.point ( point, param );

  return std::abs ( point[0] );

}

//! Compute the squre of the curve in the parametric point param
double
e_abs_curve_point ( Go::SplineCurve curve,
                    double param,
                    double sd )
{

  Go::Point point;

  curve.point ( point, param );
  double mean = point[0];
  double ratio = mean / sd;
  
  normal z;
  
  double value = sd * std::sqrt ( 2 ) / boost::math::constants::root_pi<double>() * std::exp ( - std::pow ( ratio, 2 ) / 2 ) + mean * ( 1 - 2 * cdf ( z, - ratio ) );
  
  return value;
}

//! Compute the squre of the surface in the parametric point param
int
square_surface_point ( unsigned ndim,
                       const double *x,
                       void *surface,
                       unsigned fdim,
                       double *fval )
{
  Go::SplineSurface * surf = static_cast<Go::SplineSurface*> ( surface );

  Go::Point point;

  surf->point ( point, x[0], x[1] );

  fval[0] = pow ( point[0], 2 );

  return 0; // success
}

//! Compute the absolute value of the surface in the parametric point param
int
abs_surface_point ( unsigned ndim,
                    const double *x,
                    void *surface,
                    unsigned fdim,
                    double *fval )
{
  Go::SplineSurface * surf = static_cast<Go::SplineSurface*> ( surface );

  Go::Point point;

  surf->point ( point, x[0], x[1] );

  fval[0] = std::abs ( point[0] );

  return 0; // success
}

//! Objective function for the maximization of the expected improvment
double
ObjectiveFunction ( const vector<double> &x,
                    vector<double> &grad,
                    void* ego_obj )
{

  fkrig::EgoBase* ego = reinterpret_cast<fkrig::EgoBase*> ( ego_obj );

//   MatrixXd coord_ego;

  MatrixXd coord_ego = ego->get_coord ();

  // Fill the row of the second point
  for ( size_t i = 0; i < x.size (); ++i )
    coord_ego ( 1,i ) = x[i];

  // Compute mean and variance of random variable
  double mean = ego->ComputeMean ( coord_ego );
  double sigma = std::sqrt ( ego->ComputeVariance ( coord_ego ) );
  double ratio = mean / sigma;

  // Compute the value of the expected improvment
  double value = sigma * boost::math::pdf ( ego->z_, ratio ) + mean * boost::math::cdf ( ego->z_, ratio );

  return value;
}

// struct covariance_residual_mater_5_2_fix_nugget {
//   covariance_residual_mater_5_2_fix_nugget(double x, double y)
//       : x_(x), y_(y) {}
//
//   template <typename T>
//   bool operator()(const T* const sigma_sq, const T* const phi, T* residual) const {
//     const T* const x_phi = x_ / phi;
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
