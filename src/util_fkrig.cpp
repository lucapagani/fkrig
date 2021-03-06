#include "util_fkrig.hpp"

namespace fkrig {

//! Recursive auxiliary function for AdaptiveSimpsons() function below
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
  return AdaptiveSimpsonsAux ( f, curve, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1 ) +
         AdaptiveSimpsonsAux ( f, curve, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1 );
}


//! Adaptive Simpson's Rule
double
AdaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double ), // ptr to function
                   Go::SplineCurve curve, // curve
                   double a, double b,  // interval [a,b]
                   double epsilon,  // error tolerance
                   int maxRecursionDepth )  // recursion cap
{
  double c = ( a + b ) /2, h = b - a;
  double fa = f ( curve, a ), fb = f ( curve, b ), fc = f ( curve, c );
  double S = ( h/6 ) * ( fa + 4*fc + fb );
  return AdaptiveSimpsonsAux ( f, curve, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth );
}

//! Recursive auxiliary function for AdaptiveSimpsons() function below
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
  return AdaptiveSimpsonsAux ( f, curve, sd, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1 ) +
         AdaptiveSimpsonsAux ( f, curve, sd, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1 );

}

//! Adaptive Simpson's Rule
double
AdaptiveSimpsons ( double ( *f ) ( Go::SplineCurve, double, double ), // ptr to function
                   Go::SplineCurve curve, // curve
                   double sd, // standard deviation
                   double a, double b,  // interval [a,b]
                   double epsilon,  // error tolerance
                   int maxRecursionDepth )  // recursion cap
{
  double c = ( a + b ) /2, h = b - a;
  double fa = f ( curve, a, sd ), fb = f ( curve, b, sd ), fc = f ( curve, c, sd );
  double S = ( h/6 ) * ( fa + 4*fc + fb );
  return AdaptiveSimpsonsAux ( f, curve, sd, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth );
}

//! Recursive auxiliary function for AdaptiveSimpsons() function below
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
                      int bottom )
{
  double c = ( a + b ) /2, h = b - a;
  double d = ( a + c ) /2, e = ( c + b ) /2;
  double fd = f ( curve, d, llt_sigma_folded ), fe = f ( curve, e, llt_sigma_folded );
  double Sleft = ( h/12 ) * ( fa + 4*fd + fc );
  double Sright = ( h/12 ) * ( fc + 4*fe + fb );
  double S2 = Sleft + Sright;
  if ( bottom <= 0 || fabs ( S2 - S ) <= 15*epsilon )
    return S2 + ( S2 - S ) /15;
  return AdaptiveSimpsonsAux ( f, curve, llt_sigma_folded, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1 ) +
         AdaptiveSimpsonsAux ( f, curve, llt_sigma_folded, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1 );

}

//! Adaptive Simpson's Rule
double
AdaptiveSimpsons ( double ( *f ) ( vector< shared_ptr< Go::SplineCurve > >, double, vector<MatrixXd> ), // ptr to function
                   vector< shared_ptr< Go::SplineCurve > > curve, // curve
                   vector<MatrixXd> llt_sigma_folded, // standard deviation
                   double a, double b,  // interval [a,b]
                   double epsilon,  // error tolerance
                   int maxRecursionDepth )  // recursion cap
{
  double c = ( a + b ) /2, h = b - a;
  double fa = f ( curve, a, llt_sigma_folded ), fb = f ( curve, b, llt_sigma_folded ), fc = f ( curve, c, llt_sigma_folded );
  double S = ( h/6 ) * ( fa + 4*fc + fb );
  return AdaptiveSimpsonsAux ( f, curve, llt_sigma_folded, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth );
}

//! Compute the squre of the curve in the parametric point param
double
SquareCurvePoint ( Go::SplineCurve curve,
                   double param )
{

  Go::Point point;

  curve.point ( point, param );

  return pow ( point[0], 2 );

}

//! Compute the squre of the curve in the parametric point param
double
AbsCurvePoint ( Go::SplineCurve curve,
                double param )
{

  Go::Point point;

  curve.point ( point, param );

  return std::abs ( point[0] );

}

//! Compute the expected value of the absolute value of the curve in the parametric point param
double
EAbsCurvePoint ( Go::SplineCurve curve,
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

//! Compute the variance of the absolute value of the curve in the parametric point param
double
VarAbsCurvePoint ( Go::SplineCurve curve,
                   double param,
                   double sd )
{

  Go::Point point;

  curve.point ( point, param );
  double mean = point[0];
  double ratio = mean / sd;

  normal z;

  double value = std::pow ( mean, 2 ) + std::pow ( sd, 2 ) - std::pow ( sd * std::sqrt ( 2 ) / boost::math::constants::root_pi<double>() * std::exp ( - std::pow ( ratio, 2 ) / 2 ) + mean * ( 1 - 2 * cdf ( z, - ratio ) ), 2 );

  return value;
}

//! Compute the expected value for the expected improvment
double
MeanEiCurve ( vector< shared_ptr <Go::SplineCurve> > curve_ptr,
              double param,
              vector<MatrixXd> llt_sigma_folded )
{

  // Compute the expected value of the 2d folden normal distribution
  Eigen::Vector2d mu = MeanCurve ( curve_ptr, param, llt_sigma_folded );

//   std::cout << mu << "\n";

  // Compute mu_0 - mu_1
  double value = mu ( 0 ) - mu ( 1 );

  return value;
}

//! Compute the expected value of the 2d folded normal distribution
Eigen::Vector2d
MeanCurve ( vector< shared_ptr <Go::SplineCurve> > curve_ptr,
            double param,
            vector<MatrixXd> llt_sigma_folded )
{
  // Compute the inverse of llt_sigma_folded
  MatrixXd inv_llt_sigma_folded_0 = llt_sigma_folded[0].inverse();
  MatrixXd inv_llt_sigma_folded_1 = llt_sigma_folded[1].inverse();

  // Compute the mean on the param value
  Go::Point point_0, point_1;
  curve_ptr[0]->point ( point_0, param );
  curve_ptr[1]->point ( point_1, param );

  // Store values in a vector
  Eigen::Vector2d mu ( 2 );
  mu << point_0[0], point_1[0];

  normal z;

  // Compute values of inv_llt_sigma_folded * mu
  // k1 = k2 = 1
  Eigen::Vector2d prod_0 = inv_llt_sigma_folded_0 * mu;
  // k1 = 1, k2 = -1
  mu ( 1 ) = - mu ( 1 );
  Eigen::Vector2d prod_1 = inv_llt_sigma_folded_1 * mu;
  // k1 = k2 = -1
  mu ( 0 ) = - mu ( 0 );
  Eigen::Vector2d prod_3 = inv_llt_sigma_folded_0 * mu;
  // k1 = -1, k2 = 1
  mu ( 1 ) = - mu ( 1 );
  Eigen::Vector2d prod_2 = inv_llt_sigma_folded_1 * mu;

  // Compute mean of the first entry of the vector mu
  // k1 = k2 = 1
  double mu_0_0 = point_0[0] * cdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) ) + llt_sigma_folded[0] ( 0,0 ) * pdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) );
  // k1 = 1, k2 = -1
  double mu_0_1 = point_0[0] * cdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) ) + llt_sigma_folded[1] ( 0,0 ) * pdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) );
  // k1 = -1, k2 = 1
  double mu_0_2 = - point_0[0] * cdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) ) + llt_sigma_folded[1] ( 0,0 ) * pdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) );
  // k1 = -1, k2 = -1
  double mu_0_3 = - point_0[0] * cdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) ) + llt_sigma_folded[0] ( 0,0 ) * pdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) );

  // k1 = k2 = 1
  double mu_1_0 = point_1[0] * cdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) ) + llt_sigma_folded[0] ( 1,0 ) * pdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) ) + llt_sigma_folded[0] ( 1,1 ) * pdf ( z, prod_0 ( 1 ) ) * cdf ( z, prod_0 ( 0 ) );
  // k1 = 1, k2 = -1
  double mu_1_1 = - point_1[0] * cdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) ) + llt_sigma_folded[1] ( 1,0 ) * pdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) ) + llt_sigma_folded[1] ( 1,1 ) * pdf ( z, prod_1 ( 1 ) ) * cdf ( z, prod_1 ( 0 ) );
  // k1 = -1, k2 = 1
  double mu_1_2 = point_1[0] * cdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) ) + llt_sigma_folded[1] ( 1,0 ) * pdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) ) + llt_sigma_folded[1] ( 1,1 ) * pdf ( z, prod_2 ( 1 ) ) * cdf ( z, prod_2 ( 0 ) );
  // k1 = -1, k2 = -1
  double mu_1_3 = - point_1[0] * cdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) ) + llt_sigma_folded[0] ( 1,0 ) * pdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) ) + llt_sigma_folded[0] ( 1,1 ) * pdf ( z, prod_3 ( 1 ) ) * cdf ( z, prod_3 ( 0 ) );

  // Compute mu
  Eigen::Vector2d mu_delta;
  mu_delta << mu_0_0 + mu_0_1 + mu_0_2 + mu_0_3, mu_1_0 + mu_1_1 + mu_1_2 + mu_1_3;

  return mu_delta;
}

//! Compute the variance for the expected improvment
double
VarianceEiCurve ( vector< shared_ptr <Go::SplineCurve> > curve_ptr,
                  double param,
                  vector<MatrixXd> llt_sigma_folded )
{
  // Compute the inverse of llt_sigma_folded
  MatrixXd inv_llt_sigma_folded_0 = llt_sigma_folded[0].inverse();
  MatrixXd inv_llt_sigma_folded_1 = llt_sigma_folded[1].inverse();

  // Compute the mean on the param value
  Go::Point point_0, point_1;
  curve_ptr[0]->point ( point_0, param );
  curve_ptr[1]->point ( point_1, param );

  // Store values in a vector
  Eigen::Vector2d mu_0 ( 2 );
  mu_0 << point_0[0], point_1[0];

  // Compute the expected value of the 2d folden normal distribution
  Eigen::Vector2d mu = MeanCurve ( curve_ptr, param, llt_sigma_folded );

//   std::cout << mu << "\n";

  normal z;

  // Compute values of inv_llt_sigma_folded * mu
  // k1 = k2 = 1
  Eigen::Vector2d prod_0 = inv_llt_sigma_folded_0 * mu_0;
  // k1 = 1, k2 = -1
  Eigen::Vector2d mu_1;
  mu_1 << mu_0 ( 0 ),  - mu_0 ( 1 );
  Eigen::Vector2d prod_1 = inv_llt_sigma_folded_1 * mu_1;
  // k1 = -1, k2 = 1
  Eigen::Vector2d mu_2 = - mu_1;
  Eigen::Vector2d prod_2 = inv_llt_sigma_folded_1 * mu_2;
  // k1 = k2 = -1
  Eigen::Vector2d mu_3 = - mu_0;
  Eigen::Vector2d prod_3 = inv_llt_sigma_folded_0 * mu_3;

  // Compute the vector of pdf and cdf
  Eigen::Vector2d pdf_0;
  pdf_0 << pdf ( z, prod_0 ( 0 ) ), pdf ( z, prod_0 ( 1 ) );
  Eigen::Vector2d pdf_1;
  pdf_1 << pdf ( z, prod_1 ( 0 ) ), pdf ( z, prod_1 ( 1 ) );
  Eigen::Vector2d pdf_2;
  pdf_2 << pdf ( z, prod_2 ( 0 ) ), pdf ( z, prod_2 ( 1 ) );
  Eigen::Vector2d pdf_3;
  pdf_3 << pdf ( z, prod_3 ( 0 ) ), pdf ( z, prod_3 ( 1 ) );
  Eigen::Vector2d cdf_0;
  cdf_0 << cdf ( z, prod_0 ( 0 ) ), cdf ( z, prod_0 ( 1 ) );
  Eigen::Vector2d cdf_1;
  cdf_1 << cdf ( z, prod_1 ( 0 ) ), cdf ( z, prod_1 ( 1 ) );
  Eigen::Vector2d cdf_2;
  cdf_2 << cdf ( z, prod_2 ( 0 ) ), cdf ( z, prod_2 ( 1 ) );
  Eigen::Vector2d cdf_3;
  cdf_3 << cdf ( z, prod_3 ( 0 ) ), cdf ( z, prod_3 ( 1 ) );

  // Compute the product between the cdf and pdf
  Eigen::Vector2d pf_0;
  pf_0 << pdf_0 ( 0 ) * cdf_0 ( 1 ), pdf_0 ( 1 ) * cdf_0 ( 0 );
  Eigen::Vector2d pf_1;
  pf_1 << pdf_1 ( 0 ) * cdf_1 ( 1 ), pdf_1 ( 1 ) * cdf_1 ( 0 );
  Eigen::Vector2d pf_2;
  pf_2 << pdf_2 ( 0 ) * cdf_2 ( 1 ), pdf_2 ( 1 ) * cdf_2 ( 0 );
  Eigen::Vector2d pf_3;
  pf_3 << pdf_3 ( 0 ) * cdf_3 ( 1 ), pdf_3 ( 1 ) * cdf_3 ( 0 );

  // Compute the terms I_k_i
  Eigen::Matrix2d i_0;
  i_0 << cdf_0 ( 1 ) * ( cdf_0 ( 0 ) - prod_0 ( 0 ) * pdf_0 ( 0 ) ), pdf_0 ( 0 ) * pdf_0 ( 1 ),
      pdf_0 ( 0 ) * pdf_0 ( 1 ), cdf_0 ( 0 ) * ( cdf_0 ( 1 ) - prod_0 ( 1 ) * pdf_0 ( 1 ) );
  Eigen::Matrix2d i_1;
  i_1 << cdf_1 ( 1 ) * ( cdf_1 ( 0 ) - prod_1 ( 0 ) * pdf_1 ( 0 ) ), pdf_1 ( 0 ) * pdf_1 ( 1 ),
      pdf_1 ( 0 ) * pdf_1 ( 1 ), cdf_1 ( 0 ) * ( cdf_1 ( 1 ) - prod_1 ( 1 ) * pdf_1 ( 1 ) );
  Eigen::Matrix2d i_2;
  i_2 << cdf_2 ( 1 ) * ( cdf_2 ( 0 ) - prod_2 ( 0 ) * pdf_2 ( 0 ) ), pdf_2 ( 0 ) * pdf_2 ( 1 ),
      pdf_2 ( 0 ) * pdf_2 ( 1 ), cdf_2 ( 0 ) * ( cdf_2 ( 1 ) - prod_2 ( 1 ) * pdf_2 ( 1 ) );
  Eigen::Matrix2d i_3;
  i_3 << cdf_3 ( 1 ) * ( cdf_3 ( 0 ) - prod_3 ( 0 ) * pdf_3 ( 0 ) ), pdf_3 ( 0 ) * pdf_3 ( 1 ),
      pdf_3 ( 0 ) * pdf_3 ( 1 ), cdf_3 ( 0 ) * ( cdf_3 ( 1 ) - prod_3 ( 1 ) * pdf_3 ( 1 ) );

  // Compute the first addend of the covariance matrix
  // k1 = k2 = 1
  Eigen::Matrix2d add_0_0 = mu_0 * mu_0.transpose() * cdf_0 ( 0 ) * cdf_0 ( 1 );
  // k1 = 1, k2 = -1
  Eigen::Matrix2d add_0_1 = mu_1 * mu_1.transpose() * cdf_1 ( 0 ) * cdf_1 ( 1 );
  // k1 = -1, k2 = 1
  Eigen::Matrix2d add_0_2 = mu_2 * mu_2.transpose() * cdf_2 ( 0 ) * cdf_2 ( 1 );
  // k1 = k2 = -1
  Eigen::Matrix2d add_0_3 = mu_3 * mu_3.transpose() * cdf_3 ( 0 ) * cdf_3 ( 1 );

  // Compute the second addend of the covariance matrix
  // k1 = k2 = 1
  Eigen::Matrix2d add_1_0 = llt_sigma_folded[0] * pf_0 * mu_0.transpose();
  // k1 = 1, k2 = -1
  Eigen::Matrix2d add_1_1 = llt_sigma_folded[1] * pf_1 * mu_1.transpose();
  // k1 = -1, k2 = 1
  Eigen::Matrix2d add_1_2 = llt_sigma_folded[1] * pf_2 * mu_2.transpose();
  // k1 = k2 = -1
  Eigen::Matrix2d add_1_3 = llt_sigma_folded[0] * pf_3 * mu_3.transpose();

  // Compute the forth addend of the covariance matrix
  // k1 = k2 = 1
  Eigen::Matrix2d add_3_0 = llt_sigma_folded[0] * i_0 * llt_sigma_folded[0].transpose();
  // k1 = 1, k2 = -1
  Eigen::Matrix2d add_3_1 = llt_sigma_folded[1] * i_1 * llt_sigma_folded[1].transpose();
  // k1 = -1, k2 = 1
  Eigen::Matrix2d add_3_2 = llt_sigma_folded[1] * i_2 * llt_sigma_folded[1].transpose();
  // k1 = k2 = -1
  Eigen::Matrix2d add_3_3 = llt_sigma_folded[0] * i_3 * llt_sigma_folded[0].transpose();

  // Compute the covariance matrix
  Eigen::Matrix2d cov = add_0_0 + add_1_0 + add_1_0.transpose() + add_3_0 +
                        add_0_1 + add_1_1 + add_1_1.transpose() + add_3_1 +
                        add_0_2 + add_1_2 + add_1_2.transpose() + add_3_2 +
                        add_0_3 + add_1_3 + add_1_3.transpose() + add_3_3 -
                        mu * mu.transpose();

  // Compute the variance of the expected improvment
  double variance = cov ( 0, 0 ) + cov ( 1,1 ) - cov ( 0,1 ) - cov ( 1,0 );

  return variance;

}

//! Compute the squre of the surface in the parametric point param
int
SquareSurfacePoint ( unsigned ndim,
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
AbsSurfPoint ( unsigned ndim,
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

//! Compute the expected value of the absolute value of the curve in the parametric point param
int
EAbsSurfPoint ( unsigned ndim,
                const double *x,
                void *util_surf,
                unsigned fdim,
                double *fval )
{

  std::pair< shared_ptr<Go::SplineSurface>, double>* util = static_cast< std::pair< shared_ptr<Go::SplineSurface>, double>* > ( util_surf );

  Go::Point point;

  util->first->point ( point, x[0], x[1] );

  double mean = point[0];
  double ratio = mean / util->second;

  normal z;

  fval[0] = util->second * std::sqrt ( 2 ) / boost::math::constants::root_pi<double>() * std::exp ( - std::pow ( ratio, 2 ) / 2 ) + mean * ( 1 - 2 * cdf ( z, - ratio ) );

  return 0;
}

//! Compute the variance of the absolute value of the curve in the parametric point param
int
VarAbsSurfPoint ( unsigned ndim,
                  const double *x,
                  void *util_surf,
                  unsigned fdim,
                  double *fval )
{
  std::pair< shared_ptr<Go::SplineSurface>, double>* util = static_cast< std::pair< shared_ptr<Go::SplineSurface>, double>* > ( util_surf );

  Go::Point point;

  util->first->point ( point, x[0], x[1] );
  double mean = point[0];
  double ratio = mean / util->second;

  normal z;

  fval[0] = std::pow ( mean, 2 ) + std::pow ( util->second, 2 ) - std::pow ( util->second * std::sqrt ( 2 ) / boost::math::constants::root_pi<double>() * std::exp ( - std::pow ( ratio, 2 ) / 2 ) + mean * ( 1 - 2 * cdf ( z, - ratio ) ), 2 );

  return 0;

}

//! Compute the expected value for the expected improvment
int
MeanEiSurf ( unsigned ndim,
             const double *x,
             void *util_surf,
             unsigned fdim,
             double *fval )
{
  std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> >* util = static_cast< std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> >* > ( util_surf );

  vector< shared_ptr<Go::SplineSurface> > surf_ptr;
  surf_ptr.resize ( 2 );
  ( surf_ptr[0] ).reset ( new Go::SplineSurface );
  ( surf_ptr[1] ).reset ( new Go::SplineSurface );
  surf_ptr[0] = util->first[0];
  surf_ptr[1] = util->first[1];

  Eigen::LLT<MatrixXd> llt;
  vector<MatrixXd> llt_sigma_folded;
  llt_sigma_folded.resize ( 2 );
  llt_sigma_folded[0] = util->second[0];
  llt_sigma_folded[1] = util->second[1];

  // Compute the expected value of the 2d folden normal distribution
  Eigen::Vector2d mu = MeanSurf ( surf_ptr, x, llt_sigma_folded );

//   std::cout << mu << "\n";

  // Compute mu_0 - mu_1
  fval[0] = mu ( 0 ) - mu ( 1 );

  return 0;
}

//! Compute the expected value of the 2d folded normal distribution
Eigen::Vector2d
MeanSurf ( vector< shared_ptr <Go::SplineSurface> > surf_ptr,
           const double *x,
           vector<MatrixXd> llt_sigma_folded )
{
  // Compute the inverse of llt_sigma_folded
  MatrixXd inv_llt_sigma_folded_0 = llt_sigma_folded[0].inverse();
  MatrixXd inv_llt_sigma_folded_1 = llt_sigma_folded[1].inverse();

  // Compute the mean on the param value
  Go::Point point_0, point_1;
  surf_ptr[0]->point ( point_0, x[0], x[1] );
  surf_ptr[1]->point ( point_1, x[0], x[1] );

  // Store values in a vector
  Eigen::Vector2d mu ( 2 );
  mu << point_0[0], point_1[0];

  normal z;

  // Compute values of inv_llt_sigma_folded * mu
  // k1 = k2 = 1
  Eigen::Vector2d prod_0 = inv_llt_sigma_folded_0 * mu;
  // k1 = 1, k2 = -1
  mu ( 1 ) = - mu ( 1 );
  Eigen::Vector2d prod_1 = inv_llt_sigma_folded_1 * mu;
  // k1 = k2 = -1
  mu ( 0 ) = - mu ( 0 );
  Eigen::Vector2d prod_3 = inv_llt_sigma_folded_0 * mu;
  // k1 = -1, k2 = 1
  mu ( 1 ) = - mu ( 1 );
  Eigen::Vector2d prod_2 = inv_llt_sigma_folded_1 * mu;

  // Compute mean of the first entry of the vector mu
  // k1 = k2 = 1
  double mu_0_0 = point_0[0] * cdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) ) + llt_sigma_folded[0] ( 0,0 ) * pdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) );
  // k1 = 1, k2 = -1
  double mu_0_1 = point_0[0] * cdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) ) + llt_sigma_folded[1] ( 0,0 ) * pdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) );
  // k1 = -1, k2 = 1
  double mu_0_2 = - point_0[0] * cdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) ) + llt_sigma_folded[1] ( 0,0 ) * pdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) );
  // k1 = -1, k2 = -1
  double mu_0_3 = - point_0[0] * cdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) ) + llt_sigma_folded[0] ( 0,0 ) * pdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) );

  // k1 = k2 = 1
  double mu_1_0 = point_1[0] * cdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) ) + llt_sigma_folded[0] ( 1,0 ) * pdf ( z, prod_0 ( 0 ) ) * cdf ( z, prod_0 ( 1 ) ) + llt_sigma_folded[0] ( 1,1 ) * pdf ( z, prod_0 ( 1 ) ) * cdf ( z, prod_0 ( 0 ) );
  // k1 = 1, k2 = -1
  double mu_1_1 = - point_1[0] * cdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) ) + llt_sigma_folded[1] ( 1,0 ) * pdf ( z, prod_1 ( 0 ) ) * cdf ( z, prod_1 ( 1 ) ) + llt_sigma_folded[1] ( 1,1 ) * pdf ( z, prod_1 ( 1 ) ) * cdf ( z, prod_1 ( 0 ) );
  // k1 = -1, k2 = 1
  double mu_1_2 = point_1[0] * cdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) ) + llt_sigma_folded[1] ( 1,0 ) * pdf ( z, prod_2 ( 0 ) ) * cdf ( z, prod_2 ( 1 ) ) + llt_sigma_folded[1] ( 1,1 ) * pdf ( z, prod_2 ( 1 ) ) * cdf ( z, prod_2 ( 0 ) );
  // k1 = -1, k2 = -1
  double mu_1_3 = - point_1[0] * cdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) ) + llt_sigma_folded[0] ( 1,0 ) * pdf ( z, prod_3 ( 0 ) ) * cdf ( z, prod_3 ( 1 ) ) + llt_sigma_folded[0] ( 1,1 ) * pdf ( z, prod_3 ( 1 ) ) * cdf ( z, prod_3 ( 0 ) );

  // Compute mu
  Eigen::Vector2d mu_delta;
  mu_delta << mu_0_0 + mu_0_1 + mu_0_2 + mu_0_3, mu_1_0 + mu_1_1 + mu_1_2 + mu_1_3;

  return mu_delta;
}

//! Compute the variance for the expected improvment
int
VarianceEiSurf ( unsigned ndim,
                 const double *x,
                 void *util_surf,
                 unsigned fdim,
                 double *fval )
{
  std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> >* util = static_cast< std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> >* > ( util_surf );

  vector< shared_ptr<Go::SplineSurface> > surf_ptr;
  surf_ptr.resize ( 2 );
  ( surf_ptr[0] ).reset ( new Go::SplineSurface );
  ( surf_ptr[1] ).reset ( new Go::SplineSurface );
  surf_ptr[0] = util->first[0];
  surf_ptr[1] = util->first[1];

  Eigen::LLT<MatrixXd> llt;
  vector<MatrixXd> llt_sigma_folded;
  llt_sigma_folded.resize ( 2 );
  llt_sigma_folded[0] = util->second[0];
  llt_sigma_folded[1] = util->second[1];

  // Compute the inverse of llt_sigma_folded
  MatrixXd inv_llt_sigma_folded_0 = llt_sigma_folded[0].inverse();
  MatrixXd inv_llt_sigma_folded_1 = llt_sigma_folded[1].inverse();

  // Compute the mean on the param value
  Go::Point point_0, point_1;
  surf_ptr[0]->point ( point_0, x[0], x[1] );
  surf_ptr[1]->point ( point_1, x[0], x[1] );

  // Store values in a vector
  Eigen::Vector2d mu_0 ( 2 );
  mu_0 << point_0[0], point_1[0];

  // Compute the expected value of the 2d folden normal distribution
  Eigen::Vector2d mu = MeanSurf ( surf_ptr, x, llt_sigma_folded );

//   std::cout << mu << "\n";

  normal z;

  // Compute values of inv_llt_sigma_folded * mu
  // k1 = k2 = 1
  Eigen::Vector2d prod_0 = inv_llt_sigma_folded_0 * mu_0;
  // k1 = 1, k2 = -1
  Eigen::Vector2d mu_1;
  mu_1 << mu_0 ( 0 ),  - mu_0 ( 1 );
  Eigen::Vector2d prod_1 = inv_llt_sigma_folded_1 * mu_1;
  // k1 = -1, k2 = 1
  Eigen::Vector2d mu_2 = - mu_1;
  Eigen::Vector2d prod_2 = inv_llt_sigma_folded_1 * mu_2;
  // k1 = k2 = -1
  Eigen::Vector2d mu_3 = - mu_0;
  Eigen::Vector2d prod_3 = inv_llt_sigma_folded_0 * mu_3;

  // Compute the vector of pdf and cdf
  Eigen::Vector2d pdf_0;
  pdf_0 << pdf ( z, prod_0 ( 0 ) ), pdf ( z, prod_0 ( 1 ) );
  Eigen::Vector2d pdf_1;
  pdf_1 << pdf ( z, prod_1 ( 0 ) ), pdf ( z, prod_1 ( 1 ) );
  Eigen::Vector2d pdf_2;
  pdf_2 << pdf ( z, prod_2 ( 0 ) ), pdf ( z, prod_2 ( 1 ) );
  Eigen::Vector2d pdf_3;
  pdf_3 << pdf ( z, prod_3 ( 0 ) ), pdf ( z, prod_3 ( 1 ) );
  Eigen::Vector2d cdf_0;
  cdf_0 << cdf ( z, prod_0 ( 0 ) ), cdf ( z, prod_0 ( 1 ) );
  Eigen::Vector2d cdf_1;
  cdf_1 << cdf ( z, prod_1 ( 0 ) ), cdf ( z, prod_1 ( 1 ) );
  Eigen::Vector2d cdf_2;
  cdf_2 << cdf ( z, prod_2 ( 0 ) ), cdf ( z, prod_2 ( 1 ) );
  Eigen::Vector2d cdf_3;
  cdf_3 << cdf ( z, prod_3 ( 0 ) ), cdf ( z, prod_3 ( 1 ) );

  // Compute the product between the cdf and pdf
  Eigen::Vector2d pf_0;
  pf_0 << pdf_0 ( 0 ) * cdf_0 ( 1 ), pdf_0 ( 1 ) * cdf_0 ( 0 );
  Eigen::Vector2d pf_1;
  pf_1 << pdf_1 ( 0 ) * cdf_1 ( 1 ), pdf_1 ( 1 ) * cdf_1 ( 0 );
  Eigen::Vector2d pf_2;
  pf_2 << pdf_2 ( 0 ) * cdf_2 ( 1 ), pdf_2 ( 1 ) * cdf_2 ( 0 );
  Eigen::Vector2d pf_3;
  pf_3 << pdf_3 ( 0 ) * cdf_3 ( 1 ), pdf_3 ( 1 ) * cdf_3 ( 0 );

  // Compute the terms I_k_i
  Eigen::Matrix2d i_0;
  i_0 << cdf_0 ( 1 ) * ( cdf_0 ( 0 ) - prod_0 ( 0 ) * pdf_0 ( 0 ) ), pdf_0 ( 0 ) * pdf_0 ( 1 ),
      pdf_0 ( 0 ) * pdf_0 ( 1 ), cdf_0 ( 0 ) * ( cdf_0 ( 1 ) - prod_0 ( 1 ) * pdf_0 ( 1 ) );
  Eigen::Matrix2d i_1;
  i_1 << cdf_1 ( 1 ) * ( cdf_1 ( 0 ) - prod_1 ( 0 ) * pdf_1 ( 0 ) ), pdf_1 ( 0 ) * pdf_1 ( 1 ),
      pdf_1 ( 0 ) * pdf_1 ( 1 ), cdf_1 ( 0 ) * ( cdf_1 ( 1 ) - prod_1 ( 1 ) * pdf_1 ( 1 ) );
  Eigen::Matrix2d i_2;
  i_2 << cdf_2 ( 1 ) * ( cdf_2 ( 0 ) - prod_2 ( 0 ) * pdf_2 ( 0 ) ), pdf_2 ( 0 ) * pdf_2 ( 1 ),
      pdf_2 ( 0 ) * pdf_2 ( 1 ), cdf_2 ( 0 ) * ( cdf_2 ( 1 ) - prod_2 ( 1 ) * pdf_2 ( 1 ) );
  Eigen::Matrix2d i_3;
  i_3 << cdf_3 ( 1 ) * ( cdf_3 ( 0 ) - prod_3 ( 0 ) * pdf_3 ( 0 ) ), pdf_3 ( 0 ) * pdf_3 ( 1 ),
      pdf_3 ( 0 ) * pdf_3 ( 1 ), cdf_3 ( 0 ) * ( cdf_3 ( 1 ) - prod_3 ( 1 ) * pdf_3 ( 1 ) );

  // Compute the first addend of the covariance matrix
  // k1 = k2 = 1
  Eigen::Matrix2d add_0_0 = mu_0 * mu_0.transpose() * cdf_0 ( 0 ) * cdf_0 ( 1 );
  // k1 = 1, k2 = -1
  Eigen::Matrix2d add_0_1 = mu_1 * mu_1.transpose() * cdf_1 ( 0 ) * cdf_1 ( 1 );
  // k1 = -1, k2 = 1
  Eigen::Matrix2d add_0_2 = mu_2 * mu_2.transpose() * cdf_2 ( 0 ) * cdf_2 ( 1 );
  // k1 = k2 = -1
  Eigen::Matrix2d add_0_3 = mu_3 * mu_3.transpose() * cdf_3 ( 0 ) * cdf_3 ( 1 );

  // Compute the second addend of the covariance matrix
  // k1 = k2 = 1
  Eigen::Matrix2d add_1_0 = llt_sigma_folded[0] * pf_0 * mu_0.transpose();
  // k1 = 1, k2 = -1
  Eigen::Matrix2d add_1_1 = llt_sigma_folded[1] * pf_1 * mu_1.transpose();
  // k1 = -1, k2 = 1
  Eigen::Matrix2d add_1_2 = llt_sigma_folded[1] * pf_2 * mu_2.transpose();
  // k1 = k2 = -1
  Eigen::Matrix2d add_1_3 = llt_sigma_folded[0] * pf_3 * mu_3.transpose();

  // Compute the forth addend of the covariance matrix
  // k1 = k2 = 1
  Eigen::Matrix2d add_3_0 = llt_sigma_folded[0] * i_0 * llt_sigma_folded[0].transpose();
  // k1 = 1, k2 = -1
  Eigen::Matrix2d add_3_1 = llt_sigma_folded[1] * i_1 * llt_sigma_folded[1].transpose();
  // k1 = -1, k2 = 1
  Eigen::Matrix2d add_3_2 = llt_sigma_folded[1] * i_2 * llt_sigma_folded[1].transpose();
  // k1 = k2 = -1
  Eigen::Matrix2d add_3_3 = llt_sigma_folded[0] * i_3 * llt_sigma_folded[0].transpose();

  // Compute the covariance matrix
  Eigen::Matrix2d cov = add_0_0 + add_1_0 + add_1_0.transpose() + add_3_0 +
                        add_0_1 + add_1_1 + add_1_1.transpose() + add_3_1 +
                        add_0_2 + add_1_2 + add_1_2.transpose() + add_3_2 +
                        add_0_3 + add_1_3 + add_1_3.transpose() + add_3_3 -
                        mu * mu.transpose();

  // Compute the variance of the expected improvment
  fval[0] = cov ( 0, 0 ) + cov ( 1,1 ) - cov ( 0,1 ) - cov ( 1,0 );

  return 0;

}

//! Compute the squre of the surface in the parametric point param if there is a polygonal boundary
int
SquareSurfacePointPoly ( unsigned ndim,
                         const double *x,
                         void *util_surf,
                         unsigned fdim,
                         double *fval )
{
  std::pair<Go::SplineSurface *, vector<Point> >* util = static_cast<std::pair<Go::SplineSurface*, vector<Point> >* > ( util_surf );

  if ( PnPoly ( x[0], x[1], util->second ) == true ) {
    Go::Point point;

    util->first->point ( point, x[0], x[1] );

    fval[0] = pow ( point[0], 2 );
  } else {
    fval[0] = 0.;
  }

  return 0; // success
}

//! Compute the expected value of the absolute value of the curve in the parametric point param if there is a polygonal boundary
int
AbsSurfPointPoly ( unsigned ndim,
                   const double *x,
                   void *util_surf,
                   unsigned fdim,
                   double *fval )
{

  std::pair< shared_ptr<Go::SplineSurface>, vector<Point> >* util = static_cast< std::pair< shared_ptr<Go::SplineSurface>, vector<Point> >* > ( util_surf );

  if ( PnPoly ( x[0], x[1], util->second ) == true ) {
    Go::Point point;

    util->first->point ( point, x[0], x[1] );

    fval[0] = std::abs ( point[0] );
  } else {
    fval[0] = 0.;
  }

  return 0; // success
}

//! Compute the expected value of the absolute value of the curve in the parametric point param if there is a polygonal boundary
int
EAbsSurfPointPoly ( unsigned ndim,
                    const double *x,
                    void *util_surf,
                    unsigned fdim,
                    double *fval )
{

  std::tuple< shared_ptr<Go::SplineSurface>, double, vector<Point> >* util = static_cast< std::tuple< shared_ptr<Go::SplineSurface>, double, vector<Point> >* > ( util_surf );

//   std::tuple_element<2,decltype ( util ) >::type polygon = std::get<2> ( *util );
  vector<Point> polygon = std::get<2> ( *util );

  if ( PnPoly ( x[0], x[1], polygon ) == true ) {

    Go::Point point;

    shared_ptr<Go::SplineSurface> surf = std::get<0> ( *util );
    double sd = std::get<1> ( *util );

    surf->point ( point, x[0], x[1] );

    double mean = point[0];
    double ratio = mean / sd;

    normal z;

    fval[0] = sd * std::sqrt ( 2 ) / boost::math::constants::root_pi<double>() * std::exp ( - std::pow ( ratio, 2 ) / 2 ) + mean * ( 1 - 2 * cdf ( z, - ratio ) );
  } else {
    fval[0] = 0.;
  }

  return 0;
}

//! Compute the variance of the absolute value of the curve in the parametric point param if there is a polygonal boundary
int
VarAbsSurfPointPoly ( unsigned ndim,
                      const double *x,
                      void *util_surf,
                      unsigned fdim,
                      double *fval )
{
  std::tuple< shared_ptr<Go::SplineSurface>, double, vector<Point> >* util = static_cast< std::tuple< shared_ptr<Go::SplineSurface>, double, vector<Point> >* > ( util_surf );

  vector<Point> polygon = std::get<2> ( *util );
  
  if ( PnPoly ( x[0], x[1], polygon ) == true ) {
  
    double sd = std::get<1> ( *util );
    shared_ptr<Go::SplineSurface> surf = std::get<0> ( *util );
    
    Go::Point point;

    surf->point ( point, x[0], x[1] );
    double mean = point[0];
    double ratio = mean / sd;

    normal z;

    fval[0] = std::pow ( mean, 2 ) + std::pow ( sd, 2 ) - std::pow ( sd * std::sqrt ( 2 ) / boost::math::constants::root_pi<double>() * std::exp ( - std::pow ( ratio, 2 ) / 2 ) + mean * ( 1 - 2 * cdf ( z, - ratio ) ), 2 );
  } else {
    fval[0] = 0.;
  }

  return 0;

}

//! Compute the expected value for the expected improvment if there is a polygonal boundary
int
MeanEiSurfPoly ( unsigned ndim,
                 const double *x,
                 void *util_surf,
                 unsigned fdim,
                 double *fval )
{
  std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> >* util = static_cast< std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> >* > ( util_surf );

  vector<Point> polygon = std::get<2> ( *util );
  
  if ( PnPoly ( x[0], x[1], polygon ) == true ) {
  
    vector< shared_ptr<Go::SplineSurface> > surf_ptr;
    surf_ptr.resize ( 2 );
    ( surf_ptr[0] ).reset ( new Go::SplineSurface );
    ( surf_ptr[1] ).reset ( new Go::SplineSurface );
    surf_ptr = std::get<0> ( *util );
    //   surf_ptr[0] = util [0];
    //   surf_ptr[1] = util [1];
    
    Eigen::LLT<MatrixXd> llt;
    vector<MatrixXd> llt_sigma_folded = std::get<1> ( *util );
    
    // Compute the expected value of the 2d folden normal distribution
    Eigen::Vector2d mu = MeanSurf ( surf_ptr, x, llt_sigma_folded );
    
    // Compute mu_0 - mu_1
    fval[0] = mu ( 0 ) - mu ( 1 );
    
  } else {
    fval[0] = 0.;  
  }
  
  return 0;
}

//! Compute the variance for the expected improvment if there is a polygonal boundary
int
VarianceEiSurfPoly ( unsigned ndim,
                     const double *x,
                     void *util_surf,
                     unsigned fdim,
                     double *fval )
{
  std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> >* util = static_cast< std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> >* > ( util_surf );

  vector<Point> polygon = std::get<2> ( *util );
  
  if ( PnPoly ( x[0], x[1], polygon ) == true ) {
  
    vector< shared_ptr<Go::SplineSurface> > surf_ptr;
    surf_ptr.resize ( 2 );
    ( surf_ptr[0] ).reset ( new Go::SplineSurface );
    ( surf_ptr[1] ).reset ( new Go::SplineSurface );
    surf_ptr = std::get<0> ( *util );
    //   surf_ptr[0] = util->first[0];
    //   surf_ptr[1] = util->first[1];
    
    Eigen::LLT<MatrixXd> llt;
    vector<MatrixXd> llt_sigma_folded = std::get<1> ( *util );
//     llt_sigma_folded.resize ( 2 );
//     llt_sigma_folded[0] = util->second[0];
//     llt_sigma_folded[1] = util->second[1];
    
    // Compute the inverse of llt_sigma_folded
    MatrixXd inv_llt_sigma_folded_0 = llt_sigma_folded[0].inverse();
    MatrixXd inv_llt_sigma_folded_1 = llt_sigma_folded[1].inverse();
    
    // Compute the mean on the param value
    Go::Point point_0, point_1;
    surf_ptr[0]->point ( point_0, x[0], x[1] );
    surf_ptr[1]->point ( point_1, x[0], x[1] );
    
    // Store values in a vector
    Eigen::Vector2d mu_0 ( 2 );
    mu_0 << point_0[0], point_1[0];
    
    // Compute the expected value of the 2d folden normal distribution
    Eigen::Vector2d mu = MeanSurf ( surf_ptr, x, llt_sigma_folded );
    
    normal z;
    
    // Compute values of inv_llt_sigma_folded * mu
    // k1 = k2 = 1
    Eigen::Vector2d prod_0 = inv_llt_sigma_folded_0 * mu_0;
    // k1 = 1, k2 = -1
    Eigen::Vector2d mu_1;
    mu_1 << mu_0 ( 0 ),  - mu_0 ( 1 );
    Eigen::Vector2d prod_1 = inv_llt_sigma_folded_1 * mu_1;
    // k1 = -1, k2 = 1
    Eigen::Vector2d mu_2 = - mu_1;
    Eigen::Vector2d prod_2 = inv_llt_sigma_folded_1 * mu_2;
    // k1 = k2 = -1
    Eigen::Vector2d mu_3 = - mu_0;
    Eigen::Vector2d prod_3 = inv_llt_sigma_folded_0 * mu_3;
    
    // Compute the vector of pdf and cdf
    Eigen::Vector2d pdf_0;
    pdf_0 << pdf ( z, prod_0 ( 0 ) ), pdf ( z, prod_0 ( 1 ) );
    Eigen::Vector2d pdf_1;
    pdf_1 << pdf ( z, prod_1 ( 0 ) ), pdf ( z, prod_1 ( 1 ) );
    Eigen::Vector2d pdf_2;
    pdf_2 << pdf ( z, prod_2 ( 0 ) ), pdf ( z, prod_2 ( 1 ) );
    Eigen::Vector2d pdf_3;
    pdf_3 << pdf ( z, prod_3 ( 0 ) ), pdf ( z, prod_3 ( 1 ) );
    Eigen::Vector2d cdf_0;
    cdf_0 << cdf ( z, prod_0 ( 0 ) ), cdf ( z, prod_0 ( 1 ) );
    Eigen::Vector2d cdf_1;
    cdf_1 << cdf ( z, prod_1 ( 0 ) ), cdf ( z, prod_1 ( 1 ) );
    Eigen::Vector2d cdf_2;
    cdf_2 << cdf ( z, prod_2 ( 0 ) ), cdf ( z, prod_2 ( 1 ) );
    Eigen::Vector2d cdf_3;
    cdf_3 << cdf ( z, prod_3 ( 0 ) ), cdf ( z, prod_3 ( 1 ) );
    
    // Compute the product between the cdf and pdf
    Eigen::Vector2d pf_0;
    pf_0 << pdf_0 ( 0 ) * cdf_0 ( 1 ), pdf_0 ( 1 ) * cdf_0 ( 0 );
    Eigen::Vector2d pf_1;
    pf_1 << pdf_1 ( 0 ) * cdf_1 ( 1 ), pdf_1 ( 1 ) * cdf_1 ( 0 );
    Eigen::Vector2d pf_2;
    pf_2 << pdf_2 ( 0 ) * cdf_2 ( 1 ), pdf_2 ( 1 ) * cdf_2 ( 0 );
    Eigen::Vector2d pf_3;
    pf_3 << pdf_3 ( 0 ) * cdf_3 ( 1 ), pdf_3 ( 1 ) * cdf_3 ( 0 );
    
    // Compute the terms I_k_i
    Eigen::Matrix2d i_0;
    i_0 << cdf_0 ( 1 ) * ( cdf_0 ( 0 ) - prod_0 ( 0 ) * pdf_0 ( 0 ) ), pdf_0 ( 0 ) * pdf_0 ( 1 ),
    pdf_0 ( 0 ) * pdf_0 ( 1 ), cdf_0 ( 0 ) * ( cdf_0 ( 1 ) - prod_0 ( 1 ) * pdf_0 ( 1 ) );
    Eigen::Matrix2d i_1;
    i_1 << cdf_1 ( 1 ) * ( cdf_1 ( 0 ) - prod_1 ( 0 ) * pdf_1 ( 0 ) ), pdf_1 ( 0 ) * pdf_1 ( 1 ),
    pdf_1 ( 0 ) * pdf_1 ( 1 ), cdf_1 ( 0 ) * ( cdf_1 ( 1 ) - prod_1 ( 1 ) * pdf_1 ( 1 ) );
    Eigen::Matrix2d i_2;
    i_2 << cdf_2 ( 1 ) * ( cdf_2 ( 0 ) - prod_2 ( 0 ) * pdf_2 ( 0 ) ), pdf_2 ( 0 ) * pdf_2 ( 1 ),
    pdf_2 ( 0 ) * pdf_2 ( 1 ), cdf_2 ( 0 ) * ( cdf_2 ( 1 ) - prod_2 ( 1 ) * pdf_2 ( 1 ) );
    Eigen::Matrix2d i_3;
    i_3 << cdf_3 ( 1 ) * ( cdf_3 ( 0 ) - prod_3 ( 0 ) * pdf_3 ( 0 ) ), pdf_3 ( 0 ) * pdf_3 ( 1 ),
    pdf_3 ( 0 ) * pdf_3 ( 1 ), cdf_3 ( 0 ) * ( cdf_3 ( 1 ) - prod_3 ( 1 ) * pdf_3 ( 1 ) );
    
    // Compute the first addend of the covariance matrix
    // k1 = k2 = 1
    Eigen::Matrix2d add_0_0 = mu_0 * mu_0.transpose() * cdf_0 ( 0 ) * cdf_0 ( 1 );
    // k1 = 1, k2 = -1
    Eigen::Matrix2d add_0_1 = mu_1 * mu_1.transpose() * cdf_1 ( 0 ) * cdf_1 ( 1 );
    // k1 = -1, k2 = 1
    Eigen::Matrix2d add_0_2 = mu_2 * mu_2.transpose() * cdf_2 ( 0 ) * cdf_2 ( 1 );
    // k1 = k2 = -1
    Eigen::Matrix2d add_0_3 = mu_3 * mu_3.transpose() * cdf_3 ( 0 ) * cdf_3 ( 1 );
    
    // Compute the second addend of the covariance matrix
    // k1 = k2 = 1
    Eigen::Matrix2d add_1_0 = llt_sigma_folded[0] * pf_0 * mu_0.transpose();
    // k1 = 1, k2 = -1
    Eigen::Matrix2d add_1_1 = llt_sigma_folded[1] * pf_1 * mu_1.transpose();
    // k1 = -1, k2 = 1
    Eigen::Matrix2d add_1_2 = llt_sigma_folded[1] * pf_2 * mu_2.transpose();
    // k1 = k2 = -1
    Eigen::Matrix2d add_1_3 = llt_sigma_folded[0] * pf_3 * mu_3.transpose();
    
    // Compute the forth addend of the covariance matrix
    // k1 = k2 = 1
    Eigen::Matrix2d add_3_0 = llt_sigma_folded[0] * i_0 * llt_sigma_folded[0].transpose();
    // k1 = 1, k2 = -1
    Eigen::Matrix2d add_3_1 = llt_sigma_folded[1] * i_1 * llt_sigma_folded[1].transpose();
    // k1 = -1, k2 = 1
    Eigen::Matrix2d add_3_2 = llt_sigma_folded[1] * i_2 * llt_sigma_folded[1].transpose();
    // k1 = k2 = -1
    Eigen::Matrix2d add_3_3 = llt_sigma_folded[0] * i_3 * llt_sigma_folded[0].transpose();
    
    // Compute the covariance matrix
    Eigen::Matrix2d cov = add_0_0 + add_1_0 + add_1_0.transpose() + add_3_0 +
    add_0_1 + add_1_1 + add_1_1.transpose() + add_3_1 +
    add_0_2 + add_1_2 + add_1_2.transpose() + add_3_2 +
    add_0_3 + add_1_3 + add_1_3.transpose() + add_3_3 -
    mu * mu.transpose();
    
    // Compute the variance of the expected improvment
    fval[0] = cov ( 0, 0 ) + cov ( 1,1 ) - cov ( 0,1 ) - cov ( 1,0 );
    
  } else {
    fval[0] = 0.;
  }

  return 0;

}

//! Objective function for the maximization of the expected improvment
double
ObjectiveFunction ( const vector<double> &x,
                    vector<double> &grad,
                    void* ego_obj )
{

  fkrig::EgoBase* ego = reinterpret_cast<fkrig::EgoBase*> ( ego_obj );

  RVectorXd coord_ego ( x.size() );

  // Fill the row of the second point
  for ( size_t i = 0; i < x.size (); ++i )
    coord_ego ( i ) = x[i];

  // Compute mean and variance of random variable
  double mean = ego->ComputeMean ( coord_ego );
  double sigma = std::sqrt ( ego->ComputeVariance ( coord_ego ) );

  double value = 0.;

  if ( sigma > 1e-16 ) {
    double ratio = mean / sigma;

    // Compute the value of the expected improvment
    value = sigma * boost::math::pdf ( ego->z_, ratio ) + mean * boost::math::cdf ( ego->z_, ratio );

  } else {
    value = std::max ( mean, 0. );
  }

  return value;
}

//! Objective function for the minization of the distance between the predicted function and the nominal function
double
ObjectiveFunctionMinCurve ( const vector<double> &x,
                            vector<double> &grad,
                            void* ego_obj )
{

  fkrig::EgoCurve* ego = reinterpret_cast<fkrig::EgoCurve*> ( ego_obj );

  RVectorXd coord ( x.size() );

  // Fill the row of the second point
  for ( size_t i = 0; i < x.size (); ++i )
    coord ( i ) = x[i];

  // Obtain fkrig object
  shared_ptr<CurveBase> fkrig = ego->get_fkrig();

  // Predict the curve
  Go::SplineCurve curve = fkrig->Predict ( coord );

  // Compute mean and variance of random variable
  double mean = ego->ComputeL1Mean ( curve, coord );
//   double mean = ego->ComputeL1 ( curve );

  return mean;

}

//! Objective function for the minization of the distance between the predicted function and the nominal function
double
ObjectiveFunctionMinSurf ( const vector<double> &x,
                           vector<double> &grad,
                           void* ego_obj )
{

  fkrig::EgoSurf* ego = reinterpret_cast<fkrig::EgoSurf*> ( ego_obj );

  RVectorXd coord ( x.size() );

  // Fill the row of the second point
  for ( size_t i = 0; i < x.size (); ++i )
    coord ( i ) = x[i];

  // Obtain fkrig object
  shared_ptr<SurfBase> fkrig = ego->get_fkrig();

  // Predict the curve
  Go::SplineSurface curve = fkrig->Predict ( coord );

  // Compute mean and variance of random variable
  double mean = ego->ComputeL1Mean ( curve, coord );
//   double mean = ego->ComputeL1 ( curve );

  return mean;

}

//! Check if a point is on the polygon
bool
PnPoly ( double u,
         double v,
         vector<Point> polygon )
{

  bool in_poly = false;

  for ( size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++ ) {
    if ( ( ( polygon[i].v >= v ) != ( polygon[j].v >= v ) ) && ( u <= ( polygon[j].u - polygon[i].u ) * ( v - polygon[i].v ) / ( polygon[j].v - polygon[i].v ) + polygon[i].u ) )
      in_poly = !in_poly;
  }

  return in_poly;

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
