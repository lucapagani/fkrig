#include "curve_link.hpp"
#include "util_fkrig.hpp"

using std::vector;
// using namespace fkrig::CurveLinkBase;

namespace fkrig {

// Constructor with implicit parameterization (only if the dimension of the geometric space is 2)
fkrig::CurveLink::CurveLink ( vector< vector<double> >& points,
                              shared_ptr<CurveBase> f_curve,
                              const int dim,
                              MatrixXd coord,
                              model_type model,
                              double n_max,
                              double tol )
{

  // Check if the dimension is 2
  if ( dim != 2 ) {
    std::cerr << "The implicit parameterization is possible only if the dimension of the geometric space is 2" << std::endl;
  }

  // Fill the param_ vector
  size_t n_points = 0;
  vector< vector <double> >param;
  param.resize ( points.size() );
  for ( size_t i = 0; i < param.size(); ++i ) {
    // Find the number of points of the i-th dataset
    n_points = points[i].size() / dim;
    // Resize the vector of parameters of the i-th dataset
    param[i].resize ( n_points );
    for ( size_t j = 0, k = 0; j < param[i].size(); ++j, k+=dim )
      param[i][j] = points[i][k];
  }

//   curve::curve ( points, param, dim, coord, model, n_max, tol );

}

// Constructor with parameterization as input
fkrig::CurveLink::CurveLink ( vector< vector<double> >& points,
                              vector< vector<double> >& param,
                              shared_ptr<CurveBase> f_curve,
                              const int dim,
                              MatrixXd coord,
                              model_type model,
                              double n_max,
                              double tol )
  : Curve ( points, param, dim, coord, model, n_max, tol ) //, f_curve_ ( f_curve )
{

  // Predidict in the coordinates coord
  vector<Go::SplineCurve> f_curve_low = f_curve->Predict ( coord );

  // Find the difference between the high fidelity points and the prediction with the low fidelity model
  vector< vector<double> > diff_points;
  diff_points.resize ( points.size () );

  Go::Point point;

  for ( size_t i = 0; i < param.size (); ++i ) {
    diff_points[i].reserve ( param[i].size () );
    for ( size_t j = 0; j < param[i].size (); ++j ) {
      // Predict point at param
      f_curve_low[i].point ( point, param[i][j] );
      diff_points[i].push_back ( points[i][j] - point[0] );
    }
  }

  Curve::points_ = diff_points;

  CurveBase::f_curve_ = f_curve;

}

// Compute the functional krigin parameters
void
fkrig::CurveLink::Eval()
{
 
  Curve::Eval ();

}

// TODO: prediction with dim_ > 1
// Predict the value at the design coordinates coord in the geometric value param
vector<vector<double> >
fkrig::CurveLink::Predict ( MatrixXd& coord,
                            vector<double>& param ) const
{

  vector< vector<double> > points = Curve::Predict ( coord, param );
  
  // Compute the predicted curves
  vector<Go::SplineCurve> curve = CurveBase::f_curve_->Predict ( coord );

  // Compute the prediction at a parametric points
  Go::Point int_point;

  for ( size_t i = 0; i < points.size(); ++i ) {
    points[i].resize ( param.size() * dim_ );
    for ( size_t j = 0; j < points[i].size(); ++j ) {
      curve[i].point ( int_point, param[j] );
      points[i][j] += int_point[0];
    }
  }

  return points;
}

// Predict the spline curve at the design coordinates coord
Go::SplineCurve
fkrig::CurveLink::Predict ( RVectorXd& coord ) const
{

  Go::SplineCurve temp_curve_0 = Curve::Predict ( coord );
  Go::SplineCurve temp_curve_1 = CurveBase::f_curve_->Predict ( coord );
  
  shared_ptr<Go::SplineCurve> curve = Go::GeometryTools::curveSum ( temp_curve_0, 1, temp_curve_1, 1 );

  return *curve;
}

// Predict the spline curve at the design coordinates coord
vector<Go::SplineCurve>
fkrig::CurveLink::Predict ( MatrixXd& coord ) const
{

  vector<Go::SplineCurve> pred;
  pred.reserve ( coord.rows () );

  RVectorXd coord_v;
  
  for ( size_t i = 0; i < coord.rows (); ++i ) {
    coord_v = coord.row (i);
    pred.push_back ( Predict ( coord_v ) );
  }
  
  return pred;
  
}

//! Predicted the total covarince matrix at the design coordinates coord
MatrixXd
fkrig::CurveLink::PredictCovariance ( MatrixXd& coord ) const
{

  MatrixXd cov = Curve::PredictCovariance ( coord );

  return cov;
}

// // Compute the matrix of coefficients of the term a with the Least Square method
// void
// fkrig::CurveLink::ComputeALs ()
// {
// 
//   Curve::ComputeALs ();
//   
// }
// 
// // Compute the matrix of coefficients of the term a with the Generalized Least Square method
// void
// fkrig::CurveLink::ComputeAGls ()
// {
// 
//   Curve::ComputeAGls ();
// 
// }
// 
// // Compute the matrix of coefficients of the term m with the Least Square method
// void
// fkrig::CurveLink::ComputeMLs ()
// {
// 
//   Curve::ComputeMLs ();
//   
// }
// 
// // Compute the matrix of coefficients of the term m with the Generalized Least Square method
// void
// fkrig::CurveLink::ComputeMGls ()
// {
// 
//   Curve::ComputeMGls ();
//   
// }
// 
// // Vector of coefficients of the difference between the curves and the splines of the term m
// void
// fkrig::CurveLink::ComputeDiffCoefs ()
// {
// 
//   Curve::ComputeDiffCoefs ();
// 
// }
// 
// // Compute the vector of the square of the parwise distances
// void
// fkrig::CurveLink::ComputeSqPairwiseDistances ()
// {
// 
//   Curve::ComputeSqPairwiseDistances ();
// 
// }

} // end of namespace
