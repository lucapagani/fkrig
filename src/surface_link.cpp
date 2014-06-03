#include "surface_link.hpp"
#include "util_fkrig.hpp"

using std::vector;

namespace fkrig {

// Constructor with parameterization as input
fkrig::SurfLink::SurfLink ( vector< vector<double> >& points,
                            vector< vector<double> >& param_u,
                            vector< vector<double> >& param_v,
                            shared_ptr<SurfBase> f_surf,
                            const int dim,
                            MatrixXd coord,
                            model_type model,
                            double n_max,
                            double tol )
  : Surf ( points, param_u, param_v, dim, coord, model, n_max, tol )
{

  // Predidict in the coordinates coord
  vector<Go::SplineSurface> f_surf_low = f_surf->Predict ( coord );

  // Find the difference between the high fidelity points and the prediction with the low fidelity model
  vector< vector<double> > diff_points;
  diff_points.resize ( points.size () );

  Go::Point point;

  for ( size_t i = 0; i < param_u.size (); ++i ) {
    diff_points[i].reserve ( param_u[i].size () );
    for ( size_t j = 0; j < param_u[i].size (); ++j ) {
      // Predict point at param
      f_surf_low[i].point ( point, param_u[i][j], param_v[i][j] );
      diff_points[i].push_back ( points[i][j] - point[0] );
    }
  }

  Surf::points_ = diff_points;

  SurfBase::f_surf_ = f_surf;

}

// Compute the functional krigin parameters
void
fkrig::SurfLink::Eval()
{

  Surf::Eval ();

}

// TODO: prediction with dim_ > 1
// Predict the value at the design coordinates coord in the geometric value param
vector<vector<double> >
fkrig::SurfLink::Predict ( MatrixXd& coord,
                           vector<double>& param_u,
                           vector<double>& param_v ) const
{

  vector< vector<double> > points = Surf::Predict ( coord, param_u, param_v );

  // Compute the predicted curves
  vector<Go::SplineSurface> surf = SurfBase::f_surf_->Predict ( coord );

  // Compute the prediction at a parametric points
  Go::Point int_point;

  for ( size_t i = 0; i < points.size(); ++i ) {
    points[i].resize ( param_u.size() * dim_ );
    for ( size_t j = 0; j < points[i].size(); ++j ) {
      surf[i].point ( int_point, param_u[j], param_v[j] );
      points[i][j] += int_point[0];
    }
  }

  return points;
}

// Predict the spline curve at the design coordinates coord
Go::SplineSurface
fkrig::SurfLink::Predict ( RVectorXd& coord ) const
{

  Go::SplineSurface temp_surf_0 = Surf::Predict ( coord );
  Go::SplineSurface temp_surf_1 = SurfBase::f_surf_->Predict ( coord );

  shared_ptr<Go::SplineSurface> surf = Go::GeometryTools::surfaceSum ( temp_surf_0, 1, temp_surf_1, 1 );

  return *surf;
}

// Predict the spline curve at the design coordinates coord
vector<Go::SplineSurface>
fkrig::SurfLink::Predict ( MatrixXd& coord ) const
{

  vector<Go::SplineSurface> pred;
  pred.reserve ( coord.rows () );

  RVectorXd coord_v;

  for ( size_t i = 0; i < coord.rows (); ++i ) {
    coord_v = coord.row ( i );
    pred.push_back ( Predict ( coord_v ) );
  }

  return pred;

}

//! Predicted the total covarince matrix at the design coordinates coord
MatrixXd
fkrig::SurfLink::PredictCovariance ( MatrixXd& coord ) const
{

  MatrixXd cov = Surf::PredictCovariance ( coord );

  return cov;
}

} // end of namespace
