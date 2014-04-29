#include "ego_curve.hpp"
#include "util_fkrig.hpp"

namespace fkrig {

fkrig::EgoCurve::EgoCurve ( std::shared_ptr<fkrig::CurveBase>& f_curve,
                            std::shared_ptr<Go::SplineCurve>& nominal_curve,
                            vector<double> lb,
                            vector<double> ub,
                            nlopt::algorithm glob_alg,
                            nlopt::algorithm loc_alg,
                            double tol_glob,
                            double tol_loc,
                            double max_iter_glob,
                            double max_iter_loc )
  : EgoBase ( lb, ub, glob_alg, loc_alg, tol_glob, tol_loc, max_iter_glob, max_iter_loc ), f_curve_ ( f_curve ), nominal_curve_ ( nominal_curve )
{
  // Compute range
  ComputeUniqueRange ();

  // Compute min
  ComputeMin ();
};

//! Compute the L1 distance to the nominal funciton
double
fkrig::EgoCurve::ComputeL1 ( Go::SplineCurve& curve ) const
{

  // Compute the difference between the curve and the nominal curve
  shared_ptr<Go::SplineCurve> diff = Go::GeometryTools::curveSum ( curve, 1, *nominal_curve_, -1 );

  // Compute
  double value = fkrig::adaptiveSimpsons ( fkrig::abs_curve_point, *diff, range_points_.first, range_points_.second, 1e-6, 10 );

  return value;
}

//! Compute the expected value of the L1 distance to the nominal funciton
double
fkrig::EgoCurve::ComputeL1Mean ( Go::SplineCurve& curve,
                                 MatrixXd coord ) const
{
  // Compute the difference between the curve and the nominal curve
  shared_ptr<Go::SplineCurve> diff = Go::GeometryTools::curveSum ( curve, 1, *nominal_curve_, -1 );

  // Obtain range
  std::pair<double, double> range = f_curve_->get_range();

  // Compute the standard deviation
  MatrixXd temp = f_curve_->PredictCovariance ( coord ) / ( range.second - range.first );
  // Check if is positive
  if ( temp ( 0,0 ) < 0. ) {
    temp ( 0,0 ) = 0.;
  }

  double sd = std::sqrt ( temp ( 0,0 ) );

  double value = 0.;
  // Compute
  if ( sd < 1e-12 ) {
    value = fkrig::adaptiveSimpsons ( fkrig::abs_curve_point, *diff, range_points_.first, range_points_.second, 1e-6, 10 );
  } else {
    value = fkrig::adaptiveSimpsons ( fkrig::e_abs_curve_point, *diff, sd, range_points_.first, range_points_.second, 1e-6, 10 );
  }

  return value;
}

//! Compute the value of the mean at the design coordinates coord
double
fkrig::EgoCurve::ComputeMean ( RVectorXd coord ) const
{

  double value = 0.;

  // Check if coord is equal to the minimum
  vector<bool> comp ( coord.cols(), false );
  for ( size_t i = 0; i < coord.cols(); ++i )
    comp[i] = EgoBase::coord_ego_ ( 0, i ) == coord ( i );

  if ( std::any_of ( comp.begin(), comp.end(), [] ( int i ) {
  return ( i == false );
  } ) ) {

    // Predict in the coordinate coord
    Go::SplineCurve curve = f_curve_->Predict ( coord );

    MatrixXd coord_ego ( EgoBase::coord_ego_.rows(), EgoBase::coord_ego_.cols() );

    coord_ego.row ( 0 ) = EgoBase::coord_ego_.row ( 0 );
    coord_ego.row ( 1 ) = coord.row ( 0 );

//   for ( size_t i = 0; i < coord.cols(); ++i )
//     EgoBase::coord_ego_(1,i) = coord(i);

    // Create vector of shared pointer to curves
    vector< shared_ptr<Go::SplineCurve> > curves_ptr;
    curves_ptr.resize ( 2 );
    ( curves_ptr[0] ).reset ( new Go::SplineCurve );
    ( curves_ptr[1] ).reset ( new Go::SplineCurve );
    curves_ptr[0] = Go::GeometryTools::curveSum ( curve_min_, 1, *nominal_curve_, -1 );
    curves_ptr[1] = Go::GeometryTools::curveSum ( curve, 1, *nominal_curve_, -1 );

    // Compute the covariance matrix between the points of coord (s1 = s2 = 1)
    MatrixXd sigma_folded ( 2, 2 );
    sigma_folded = f_curve_->PredictCovariance ( coord_ego );

    // Obtain range
    std::pair<double, double> range = f_curve_->get_range();

    // Compute the covariance for each s
    sigma_folded /= ( range.second - range.first );

    if ( sigma_folded ( 0,0 ) < 1e-12 && sigma_folded ( 1,1 ) < 1e-12 ) {

      value = fkrig::adaptiveSimpsons ( fkrig::abs_curve_point, * ( curves_ptr[0] ), range_points_.first, range_points_.second, 1e-6, 10 );
      value -= fkrig::adaptiveSimpsons ( fkrig::abs_curve_point, * ( curves_ptr[1] ), range_points_.first, range_points_.second, 1e-6, 10 );

    } else {

      Eigen::LLT<MatrixXd> llt;
      vector<MatrixXd> llt_sigma_folded;
      llt_sigma_folded.resize ( 2 );
      llt_sigma_folded[0].resize ( coord.rows(), coord.rows() );
      llt_sigma_folded[1].resize ( coord.rows(), coord.rows() );

      llt.compute ( sigma_folded );
      llt_sigma_folded[0] = llt.matrixL();

      // Compute the covariance matrix with s1 = 1, s2 = -1
      sigma_folded ( 0,1 ) = - sigma_folded ( 0,1 );
      sigma_folded ( 1,0 ) = - sigma_folded ( 1,0 );

      llt.compute ( sigma_folded );
      llt_sigma_folded[1] = llt.matrixL();

      value = fkrig::adaptiveSimpsons ( &fkrig::MeanEiCurve, curves_ptr, llt_sigma_folded, range_points_.first, range_points_.second, 1e-6, 10 );

    }

  }

  return value;
}

//! Compute the value of the variance at the design coordinates coord
double
fkrig::EgoCurve::ComputeVariance ( RVectorXd coord ) const
{

  double value = 0.;

  // Check if coord is equal to the minimum
  vector<bool> comp ( coord.cols(), false );
  for ( size_t i = 0; i < coord.cols(); ++i )
    comp[i] = EgoBase::coord_ego_ ( 0, i ) == coord ( i );

  if ( std::any_of ( comp.begin(), comp.end(), [] ( int i ) {
  return ( i == false );
  } ) ) {

    // Predict in the coordinate coord
    Go::SplineCurve curve = f_curve_->Predict ( coord );

    MatrixXd coord_ego ( EgoBase::coord_ego_.rows(), EgoBase::coord_ego_.cols() );

    coord_ego.row ( 0 ) = EgoBase::coord_ego_.row ( 0 );
    coord_ego.row ( 1 ) = coord.row ( 0 );

//   for ( size_t i = 0; i < coord.cols(); ++i )
//     EgoBase::coord_ego_(1,i) = coord(i);

    // Create vector of shared pointer to curves
    vector< shared_ptr<Go::SplineCurve> > curves_ptr;
    curves_ptr.resize ( 2 );
    ( curves_ptr[0] ).reset ( new Go::SplineCurve );
    ( curves_ptr[1] ).reset ( new Go::SplineCurve );
    curves_ptr[0] = Go::GeometryTools::curveSum ( curve_min_, 1, *nominal_curve_, -1 );
    curves_ptr[1] = Go::GeometryTools::curveSum ( curve, 1, *nominal_curve_, -1 );

    // Compute the covariance matrix between the points of coord (s1 = s2 = 1)
    MatrixXd sigma_folded ( 2, 2 );
    sigma_folded = f_curve_->PredictCovariance ( coord_ego );

    // Obtain range
    std::pair<double, double> range = f_curve_->get_range();

    // Compute the covariance for each s
    sigma_folded /= ( range.second - range.first );

    if ( sigma_folded ( 0,0 ) >= 1e-12 || sigma_folded ( 1,1 ) >= 1e-12 ) {

      Eigen::LLT<MatrixXd> llt;
      vector<MatrixXd> llt_sigma_folded;
      llt_sigma_folded.resize ( 2 );
      llt_sigma_folded[0].resize ( coord.rows(), coord.rows() );
      llt_sigma_folded[1].resize ( coord.rows(), coord.rows() );

      llt.compute ( sigma_folded );
      llt_sigma_folded[0] = llt.matrixL();

      // Compute the covariance matrix with s1 = 1, s2 = -1
      sigma_folded ( 0,1 ) = - sigma_folded ( 0,1 );
      sigma_folded ( 1,0 ) = - sigma_folded ( 1,0 );

      llt.compute ( sigma_folded );
      llt_sigma_folded[1] = llt.matrixL();

      value = fkrig::adaptiveSimpsons ( &fkrig::VarianceEiCurve, curves_ptr, llt_sigma_folded, range_points_.first, range_points_.second, 1e-6, 10 );

    }

  }

  return value;  
  
}

//! Compute the expected improvment in location coord
double
fkrig::EgoCurve::ComputeFunction ( RVectorXd coord )
{
  return 0.;
}

//! Compute the L1 distance to the nominal funciton
void
fkrig::EgoCurve::ComputeMin ()
{
  // Obtain the matrix of coordinates
  MatrixXd coord = f_curve_->get_coord();

  // Predict the curves in the design locations
  vector<Go::SplineCurve> curves = f_curve_->Predict ( coord );

  // Compute the expected value of the L1 distance between the predicted curves and the nominal curve
  vector<double> distance ( coord.rows (), 0. );
  for ( size_t i = 0; i < distance.size (); ++i )
    distance[i] = ComputeL1Mean ( curves[i], coord.row ( i ) );

  // Find the index of the min
  EgoBase::index_min_ = ( std::min_element ( distance.begin (), distance.end () ) - distance.begin () );

  // Save the geometric coordinate of the minimum
  EgoBase::coord_ego_.row ( 0 ) = coord.row ( EgoBase::index_min_ );
  // Save the predicted curve with miminum expected L1 distance
  curve_min_ = curves[EgoBase::index_min_];
}

//! Compute the range of the knots vector for the curves
void
fkrig::EgoCurve::ComputeUniqueRange ()
{

  // Obtain the range of the curves
  std::pair<double, double> range_curves = f_curve_->get_range();
  // Obtain the range of the nominal curve
  double min_par = nominal_curve_->startparam();
  double max_par = nominal_curve_->endparam();

  // The range is computed as the maximum of the minimum value and the minimum of the maximum values
  range_points_ = std::make_pair ( std::max ( range_curves.first, min_par ), std::min ( range_curves.second, max_par ) );
}

} // End of namespace
