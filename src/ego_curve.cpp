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
  double value = fkrig::AdaptiveSimpsons ( fkrig::AbsCurvePoint, *diff, range_points_.first, range_points_.second, 1e-6, 10 );

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
    value = fkrig::AdaptiveSimpsons ( fkrig::AbsCurvePoint, *diff, range_points_.first, range_points_.second, 1e-6, 10 );
  } else {
    value = fkrig::AdaptiveSimpsons ( fkrig::EAbsCurvePoint, *diff, sd, range_points_.first, range_points_.second, 1e-6, 10 );
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

    MatrixXd coord_ego ( 2, EgoBase::coord_ego_.cols() );

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

    if ( ( sigma_folded ( 0,0 ) < 1e-12 && sigma_folded ( 1,1 ) < 1e-12 ) || sigma_folded.determinant () < 1e-6 ) {

      value = fkrig::AdaptiveSimpsons ( fkrig::AbsCurvePoint, * ( curves_ptr[0] ), range_points_.first, range_points_.second, 1e-6, 10 );
      value -= fkrig::AdaptiveSimpsons ( fkrig::AbsCurvePoint, * ( curves_ptr[1] ), range_points_.first, range_points_.second, 1e-6, 10 );

    } else if ( sigma_folded ( 0,0 ) < 1e-12 ) {

      value = fkrig::AdaptiveSimpsons ( fkrig::AbsCurvePoint, * ( curves_ptr[0] ), range_points_.first, range_points_.second, 1e-6, 10 );
      value -= fkrig::AdaptiveSimpsons ( fkrig::EAbsCurvePoint, * ( curves_ptr[1] ), sigma_folded ( 1,1 ), range_points_.first, range_points_.second, 1e-6, 10 );

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

      value = fkrig::AdaptiveSimpsons ( &fkrig::MeanEiCurve, curves_ptr, llt_sigma_folded, range_points_.first, range_points_.second, 1e-6, 10 );

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

    MatrixXd coord_ego ( 2, EgoBase::coord_ego_.cols() );

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

//     std::cout << sigma_folded.determinant() << "\n";

    if ( ( sigma_folded ( 0,0 ) >= 1e-12 && sigma_folded ( 1,1 ) >= 1e-12 ) && sigma_folded.determinant () >= 1e-6 ) {

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

      value = fkrig::AdaptiveSimpsons ( &fkrig::VarianceEiCurve, curves_ptr, llt_sigma_folded, range_points_.first, range_points_.second, 1e-6, 10 );

//     } else if ( sigma_folded ( 1,1 ) >= 1e-12 && sigma_folded.determinant () >= 1e-6 ) {
    } else if ( sigma_folded ( 0,0 ) >= 1e-12 || sigma_folded ( 1,1 ) >= 1e-12 ){

      value = fkrig::AdaptiveSimpsons ( &fkrig::VarAbsCurvePoint, * ( curves_ptr[0] ), sigma_folded ( 1,1 ), range_points_.first, range_points_.second, 1e-6, 10 );
      value += fkrig::AdaptiveSimpsons ( &fkrig::VarAbsCurvePoint, * ( curves_ptr[1] ), sigma_folded ( 1,1 ), range_points_.first, range_points_.second, 1e-6, 10 );

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

  // Check if it is a link model
  bool stop = false;
  shared_ptr<CurveBase> temp_ptr;
//   temp_ptr = f_curve_;
  MatrixXd temp_coord;
  temp_ptr = f_curve_->get_f_curve ();
  
  while ( stop == false ) {
    if ( ( temp_ptr != NULL ) == true ) {
      // Obtain the coordinates of the nested object
      temp_coord = temp_ptr->get_coord ();
      // Resize matrix
      coord.conservativeResize ( coord.rows () + temp_coord.rows (), Eigen::NoChange );
      // Bind columns
      coord.block ( coord.rows () - temp_coord.rows (), 0, temp_coord.rows (), 2 ) = temp_coord;
      // Obtain the nested object
      temp_ptr = temp_ptr->get_f_curve ();     
    } else {
      stop = true;
    }
  }
  
  // Predict the curves in the design locations
  vector<Go::SplineCurve> curves = f_curve_->Predict ( coord );

  // Compute the expected value of the L1 distance between the predicted curves and the nominal curve
  vector<double> distance ( coord.rows (), 0. );
  for ( size_t i = 0; i < distance.size (); ++i )
    distance[i] = ComputeL1Mean ( curves[i], coord.row ( i ) );

  // Find the index of the min
  EgoBase::index_min_ = ( std::min_element ( distance.begin (), distance.end () ) - distance.begin () );

  // Save the geometric coordinate of the minimum
  EgoBase::coord_ego_ = coord.row ( EgoBase::index_min_ );
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

//! Find the design coordinate that minimize the distance between the predicted surface and the nominal surface
void
fkrig::EgoCurve::ComputeMinDist ()
{
  
  // Create a global optimization object
  nlopt::opt opt_glob ( EgoBase::glob_alg_, EgoBase::lb_.size() );
  // Create a local optimization object
  nlopt::opt opt_loc ( EgoBase::loc_alg_, EgoBase::lb_.size() );

  nlopt::vfunc f = &fkrig::ObjectiveFunctionMinCurve;

  // Set bounds
  opt_glob.set_lower_bounds ( EgoBase::lb_ );
  opt_glob.set_upper_bounds ( EgoBase::ub_ );
  opt_loc.set_lower_bounds ( EgoBase::lb_ );
  opt_loc.set_upper_bounds ( EgoBase::ub_ );

  // Set the objective function
  opt_glob.set_min_objective ( f, this );  
  
  // Set the relative x tollerance
  opt_glob.set_xtol_rel ( EgoBase::x_tol_glob_ );
  opt_loc.set_xtol_rel ( EgoBase::x_tol_loc_ );

  // Set the maximum number of iterations
  opt_glob.set_maxeval ( EgoBase::max_iter_glob_ );
  opt_loc.set_maxeval ( EgoBase::max_iter_loc_ );

  opt_glob.set_local_optimizer ( opt_loc );

  // Chose a starting point
  std::vector<double> x0 ( EgoBase::lb_.size(), 0. );
  for ( size_t i = 0; i < x0.size(); ++i )
    x0[i] = ( EgoBase::lb_[i] + EgoBase::ub_[i] ) / 2;

  // Preform the optimization
  EgoBase::result_min_ = opt_glob.optimize ( x0, EgoBase::value_min_ );

  EgoBase::x_min_ = x0;  
  
}

} // End of namespace
