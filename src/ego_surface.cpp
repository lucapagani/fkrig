#include "ego_surface.hpp"
#include "util_fkrig.hpp"

namespace fkrig {

fkrig::EgoSurf::EgoSurf ( std::shared_ptr<fkrig::SurfBase>& f_surf,
                          std::shared_ptr<Go::SplineSurface>& nominal_surf,
                          vector<double> lb,
                          vector<double> ub,
                          nlopt::algorithm glob_alg,
                          nlopt::algorithm loc_alg,
                          double tol_glob,
                          double tol_loc,
                          double max_iter_glob,
                          double max_iter_loc )
  : EgoBase ( lb, ub, glob_alg, loc_alg, tol_glob, tol_loc, max_iter_glob, max_iter_loc ), f_surf_ ( f_surf ), nominal_surf_ ( nominal_surf ), polygon_ ( f_surf->get_polygon () )
{
  // Compute range
  ComputeUniqueRange ();

  // Compute min
  ComputeMin ();
};

//! Compute the L1 distance to the nominal function
double
fkrig::EgoSurf::ComputeL1 ( Go::SplineSurface& surf ) const
{

  // Compute the difference between the curve and the nominal curve
  shared_ptr<Go::SplineSurface> diff = Go::GeometryTools::surfaceSum ( surf, 1, *nominal_surf_, -1 );

  // Compute
  double value ( 0. ), err ( 0. );
  double range_min[2] = { EgoSurf::range_points_u_.first, EgoSurf::range_points_v_.first };
  double range_max[2] = { EgoSurf::range_points_u_.second, EgoSurf::range_points_v_.second };
  if ( polygon_.empty () == true ) {
    hcubature ( 1, fkrig::AbsSurfPoint, diff.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
  } else {
    std::pair< shared_ptr<Go::SplineSurface>, vector<Point> >* util_surf ( new std::pair< shared_ptr<Go::SplineSurface>, vector<Point> > ( diff, polygon_ ) );
    hcubature ( 1, fkrig::AbsSurfPointPoly, util_surf, 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
    delete util_surf;
  }

  return value;
}

//! Compute the expected value of the L1 distance to the nominal function
double
fkrig::EgoSurf::ComputeL1Mean ( Go::SplineSurface& surf,
                                MatrixXd coord ) const
{
  // Compute the difference between the curve and the nominal curve
  shared_ptr<Go::SplineSurface> diff = Go::GeometryTools::surfaceSum ( surf, 1, *nominal_surf_, -1 );

  // Obtain area of the geometric domain
  double area = f_surf_->get_domain_range ();

  // Compute the standard deviation
  MatrixXd temp = f_surf_->PredictCovariance ( coord ) / area;
  // Check if is positive
  if ( temp ( 0,0 ) < 0. ) {
    temp ( 0,0 ) = 0.;
  }

  double sd = std::sqrt ( temp ( 0,0 ) );

  // Compute
  double value ( 0. ), err ( 0. );
  double range_min[2] = { EgoSurf::range_points_u_.first, EgoSurf::range_points_v_.first };
  double range_max[2] = { EgoSurf::range_points_u_.second, EgoSurf::range_points_v_.second };
  if ( sd < 1e-12 ) {
    if ( polygon_.empty () == true ) {
      hcubature ( 1, fkrig::AbsSurfPoint, diff.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
    } else {
      std::pair< shared_ptr<Go::SplineSurface>, vector<Point> >* util_surf ( new std::pair< shared_ptr<Go::SplineSurface>, vector<Point> > ( diff, polygon_ ) );
      hcubature ( 1, fkrig::AbsSurfPointPoly, util_surf, 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
      delete util_surf;
    }
  } else {
    if ( polygon_.empty () == true ) {
      shared_ptr< std::pair< shared_ptr<Go::SplineSurface>, double > > util ( new std::pair< shared_ptr<Go::SplineSurface>, double > ( diff, sd ) );
      hcubature ( 1, fkrig::EAbsSurfPoint, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
    } else {
      shared_ptr< std::tuple< shared_ptr<Go::SplineSurface>, double, vector<Point> > > util ( new std::tuple< shared_ptr<Go::SplineSurface>, double, vector<Point> > ( diff, sd, polygon_ ) );
      hcubature ( 1, fkrig::EAbsSurfPoint, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
    }
  }

  return value;
}

//! Compute the value of the mean at the design coordinates coord
double
fkrig::EgoSurf::ComputeMean ( RVectorXd coord ) const
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
    Go::SplineSurface surf = f_surf_->Predict ( coord );

    MatrixXd coord_ego ( 2, EgoBase::coord_ego_.cols() );

    coord_ego.row ( 0 ) = EgoBase::coord_ego_.row ( 0 );
    coord_ego.row ( 1 ) = coord.row ( 0 );

//   for ( size_t i = 0; i < coord.cols(); ++i )
//     EgoBase::coord_ego_(1,i) = coord(i);

    // Create vector of shared pointer to surfaces
    vector< shared_ptr<Go::SplineSurface> > surf_ptr;
    surf_ptr.resize ( 2 );
    ( surf_ptr[0] ).reset ( new Go::SplineSurface );
    ( surf_ptr[1] ).reset ( new Go::SplineSurface );
    surf_ptr[0] = Go::GeometryTools::surfaceSum ( surf_min_, 1, *nominal_surf_, -1 );
    surf_ptr[1] = Go::GeometryTools::surfaceSum ( surf, 1, *nominal_surf_, -1 );

    // Compute the covariance matrix between the points of coord (s1 = s2 = 1)
    MatrixXd sigma_folded ( 2, 2 );
    sigma_folded = f_surf_->PredictCovariance ( coord_ego );

    // Obtain area of the geometric domain
    double area = f_surf_->get_domain_range ();

    // Compute the covariance for each s
    sigma_folded /= area;

    // Initialize values for integration
    double value_temp ( 0. ), err ( 0. );
    double range_min[2] = { EgoSurf::range_points_u_.first, EgoSurf::range_points_v_.first };
    double range_max[2] = { EgoSurf::range_points_u_.second, EgoSurf::range_points_v_.second };

    if ( ( sigma_folded ( 0,0 ) < 1e-12 && sigma_folded ( 1,1 ) < 1e-12 ) || sigma_folded.determinant () < 1e-6 ) {

      if ( polygon_.empty () == true ) {
        hcubature ( 1, fkrig::AbsSurfPoint, ( surf_ptr[0] ).get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
        hcubature ( 1, fkrig::AbsSurfPoint, ( surf_ptr[1] ).get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value_temp, &err );
      } else {
        shared_ptr< std::pair< shared_ptr<Go::SplineSurface>, vector<Point> > > util_surf ( new std::pair< shared_ptr<Go::SplineSurface>, vector<Point> > ( surf_ptr[0], polygon_ ) );
        hcubature ( 1, fkrig::AbsSurfPointPoly, util_surf.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );

        util_surf.reset ( new std::pair< shared_ptr<Go::SplineSurface>, vector<Point> > ( surf_ptr[1], polygon_ ) );
        hcubature ( 1, fkrig::AbsSurfPointPoly, util_surf.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
      }

      value -= value_temp;

    } else if ( sigma_folded ( 0,0 ) < 1e-12 ) {

      if ( polygon_.empty () == true ) {
        hcubature ( 1, fkrig::AbsSurfPoint, ( surf_ptr[0] ).get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );

        shared_ptr< std::pair< shared_ptr< Go::SplineSurface >, double > > util ( new std::pair< shared_ptr< Go::SplineSurface >, double > ( surf_ptr[1], sigma_folded ( 1, 1 ) ) );

        hcubature ( 1, fkrig::EAbsSurfPoint, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value_temp, &err );
      } else {
        shared_ptr< std::pair< shared_ptr<Go::SplineSurface>, vector<Point> > > util_surf_0 ( new std::pair< shared_ptr<Go::SplineSurface>, vector<Point> > ( surf_ptr[0], polygon_ ) );
        hcubature ( 1, fkrig::AbsSurfPointPoly, util_surf_0.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );

        shared_ptr< std::tuple< shared_ptr<Go::SplineSurface>, double, vector<Point> > > util_surf_1 ( new std::tuple< shared_ptr<Go::SplineSurface>, double,  vector<Point> > ( surf_ptr[1], sigma_folded ( 1, 1 ), polygon_ ) );
        hcubature ( 1, fkrig::EAbsSurfPointPoly, util_surf_1.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
      }

      value -= value_temp;

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

      if ( polygon_.empty () == true ) {
        shared_ptr< std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> > > util ( new std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> > ( surf_ptr, llt_sigma_folded ) );

        hcubature ( 1, fkrig::MeanEiSurf, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
      } else {
        shared_ptr< std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> > > util ( new std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> > ( surf_ptr, llt_sigma_folded, polygon_ ) );

        hcubature ( 1, fkrig::MeanEiSurfPoly, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
      }
    }

  }

  return value;
}

//! Compute the value of the variance at the design coordinates coord
double
fkrig::EgoSurf::ComputeVariance ( RVectorXd coord ) const
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
    Go::SplineSurface surf = f_surf_->Predict ( coord );

    MatrixXd coord_ego ( 2, EgoBase::coord_ego_.cols() );

    coord_ego.row ( 0 ) = EgoBase::coord_ego_.row ( 0 );
    coord_ego.row ( 1 ) = coord.row ( 0 );

//   for ( size_t i = 0; i < coord.cols(); ++i )
//     EgoBase::coord_ego_(1,i) = coord(i);

    // Create vector of shared pointer to curves
    vector< shared_ptr<Go::SplineSurface> > surf_ptr;
    surf_ptr.resize ( 2 );
    ( surf_ptr[0] ).reset ( new Go::SplineSurface );
    ( surf_ptr[1] ).reset ( new Go::SplineSurface );
    surf_ptr[0] = Go::GeometryTools::surfaceSum ( surf_min_, 1, *nominal_surf_, -1 );
    surf_ptr[1] = Go::GeometryTools::surfaceSum ( surf, 1, *nominal_surf_, -1 );

    // Compute the covariance matrix between the points of coord (s1 = s2 = 1)
    MatrixXd sigma_folded ( 2, 2 );
    sigma_folded = f_surf_->PredictCovariance ( coord_ego );

    // Obtain area of the geometric domain
    double area = f_surf_->get_domain_range ();

    // Compute the covariance for each s
    sigma_folded /= area;

    // Initialize values for integration
    double err ( 0. );
    double range_min[2] = { EgoSurf::range_points_u_.first, EgoSurf::range_points_v_.first };
    double range_max[2] = { EgoSurf::range_points_u_.second, EgoSurf::range_points_v_.second };

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

      if ( polygon_.empty () == true ) {
        shared_ptr< std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> > > util ( new std::pair< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd> > ( surf_ptr, llt_sigma_folded ) );

        hcubature ( 1, fkrig::VarianceEiSurf, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
      } else {
         shared_ptr< std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> > > util ( new std::tuple< vector< shared_ptr< Go::SplineSurface > >, vector<MatrixXd>, vector<Point> > ( surf_ptr, llt_sigma_folded, polygon_ ) );

        hcubature ( 1, fkrig::VarianceEiSurfPoly, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );       
      }
      
//     } else if ( sigma_folded ( 1,1 ) >= 1e-12 && sigma_folded.determinant () >= 1e-6 ) {
    } else if ( sigma_folded ( 0,0 ) >= 1e-12 || sigma_folded ( 1,1 ) >= 1e-12 ) {

      double value_temp ( 0. );

      if ( polygon_.empty () == true ) {
        shared_ptr< std::pair< shared_ptr< Go::SplineSurface >, double > > util ( new std::pair< shared_ptr< Go::SplineSurface >, double > ( surf_ptr[0], sigma_folded ( 0, 0 ) ) );
        
        hcubature ( 1, fkrig::VarAbsSurfPoint, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
        
        util.reset ( new std::pair< shared_ptr< Go::SplineSurface >, double > ( surf_ptr[1], sigma_folded ( 1, 1 ) ) );
        
        hcubature ( 1, fkrig::VarAbsSurfPoint, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value_temp, &err );
      } else {
        shared_ptr< std::tuple< shared_ptr< Go::SplineSurface >, double, vector<Point> > > util ( new std::tuple< shared_ptr< Go::SplineSurface >, double, vector<Point> > ( surf_ptr[0], sigma_folded ( 0, 0 ), polygon_ ) );
        
        hcubature ( 1, fkrig::VarAbsSurfPointPoly, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );
        
        util.reset ( new std::tuple< shared_ptr< Go::SplineSurface>, double, vector<Point> > ( surf_ptr[1], sigma_folded ( 1, 1 ), polygon_ ) );
        
        hcubature ( 1, fkrig::VarAbsSurfPointPoly, util.get (), 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value_temp, &err );        
      }
      value += value_temp;
    }
  }

  return value;

}

//! Compute the expected improvment in location coord
double
fkrig::EgoSurf::ComputeFunction ( RVectorXd coord )
{
  return 0.;
}

//! Compute the L1 distance to the nominal function
void
fkrig::EgoSurf::ComputeMin ()
{
  // Obtain the matrix of coordinates
  MatrixXd coord = f_surf_->get_coord();

  // Check if it is a link model
  bool stop = false;
  shared_ptr<SurfBase> temp_ptr;
  MatrixXd temp_coord;
  temp_ptr = f_surf_->get_f_surf ();
  
  while ( stop == false ) {
    if ( ( temp_ptr != NULL ) == true ) {
      // Obtain the coordinates of the nested object
      temp_coord = temp_ptr->get_coord ();
      // Resize matrix
      coord.conservativeResize ( coord.rows () + temp_coord.rows (), Eigen::NoChange );
      // Bind columns
      coord.block ( coord.rows () - temp_coord.rows (), 0, temp_coord.rows (), 2 ) = temp_coord;
      // Obtain the nested object
      temp_ptr = temp_ptr->get_f_surf ();     
    } else {
      stop = true;
    }
  }  
  
  // Predict the curves in the design locations
  vector<Go::SplineSurface> surfaces = f_surf_->Predict ( coord );

  // Compute the expected value of the L1 distance between the predicted curves and the nominal curve
  vector<double> distance ( coord.rows (), 0. );
  for ( size_t i = 0; i < distance.size (); ++i )
    distance[i] = ComputeL1Mean ( surfaces[i], coord.row ( i ) );

  // Find the index of the min
  EgoBase::index_min_ = ( std::min_element ( distance.begin (), distance.end () ) - distance.begin () );

  // Save the geometric coordinate of the minimum
  EgoBase::coord_ego_ = coord.row ( EgoBase::index_min_ );
  // Save the predicted surface with miminum expected L1 distance
  surf_min_ = surfaces[EgoBase::index_min_];
}

//! Compute the range of the knots vector for the surfaces
void
fkrig::EgoSurf::ComputeUniqueRange ()
{
  // Obtain the range of the surface
  std::pair<double, double> range_u = f_surf_->get_range_u ();
  std::pair<double, double> range_v = f_surf_->get_range_v ();

  // Obtain the range of the nominal surface
  double min_par_u = nominal_surf_->startparam_u();
  double max_par_u = nominal_surf_->endparam_u();
  double min_par_v = nominal_surf_->startparam_v();
  double max_par_v = nominal_surf_->endparam_v();

  // The range is computed as the maximum of the minimum value and the minimum of the maximum values
  range_points_u_ = std::make_pair ( std::max ( range_u.first, min_par_u ), std::min ( range_u.second, max_par_u ) );
  range_points_v_ = std::make_pair ( std::max ( range_v.first, min_par_v ), std::min ( range_v.second, max_par_v ) );
}

//! Find the design coordinate that minimize the distance between the predicted surface and the nominal surface
void
fkrig::EgoSurf::ComputeMinDist ()
{

  // Create a global optimization object
  nlopt::opt opt_glob ( EgoBase::glob_alg_, EgoBase::lb_.size() );
  // Create a local optimization object
  nlopt::opt opt_loc ( EgoBase::loc_alg_, EgoBase::lb_.size() );

  nlopt::vfunc f = &fkrig::ObjectiveFunctionMinSurf;

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

//! Set the polygonal boundary
void
fkrig::EgoSurf::set_polygon ( vector<Point> polygon )
{

  polygon_ = polygon;

};

} // End of namespace
