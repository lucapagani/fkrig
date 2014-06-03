#include "surface_base.hpp"
#include "util_fkrig.hpp"

using std::vector ;

namespace fkrig {

//! Constructor with parameterization as input
fkrig::SurfBase::SurfBase ( vector< vector<double> >& points,
                            vector< vector<double> >& param_u,
                            vector< vector<double> >& param_v,
                            int dim,
                            MatrixXd coord,
                            model_type model,
                            double n_max,
                            double tol )
  : points_ ( points ), param_u_ ( param_u ), param_v_ ( param_v ), dim_ ( dim ), coord_ ( coord ), model_ ( model ), n_max_ ( n_max ), tol_ ( tol ), polygon_ ()
{

  surf_ptr_.resize ( points_.size() );

  // Construct the model matrix
  MakeModelMatrix ();

  // Initialization vector dist_
  dist_.resize ( coord_.rows() * ( coord_.rows() - 1 ) / 2 );
  // Initialization vector par_sq_dist_
  par_sq_dist_.resize ( coord_.rows() * ( coord_.rows() - 1 ) / 2 );

  Eigen::Matrix<double, 1, Eigen::Dynamic> temp ( 1, coord_.cols() );
  int count = 0;

  // Compute pairwise distances between points in the coordinate space
  for ( size_t i = 0; i < coord_.rows() - 1; ++i )
    for ( size_t j = i + 1; j < coord_.rows(); ++j, ++count ) {
      temp = coord_.row ( i ) - coord_.row ( j );
      dist_ ( count ) = std::sqrt ( temp * temp.adjoint () );
    }

  // Make unique range
  MakeUniqueRange ();

}

//! Approximate the input data with MBA algorithm
void
fkrig::SurfBase::Approximate ( int h,
                               int m0,
                               int n0,
                               int smoothing_steps )
{
  // Compute the approximate surfaces
  for ( size_t i = 0; i < points_.size(); ++i )
    surf_ptr_[i] = ComputeMba ( points_[i], param_u_[i], param_v_[i], h, m0, n0, smoothing_steps );

  // Make the same knots vector
  Go::GeometryTools::unifySurfaceSplineSpace ( surf_ptr_, 1e-12 );
}


//TODO: compute the mse and choose the best
//! Approximate the input data with MBA algorithm
void
fkrig::SurfBase::Approximate ( std::pair<int, int> h_range,
                               int m0,
                               int n0,
                               int smoothing_iterations )
{
  Approximate ( h_range.first, m0, n0, smoothing_iterations );
}

// TODO: smoothing spline
//! Smoothing the input data
void
fkrig::SurfBase::Smoothing ( int num_coefs )
{

  Approximate();

}

//! Set covariance type
void
fkrig::SurfBase::set_covariance ( double phi,
                                  double sigma_sq,
                                  cov_type type,
                                  bool fix_nugget,
                                  double nugget )
{

  cov_.reset ( new fkrig::Covariance ( phi, sigma_sq, coord_.rows(), type, fix_nugget, nugget ) );

}

// Set covariance type
void
fkrig::SurfBase::set_covariance ( cov_type type,
                                  bool fix_nugget,
                                  double nugget )
{

  cov_.reset ( new fkrig::Covariance ( coord_.rows(), type, fix_nugget, nugget ) );

}

// Set covariance matrix
void
fkrig::SurfBase::set_covaraince ( std::unique_ptr<fkrig::Covariance> cov )
{

  *cov_ = *cov;

};

//! Set the polygonal boundary
void
fkrig::SurfBase::set_polygon ( vector<Point> polygon )
{

  polygon_ = polygon;

};

// Predict the value at the design coordinates coord in nu equally spaced geomtric points
vector<vector<double> >
fkrig::SurfBase::Predict ( MatrixXd& coord,
                           int n_u,
                           int n_v ) const
{

  // Create vector of n_u eqaully space value in the param space
  vector<double> param_u ( n_u * n_v, 0. ), param_v ( n_u * n_v, 0. );
  double delta_u = ( range_points_u_.second - range_points_u_.first ) / ( n_u - 1 );
  double delta_v = ( range_points_v_.second - range_points_v_.first ) / ( n_v - 1 );

  size_t count = 0;
  for ( size_t i = 0; i < n_u; ++i )
    for ( size_t j = 0; j < n_v; ++j, ++count ) {
      param_u[count] = range_points_u_.first + delta_u * i;
      param_v[count] = range_points_v_.first + delta_v * j;
    }

  // Predict the value
  vector< vector<double> > points = Predict ( coord, param_u, param_v );

  return points;

}

// Compute the model matrix
void
fkrig::SurfBase::MakeModelMatrix ()
{

  Eigen::FullPivLU<MatrixXd> lu_decomp;

  switch ( model_ ) {
  case 0: {
    F_.resize ( coord_.rows(), 1 );
    F_.setOnes();
    break;
  }
  case 1: {
    F_.resize ( coord_.rows(), coord_.cols() + 1 );
    F_.col ( 0 ).setOnes();

    // Linear terms
    for ( size_t i = 0; i < coord_.cols(); ++i )
      F_.col ( i + 1 ) = coord_.col ( i );

    // Check rank
    lu_decomp.compute ( F_ );
    if ( lu_decomp.rank() < F_.cols() ) {
      std::cerr << "The model matrix has not full rank" << std::endl ;
      exit ( 1 );
    }

    break;
  }
  case 2: {
    F_.resize ( coord_.rows(), 2 * coord_.cols() + coord_.cols() * ( coord_.cols() - 1 ) / 2 + 1 );
    F_.col ( 0 ).setOnes();

    // Linear terms
    for ( size_t i = 0; i < coord_.cols(); ++i )
      F_.col ( i + 1 ) = coord_.col ( i );

    // Cross terms
    size_t count_2 = coord_.cols() + 1;
    for ( size_t i = 0; i < coord_.cols() - 1; ++i )
      for ( size_t j = i + 1; j < coord_.cols(); ++j, ++count_2 ) {
        F_.col ( count_2 ) = coord_.col ( i ).cwiseProduct ( coord_.col ( j ) );
      }

    // Quadratic terms
    for ( size_t i = 0; i < coord_.cols(); ++i, ++count_2 )
      F_.col ( count_2 ) = coord_.col ( i ).cwiseProduct ( coord_.col ( i ) );

    // Check rank
    lu_decomp.compute ( F_ );
    if ( lu_decomp.rank() < F_.cols() ) {
      std::cerr << "The model matrix has not full rank" << std::endl ;
      exit ( 1 );
    }

    break;
  }
  case 3: {
    F_.resize ( coord_.rows(), coord_.cols() * 2 + 1 );
    F_.col ( 0 ).setOnes();

    // Linear terms
    for ( size_t i = 0; i < coord_.cols(); ++i )
      F_.col ( i + 1 ) = coord_.col ( i );

    // Quadratic terms
    size_t count_3 = coord_.cols() + 1;

    for ( size_t i = 0; i < coord_.cols(); ++i, ++count_3 )
      F_.col ( count_3 ) = coord_.col ( i ).cwiseProduct ( coord_.col ( i ) );

    // Check rank
    lu_decomp.compute ( F_ );
    if ( lu_decomp.rank() < F_.cols() ) {
      std::cerr << "The model matrix has not full rank" << std::endl ;
      exit ( 1 );
    }

    break;
  }
  default:
    std::cerr << "Non supproted model type" << std::endl;
    exit ( 1 );
  }
}

//! Compute the model matrix based on the matrix coord
MatrixXd
fkrig::SurfBase::MakeModelMatrix ( MatrixXd coord ) const
{

  MatrixXd F;

  switch ( model_ ) {
  case 0: {
    F.resize ( coord.rows(), 1 );
    F.setOnes();
    break;
  }
  case 1: {
    F.resize ( coord.rows(), coord.cols() + 1 );
    F.col ( 0 ).setOnes();

    // Linear terms
    for ( size_t i = 0; i < coord.cols(); ++i )
      F.col ( i + 1 ) = coord.col ( i );

    break;
  }
  case 2: {
    F.resize ( coord.rows(), 2 * coord.cols() + coord.cols() * ( coord.cols() - 1 ) / 2 + 1 );
    F.col ( 0 ).setOnes();

    // Linear terms
    for ( size_t i = 0; i < coord.cols(); ++i )
      F.col ( i + 1 ) = coord.col ( i );

    // Cross terms
    size_t count_2 = coord.cols() + 1;
    for ( size_t i = 0; i < coord.cols() - 1; ++i )
      for ( size_t j = i + 1; j < coord.cols(); ++j, ++count_2 ) {
        F.col ( count_2 ) = coord.col ( i ).cwiseProduct ( coord.col ( j ) );
      }

    // Quadratic terms
    for ( size_t i = 0; i < coord.cols(); ++i, ++count_2 )
      F.col ( count_2 ) = coord.col ( i ).cwiseProduct ( coord.col ( i ) );

    break;
  }
  case 3: {
    F.resize ( coord.rows(), coord.cols() * 2 + 1 );
    F.col ( 0 ).setOnes();

    // Linear terms
    for ( size_t i = 0; i < coord.cols(); ++i )
      F.col ( i + 1 ) = coord.col ( i );

    // Quadratic terms
    size_t count_3 = coord.cols() + 1;

    for ( size_t i = 0; i < coord.cols(); ++i, ++count_3 )
      F.col ( count_3 ) = coord.col ( i ).cwiseProduct ( coord.col ( i ) );

    break;
  }
  default:
    std::cerr << "Non supproted model type" << std::endl;
    exit ( 1 );
  }

  return F;
}

//! Make the same knots vector for the curves
void
fkrig::SurfBase::MakeUniqueRange ()
{
  // Compute range of the surface
  std::pair<vector<double>, vector<double>> min_element, max_element;
  min_element.first.resize ( param_u_.size(), 0 ), min_element.second.resize ( param_v_.size(), 0 );
  max_element.first.resize ( param_u_.size(), 0 ), max_element.second.resize ( param_v_.size(), 0 );

  for ( size_t i = 0; i < min_element.first.size(); ++i ) {
    min_element.first[i] = *std::min_element ( param_u_[i].begin(), param_u_[i].end() );
    max_element.first[i] = *std::max_element ( param_u_[i].begin(), param_u_[i].end() );
  }

  for ( size_t i = 0; i < min_element.second.size(); ++i ) {
    min_element.second[i] = *std::min_element ( param_v_[i].begin(), param_v_[i].end() );
    max_element.second[i] = *std::max_element ( param_v_[i].begin(), param_v_[i].end() );
  }

  // The range is computed as the maximum of the minimum value and the minimum of the maximum values
  range_points_u_ = std::make_pair ( *std::max_element ( min_element.first.begin(), min_element.first.end() ), *std::min_element ( max_element.first.begin(), max_element.first.end() ) );
  range_points_v_ = std::make_pair ( *std::max_element ( min_element.second.begin(), min_element.second.end() ), *std::min_element ( max_element.second.begin(), max_element.second.end() ) );

  // Check if the value of u and v are uniques
  if ( range_points_u_.first != *std::min_element ( min_element.first.begin(), min_element.first.end() ) || range_points_u_.second != *std::max_element ( max_element.first.begin(), max_element.first.end() ) || range_points_v_.first != *std::min_element ( min_element.second.begin(), min_element.second.end() ) || range_points_v_.second != *std::max_element ( max_element.second.begin(), max_element.second.end() ) ) {

    // Create an empty MBA
    std::unique_ptr<MBA> mba ( new MBA );
    UCBspl::SplineSurface surf;
    int pos;
    vector<double>::iterator it;
    boost::shared_ptr<dVec> u_arr ( new dVec ), v_arr ( new dVec ), z_arr ( new dVec );
    // Number of points for each direction
    size_t n_points_v = round ( ( range_points_v_.second - range_points_v_.first ) / dist_.minCoeff() );
    double dv = ( range_points_v_.second - range_points_v_.first ) / ( n_points_v - 1 );
    size_t n_points_u = round ( ( range_points_u_.second - range_points_u_.first ) / dist_.minCoeff() );
    double du = ( range_points_u_.second - range_points_u_.first ) / ( n_points_u - 1 );

    for ( size_t i = 0; i < points_.size(); ++i ) {

      if ( param_u_[i][0] != range_points_u_.first || * ( param_u_[i].end() - 1 ) != range_points_u_.second ) {
        // Define surface
        *u_arr = param_u_[i];
        *v_arr = param_v_[i];
        *z_arr = points_[i];

        mba.reset ( new MBA ( u_arr, v_arr, z_arr ) );

        // Compute approximation
        mba->MBAalg ( 1, 1, 10 );

        // Obtain surface object
        surf = mba->getSplineSurface();

        // Find the first value of the knots vector less or equal to the minimu value of the range
        it = std::lower_bound ( param_u_[i].begin(), param_u_[i].end(), range_points_u_.first );
        --it;
        pos = it - param_u_[i].begin ();
        // Delete the values out of range
        param_u_[i].erase ( param_u_[i].begin (), param_u_[i].begin () + pos + 1 );
        param_v_[i].erase ( param_v_[i].begin (), param_v_[i].begin () + pos + 1 );
        // Delete the points out the range
        points_[i].erase ( points_[i].begin(), points_[i].begin() + dim_ + pos * dim_ );

        for ( size_t i = 0; i < n_points_v; ++i ) {
          // Add the first param value
          param_u_[i].insert ( param_u_[i].begin (), range_points_u_.first );
          param_v_[i].insert ( param_v_[i].begin (), range_points_v_.first + i * dv );
          // Use the surface to predict the first value
          points_[i].insert ( points_[i].begin (), surf.f ( range_points_u_.first, range_points_v_.first + i * dv ) );
        }

        // Find the first value of the knots vector less or equal to the minimu value of the range
        it = std::upper_bound ( param_u_[i].begin(), param_u_[i].end(), range_points_u_.second );
        // Delete the values out of range
        pos = it - param_u_[i].begin ();
        param_u_[i].erase ( param_u_[i].begin() + pos, param_u_[i].end () );
        param_v_[i].erase ( param_v_[i].begin() + pos, param_u_[i].end () );
        // Delete the points out the range
        points_[i].erase ( points_[i].begin() + dim_ - 1 + pos * dim_, points_[i].end () );

        for ( size_t i = 0; i < n_points_v; ++i ) {
          // Add the first param value
          param_u_[i].push_back ( range_points_u_.second );
          param_v_[i].push_back ( range_points_v_.first + i * dv );
          // Use the surface to predict the first value
          points_[i].push_back ( surf.f ( range_points_u_.second, range_points_v_.first + i * dv ) );
        }

      }

      if ( param_v_[i][0] != range_points_v_.first || * ( param_v_[i].end() - 1 ) != range_points_v_.second ) {
        // Define surface
        *u_arr = param_u_[i];
        *v_arr = param_v_[i];
        *z_arr = points_[i];

        mba.reset ( new MBA ( u_arr, v_arr, z_arr ) );

        // Compute approximation
        mba->MBAalg ( 1, 1, 10 );

        // Obtain surface object
        surf = mba->getSplineSurface();

        // Find the first value of the knots vector less or equal to the minimu value of the range
        it = std::lower_bound ( param_v_[i].begin(), param_v_[i].end(), range_points_v_.first );
        --it;
        pos = it - param_v_[i].begin ();
        // Delete the values out of range
        param_u_[i].erase ( param_u_[i].begin (), param_u_[i].begin () + pos + 1 );
        param_v_[i].erase ( param_v_[i].begin (), param_v_[i].begin () + pos + 1 );
        // Delete the points out the range
        points_[i].erase ( points_[i].begin(), points_[i].begin() + dim_ + pos * dim_ );

        for ( size_t i = 0; i < n_points_u; ++i ) {
          // Add the first param value
          param_u_[i].insert ( param_u_[i].begin (), range_points_u_.first + i * du );
          param_v_[i].insert ( param_v_[i].begin (), range_points_v_.first );
          // Use the surface to predict the first value
          points_[i].insert ( points_[i].begin (), surf.f ( range_points_u_.first + i * du, range_points_v_.first ) );
        }

        // Find the first value of the knots vector less or equal to the minimu value of the range
        it = std::upper_bound ( param_v_[i].begin(), param_v_[i].end(), range_points_v_.second );
        // Delete the values out of range
        pos = it - param_v_[i].begin ();
        param_u_[i].erase ( param_u_[i].begin() + pos, param_u_[i].end () );
        param_v_[i].erase ( param_v_[i].begin() + pos, param_u_[i].end () );
        // Delete the points out the range
        points_[i].erase ( points_[i].begin() + dim_ - 1 + pos * dim_, points_[i].end () );

        for ( size_t i = 0; i < n_points_u; ++i ) {
          // Add the first param value
          param_u_[i].push_back ( range_points_u_.first + i * du );
          param_v_[i].push_back ( range_points_v_.second );
          // Use the surface to predict the first value
          points_[i].push_back ( surf.f ( range_points_u_.first + i * du, range_points_v_.second ) );
        }

      }

    }

  }


//   // Check if the range of the u value is unique
//   if ( range_points_u_.first != *std::min_element ( min_element.first.begin(), min_element.first.end() ) || range_points_u_.second != *std::max_element ( max_element.first.begin(), max_element.first.end() ) ) {
//     // Create an empty spline interpolator
//     std::unique_ptr<MBA> mba ( new MBA );
//     UCBspl::SplineSurface surf;
//     int pos;
//     vector<double>::iterator it;
//     boost::shared_ptr<dVec> u_arr ( new dVec ), v_arr ( new dVec ), z_arr ( new dVec );
//     // Number of points for each direction
//     size_t n_points_v = round ( ( range_points_v_.second - range_points_v_.first ) / dist_.minCoeff() );
//     double dv = ( range_points_v_.second - range_points_v_.first ) / ( n_points_v - 1 );
//
//     for ( size_t i = 0; i < points_.size(); ++i ) {
//
//       if ( param_u_[i][0] != range_points_u_.first || * ( param_u_[i].end() - 1 ) != range_points_u_.second ) {
//         // Define surface
//         *u_arr = param_u_[i];
//         *v_arr = param_v_[i];
//         *z_arr = points_[i];
//
//         mba.reset ( new MBA ( u_arr, v_arr, z_arr ) );
//
//         // Compute approximation
//         mba->MBAalg ( 1, 1, 10 );
//
//         // Obtain surface object
//         surf = mba->getSplineSurface();
//
//         // Find the first value of the knots vector less or equal to the minimu value of the range
//         it = std::lower_bound ( param_u_[i].begin(), param_u_[i].end(), range_points_u_.first );
//         --it;
//         pos = it - param_u_[i].begin ();
//         // Delete the values out of range
//         param_u_[i].erase ( param_u_[i].begin (), param_u_[i].begin () + pos + 1 );
//         param_v_[i].erase ( param_v_[i].begin (), param_v_[i].begin () + pos + 1 );
//         // Delete the points out the range
//         points_[i].erase ( points_[i].begin(), points_[i].begin() + dim_ + pos * dim_ );
//
//         for ( size_t i = 0; i < n_points_v; ++i ) {
//           // Add the first param value
//           param_u_[i].insert ( param_u_[i].begin (), range_points_u_.first );
//           param_v_[i].insert ( param_v_[i].begin (), range_points_v_.first + i * dv );
//           // Use the surface to predict the first value
//           points_[i].insert ( points_[i].begin (), surf.f ( range_points_u_.first, range_points_v_.first + i * dv ) );
//         }
//
//         // Find the first value of the knots vector less or equal to the minimu value of the range
//         it = std::upper_bound ( param_u_[i].begin(), param_u_[i].end(), range_points_u_.second );
//         // Delete the values out of range
//         pos = it - param_u_[i].begin ();
//         param_u_[i].erase ( param_u_[i].begin() + pos, param_u_[i].end () );
//         param_v_[i].erase ( param_v_[i].begin() + pos, param_u_[i].end () );
//         // Delete the points out the range
//         points_[i].erase ( points_[i].begin() + dim_ - 1 + pos * dim_, points_[i].end () );
//
//         for ( size_t i = 0; i < n_points_v; ++i ) {
//           // Add the first param value
//           param_u_[i].push_back ( range_points_u_.second );
//           param_v_[i].push_back ( range_points_v_.first + i * dv );
//           // Use the surface to predict the first value
//           points_[i].push_back ( surf.f ( range_points_u_.second, range_points_v_.first + i * dv ) );
//         }
//
//       }
//
//     }
//
//   }
//
//   // Check if the range of the u value is unique
//   if ( range_points_v_.first != *std::min_element ( min_element.second.begin(), min_element.second.end() ) || range_points_v_.second != *std::max_element ( max_element.second.begin(), max_element.second.end() ) ) {
//     // Create an empty spline interpolator
//     std::unique_ptr<MBA> mba ( new MBA );
//     UCBspl::SplineSurface surf;
//     int pos;
//     vector<double>::iterator it;
//     boost::shared_ptr<dVec> u_arr ( new dVec ), v_arr ( new dVec ), z_arr ( new dVec );
//     // Number of points for each direction
//     size_t n_points_u = round ( ( range_points_u_.second - range_points_u_.first ) / dist_.minCoeff() );
//     double du = ( range_points_u_.second - range_points_u_.first ) / ( n_points_u - 1 );
//
//     for ( size_t i = 0; i < points_.size(); ++i ) {
//
//       if ( param_v_[i][0] != range_points_v_.first || * ( param_v_[i].end() - 1 ) != range_points_v_.second ) {
//         // Define surface
//         *u_arr = param_u_[i];
//         *v_arr = param_v_[i];
//         *z_arr = points_[i];
//
//         mba.reset ( new MBA ( u_arr, v_arr, z_arr ) );
//
//         // Compute approximation
//         mba->MBAalg ( 1, 1, 10 );
//
//         // Obtain surface object
//         surf = mba->getSplineSurface();
//
//         // Find the first value of the knots vector less or equal to the minimu value of the range
//         it = std::lower_bound ( param_v_[i].begin(), param_v_[i].end(), range_points_v_.first );
//         --it;
//         pos = it - param_v_[i].begin ();
//         // Delete the values out of range
//         param_u_[i].erase ( param_u_[i].begin (), param_u_[i].begin () + pos + 1 );
//         param_v_[i].erase ( param_v_[i].begin (), param_v_[i].begin () + pos + 1 );
//         // Delete the points out the range
//         points_[i].erase ( points_[i].begin(), points_[i].begin() + dim_ + pos * dim_ );
//
//         for ( size_t i = 0; i < n_points_u; ++i ) {
//           // Add the first param value
//           param_u_[i].insert ( param_u_[i].begin (), range_points_u_.first + i * du );
//           param_v_[i].insert ( param_v_[i].begin (), range_points_v_.first );
//           // Use the surface to predict the first value
//           points_[i].insert ( points_[i].begin (), surf.f ( range_points_u_.first + i * du, range_points_v_.first ) );
//         }
//
//         // Find the first value of the knots vector less or equal to the minimu value of the range
//         it = std::upper_bound ( param_v_[i].begin(), param_v_[i].end(), range_points_v_.second );
//         // Delete the values out of range
//         pos = it - param_v_[i].begin ();
//         param_u_[i].erase ( param_u_[i].begin() + pos, param_u_[i].end () );
//         param_v_[i].erase ( param_v_[i].begin() + pos, param_u_[i].end () );
//         // Delete the points out the range
//         points_[i].erase ( points_[i].begin() + dim_ - 1 + pos * dim_, points_[i].end () );
//
//         for ( size_t i = 0; i < n_points_u; ++i ) {
//           // Add the first param value
//           param_u_[i].push_back ( range_points_u_.first + i * du );
//           param_v_[i].push_back ( range_points_v_.second );
//           // Use the surface to predict the first value
//           points_[i].push_back ( surf.f ( range_points_u_.first + i * du, range_points_v_.second ) );
//         }
//
//       }
//
//     }
//
//   }

}

//! Compute the parameters of the covariance matrix
void
fkrig::SurfBase::ComputeCovPar ()
{

  // Compute the vector or square pairwise distances
  ComputeSqPairwiseDistances ();

  // Estimate the parameters of the covariance matrix
  cov_->Estimate ( dist_, par_sq_dist_ );

}

//! Compute the MBA approximation and return a shared pointer to a SplineSurface
shared_ptr< Go::SplineSurface >
fkrig::SurfBase::ComputeMba ( vector<double> points,
                              vector<double> param_u,
                              vector<double> param_v,
                              int h,
                              int m0,
                              int n0,
                              int smoothing_iterations )
{

  // Create vector for MBA approximation
  boost::shared_ptr< vector<double> > u_arr ( new vector<double> ), v_arr ( new vector<double> ), z_arr ( new vector <double> );
  // Fill vector
  *u_arr = param_u;
  *v_arr = param_v;
  *z_arr = points;

  // Create an MBA object
  MBA mba ( u_arr, v_arr, z_arr );

  // Perform the approximation
  mba.MBAalg ( m0, n0, h, smoothing_iterations );

  UCBspl::SplineSurface surf = mba.getSplineSurface();

  // Obtain the matrix of coefficients
  boost::shared_ptr< GenMatrixType > phi ( new GenMatrixType );
  phi = mba.PHI();

  int noX = phi->noX(), noY = phi->noY();
//   double temp = sqrt ( param_u.size() * param_v.size() / m0 / n0 );
//   int noX = static_cast<int> ( m0 * temp ), noY = static_cast<int> ( n0 * temp );

  // Create an empty spline interpolator
  Go::SplineInterpolator interpol;
  interpol.setFreeConditions();

  // Create vector of parameters and predicted points on a regular grid noX x noY
  std::vector<double> u ( noX, 0. ), v ( noY, 0. ), z ( u.size() * v.size(), 0. );

  std::vector<double> knots_x ( noX, 0. );
  double umin = surf.umin(), umax = surf.umax();
  double dx = ( umax - umin ) / ( noX - 1 );

  std::vector<double> knots_y ( noY, 0. );
  double vmin = surf.vmin(), vmax = surf.vmax();
  double dy = ( vmax - vmin ) / ( noY - 1 );

  for ( size_t i = 0; i < u.size(); ++i )
    u[i] = umin + dx * i;

  for ( size_t i = 0; i < v.size(); ++i )
    v[i] = vmin + dy * i;

  int count_z = 0;

  for ( size_t i = 0; i < v.size(); ++i )
    for ( size_t j = 0; j < u.size(); ++j, ++count_z ) {
      z[count_z] = surf.f ( u[j], v[i] );
//       surf.derivatives( u[j], u[i], z[++count_z], z[++count_z] );
    }

  // Approximate the MBA surface with a SplineSurface
  shared_ptr<Go::SplineSurface> sisl_surf ( new Go::SplineSurface );
  sisl_surf->interpolate ( interpol, interpol, noX, noY, 1, &u[0], &v[0], &z[0] );

  return sisl_surf;
}


//! Return the area of the parametric space
double
fkrig::SurfBase::get_domain_range () const
{

  double value = 0.;

  if ( polygon_.empty() ) {
    value = ( range_points_u_.second - range_points_u_.first ) * ( range_points_v_.second - range_points_v_.first );
  } else {
    //  Public-domain function by Darel Rex Finley, 2006
    size_t j = polygon_.size () - 1;
    for ( size_t i = 0; i < polygon_.size (); ++i ) {
      value += ( polygon_[j].u + polygon_[i].u ) * ( polygon_[j].v - polygon_[i].v );
      j = i;
    }
    value *= -.5;
  }

  return value;
};

} // end of namespace
