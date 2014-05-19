#include "curve_base.hpp"
#include "util_fkrig.hpp"

using std::vector;

namespace fkrig {

// Constructor with parameterization as input
fkrig::CurveBase::CurveBase ( vector< vector<double> >& points,
                              vector< vector<double> >& param,
                              int dim,
                              MatrixXd coord,
                              model_type model,
                              double n_max,
                              double tol )
  : points_ ( points ), param_ ( param ), dim_ ( dim ), coord_ ( coord ), model_ ( model ), n_max_ ( n_max ), tol_ ( tol )
{

  curve_ptr_.reserve ( points_.size() );

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

// Interpolate the iput data
void
fkrig::CurveBase::Interpolate ()
{

  // Create an empty spline interpolator
  Go::SplineInterpolator interpol;
  // set free end condition
  interpol.setFreeConditions();
  // Create vector of empty curves
  vector< Go::SplineCurve > curve ;
  curve.resize ( points_.size() );

  // Compute the interpolation curves
  for ( size_t i = 0; i < points_.size(); ++i ) {
    curve[i].interpolate ( interpol, param_[i].size(), dim_, & ( param_[i] ) [0], & ( points_[i] ) [0] );
    curve_ptr_.push_back ( std::make_shared<Go::SplineCurve> ( curve[i] ) );
  }

  // Make the same knots vector
  Go::GeometryTools::unifyCurveSplineSpace ( curve_ptr_, 1e-12 );

}

// TODO: cv or similar to choose num_coefs
// Smoothing the input data
void
fkrig::CurveBase::Smoothing ( int num_coefs )
{
  // Create empty spline approximator
  Go::SplineApproximator approx;
  // Set number of coefficients (number of ponts of the control polygon)
  approx.setNumCoefs ( num_coefs );
  // Create vector of empty curves
  vector< Go::SplineCurve > curve ;
  curve.resize ( points_.size() );

  // Compute the interpolation curves
  for ( size_t i = 0; i < points_.size(); ++i ) {
    curve[i].interpolate ( approx, param_[i].size(), dim_, & ( param_[i] ) [0], & ( points_[i] ) [0] );
    curve_ptr_.push_back ( std::make_shared<Go::SplineCurve> ( curve[i] ) );
  }

  // Make the same knots vector
  Go::GeometryTools::unifyCurveSplineSpace ( curve_ptr_, 1e-12 );

}

// Set covariance type
void
fkrig::CurveBase::set_covariance ( double phi,
                                   double sigma_sq,
                                   cov_type type,
                                   bool fix_nugget,
                                   double nugget )
{

  cov_.reset ( new fkrig::Covariance ( phi, sigma_sq, coord_.rows(), type, fix_nugget, nugget ) );

}

// Set covariance type
void
fkrig::CurveBase::set_covariance ( cov_type type,
                                   bool fix_nugget,
                                   double nugget )
{

  cov_.reset ( new fkrig::Covariance ( coord_.rows(), type, fix_nugget, nugget ) );

}

// Set covariance matrix
void
fkrig::CurveBase::set_covaraince ( std::unique_ptr<fkrig::Covariance> cov )
{

  *cov_ = *cov;

};

// Predict the value at the design coordinates coord in nu equally spaced geomtric points
vector<vector<double> >
fkrig::CurveBase::Predict ( MatrixXd& coord,
                            int n_u ) const
{

  // Create vector of n_u eqaully space value in the param space
  vector<double> param ( n_u, 0. );
  double delta = ( range_points_.second - range_points_.first ) / ( n_u - 1 );

  param[0] = range_points_.first;
  for ( size_t i = 1; i < param.size(); ++i )
    param[i] = param[i - 1] + delta;

  // Predict the value
  vector< vector<double> > points = this->Predict ( coord, param );

  return points;

}

// Compute the model matrix
void
fkrig::CurveBase::MakeModelMatrix ()
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

// Compute the model matrix based on the matrix coord
MatrixXd
fkrig::CurveBase::MakeModelMatrix ( MatrixXd coord ) const
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

// Make the same knots vector for the curves
void
fkrig::CurveBase::MakeUniqueRange ()
{
  // Compute range of the curves
  vector<double> min_element ( param_.size(), 0. ), max_element ( param_.size(), 0. );

  for ( size_t i = 0; i < min_element.size(); ++i ) {
    min_element[i] = *std::min_element ( param_[i].begin(), param_[i].end() );
    max_element[i] = *std::max_element ( param_[i].begin(), param_[i].end() );
  }

  // The range is computed as the maximum of the minimum value and the minimum of the maximum values
  range_points_ = std::make_pair ( *std::max_element ( min_element.begin(), min_element.end() ), *std::min_element ( max_element.begin(), max_element.end() ) );

  // Check if min value is unique
  if ( range_points_.first != *std::min_element ( min_element.begin(), min_element.end() ) ) {
    // Create an empty spline interpolator
    Go::SplineInterpolator interpol;
    // set free end condition
    interpol.setFreeConditions();
    // Create vector of empty curves
    Go::SplineCurve curve ;
    // Create an empty point
    Go::Point point;
    double value = 0;
    int pos;
    vector<double>::iterator it;

    for ( size_t i = 0; i < points_.size(); ++i ) {

      if ( param_[i][0] != range_points_.first ) {
        // Compute the interpolation curves
        curve.interpolate ( interpol, param_[i].size(), dim_, & ( param_[i] ) [0], & ( points_[i] ) [0] );
        // Find the first value of the knots vector less or equal to the minimu value of the range
        it = std::lower_bound ( param_[i].begin(), param_[i].end(), range_points_.first );
        --it;
        pos = it - param_[i].begin ();
        // Delete the values out of range
        param_[i].erase ( param_[i].begin (), param_[i].begin () + pos + 1 );
        // Delete the points out the range
        points_[i].erase ( points_[i].begin(), points_[i].begin() + dim_ + pos * dim_ );
        // Add the first param value
        param_[i].insert ( param_[i].begin (), range_points_.first );
        // Use the curve to predict the first value
        curve.point ( point, range_points_.first );
        for ( int j = dim_ - 1; j >= 0; --j ) {
          value = point[j];
          // Add the values to points_
          points_[i].insert ( points_[i].begin (), value );
        }
      }

    }

  }

  // Check if max value is unique
  if ( range_points_.second != *std::max_element ( max_element.begin(), max_element.end() ) ) {
    // Create an empty spline interpolator
    Go::SplineInterpolator interpol;
    // set free end condition
    interpol.setFreeConditions();
    // Create vector of empty curves
    Go::SplineCurve curve ;
    // Create an empty point
    Go::Point point;
    double value;
    int pos;
    vector<double>::iterator it;

    for ( size_t i = 0; i < points_.size(); ++i ) {

      if ( * ( param_[i].end() - 1 ) != range_points_.second ) {
        // Compute the interpolation curves
        curve.interpolate ( interpol, param_[i].size(), dim_, & ( param_[i] ) [0], & ( points_[i] ) [0] );
        // Find the first value of the knots vector less or equal to the minimu value of the range
        it = std::upper_bound ( param_[i].begin(), param_[i].end(), range_points_.second );
        // Delete the values out of range
        pos = it - param_[i].begin ();
        param_[i].erase ( param_[i].begin() + pos, param_[i].end () );
        // Delete the points out the range
        points_[i].erase ( points_[i].begin() + dim_ - 1 + pos * dim_, points_[i].end () );
        // Add the first param value
        param_[i].push_back ( range_points_.second );
        // Use the curve to predict the first value
        curve.point ( point, range_points_.second );
        for ( size_t j = 0; j < dim_; ++j ) {
          value = point[j];
          // Add the values to points_
          points_[j].push_back ( value );
        }

      }

    }

  }

}

// Compute the parameters of the covariance matrix
void
fkrig::CurveBase::ComputeCovPar ()
{

  // Compute the vector or square pairwise distances
  ComputeSqPairwiseDistances ();

  // Estimate the parameters of the covariance matrix
  cov_->Estimate ( dist_, par_sq_dist_ );

}

} // end of namespace
