#include "curve.hpp"
#include "util_fkrig.hpp"

using std::vector;
// using namespace fkrig::CurveBase;

namespace fkrig {

// Constructor with implicit parameterization (only if the dimension of the geometric space is 2)
fkrig::Curve::Curve ( vector< vector<double> >& points,
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
fkrig::Curve::Curve ( vector< vector<double> >& points,
                      vector< vector<double> >& param,
                      const int dim,
                      MatrixXd coord,
                      model_type model,
                      double n_max,
                      double tol )
  : CurveBase ( points, param, dim, coord, model, n_max, tol ) {
    
    CurveBase::f_curve_ = NULL;
    
  }

// // Interpolate the iput data
// void
// fkrig::Curve::Interpolate ()
// {
//
//   // Create an empty spline interpolator
//   Go::SplineInterpolator interpol;
//   // set free end condition
//   interpol.setFreeConditions();
//   // Create vector of empty curves
//   vector< Go::SplineCurve > curve ;
//   curve.resize ( points_.size() );
//
//   // Compute the interpolation curves
//   for ( size_t i = 0; i < points_.size(); ++i ) {
//     curve[i].interpolate ( interpol, param_[i].size(), dim_, & ( param_[i] ) [0], & ( points_[i] ) [0] );
//     curve_ptr_.push_back ( std::make_shared<Go::SplineCurve> ( curve[i] ) );
//   }
//
//   // Make the same knots vector
//   Go::GeometryTools::unifyCurveSplineSpace ( curve_ptr_, 1e-12 );
//
// }

// // TODO: cv or similar to choose num_coefs
// // Smoothing the input data
// void
// fkrig::Curve::Smoothing ( int num_coefs )
// {
//   // Create empty spline approximator
//   Go::SplineApproximator approx;
//   // Set number of coefficients (number of ponts of the control polygon)
//   approx.setNumCoefs ( num_coefs );
//   // Create vector of empty curves
//   vector< Go::SplineCurve > curve ;
//   curve.resize ( points_.size() );
//
//   // Compute the interpolation curves
//   for ( size_t i = 0; i < points_.size(); ++i ) {
//     curve[i].interpolate ( approx, param_[i].size(), dim_, & ( param_[i] ) [0], & ( points_[i] ) [0] );
//     curve_ptr_.push_back ( std::make_shared<Go::SplineCurve> ( curve[i] ) );
//   }
//
//   // Make the same knots vector
//   Go::GeometryTools::unifyCurveSplineSpace ( curve_ptr_, 1e-12 );
//
// }

// Compute the functional krigin parameters
void
fkrig::Curve::Eval()
{
  // Compute the coefficients of the term a with LS method
  ComputeALs ();

  // Store the coefficients of the term a in a temp variable
  MatrixXd temp_a ( a_coefs_.rows(), a_coefs_.cols() );
  temp_a = a_coefs_;
  
  // Compute the coefficients of the term m with LS method
  ComputeMLs ();

  // Compute the coefficients of the difference between the spline
  ComputeDiffCoefs ();

  // Copmute the parameters of the covariance matrix
  CurveBase::ComputeCovPar ();
  // Store the coefficients of the covariance matrix in a temp variable
  vector<double> temp_cov_par ( 4, 0. );
  temp_cov_par[0] = CurveBase::cov_->get_nugget(), temp_cov_par[1] = CurveBase::cov_->get_sigma_sq(), temp_cov_par[2] = CurveBase::cov_->get_phi(), temp_cov_par[3] = CurveBase::cov_->get_kappa();

  CurveBase::conv_ = false;
  CurveBase::n_iter_ = 1;

  while ( CurveBase::conv_ == false && CurveBase::n_iter_ < CurveBase::n_max_ ) {

    // Obtain the covariance matrix
    CurveBase::sigma_ = CurveBase::cov_->get_cov();

    // Compute the coefficients of the term a with GLS method
    ComputeAGls ();

    // Compute the coefficients of the term m with GLS method
    ComputeMGls ();

    // Compute the coefficients of the difference between the spline
    ComputeDiffCoefs ();

    // Copmute the parameters of the covariance matrix
    CurveBase::ComputeCovPar ();

    // Compute the absolute value of the difference between the coefficients
    vector<double> diff_coefs ( a_coefs_.rows() * a_coefs_.cols() + temp_cov_par.size(), 0. );

    size_t count = 0;

    // Difference between the coefficients of the curves
    for ( size_t i = 0; i < a_coefs_.rows(); ++i )
      for ( size_t j = 0; j < a_coefs_.cols(); ++j, ++count )
        diff_coefs[count] = std::abs ( a_coefs_ ( i, j ) - temp_a ( i, j ) );

    // Difference between tha values of the covariance parameters
    diff_coefs[count] = std::abs ( CurveBase::cov_->get_nugget() - temp_cov_par[0] ), diff_coefs[count + 1] = std::abs ( CurveBase::cov_->get_sigma_sq() - temp_cov_par[1] ), diff_coefs[count + 2] = std::abs ( CurveBase::cov_->get_phi() - temp_cov_par[2] ), diff_coefs[count + 3] = std::abs ( CurveBase::cov_->get_kappa() - temp_cov_par[3] );

    // Check convergence
    if ( std::any_of ( diff_coefs.begin(), diff_coefs.end(), [ this ] ( double i ) {
    return i > this->tol_;
  } ) ) {

      temp_a = a_coefs_;
      temp_cov_par[0] = CurveBase::cov_->get_nugget(), temp_cov_par[1] = CurveBase::cov_->get_sigma_sq(), temp_cov_par[2] = CurveBase::cov_->get_phi(), temp_cov_par[3] = CurveBase::cov_->get_kappa();

      ++CurveBase::n_iter_;

    }
    else {
      CurveBase::conv_ = true;
    }

  }

  // Obtain the covariance matrix at the last iteration
  CurveBase::sigma_ = CurveBase::cov_->get_cov ();

  // Compute the coefficients of the term a with GLS method
  ComputeAGls ();

  // Compute the coefficients of the term m with GLS method
  ComputeMGls ();

  // Compute the coefficients of the difference between the spline
  ComputeDiffCoefs ();

  // Compute the coefficients of inv(S) (z - F a)
  Eigen::LLT<MatrixXd> llt;
  llt.compute ( CurveBase::sigma_ );
  dx_coefs_.resize ( CurveBase::sigma_.rows(), m_diff_coefs_.cols() );
  dx_coefs_ = llt.solve ( m_diff_coefs_ );

}

// TODO: prediction with dim_ > 1
// Predict the value at the design coordinates coord in the geometric value param
vector<vector<double> >
fkrig::Curve::Predict ( MatrixXd& coord,
                        vector<double>& param ) const
{

  // Compute the predicted curves
  vector<Go::SplineCurve> curve = this->Predict ( coord );
  vector<vector<double> > points;
  points.resize ( curve.size() );

  // Compute the prediction at a parametric points
  Go::Point int_point;

  for ( size_t i = 0; i < points.size(); ++i ) {
    points[i].resize ( param.size() * dim_ );
    for ( size_t j = 0; j < points[i].size(); ++j ) {
      curve[i].point ( int_point, param[j] );
      points[i][j] = int_point[0];
    }
  }

  return points;
}

// Predict the spline curve at the design coordinates coord
Go::SplineCurve
fkrig::Curve::Predict ( RVectorXd& coord ) const
{

  // Check matrix coord
  if ( CurveBase::coord_.cols() != coord.cols() ) {
    std::cerr << "Matrix of coordinates must have the same number of columns" << std::endl;
  }

  // Compute the matrix of pariwise distance between coord and coord_
  VectorXd u ( CurveBase::coord_.rows() );
  Eigen::Matrix<double, 1, Eigen::Dynamic> temp ( CurveBase::coord_.cols() );

  for ( size_t i = 0; i < CurveBase::coord_.rows(); ++i ) {
    temp = coord.row ( 0 ) - CurveBase::coord_.row ( i );
    u ( i ) = std::sqrt ( temp * temp.adjoint() );
  }

  // Compute the covariance matrix between coord and coord_
  MatrixXd cov = CurveBase::cov_->ComputeCov ( u, 1, CurveBase::coord_.rows() );

  // Compute the coefficients of the curves of the term S0 inv(S) (z - F a)
  MatrixXd dx_coefs = cov * dx_coefs_;

  // Compute the coefficients of the curves of the term F_0 a
  MatrixXd F_0 = CurveBase::MakeModelMatrix ( coord ) ;

  // Obtain the vector of the coefficients of the curves and store them in a matrix
  MatrixXd curve_coefs ( CurveBase::curve_ptr_.size(), CurveBase::curve_ptr_[0]->numCoefs() );

  // Fill the matrix
  std::vector< double >::const_iterator i_coefs_begin, i_coefs_end;
  for ( size_t i = 0; i < CurveBase::curve_ptr_.size(); ++i ) {
    i_coefs_begin = CurveBase::curve_ptr_[i]->coefs_begin () + CurveBase::dim_ - 1;
    i_coefs_end = CurveBase::curve_ptr_[i]->coefs_end ();
    for ( size_t j = 0 ; i_coefs_begin < i_coefs_end; i_coefs_begin += CurveBase::dim_, ++j )
      curve_coefs ( i, j ) = *i_coefs_begin;
  }

  MatrixXd sx_coefs = F_0 * a_coefs_ * curve_coefs;

  // Compute the coefficients of the predicted curve
  vector<double> coefs;

  coefs.resize ( sx_coefs.cols() );
  for ( size_t i = 0; i < sx_coefs.cols(); ++i )
    coefs[i] = sx_coefs ( i ) + dx_coefs ( i );

  // Compute the predicted curves
  Go::SplineCurve pred;

  size_t n_coefs = CurveBase::curve_ptr_[0]->numCoefs (), order = CurveBase::curve_ptr_[0]->order ();
  vector<double>::const_iterator it_begin = CurveBase::curve_ptr_[0]->knotsBegin ();

  pred = Go::SplineCurve ( n_coefs, order, it_begin, coefs.begin(), CurveBase::dim_ );

  return pred;
}

// Predict the spline curve at the design coordinates coord
vector<Go::SplineCurve>
fkrig::Curve::Predict ( MatrixXd& coord ) const
{

  // Check matrix coord
  if ( CurveBase::coord_.cols() != coord.cols() ) {
    std::cerr << "Matrix of coordinates must have the same number of columns" << std::endl;
  }

  // Compute the matrix of pariwise distance between coord and coord_
  VectorXd u ( CurveBase::coord_.rows() * coord.rows() );
  Eigen::Matrix<double, 1, Eigen::Dynamic> temp ( CurveBase::coord_.cols() );

  size_t count = 0;

  for ( size_t i = 0; i < coord.rows() ; ++i )
    for ( size_t j = 0; j < CurveBase::coord_.rows(); ++j, ++count ) {
      temp = coord.row ( i ) - CurveBase::coord_.row ( j );
      u ( count ) = std::sqrt ( temp * temp.adjoint() );
    }

  // Compute the covariance matrix between coord and coord_
  MatrixXd cov = CurveBase::cov_->ComputeCov ( u, coord.rows(), CurveBase::coord_.rows() );

  // Compute the coefficients of the curves of the term S0 inv(S) (z - F a)
  MatrixXd dx_coefs = cov * dx_coefs_;

  // Compute the coefficients of the curves of the term F_0 a
  MatrixXd F_0 = CurveBase::MakeModelMatrix ( coord ) ;

  // Obtain the vector of the coefficients of the curves and store them in a matrix
  MatrixXd curve_coefs ( CurveBase::curve_ptr_.size(), CurveBase::curve_ptr_[0]->numCoefs() );

  // Fill the matrix
  std::vector< double >::const_iterator i_coefs_begin, i_coefs_end;
  for ( size_t i = 0; i < CurveBase::curve_ptr_.size(); ++i ) {
    i_coefs_begin = CurveBase::curve_ptr_[i]->coefs_begin () + CurveBase::dim_ - 1;
    i_coefs_end = CurveBase::curve_ptr_[i]->coefs_end ();
    for ( size_t j = 0 ; i_coefs_begin < i_coefs_end; i_coefs_begin += CurveBase::dim_, ++j )
      curve_coefs ( i, j ) = *i_coefs_begin;
  }

  MatrixXd sx_coefs = F_0 * a_coefs_ * curve_coefs;

  // Compute the coefficients of the predicted curves
  vector< vector<double> > coefs;
  coefs.resize ( sx_coefs.rows() );

  for ( size_t i = 0; i < sx_coefs.rows(); ++i ) {
    coefs[i].resize ( sx_coefs.cols() );
    for ( size_t j = 0; j < sx_coefs.cols(); ++j )
      coefs[i][j] = sx_coefs ( i,j ) + dx_coefs ( i,j );
  }

  // Compute the predicted curves
  vector<Go::SplineCurve> pred;
  pred.resize ( coord.rows() );

  size_t n_coefs = CurveBase::curve_ptr_[0]->numCoefs (), order = CurveBase::curve_ptr_[0]->order ();
  vector<double>::const_iterator it_begin = CurveBase::curve_ptr_[0]->knotsBegin ();

  for ( size_t i = 0; i < pred.size(); ++i )
    pred[i] = Go::SplineCurve ( n_coefs, order, it_begin, coefs[i].begin(), CurveBase::dim_ );

  return pred;
}

//! Predicted the total covarince matrix at the design coordinates coord
MatrixXd
fkrig::Curve::PredictCovariance ( MatrixXd& coord ) const
{

  // Check matrix coord
  if ( CurveBase::coord_.cols() != coord.cols() ) {
    std::cerr << "Matrix of coordinates must have the same number of columns" << std::endl;
  }

  // Compute the matrix of pariwise distance between the rows of coord
  VectorXd u ( coord.rows() * coord.rows() );
  Eigen::Matrix<double, 1, Eigen::Dynamic> temp ( coord.cols() );

  size_t count = 0;

  for ( size_t i = 0; i < coord.rows() ; ++i )
    for ( size_t j = 0; j < coord.rows(); ++j, ++count ) {
      temp = coord.row ( i ) - coord.row ( j );
      u ( count ) = std::sqrt ( temp * temp.adjoint() );
    }

  // Compute the covariance matrix between the rows of coord
  MatrixXd cov_0 = CurveBase::cov_->ComputeCov ( u, coord.rows(), coord.rows() );

  // Compute the matrix of pariwise distance between coord and coord_
  u.resize ( CurveBase::coord_.rows() * coord.rows() );
  temp.resize ( CurveBase::coord_.cols() );

  count = 0;

  for ( size_t i = 0; i < coord.rows() ; ++i )
    for ( size_t j = 0; j < CurveBase::coord_.rows(); ++j, ++count ) {
      temp = coord.row ( i ) - CurveBase::coord_.row ( j );
      u ( count ) = std::sqrt ( temp * temp.adjoint() );
    }

  // Compute the covariance matrix between coord and coord_
  MatrixXd cov_1 = CurveBase::cov_->ComputeCov ( u, coord.rows(), CurveBase::coord_.rows() );

  // Compute the coefficients of the curves of the term F_0 a
  MatrixXd F_0 = CurveBase::MakeModelMatrix ( coord ) ;

  // Compute cholesky decomposition of sigma_
  Eigen::LLT<MatrixXd> llt_s, llt_f;
  llt_s.compute ( CurveBase::sigma_ );

  // Compute inv( simga_ ) F
  MatrixXd sf = llt_s.solve ( CurveBase::F_ );

  // Compute cholesky decomposition of F' inv( sigma_ ) F
  llt_f.compute ( CurveBase::F_.transpose () * sf );

  // Compute F_0 - cov_1' inv( sigma_ ) F_
//   MatrixXd temp_1 = F_0 - cov_1.transpose () * sf;
  MatrixXd temp_1 = F_0 - cov_1 * sf;

  MatrixXd cov = cov_0 - cov_1 * llt_s.solve ( cov_1.transpose () ); // + temp_1 * llt_f.solve ( temp_1.transpose () );

  return cov;
}

// Compute the matrix of coefficients of the term a with the Least Square method
void
fkrig::Curve::ComputeALs ()
{
  // Resize a_coefs_
  a_coefs_.resize ( CurveBase::F_.cols(), CurveBase::F_.rows() );

  if ( CurveBase::model_ == 0 ) {
    // Compute F'F
    double value = CurveBase::F_.col ( 0 ).adjoint () * CurveBase::F_.col ( 0 );
    // Compute inv ( F' * F ) * F'
    a_coefs_ = CurveBase::F_.transpose() / value ;
  } else {
    // Compute F'F
    MatrixXd F_F = CurveBase::F_.transpose() * CurveBase::F_;

    // Compute cholesky decomposition of F'F
    Eigen::LLT<MatrixXd> llt;
    llt.compute ( F_F ) ;

    // Compute inv ( F' * F ) * F'
    a_coefs_ = llt.solve ( CurveBase::F_.transpose () );
  }
}

// Compute the matrix of coefficients of the term a with the Generalized Least Square method
void
fkrig::Curve::ComputeAGls ()
{

  // Compute cholesky decomposition of sigma_
  Eigen::LLT<MatrixXd> llt;
  llt.compute ( CurveBase::sigma_ );
  // Compute inv(F) F
  MatrixXd S_F = llt.solve ( F_ );

  if ( CurveBase::model_ == 0 ) {
    // Compute F' inv(S) F
    double value = CurveBase::F_.col ( 0 ).adjoint () * S_F.col ( 0 );
    // Compute inv(F' inv(S) F) F' inv(S)
    a_coefs_ = S_F.transpose () / value;
  } else {
    // Compute F' inv(S) F
    MatrixXd F_S_F;
    F_S_F = CurveBase::F_.adjoint () * S_F;

    // Compute cholesky decomposition of F' inv(S) F
    Eigen::LLT<MatrixXd> llt2;
    llt2.compute ( F_S_F ) ;

    // Compute inv(F' inv(S) F) F' inv(S)
    a_coefs_ = llt2.solve ( S_F.transpose () );
  }

}

// Compute the matrix of coefficients of the term m with the Least Square method
void
fkrig::Curve::ComputeMLs ()
{
  // Compute F a
  m_coefs_ = CurveBase::F_ * a_coefs_;
}

// Compute the matrix of coefficients of the term m with the Generalized Least Square method
void
fkrig::Curve::ComputeMGls ()
{
  // Compute F a
  m_coefs_ = CurveBase::F_ * a_coefs_;
}

// Vector of coefficients of the difference between the curves and the splines of the term m
void
fkrig::Curve::ComputeDiffCoefs ()
{

  // Obtain the vector of the coefficients of the curves and store them in a matrix
  MatrixXd curve_coefs ( CurveBase::curve_ptr_.size(), CurveBase::curve_ptr_[0]->numCoefs() );
  MatrixXd m_curve_coefs ( curve_coefs.rows(), curve_coefs.cols() );

  // Fill the matrix
  std::vector< double >::const_iterator i_coefs_begin, i_coefs_end;
  for ( size_t i = 0; i < CurveBase::curve_ptr_.size(); ++i ) {
    i_coefs_begin = CurveBase::curve_ptr_[i]->coefs_begin () + CurveBase::dim_ - 1;
    i_coefs_end = CurveBase::curve_ptr_[i]->coefs_end ();
    for ( size_t j = 0 ; i_coefs_begin < i_coefs_end; i_coefs_begin += CurveBase::dim_, ++j )
      curve_coefs ( i, j ) = *i_coefs_begin;

  }

  // Compute the coefficients of the splines of the term m
  m_curve_coefs = m_coefs_ * curve_coefs;

  // Matrix of coefficients of the difference between the curves and the splines of the term m
  m_diff_coefs_ = curve_coefs - m_curve_coefs;

}

// Compute the vector of the square of the parwise distances
void
fkrig::Curve::ComputeSqPairwiseDistances ()
{

  vector<double> diff_i_j ( m_diff_coefs_.cols(), 0. );
  size_t n_coefs = CurveBase::curve_ptr_[0]->numCoefs (), order = CurveBase::curve_ptr_[0]->order ();
  vector<double>::const_iterator it_begin = CurveBase::curve_ptr_[0]->knotsBegin ();
  Go::SplineCurve curve;
//   shared_ptr<Go::SplineCurve> curve;
  int it = 0;

  double value = 0.;

  for ( size_t i = 0; i < m_diff_coefs_.rows() - 1; ++i )
    for ( size_t j = i + 1; j < m_diff_coefs_.rows(); ++j, ++it ) {
      // Compute the pairwise difference of the coefficients
      for ( size_t k = 0; k < diff_i_j.size(); ++k )
        diff_i_j[k] = m_diff_coefs_ ( i, k ) - m_diff_coefs_ ( j, k );

      // Create the curve of the pairwise difference
      curve = Go::SplineCurve ( n_coefs, order, it_begin, diff_i_j.begin(), CurveBase::dim_ );

//       // Create the curve of the pairwise difference
//       curve.reset ( new Go::SplineCurve ( n_coefs, order, it_begin, diff_i_j.begin(), CurveBase::dim_ ) );
//       
//       // Compute the difference between the curves and the difference between the mean term
//       curve = Go::GeometryTools::curveSum ( *( CurveBase::curve_ptr_[i] ), 1, *curve, 1 );
//       curve = Go::GeometryTools::curveSum ( *( CurveBase::curve_ptr_[j] ), -1, *curve, 1 );
      
      // Fill the vector of the differeces between the observed spline and the mean terms
      value = fkrig::AdaptiveSimpsons ( fkrig::SquareCurvePoint, curve, CurveBase::range_points_.first, CurveBase::range_points_.second, 1e-6, 10 );
      CurveBase::par_sq_dist_[it] = value;
    }

}

} // end of namespace
