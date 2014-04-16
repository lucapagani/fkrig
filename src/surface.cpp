#include "surface.hpp"
#include "util_fkrig.hpp"

using std::vector;
// using namespace fkrig::SurfBase;

namespace fkrig {

//! Constructor with parameterization as input
fkrig::Surf::Surf ( vector< vector<double> >& points,
                    vector< vector<double> >& param_u,
                    vector< vector<double> >& param_v,
                    const int dim,
                    MatrixXd coord,
                    model_type model,
                    double n_max,
                    double tol )
  : SurfBase ( points, param_u, param_v, dim, coord, model, n_max, tol ) {}


//! Compute the functional krigin parameters
void
fkrig::Surf::Eval()
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
  SurfBase::ComputeCovPar ();
  // Store the coefficients of the covariance matrix in a temp variable
  vector<double> temp_cov_par ( 4, 0. );
  temp_cov_par[0] = SurfBase::cov_->get_nugget(), temp_cov_par[1] = SurfBase::cov_->get_sigma_sq(), temp_cov_par[2] = SurfBase::cov_->get_phi(), temp_cov_par[3] = SurfBase::cov_->get_kappa();

  SurfBase::conv_ = false;
  SurfBase::n_iter_ = 1;

  while ( SurfBase::conv_ == false && SurfBase::n_iter_ < SurfBase::n_max_ ) {

    // Obtain the covariance matrix
    SurfBase::sigma_ = SurfBase::cov_->get_cov();

    // Compute the coefficients of the term a with GLS method
    ComputeAGls ();

    // Compute the coefficients of the term m with GLS method
    ComputeMGls ();

    // Compute the coefficients of the difference between the spline
    ComputeDiffCoefs ();

    // Copmute the parameters of the covariance matrix
    SurfBase::ComputeCovPar ();

    // Compute the absolute value of the difference between the coefficients
    vector<double> diff_coefs ( a_coefs_.rows() * a_coefs_.cols() + temp_cov_par.size(), 0. );

    size_t count = 0;

    // Difference between the coefficients of the curves
    for ( size_t i = 0; i < a_coefs_.rows(); ++i )
      for ( size_t j = 0; j < a_coefs_.cols(); ++j, ++count )
        diff_coefs[count] = std::abs ( a_coefs_ ( i, j ) - temp_a ( i, j ) );

    // Difference between tha values of the covariance parameters
    diff_coefs[count] = std::abs ( SurfBase::cov_->get_nugget() - temp_cov_par[0] ), diff_coefs[count + 1] = std::abs ( SurfBase::cov_->get_sigma_sq() - temp_cov_par[1] ), diff_coefs[count + 2] = std::abs ( SurfBase::cov_->get_phi() - temp_cov_par[2] ), diff_coefs[count + 3] = std::abs ( SurfBase::cov_->get_kappa() - temp_cov_par[3] );

    // Check convergence
    if ( std::any_of ( diff_coefs.begin(), diff_coefs.end(), [ this ] ( double i ) {
    return i > this->tol_;
  } ) ) {

      temp_a = a_coefs_;
      temp_cov_par[0] = SurfBase::cov_->get_nugget(), temp_cov_par[1] = SurfBase::cov_->get_sigma_sq(), temp_cov_par[2] = SurfBase::cov_->get_phi(), temp_cov_par[3] = SurfBase::cov_->get_kappa();

      ++SurfBase::n_iter_;

    }
    else {
      SurfBase::conv_ = true;
    }

  }

  // Obtain the covariance matrix at the last iteration
  SurfBase::sigma_ = SurfBase::cov_->get_cov ();

  // Compute the coefficients of the term a with GLS method
  ComputeAGls ();

  // Compute the coefficients of the term m with GLS method
  ComputeMGls ();

  // Compute the coefficients of the difference between the spline
  ComputeDiffCoefs ();

  // Compute the coefficients of inv(S) (z - F a)
  Eigen::LLT<MatrixXd> llt;
  llt.compute ( SurfBase::sigma_ );
  dx_coefs_.resize ( SurfBase::sigma_.rows(), m_diff_coefs_.cols() );
  dx_coefs_ = llt.solve ( m_diff_coefs_ );

}

// Predict the value at the design coordinates coord in the geometric value param
vector<vector<double> >
fkrig::Surf::Predict ( MatrixXd& coord,
                       vector<double>& param_u,
                       vector<double>& param_v ) const
{

  // Compute the predicted curves
  vector<Go::SplineSurface> surf = this->Predict ( coord );
  vector<vector<double> > points;
  points.resize ( surf.size() );

  // Compute the prediction at a parametric points
  Go::Point int_point;

  for ( size_t i = 0; i < points.size(); ++i ) {
    points[i].resize ( param_u.size() * dim_ );
    for ( size_t j = 0; j < param_u.size(); ++j ) {
      surf[i].point ( int_point, param_u[j], param_v[j] );
      points[i][j] = int_point[0];
    }
  }

  return points;
}

// Predict the spline curve at the design coordinates coord
vector<Go::SplineSurface>
fkrig::Surf::Predict ( MatrixXd& coord ) const
{

  // Check matrix coord
  if ( SurfBase::coord_.cols() != coord.cols() ) {
    std::cerr << "Matrix of coordinates must have the same number of columns" << std::endl;
  }

  // Compute the matrix of pariwise distance between coord and coord_
  VectorXd u ( SurfBase::coord_.rows() * coord.rows() );
  Eigen::Matrix<double, 1, Eigen::Dynamic> temp ( SurfBase::coord_.cols() );

  size_t count = 0;

  for ( size_t i = 0; i < coord.rows() ; ++i )
    for ( size_t j = 0; j < SurfBase::coord_.rows(); ++j, ++count ) {
      temp = coord.row ( i ) - SurfBase::coord_.row ( j );
      u ( count ) = std::sqrt ( temp * temp.adjoint() );
    }

  // Compute the covariance matrix between coord and coord_
  MatrixXd cov = SurfBase::cov_->ComputeCov ( u, coord.rows(), SurfBase::coord_.rows() );

  // Compute the coefficients of the curves of the term S0 inv(S) (z - F a)
  MatrixXd dx_coefs = cov * dx_coefs_;

  // Compute the coefficients of the curves of the term F_0 a
  MatrixXd F_0 = SurfBase::MakeModelMatrix ( coord ) ;

  // Obtain the vector of the coefficients of the surfaces and store them in a matrix
  MatrixXd surf_coefs ( SurfBase::surf_ptr_.size(), SurfBase::surf_ptr_[0]->numCoefs_u() * SurfBase::surf_ptr_[0]->numCoefs_v() );

  // Fill the matrix
  std::vector< double >::const_iterator i_coefs_begin, i_coefs_end;
  for ( size_t i = 0; i < SurfBase::surf_ptr_.size(); ++i ) {
    i_coefs_begin = SurfBase::surf_ptr_[i]->coefs_begin () + SurfBase::dim_ - 1;
    i_coefs_end = SurfBase::surf_ptr_[i]->coefs_end ();
    for ( size_t j = 0 ; i_coefs_begin < i_coefs_end; i_coefs_begin += SurfBase::dim_, ++j )
      surf_coefs ( i, j ) = *i_coefs_begin;
  }

  MatrixXd sx_coefs = F_0 * a_coefs_ * surf_coefs;

  // Compute the coefficients of the predicted surfaces
  vector< vector<double> > coefs;
  coefs.resize ( sx_coefs.rows() );

  for ( size_t i = 0; i < sx_coefs.rows(); ++i ) {
    coefs[i].resize ( sx_coefs.cols() );
    for ( size_t j = 0; j < sx_coefs.cols(); ++j )
      coefs[i][j] = sx_coefs ( i,j ) + dx_coefs ( i,j );
  }

  // Compute the predicted surfaces
  vector<Go::SplineSurface> pred;
  pred.resize ( coord.rows() );

  const Go::BsplineBasis basis_u = SurfBase::surf_ptr_[0]->basis_u(), basis_v = SurfBase::surf_ptr_[0]->basis_v();

  for ( size_t i = 0; i < pred.size(); ++i )
    pred[i] = Go::SplineSurface ( basis_u, basis_v, coefs[i].begin(), SurfBase::dim_ );

  return pred;
}

//! Predicted the total covarince matrix at the design coordinates coord
MatrixXd
fkrig::Surf::PredictCovariance ( MatrixXd& coord ) const
{

  // Check matrix coord
  if ( SurfBase::coord_.cols() != coord.cols() ) {
    std::cerr << "Matrix of coordinates must have the same number of columns" << std::endl;
  }

  // Compute the matrix of pariwise distance between the rows of coord
  VectorXd u ( coord.rows() * coord.rows() );
  Eigen::Matrix<double, 1, Eigen::Dynamic> temp ( coord.cols() );  
  
  size_t count = 0;

  for ( size_t i = 0; i < coord.rows() ; ++i )
    for ( size_t j = 0; j < SurfBase::coord_.rows(); ++j, ++count ) {
      temp = coord.row ( i ) - SurfBase::coord_.row ( j );
      u ( count ) = std::sqrt ( temp * temp.adjoint() );
    }  
  
  // Compute the covariance matrix between the rows of coord
  MatrixXd cov_0 = SurfBase::cov_->ComputeCov ( u, coord.rows(), coord.rows() );  
  
  // Compute the matrix of pariwise distance between coord and coord_
  u.resize( SurfBase::coord_.rows() * coord.rows() );
  temp.resize ( SurfBase::coord_.cols() );

  count = 0;

  for ( size_t i = 0; i < coord.rows() ; ++i )
    for ( size_t j = 0; j < SurfBase::coord_.rows(); ++j, ++count ) {
      temp = coord.row ( i ) - SurfBase::coord_.row ( j );
      u ( count ) = std::sqrt ( temp * temp.adjoint() );
    }

  // Compute the covariance matrix between coord and coord_
  MatrixXd cov_1 = SurfBase::cov_->ComputeCov ( u, coord.rows(), SurfBase::coord_.rows() );

  // Compute the coefficients of the curves of the term F_0 a
  MatrixXd F_0 = SurfBase::MakeModelMatrix ( coord ) ;  
  
  // Compute cholesky decomposition of sigma_
  Eigen::LLT<MatrixXd> llt_s, llt_f;
  llt_s.compute ( SurfBase::sigma_ );
 
  // Compute inv( simga_ ) F
  MatrixXd sf = llt_s.solve ( SurfBase::F_ );
  
  // Compute cholesky decomposition of F' sigma_ F
  llt_f.compute ( SurfBase::F_.transpose () * sf );
  
  // Compute F_0 - cov_1' inv( sigma_ ) F_
  MatrixXd temp_1 = F_0 - cov_1.transpose () * sf;
  
  MatrixXd cov = cov_0 - cov_1.transpose () * llt_s.solve ( cov_1 ) + temp_1 * llt_f.solve ( temp_1.transpose () );
  
  return cov;
}

//! Compute the matrix of coefficients of the term a with the Least Square method
void
fkrig::Surf::ComputeALs ()
{
  // Resize a_coefs_
  a_coefs_.resize ( SurfBase::F_.cols(), SurfBase::F_.rows() );

  if ( SurfBase::model_ == 0 ) {
    // Compute F'F
    double value = SurfBase::F_.col ( 0 ).adjoint () * SurfBase::F_.col ( 0 );
    // Compute inv ( F' * F ) * F'
    a_coefs_ = SurfBase::F_.transpose() / value ;
  } else {
    // Compute F'F
    MatrixXd F_F = SurfBase::F_.transpose() * SurfBase::F_;

    // Compute cholesky decomposition of F'F
    Eigen::LLT<MatrixXd> llt;
    llt.compute ( F_F ) ;

    // Compute inv ( F' * F ) * F'
    a_coefs_ = llt.solve ( SurfBase::F_.transpose () );
  }
}

// Compute the matrix of coefficients of the term a with the Generalized Least Square method
void
fkrig::Surf::ComputeAGls ()
{

  // Compute cholesky decomposition of sigma_
  Eigen::LLT<MatrixXd> llt;
  llt.compute ( SurfBase::sigma_ );
  // Compute inv(F) F
  MatrixXd S_F = llt.solve ( SurfBase::F_ );

  if ( SurfBase::model_ == 0 ) {
    // Compute F' inv(S) F
    double value = SurfBase::F_.col ( 0 ).adjoint () * S_F.col ( 0 );
    // Compute inv(F' inv(S) F) F' inv(S)
    a_coefs_ = S_F.transpose () / value;
  } else {
    // Compute F' inv(S) F
    MatrixXd F_S_F;
    F_S_F = SurfBase::F_.adjoint () * S_F;

    // Compute cholesky decomposition of F' inv(S) F
    Eigen::LLT<MatrixXd> llt2;
    llt2.compute ( F_S_F ) ;

    // Compute inv(F' inv(S) F) F' inv(S)
    a_coefs_ = llt2.solve ( S_F.transpose () );
  }

}

// Compute the matrix of coefficients of the term m with the Least Square method
void
fkrig::Surf::ComputeMLs ()
{
  // Compute F a
  m_coefs_ = SurfBase::F_ * a_coefs_;
}

// Compute the matrix of coefficients of the term m with the Generalized Least Square method
void
fkrig::Surf::ComputeMGls ()
{
  // Compute F a
  m_coefs_ = SurfBase::F_ * a_coefs_;
}

// Vector of coefficients of the difference between the curves and the splines of the term m
void
fkrig::Surf::ComputeDiffCoefs ()
{

  // Obtain the vector of the coefficients of the curves and store them in a matrix
  MatrixXd surf_coefs ( SurfBase::surf_ptr_.size(), SurfBase::surf_ptr_[0]->numCoefs_u() * SurfBase::surf_ptr_[0]->numCoefs_v() );
  MatrixXd m_surf_coefs ( surf_coefs.rows(), surf_coefs.cols() );

  // Fill the matrix
  std::vector< double >::const_iterator i_coefs_begin, i_coefs_end;
  for ( size_t i = 0; i < SurfBase::surf_ptr_.size(); ++i ) {
    i_coefs_begin = SurfBase::surf_ptr_[i]->coefs_begin () + SurfBase::dim_ - 1;
    i_coefs_end = SurfBase::surf_ptr_[i]->coefs_end ();
    for ( size_t j = 0 ; i_coefs_begin < i_coefs_end; i_coefs_begin += SurfBase::dim_, ++j )
      surf_coefs ( i, j ) = *i_coefs_begin;

  }

  // Compute the coefficients of the splines of the term m
  m_surf_coefs = m_coefs_ * surf_coefs;

  // Matrix of coefficients of the difference between the curves and the splines of the term m
  m_diff_coefs_ = surf_coefs - m_surf_coefs;

}

//! Compute the vector of the square of the parwise distances
void
fkrig::Surf::ComputeSqPairwiseDistances ()
{

  vector<double> diff_i_j ( m_diff_coefs_.cols(), 0. );
  Go::BsplineBasis basis_u = SurfBase::surf_ptr_[0]->basis_u (), basis_v = SurfBase::surf_ptr_[0]->basis_v ();
//   vector<double> knots_u = basis_u.getKnots(), knots_v = basis_v.getKnots();
//   vector<double>::const_iterator it_u_begin = knots_u.begin (), it_v_begin = knots_v.begin();
  Go::SplineSurface* surf ( new Go::SplineSurface );
  int it = 0;

  double value = 0., err = 0.;
  double range_min[2] = { SurfBase::range_points_u_.first, SurfBase::range_points_v_.first };
  double range_max[2] = { SurfBase::range_points_u_.second, SurfBase::range_points_v_.second };

  for ( size_t i = 0; i < m_diff_coefs_.rows() - 1; ++i )
    for ( size_t j = i + 1; j < m_diff_coefs_.rows(); ++j, ++it ) {
      // Compute the pairwise difference of the coefficients
      for ( size_t k = 0; k < diff_i_j.size(); ++k )
        diff_i_j[k] = m_diff_coefs_ ( i, k ) - m_diff_coefs_ ( j, k );

      // Create the curve of the pairwise difference
      *surf = Go::SplineSurface ( basis_u, basis_v, diff_i_j.begin(), SurfBase::dim_ );

      // Fill the vector of the differeces between the observed spline and the mean terms
      //TODO 2d integration in a polygon
      hcubature ( 1, fkrig::square_surface_point, surf, 2, range_min, range_max, 0, 0, 1e-4, ERROR_INDIVIDUAL, &value, &err );

      SurfBase::par_sq_dist_[it] = value;
    }

  delete surf;
}

} // end of namespace
