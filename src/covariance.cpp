#include "covariance.hpp"

namespace fkrig {

fkrig::Covariance::Covariance ( double phi,
                                double sigma,
                                int n,
                                cov_type type,
                                bool fix_nugget,
                                double nugget )
  : phi_ ( phi ), sigma_sq_ ( sigma ), n_ ( n ), type_ ( type ), fix_nugget_ ( fix_nugget ), nugget_ ( nugget )
{

  corr_matr_.resize ( n_, n_ );
  corr_matr_.diagonal().setOnes();

  cov_matr_.resize ( n_, n_ );

  corr_vec_.resize ( n_ * ( n_ - 1 ) / 2 );

  kappa_ = 0.;
}

// Constructor
fkrig::Covariance::Covariance ( int n,
                                cov_type type,
                                bool fix_nugget,
                                double nugget )
  : phi_ ( 0. ), sigma_sq_ ( 0. ), n_ ( n ), type_ ( type ), fix_nugget_ ( fix_nugget ), nugget_ ( nugget )
{

  corr_matr_.resize ( n_, n_ );
  corr_matr_.diagonal().setOnes();

  cov_matr_.resize ( n_, n_ );

  corr_vec_.resize ( n_ * ( n_ - 1 ) / 2 );

  kappa_ = 0.;
}

// Compute the covariance matrix
void
fkrig::Covariance::Eval ( VectorXd& u )
{
  // Check the size of the vector u
  if ( u.rows() != n_ * ( n_ - 1 ) / 2 ) {
    std::cerr << "Length of vector don't match the dimension of the matrix" << std::endl ;
    exit ( 1 );
  }

  // Compute correlation vector
  ComputeCorr ( u );

  // Fill the correlation matrix
  size_t count = 0;

  for ( size_t i = 0; i < corr_matr_.rows(); ++i )
    for ( size_t j = i + 1; j < corr_matr_.cols(); ++j, ++count )
      corr_matr_ ( i, j ) = corr_matr_ ( j, i ) = corr_vec_ ( count );

  // Fill the covariance matrix
  cov_matr_ = sigma_sq_ * corr_matr_;
  
  if ( nugget_ != 0. ) {
    for ( size_t i = 0; i < cov_matr_.rows(); ++i )
      cov_matr_(i, i) += nugget_;
  }

}

// Compute the correlation vector
void
fkrig::Covariance::ComputeCorr ( VectorXd& u )
{

  VectorXd u_phi = u / phi_;

  switch ( type_ ) {
  case 0: {
    for ( size_t i = 0; i < u_phi.rows(); ++i )
      u_phi[i] > 0. ? corr_vec_[i] = ( 1. + u_phi[i] + pow ( u_phi[i], 2 ) / 3. ) * exp ( - u_phi[i] ) : corr_vec_[i] = 1.;
//       u_phi[i] > 0. ? corr_vec_[i] = ( 1. + sqrt ( 5. ) * u_phi[i] + 5. / 3. * pow ( u_phi[i], 2 ) ) * exp ( - sqrt ( 5. ) * u_phi[i] ) : corr_vec_[i] = 1.;

    break;
  }
  case 1: {
    for ( size_t i = 0; i < u_phi.rows(); ++i )
      u_phi[i] > 0. ? corr_vec_[i] = ( 1. + u_phi[i] ) * exp ( - u_phi[i] ) : corr_vec_[i] = 1.;
//       u_phi[i] > 0. ? corr_vec_[i] = ( 1. + sqrt ( 3. ) * u_phi[i] ) * exp ( - sqrt ( 3. ) * u_phi[i] ) : corr_vec_[i] = 1.;

    break;
  }
  case 2: {
    for ( size_t i = 0; i < u.rows(); ++i )
      corr_vec_[i] = exp ( - pow ( u_phi[i] , kappa_ ) );

    break;
  }
  case 3: {
    for ( size_t i = 0; i < u.rows(); ++i )
      corr_vec_[i] = exp ( - u_phi[i] );

    break;
  }
  case 4: {
    for ( size_t i = 0; i < u.rows(); ++i )
      corr_vec_[i] = exp ( - pow ( u_phi[i], 2 ) ) ;

    break;
  }
  default:
    std::cerr << "Non supproted covariance type" << std::endl;
    exit ( 1 );
  }
}

// Compute the correalation at distance u
double
fkrig::Covariance::ComputeCorr ( const double& u ) const
{

  double value = 0.;
  double u_phi = u / phi_;

  switch ( type_ ) {
  case 0: {
    u_phi > 0. ? value = ( 1. + u_phi + pow ( u_phi, 2 ) / 3. ) * exp ( - u_phi ) : value = 1.;

    break;
  }
  case 1: {
    u_phi > 0. ? value = ( 1. + u_phi ) * exp ( - u_phi ) : value = 1.;

    break;
  }
  case 2: {
    value = exp ( - pow ( u_phi, kappa_ ) );

    break;
  }
  case 3: {
    value = exp ( - u_phi );

    break;
  }
  case 4: {
    value = exp ( - pow ( u_phi , 2 ) );

    break;
  }
  default:
    std::cerr << "Non supproted covariance type" << std::endl;
    exit ( 1 );
  }

  return value;
}

// Compute the correlation matrix ad distance u and rearrange the value in a row-based matrix n_row x n_col
MatrixXd
fkrig::Covariance::ComputeCorr ( VectorXd& u,
                                 int n_row,
                                 int n_col ) const
{

  // Check the size of vector u
  if ( u.rows() != n_row * n_col ) {
    std::cerr << "Length of vector don't match the dimension of the matrix" << std::endl ;
    exit ( 1 );
  }

  MatrixXd corr ( n_row, n_col );

  size_t count = 0;

  for ( size_t i = 0; i < corr.rows(); ++i )
    for ( size_t j = 0; j < corr.cols(); ++j, ++count )
      corr ( i, j ) = ComputeCorr ( u ( count ) );

  return corr;

}

// Compute the covariance matrix ad distance u and rearrange the value in a row-based matrix n_row x n_col
MatrixXd
fkrig::Covariance::ComputeCov ( VectorXd& u,
                                int n_row,
                                int n_col ) const
{
  MatrixXd cov = sigma_sq_ * ComputeCorr ( u, n_row, n_col );

  return cov;
}

//! Functor with MATERN_5_2 if fix_nugget is true
struct covariance_residual_matern_5_2_fix_nugget {
  covariance_residual_matern_5_2_fix_nugget ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, T* residual ) const {
    x_ > 0. ? residual[0] = T ( y_ ) - exp ( lsigma_sq[0] ) * ( 1. - ( 1. + x_ / exp ( lphi[0] ) + pow ( x_ / exp ( lphi[0] ), 2 ) / 3. ) * exp ( - x_ / exp ( lphi[0] ) ) ) : residual[0] = T ( y_ );
//     x_ > 0. ? residual[0] = T ( y_ ) - sigma_sq[0] * ( 1. - ( 1. + sqrt ( 5. ) * x_ / exp ( lphi[0] ) + 5. / 3. * pow ( x_ / exp ( lphi[0] ), 2 ) ) * exp ( - sqrt ( 5. ) * x_ / exp ( lphi[0] ) ) ) : residual[0] = T ( y_ ) ;
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with MATERN_5_2 if fix_nugget is false
struct covariance_residual_matern_5_2 {
  covariance_residual_matern_5_2 ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, const T* const nugget, T* residual ) const {
    x_ > 0. ? residual[0] = T ( y_ ) - nugget [0] - exp ( lsigma_sq[0] ) * ( 1. - ( 1. + x_ / exp ( lphi[0] ) + pow ( x_ / exp ( lphi[0] ), 2 ) / 3. ) * exp ( - x_ / exp ( lphi[0] ) ) ) : residual[0] = T ( y_ ) - nugget [0];
//     x_ > 0. ? residual[0] = T ( y_ ) - nugget[0] - sigma_sq[0] * ( 1. - ( 1. + sqrt ( 5. ) * x_ / exp ( lphi[0] ) + 5. / 3. * pow ( x_ / exp ( lphi[0] ), 2 ) ) * exp ( - sqrt ( 5. ) * x_ / exp ( lphi[0] ) ) ) : residual[0] = T ( y_ ) - nugget;
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with MATERN_3_2 if fix_nugget is true
struct covariance_residual_matern_3_2_fix_nugget {
  covariance_residual_matern_3_2_fix_nugget ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, T* residual ) const {
    x_ > 0. ? residual[0] = T ( y_ ) - exp ( lsigma_sq[0] ) * ( 1. - ( 1. +  x_ / exp ( lphi[0] ) ) * exp ( - x_ / exp ( lphi[0] ) ) ) : residual[0] = T ( y_ );
//     x_ > 0. ? residual[0] = T ( y_ ) - sigma_sq[0] * ( 1. - ( 1. + sqrt ( 3. ) * x_ / phi[0] ) * exp ( - sqrt ( 3. ) * x_ / phi[0] ) ) : residual[0] = T ( y_ ) ;
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with MATERN_3_2 if fix_nugget is false
struct covariance_residual_matern_3_2 {
  covariance_residual_matern_3_2 ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, const T* const nugget, T* residual ) const {
    x_ > 0. ? residual[0] = T ( y_ ) - nugget[0] - exp ( lsigma_sq[0] ) * ( 1. - ( 1. +  x_ / exp ( lphi[0] ) ) * exp ( - x_ / exp ( lphi[0] ) ) ) : residual[0] = T ( y_ ) - nugget[0];
//     x_ > 0. ? residual[0] = T ( y_ ) - nugget - sigma_sq[0] * ( 1. - ( 1. + sqrt ( 3. ) * x_ / phi[0] ) * exp ( - sqrt ( 3. ) * x_ / phi[0] ) ) : residual[0] = T ( y_ ) - nugget ;
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with POW_EXP if fix_nugget is true
struct covariance_residual_pow_exp_fix_nugget {
  covariance_residual_pow_exp_fix_nugget ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, const T* const tkappa, T* residual ) const {
    residual[0] = T ( y_ ) - exp ( lsigma_sq[0] ) * ( 1. - exp ( - pow ( x_ / exp ( lphi[0] ), 2. * exp ( tkappa[0] ) / ( 1. + exp ( tkappa[0] ) ) ) ) );
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with POW_EXP if fix_nugget is false
struct covariance_residual_pow_exp {
  covariance_residual_pow_exp ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, const T* const tkappa, const T* const nugget, T* residual ) const {
    residual[0] = T ( y_ ) - nugget[0] - exp ( lsigma_sq[0] ) * ( 1. - exp ( - pow ( x_ / exp ( lphi[0] ), 2. * exp ( tkappa[0] ) / ( 1. + exp ( tkappa[0] ) ) ) ) );
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with EXP if fix_nugget is true
struct covariance_residual_exp_fix_nugget {
  covariance_residual_exp_fix_nugget ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, T* residual ) const {
    residual[0] = T ( y_ ) - exp ( lsigma_sq[0] ) * ( 1. - exp ( - x_ / exp ( lphi[0] ) ) );
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with EXP if fix_nugget is false
struct covariance_residual_exp {
  covariance_residual_exp ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, const T* const nugget, T* residual ) const {
    residual[0] = T ( y_ ) - nugget[0] - exp ( lsigma_sq[0] ) * ( 1. - exp ( - x_ / exp ( lphi[0] ) ) );
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with GAUSS if fix_nugget is true
struct covariance_residual_gauss_fix_nugget {
  covariance_residual_gauss_fix_nugget ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, T* residual ) const {
    residual[0] = T ( y_ ) - exp ( lsigma_sq[0] ) * ( 1. - exp ( - pow ( x_ / exp ( lphi[0] ), 2 ) ) );
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

//! Functor with GAUSS if fix_nugget is false
struct covariance_residual_gauss {
  covariance_residual_gauss ( double x, double y )
    : x_ ( x ), y_ ( y ) {}

  template <typename T>
  bool operator() ( const T* const lsigma_sq, const T* const lphi, const T* const nugget, T* residual ) const {
    residual[0] = T ( y_ ) - nugget[0] - exp ( lsigma_sq[0] ) * ( 1. - exp ( - pow ( x_ / exp ( lphi[0] ), 2 ) ) );
    return true;
  }

private:
  // Observations for a sample.
  const double x_;
  const double y_;
};

// Estimate the parameter of the covariance function
void
fkrig::Covariance::Estimate ( VectorXd& u,
                              VectorXd& v )
{

  // Create optimization problem
  ceres::Problem problem;
  // Initial value for the optimization
  vector<double> init_sigma_sq = { v.maxCoeff () / 2., v.maxCoeff () * 3. / 4., v.maxCoeff () };
  vector<double> init_phi = { .4 / 5. * u.maxCoeff (), .8 / 5. * u.maxCoeff (), 1.6 / 5. * u.maxCoeff (), 2.4 / 5. * u.maxCoeff (), 3.2 / 5. * u.maxCoeff (), .8 * u.maxCoeff () };

  double lsigma_sq = 0., lphi = 0.;

  // Compute the value of the obliective function on a grid
  if ( fix_nugget_ == true ) {

    vector<double> f_0 ( init_sigma_sq.size() * init_phi.size(), 0. );

    vector<double>::iterator it = f_0.begin();

    for ( vector<double>::const_iterator it_j = init_phi.begin(); it_j != init_phi.end(); ++it_j, ++it ) {
      phi_ = *it_j;
      for ( vector<double>::const_iterator it_i = init_sigma_sq.begin(); it_i != init_sigma_sq.end(); ++it_i )
        for ( size_t i = 0; i < u.rows(); ++i )
          *it += pow ( v ( i ) - nugget_ + *it_i * ( 1. - ComputeCorr ( u ( i ) ) ), 2 );
    }

    // Find the minimum element of the vector
    vector<double>::iterator it_min = std::min_element ( f_0.begin(), f_0.end() );
    size_t it_min_i = std::distance ( f_0.begin(), it_min );

    int index_j = static_cast<int> ( it_min_i / init_sigma_sq.size() );
    int index_i = std::fmod ( it_min_i, init_sigma_sq.size() );

    lsigma_sq = log ( init_sigma_sq[index_i] );
    lphi = log ( init_phi[index_j] );

  } else {
    // Initial value for nugget
    vector<double> init_nugget = { v.maxCoeff () / 10., v.maxCoeff () / 4., v.maxCoeff () / 2. };

    vector<double> f_0 ( init_sigma_sq.size() * init_phi.size() * init_nugget.size(), 0. );

    vector<double>::iterator it = f_0.begin();

    for ( vector<double>::const_iterator it_j = init_phi.begin(); it_j != init_phi.end(); ++it_j ) {
      phi_ = *it_j;
      for ( vector<double>::const_iterator it_i = init_sigma_sq.begin(); it_i != init_sigma_sq.end(); ++it_i )
        for ( vector<double>::const_iterator it_k = init_nugget.begin(); it_k != init_nugget.end(); ++it_k, ++it )
          for ( size_t i = 0; i < u.rows(); ++i )
            *it += pow ( v ( i ) - *it_k + *it_i * ( 1. - ComputeCorr ( u ( i ) ) ), 2 );
    }

    // Find the minimum element of the vector
    vector<double>::iterator it_min = std::min_element ( f_0.begin(), f_0.end() );
    size_t it_min_i = std::distance ( f_0.begin(), it_min );

    int index_j = static_cast<int> ( it_min_i / init_sigma_sq.size() / init_nugget.size() );
    int index_i = std::fmod ( static_cast<int> ( it_min_i / init_nugget.size() ), init_sigma_sq.size() );
    int index_k = std::fmod ( it_min_i, init_nugget.size() );

    lsigma_sq = log ( init_sigma_sq[index_i] );
    lphi = log ( init_phi[index_j] );
    nugget_ = init_nugget[index_k];

  }

  // Set the optimization options
  ceres::Solver::Options options;
  // Maximum number of iteration
  options.max_num_iterations = 1e2;
  // Set linear solver type
  options.linear_solver_type = ceres::DENSE_QR;
  // Set no verbose output
  options.minimizer_progress_to_stdout = false;

  // Set optimization problem
  if ( fix_nugget_ == true ) {

    // Remove the nugget for the parameters estimation
    VectorXd temp ( v.rows() );
    for ( size_t i = 0; i < v.rows(); ++i )
      temp[i] = v[i] - nugget_;

    v = temp;

    switch ( type_ ) {
    case 0: {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_matern_5_2_fix_nugget, 1, 1, 1> ( new covariance_residual_matern_5_2_fix_nugget ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi );

      break;
    }
    case 1: {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_matern_3_2_fix_nugget, 1, 1, 1> ( new covariance_residual_matern_3_2_fix_nugget ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi );

      break;
    }
    case 2: {
      // Initialization of the kappa parameter to 1.5
      double tkappa = - log ( 1.5 / ( 2 - 1.5 ) );

      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_pow_exp_fix_nugget, 1, 1, 1, 1> ( new covariance_residual_pow_exp_fix_nugget ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi, &tkappa );

      // Set value of kappa_
      kappa_ = 2 * exp ( tkappa ) / ( 1 + exp ( tkappa ) );

      break;
    }
    case 3: {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_exp_fix_nugget, 1, 1, 1> ( new covariance_residual_exp_fix_nugget ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi );

      break;
    }
    case 4 : {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_gauss_fix_nugget, 1, 1, 1> ( new covariance_residual_gauss_fix_nugget ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi );

      break;
    }
    }
  } else {

    switch ( type_ ) {
    case 0: {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_matern_5_2, 1, 1, 1, 1> ( new covariance_residual_matern_5_2 ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi, &nugget_ );

      break;
    }
    case 1: {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_matern_3_2, 1, 1, 1, 1> ( new covariance_residual_matern_3_2 ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi, &nugget_ );

      break;
    }
    case 2: {
      // Initialization of the kappa parameter to 1.5
      double tkappa = log ( 2 / ( 2 - 1.5 ) );

      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_pow_exp, 1, 1, 1, 1, 1> ( new covariance_residual_pow_exp ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi, &tkappa, &nugget_ );

      // Set value of kappa_
      kappa_ = 2 * exp ( tkappa ) / ( 1 + exp ( tkappa ) );

      break;
    }
    case 3: {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_exp, 1, 1, 1, 1> ( new covariance_residual_exp ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi, &nugget_ );

      break;
    }
    case 4 : {
      for ( size_t i = 0; i < u.rows(); ++i )
        problem.AddResidualBlock ( new ceres::AutoDiffCostFunction<covariance_residual_gauss, 1, 1, 1, 1> ( new covariance_residual_gauss ( u[i], v[i] ) ), NULL, &lsigma_sq, &lphi, &nugget_ );

      break;
    }
    }

  }

  // Perform the optimization
  ceres::Solve ( options, &problem, &summary_ );

  // Save the parameters find
  sigma_sq_ = exp ( lsigma_sq );
  phi_ = exp ( lphi );

  // Compute the correlation and covariance matrix
  fkrig::Covariance::Eval ( u );

  if ( nugget_ < 0. ) {
    std::cout <<  "\33[0;31m" << "Nugget less than 0 -> setting nugget to 0 and re-run the algorithm\n" << "\33[0m" ;
    nugget_ = 0.;
    fix_nugget_ = true;
    Estimate ( u, v );
  }
//   std::cout << summary_.BriefReport() << "\n";
//   std::cout << "Final   sigma_sq: " << sigma_sq << " phi: " << exp ( lphi ) << "\n";

}

// template<>
// void
// fkrig::fkrig::Covariance::compute_corr <fkrig::MATERN_5_2> ( VectorXd u )
// {
//
//   VectorXd u_phi = u / phi_;
//
//   for ( size_t i = 0; i < u_phi.rows(); ++i )
//     u_phi[i] > 0 ? corr_vec_[i] = ( 1 + sqrt ( 5 ) * u_phi[i] + 5 / 3 * pow ( u_phi[i], 2 ) ) * exp ( - sqrt ( 5 ) * u_phi[i] ) : corr_vec_[i] = 1;
// //     u_phi[i] > 0 ? corr_vec_[i] = pow ( 2, ( - kappa_ + 1 ) ) / std::tgamma( kappa_ ) * pow ( u_phi[i], kappa_ ) * boost::math::cyl_bessel_k<double, double>( u_phi[i], kappa_ ) : corr_vec_[i] = 1;
//
// }
//
// template<>
// void
// fkrig::fkrig::Covariance::compute_corr <fkrig::MATER_3_2> ( VectorXd u )
// {
//
//   VectorXd u_phi = u / phi_;
//
//   for ( size_t i = 0; i < u_phi.rows(); ++i )
//     u_phi[i] > 0 ? corr_vec_[i] = ( 1 + sqrt ( 3 ) * u_phi[i] ) * exp ( sqrt ( 3 ) * u_phi[i] ) : corr_vec_[i] = 1;
// }
//
// template<>
// void
// fkrig::fkrig::Covariance::compute_corr <fkrig::POW_EXP> ( VectorXd u )
// {
//   VectorXd u_phi = u / phi_;
//
//   for ( size_t i = 0; i < u.rows(); ++i )
//     corr_vec_[i] = exp ( pow ( abs ( u_phi[i] ), kappa_ ) );
// }
//
// template<>
// void
// fkrig::fkrig::Covariance::compute_corr <fkrig::EXP> ( VectorXd u )
// {
//   VectorXd u_phi = u / phi_;
//
//   for ( size_t i = 0; i < u.rows(); ++i )
//     corr_vec_[i] = exp ( abs ( u_phi[i] ) );
// }
//
// template<>
// void
// fkrig::fkrig::Covariance::compute_corr <fkrig::GAUSS> ( VectorXd u )
// {
//   VectorXd u_phi = u / phi_;
//
//   for ( size_t i = 0; i < u.rows(); ++i )
//     corr_vec_[i] = exp ( pow ( abs ( u_phi[i] ), 2 ) );
// }

} // End of namespace
