#include "ego_base.hpp"
#include "util_fkrig.hpp"

// using boost::math::normal; // typedef provides default type is double.

namespace fkrig {

fkrig::EgoBase::EgoBase ( vector<double> lb,
                          vector<double> ub,
                          nlopt::algorithm glob_alg,
                          nlopt::algorithm loc_alg,
                          double tol_glob,
                          double tol_loc,
                          double max_iter_glob,
                          double max_iter_loc )
  : lb_ ( lb ), ub_ ( ub ), glob_alg_ ( glob_alg ), loc_alg_ ( loc_alg ), x_tol_glob_ ( tol_glob ), x_tol_loc_ ( tol_loc ), max_iter_glob_ ( max_iter_glob ), max_iter_loc_ ( max_iter_loc )
{
  if ( lb.size() != ub.size () ) {
    std::cerr << "Length of lower bound and upper bound must be the same." << std::endl;
    exit ( 1 );
  }

  // Resize matrix
  coord_ego_.resize ( 2, lb.size () );
}

//! Find the design coordinate that maximize the expected improvment
void
fkrig::EgoBase::Compute ()
{
  
  // Find the function with the minimum value
  ComputeMin();

  // Create a global optimization object
  nlopt::opt opt_glob ( glob_alg_, lb_.size() );
  // Create a local optimization object
  nlopt::opt opt_loc ( loc_alg_, lb_.size() );
  
  nlopt::vfunc f = &fkrig::ObjectiveFunction;
  
  // Set the objective function
  opt_glob.set_max_objective ( f, this );

  // Set the relative x tollerance
  opt_glob.set_xtol_rel ( x_tol_glob_ );
  opt_loc.set_xtol_rel (x_tol_loc_ );
  
  // Set the maximum number of iterations
  opt_glob.set_maxeval ( max_iter_glob_ );
  opt_loc.set_maxeval ( max_iter_loc_ );
  
  opt_glob.set_local_optimizer( opt_loc );
  
  // Chose a starting point
  std::vector<double> x0 ( lb_.size(), 0. );
  for ( size_t i = 0; i < x0.size(); ++i )
    x0[i] = ( lb_[i] + ub_[i] ) / 2; 

  // Preform the optimization
  result_ = opt_glob.optimize ( x0, value_ );

}

//! Compute the expected improvment in location coord
// double
// fkrig::EgoBase::ComputeFunction ( RVectorXd coord )
// {
//   // Compute mean and variance of random variable
//   double mean = ComputeMean ( coord_ego_ );
//   double sigma = std::sqrt ( ComputeVariance ( coord_ego_ ) );
//   double ratio = mean / sigma;
// 
//   // Compute the value of the expected improvment
//   double value = sigma * boost::math::pdf ( z_, ratio ) + mean * boost::math::cdf ( z_, ratio );
// 
//   return value;
// }

//! Objective function for the maximization of the expected improvment
// double
// fkrig::EgoBase::ObjectiveFunction ( const vector<double> &x,
//                                     const vector<double> &grad,
//                                     void* param )
// {
// 
//   // Fill the row of the second point
//   for ( size_t i = 0; i < x.size (); ++i )
//     coord_ego_ ( 1,i ) = x[i];
// 
//   // Compute mean and variance of random variable
//   double mean = ComputeMean ( coord_ego_ );
//   double sigma = std::sqrt ( ComputeVariance ( coord_ego_ ) );
//   double ratio = mean / sigma;
// 
//   // Compute the value of the expected improvment
//   double value = sigma * boost::math::pdf ( z_, ratio ) + mean * boost::math::cdf ( z_, ratio );
// 
//   return value;
// }

} // End of namespace
