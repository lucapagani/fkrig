#ifndef H_EGO__
#define H_EGO__

#include "surface.hpp"
// #include "curve.hpp"
#include <nlopt.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RVectorXd;

using std::vector;

namespace fkrig {

class EgoBase {

public:

  //! Default constructor
  EgoBase () = default;

//   /*! Contructor
//    *
//    *  @param chi functional variable
//    *  @param nominal nominal curve or surface
//    */
//   EgoBase ( std::unique_ptr<T_fkrig> chi,
//         std::unique_ptr<T_go> nominal );

  /*! Consrtuctor
   *
   *  @param lb vector with the lower bound of the optimization
   *  @param ub vector with the upper bound of the optimazation
   *  @param glob_alg global optimization algorithm
   *  @param loc_alg local optimization algorithm
   *  @param tol_glob relative tollerance for the global optimization algorithm
   *  @param tol_loc relative tollerance for the local optimization algorithm
   *  @param max_iter_glob maximum iteration for the global optimization algorithm
   *  @param max_iter_loc maximum iteration for the local optimization algorithm
   */
  EgoBase ( vector<double> lb,
            vector<double> ub,
//             unsigned int n,
            nlopt::algorithm glob_alg = nlopt::GN_DIRECT_L,
            nlopt::algorithm loc_alg = nlopt::LN_NELDERMEAD,
            double tol_glob = 1e-4,
            double tol_loc = 1e-5,
            double max_iter_glob = 1e3,
            double max_iter_loc = 1e3 );

  //! Destructor
  virtual ~EgoBase () {};

  /*! @brief Compute the L1 distance to the nominal funciton
   *
   *  Compute the L1 distance between the prediction of the i-th observed curve (surface) and the nominal curve (surface)
   *
   *  @param curve spline curve
   */
  virtual double
  ComputeL1 ( Go::SplineCurve& curve ) const = 0;

  /*! @brief Compute the expected value of the L1 distance to the nominal funciton
   *
   *  Compute the L1 distance of the expected value between the prediction of the curve (surface) at the geometric coordinate coord and the nominal curve (surface)
   *
   *  @param curve spline curve
   */
  virtual double
  ComputeL1Mean ( Go::SplineCurve& curve,
                  MatrixXd coord ) const = 0 ;

  /*! @brief Compute the L1 distance to the nominal funciton
   *
   * Compute the L1 distance between the prediction of the observed curves (surfaces) and the nominal curve (surface) and save the index of the minimum distance
   */
  virtual void
  ComputeMin () = 0;

  /*! Compute the expected improvment in location coord
   *
   *  @param coord row vector with a design coordinate
   */
  virtual double
  ComputeFunction ( RVectorXd coord ) = 0;

  //! Find the design coordinate that maximize the expected improvment
  void
  Compute ();

//   /*! Objective function for the maximization of the expected improvment
//    *
//    *  @param x value of the coordinates in the disegn space
//    *  @param grad value of the gradient
//    *  @param param vector with value of the mean and standard deviazion of the random variable
//    */
//   friend double
//   ObjectiveFunction ( const vector<double> &x,
//                       const vector<double> &grad,
//                       void* param );

  /*! Compute the value of the mean at the design coordinates coord
   *
   *  @param coord value of the coordinate in the design space
   */
  virtual double
  ComputeMean ( MatrixXd coord ) const = 0;

  /*! Compute the value of the variance at the design coordinates coord
   *
   *  @param coord value of the coordinate in the design space
   */
  virtual double
  ComputeVariance ( MatrixXd coord ) const = 0;
  
  //! Return the index of the closest curve (surface) to the nominal curve (surface)
  size_t
  get_index_min () {
    return index_min_;
  };  

  //! Return the coordinate of the maximum of the expected improvment
  RVectorXd
  get_max_ei () {
    return ei_max_;
  };

  //! Return the coord_ego_ matrix
  MatrixXd
  get_coord () {
    return coord_ego_;
  };

  //! Standard normal distribution
  boost::math::normal z_;

protected:

  //! Compute the range of the knots vector for the curves (surfaces)
  virtual void
  ComputeUniqueRange () = 0;

  //! Lower bound for the optimization
  vector<double> lb_;
  //! Upper bound for the optimization
  vector<double> ub_;
  //! Global optimization algorithm
  nlopt::algorithm glob_alg_;
  //! Local optimization algorithm
  nlopt::algorithm loc_alg_;
  //! Relative x tol for the global optimization
  double x_tol_glob_;
  //! Relative x tol for the local optimization
  double x_tol_loc_;
  //! Maximum number of iterazion allowed for the global optimization
  unsigned int max_iter_glob_;
  //! Maximum number of iterazion allowed for the local optimization
  unsigned int max_iter_loc_;
  //! Index of the minimum
  size_t index_min_;
  //! Coordinates of the maximum of the expected inprovment
  RVectorXd ei_max_;
  //! Object for the non linear optimization
  nlopt::opt opt_;
  //! Optimum value in the geometrical space with the global optimization algorithm
  nlopt::result result_;
//   //! Optimum value in the geometrical space with the local optimization algorithm
//   nlopt::result result_loc_;
  //! Optimum value of the objective function with the global optimization algorithm
  double value_;
//   //! Optimum value of the objective function with the local optimization algorithm
//   double value_local_;
  //! Matrix with the coordinates of the minimum (first row) and a new point (second row)
  MatrixXd coord_ego_;

};

} // End of namespace
#endif

