#ifndef H_EGO_SURFACE__
#define H_EGO_SURFACE__

#include "ego_base.hpp"
#include "surface.hpp"

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RVectorXd;

namespace fkrig {

class EgoSurf: public EgoBase {

public:

  //! Default constructor
  EgoSurf () = default;

  /*! Consrtuctor
   *
   *  @param f_surf functional kriging surface object
   *  @param nominal_surf nominal surface
   *  @param lb vector with the lower bound of the optimization
   *  @param ub vector with the upper bound of the optimazation
   *  @param glob_alg global optimization algorithm
   *  @param loc_alg local optimization algorithm
   *  @param tol_glob relative tollerance for the global optimization algorithm
   *  @param tol_loc relative tollerance for the local optimization algorithm
   *  @param max_iter_glob maximum iteration for the global optimization algorithm
   *  @param max_iter_loc maximum iteration for the local optimization algorithm
   */
  EgoSurf ( std::shared_ptr<fkrig::SurfBase>& f_surf,
            std::shared_ptr<Go::SplineSurface>& nominal_surf,
            vector<double> lb,
            vector<double> ub,
            nlopt::algorithm glob_alg = nlopt::GN_DIRECT_L,
            nlopt::algorithm loc_alg = nlopt::LN_NELDERMEAD,
            double tol_glob = 1e-4,
            double tol_loc = 1e-5,
            double max_iter_glob = 1e3,
            double max_iter_loc = 1e3 );

  //! Destructor
  ~EgoSurf () {};

  /*! @brief Compute the L1 distance to the nominal function
   *
   *  Compute the L1 distance between the prediction of the surface at the geometric coordinate coord and the nominal surface
   *
   *  @param surf spline surface
   */
  double
  ComputeL1 ( Go::SplineCurve& curve ) const {
    return 0.;
  };

  /*! @brief Compute the L1 distance to the nominal function
   *
   *  Compute the L1 distance between the prediction of the surface at the geometric coordinate coord and the nominal surface
   *
   *  @param surf spline surface
   */
  double
  ComputeL1 ( Go::SplineSurface& surf ) const;

  /*! @brief Compute the expected value of the L1 distance to the nominal function
   *
   *  Compute the L1 distance of the expected value between the prediction of the surface at the geometric coordinate coord and the nominal surface
   *
   *  @param surf spline surface
   */
  double
  ComputeL1Mean ( Go::SplineCurve& curve,
                  MatrixXd coord ) const {
    return 0.;
  };

  /*! @brief Compute the expected value of the L1 distance to the nominal function
   *
   *  Compute the L1 distance of the expected value between the prediction of the surface at the geometric coordinate coord and the nominal surface
   *
   *  @param surf spline surface
   */
  double
  ComputeL1Mean ( Go::SplineSurface& surf,
                  MatrixXd coord ) const;

  /*! @brief Compute the L1 distance to the nominal function
   *
   * Compute the L1 distance between the prediction of the observed surfaces and the nominal surface and save the index of the minimum distance
   */
  void
  ComputeMin ();

  /*! Compute the value of the mean at the design coordinates coord
   *
   *  @param coord value of the coordinate in the design space
   */
  double
  ComputeMean ( RVectorXd coord ) const;

  /*! Compute the value of the variance at the design coordinates coord
   *
   *  @param coord value of the coordinate in the design space
   */
  double
  ComputeVariance ( RVectorXd coord ) const;

  /*! Compute the expected improvment in location coord
   *
   *  @param coord row vector with a design coordinate
   */
  double
  ComputeFunction ( RVectorXd coord );
  
  //! Find the design coordinate that minimize the distance between the predicted surface and the nominal surface
  void
  ComputeMinDist ();

  //! Return fkrig object
  std::shared_ptr<SurfBase>
  get_fkrig () const {
    return f_surf_;
  };

//   /*! Objective function for the maximization of the expected improvment
//    *
//    *  @param x value of the coordinates in the disegn space
//    *  @param grad value of the gradient
//    *  @param param vector with value of the mean and standard deviazion of the random variable
//    */
//   static double
//   ObjectiveFunction ( const vector<double> &x,
//                       const vector<double> &grad,
//                       void* param );

private:

  //! Compute the range of the knots vector for the curves
  void
  ComputeUniqueRange ();

  //! Functional krigin object
  std::shared_ptr<SurfBase> f_surf_;
  //! Nominal curve
  std::shared_ptr<Go::SplineSurface> nominal_surf_;
  //! Range of the surfaces in the u direction
  std::pair<double, double> range_points_u_;
  //! Range of the surfaces in the v direction
  std::pair<double, double> range_points_v_;
  //! Curve with the mimum expected L1 distance
  Go::SplineSurface surf_min_;
//   //! Index of the minimum
//   size_t index_min_;

};

} // End of namespace
#endif
