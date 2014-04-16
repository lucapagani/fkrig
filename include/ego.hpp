#ifndef H_EGO__
#define H_EGO__

#include "surface.hpp"
#include "curve.hpp"

// typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RVectorXd;

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

  //! Destructor
  virtual ~Ego () {};

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

  //! Return the coordinate of the maximum of the expected improvment
  RVectorXd
  get_max_ei () {
    return ei_max_;
  };

protected:

  //! Index of the minimum
  size_t index_min_;
  //! Coordinates of the maximum of the expected inprovment
  RVectorXd ei_max_;

};

} // End of namespace
#endif
