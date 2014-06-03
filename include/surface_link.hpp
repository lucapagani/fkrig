#ifndef H_FKRIG_SURFACE_LINK__
#define H_FKRIG_SURFACE_LINK__

#include "surface.hpp"

namespace fkrig {

class SurfLink: public Surf {

public:

  //! Default constructor
  SurfLink () = default;

  //! Desctructor
  ~SurfLink() {};

  /*! @brief Constructor with parameterization as input
   *
   *  Constructor of the surface with points and parameterization of each point as input
   *
   *  @param points a vector with the vectors of the observed points
   *  @param param_u a vector with the vectors of the parameterization of the points in the u direction
   *  @param param_v a vector with the vectors of the parameterization of the points in the v direction
   *  @param dim the dimension of the points (only dim=1 supported)
   *  @param coord matrix of the coordinates of the design space
   *  @param model the type of the model matrix
   *  @param n_max number of maximum iterazions of the algorithm the estimate the functional kriging parameters
   *  @param tol tollerance of the diffenrence of the estimated parameters between two different steps of the algorithm
   */
  SurfLink ( vector< vector<double> >& points,
             vector< vector<double> >& param_u,
             vector< vector<double> >& param_v,
             shared_ptr<SurfBase> f_surf,
             const int dim,
             MatrixXd coord,
             model_type model = CONSTANT,
             double n_max = 1e3,
             double tol = 1e-6 );

  //! Compute the functional krigin parameters
  void
  Eval ();

  /*! Predict the value at the design coordinates coord in the geometric value param
   *  @param coord matrix of coordinates of the new design locations
   *  @param param_u vector of parameter values in the u direction for each new design location
   *  @param param_v vector of parameter values in the v direction for each new design location
   */
  vector<vector<double> >
  Predict ( MatrixXd& coord,
            vector<double>& param_u,
            vector<double>& param_v ) const;

  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord row vector of coordinates of the new desing locations
   */
  Go::SplineSurface
  Predict ( RVectorXd& coord ) const;

  /*! Predict the spline curve at the design coordinates coord
   *
   * @param coord matrix of coordinates of the new desing locations
   */
  vector<Go::SplineSurface>
  Predict ( MatrixXd& coord ) const;

  /*! Predicted the total covarince matrix at the design coordinates coord
   *
   *  @param coord matrix of coordinates of the new design locations
   */
  MatrixXd
  PredictCovariance ( MatrixXd& coord ) const;

};

} //! end of namespace
#endif
