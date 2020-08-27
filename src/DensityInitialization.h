#ifndef _DENSITY_INITIALIZATION_HPP_
#define _DENSITY_INITIALIZATION_HPP_

// This file contains the initialization procedure for the Density Estimation problem 

/*!  @brief An abstract base class dealing with the density initialization feature.
*/
template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DensityInitialization{
  protected:
    // A member to access data problem methods
    const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dataProblem_;

  public:
    //! A constructor.
    DensityInitialization(const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dp): dataProblem_(dp){};
    //! A destructor.
    virtual ~DensityInitialization(){};
    //! A pure virtual method to compute density initialization.
    virtual const VectorXr* chooseInitialization(Real lambda) const = 0;
};

/*!  @brief A class dealing with the user's density initialization.
*/
template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class UserInitialization : public DensityInitialization<Integrator, Integrator_noPoly, ORDER, mydim, ndim>{
  private:
    // A VectorXr which contains the user initialization
    VectorXr initialization;

  public:
    //! Delegating constructor
    UserInitialization(const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dp);
    //! An overridden method to compute density initialization when the user gives it.
    const VectorXr* chooseInitialization(Real lambda) const override;
};

/*!  @brief A class dealing with the density initialization given by a discretized heat diffusion process.
    It implements some methods useful to perform the discretized heat diffusion process and to choose the best initialization.
*/
template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class HeatProcess : public DensityInitialization<Integrator, Integrator_noPoly, ORDER, mydim, ndim>{
  private:
    // A member to access functional methods
    const FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& funcProblem_;
    // A vector of vectors containing all the possibile initial densities given by the heat diffusion process
    std::vector<VectorXr> init_proposals_;
    // patch_areas_: for each mesh node it saves the sum of the area of each triangle that has that node
    VectorXr patch_areas_;
    //  Parameters useful for the initialization of the density
    UInt niter_;
    Real alpha_;
    // Parameter to modify the initial density
    const Real epsilon_ = 1e-10;
    // Loglikelihood for each possible initial density
    VectorXr llik_;
    // Penalization term for each possible initial density
    VectorXr penTerm_;

    //! A method to compute the patch_areas_.
    void computePatchAreas();
    //! A method to compute the density exploting only the data.
    VectorXr computeDensityOnlyData();
    //! A method that provides a set of starting densities.
    void computeStartingDensities();

  public:
    //! A Constructor.
    HeatProcess(const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dp,
      const FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& fp);
    //! An overridden method to compute density initialization when it needes to be choosen among the proposals given by a discretized heat diffusion process.
    const VectorXr* chooseInitialization(Real lambda) const override;
};

#include "DensityInitialization_imp.h"
#endif
