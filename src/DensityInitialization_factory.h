#ifndef _DENS_INIT_FACTORY_HPP_
#define _DENS_INIT_FACTORY_HPP_

#include <memory>

//! brief@ A Factory class: a class for the choice of the density initialization.
template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DensityInitialization_factory
{
	public:
	//! A method that builds a pointer to the right object for the initialization choice.
	static std::unique_ptr<DensityInitialization<Integrator, Integrator_noPoly, ORDER,  mydim,  ndim>>
  createInitializationSolver(const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& fp){

      if(!dp.isFvecEmpty())
    		return make_unique<UserInitialization<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp);
    	else
    		return make_unique<HeatProcess<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp);

  }

};

#endif
