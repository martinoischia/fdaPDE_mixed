#ifndef _PREPROCESS_FACTORY_HPP_
#define _PREPROCESS_FACTORY_HPP_

#include <memory>

//! @brief A Factory class: a class for the choice of the cross-validation method.
template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class Preprocess_factory
{
	public:
	//! A method that builds a pointer to the right object for the cross-validation method choice, taking as parameters a string and others objects needed for constructor.
	static std::unique_ptr<Preprocess<Integrator, Integrator_noPoly, ORDER,  mydim,  ndim>>
  createPreprocessSolver(const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& fp,
    std::shared_ptr<MinimizationAlgorithm<Integrator, Integrator_noPoly, ORDER, mydim, ndim>> ma, const std::string& p){

			if(p=="RightCV")
				return make_unique<RightCrossValidation<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, ma);
			else if(p=="SimplifiedCV")
				return make_unique<SimplifiedCrossValidation<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, ma);
			else if(p=="NoCrossValidation")
      	return make_unique<NoCrossValidation<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp);
			else
				return make_unique<RightCrossValidation<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, ma);

  }

};

#endif
