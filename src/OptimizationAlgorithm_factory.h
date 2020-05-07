#ifndef _OMP_ALGO_FACTORY_HPP
#define _OMP_ALGO_FACTORY_HPP

#include <memory>

//!brief@ A Factory class: a class for the choice of the step mehod for the optimization algorithm.
template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class MinimizationAlgorithm_factory
{
	public:
		//! A method that builds a pointer to the right object for the step choice, taking as parameters a string and others objects needed for constructor.
	static std::shared_ptr<MinimizationAlgorithm<Integrator, Integrator_noPoly, ORDER,  mydim,  ndim>>
  createStepSolver(const DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& dp,
    const FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>& fp,
		const std::string& d, const std::string& s)
	{
		if(s == "Fixed_Step") return std::make_shared<FixedStep<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, d);

    else if(s == "Backtracking_Method") return std::make_shared<BacktrackingMethod<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, d);

    else if(s == "Wolfe_Method") return std::make_shared<WolfeMethod<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp, d);

		else{

      #ifdef R_VERSION_
      Rprintf("Unknown step option - using fixed step\n");
      #else
      std::cout<<"Unknown step option - using fixed step"<<std::endl;
      #endif

			return std::make_shared<FixedStep<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dp, fp,  std::move(d));
		}
  }

};

#endif
