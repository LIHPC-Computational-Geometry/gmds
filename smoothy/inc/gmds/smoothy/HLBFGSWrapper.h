/*----------------------------------------------------------------------------*/
#ifndef __HLBFGS_WRAPPER_H__
#define __HLBFGS_WRAPPER_H__
/*----------------------------------------------------------------------------*/
#include <cassert>
#include <functional>
#include <vector>
/*----------------------------------------------------------------------------*/
#include "HLBFGS.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
struct HLBFGSWrapper {
	typedef std::function<void(const std::vector<double> &x, double &f, std::vector<double> &g)> gradient_eval;

	HLBFGSWrapper(gradient_eval func) : gradient(func) {}
	void run(std::vector<double> &sol);

	const gradient_eval gradient;
	int maxiter = 10000;     // Maximum number of quasi-Newton updates
	double gtol = 1e-10;     // The iteration will stop when ||g||/max(1,||x||) <= gtol
	bool verbose = true;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     //__HLBFGS_WRAPPER_H__
/*----------------------------------------------------------------------------*/
