/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/HLBFGSWrapper.h>
/*----------------------------------------------------------------------------*/
#if WIN32
// disable int to size_t warning
#	pragma warning(disable : 4267)
#endif
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
static const std::function<void(int N, double *x, double *prev_x, double *f, double *g)> *local_gradient = nullptr;
/*----------------------------------------------------------------------------*/
static void
static_gradient(int N, double *x, double *prev_x, double *f, double *g)
{
	(*local_gradient)(N, x, prev_x, f, g);
}
/*----------------------------------------------------------------------------*/
static void static_newiter_callback(int, int, double *, double *, double *, double *)
{}
/*----------------------------------------------------------------------------*/
void HLBFGSWrapper::run(std::vector<double> &sol)
{
	int hlbfgs_info[20] = {0};
	double parameter[20] = {0};
	INIT_HLBFGS(parameter, hlbfgs_info);
	hlbfgs_info[3] = 1;           // b_m1qn3_ ? 1 : 0; // determines whether we use m1qn3
	hlbfgs_info[4] = maxiter;     // max iterations
	hlbfgs_info[5] = verbose;     // verbose

	parameter[5] = gtol;
	parameter[6] = gtol;

	std::function<void(int N, double *x, double *prev_x, double *f, double *g)> wrapfunc = [&](int N, double *x, double *, double *f, double *g) {
		std::vector<double> array_x(N), array_g(N);
		for (int i = 0; i < N; i++)
			array_x[i] = x[i];

		gradient(array_x, *f, array_g);

		for (int i = 0; i < N; i++)
			g[i] = array_g[i];
	};

	local_gradient = &wrapfunc;
	HLBFGS((int) sol.size(), (int) 5, sol.data(), static_gradient, nullptr, HLBFGS_UPDATE_Hessian, static_newiter_callback, parameter, hlbfgs_info);
}
/*----------------------------------------------------------------------------*/
