/*----------------------------------------------------------------------------*/
#include <pybind11/pybind11.h>
/*----------------------------------------------------------------------------*/
namespace py = pybind11;
/*----------------------------------------------------------------------------*/
void bind_math(py::module &);
void bind_mesh(py::module &);
void bind_geometry(py::module &);
#if ENABLE_BLOCKING
void bind_blocking(py::module &);
#endif
/*----------------------------------------------------------------------------*/
// ig is the submodule name
PYBIND11_MODULE(gmds, m)
{
	// current version
	m.attr("__version__") = "dev";
	// module description
	m.doc() = "gmds is the python binding to a set of features provided by the gmds project";

	auto sub_math = m.def_submodule("math", "math functionalities");
	bind_math(sub_math);
	auto sub_mesh = m.def_submodule("mesh", "mesh kernel");
	bind_mesh(sub_mesh);
	auto sub_geom = m.def_submodule("cad", "cad interface");
	bind_geometry(sub_geom);
#if ENABLE_BLOCKING
	auto sub_block = m.def_submodule("blocking", "blocking kernel");
	bind_blocking(sub_block);
#endif
}
