/*----------------------------------------------------------------------------*/
#include "gmds/cadfac/FACManager.h"
#include "gmds/cadfac/FACPoint.h"
#include "gmds/utils/CommonTypes.h"
/*----------------------------------------------------------------------------*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/*----------------------------------------------------------------------------*/
namespace py = pybind11;
/*----------------------------------------------------------------------------*/
void bind_geometry(py::module &m){
	py::enum_<gmds::cad::GeomCurve::CurvatureInfo>(m, "CurvatureInfo")
	   .value("CURVATURE_FLAT", gmds::cad::GeomCurve::Flat)
	   .value("CURVATURE_CONVEX", gmds::cad::GeomCurve::Convex)
	   .value("CURVATURE_CONCAVE", gmds::cad::GeomCurve::Concave)
	   .value("CURVATURE_SMOOTH_CONVEX", gmds::cad::GeomCurve::SmoothConvex)
	   .value("CURVATURE_SMOOTH_CONCAV", gmds::cad::GeomCurve::SmoothConcave)
	   .export_values();

	py::class_<gmds::cad::GeomManager>(m, "GeomManager");

	py::class_<gmds::cad::FACManager, gmds::cad::GeomManager>(m, "FACManager")
	   .def(py::init<>())
	   .def("get_nb_volumes", &gmds::cad::FACManager::getNbVolumes)
	   .def("get_nb_surfaces", &gmds::cad::FACManager::getNbSurfaces)
	   .def("get_nb_curves", &gmds::cad::FACManager::getNbCurves)
	   .def("get_nb_points", &gmds::cad::FACManager::getNbPoints)
	   .def("get_volumes", static_cast<std::vector<gmds::cad::GeomVolume *> (gmds::cad::FACManager::*)() const>(&gmds::cad::FACManager::getVolumes),
	        py::return_value_policy::reference)
	   .def("get_surfaces", static_cast<std::vector<gmds::cad::GeomSurface *> (gmds::cad::FACManager::*)() const>(&gmds::cad::FACManager::getSurfaces),
	        py::return_value_policy::reference)
	   .def("get_curves", static_cast<std::vector<gmds::cad::GeomCurve *> (gmds::cad::FACManager::*)() const>(&gmds::cad::FACManager::getCurves),
	        py::return_value_policy::reference)
	   .def("get_points", static_cast<std::vector<gmds::cad::GeomPoint *> (gmds::cad::FACManager::*)() const>(&gmds::cad::FACManager::getPoints),
	        py::return_value_policy::reference)
	   .def("get_volume", &gmds::cad::FACManager::getVolume)
	   .def("get_surface", &gmds::cad::FACManager::getSurface)
	   .def("get_curve", &gmds::cad::FACManager::getCurve)
	   .def("get_point", &gmds::cad::FACManager::getPoint)
	   .def("init_from_3d_mesh", &gmds::cad::FACManager::initFrom3DMesh);

	py::class_<gmds::cad::GeomPoint>(m, "GeomPoint");
	py::class_<gmds::cad::FACPoint, gmds::cad::GeomPoint>(m, "FACPoint")
	   .def("x", &gmds::cad::FACPoint::X)
	   .def("y", &gmds::cad::FACPoint::Y)
	   .def("z", &gmds::cad::FACPoint::Z)
	   .def("id", &gmds::cad::FACPoint::id)
	   .def("bbox", &gmds::cad::FACPoint::BBox)
	   .def("curves", &gmds::cad::FACPoint::curves, py::return_value_policy::reference)
	   .def("surfaces", &gmds::cad::FACPoint::surfaces, py::return_value_policy::reference)
	   .def("volumes", &gmds::cad::FACPoint::volumes, py::return_value_policy::reference);

	py::class_<gmds::cad::GeomCurve>(m, "GeomCurve");
	py::class_<gmds::cad::FACCurve, gmds::cad::GeomCurve>(m, "FACCurve")
	   .def("id", &gmds::cad::FACCurve::id)
	   .def("length", &gmds::cad::FACCurve::length)
	   .def("tangent", &gmds::cad::FACCurve::computeTangent)
	   .def("closest_point", &gmds::cad::FACCurve::closestPoint)
	   .def("bbox", &gmds::cad::FACCurve::BBox)
	   .def("is_loop", &gmds::cad::FACCurve::isALoop)
	   .def("curvature_info", &gmds::cad::FACCurve::getCurvatureInfo)
	   .def("dihedral_angle", &gmds::cad::FACCurve::computeDihedralAngle)
	   .def("points", &gmds::cad::FACCurve::points, py::return_value_policy::reference)
	   .def("surfaces", &gmds::cad::FACCurve::surfaces, py::return_value_policy::reference)
	   .def("volumes", &gmds::cad::FACCurve::volumes, py::return_value_policy::reference);

	py::class_<gmds::cad::GeomSurface>(m, "GeomSurface");
	py::class_<gmds::cad::FACSurface, gmds::cad::GeomSurface>(m, "FACSurface")
	   .def("id", &gmds::cad::FACSurface::id)
	   .def("area", &gmds::cad::FACSurface::computeArea)
	   .def("center", &gmds::cad::FACSurface::center)
	   .def("closest_point", &gmds::cad::FACSurface::closestPoint)
	   .def("bbox", &gmds::cad::FACSurface::BBox)
	   .def("points", &gmds::cad::FACSurface::points, py::return_value_policy::reference)
	   .def("curves", &gmds::cad::FACSurface::curves, py::return_value_policy::reference)
	   .def("volumes", &gmds::cad::FACSurface::volumes, py::return_value_policy::reference);

	py::class_<gmds::cad::GeomVolume>(m, "GeomVolume");
	py::class_<gmds::cad::FACVolume, gmds::cad::GeomVolume>(m, "FACVolume")
	   .def("id", &gmds::cad::FACVolume::id)
	   .def("bbox", &gmds::cad::FACVolume::BBox)
	   .def("points", &gmds::cad::FACVolume::points, py::return_value_policy::reference)
	   .def("curves", &gmds::cad::FACVolume::curves, py::return_value_policy::reference)
	   .def("surfaces", &gmds::cad::FACVolume::surfaces, py::return_value_policy::reference);

}