/*----------------------------------------------------------------------------*/
#include <pybind11/pybind11.h>
#include <gmds/ig/Mesh.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
namespace py = pybind11;
//pygmds is the name of the module we export!!!
PYBIND11_MODULE(pygmds, m) {


    py::enum_<gmds::ECellType>(m, "CellType")
            .value("GMDS_NODE", gmds::ECellType::GMDS_NODE)
            .value("GMDS_EDGE", gmds::ECellType::GMDS_EDGE)
            .value("GMDS_FACE", gmds::ECellType::GMDS_FACE)
            .value("GMDS_REGION", gmds::ECellType::GMDS_REGION)
            .export_values();

    py::enum_<gmds::EMeshDefinition>(m, "ModelDefinition", py::arithmetic())
            .value("DIM3", gmds::EMeshDefinition::DIM3)
            .value("N", gmds::EMeshDefinition::N)
            .value("E", gmds::EMeshDefinition::E)
            .value("F", gmds::EMeshDefinition::F)
            .value("R", gmds::EMeshDefinition::R)
            .value("N2E", gmds::EMeshDefinition::N2E)
            .value("N2F", gmds::EMeshDefinition::N2F)
            .value("N2R", gmds::EMeshDefinition::N2R)
            .value("E2N", gmds::EMeshDefinition::E2N)
            .value("E2F", gmds::EMeshDefinition::E2F)
            .value("E2R", gmds::EMeshDefinition::E2R)
            .value("F2N", gmds::EMeshDefinition::F2N)
            .value("F2E", gmds::EMeshDefinition::F2E)
            .value("F2R", gmds::EMeshDefinition::F2R)
            .value("R2N", gmds::EMeshDefinition::R2N)
            .value("R2E", gmds::EMeshDefinition::R2E)
            .value("R2F", gmds::EMeshDefinition::R2F)
            .export_values();

    py::class_<gmds::MeshModel>(m, "MeshModel")
            .def(py::init<const int &>());

    py::class_<gmds::Variable<int>>(m, "VariableInt")
            .def(py::init<const std::string&>())
            .def("set", &gmds::Variable<int>::set)
            .def("value", static_cast<int& (gmds::Variable<int>::*)(const int&)>(&gmds::Variable<int>::value));



    py::class_<gmds::Mesh>(m, "Mesh")
            .def(py::init<const gmds::MeshModel&>())
            .def("getNbNodes", &gmds::Mesh::getNbNodes)
            .def("getNbEdges", &gmds::Mesh::getNbEdges)
            .def("getNbFaces", &gmds::Mesh::getNbFaces)
            .def("getNbRegions", &gmds::Mesh::getNbRegions)
            .def("getNode", &gmds::Mesh::get<gmds::Node>)
            .def("getEdge", &gmds::Mesh::get<gmds::Edge>)
            .def("getFace", &gmds::Mesh::get<gmds::Face>)
            .def("getRegion", &gmds::Mesh::get<gmds::Region>)
            .def("newIntVariableNode", [](gmds::Mesh& AThis, std::string &AName) {
                     return AThis.getOrCreateVariable<int, gmds::GMDS_NODE>("AName");
                 });

    py::class_<gmds::IMeshIOService>(m, "AbstractIOService");

    py::class_<gmds::IGMeshIOService, gmds::IMeshIOService>(m, "IOService")
            .def(py::init<gmds::Mesh*>());


    py::class_<gmds::VTKReader>(m, "VTKReader")
            .def(py::init<gmds::IMeshIOService*>())
            .def("read", &gmds::VTKReader::read)
            .def("setCellOptions",&gmds::VTKReader::setCellOptions)
            .def("setDataOptions",&gmds::VTKReader::setDataOptions);

    m.attr("__version__") = "dev";
}
