
/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
int GeomMeshLinker::m_global_link_id = 0;
/*----------------------------------------------------------------------------*/
GeomMeshLinker::GeomMeshLinker() : m_link_id(m_global_link_id++), m_mesh(nullptr), m_geometry(nullptr) {}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::GeomMeshLinker(Mesh *AMesh, GeomManager *AGeometry) : m_link_id(m_global_link_id++), m_mesh(nullptr), m_geometry(nullptr)
{
	setMesh(AMesh);
	setGeometry(AGeometry);
}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::~GeomMeshLinker()
{

	clear();
}

/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::clear()
{
	if (m_mesh != nullptr) {
		m_mesh->deleteVariable(GMDS_NODE, m_node_classification_dim);
		m_mesh->deleteVariable(GMDS_NODE, m_node_classification_id);
		m_mesh->deleteVariable(GMDS_EDGE, m_node_classification_dim);
		m_mesh->deleteVariable(GMDS_EDGE, m_edge_classification_id);
		m_mesh->deleteVariable(GMDS_FACE, m_face_classification_dim);
		m_mesh->deleteVariable(GMDS_FACE, m_face_classification_id);
	}
	m_mesh = nullptr;
	m_geometry = nullptr;
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::setMesh(Mesh *AMesh)
{
	if (m_mesh != nullptr) {
		clear();
	}

	m_mesh = AMesh;

	m_node_classification_dim = m_mesh->newVariable<eLink, GMDS_NODE>("geom_link_dim_" + std::to_string(m_link_id));
	m_node_classification_id = m_mesh->newVariable<int, GMDS_NODE>("geom_link_id_" + std::to_string(m_link_id));

	m_edge_classification_dim = m_mesh->newVariable<eLink, GMDS_EDGE>("geom_link_dim_" + std::to_string(m_link_id));
	m_edge_classification_id = m_mesh->newVariable<int, GMDS_EDGE>("geom_link_id_" + std::to_string(m_link_id));

	m_face_classification_dim = m_mesh->newVariable<eLink, GMDS_FACE>("geom_link_dim_" + std::to_string(m_link_id));
	m_face_classification_id = m_mesh->newVariable<int, GMDS_FACE>("geom_link_id_" + std::to_string(m_link_id));
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::setGeometry(GeomManager *AGeometry)
{
	m_geometry = AGeometry;
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkToPoint(const Node &AN, const int AGeomId)
{
	linkNodeToPoint(AN.id(), AGeomId);
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkNodeToPoint(const TCellID &AN, const int AGeomId)
{
	(*m_node_classification_dim)[AN] = LinkPoint;
	(*m_node_classification_id)[AN] = AGeomId;
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkToCurve(const Edge &AE, const int AGeomId)
{
	linkEdgeToCurve(AE.id(), AGeomId);
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkEdgeToCurve(const TCellID &AE, const int AGeomId)
{
	(*m_edge_classification_dim)[AE] = LinkCurve;
	(*m_edge_classification_id)[AE] = AGeomId;
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkToCurve(const Node &AN, const int AGeomId)
{
	linkNodeToCurve(AN.id(), AGeomId);
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkNodeToCurve(const TCellID &AN, const int AGeomId)
{
	(*m_node_classification_dim)[AN] = LinkCurve;
	(*m_node_classification_id)[AN] = AGeomId;
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkToSurface(const Node &AN, const int AGeomId)
{
	linkNodeToSurface(AN.id(), AGeomId);
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkNodeToSurface(const TCellID &AN, const int AGeomId)
{
	(*m_node_classification_dim)[AN] = LinkSurface;
	(*m_node_classification_id)[AN] = AGeomId;
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkToSurface(const Edge &AE, const int AGeomId)
{
	linkEdgeToSurface(AE.id(), AGeomId);
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkEdgeToSurface(const TCellID &AE, const int AGeomId)
{
	(*m_edge_classification_dim)[AE] = LinkSurface;
	(*m_edge_classification_id)[AE] = AGeomId;
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkToSurface(const Face &AF, const int AGeomId)
{
	linkFaceToSurface(AF.id(), AGeomId);
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::linkFaceToSurface(const TCellID &AF, const int AGeomId)
{
	(*m_face_classification_dim)[AF] = LinkSurface;
	(*m_face_classification_id)[AF] = AGeomId;
}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::eLink
GeomMeshLinker::getGeomDim(const Node &AN)
{
	return getGeomDim<Node>(AN.id());
}
/*----------------------------------------------------------------------------*/
int
GeomMeshLinker::getGeomId(const Node &AN)
{
	return getGeomId<Node>(AN.id());
}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::eLink
GeomMeshLinker::getGeomDim(const Edge &AE)
{
	return getGeomDim<Edge>(AE.id());
}
/*----------------------------------------------------------------------------*/
int
GeomMeshLinker::getGeomId(const Edge &AE)
{
	return getGeomId<Edge>(AE.id());
}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::eLink
GeomMeshLinker::getGeomDim(const Face &AF)
{
	return getGeomDim<Face>(AF.id());
}
/*----------------------------------------------------------------------------*/
int
GeomMeshLinker::getGeomId(const Face &AF)
{
	return getGeomId<Face>(AF.id());
}
/*----------------------------------------------------------------------------*/
std::pair<GeomMeshLinker::eLink, int>
GeomMeshLinker::getGeomInfo(const Node &AN)
{
	return getGeomInfo<Node>(AN.id());
}
/*----------------------------------------------------------------------------*/
std::pair<GeomMeshLinker::eLink, int>
GeomMeshLinker::getGeomInfo(const Edge &AN)
{
	return getGeomInfo<Edge>(AN.id());
}
/*----------------------------------------------------------------------------*/
std::pair<GeomMeshLinker::eLink, int>
GeomMeshLinker::getGeomInfo(const Face &AN)
{
	return getGeomInfo<Face>(AN.id());
}
/*----------------------------------------------------------------------------*/
void
GeomMeshLinker::writeVTKDebugMesh(const std::string &AFileName)
{
	Mesh m(MeshModel(DIM3 | F | E | N | F2N | E2N));
	std::map<TCellID, TCellID> n2n;

	Variable<int> *var_node_id = m.newVariable<int, GMDS_NODE>("REF_ID");
	Variable<int> *var_edge_id = m.newVariable<int, GMDS_EDGE>("REF_ID");
	Variable<int> *var_face_id = m.newVariable<int, GMDS_FACE>("REF_ID");

	Variable<int> *var_node_geom_dim = m.newVariable<int, GMDS_NODE>("classif_dim");
	Variable<int> *var_node_geom_id = m.newVariable<int, GMDS_NODE>("classif_id");

	Variable<int> *var_edge_geom_dim = m.newVariable<int, GMDS_EDGE>("classif_dim");
	Variable<int> *var_edge_geom_id = m.newVariable<int, GMDS_EDGE>("classif_id");

	Variable<int> *var_face_geom_dim = m.newVariable<int, GMDS_FACE>("classif_dim");
	Variable<int> *var_face_geom_id = m.newVariable<int, GMDS_FACE>("classif_id");

	for (auto n_id : m_mesh->nodes()) {
		if (m_node_classification_dim->value(n_id) == GeomMeshLinker::LinkSurface || m_node_classification_dim->value(n_id) == GeomMeshLinker::LinkCurve
		    || m_node_classification_dim->value(n_id) == GeomMeshLinker::LinkPoint) {
			// so a node classified on the boundary
			Node from_node = m_mesh->get<Node>(n_id);
			Node to_node = m.newNode(from_node.point());
			n2n[n_id] = to_node.id();
			var_node_id->value(to_node.id()) = n_id;

			var_node_geom_dim->value(to_node.id()) = m_node_classification_dim->value(n_id) - 1;
			var_node_geom_id->value(to_node.id()) = m_node_classification_id->value(n_id);
		}
	}
	for (auto e_id : m_mesh->edges()) {
		if (m_edge_classification_dim->value(e_id) < 4 && m_edge_classification_dim->value(e_id) > 1) {
			// so an edge classified on the boundary
			Edge from_edge = m_mesh->get<Edge>(e_id);
			std::vector<TCellID> end_points = from_edge.getIDs<Node>();
			Edge to_edge = m.newEdge(n2n[end_points[0]], n2n[end_points[1]]);
			var_edge_id->value(to_edge.id()) = e_id;
			var_edge_geom_dim->value(to_edge.id()) = m_edge_classification_dim->value(e_id) - 1;
			var_edge_geom_id->value(to_edge.id()) = m_edge_classification_id->value(e_id);
		}
	}
	for (auto f_id : m_mesh->faces()) {
		if (m_face_classification_dim->value(f_id) == GeomMeshLinker::LinkSurface) {
			// so an edge classified on the boundary
			Face from_face = m_mesh->get<Face>(f_id);
			std::vector<TCellID> from_node_ids = from_face.getIDs<Node>();
			std::vector<TCellID> to_node_ids;
			to_node_ids.resize(from_node_ids.size());
			for (auto i = 0; i < from_node_ids.size(); i++)
				to_node_ids[i] = n2n[from_node_ids[i]];
			Face to_face = m.newFace(to_node_ids);
			var_face_id->value(to_face.id()) = f_id;
			var_face_geom_dim->value(to_face.id()) = m_face_classification_dim->value(f_id) - 1;
			var_face_geom_id->value(to_face.id()) = m_face_classification_id->value(f_id);
		}
	}

	IGMeshIOService ios(&m);

	VTKWriter vtkWriter(&ios);
	vtkWriter.setCellOptions(N | F | E);
	vtkWriter.setDataOptions(N | E | F);
	vtkWriter.write(AFileName);
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API GeomMeshLinker::eLink
GeomMeshLinker::getGeomDim<Node>(const TCellID &AN)
{
	return (*m_node_classification_dim)[AN];
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API int
GeomMeshLinker::getGeomId<Node>(const TCellID &AN)
{
	return (*m_node_classification_id)[AN];
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API GeomMeshLinker::eLink
GeomMeshLinker::getGeomDim<Edge>(const TCellID &AN)
{
	return (*m_edge_classification_dim)[AN];
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API int
GeomMeshLinker::getGeomId<Edge>(const TCellID &AN)
{
	return (*m_edge_classification_id)[AN];
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API GeomMeshLinker::eLink
GeomMeshLinker::getGeomDim<Face>(const TCellID &AN)
{
	return (*m_face_classification_dim)[AN];
}
/*----------------------------------------------------------------------------*/
template<> int
GeomMeshLinker::getGeomId<Face>(const TCellID &AN)
{
	return (*m_face_classification_id)[AN];
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API std::pair<GeomMeshLinker::eLink, int>
GeomMeshLinker::getGeomInfo<Node>(const TCellID &AN)
{
	return std::make_pair((*m_node_classification_dim)[AN], (*m_node_classification_id)[AN]);
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API std::pair<GeomMeshLinker::eLink, int>
GeomMeshLinker::getGeomInfo<Edge>(const TCellID &AE)
{
	return std::make_pair((*m_edge_classification_dim)[AE], (*m_edge_classification_id)[AE]);
}
/*----------------------------------------------------------------------------*/
template<> GMDSCad_API std::pair<GeomMeshLinker::eLink, int>
GeomMeshLinker::getGeomInfo<Face>(const TCellID &AF)
{
	return std::make_pair((*m_edge_classification_dim)[AF], (*m_edge_classification_id)[AF]);
}
/*----------------------------------------------------------------------------*/
