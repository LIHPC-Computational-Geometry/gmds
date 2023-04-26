/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CurvedBlocking.h"
#include <CGAL/Linear_cell_complex_base.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::blocking;
/*----------------------------------------------------------------------------*/
int CurvedBlocking::m_counter_nodes = 0;
int CurvedBlocking::m_counter_edges = 0;
int CurvedBlocking::m_counter_faces = 0;
int CurvedBlocking::m_counter_blocks = 0;
/*----------------------------------------------------------------------------*/
CurvedBlocking::CurvedBlocking(cad::GeomManager *AGeomModel, bool AInitAsBoundingBox):
m_geom_model(AGeomModel)
{
	if(AInitAsBoundingBox){
		TCoord min[3]={MAXFLOAT,MAXFLOAT,MAXFLOAT};
		TCoord max[3]={-MAXFLOAT,-MAXFLOAT,-MAXFLOAT};
		std::vector<cad::GeomVolume*> vols;
		m_geom_model->getVolumes(vols);
		for(auto v:vols){
			TCoord v_min[3], v_max[3];
			v->computeBoundingBox(v_min,v_max);
			for(auto i=0;i<2;i++)
				if (v_min[i]<min[i]) min[i]=v_min[i];
			for(auto i=0;i<2;i++)
				if (v_max[i]>max[i]) max[i]=v_max[i];
		}
		math::Point p1(min[0],min[1],min[2]);
		math::Point p2(min[0],max[1],min[2]);
		math::Point p3(max[0],max[1],min[2]);
		math::Point p4(max[0],min[1],min[2]);
		math::Point p5(min[0],min[1],max[2]);
		math::Point p6(min[0],max[1],max[2]);
		math::Point p7(max[0],max[1],max[2]);
		math::Point p8(max[0],min[1],max[2]);
		createBlock(p1,p2,p3,p4,p5,p6,p7,p8);
	}

}
/*----------------------------------------------------------------------------*/
CurvedBlocking::~CurvedBlocking() {}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Node
CurvedBlocking::createNode(const int AGeomDim, const int AGeomId, math::Point &APoint)
{
	return m_gmap.create_attribute<0>(NodeInfo(m_counter_nodes++, AGeomDim, AGeomId, APoint));
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Edge
CurvedBlocking::createEdge(const int AGeomDim, const int AGeomId)
{
	return m_gmap.create_attribute<1>(CellInfo(1, m_counter_edges++, AGeomDim, AGeomId));
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Face
CurvedBlocking::createFace(const int AGeomDim, const int AGeomId)
{
	return m_gmap.create_attribute<2>(CellInfo(2, m_counter_faces++, AGeomDim, AGeomId));
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Block
CurvedBlocking::createBlock(const int AGeomDim, const int AGeomId)
{
	return m_gmap.create_attribute<3>(CellInfo(3, m_counter_blocks++, AGeomDim, AGeomId));
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Face>
CurvedBlocking::get_faces_of_block(const CurvedBlocking::Block AB)
{
	Dart3 d1 = AB->dart();
	std::vector<CurvedBlocking::Face> faces;
	faces.resize(6);

	faces[0] = m_gmap.attribute<2>(d1);
	faces[1] = m_gmap.attribute<2>(m_gmap.alpha<2, 1, 0, 1, 2>(d1));
	faces[2] = m_gmap.attribute<2>(m_gmap.alpha<2>(d1));
	faces[3] = m_gmap.attribute<2>(m_gmap.alpha<1, 0, 1, 2>(d1));
	faces[4] = m_gmap.attribute<2>(m_gmap.alpha<1, 2>(d1));
	faces[5] = m_gmap.attribute<2>(m_gmap.alpha<0, 1, 2>(d1));
	return faces;
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_block(const CurvedBlocking::Block AB)
{
	Dart3 d1 = AB->dart();
	std::vector<CurvedBlocking::Node> nodes;
	nodes.resize(8);
	nodes[0] = m_gmap.attribute<0>(d1);
	nodes[1] = m_gmap.attribute<0>(m_gmap.alpha<0, 1>(d1));
	nodes[2] = m_gmap.attribute<0>(m_gmap.alpha<0, 1, 0, 1>(d1));
	nodes[3] = m_gmap.attribute<0>(m_gmap.alpha<1, 0>(d1));

	Dart3 d2 = m_gmap.alpha<2, 1, 0, 1, 2>(d1);
	nodes[4] = m_gmap.attribute<0>(d2);
	nodes[5] = m_gmap.attribute<0>(m_gmap.alpha<0, 1>(d2));
	nodes[6] = m_gmap.attribute<0>(m_gmap.alpha<0, 1, 0, 1>(d2));
	nodes[7] = m_gmap.attribute<0>(m_gmap.alpha<1, 0>(d2));
	return nodes;
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_face(const CurvedBlocking::Face AF)
{
	Dart3 d1 = AF->dart();
	std::vector<CurvedBlocking::Node> nodes;
	Dart3 d = d1;
	do {
		nodes.push_back(m_gmap.attribute<0>(d));
		d = m_gmap.alpha<0, 1>(d);
	} while (d != d1);
	return nodes;
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_edge(const CurvedBlocking::Edge AE)
{
	std::vector<CurvedBlocking::Node> nodes;
	nodes.reserve(2);
	Dart3 d= AE->dart();
	nodes[0] = m_gmap.attribute<0>(d);
	nodes[1] = m_gmap.attribute<0>(m_gmap.alpha<0>(d));
	return nodes;
}
/*----------------------------------------------------------------------------*/
math::Point
CurvedBlocking::get_center_of_face(const Face AF)
{
	std::vector<Node> nodes = get_nodes_of_face(AF);
	math::Point center(0, 0, 0);
	for (auto n : nodes) {
		center = center + n->info().point;
	}
	return (1.0 / nodes.size()) * center;
}
/*----------------------------------------------------------------------------*/
math::Point
CurvedBlocking::get_center_of_block(const Block AB)
{
	std::vector<Node> nodes = get_nodes_of_block(AB);
	math::Point center(0, 0, 0);
	for (auto n : nodes) {
		center = center + n->info().point;
	}
	return (1.0 / nodes.size()) * center;
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Block
CurvedBlocking::createBlock(
   math::Point &AP1, math::Point &AP2, math::Point &AP3, math::Point &AP4, math::Point &AP5, math::Point &AP6, math::Point &AP7, math::Point &AP8)
{
	// Initialize attribute for the created hexahedron
	Dart3 d1 = m_gmap.make_combinatorial_hexahedron();

	// 0D attribute
	m_gmap.set_attribute<0>(d1, createNode(1, 4, AP1));
	m_gmap.set_attribute<0>(m_gmap.alpha<0>(d1), createNode(1, 4, AP2));
	m_gmap.set_attribute<0>(m_gmap.alpha<0, 1, 0>(d1), createNode(1, 4, AP3));
	m_gmap.set_attribute<0>(m_gmap.alpha<1, 0>(d1), createNode(1, 4, AP4));
	m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2>(d1), createNode(1, 4, AP5));
	m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2, 0, 1>(d1), createNode(1, 4, AP6));
	m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2, 0, 1, 0, 1>(d1), createNode(1, 4, AP7));
	m_gmap.set_attribute<0>(m_gmap.alpha<2, 1, 0, 1, 2, 1, 0>(d1), createNode(1, 4, AP8));

	// go through all the edges
	for (auto it = m_gmap.one_dart_per_incident_cell<1, 3>(d1).begin(), itend = m_gmap.one_dart_per_incident_cell<1, 3>(d1).end(); it != itend; ++it) {
		m_gmap.set_attribute<1>(it, createEdge(4, NullID));
	}

	// go through all the faces
	for (auto it = m_gmap.one_dart_per_incident_cell<2, 3>(d1).begin(), itend = m_gmap.one_dart_per_incident_cell<2, 3>(d1).end(); it != itend; ++it) {
		m_gmap.set_attribute<2>(it, createFace(4, NullID));
	}
	Block b = createBlock(4, NullID);
	m_gmap.set_attribute<3>(d1, b);
	return b;
}
/*----------------------------------------------------------------------------*/
bool
CurvedBlocking::isValidTopology() const
{
	return m_gmap.is_valid();
}
/*----------------------------------------------------------------------------*/
std::string
CurvedBlocking::info() const
{
	std::ostringstream mess;
	mess << "Blocking Info: " << std::endl;
	m_gmap.display_characteristics(mess);
	mess << ", validity=" << (isValidTopology() ? "true" : "false");
	mess << "\nOd-attributes: ";
	for (auto it = m_gmap.attributes<0>().begin(), itend = m_gmap.attributes<0>().end(); it != itend; ++it) {
		auto at = m_gmap.info_of_attribute<0>(it);
		mess << "[" << at.point << ",(dim:" << at.geom_dim << ",id:" << at.geom_id << ")]; ";
	}
	mess << "\n1d-attributes: ";
	for (auto it = m_gmap.attributes<1>().begin(), itend = m_gmap.attributes<1>().end(); it != itend; ++it) {
		auto at = m_gmap.info_of_attribute<1>(it);
		mess << "(dim:" << at.geom_dim << ",id:" << at.geom_id << "); ";
	}
	mess << "\n2d-attributes: ";
	for (auto it = m_gmap.attributes<2>().begin(), itend = m_gmap.attributes<2>().end(); it != itend; ++it) {
		auto at = m_gmap.info_of_attribute<2>(it);
		mess << "(dim:" << at.geom_dim << ",id:" << at.geom_id << "); ";
	}
	mess << "\n3d-attributes: ";
	for (auto it = m_gmap.attributes<3>().begin(), itend = m_gmap.attributes<3>().end(); it != itend; ++it) {
		auto at = m_gmap.info_of_attribute<3>(it);
		mess << "(dim:" << at.geom_dim << ",id:" << at.geom_id << "); ";
	}
	return mess.str();
}
/*----------------------------------------------------------------------------*/
void CurvedBlocking::convert_to_mesh(Mesh &AMesh)
{
	MeshModel model = AMesh.getModel();
	if(!model.has(N) || !model.has(E) ||!model.has(F) ||!model.has(R) ||
	    !model.has(E2N) ||!model.has(F2N) ||!model.has(R2N))
		throw GMDSException("Wrong mesh model for block->mesh conversion");

	AMesh.clear();

	Variable<int>* var_node_topo_id  = AMesh.getOrCreateVariable<int, GMDS_NODE>("blocking_topo_id");
	Variable<int>* var_node_topo_dim = AMesh.getOrCreateVariable<int, GMDS_NODE>("blocking_topo_dim");
	Variable<int>* var_node_geom_id  = AMesh.getOrCreateVariable<int, GMDS_NODE>("blocking_geom_id");
	Variable<int>* var_node_geom_dim = AMesh.getOrCreateVariable<int, GMDS_NODE>("blocking_geom_dim");

	Variable<int>* var_edge_topo_id  = AMesh.getOrCreateVariable<int, GMDS_EDGE>("blocking_topo_id");
	Variable<int>* var_edge_topo_dim = AMesh.getOrCreateVariable<int, GMDS_EDGE>("blocking_topo_dim");
	Variable<int>* var_edge_geom_id  = AMesh.getOrCreateVariable<int, GMDS_EDGE>("blocking_geom_id");
	Variable<int>* var_edge_geom_dim = AMesh.getOrCreateVariable<int, GMDS_EDGE>("blocking_geom_dim");

	Variable<int>* var_face_topo_id  = AMesh.getOrCreateVariable<int, GMDS_FACE>("blocking_topo_id");
	Variable<int>* var_face_topo_dim = AMesh.getOrCreateVariable<int, GMDS_FACE>("blocking_topo_dim");
	Variable<int>* var_face_geom_id  = AMesh.getOrCreateVariable<int, GMDS_FACE>("blocking_geom_id");
	Variable<int>* var_face_geom_dim = AMesh.getOrCreateVariable<int, GMDS_FACE>("blocking_geom_dim");

	Variable<int>* var_region_topo_id  = AMesh.getOrCreateVariable<int, GMDS_REGION>("blocking_topo_id");
	Variable<int>* var_region_topo_dim = AMesh.getOrCreateVariable<int, GMDS_REGION>("blocking_topo_dim");
	Variable<int>* var_region_geom_id  = AMesh.getOrCreateVariable<int, GMDS_REGION>("blocking_geom_id");
	Variable<int>* var_region_geom_dim = AMesh.getOrCreateVariable<int, GMDS_REGION>("blocking_geom_dim");

	//mapping from blocking node ids to mesh node ids
	std::map<int,TCellID> n2n;

	//nodes
	for (auto it = m_gmap.attributes<0>().begin(),
	          itend = m_gmap.attributes<0>().end(); it != itend; ++it)
	{
		auto att = m_gmap.info_of_attribute<0>(it);
		gmds::Node n = AMesh.newNode(att.point);
		var_node_topo_id->set(n.id(),att.topo_id);
		var_node_topo_dim->set(n.id(),att.topo_dim);
		var_node_geom_id->set(n.id(),att.geom_id);
		var_node_geom_dim->set(n.id(),att.geom_dim);
		n2n[att.topo_id]=n.id();
	}
	//edges
	for (auto it = m_gmap.attributes<1>().begin(),
	          itend = m_gmap.attributes<1>().end(); it != itend; ++it)
	{
		auto att = m_gmap.info_of_attribute<1>(it);
		std::vector<Node> cell_nodes = get_nodes_of_edge(it);
		gmds::Edge e = AMesh.newEdge(n2n[cell_nodes[0]->info().topo_id],
		                             n2n[cell_nodes[1]->info().topo_id]);
		var_edge_topo_id->set(e.id(),att.topo_id);
		var_edge_topo_dim->set(e.id(),att.topo_dim);
		var_edge_geom_id->set(e.id(),att.geom_id);
		var_edge_geom_dim->set(e.id(),att.geom_dim);
	}
	//faces
	for (auto it = m_gmap.attributes<2>().begin(),
	          itend = m_gmap.attributes<2>().end(); it != itend; ++it)
	{
		auto att = m_gmap.info_of_attribute<2>(it);
		std::vector<Node> cell_nodes = get_nodes_of_face(it);
		if(cell_nodes.size()!=4)
			throw GMDSException("Only quad blocking faces can be converted into mesh");

		gmds::Face f = AMesh.newQuad(n2n[cell_nodes[0]->info().topo_id],
		                             n2n[cell_nodes[1]->info().topo_id],
		                             n2n[cell_nodes[2]->info().topo_id],
		                             n2n[cell_nodes[3]->info().topo_id]);

		var_face_topo_id->set(f.id(),att.topo_id);
		var_face_topo_dim->set(f.id(),att.topo_dim);
		var_face_geom_id->set(f.id(),att.geom_id);
		var_face_geom_dim->set(f.id(),att.geom_dim);
	}
	//blocks
	for (auto it = m_gmap.attributes<3>().begin(),
	          itend = m_gmap.attributes<3>().end(); it != itend; ++it)
	{
		auto att = m_gmap.info_of_attribute<3>(it);
		std::vector<Node> cell_nodes = get_nodes_of_block(it);
		if(cell_nodes.size()!=8)
			throw GMDSException("Only hex blocks can be converted into mesh");

		gmds::Region r = AMesh.newHex(n2n[cell_nodes[0]->info().topo_id],
		                             n2n[cell_nodes[1]->info().topo_id],
		                             n2n[cell_nodes[2]->info().topo_id],
		                             n2n[cell_nodes[3]->info().topo_id],
		                               n2n[cell_nodes[4]->info().topo_id],
		                               n2n[cell_nodes[5]->info().topo_id],
		                               n2n[cell_nodes[6]->info().topo_id],
		                               n2n[cell_nodes[7]->info().topo_id]);

		var_region_topo_id->set(r.id(),att.topo_id);
		var_region_topo_dim->set(r.id(),att.topo_dim);
		var_region_geom_id->set(r.id(),att.geom_id);
		var_region_geom_dim->set(r.id(),att.geom_dim);
	}
}