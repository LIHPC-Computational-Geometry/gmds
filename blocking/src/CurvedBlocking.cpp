/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CurvedBlocking.h"
#include <CGAL/Linear_cell_complex_base.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::blocking;
/*----------------------------------------------------------------------------*/
int CurvedBlocking::m_counter_nodes=0;
int CurvedBlocking::m_counter_edges=0;
int CurvedBlocking::m_counter_faces=0;
int CurvedBlocking::m_counter_blocks=0;
/*----------------------------------------------------------------------------*/

CurvedBlocking::CurvedBlocking()
{}
/*----------------------------------------------------------------------------*/
CurvedBlocking::~CurvedBlocking()
{}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Node
CurvedBlocking::createNode(const int AGeomDim, const int AGeomId,
                           math::Point& APoint){
	return m_gmap.create_attribute<0>(NodeInfo(m_counter_nodes++,AGeomDim,AGeomId,APoint));
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Edge
CurvedBlocking::createEdge(const int AGeomDim, const int AGeomId) {
	return m_gmap.create_attribute<1>(CellInfo(1,m_counter_edges++,AGeomDim,AGeomId));
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Face
CurvedBlocking::createFace(const int AGeomDim, const int AGeomId) {
	return m_gmap.create_attribute<2>(CellInfo(2,m_counter_faces++,AGeomDim,AGeomId));
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Block
CurvedBlocking::createBlock(const int AGeomDim, const int AGeomId) {
	return m_gmap.create_attribute<3>(CellInfo(3,m_counter_blocks++,AGeomDim,AGeomId));
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Face>
CurvedBlocking::get_faces_of_block(const CurvedBlocking::Block AB){
	DartDescriptor d1 = AB->dart();
	std::vector<CurvedBlocking::Face> faces;
	faces.resize(6);

	faces[0]= m_gmap.attribute<2>(d1);
	faces[1]= m_gmap.attribute<2>(m_gmap.alpha<2,1,0,1,2>(d1));
	faces[2]= m_gmap.attribute<2>(m_gmap.alpha<2>(d1));
	faces[3]= m_gmap.attribute<2>(m_gmap.alpha<1,0,1,2>(d1));
	faces[4]= m_gmap.attribute<2>(m_gmap.alpha<1,2>(d1));
	faces[5]= m_gmap.attribute<2>(m_gmap.alpha<0,1,2>(d1));
	return faces;
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_block(const CurvedBlocking::Block AB){
	DartDescriptor d1 = AB->dart();
	std::vector<CurvedBlocking::Node> nodes;
	nodes.resize(8);
	nodes[0]= m_gmap.attribute<0>(d1);
	nodes[1]= m_gmap.attribute<0>(m_gmap.alpha<0,1>(d1));
	nodes[2]= m_gmap.attribute<0>(m_gmap.alpha<0,1,0,1>(d1));
	nodes[3]= m_gmap.attribute<0>(m_gmap.alpha<1,0>(d1));

	DartDescriptor d2 = m_gmap.alpha<2,1,0,1,2>(d1);
	nodes[4]= m_gmap.attribute<0>(d2);
	nodes[5]= m_gmap.attribute<0>(m_gmap.alpha<0,1>(d2));
	nodes[6]= m_gmap.attribute<0>(m_gmap.alpha<0,1,0,1>(d2));
	nodes[7]= m_gmap.attribute<0>(m_gmap.alpha<1,0>(d2));
	return nodes;
}
/*----------------------------------------------------------------------------*/
std::vector<CurvedBlocking::Node>
CurvedBlocking::get_nodes_of_face(const CurvedBlocking::Face AF){
	DartDescriptor d1 = AF->dart();
	std::vector<CurvedBlocking::Node> nodes;
	DartDescriptor d=d1;
	do{
		nodes.push_back(m_gmap.attribute<0>(d));
		d= m_gmap.alpha<0,1>(d);
	} while (d!=d1);
	return nodes;
}
/*----------------------------------------------------------------------------*/
math::Point CurvedBlocking::get_center_of_face(const Face AF){
	std::vector<Node> nodes = get_nodes_of_face(AF);
	math::Point center(0,0,0);
	for(auto n:nodes){
		center = center+n->info().point;
	}
	return (1.0/nodes.size())*center;
}
/*----------------------------------------------------------------------------*/
math::Point CurvedBlocking::get_center_of_block(const Block AB){
	std::vector<Node> nodes = get_nodes_of_block(AB);
	math::Point center(0,0,0);
	for(auto n:nodes){
		center = center+n->info().point;
	}
	return (1.0/nodes.size())*center;
}
/*----------------------------------------------------------------------------*/
CurvedBlocking::Block
CurvedBlocking::createBlock(math::Point& AP1, math::Point& AP2,
                            math::Point& AP3, math::Point& AP4,
                            math::Point& AP5, math::Point& AP6,
                            math::Point& AP7, math::Point& AP8)
{
	// Initialize attribute for the created hexahedron
	DartDescriptor d1 = m_gmap.make_combinatorial_hexahedron();

	//0D attribute
	m_gmap.set_attribute<0>(d1, createNode(1,4,AP1));
	m_gmap.set_attribute<0>(m_gmap.alpha<0>(d1),createNode(1,4,AP2));
	m_gmap.set_attribute<0>(m_gmap.alpha<0,1,0>(d1),
	                             createNode(1,4,AP3));
	m_gmap.set_attribute<0>(m_gmap.alpha<1,0>(d1),createNode(1,4,AP4));
	m_gmap.set_attribute<0>(m_gmap.alpha<2,1,0,1,2>(d1),createNode(1,4,AP5));
	m_gmap.set_attribute<0>(m_gmap.alpha<2,1,0,1,2, 0,1>(d1),createNode(1,4,AP6));
	m_gmap.set_attribute<0>(m_gmap.alpha<2,1,0,1,2, 0,1,0,1>(d1),createNode(1,4,AP7));
	m_gmap.set_attribute<0>(m_gmap.alpha<2,1,0,1,2, 1,0>(d1),createNode(1,4,AP8));

	//go through all the edges
	for (auto it=m_gmap.one_dart_per_incident_cell<1,3>(d1).begin(),
	          itend=m_gmap.one_dart_per_incident_cell<1,3>(d1).end();
	     it!=itend; ++it){
		m_gmap.set_attribute<1>(it, createEdge(4,NullID));
	}

	//go through all the faces
	for (auto it=m_gmap.one_dart_per_incident_cell<2,3>(d1).begin(),
	          itend=m_gmap.one_dart_per_incident_cell<2,3>(d1).end();
	     it!=itend; ++it){
		m_gmap.set_attribute<2>(it, createFace(4,NullID));
	}
	Block b = createBlock(4,NullID);
	m_gmap.set_attribute<3>(d1, b);
	return b;
}
/*----------------------------------------------------------------------------*/
std::string CurvedBlocking::info() const{
	std::ostringstream mess;
	mess<<"Blocking Info: "<<std::endl;
	m_gmap.display_characteristics(mess);
	mess<<", validity="<<(isValidTopology()?"true":"false");
	mess<<"\nOd-attributes: ";
	for (auto it=m_gmap.attributes<0>().begin(),
	          itend=m_gmap.attributes<0>().end(); it!=itend; ++it)
	{
		auto at = m_gmap.info_of_attribute<0>(it);
		mess<<"["<<at.point<<",(dim:"<<at.geom_dim<<",id:"<<at.geom_id<<")]; ";
	}
	mess<<"\n1d-attributes: ";
	for (auto it=m_gmap.attributes<1>().begin(),
	          itend=m_gmap.attributes<1>().end(); it!=itend; ++it)
	{
		auto at = m_gmap.info_of_attribute<1>(it);
		mess<<"(dim:"<<at.geom_dim<<",id:"<<at.geom_id<<"); ";
	}
	mess<<"\n2d-attributes: ";
	for (auto it=m_gmap.attributes<2>().begin(),
	          itend=m_gmap.attributes<2>().end(); it!=itend; ++it)
	{ 	auto at = m_gmap.info_of_attribute<2>(it);
		mess<<"(dim:"<<at.geom_dim<<",id:"<<at.geom_id<<"); "; }
	mess<<"\n3d-attributes: ";
	for (auto it=m_gmap.attributes<3>().begin(),
	          itend=m_gmap.attributes<3>().end(); it!=itend; ++it)
	{ 	auto at = m_gmap.info_of_attribute<3>(it);
		mess<<"(dim:"<<at.geom_dim<<",id:"<<at.geom_id<<"); ";}
	return mess.str();
}
