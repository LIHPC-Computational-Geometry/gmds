/*----------------------------------------------------------------------------*/
/*
 * Node.cpp
 *
 *  Created on: 5 f�vr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Node.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/EdgeContainer.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
Node::
Node()
: Cell(nullptr,GMDS_NODE,NullID)
{
	m_nodes_container = nullptr;
}
/*----------------------------------------------------------------------------*/
Node::
Node(Mesh* AMesh,  const TCellID& AID,
		const TCoord& AX, const TCoord& AY, const TCoord& AZ)
: Cell(AMesh,GMDS_NODE,AID)
{
	if(AMesh!=nullptr){
		m_nodes_container  = AMesh->m_nodes_container;
		m_nodes_container->m_node_coords[m_id].setXYZ(AX,AY,AZ);
	}
	else
		m_nodes_container=nullptr;
}
/*----------------------------------------------------------------------------*/
Node::
Node(Mesh* AMesh, const TCellID& AID,
		const math::Point& APt)
: Cell(AMesh,GMDS_NODE,AID)
{
	if(AMesh!=nullptr){
		m_nodes_container  = AMesh->m_nodes_container;
		m_nodes_container->m_node_coords[m_id]=APt;
	}
	else
		m_nodes_container=nullptr;
}
/*----------------------------------------------------------------------------*/
Node::
Node(const Node& ANode)
: Cell(ANode.m_owner,GMDS_NODE,ANode.m_id)
{
	if(m_owner!=nullptr)
		m_nodes_container  = m_owner->m_nodes_container;
	else
		m_nodes_container = nullptr;
}
/*----------------------------------------------------------------------------*/
void Node::operator=(const Node& ANode)
{
	m_owner = ANode.m_owner;
	m_id = ANode.m_id;
	if(m_owner!=nullptr)
		m_nodes_container  = m_owner->m_nodes_container;
	else
		m_nodes_container = nullptr;
}
/*----------------------------------------------------------------------------*/
bool Node::operator==(const Node& ANode) const
{
	return (m_owner == ANode.m_owner && m_id==ANode.m_id);
}
/*----------------------------------------------------------------------------*/
bool Node::operator!=(const Node& ANode) const
{
        return (!(*this == ANode));
}
/*----------------------------------------------------------------------------*/
Node::~Node()
= default;
/*----------------------------------------------------------------------------*/
TInt Node::nbNodes() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2N)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Node::nbEdges() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2E)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Node::nbFaces() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2F)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Node::nbRegions() const
{
	TabCellID<size_undef> cells = (*m_nodes_container->m_N2R)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
gmds::math::Point Node::center() const
{
        return m_nodes_container->m_node_coords[m_id];
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2N)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2E)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2F)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGet(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2R)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetNodeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");

	(*m_nodes_container->m_N2N)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetEdgeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetFaceIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetRegionIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2N)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2E)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2F)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAll(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_nodes_container->m_N2R)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2N))
                throw GMDSException("N2N adjacency is not supported by the mesh model");

        (*m_nodes_container->m_N2N)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2E))
                throw GMDSException("N2E adjacency is not supported by the mesh model");
        (*m_nodes_container->m_N2E)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2F))
                throw GMDSException("N2F adjacency is not supported by the mesh model");
        (*m_nodes_container->m_N2F)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(N2R))
                throw GMDSException("N2R adjacency is not supported by the mesh model");
        (*m_nodes_container->m_N2R)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetNodeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetEdgeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetFaceIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(N2F))
			throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateSetRegionIDs(const std::vector<TCellID>& ACells)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Node::delegateNodeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateEdgeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateFaceAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateRegionAdd(TCellID AID)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateNodeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateEdgeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateFaceRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateRegionRemove(TCellID AID)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Node::delegateNodeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2N))
		throw GMDSException("N2N adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2N)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Node::delegateEdgeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2E))
		throw GMDSException("N2E adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2E)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Node::delegateFaceReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2F))
		throw GMDSException("N2F adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2F)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Node::delegateRegionReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(N2R))
		throw GMDSException("N2R adjacency is not supported by the mesh model");
	(*m_nodes_container->m_N2R)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
math::Point Node::point() const
{
	return  m_nodes_container->m_node_coords[m_id];
}
/*----------------------------------------------------------------------------*/
void Node::setPoint(const math::Point& APnt)
{
	m_nodes_container->m_node_coords[m_id] = APnt;
}
/*----------------------------------------------------------------------------*/
TCoord Node::X() const
{
	return m_nodes_container->m_node_coords[m_id].X();
}
/*----------------------------------------------------------------------------*/
TCoord Node::Y() const
{
	return m_nodes_container->m_node_coords[m_id].Y();
}
/*----------------------------------------------------------------------------*/
TCoord Node::Z() const
{
	return m_nodes_container->m_node_coords[m_id].Z();
}
/*----------------------------------------------------------------------------*/
TCoord& Node::X()
{
	return m_nodes_container->m_node_coords[m_id].X();
}
/*----------------------------------------------------------------------------*/
TCoord& Node::Y()
{
	return m_nodes_container->m_node_coords[m_id].Y();
}
/*----------------------------------------------------------------------------*/
TCoord& Node::Z()
{
	return m_nodes_container->m_node_coords[m_id].Z();
}
/*----------------------------------------------------------------------------*/
void Node::setX(const TCoord AVal)
{
	m_nodes_container->m_node_coords[m_id].setX(AVal);
}
/*----------------------------------------------------------------------------*/
void Node::setY(const TCoord AVal)
{
	m_nodes_container->m_node_coords[m_id].setY(AVal);
}
/*----------------------------------------------------------------------------*/
void Node::setZ(const TCoord AVal)
{
	m_nodes_container->m_node_coords[m_id].setZ(AVal);
}
/*----------------------------------------------------------------------------*/
void Node::setXYZ(const TCoord AX, const TCoord AY,const TCoord AZ)
{
	m_nodes_container->m_node_coords[m_id].setXYZ(AX,AY,AZ);
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStream, const Node& AN)
{
        AStream << "Node " << AN.id() << " ("
                << AN.point().X() << ", "
                << AN.point().Y() << ", "
                << AN.point().Z() << ")";
        return AStream;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
