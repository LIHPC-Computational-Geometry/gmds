/*----------------------------------------------------------------------------*/
/*
 * Edge.cpp
 *
 *  Created on: 19 may 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Edge.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/EdgeContainer.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Edge::
Edge()
: Cell(nullptr,GMDS_EDGE,NullID)
{
	m_edges_container = nullptr;
}
/*----------------------------------------------------------------------------*/
Edge::
Edge(Mesh* AMesh, const TCellID& AID)
: Cell(AMesh,GMDS_EDGE,AID)
{
	if(AMesh!=nullptr)
		m_edges_container  = AMesh->m_edges_container;
	else
		m_edges_container  = nullptr;
}
/*----------------------------------------------------------------------------*/
Edge::
Edge(const Edge& AEdge)
: Cell(AEdge.m_owner,GMDS_EDGE,AEdge.m_id)
{
	if(m_owner!=nullptr)
		m_edges_container = m_owner->m_edges_container;
	else
		m_edges_container = nullptr;
}
/*----------------------------------------------------------------------------*/
Edge::~Edge()
{}
/*----------------------------------------------------------------------------*/
bool Edge::operator==(const Edge& AEdge) const
{
        return (m_owner == AEdge.m_owner && m_id==AEdge.m_id);
}
/*----------------------------------------------------------------------------*/
bool Edge::operator!=(const Edge& AEdge) const
{
        return (!(*this == AEdge));
}
/*----------------------------------------------------------------------------*/
TInt Edge::nbNodes() const
{
	return 2;
}
/*----------------------------------------------------------------------------*/
TInt Edge::nbEdges() const
{
	TabCellID<size_undef> cells = (*m_edges_container->m_E2E)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Edge::nbFaces() const
{
	TabCellID<size_undef> cells = (*m_edges_container->m_E2F)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TInt Edge::nbRegions() const
{
	TabCellID<size_undef> cells = (*m_edges_container->m_E2R)[m_id];
	return cells.size();
}
/*----------------------------------------------------------------------------*/
TCoord Edge::length() const
{
	std::vector<Node> nodes;
	get<Node>(nodes);
	return math::Vector3d(nodes[0].point(), nodes[1].point()).norm();
}
/*----------------------------------------------------------------------------*/
math::Point Edge::center() const
{
	std::vector<Node> nodes;
	get<Node>(nodes);
	return 0.5* nodes[0].point() + 0.5 * nodes[1].point();
}
/*----------------------------------------------------------------------------*/
math::Segment Edge::segment() const
{
    std::vector<Node> nodes;
    get<Node>(nodes);
    return math::Segment(nodes[0].point(), nodes[1].point());
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");

	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2N)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2E)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2F)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGet(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2R)[m_id].values(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetNodeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetEdgeIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2E)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetFaceIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2F)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetRegionIDs(std::vector<TCellID>&   ACells) const
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id].values(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Node>&   ACells) const
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2N)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Edge>&   ACells) const
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2N)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_edges_container->buildEdge(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Face>&   ACells) const
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2F)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAll(std::vector<Region>& ACells) const
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	std::vector<TCellID> cellIDs;
	(*m_edges_container->m_E2R)[m_id].allValues(cellIDs);
	ACells.resize(cellIDs.size());
	for(unsigned int i=0;i<cellIDs.size();i++){
		ACells[i] =m_owner->m_regions_container->buildRegion(cellIDs[i]);
	}
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2N))
                throw GMDSException("E2N adjacency is not supported by the mesh model");
                
        (*m_edges_container->m_E2N)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2E))
                throw GMDSException("E2E adjacency is not supported by the mesh model");
        (*m_edges_container->m_E2E)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2F))
                throw GMDSException("E2F adjacency is not supported by the mesh model");
        (*m_edges_container->m_E2F)[m_id].allValues(ACells);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(E2R))
                throw GMDSException("E2R adjacency is not supported by the mesh model");
        (*m_edges_container->m_E2R)[m_id].allValues(ACells);
}

/*----------------------------------------------------------------------------*/
void Edge::delegateSetNodeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateSetEdgeIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2E)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateSetFaceIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2F)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateSetRegionIDs(const std::vector<TCellID>&   ACells)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id]=ACells;
}
/*----------------------------------------------------------------------------*/
void Edge::delegateNodeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateEdgeAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2E)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateFaceAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2F)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateRegionAdd(TCellID AID)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id].add(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateNodeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2N)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateEdgeRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2E)[m_id].del(AID);

}
/*----------------------------------------------------------------------------*/
void Edge::delegateFaceRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2F)[m_id].del(AID);

}
/*----------------------------------------------------------------------------*/
void Edge::delegateRegionRemove(TCellID AID)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");
	(*m_edges_container->m_E2R)[m_id].del(AID);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateNodeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2N))
		throw GMDSException("E2N adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2N)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateEdgeReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2E))
		throw GMDSException("E2E adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2E)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateFaceReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2F))
		throw GMDSException("E2F adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2F)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void Edge::delegateRegionReplace(TCellID AID1, TCellID AID2)
{
	if(!m_owner->m_model.has(E2R))
		throw GMDSException("E2R adjacency is not supported by the mesh model");

	(*m_edges_container->m_E2R)[m_id].replace(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
std::ostream & operator << (std::ostream & AStream, const Edge & AN)
{
	AStream<<"Edge "<< AN.id();
	return AStream;
}
/*----------------------------------------------------------------------------*/
Node Edge::getOppositeNode(const Node &ANode) const
{
	 return  m_owner->get<Node>(getOppositeNodeId(ANode.id()));
}
/*----------------------------------------------------------------------------*/
Node Edge::getOppositeNode(const TCellID &ANodeID) const
{
	return m_owner->get<Node>(getOppositeNodeId(ANodeID));

}
/*----------------------------------------------------------------------------*/
TCellID Edge::getOppositeNodeId(const Node &ANode) const
{
	return getOppositeNodeId(ANode.id());
}
/*----------------------------------------------------------------------------*/
TCellID Edge::getOppositeNodeId(const TCellID &AID) const
{
	std::vector<TCellID> node_ids = getIDs<Node>();
	return (node_ids[0]==AID)?node_ids[1]:node_ids[0];
}