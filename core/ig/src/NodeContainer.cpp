/*----------------------------------------------------------------------------*/
/*
 * NodeContainer.cpp
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/NodeContainer.h>
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
NodeContainer::NodeContainer( Mesh* AMesh)
:m_mesh(AMesh),m_model(AMesh->getModel())
{
	if(m_model.has(N2N))
		m_N2N = new IndexedVector<TabCellID<size_undef> >();
	else
		m_N2N = nullptr;
	if(m_model.has(N2E))
		m_N2E = new IndexedVector<TabCellID<size_undef> >();
	else
		m_N2E = nullptr;
	if(m_model.has(N2F))
		m_N2F = new IndexedVector<TabCellID<size_undef> >();
	else
		m_N2F = nullptr;
	if(m_model.has(N2R))
		m_N2R = new IndexedVector<TabCellID<size_undef> >();
	else
		m_N2R = nullptr;
}
/*----------------------------------------------------------------------------*/
NodeContainer::~NodeContainer()
{

		delete m_N2N;

		delete m_N2E;

		delete m_N2F;

		delete m_N2R;

    m_N2N=nullptr;
    m_N2E=nullptr;
    m_N2F=nullptr;
    m_N2R=nullptr;
}
/*----------------------------------------------------------------------------*/
Node NodeContainer::add(const TCoord& AX, const TCoord& AY, const TCoord& AZ)
{
	TInt index = m_node_ids.getFreeIndex();
	return add(AX,AY,AZ,index);
}
/*----------------------------------------------------------------------------*/
Node NodeContainer::add(const TCoord& AX, const TCoord& AY, const TCoord& AZ,
		const TCellID& AGID)
{
	TCellID index = AGID;
	m_node_ids.assign(index);

	math::Point p(AX,AY,AZ);
	m_node_coords.assign(p, index);

	if(m_model.has(N2N)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		m_N2N->assign(t,index);
	}
	if(m_model.has(N2E)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		m_N2E->assign(t,index);
	}
	if(m_model.has(N2F)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		m_N2F->assign(t,index);
	}
	if(m_model.has(N2R)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		m_N2R->assign(t,index);
	}

	Node n(m_mesh,index,p);
	return n;
}
/*----------------------------------------------------------------------------*/
void NodeContainer::addConnectivityContainers(const TInt ADim)
{
	/** WARNING
	 * When wa add a new adjacency container, we must define its size. If the
	 * cell was already used, then X2N exists and is filled, so we resize to
	 * X2N->size()
	 */
	if(ADim==0){
		if(m_N2N==nullptr){
			m_N2N = new IndexedVector<TabCellID<size_undef> >(
				m_node_ids.capacity());
			m_node_ids.update();
		}
	}
	else if (ADim==1){
		if (m_N2E==nullptr){
                        m_N2E = new IndexedVector<TabCellID<size_undef> >(
                                m_node_ids.capacity());
			m_node_ids.update();
		}
	}
	else if (ADim==2){
		if (m_N2F==nullptr){
                        m_N2F = new IndexedVector<TabCellID<size_undef> >(
                                m_node_ids.capacity());
			m_node_ids.update();
		}
	}
	else if (ADim==3){
		if(m_N2R==nullptr){
                        m_N2R = new IndexedVector<TabCellID<size_undef> >(
                                m_node_ids.capacity());
			m_node_ids.update();
		}
	}
}
/*----------------------------------------------------------------------------*/
void NodeContainer::removeConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_N2N!=nullptr){
			delete m_N2N;
			m_N2N=nullptr;
		}
	}
	else if (ADim==1){
		if(m_N2E!=nullptr){
			delete m_N2E;
			m_N2E=nullptr;
		}
	}
	else if (ADim==2){
		if(m_N2F!=nullptr){
			delete m_N2F;
			m_N2F=nullptr;
		}
	}
	else if (ADim==3)
	{
		if(m_N2R!=nullptr){
			delete m_N2R;
			m_N2R=nullptr;
		}
	}
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getNodesData(const TCellID& AID, int& ANbNodes) const
{
	if(m_N2N==nullptr)
		ANbNodes = 0;
	else
		ANbNodes = (*m_N2N)[AID].size();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getEdgesData(const TCellID& AID, int& ANbEdges) const
{
	if(m_N2E==nullptr)
		ANbEdges = 0;
	else
		ANbEdges = (*m_N2E)[AID].size();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getFacesData(const TCellID& AID, int& ANbFaces) const
{
	if(m_N2F==nullptr)
		ANbFaces = 0;
	else
		ANbFaces = (*m_N2F)[AID].size();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::getRegionsData(const TCellID& AID, int& ANbRegions) const
{
	if(m_N2R==nullptr)
		ANbRegions = 0;
	else
		ANbRegions = (*m_N2R)[AID].size();
}
/*----------------------------------------------------------------------------*/
Node NodeContainer::buildNode(const TInt index) const
{
	math::Point p = m_node_coords[index];
	Node n(m_mesh,index,p);
	return n;
}
/*----------------------------------------------------------------------------*/
bool NodeContainer::has(const TCellID& AID) const
{
	return m_node_ids[AID];
}

/*----------------------------------------------------------------------------*/
void NodeContainer::update()
{
	m_node_ids.update();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::clear()
{
	m_node_ids.clear();
	m_node_coords.clear();
	if(m_N2N)
		m_N2N->clear();
	if(m_N2E)
		m_N2E->clear();
	if(m_N2F)
		m_N2F->clear();
	if(m_N2R)
		m_N2R->clear();
}
/*----------------------------------------------------------------------------*/
void NodeContainer::resize(const TInt ASize)
{
	m_node_ids.resize(ASize);
	m_node_coords.resize(ASize);
}
/*----------------------------------------------------------------------------*/
void NodeContainer::serialize(std::ostream& AStr)
{
	m_node_ids.serialize(AStr);
	m_node_coords.serialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.serialize(AStr);
#endif //GMDS_PARALLEL

	if(m_N2N)
		m_N2N->serialize(AStr);
	if(m_N2E)
		m_N2E->serialize(AStr);
	if(m_N2F)
		m_N2F->serialize(AStr);
	if(m_N2R)
		m_N2R->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void NodeContainer::unserialize(std::istream& AStr)
{
	clear();
	m_node_ids.unserialize(AStr);
	m_node_coords.unserialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.unserialize(AStr);
#endif //GMDS_PARALLEL

	if(m_N2N)
		m_N2N->unserialize(AStr);
	if(m_N2E)
		m_N2E->unserialize(AStr);
	if(m_N2F)
		m_N2F->unserialize(AStr);
	if(m_N2R)
		m_N2R->unserialize(AStr);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
