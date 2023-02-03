/*----------------------------------------------------------------------------*/
/*
 * EdgeContainer.cpp
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/EdgeContainer.h>
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
EdgeContainer::EdgeContainer( Mesh* AMesh)
:m_mesh(AMesh),m_model(AMesh->getModel()),
 m_E2N(nullptr),m_E2E(nullptr),m_E2F(nullptr),m_E2R(nullptr)
{
	if(m_model.has(E2N))
	{
		m_E2N = new IndexedVector<TabCellID<2> >();
	}
	if(m_model.has(E2E)){
		m_E2E = new IndexedVector<TabCellID<size_undef> >();
	}
	if(m_model.has(E2F)){
		m_E2F = new IndexedVector<TabCellID<size_undef> >();
	}
	if(m_model.has(E2R)){
		m_E2R = new IndexedVector<TabCellID<size_undef> >();
	}
}
/*----------------------------------------------------------------------------*/
EdgeContainer::~EdgeContainer()
{

		delete m_E2N;

		delete m_E2E;

		delete m_E2F;

		delete m_E2R;
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::addConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_E2N==nullptr){
			m_E2N = new IndexedVector<TabCellID<2> >(m_edge_ids.capacity());
			m_edge_ids.update();
		}
	}
	else if (ADim==1){
		if (m_E2E==nullptr){
			m_E2E = new IndexedVector<TabCellID<size_undef> >(m_edge_ids.capacity());
			m_edge_ids.update();
		}
	}
	else if (ADim==2){
		if (m_E2F==nullptr){
			m_E2F = new IndexedVector<TabCellID<size_undef> >(m_edge_ids.capacity());
			m_edge_ids.update();
		}
	}
	else if (ADim==3){
		if(m_E2R==nullptr){
			m_E2R = new IndexedVector<TabCellID<size_undef> >(m_edge_ids.capacity());
			m_edge_ids.update();
		}
	}
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::removeConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_E2N!=nullptr){
			delete m_E2N;
			m_E2N=nullptr;
		}
	}
	else if (ADim==1){
		if(m_E2E!=nullptr){
			delete m_E2E;
			m_E2E=nullptr;
		}
	}
	else if (ADim==2){
		if(m_E2F!=nullptr){
			delete m_E2F;
			m_E2F=nullptr;
		}
	}
	else if (ADim==3)
	{
		if(m_E2R!=nullptr){
			delete m_E2R;
			m_E2R=nullptr;
		}
	}
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getNodesData(const TCellID& AID, TInt& ANbNodes) const
{
	if(m_E2N==nullptr)
		ANbNodes = 0;
	else
		ANbNodes = 2;
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getEdgesData(const TCellID& AID, TInt& ANbEdges) const
{
	if(m_E2E==nullptr)
		ANbEdges = 0;
	else
		ANbEdges = static_cast<int>((*m_E2E)[AID].size());
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getFacesData(const TCellID& AID, TInt& ANbFaces) const
{
	if(m_E2F==nullptr)
		ANbFaces = 0;
	else
		ANbFaces = (*m_E2F)[AID].size();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::getRegionsData(const TCellID& AID, TInt& ANbRegions) const
{
	if(m_E2R==nullptr)
		ANbRegions = 0;
	else
		ANbRegions = (*m_E2R)[AID].size();
}

/*----------------------------------------------------------------------------*/
Edge EdgeContainer::add(const TCellID& AN1,const TCellID& AN2)
{
	//===========================================
	// STEP 1 - we value the edgeindex
	//===========================================
	TCellID index = m_edge_ids.getFreeIndex();
	return add(AN1,AN2,index);
}
/*----------------------------------------------------------------------------*/
bool EdgeContainer::has(const TCellID& AID) const
{
	return m_edge_ids[AID];
}
/*----------------------------------------------------------------------------*/
Edge EdgeContainer::
add(const TCellID& AN1,const TCellID& AN2, const TCellID& AID)
{
	TCellID index = AID;
	m_edge_ids.assign(index);
	//============================================
	// STEP 2 - we fill in F2N
	//============================================
	if(m_model.has(E2N)){
		//here only nodes are updated
		TabCellID<2> t;
		t.add(AN1);
		t.add(AN2);
		m_E2N->assign(t,index);
	}
	if(m_model.has(E2E)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		m_E2E->assign(t,index);
	}
	if(m_model.has(E2F)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		m_E2F->assign(t,index);
	}
	if(m_model.has(E2R)){
		//here only nodes are updated
		TabCellID<size_undef> t;
		m_E2R->assign(t,index);
	}

	//============================================
	// STEP 3 - we create an edge object
	//============================================
	Edge e(m_mesh,index);
	return e;
}
/*----------------------------------------------------------------------------*/
Edge EdgeContainer::buildEdge(const TInt index) const
{
	Edge e(m_mesh,index);
	return e;
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::clear()
{
	m_edge_ids.clear();
	if(m_E2N)
		m_E2N->clear();
	if(m_E2E)
		m_E2E->clear();
	if(m_E2F)
		m_E2F->clear();
	if(m_E2R)
		m_E2R->clear();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::resize(const TInt ASize)
{
	m_edge_ids.resize(ASize);
}

/*----------------------------------------------------------------------------*/
void EdgeContainer::update()
{
	m_edge_ids.update();
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::serialize(std::ostream& AStr)
{
	m_edge_ids.serialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.serialize(AStr);
#endif //GMDS_PARALLEL

	if(m_E2N)
		m_E2N->serialize(AStr);
	if(m_E2E)
		m_E2E->serialize(AStr);
	if(m_E2F)
		m_E2F->serialize(AStr);
	if(m_E2R)
		m_E2R->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void EdgeContainer::unserialize(std::istream& AStr)
{
	m_edge_ids.unserialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.unserialize(AStr);
#endif //GMDS_PARALLEL

	if(m_E2N)
		m_E2N->unserialize(AStr);
	if(m_E2E)
		m_E2E->unserialize(AStr);
	if(m_E2F)
		m_E2F->unserialize(AStr);
	if(m_E2R)
		m_E2R->unserialize(AStr);

}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
