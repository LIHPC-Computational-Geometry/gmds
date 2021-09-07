/*----------------------------------------------------------------------------*/
/*
 * FaceContainer.cpp
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/FaceContainer.h>
#include <gmds/ig/Mesh.h>

/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
FaceContainer::FaceContainer( Mesh* AMesh)
 :m_mesh(AMesh),m_model(AMesh->getModel()),
 m_T2N(0),m_T2E(0),m_T2F(0),m_T2R(0),
 m_Q2N(0),m_Q2E(0),m_Q2F(0),m_Q2R(0),
 m_P2N(0),m_P2E(0),m_P2F(0),m_P2R(0),
 m_triangles(0), m_quads(0), m_polygons(0)
{

	setConnectivityContainers();
}
/*----------------------------------------------------------------------------*/
FaceContainer::~FaceContainer()
{
	//TRIANGLES
	if(m_T2N)
		delete m_T2N;
	if(m_T2E)
		delete m_T2E;
	if(m_T2F)
		delete m_T2F;
	if(m_T2R)
		delete m_T2R;

	//QUADS
	if(m_Q2N)
		delete m_Q2N;
	if(m_Q2E)
		delete m_Q2E;
	if(m_Q2F)
		delete m_Q2F;
	if(m_Q2R)
		delete m_Q2R;

	//POLYGONS
	if(m_P2N)
		delete m_P2N;
	if(m_P2E)
		delete m_P2E;
	if(m_P2F)
		delete m_P2F;
	if(m_P2R)
		delete m_P2R;

	if(m_triangles)
		delete m_triangles;
	if(m_quads)
		delete m_quads;
	if(m_polygons)
		delete m_polygons;
}
/*----------------------------------------------------------------------------*/
void FaceContainer::setConnectivityContainers()
{
	if(m_model.has(F2N)){
		m_T2N = new SmartVector<TabCellID<3> >();
		m_Q2N = new SmartVector<TabCellID<4> >();
		m_P2N = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(F2E)){
		m_T2E = new SmartVector<TabCellID<3> >();
		m_Q2E = new SmartVector<TabCellID<4> >();
		m_P2E = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(F2F)){
		m_T2F = new SmartVector<TabCellID<3> >();
		m_Q2F = new SmartVector<TabCellID<4> >();
		m_P2F = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(F2R)){
		m_T2R = new SmartVector<TabCellID<2> >();
		m_Q2R = new SmartVector<TabCellID<2> >();
		m_P2R = new SmartVector<TabCellID<2> >();
	}
	m_triangles = new TAccessor(this,this->m_mesh->getModel());
	m_quads     = new QAccessor(this,this->m_mesh->getModel());
	m_polygons  = new PAccessor(this,this->m_mesh->getModel());
}
/*----------------------------------------------------------------------------*/
void FaceContainer::addConnectivityContainers(const TInt ADim)
{
	/** WARNING
	 * When wa add a new adjacency container, we must define its size. If the
	 * cell was already used, then X2N exists and is filled, so we resize to
	 * X2N->size()
	 */
	if(ADim==0){
		if(m_T2N==0){
			m_T2N = new SmartVector<TabCellID<3> >();
			m_triangles->m_N = m_T2N;
			if(m_triangles->adj_N!=0)
				delete m_triangles->adj_N;
			m_triangles->adj_N = new AdjUpdate<3>(m_T2N);
		}
		if(m_Q2N==0){
			m_Q2N = new SmartVector<TabCellID<4> >();
			m_quads->m_N = m_Q2N;
			if(m_quads->adj_N!=0)
				delete m_quads->adj_N;
			m_quads->adj_N = new AdjUpdate<4>(m_Q2N);
		}
		if(m_P2N==0){
			m_P2N = new SmartVector<TabCellID<size_undef> >();
			m_polygons->m_N = m_P2N;
			if(m_polygons->adj_N!=0)
				delete m_polygons->adj_N;
			m_polygons->adj_N = new AdjUpdate<size_undef>(m_P2N);
		}
	}
	else if (ADim==1){
		if (m_T2E==0){
			m_T2E = new SmartVector<TabCellID<3> >(m_T2N->getMarks());
			m_triangles->m_E = m_T2E;
			if(m_triangles->adj_E!=0)
				delete m_triangles->adj_E;
			m_triangles->adj_E = new AdjUpdate<3>(m_T2E);
		}
		if (m_Q2E==0){
			m_Q2E = new SmartVector<TabCellID<4> >(m_Q2N->getMarks());
			m_quads->m_E = m_Q2E;
			if(m_quads->adj_E!=0)
				delete m_quads->adj_E;
			m_quads->adj_E = new AdjUpdate<4>(m_Q2E);
		}
		if (m_P2E==0){
			m_P2E = new SmartVector<TabCellID<size_undef> >(m_P2N->getMarks());
			m_polygons->m_E = m_P2E;
			if(m_polygons->adj_E!=0)
				delete m_polygons->adj_E;
			m_polygons->adj_E = new AdjUpdate<size_undef>(m_P2E);
		}
	}
	else if (ADim==2){
		if (m_T2F==0){
			m_T2F = new SmartVector<TabCellID<3> >(m_T2N->getMarks());
			m_triangles->m_F = m_T2F;
			if(m_triangles->adj_F!=0)
				delete m_triangles->adj_F;
			m_triangles->adj_F = new AdjUpdate<3>(m_T2F);
		}
		if (m_Q2F==0){
			m_Q2F = new SmartVector<TabCellID<4> >(m_Q2N->getMarks());
			m_quads->m_F = m_Q2F;
			if(m_quads->adj_F!=0)
				delete m_quads->adj_F;
			m_quads->adj_F = new AdjUpdate<4>(m_Q2F);
		}
		if (m_P2F==0){
			m_P2F = new SmartVector<TabCellID<size_undef> >(m_P2N->getMarks());
			m_polygons->m_F = m_P2F;
			if(m_polygons->adj_F!=0)
				delete m_polygons->adj_F;
			m_polygons->adj_F = new AdjUpdate<size_undef>(m_P2F);
		}
	}
	else if (ADim==3)
	{
		if(m_T2R==0){
			m_T2R = new SmartVector<TabCellID<2> >(m_T2N->getMarks());
			m_triangles->m_R = m_T2R;
			if(m_triangles->adj_R!=0)
				delete m_triangles->adj_R;
			m_triangles->adj_R = new AdjUpdate<2>(m_T2R);
		}
		if(m_Q2R==0){
			m_Q2R = new SmartVector<TabCellID<2> >(m_Q2N->getMarks());
			m_quads->m_R = m_Q2R;
			if(m_quads->adj_R!=0)
				delete m_quads->adj_R;
			m_quads->adj_R = new AdjUpdate<2>(m_Q2R);
		}
		if(m_P2R==0){
			m_P2R = new SmartVector<TabCellID<2> >(m_P2N->getMarks());
			m_polygons->m_R = m_P2R;
			if(m_polygons->adj_R!=0)
				delete m_polygons->adj_R;
			m_polygons->adj_R = new AdjUpdate<2>(m_P2R);
		}

	}
}
/*----------------------------------------------------------------------------*/
void FaceContainer::removeConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_T2N!=0){
			delete m_T2N;
			m_T2N=0;
                        if(m_triangles->adj_N != 0) {
                                delete m_triangles->adj_N;
                                m_triangles->adj_N = 0;
                        }
		}
		if(m_Q2N!=0){
			delete m_Q2N;
			m_Q2N=0;
                        if(m_quads->adj_N != 0) {
                                delete m_quads->adj_N;
                                m_quads->adj_N = 0;
                        }
		}
		if(m_P2N!=0){
			delete m_P2N;
			m_P2N=0;
                        if(m_polygons->adj_N != 0) {
                                delete m_polygons->adj_N;
                                m_polygons->adj_N = 0;
                        }
		}
	}
	else if (ADim==1){
		if(m_T2E!=0){
			delete m_T2E;
			m_T2E=0;
                        if(m_triangles->adj_E != 0) {
                                delete m_triangles->adj_E;
                                m_triangles->adj_E = 0;
                        }
		}
		if(m_Q2E!=0){
			delete m_Q2E;
			m_Q2E=0;
			if(m_quads->adj_E != 0) {
				delete m_quads->adj_E;
				m_quads->adj_E = 0;
			}
		}
		if(m_P2E!=0){
			delete m_P2E;
			m_P2E=0;
                        if(m_polygons->adj_E != 0) {
                                delete m_polygons->adj_E;
                                m_polygons->adj_E = 0;
                        }
		}

	}
	else if (ADim==2){
		if(m_T2F!=0){
			delete m_T2F;
			m_T2F=0;
                        if(m_triangles->adj_F != 0) {
                                delete m_triangles->adj_F;
                                m_triangles->adj_F = 0;
                        }
		}
		if(m_Q2F!=0){
			delete m_Q2F;
			m_Q2F=0;
                        if(m_quads->adj_F != 0) {
                                delete m_quads->adj_F;
                                m_quads->adj_F = 0;
                        }
		}
		if(m_P2F!=0){
			delete m_P2F;
			m_P2F=0;
                        if(m_polygons->adj_F != 0) {
                                delete m_polygons->adj_F;
                                m_polygons->adj_F = 0;
                        }
		}
	}
	else if (ADim==3)
	{
		if(m_T2R!=0){
			delete m_T2R;
			m_T2R=0;
                        if(m_triangles->adj_R != 0) {
                                delete m_triangles->adj_R;
                                m_triangles->adj_R = 0;
                        }
		}
		if(m_Q2R!=0){
			delete m_Q2R;
			m_Q2R=0;
                        if(m_quads->adj_R != 0) {
                                delete m_quads->adj_R;
                                m_quads->adj_R = 0;
                        }
		}
		if(m_P2R!=0){
			delete m_P2R;
			m_P2R=0;
                        if(m_polygons->adj_R != 0) {
                                delete m_polygons->adj_R;
                                m_polygons->adj_R = 0;
                        }
		}
	}
}
/*----------------------------------------------------------------------------*/
Face FaceContainer::
addTriangle()
{
	return addTriangle(NullID,NullID,NullID);
}
/*----------------------------------------------------------------------------*/
Face FaceContainer::
addTriangle(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3)
{
	//===========================================
	// STEP 1 - we value the face index
	//===========================================
	TInt index_face = m_face_ids.getFreeIndex();
	return addTriangle(AN1,AN2,AN3,index_face);
}
/*----------------------------------------------------------------------------*/
Face FaceContainer::
addTriangle(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
		const TCellID& AID)
{
	TCellID index_face = AID;

	m_face_ids.assign(index_face);
	//===========================================
	// STEP 2 - we value the triangle index
	//===========================================
	int index_type= m_triangles->getID();
	//============================================
	// STEP 3 - we build and assign the face infos
	//============================================
	FaceInfo info(GMDS_TRIANGLE,index_type);
	m_face_types.assign(info, index_face);

	//============================================
	// STEP 4 - we fill in F2N
	//============================================
	if(m_model.has(F2N)){
		//here only nodes are updated
		TabCellID<3> t;
		t.add(AN1);
		t.add(AN2);
		t.add(AN3);
		m_T2N->assign(t,index_type);
	}
	if(m_model.has(F2E)){
		TabCellID<3> t;
		m_T2E->assign(t,index_type);
	}
	if(m_model.has(F2F)){
		TabCellID<3> t;
		m_T2F->assign(t,index_type);
	}
	if(m_model.has(F2R)){
		TabCellID<2> t;
		m_T2R->assign(t,index_type);
	}


	//===============================================
	// STEP 5 - Temporary (or handle) face is created
	//===============================================
	Face f(m_mesh,GMDS_TRIANGLE,index_face);

	return f;
}
/*----------------------------------------------------------------------------*/
Face FaceContainer::
addQuad(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
		const TCellID& AN4)
{
	//===========================================
	// STEP 1 - we value the face index
	//===========================================
	TInt index_face = m_face_ids.getFreeIndex();
	return addQuad(AN1,AN2,AN3,AN4,index_face);

}
/*----------------------------------------------------------------------------*/
Face FaceContainer::
addQuad(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
		const TCellID& AN4, const TCellID& AID)
{
	TInt index_face = AID;
	m_face_ids.assign(index_face);
	//value the quad index
	int index_type= m_quads->getID();
	// Face Info creation
	FaceInfo info(GMDS_QUAD,index_type);
	m_face_types.assign(info, index_face);

	//nodes info must be updated


	if(m_model.has(F2N)){
		//here only nodes are updated
		TabCellID<4> t;
		t.add(AN1);
		t.add(AN2);
		t.add(AN3);
		t.add(AN4);
		m_Q2N->assign(t,index_type);
	}
	if(m_model.has(F2E)){
		TabCellID<4> t;
		m_Q2E->assign(t,index_type);
	}
	if(m_model.has(F2F)){
		TabCellID<4> t;
		m_Q2F->assign(t,index_type);
	}
	if(m_model.has(F2R)){
		TabCellID<2> t;
		m_Q2R->assign(t,index_type);
	}

	Face f(m_mesh,GMDS_QUAD,index_face);
	return f;
}
/*----------------------------------------------------------------------------*/
Face FaceContainer::
addPolygon(const std::vector<TCellID>& ANodes)
{
	//===========================================
	// STEP 1 - we value the face index
	//===========================================
	TCellID index_face = m_face_ids.getFreeIndex();
	return addPolygon(ANodes,index_face);
}
/*----------------------------------------------------------------------------*/
Face FaceContainer::
addPolygon(const std::vector<TCellID>& ANodes, const TCellID& AID)
{
	TCellID index_face = AID;
	m_face_ids.assign(index_face);
	//===========================================
	// STEP 2 - we value the polygon index
	//===========================================
	int index_type= m_polygons->getID();
	//============================================
	// STEP 3 - we build and assign the face infos
	//============================================
	FaceInfo info(GMDS_POLYGON,index_type);
	m_face_types.assign(info, index_face);

	//============================================
	// STEP 4 - we fill in F2N
	//============================================
	if(m_model.has(F2N)){
		//here only nodes are updated

		TabCellID<size_undef> t;
		t=ANodes;
		m_P2N->assign(t,index_type);
	}
	if(m_model.has(F2E)){
		TabCellID<size_undef> t;
		m_P2E->assign(t,index_type);
	}
	if(m_model.has(F2F)){
		TabCellID<size_undef> t;
		m_P2F->assign(t,index_type);
	}
	if(m_model.has(F2R)){
		TabCellID<2> t;
		m_P2R->assign(t,index_type);
	}

	//===============================================
	// STEP 5 - Temporary (or handle) face is created
	//===============================================
	Face f(m_mesh,GMDS_POLYGON,index_face);

	return f;
}
/*----------------------------------------------------------------------------*/
Face FaceContainer::buildFace(const TInt AIndex) const
{
	ECellType type = m_face_types[AIndex].type;
	Face f(m_mesh,type,AIndex);
	return f;
}
/*----------------------------------------------------------------------------*/
void FaceContainer::
getNodesData(const TCellID& AID, int& ANbNodes) const
{
	FaceInfo info = m_face_types[AID];
	if(info.type==GMDS_TRIANGLE)
	{
		if( m_triangles->adj_N==0)
			ANbNodes=0;
		else
			ANbNodes = m_triangles->adj_N->getSizeOfAnElement();
	}
	else if(info.type==GMDS_QUAD)
	{
		if( m_quads->adj_N==0)
			ANbNodes=0;
		else
			ANbNodes = m_quads->adj_N->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYGON)
	{
		ANbNodes = (*m_polygons->m_N)[info.type_id].size();
	}
}
/*----------------------------------------------------------------------------*/
void FaceContainer::
getEdgesData(const TCellID& AID, int& ANbEdges) const
{
	FaceInfo info = m_face_types[AID];
	if(info.type==GMDS_TRIANGLE)
	{
		if( m_triangles->adj_E==0)
			ANbEdges=0;
		else
			ANbEdges = m_triangles->adj_E->getSizeOfAnElement();
	}
	else if(info.type==GMDS_QUAD)
	{
		if( m_quads->adj_E==0)
			ANbEdges=0;
		else
			ANbEdges = m_quads->adj_E->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYGON)
	{
		ANbEdges = (*m_polygons->m_E)[info.type_id].size();
	}
}
/*----------------------------------------------------------------------------*/
void FaceContainer::
getFacesData(const TCellID& AID, int& ANbFaces) const
{
	FaceInfo info = m_face_types[AID];
	if(info.type==GMDS_TRIANGLE)
	{
		if( m_triangles->adj_F==0)
			ANbFaces=0;
		else
		ANbFaces = m_triangles->adj_F->getSizeOfAnElement();
	}
	else if(info.type==GMDS_QUAD)
	{
		if( m_quads->adj_F==0)
			ANbFaces=0;
		else
			ANbFaces = m_quads->adj_F->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYGON)
	{
		ANbFaces = (*m_polygons->m_F)[info.type_id].size();
	}
}
/*----------------------------------------------------------------------------*/
void FaceContainer::
getRegionsData(const TCellID& AID, int& ANbRegions) const
{
	FaceInfo info = m_face_types[AID];
	if(info.type==GMDS_TRIANGLE)
	{
		if( m_triangles->adj_R==0)
			ANbRegions=0;
		else
			ANbRegions = m_triangles->adj_R->getSizeOfAnElement();
	}
	else if(info.type==GMDS_QUAD)
	{
		if( m_quads->adj_R==0)
			ANbRegions=0;
		else
			ANbRegions = m_quads->adj_R->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYGON)
	{
		ANbRegions = 2;
	}
}
/*----------------------------------------------------------------------------*/
void FaceContainer::remove(TInt AIndex)
{
	m_face_ids.unselect(AIndex);
	if(m_face_types[AIndex].type==GMDS_TRIANGLE)
	{
		TCellID type_index = m_face_types[AIndex].type_id;
		if(m_model.has(F2N))
			m_T2N->remove(type_index);
		if(m_model.has(F2E))
			m_T2E->remove(type_index);
		if(m_model.has(F2F))
			m_T2F->remove(type_index);
		if(m_model.has(F2R))
			m_T2R->remove(type_index);

	}
	else if(m_face_types[AIndex].type==GMDS_QUAD)
	{
		TCellID type_index = m_face_types[AIndex].type_id;
		if(m_model.has(F2N))
			m_Q2N->remove(type_index);
		if(m_model.has(F2E))
			m_Q2E->remove(type_index);
		if(m_model.has(F2F))
			m_Q2F->remove(type_index);
		if(m_model.has(F2R))
			m_Q2R->remove(type_index);

	}
	else //POLYGON
	{
		TCellID type_index = m_face_types[AIndex].type_id;
		if(m_model.has(F2N))
			m_P2N->remove(type_index);
		if(m_model.has(F2E))
			m_P2E->remove(type_index);
		if(m_model.has(F2F))
			m_P2F->remove(type_index);
		if(m_model.has(F2R))
			m_P2R->remove(type_index);
	}
}
/*----------------------------------------------------------------------------*/
bool FaceContainer::has(const TCellID& AID) const
{
	return m_face_ids[AID];
}

/*----------------------------------------------------------------------------*/
void FaceContainer::clear()
{
	m_face_ids.clear();
	m_face_types.clear();

	if(m_T2N)
		m_T2N->clear();
	if(m_T2E)
		m_T2E->clear();
	if(m_T2F)
		m_T2F->clear();
	if(m_T2R)
		m_T2R->clear();

	if(m_Q2N)
		m_Q2N->clear();
	if(m_Q2E)
		m_Q2E->clear();
	if(m_Q2F)
		m_Q2F->clear();
	if(m_Q2R)
		m_Q2R->clear();

	if(m_P2N)
		m_P2N->clear();
	if(m_P2E)
		m_P2E->clear();
	if(m_P2F)
		m_P2F->clear();
	if(m_P2R)
		m_P2R->clear();

}
/*----------------------------------------------------------------------------*/
void FaceContainer::resize(const TInt ASize)
{
	m_face_ids.resize(ASize);
	m_face_types.resize(ASize);
}
/*----------------------------------------------------------------------------*/
void FaceContainer::update()
{
	m_face_ids.update();
}
/*----------------------------------------------------------------------------*/
void FaceContainer::serialize(std::ostream& AStr)
{
	m_face_ids.serialize(AStr);
	m_face_types.serialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.serialize(AStr);
#endif //GMDS_PARALLEL

	if(m_T2N)
		m_T2N->serialize(AStr);
	if(m_T2E)
		m_T2E->serialize(AStr);
	if(m_T2F)
		m_T2F->serialize(AStr);
	if(m_T2R)
		m_T2R->serialize(AStr);

	if(m_Q2N)
		m_Q2N->serialize(AStr);
	if(m_Q2E)
		m_Q2E->serialize(AStr);
	if(m_Q2F)
		m_Q2F->serialize(AStr);
	if(m_Q2R)
		m_Q2R->serialize(AStr);

	if(m_P2N)
		m_P2N->serialize(AStr);
	if(m_P2E)
		m_P2E->serialize(AStr);
	if(m_P2F)
		m_P2F->serialize(AStr);
	if(m_P2R)
		m_P2R->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void FaceContainer::unserialize(std::istream& AStr)
{
	m_face_ids.unserialize(AStr);
	m_face_types.unserialize(AStr);
#ifdef GMDS_PARALLEL
	m_distributed_data.unserialize(AStr);
#endif //GMDS_PARALLEL

	if(m_T2N)
		m_T2N->unserialize(AStr);
	if(m_T2E)
		m_T2E->unserialize(AStr);
	if(m_T2F)
		m_T2F->unserialize(AStr);
	if(m_T2R)
		m_T2R->unserialize(AStr);

	if(m_Q2N)
		m_Q2N->unserialize(AStr);
	if(m_Q2E)
		m_Q2E->unserialize(AStr);
	if(m_Q2F)
		m_Q2F->unserialize(AStr);
	if(m_Q2R)
		m_Q2R->unserialize(AStr);

	if(m_P2N)
		m_P2N->unserialize(AStr);
	if(m_P2E)
		m_P2E->unserialize(AStr);
	if(m_P2F)
		m_P2F->unserialize(AStr);
	if(m_P2R)
		m_P2R->unserialize(AStr);
}

/*----------------------------------------------------------------------------*/
FaceContainer::TAccessor::TAccessor(FaceContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(F2N))?AOwner->m_T2N:0;
	m_E =(AOwner->m_model.has(F2E))?AOwner->m_T2E:0;
	m_F =(AOwner->m_model.has(F2F))?AOwner->m_T2F:0;
	m_R =(AOwner->m_model.has(F2R))?AOwner->m_T2R:0;

	adj_N = (AOwner->m_model.has(F2N))? new AdjUpdate<3>(m_N):0;
	adj_E = (AOwner->m_model.has(F2E))? new AdjUpdate<3>(m_E):0;
	adj_F = (AOwner->m_model.has(F2F))? new AdjUpdate<3>(m_F):0;
	adj_R = (AOwner->m_model.has(F2R))? new AdjUpdate<2>(m_R):0;
}
/*----------------------------------------------------------------------------*/
FaceContainer::TAccessor::~TAccessor()
{
	if (adj_N)
		delete adj_N;
	if(adj_E)
		delete adj_E;
	if(adj_F)
		delete adj_F;
	if(adj_R)
		delete adj_R;
}

/*----------------------------------------------------------------------------*/
TInt FaceContainer::TAccessor::getID()
{
	TInt i = NullID;
	if (adj_N)
		i=adj_N->select();
	if(adj_E)
		i=adj_E->select();
	if(adj_F)
		i=adj_F->select();
	if(adj_R)
		i=adj_R->select();
	return i;
}
/*----------------------------------------------------------------------------*/
FaceContainer::QAccessor::QAccessor(FaceContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(F2N))?AOwner->m_Q2N:0;
	m_E =(AOwner->m_model.has(F2E))?AOwner->m_Q2E:0;
	m_F =(AOwner->m_model.has(F2F))?AOwner->m_Q2F:0;
	m_R =(AOwner->m_model.has(F2R))?AOwner->m_Q2R:0;

	adj_N = (AOwner->m_model.has(F2N))? new AdjUpdate<4>(m_N):0;
	adj_E = (AOwner->m_model.has(F2E))? new AdjUpdate<4>(m_E):0;
	adj_F = (AOwner->m_model.has(F2F))? new AdjUpdate<4>(m_F):0;
	adj_R = (AOwner->m_model.has(F2R))? new AdjUpdate<2>(m_R):0;
}
/*----------------------------------------------------------------------------*/
FaceContainer::QAccessor::~QAccessor()
{
	if (adj_N)
		delete adj_N;
	if(adj_E)
		delete adj_E;
	if(adj_F)
		delete adj_F;
	if(adj_R)
		delete adj_R;
}
/*----------------------------------------------------------------------------*/
TInt FaceContainer::QAccessor::getID()
{
	TInt i = NullID;
	if (adj_N)
		i=adj_N->select();
	if(adj_E)
		i=adj_E->select();
	if(adj_F)
		i=adj_F->select();
	if(adj_R)
		i=adj_R->select();
	return i;
}
/*----------------------------------------------------------------------------*/
FaceContainer::PAccessor::PAccessor(FaceContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(F2N))?AOwner->m_P2N:0;
	m_E =(AOwner->m_model.has(F2E))?AOwner->m_P2E:0;
	m_F =(AOwner->m_model.has(F2F))?AOwner->m_P2F:0;
	m_R =(AOwner->m_model.has(F2R))?AOwner->m_P2R:0;

	adj_N = (AOwner->m_model.has(F2N))? new AdjUpdate<size_undef>(m_N):0;
	adj_E = (AOwner->m_model.has(F2E))? new AdjUpdate<size_undef>(m_E):0;
	adj_F = (AOwner->m_model.has(F2F))? new AdjUpdate<size_undef>(m_F):0;
	adj_R = (AOwner->m_model.has(F2R))? new AdjUpdate<2>(m_R):0;
}
/*----------------------------------------------------------------------------*/
FaceContainer::PAccessor::~PAccessor()
{
	if (adj_N)
		delete adj_N;
	if(adj_E)
		delete adj_E;
	if(adj_F)
		delete adj_F;
	if(adj_R)
		delete adj_R;
}

/*----------------------------------------------------------------------------*/
TInt FaceContainer::PAccessor::getID()
{
	TInt i = NullID;
	if (adj_N)
		i=adj_N->select();
	if(adj_E)
		i=adj_E->select();
	if(adj_F)
		i=adj_F->select();
	if(adj_R)
		i=adj_R->select();
	return i;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
