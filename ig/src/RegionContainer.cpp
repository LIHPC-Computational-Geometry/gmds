/*----------------------------------------------------------------------------*/
/*
 * RegionContainer.cpp
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/RegionContainer.h>
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
RegionContainer::RegionContainer( Mesh* AMesh)
 :m_mesh(AMesh),m_model(AMesh->getModel()),
 m_T2N(0),m_T2E(0),m_T2F(0),m_T2R(0),
 m_PY2N(0),m_PY2E(0),m_PY2F(0),m_PY2R(0),
 m_PR2N(0),m_PR2E(0),m_PR2F(0),m_PR2R(0),
 m_H2N(0),m_H2E(0),m_H2F(0),m_H2R(0),
 m_P2N(0),m_P2E(0),m_P2F(0),m_P2R(0),
 m_tet(0), m_pyra(0), m_hex(0), m_poly(0)
{

	setConnectivityContainers();
}
/*----------------------------------------------------------------------------*/
RegionContainer::~RegionContainer()
{
	//TETRA
	if(m_T2N)
		delete m_T2N;
	if(m_T2E)
		delete m_T2E;
	if(m_T2F)
		delete m_T2F;
	if(m_T2R)
		delete m_T2R;

	//PYRAMID
	if(m_PY2N)
		delete m_PY2N;
	if(m_PY2E)
		delete m_PY2E;
	if(m_PY2F)
		delete m_PY2F;
	if(m_PY2R)
		delete m_PY2R;

	//PRISM3
	if(m_PR2N)
		delete m_PR2N;
	if(m_PR2E)
		delete m_PR2E;
	if(m_PR2F)
		delete m_PR2F;
	if(m_PR2R)
		delete m_PR2R;

	//HEX
	if(m_H2N)
		delete m_H2N;
	if(m_H2E)
		delete m_H2E;
	if(m_H2F)
		delete m_H2F;
	if(m_H2R)
		delete m_H2R;

	//POLY
	if(m_P2N)
		delete m_P2N;
	if(m_P2E)
		delete m_P2E;
	if(m_P2F)
		delete m_P2F;
	if(m_P2R)
		delete m_P2R;

	if(m_tet)
		delete m_tet;
	if(m_pyra)
		delete m_pyra;
	if(m_prism3)
		delete m_prism3;
	if(m_hex)
		delete m_hex;
	if(m_poly)
		delete m_poly;
}
/*----------------------------------------------------------------------------*/
void RegionContainer::setConnectivityContainers()
{
	if(m_model.has(R2N)){
		m_T2N = new SmartVector<TabCellID<4> >();
		m_PY2N = new SmartVector<TabCellID<5> >();
		m_PR2N = new SmartVector<TabCellID<6> >();
		m_H2N = new SmartVector<TabCellID<8> >();
		m_P2N = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(R2E)){
		m_T2E = new SmartVector<TabCellID<6> >();
		m_PY2E = new SmartVector<TabCellID<8> >();
		m_PR2E = new SmartVector<TabCellID<9> >();
		m_H2E = new SmartVector<TabCellID<12> >();
		m_P2E = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(R2F)){
		m_T2F = new SmartVector<TabCellID<4> >();
		m_PY2F = new SmartVector<TabCellID<5> >();
		m_PR2F = new SmartVector<TabCellID<5> >();
		m_H2F = new SmartVector<TabCellID<6> >();
		m_P2F = new SmartVector<TabCellID<size_undef> >();
	}
	if(m_model.has(R2R)){
		m_T2R = new SmartVector<TabCellID<4> >();
		m_PY2R = new SmartVector<TabCellID<5> >();
		m_PR2R = new SmartVector<TabCellID<5> >();
		m_H2R = new SmartVector<TabCellID<6> >();
		m_P2R = new SmartVector<TabCellID<size_undef> >();
	}
	m_tet    = new TAccessor(this,this->m_mesh->getModel());
	m_pyra   = new PyAccessor(this,this->m_mesh->getModel());
	m_prism3 = new PrAccessor(this,this->m_mesh->getModel());
	m_hex    = new HAccessor(this,this->m_mesh->getModel());
	m_poly   = new PAccessor(this,this->m_mesh->getModel());
}
/*----------------------------------------------------------------------------*/
void RegionContainer::addConnectivityContainers(const TInt ADim)
{
	/** WARNING
	 * When wa add a new adjacency container, we must define its size. If the
	 * cell was already used, then X2N exists and is filled, so we resize to
	 * X2N->size()
	 */
	if(ADim==0){
		if(m_T2N==0){
			m_T2N = new SmartVector<TabCellID<4> >();
			m_tet->m_N = m_T2N;
			if(m_tet->adj_N!=0)
				delete m_tet->adj_N;
			m_tet->adj_N = new AdjUpdate<4>(m_T2N);
		}
		if(m_PY2N==0){
			m_PY2N = new SmartVector<TabCellID<5> >();
			m_pyra->m_N = m_PY2N;
			if(m_pyra->adj_N!=0)
				delete m_pyra->adj_N;
			m_pyra->adj_N = new AdjUpdate<5>(m_PY2N);
		}
		if(m_PR2N==0){
			m_PR2N = new SmartVector<TabCellID<6> >();
			m_prism3->m_N = m_PR2N;
			if(m_prism3->adj_N!=0)
				delete m_prism3->adj_N;
			m_prism3->adj_N = new AdjUpdate<6>(m_PR2N);
		}
		if(m_H2N==0){
			m_H2N = new SmartVector<TabCellID<8> >();
			m_hex->m_N = m_H2N;
			if(m_hex->adj_N!=0)
				delete m_hex->adj_N;
			m_hex->adj_N = new AdjUpdate<8>(m_H2N);
		}
		if(m_P2N==0){
			m_P2N = new SmartVector<TabCellID<size_undef> >();
			m_poly->m_N = m_P2N;
			if(m_poly->adj_N!=0)
				delete m_poly->adj_N;
			m_poly->adj_N = new AdjUpdate<size_undef>(m_P2N);
		}
	}
	else if (ADim==1){
		if (m_T2E==0){
			m_T2E = new SmartVector<TabCellID<6> >(m_T2N->getMarks());
			m_tet->m_E = m_T2E;
			if(m_tet->adj_E!=0)
				delete m_tet->adj_E;
			m_tet->adj_E = new AdjUpdate<6>(m_T2E);
		}
		if(m_PY2E==0){
			m_PY2E = new SmartVector<TabCellID<8> >(m_PY2N->getMarks());
			m_pyra->m_E = m_PY2E;
			if(m_pyra->adj_E!=0)
				delete m_pyra->adj_E;
			m_pyra->adj_E = new AdjUpdate<8>(m_PY2E);
		}

		if(m_PR2E==0){
			m_PR2E = new SmartVector<TabCellID<9> >(m_PR2N->getMarks());
			m_prism3->m_E = m_PR2E;
			if(m_prism3->adj_E!=0)
				delete m_prism3->adj_E;
			m_prism3->adj_E = new AdjUpdate<9>(m_PR2E);
		}

		if (m_H2E==0){
			m_H2E = new SmartVector<TabCellID<12> >(m_H2N->getMarks());
			m_hex->m_E = m_H2E;
			if(m_hex->adj_E!=0)
				delete m_hex->adj_E;
			m_hex->adj_E = new AdjUpdate<12>(m_H2E);
		}
		if (m_P2E==0){
			m_P2E = new SmartVector<TabCellID<size_undef> >(m_P2N->getMarks());
			m_poly->m_E = m_P2E;
			if(m_poly->adj_E!=0)
				delete m_poly->adj_E;
			m_poly->adj_E = new AdjUpdate<size_undef>(m_P2E);
		}
	}
	else if (ADim==2){
		if (m_T2F==0){
			m_T2F = new SmartVector<TabCellID<4> >(m_T2N->getMarks());
			m_tet->m_F = m_T2F;
			if(m_tet->adj_F!=0)
				delete m_tet->adj_F;
			m_tet->adj_F = new AdjUpdate<4>(m_T2F);
		}
		if(m_PY2F==0){
			m_PY2F = new SmartVector<TabCellID<5> >(m_PY2N->getMarks());
			m_pyra->m_F = m_PY2F;
			if(m_pyra->adj_F!=0)
			delete m_pyra->adj_F;
			m_pyra->adj_F = new AdjUpdate<5>(m_PY2F);
		}
		if(m_PR2F==0){
			m_PR2F = new SmartVector<TabCellID<5> >(m_PR2N->getMarks());
			m_prism3->m_F = m_PR2F;
			if(m_prism3->adj_F!=0)
			delete m_prism3->adj_F;
			m_prism3->adj_F = new AdjUpdate<5>(m_PR2F);
		}
		if (m_H2F==0){
			m_H2F = new SmartVector<TabCellID<6> >(m_H2N->getMarks());
			m_hex->m_F = m_H2F;
			if(m_hex->adj_F!=0)
				delete m_hex->adj_F;
			m_hex->adj_F = new AdjUpdate<6>(m_H2F);
		}
		if (m_P2F==0){
			m_P2F = new SmartVector<TabCellID<size_undef> >(m_P2N->getMarks());
			m_P2F->resize(m_P2N->size());
			m_poly->m_F = m_P2F;
			if(m_poly->adj_F!=0)
				delete m_poly->adj_F;
			m_poly->adj_F = new AdjUpdate<size_undef>(m_P2F);
		}
	}
	else if (ADim==3)
	{
		if(m_T2R==0){
			m_T2R = new SmartVector<TabCellID<4> >(m_T2N->getMarks());
			m_tet->m_R = m_T2R;
			if(m_tet->adj_R!=0)
				delete m_tet->adj_R;
			m_tet->adj_R = new AdjUpdate<4>(m_T2R);
		}
		if(m_PY2R==0){
			m_PY2R = new SmartVector<TabCellID<5> >(m_PY2N->getMarks());
			m_pyra->m_R = m_PY2R;
			if(m_pyra->adj_R!=0)
				delete m_pyra->adj_R;
			m_pyra->adj_R = new AdjUpdate<5>(m_PY2R);
		}
		if(m_PR2R==0){
			m_PR2R = new SmartVector<TabCellID<5> >(m_PR2N->getMarks());
			m_prism3->m_R = m_PR2R;
			if(m_prism3->adj_R!=0)
				delete m_prism3->adj_R;
			m_prism3->adj_R = new AdjUpdate<5>(m_PR2R);
		}

		if(m_H2R==0){
			m_H2R = new SmartVector<TabCellID<6> >(m_H2N->getMarks());
			m_hex->m_R = m_H2R;
			if(m_hex->adj_R!=0)
				delete m_hex->adj_R;
			m_hex->adj_R = new AdjUpdate<6>(m_H2R);
		}
		if(m_P2R==0){
			m_P2R = new SmartVector<TabCellID<size_undef> >(m_P2N->getMarks());
			m_poly->m_R = m_P2R;
			if(m_poly->adj_R!=0)
				delete m_poly->adj_R;
			m_poly->adj_R = new AdjUpdate<size_undef>(m_P2R);
		}
	}
}
/*----------------------------------------------------------------------------*/
void RegionContainer::removeConnectivityContainers(const TInt ADim)
{
	if(ADim==0){
		if(m_T2N!=0){
			delete m_T2N;
			m_T2N=0;
                        if(m_tet->adj_N!=0) {
                                delete m_tet->adj_N;
				m_tet->adj_N = 0;
			}
		}
		if(m_PY2N!=0){
			delete m_PY2N;
			m_PY2N=0;
                        if(m_pyra->adj_N!=0) {
                                delete m_pyra->adj_N;
                                m_pyra->adj_N = 0;
                        }
		}
		if(m_PR2N!=0){
			delete m_PR2N;
			m_PR2N=0;
                        if(m_prism3->adj_N!=0) {
                                delete m_prism3->adj_N;
                                m_prism3->adj_N = 0;
                        }
		}
		if(m_H2N!=0){
			delete m_H2N;
			m_H2N=0;
                        if(m_hex->adj_N!=0) {
                                delete m_hex->adj_N;
                                m_hex->adj_N = 0;
                        }
		}
		if(m_P2N!=0){
			delete m_P2N;
			m_P2N=0;
                        if(m_poly->adj_N!=0) {
                                delete m_poly->adj_N;
                                m_poly->adj_N = 0;
                        }
		}
	}
	else if (ADim==1){
		if(m_T2E!=0){
			delete m_T2E;
			m_T2E=0;
                        if(m_tet->adj_E!=0) {
                                delete m_tet->adj_E;
                                m_tet->adj_E = 0;
                        }
		}
		if(m_PY2E!=0){
			delete m_PY2E;
			m_PY2E=0;
                        if(m_pyra->adj_E!=0) {
                                delete m_pyra->adj_E;
                                m_pyra->adj_E = 0;
                        }
		}

		if(m_PR2E!=0){
			delete m_PR2E;
			m_PR2E=0;
                        if(m_prism3->adj_E!=0) {
                                delete m_prism3->adj_E;
                                m_prism3->adj_E = 0;
                        }
		}

		if(m_H2E!=0){
			delete m_H2E;
			m_H2E=0;
                        if(m_hex->adj_E!=0) {
                                delete m_hex->adj_E;
                                m_hex->adj_E = 0;
                        }
		}
		if(m_P2E!=0){
			delete m_P2E;
			m_P2E=0;
                        if(m_poly->adj_E!=0) {
                                delete m_poly->adj_E;
                                m_poly->adj_E = 0;
                        }
		}
	}
	else if (ADim==2){
		if(m_T2F!=0){
			delete m_T2F;
			m_T2F=0;
                        if(m_tet->adj_F!=0) {
                                delete m_tet->adj_F;
                                m_tet->adj_F = 0;
                        }
		}
		if(m_PY2F!=0){
			delete m_PY2F;
			m_PY2F=0;
                        if(m_pyra->adj_F!=0) {
                                delete m_pyra->adj_F;
                                m_pyra->adj_F = 0;
                        }
		}
		if(m_PR2F!=0){
			delete m_PR2F;
			m_PR2F=0;
                        if(m_prism3->adj_F!=0) {
                                delete m_prism3->adj_F;
                                m_prism3->adj_F = 0;
                        }
		}

		if(m_H2F!=0){
			delete m_H2F;
			m_H2F=0;
                        if(m_hex->adj_F!=0) {
                                delete m_hex->adj_F;
                                m_hex->adj_F = 0;
                        }
		}
		if(m_P2F!=0){
			delete m_P2F;
			m_P2F=0;
                        if(m_poly->adj_F!=0) {
                                delete m_poly->adj_F;
                                m_poly->adj_F = 0;
                        }
		}
	}
	else if (ADim==3)
	{
		if(m_T2R!=0){
			delete m_T2R;
			m_T2R=0;
                        if(m_tet->adj_R!=0) {
                                delete m_tet->adj_R;
                                m_tet->adj_R = 0;
                        }
		}
		if(m_PY2R!=0){
			delete m_PY2R;
			m_PY2R=0;
                        if(m_pyra->adj_R!=0) {
                                delete m_pyra->adj_R;
                                m_pyra->adj_R = 0;
                        }
		}
		if(m_PR2R!=0){
			delete m_PR2R;
			m_PR2R=0;
                        if(m_prism3->adj_R!=0) {
                                delete m_prism3->adj_R;
                                m_prism3->adj_R = 0;
                        }
		}

		if(m_H2R!=0){
			delete m_H2R;
			m_H2R=0;
                        if(m_hex->adj_R!=0) {
                                delete m_hex->adj_R;
                                m_hex->adj_R = 0;
                        }
		}
		if(m_P2R!=0){
			delete m_P2R;
			m_P2R=0;
                        if(m_poly->adj_R!=0) {
                                delete m_poly->adj_R;
                                m_poly->adj_R = 0;
                        }
		}
	}
}
/*----------------------------------------------------------------------------*/
bool RegionContainer::has(const TCellID& AID) const
{
	return m_region_ids[AID];
}

/*----------------------------------------------------------------------------*/
Region RegionContainer::
addTet(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,const TCellID& AN4)
{
	//===========================================
	// STEP 1 - we value the region index
	//===========================================
	TInt index_region = m_region_ids.getFreeIndex();
	return addTet(AN1,AN2,AN3,AN4,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::
addTet(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,const TCellID& AN4,
		const TCellID& AID)
{
	TCellID index_region = AID;
	m_region_ids.assign(index_region);
	//===========================================
	// STEP 2 - we value the tet index
	//===========================================
	int index_type= m_tet->getID();
	//============================================
	// STEP 3 - we build and assign the reg. infos
	//============================================
	RegionInfo info(GMDS_TETRA,index_type);
	m_region_types.assign(info, index_region);
	//============================================
	// STEP 4 - we fill in R2N
	//============================================
	if(m_model.has(R2N)){
		//here only nodes are updated
		TabCellID<4> t;
		t.add(AN1);
		t.add(AN2);
		t.add(AN3);
		t.add(AN4);
		m_T2N->assign(t,index_type);
	}
	if(m_model.has(R2E)){
		TabCellID<6> t;
		m_T2E->assign(t,index_type);
	}
	if(m_model.has(R2F)){
		TabCellID<4> t;
		m_T2F->assign(t,index_type);
	}
	if(m_model.has(R2R)){
		TabCellID<4> t;
		m_T2R->assign(t,index_type);
	}

	//===============================================
	// STEP 5 - Temporary (or handle) face is created
	//===============================================
	return Region(m_mesh,GMDS_TETRA,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::
addPyramid(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
		const TCellID& AN4, const TCellID& AN5)
{
	//===========================================
	// STEP 1 - we value the region index
	//===========================================
	TInt index_region = m_region_ids.getFreeIndex();
	return addPyramid(AN1,AN2,AN3,AN4,AN5,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::
addPyramid(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
		const TCellID& AN4, const TCellID& AN5, const TCellID& AID)
{
	//===========================================
	// STEP 1 - we value the region index
	//===========================================
	TInt index_region = AID;

	m_region_ids.assign(index_region);
	//===========================================
	// STEP 2 - we value the pyramid index
	//===========================================
	int index_type= m_pyra->getID();
	//============================================
	// STEP 3 - we build and assign the reg. infos
	//============================================
	RegionInfo info(GMDS_PYRAMID,index_type);
	m_region_types.assign(info, index_region);
	//============================================
	// STEP 4 - we fill in R2N
	//============================================
	if(m_model.has(R2N)){
		//here only nodes are updated
		TabCellID<5> t;
		t.add(AN1);
		t.add(AN2);
		t.add(AN3);
		t.add(AN4);
		t.add(AN5);
		m_PY2N->assign(t,index_type);
	}
	if(m_model.has(R2E)){
		TabCellID<8> t;
		m_PY2E->assign(t,index_type);
	}
	if(m_model.has(R2F)){
		TabCellID<5> t;
		m_PY2F->assign(t,index_type);
	}
	if(m_model.has(R2R)){
		TabCellID<5> t;
		m_PY2R->assign(t,index_type);
	}

	//===============================================
	// STEP 5 - Temporary (or handle) face is created
	//===============================================
	return Region(m_mesh,GMDS_PYRAMID,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::
addPrism3(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
		const TCellID& AN4, const TCellID& AN5, const TCellID& AN6)
{
	//===========================================
	// STEP 1 - we value the region index
	//===========================================
	TInt index_region = m_region_ids.getFreeIndex();
	return addPrism3(AN1,AN2,AN3,AN4,AN5,AN6,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::
addPrism3(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
		const TCellID& AN4, const TCellID& AN5, const TCellID& AN6,
		const TCellID& AID)
{
	//===========================================
	// STEP 1 - we value the region index
	//===========================================
	TInt index_region = AID;

	m_region_ids.assign(index_region);
	//===========================================
	// STEP 2 - we value the pyramid index
	//===========================================
	int index_type= m_prism3->getID();
	//============================================
	// STEP 3 - we build and assign the reg. infos
	//============================================
	RegionInfo info(GMDS_PRISM3,index_type);
	m_region_types.assign(info, index_region);
	//============================================
	// STEP 4 - we fill in R2N
	//============================================
	if(m_model.has(R2N)){
		//here only nodes are updated
		TabCellID<6> t;
		t.add(AN1);
		t.add(AN2);
		t.add(AN3);
		t.add(AN4);
		t.add(AN5);
		t.add(AN6);
		m_PR2N->assign(t,index_type);
	}
	if(m_model.has(R2E)){
		TabCellID<9> t;
		m_PR2E->assign(t,index_type);
	}
	if(m_model.has(R2F)){
		TabCellID<5> t;
		m_PR2F->assign(t,index_type);
	}
	if(m_model.has(R2R)){
		TabCellID<5> t;
		m_PR2R->assign(t,index_type);
	}

	//===============================================
	// STEP 5 - Temporary (or handle) face is created
	//===============================================
	return Region(m_mesh,GMDS_PRISM3,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::
addHex(const TCellID& AN1, const TCellID& AN2, const TCellID& AN3,
			const TCellID& AN4, const TCellID& AN5, const TCellID& AN6,
			const TCellID& AN7, const TCellID& AN8)
{
	//===========================================
	// STEP 1 - we value the region index
	//===========================================
	TInt index_region = m_region_ids.getFreeIndex();
	return addHex(AN1,AN2,AN3,AN4,AN5,AN6,AN7,AN8,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::
addHex(const TCellID& AN1, const TCellID& AN2, const TCellID& AN3,
			const TCellID& AN4, const TCellID& AN5, const TCellID& AN6,
			const TCellID& AN7, const TCellID& AN8, const TCellID& AID)
{
	//===========================================
	// STEP 1 - we value the region index
	//===========================================
	TInt index_region = AID;

	m_region_ids.assign(index_region);
	//===========================================
	// STEP 2 - we value the hex index
	//===========================================
	int index_type= m_hex->getID();

	//============================================
	// STEP 3 - we build and assign the reg. infos
	//============================================
	RegionInfo info(GMDS_HEX,index_type);
	m_region_types.assign(info, index_region);

	//============================================
	// STEP 4 - we fill in R2N
	//============================================
	if(m_model.has(R2N)){
		//here only nodes are updated
		TabCellID<8> t;
		t.add(AN1);
		t.add(AN2);
		t.add(AN3);
		t.add(AN4);
		t.add(AN5);
		t.add(AN6);
		t.add(AN7);
		t.add(AN8);
		m_H2N->assign(t,index_type);
	}
	if(m_model.has(R2E)){
		TabCellID<12> t;
		m_H2E->assign(t,index_type);
	}
	if(m_model.has(R2F)){
		TabCellID<6> t;
		m_H2F->assign(t,index_type);
	}
	if(m_model.has(R2R)){
		TabCellID<6> t;
		m_H2R->assign(t,index_type);
	}
	//===============================================
	// STEP 5 - Temporary (or handle) face is created
	//===============================================
	return Region(m_mesh,GMDS_HEX,index_region);
}
/*----------------------------------------------------------------------------*/
Region RegionContainer::buildRegion(const TInt AIndex) const
{
	ECellType type = m_region_types[AIndex].type;
	return Region(m_mesh,type,AIndex);
}
/*----------------------------------------------------------------------------*/
void RegionContainer::
getNodesData(const TCellID& AID, int& ANbNodes) const
{
	RegionInfo info = m_region_types[AID];
	if(info.type==GMDS_TETRA)
	{
		if(m_tet->adj_N==0)
			ANbNodes=0;
		else
			ANbNodes = m_tet->adj_N->getSizeOfAnElement();
	}
	else if(info.type==GMDS_HEX)
	{
		if(m_hex->adj_N==0)
			ANbNodes=0;
		else
			ANbNodes = m_hex->adj_N->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PYRAMID)
	{
		if(m_pyra->adj_N==0)
			ANbNodes=0;
		else
			ANbNodes = m_pyra->adj_N->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PRISM3)
	{
		if(m_prism3->adj_N==0)
			ANbNodes=0;
		else
			ANbNodes = m_prism3->adj_N->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYHEDRA)
	{
		ANbNodes = (*m_poly->m_N)[info.type_id].size();
	}
}
/*----------------------------------------------------------------------------*/
void RegionContainer::
getEdgesData(const TCellID& AID, int& ANbEdges) const
{
	RegionInfo info = m_region_types[AID];
	if(info.type==GMDS_TETRA)
	{
		if(m_tet->adj_E==0)
			ANbEdges=0;
		else
			ANbEdges = m_tet->adj_E->getSizeOfAnElement();
	}
	else if(info.type==GMDS_HEX)
	{
		if(m_hex->adj_E==0)
			ANbEdges=0;
		else
			ANbEdges = m_hex->adj_E->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PYRAMID)
	{
		if(m_pyra->adj_E==0)
			ANbEdges=0;
		else
			ANbEdges = m_pyra->adj_E->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PRISM3)
	{
		if(m_prism3->adj_E==0)
			ANbEdges=0;
		else
			ANbEdges = m_prism3->adj_E->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYHEDRA)
	{
		ANbEdges = (*m_poly->m_E)[info.type_id].size();
	}
}
/*----------------------------------------------------------------------------*/
void RegionContainer::
getFacesData(const TCellID& AID, int& ANbFaces) const
{
	RegionInfo info = m_region_types[AID];
	if(info.type==GMDS_TETRA)
	{
		if(m_tet->adj_F==0)
			ANbFaces=0;
		else
			ANbFaces = m_tet->adj_F->getSizeOfAnElement();
	}
	else if(info.type==GMDS_HEX)
	{
		if(m_hex->adj_F==0)
			ANbFaces=0;
		else
			ANbFaces = m_hex->adj_F->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PYRAMID)
	{
		if(m_pyra->adj_F==0)
			ANbFaces=0;
		else
			ANbFaces = m_pyra->adj_F->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PRISM3)
	{
		if(m_prism3->adj_F==0)
			ANbFaces=0;
		else
			ANbFaces = m_prism3->adj_F->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYHEDRA)
	{
		ANbFaces = (*m_poly->m_F)[info.type_id].size();
	}
}
/*----------------------------------------------------------------------------*/
void RegionContainer::
getRegionsData(const TCellID& AID, int& ANbReg) const
{
	RegionInfo info = m_region_types[AID];
	if(info.type==GMDS_TETRA)
	{
		if(m_tet->adj_R==0)
			ANbReg=0;
		else
			ANbReg = m_tet->adj_R->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PYRAMID)
	{
		if(m_pyra->adj_R==0)
			ANbReg=0;
		else
			ANbReg = m_pyra->adj_R->getSizeOfAnElement();
	}
	else if(info.type==GMDS_PRISM3)
	{
		if(m_prism3->adj_R==0)
			ANbReg=0;
		else
			ANbReg = m_prism3->adj_R->getSizeOfAnElement();
	}
	else if(info.type==GMDS_HEX)
	{
		if(m_hex->adj_R==0)
			ANbReg=0;
		else
			ANbReg = m_hex->adj_R->getSizeOfAnElement();
	}
	else if(info.type==GMDS_POLYHEDRA)
	{
		ANbReg = (*m_poly->m_R)[info.type_id].size();
	}
}
/*----------------------------------------------------------------------------*/
void RegionContainer::remove(TInt AIndex)
{
	m_region_ids.unselect(AIndex);
	if(m_region_types[AIndex].type==GMDS_TETRA)
	{
		TCellID type_index = m_region_types[AIndex].type_id;
		if(m_model.has(R2N))
			m_T2N->remove(type_index);
		if(m_model.has(R2E))
			m_T2E->remove(type_index);
		if(m_model.has(R2F))
			m_T2F->remove(type_index);
		if(m_model.has(R2R))
			m_T2R->remove(type_index);

	}
	else if(m_region_types[AIndex].type==GMDS_PYRAMID)
	{
		TCellID type_index = m_region_types[AIndex].type_id;
		if(m_model.has(R2N))
			m_PY2N->remove(type_index);
		if(m_model.has(R2E))
			m_PY2E->remove(type_index);
		if(m_model.has(R2F))
			m_PY2F->remove(type_index);
		if(m_model.has(R2R))
			m_PY2R->remove(type_index);

	}
	else if(m_region_types[AIndex].type==GMDS_PRISM3)
	{
		TCellID type_index = m_region_types[AIndex].type_id;
		if(m_model.has(R2N))
			m_PR2N->remove(type_index);
		if(m_model.has(R2E))
			m_PR2E->remove(type_index);
		if(m_model.has(R2F))
			m_PR2F->remove(type_index);
		if(m_model.has(R2R))
			m_PR2R->remove(type_index);

	}
	else if(m_region_types[AIndex].type==GMDS_HEX)
	{
		TCellID type_index = m_region_types[AIndex].type_id;
		if(m_model.has(R2N))
			m_H2N->remove(type_index);
		if(m_model.has(R2E))
			m_H2E->remove(type_index);
		if(m_model.has(R2F))
			m_H2F->remove(type_index);
		if(m_model.has(R2R))
			m_H2R->remove(type_index);

	}
	else //POLYHEDRON
	{
		TCellID type_index = m_region_types[AIndex].type_id;
		if(m_model.has(R2N))
			m_P2N->remove(type_index);
		if(m_model.has(R2E))
			m_P2E->remove(type_index);
		if(m_model.has(R2F))
			m_P2F->remove(type_index);
		if(m_model.has(R2R))
			m_P2R->remove(type_index);
	}
}

/*----------------------------------------------------------------------------*/
void RegionContainer::clear()
{
	m_region_ids.clear();
	m_region_types.clear();

	if(m_T2N)
		m_T2N->clear();
	if(m_T2E)
		m_T2E->clear();
	if(m_T2F)
		m_T2F->clear();
	if(m_T2R)
		m_T2R->clear();

	if(m_H2N)
		m_H2N->clear();
	if(m_H2E)
		m_H2E->clear();
	if(m_H2F)
		m_H2F->clear();
	if(m_H2R)
		m_H2R->clear();

	if(m_PY2N)
		m_PY2N->clear();
	if(m_PY2E)
		m_PY2E->clear();
	if(m_PY2F)
		m_PY2F->clear();
	if(m_PY2R)
		m_PY2R->clear();

	if(m_PR2N)
		m_PR2N->clear();
	if(m_PR2E)
		m_PR2E->clear();
	if(m_PR2F)
		m_PR2F->clear();
	if(m_PR2R)
		m_PR2R->clear();

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
void RegionContainer::resize(const TInt ASize)
{
	m_region_ids.resize(ASize);
	m_region_types.resize(ASize);
}
/*----------------------------------------------------------------------------*/
void RegionContainer::update()
{
	m_region_ids.update();
}
/*----------------------------------------------------------------------------*/
void RegionContainer::serialize(std::ostream& AStr)
{
	m_region_ids.serialize(AStr);
	m_region_types.serialize(AStr);
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

	if(m_PY2N)
		m_PY2N->serialize(AStr);
	if(m_PY2E)
		m_PY2E->serialize(AStr);
	if(m_PY2F)
		m_PY2F->serialize(AStr);
	if(m_PY2R)
		m_PY2R->serialize(AStr);

	if(m_PR2N)
		m_PR2N->serialize(AStr);
	if(m_PR2E)
		m_PR2E->serialize(AStr);
	if(m_PR2F)
		m_PR2F->serialize(AStr);
	if(m_PR2R)
		m_PR2R->serialize(AStr);

	if(m_H2N)
		m_H2N->serialize(AStr);
	if(m_H2E)
		m_H2E->serialize(AStr);
	if(m_H2F)
		m_H2F->serialize(AStr);
	if(m_H2R)
		m_H2R->serialize(AStr);

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
void RegionContainer::unserialize(std::istream& AStr)
{
	m_region_ids.unserialize(AStr);
	m_region_types.unserialize(AStr);
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

	if(m_PY2N)
		m_PY2N->unserialize(AStr);
	if(m_PY2E)
		m_PY2E->unserialize(AStr);
	if(m_PY2F)
		m_PY2F->unserialize(AStr);
	if(m_PY2R)
		m_PY2R->unserialize(AStr);

	if(m_PR2N)
		m_PR2N->unserialize(AStr);
	if(m_PR2E)
		m_PR2E->unserialize(AStr);
	if(m_PR2F)
		m_PR2F->unserialize(AStr);
	if(m_PR2R)
		m_PR2R->unserialize(AStr);

	if(m_H2N)
		m_H2N->unserialize(AStr);
	if(m_H2E)
		m_H2E->unserialize(AStr);
	if(m_H2F)
		m_H2F->unserialize(AStr);
	if(m_H2R)
		m_H2R->unserialize(AStr);

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
RegionContainer::TAccessor::TAccessor(RegionContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(R2N))?AOwner->m_T2N:0;
	m_E =(AOwner->m_model.has(R2E))?AOwner->m_T2E:0;
	m_F =(AOwner->m_model.has(R2F))?AOwner->m_T2F:0;
	m_R =(AOwner->m_model.has(R2R))?AOwner->m_T2R:0;

	adj_N = (AOwner->m_model.has(R2N))? new AdjUpdate<4>(m_N):0;
	adj_E = (AOwner->m_model.has(R2E))? new AdjUpdate<6>(m_E):0;
	adj_F = (AOwner->m_model.has(R2F))? new AdjUpdate<4>(m_F):0;
	adj_R = (AOwner->m_model.has(R2R))? new AdjUpdate<4>(m_R):0;
}
/*----------------------------------------------------------------------------*/
RegionContainer::TAccessor::~TAccessor()
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
TInt RegionContainer::TAccessor::getID()
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
RegionContainer::PyAccessor::PyAccessor(RegionContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(R2N))?AOwner->m_PY2N:0;
	m_E =(AOwner->m_model.has(R2E))?AOwner->m_PY2E:0;
	m_F =(AOwner->m_model.has(R2F))?AOwner->m_PY2F:0;
	m_R =(AOwner->m_model.has(R2R))?AOwner->m_PY2R:0;

	adj_N = (AOwner->m_model.has(R2N))? new AdjUpdate<5>(m_N):0;
	adj_E = (AOwner->m_model.has(R2E))? new AdjUpdate<8>(m_E):0;
	adj_F = (AOwner->m_model.has(R2F))? new AdjUpdate<5>(m_F):0;
	adj_R = (AOwner->m_model.has(R2R))? new AdjUpdate<5>(m_R):0;
}
/*----------------------------------------------------------------------------*/
RegionContainer::PyAccessor::~PyAccessor()
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
TInt RegionContainer::PyAccessor::getID()
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
RegionContainer::PrAccessor::PrAccessor(RegionContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(R2N))?AOwner->m_PR2N:0;
	m_E =(AOwner->m_model.has(R2E))?AOwner->m_PR2E:0;
	m_F =(AOwner->m_model.has(R2F))?AOwner->m_PR2F:0;
	m_R =(AOwner->m_model.has(R2R))?AOwner->m_PR2R:0;

	adj_N = (AOwner->m_model.has(R2N))? new AdjUpdate<6>(m_N):0;
	adj_E = (AOwner->m_model.has(R2E))? new AdjUpdate<9>(m_E):0;
	adj_F = (AOwner->m_model.has(R2F))? new AdjUpdate<5>(m_F):0;
	adj_R = (AOwner->m_model.has(R2R))? new AdjUpdate<5>(m_R):0;
}
/*----------------------------------------------------------------------------*/
RegionContainer::PrAccessor::~PrAccessor()
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
TInt RegionContainer::PrAccessor::getID()
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
RegionContainer::HAccessor::HAccessor(RegionContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(R2N))?AOwner->m_H2N:0;
	m_E =(AOwner->m_model.has(R2E))?AOwner->m_H2E:0;
	m_F =(AOwner->m_model.has(R2F))?AOwner->m_H2F:0;
	m_R =(AOwner->m_model.has(R2R))?AOwner->m_H2R:0;

	adj_N = (AOwner->m_model.has(R2N))? new AdjUpdate<8>(m_N):0;
	adj_E = (AOwner->m_model.has(R2E))? new AdjUpdate<12>(m_E):0;
	adj_F = (AOwner->m_model.has(R2F))? new AdjUpdate<6>(m_F):0;
	adj_R = (AOwner->m_model.has(R2R))? new AdjUpdate<6>(m_R):0;
}
/*----------------------------------------------------------------------------*/
RegionContainer::HAccessor::~HAccessor()
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
TInt RegionContainer::HAccessor::getID()
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
RegionContainer::PAccessor::PAccessor(RegionContainer* AOwner, const MeshModel& AModel)
{
	m_owner = AOwner;
	m_N =(AOwner->m_model.has(R2N))?AOwner->m_P2N:0;
	m_E =(AOwner->m_model.has(R2E))?AOwner->m_P2E:0;
	m_F =(AOwner->m_model.has(R2F))?AOwner->m_P2F:0;
	m_R =(AOwner->m_model.has(R2R))?AOwner->m_P2R:0;

	adj_N = (AOwner->m_model.has(R2N))? new AdjUpdate<size_undef>(m_N):0;
	adj_E = (AOwner->m_model.has(R2E))? new AdjUpdate<size_undef>(m_E):0;
	adj_F = (AOwner->m_model.has(R2F))? new AdjUpdate<size_undef>(m_F):0;
	adj_R = (AOwner->m_model.has(R2R))? new AdjUpdate<size_undef>(m_R):0;
}
/*----------------------------------------------------------------------------*/
RegionContainer::PAccessor::~PAccessor()
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
TInt RegionContainer::PAccessor::getID()
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
