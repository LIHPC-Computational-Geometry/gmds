/*----------------------------------------------------------------------------*/
/*
 * Mesh.cpp
 *
 *  Created on: 5 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
//#include <gmds/ig/MeshDoctor.h>
//#include <gmds/CAD/GeomEntity.h>
#include <gmds/ig/CellGroup.h>

/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
template<> TInt
Mesh::newMark<Node>()
{
	if (m_nbUsedMarks_nodes == 31)
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_nodes == m_maxNbUsedMarks_nodes) m_maxNbUsedMarks_nodes = m_nbUsedMarks_nodes + 1;
#endif     // _DEBUG

	TInt mark = m_marks_nodes[m_nbUsedMarks_nodes++];
	m_usedMarks_nodes.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> TInt
Mesh::newMark<Edge>()
{
	if (m_nbUsedMarks_edges == 31)
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_edges == m_maxNbUsedMarks_edges) m_maxNbUsedMarks_edges = m_nbUsedMarks_edges + 1;
#endif     // _DEBUG

	TInt mark = m_marks_edges[m_nbUsedMarks_edges++];
	m_usedMarks_edges.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> TInt
Mesh::newMark<Face>()
{
	if (m_nbUsedMarks_faces == 31)
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_faces == m_maxNbUsedMarks_faces) m_maxNbUsedMarks_faces = m_nbUsedMarks_faces + 1;
#endif     // _DEBUG

	TInt mark = m_marks_faces[m_nbUsedMarks_faces++];
	m_usedMarks_faces.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> TInt
Mesh::newMark<Region>()
{
	if (m_nbUsedMarks_regions == 31)
		// not enough marks
		throw GMDSException("Limit of Boolean marks reached");
#ifdef _DEBUG

	if (m_nbUsedMarks_regions == m_maxNbUsedMarks_regions) m_maxNbUsedMarks_regions = m_nbUsedMarks_regions + 1;
#endif     // _DEBUG

	TInt mark = m_marks_regions[m_nbUsedMarks_regions++];
	m_usedMarks_regions.set(mark, true);
	return mark;
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::negateMaskMark<Node>(const TInt AMarkNumber)
{
	m_maskMarks_nodes.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::negateMaskMark<Edge>(const TInt AMarkNumber)
{
	m_maskMarks_edges.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::negateMaskMark<Face>(const TInt AMarkNumber)
{
	m_maskMarks_faces.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::negateMaskMark<Region>(const TInt AMarkNumber)
{
	m_maskMarks_regions.flip(AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool
Mesh::isMarked<Node>(const TCellID &ACellID, TInt AMarkNumber) const
{
	return (*m_marks[0])[ACellID][AMarkNumber] != m_maskMarks_nodes[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool
Mesh::isMarked(const Node &ACell, TInt AMarkNumber) const
{
	return isMarked<Node>(ACell.id(), AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool
Mesh::isMarked<Edge>(const TCellID &ACellID, TInt AMarkNumber) const
{
	return (*m_marks[1])[ACellID][AMarkNumber] != m_maskMarks_edges[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool
Mesh::isMarked(const Edge &ACell, TInt AMarkNumber) const
{
	return isMarked<Edge>(ACell.id(), AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool
Mesh::isMarked<Face>(const TCellID &ACellID, TInt AMarkNumber) const
{
	return (*m_marks[2])[ACellID][AMarkNumber] != m_maskMarks_faces[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool
Mesh::isMarked(const Face &ACell, TInt AMarkNumber) const
{
	return isMarked<Face>(ACell.id(), AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> bool
Mesh::isMarked<Region>(const TCellID &ACellID, TInt AMarkNumber) const
{
	return (*m_marks[3])[ACellID][AMarkNumber] != m_maskMarks_regions[AMarkNumber];
}
/*----------------------------------------------------------------------------*/
bool
Mesh::isMarked(const Region &ACell, TInt AMarkNumber) const
{
	return isMarked<Region>(ACell.id(), AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::markTo<Node>(const TCellID &ACellID, TInt AMarkNumber, bool AState)
{
	(*m_marks[0])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_nodes[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void
Mesh::markTo(const Node &ACell, TInt AMarkNumber, bool AState)
{
	markTo<Node>(ACell.id(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::markTo<Edge>(const TCellID &ACellID, TInt AMarkNumber, bool AState)
{
	(*m_marks[1])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_edges[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void
Mesh::markTo(const Edge &ACell, TInt AMarkNumber, bool AState)
{
	markTo<Edge>(ACell.id(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::markTo<Face>(const TCellID &ACellID, TInt AMarkNumber, bool AState)
{
	(*m_marks[2])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_faces[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void
Mesh::markTo(const Face &ACell, TInt AMarkNumber, bool AState)
{
	markTo<Face>(ACell.id(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::markTo<Region>(const TCellID &ACellID, TInt AMarkNumber, bool AState)
{
	(*m_marks[3])[ACellID].set(AMarkNumber, AState ^ m_maskMarks_regions[AMarkNumber]);
}
/*----------------------------------------------------------------------------*/
void
Mesh::markTo(const Region &ACell, TInt AMarkNumber, bool AState)
{
	markTo<Region>(ACell.id(), AMarkNumber, AState);
}
/*----------------------------------------------------------------------------*/
void
Mesh::mark(const Node &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::mark<Node>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Node>(ACellID, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
void
Mesh::mark(const Edge &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::mark<Edge>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Edge>(ACellID, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
void
Mesh::mark(const Face &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::mark<Face>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Face>(ACellID, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
void
Mesh::mark(const Region &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::mark<Region>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Region>(ACellID, AMarkNumber, true);
}
/*----------------------------------------------------------------------------*/
void
Mesh::unmark(const Node &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmark<Node>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Node>(ACellID, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
void
Mesh::unmark(const Edge &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmark<Edge>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Edge>(ACellID, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
void
Mesh::unmark(const Face &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmark<Face>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Face>(ACellID, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
void
Mesh::unmark(const Region &ACell, TInt AMarkNumber)
{
	markTo(ACell, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmark<Region>(const TCellID &ACellID, TInt AMarkNumber)
{
	markTo<Region>(ACellID, AMarkNumber, false);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmarkAll<Node>(const TInt AMarkNumber)
{
	for (auto id : nodes()) {
		unmark<Node>(id, AMarkNumber);
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmarkAll<Edge>(const TInt AMarkNumber)
{
	for (auto id : edges())
		unmark<Edge>(id, AMarkNumber);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmarkAll<Face>(const TInt AMarkNumber)
{
	for (auto id : faces()) {
		unmark<Face>(id, AMarkNumber);
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::unmarkAll<Region>(const TInt AMarkNumber)
{
	for (auto id : regions()) {
		unmark<Region>(id, AMarkNumber);
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::freeMark<Node>(TInt AMarkNumber)
{
	if (m_usedMarks_nodes.value(AMarkNumber)) {
		m_usedMarks_nodes.set(AMarkNumber, false);
		m_marks_nodes[--m_nbUsedMarks_nodes] = AMarkNumber;
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::freeMark<Edge>(TInt AMarkNumber)
{
	if (m_usedMarks_edges.value(AMarkNumber) == true) {
		m_usedMarks_edges.set(AMarkNumber, false);
		m_marks_edges[--m_nbUsedMarks_edges] = AMarkNumber;
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::freeMark<Face>(TInt AMarkNumber)
{
	if (m_usedMarks_faces.value(AMarkNumber) == true) {
		m_usedMarks_faces.set(AMarkNumber, false);
		m_marks_faces[--m_nbUsedMarks_faces] = AMarkNumber;
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::freeMark<Region>(TInt AMarkNumber)
{
	if (m_usedMarks_regions.value(AMarkNumber) == true) {
		m_usedMarks_regions.set(AMarkNumber, false);
		m_marks_regions[--m_nbUsedMarks_regions] = AMarkNumber;
	}
}
/*----------------------------------------------------------------------------*/
template<> Node
Mesh::get<Node>(const TCellID &AID) const
{
	return m_nodes_container->buildNode(AID);
}
/*----------------------------------------------------------------------------*/
template<> Edge
Mesh::get<Edge>(const TCellID &AID) const
{
	return m_edges_container->buildEdge(AID);
}
/*----------------------------------------------------------------------------*/
template<> Face
Mesh::get<Face>(const TCellID &AID) const
{
	return m_faces_container->buildFace(AID);
}
/*----------------------------------------------------------------------------*/
template<> Region
Mesh::get<Region>(const TCellID &AID) const
{
	return m_regions_container->buildRegion(AID);
}

/*----------------------------------------------------------------------------*/
template<> bool
Mesh::has<Node>(const TCellID &AID) const
{
	return m_nodes_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> bool
Mesh::has<Edge>(const TCellID &AID) const
{
	return m_edges_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> bool
Mesh::has<Face>(const TCellID &AID) const
{
	return m_faces_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> bool
Mesh::has<Region>(const TCellID &AID) const
{
	return m_regions_container->has(AID);
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::getAll<Node>(std::vector<Node> &AVec) const
{
	TInt nb_cells = getNbNodes();
	AVec.clear();
	AVec.resize(nb_cells);
	int i = 0;
	for (auto id : nodes()) {
		AVec[i++] = get<Node>(id);
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::getAll<Edge>(std::vector<Edge> &AVec) const
{
	TInt nb_cells = getNbEdges();
	AVec.clear();
	AVec.resize(nb_cells);
	int i = 0;
	for (auto id : edges()) {
		AVec[i++] = get<Edge>(id);
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::getAll<Face>(std::vector<Face> &AVec) const
{
	TInt nb_cells = getNbFaces();
	AVec.clear();
	AVec.resize(nb_cells);
	int i = 0;
	for (auto id : faces()) {
		AVec[i++] = get<Face>(id);
	}
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::getAll<Region>(std::vector<Region> &AVec) const
{
	TInt nb_cells = getNbRegions();
	AVec.clear();
	AVec.resize(nb_cells);
	int i = 0;
	for (auto id : regions()) {
		AVec[i++] = get<Region>(id);
	}
}
/*----------------------------------------------------------------------------*/
template<> Marks32
Mesh::getMarks<Node>(const Node &ACell) const
{
	return (*m_marks[0])[ACell.id()];
}
/*----------------------------------------------------------------------------*/
template<> Marks32
Mesh::getMarks<Edge>(const Edge &ACell) const
{
	return (*m_marks[1])[ACell.id()];
}
/*----------------------------------------------------------------------------*/
template<> Marks32
Mesh::getMarks<Face>(const Face &ACell) const
{
	return (*m_marks[2])[ACell.id()];
}
/*----------------------------------------------------------------------------*/
template<> Marks32
Mesh::getMarks<Region>(const Region &ACell) const
{
	return (*m_marks[3])[ACell.id()];
}
/*----------------------------------------------------------------------------*/
Mesh::Mesh(MeshModel model) : m_model(model)
{
	m_nodes_container = new NodeContainer(this);
	m_edges_container = new EdgeContainer(this);
	m_faces_container = new FaceContainer(this);
	m_regions_container = new RegionContainer(this);

	m_node_variable_manager = new VariableManager();
	m_edge_variable_manager = new VariableManager();
	m_face_variable_manager = new VariableManager();
	m_region_variable_manager = new VariableManager();
	/* init all the bits to false*/
	m_maskMarks_nodes.reset();
	m_maskMarks_edges.reset();
	m_maskMarks_faces.reset();
	m_maskMarks_regions.reset();
	m_usedMarks_nodes.reset();
	m_usedMarks_edges.reset();
	m_usedMarks_faces.reset();
	m_usedMarks_regions.reset();

	m_nbUsedMarks_nodes = 0;
	m_nbUsedMarks_edges = 0;
	m_nbUsedMarks_faces = 0;
	m_nbUsedMarks_regions = 0;
#ifdef _DEBUG
	m_maxNbUsedMarks_nodes = 0;
	m_maxNbUsedMarks_edges = 0;
	m_maxNbUsedMarks_faces = 0;
	m_maxNbUsedMarks_regions = 0;
#endif     // _DEBUG

	m_marks[0]= newVariable<Marks32, GMDS_NODE>("mark");
	m_marks[1]= newVariable<Marks32, GMDS_EDGE>("mark");
	m_marks[2]= newVariable<Marks32, GMDS_FACE>("mark");
	m_marks[3]= newVariable<Marks32, GMDS_REGION>("mark");
	for (auto i = 0; i < 32; ++i)
	{
		m_marks_nodes[i] = i;
		m_marks_edges[i] = i;
		m_marks_faces[i] = i;
		m_marks_regions[i] = i;
	}
}
/*----------------------------------------------------------------------------*/
Mesh::~Mesh()
{

	delete m_nodes_container;
	delete m_edges_container;
	delete m_faces_container;
	delete m_regions_container;
	delete m_node_variable_manager;
	delete m_edge_variable_manager;
	delete m_face_variable_manager;
	delete m_region_variable_manager;

	for (auto g : m_clouds) {
		delete g;
	}
	for (auto g : m_lines) {
		delete g;
	}
	for (auto g : m_surfaces) {
		delete g;
	}
	for (auto g : m_volumes) {
		delete g;
	}
}
/*----------------------------------------------------------------------------*/
MeshModel
Mesh::getModel() const
{
	return m_model;
}
/*----------------------------------------------------------------------------*/
void
Mesh::clear()
{
	if (m_nodes_container) m_nodes_container->clear();
	if (m_edges_container) m_edges_container->clear();
	if (m_faces_container) m_faces_container->clear();
	if (m_regions_container) m_regions_container->clear();

	m_node_variable_manager->clearVariables();
	m_edge_variable_manager->clearVariables();
	m_face_variable_manager->clearVariables();
	m_region_variable_manager->clearVariables();

	m_clouds.clear();
	m_lines.clear();
	m_surfaces.clear();
	m_volumes.clear();
}
/*----------------------------------------------------------------------------*/
void
Mesh::changeModel(const MeshModel &AModel, const bool &ACallDoctor)
{
	if (ACallDoctor) {
		changeModelWithDoctor(AModel);
	}
	else {
		changeModelWithoutDoctor(AModel);
	}
}
/*----------------------------------------------------------------------------*/
void
Mesh::changeModelWithoutDoctor(const MeshModel &AModel)
{
	MeshModel model_old = m_model;
	MeshModel model_new = AModel;

	/* the reference mesh model must be changed now since base value classes
	 * (as the Node or Face ones for instance), will check some model
	 * properties before applying an operation.
	 */
	m_model = AModel;
	m_nodes_container->setModel(m_model);
	m_edges_container->setModel(m_model);
	m_faces_container->setModel(m_model);
	m_regions_container->setModel(m_model);
	// We look for the cells and relationships that differ between the 2 models
	MeshModel model_oldBUTnew = MeshModel::exclusion(model_old, model_new);
	MeshModel model_newBUTold = MeshModel::exclusion(model_new, model_old);

	/* Containers for cells of any dimension remain available all the time
	 * (see m_X_container), but adjacency and incidence relations X2Y must be
	 * added or removed
	 */

	/* We add the necesseray containers
	 */
	if (model_newBUTold.has(N2N)) m_nodes_container->addConnectivityContainers(0);
	if (model_newBUTold.has(N2E)) m_nodes_container->addConnectivityContainers(1);
	if (model_newBUTold.has(N2F)) m_nodes_container->addConnectivityContainers(2);
	if (model_newBUTold.has(N2R)) m_nodes_container->addConnectivityContainers(3);

	if (model_newBUTold.has(E2N)) m_edges_container->addConnectivityContainers(0);
	if (model_newBUTold.has(E2E)) m_edges_container->addConnectivityContainers(1);
	if (model_newBUTold.has(E2F)) m_edges_container->addConnectivityContainers(2);
	if (model_newBUTold.has(E2R)) m_edges_container->addConnectivityContainers(3);

	if (model_newBUTold.has(F2N)) m_faces_container->addConnectivityContainers(0);
	if (model_newBUTold.has(F2E)) m_faces_container->addConnectivityContainers(1);
	if (model_newBUTold.has(F2F)) m_faces_container->addConnectivityContainers(2);
	if (model_newBUTold.has(F2R)) m_faces_container->addConnectivityContainers(3);

	if (model_newBUTold.has(R2N)) m_regions_container->addConnectivityContainers(0);
	if (model_newBUTold.has(R2E)) m_regions_container->addConnectivityContainers(1);
	if (model_newBUTold.has(R2F)) m_regions_container->addConnectivityContainers(2);
	if (model_newBUTold.has(R2R)) m_regions_container->addConnectivityContainers(3);

	// REMOVING
	if (model_oldBUTnew.has(N)) {
		m_nodes_container->clear();
	}
	if (model_oldBUTnew.has(N2N)) m_nodes_container->removeConnectivityContainers(0);
	if (model_oldBUTnew.has(N2E)) m_nodes_container->removeConnectivityContainers(1);
	if (model_oldBUTnew.has(N2F)) m_nodes_container->removeConnectivityContainers(2);
	if (model_oldBUTnew.has(N2R)) m_nodes_container->removeConnectivityContainers(3);

	if (model_oldBUTnew.has(E)) {
		m_edges_container->clear();
	}
	if (model_oldBUTnew.has(E2N)) m_edges_container->removeConnectivityContainers(0);
	if (model_oldBUTnew.has(E2E)) m_edges_container->removeConnectivityContainers(1);
	if (model_oldBUTnew.has(E2F)) m_edges_container->removeConnectivityContainers(2);
	if (model_oldBUTnew.has(E2R)) m_edges_container->removeConnectivityContainers(3);

	if (model_oldBUTnew.has(F)) {
		m_faces_container->clear();
	}
	if (model_oldBUTnew.has(F2N)) m_faces_container->removeConnectivityContainers(0);
	if (model_oldBUTnew.has(F2E)) m_faces_container->removeConnectivityContainers(1);
	if (model_oldBUTnew.has(F2F)) m_faces_container->removeConnectivityContainers(2);
	if (model_oldBUTnew.has(F2R)) m_faces_container->removeConnectivityContainers(3);

	if (model_oldBUTnew.has(R)) {
		m_regions_container->clear();
	}
	if (model_oldBUTnew.has(R2N)) m_regions_container->removeConnectivityContainers(0);
	if (model_oldBUTnew.has(R2E)) m_regions_container->removeConnectivityContainers(1);
	if (model_oldBUTnew.has(R2F)) m_regions_container->removeConnectivityContainers(2);
	if (model_oldBUTnew.has(R2R)) m_regions_container->removeConnectivityContainers(3);
}
/*----------------------------------------------------------------------------*/
void
Mesh::changeModelWithDoctor(const MeshModel &AModel)
{
	//	MeshModel model_old = m_model;
	//	MeshModel model_new = AModel;
	//
	//	/* the reference mesh model must be changed now since base value classes
	//	 * (as the Node or Face ones for instance), will check some model
	//	 * properties before applying an operation.
	//	 */
	//	m_model = AModel;
	//	m_nodes_container->setModel(m_model);
	//	m_edges_container->setModel(m_model);
	//	m_faces_container->setModel(m_model);
	//	m_regions_container->setModel(m_model);
	//	// We look for the cells and relationships that differ between the 2 models
	//	MeshModel model_oldBUTnew = MeshModel::exclusion	(model_old, model_new);
	//	MeshModel model_newBUTold = MeshModel::exclusion	(model_new, model_old);
	//
	//	/* Containers for cells of any dimension remain available all the time
	//	 * (see m_X_container), but adjacency and incidence relations X2Y must be
	//	 * added or removed
	//	*/
	//	MeshDoctor doc(this);
	//
	//	/* We add the necesseray containers
	//	 */
	//	if(model_newBUTold.has(N2N))
	//		m_nodes_container->addConnectivityContainers(0);
	//	if(model_newBUTold.has(N2E))
	//		m_nodes_container->addConnectivityContainers(1);
	//	if(model_newBUTold.has(N2F))
	//		m_nodes_container->addConnectivityContainers(2);
	//	if(model_newBUTold.has(N2R))
	//		m_nodes_container->addConnectivityContainers(3);
	//
	//	if(model_newBUTold.has(E2N))
	//		m_edges_container->addConnectivityContainers(0);
	//	if(model_newBUTold.has(E2E))
	//		m_edges_container->addConnectivityContainers(1);
	//	if(model_newBUTold.has(E2F))
	//		m_edges_container->addConnectivityContainers(2);
	//	if(model_newBUTold.has(E2R))
	//		m_edges_container->addConnectivityContainers(3);
	//
	//	if(model_newBUTold.has(F2N))
	//		m_faces_container->addConnectivityContainers(0);
	//	if(model_newBUTold.has(F2E))
	//		m_faces_container->addConnectivityContainers(1);
	//	if(model_newBUTold.has(F2F))
	//		m_faces_container->addConnectivityContainers(2);
	//	if(model_newBUTold.has(F2R))
	//		m_faces_container->addConnectivityContainers(3);
	//
	//	if(model_newBUTold.has(R2N))
	//		m_regions_container->addConnectivityContainers(0);
	//	if(model_newBUTold.has(R2E))
	//		m_regions_container->addConnectivityContainers(1);
	//	if(model_newBUTold.has(R2F))
	//		m_regions_container->addConnectivityContainers(2);
	//	if(model_newBUTold.has(R2R))
	//		m_regions_container->addConnectivityContainers(3);
	//
	//	// Now missing cells are created under the assumption that the X2N
	//	// adjacencies exists for X= R, F or E
	//	if(model_newBUTold.has(F))
	//		doc.buildF();
	//	if(model_newBUTold.has(E))
	//		doc.buildE();
	//
	//	//Data in newBUTold must be added, then those in old_butNew will be removed
	//	//The addition is done beginning with downard ajacencies then upward.
	//
	//	model_old = model_new;
	//	if(model_newBUTold.has(F2E))
	//		doc.buildF2E(model_old);
	//
	//	if(model_newBUTold.has(R2N))
	//		doc.buildR2N(model_old);
	//	if(model_newBUTold.has(R2E))
	//		doc.buildR2E(model_old);
	//	if(model_newBUTold.has(R2F))
	//		doc.buildR2F(model_old);
	//
	//	if(model_newBUTold.has(N2N))
	//		doc.buildN2N(model_old);
	//	if(model_newBUTold.has(N2E))
	//		doc.buildN2E(model_old);
	//	if(model_newBUTold.has(N2F))
	//		doc.buildN2F(model_old);
	//	if(model_newBUTold.has(N2R))
	//		doc.buildN2R(model_old);
	//
	//	if(model_newBUTold.has(E2E))
	//		doc.buildE2E(model_old);
	//	if(model_newBUTold.has(E2F))
	//		doc.buildE2F(model_old);
	//	if(model_newBUTold.has(E2R))
	//		doc.buildE2R(model_old);
	//
	//	if(model_newBUTold.has(F2F))
	//		doc.buildF2F(model_old);
	//	if(model_newBUTold.has(F2R))
	//		doc.buildF2R(model_old);
	//
	//	if(model_newBUTold.has(R2R))
	//		doc.buildR2R(model_old);
	//
	//
	//
	//
	//	//REMOVING
	//	if(model_oldBUTnew.has(N)) {
	//		m_nodes_container->clear();
	//	}
	//	if(model_oldBUTnew.has(N2N))
	//		m_nodes_container->removeConnectivityContainers(0);
	//	if(model_oldBUTnew.has(N2E))
	//		m_nodes_container->removeConnectivityContainers(1);
	//	if(model_oldBUTnew.has(N2F))
	//		m_nodes_container->removeConnectivityContainers(2);
	//	if(model_oldBUTnew.has(N2R))
	//		m_nodes_container->removeConnectivityContainers(3);
	//
	//	if(model_oldBUTnew.has(E)) {
	//                m_edges_container->clear();
	//        }
	//	if(model_oldBUTnew.has(E2N))
	//		m_edges_container->removeConnectivityContainers(0);
	//	if(model_oldBUTnew.has(E2E))
	//		m_edges_container->removeConnectivityContainers(1);
	//	if(model_oldBUTnew.has(E2F))
	//		m_edges_container->removeConnectivityContainers(2);
	//	if(model_oldBUTnew.has(E2R))
	//		m_edges_container->removeConnectivityContainers(3);
	//
	//	if(model_oldBUTnew.has(F)) {
	//                m_faces_container->clear();
	//        }
	//	if(model_oldBUTnew.has(F2N))
	//		m_faces_container->removeConnectivityContainers(0);
	//	if(model_oldBUTnew.has(F2E))
	//		m_faces_container->removeConnectivityContainers(1);
	//	if(model_oldBUTnew.has(F2F))
	//		m_faces_container->removeConnectivityContainers(2);
	//	if(model_oldBUTnew.has(F2R))
	//		m_faces_container->removeConnectivityContainers(3);
	//
	//	if(model_oldBUTnew.has(R)) {
	//                m_regions_container->clear();
	//        }
	//	if(model_oldBUTnew.has(R2N))
	//		m_regions_container->removeConnectivityContainers(0);
	//	if(model_oldBUTnew.has(R2E))
	//		m_regions_container->removeConnectivityContainers(1);
	//	if(model_oldBUTnew.has(R2F))
	//		m_regions_container->removeConnectivityContainers(2);
	//	if(model_oldBUTnew.has(R2R))
	//		m_regions_container->removeConnectivityContainers(3);
}
/*----------------------------------------------------------------------------*/
Node
Mesh::newNode(const TCoord &AX, const TCoord &AY, const TCoord AZ)
{
	Node n = m_nodes_container->add(AX, AY, AZ);
	m_node_variable_manager->addEntry(n.id());
	return n;
}
/*----------------------------------------------------------------------------*/
Node
Mesh::newNode(const math::Point &APnt)
{
	return newNode(APnt.X(), APnt.Y(), APnt.Z());
}
/*----------------------------------------------------------------------------*/
Edge
Mesh::newEdge(const Node &AN1, const Node &AN2)
{
	return newEdge(AN1.id(), AN2.id());
}
/*----------------------------------------------------------------------------*/
Edge
Mesh::newEdge(const TCellID &AN1, const TCellID &AN2)
{
	Edge e = m_edges_container->add(AN1, AN2);
	m_edge_variable_manager->addEntry(e.id());
	return e;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newTriangle(const Node &AN1, const Node &AN2, const Node &AN3)
{
	return newTriangle(AN1.id(), AN2.id(), AN3.id());
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newTriangle(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3)
{
	Face f = m_faces_container->addTriangle(AN1, AN2, AN3);
	m_face_variable_manager->addEntry(f.id());
	return f;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newQuad(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4)
{
	return newQuad(AN1.id(), AN2.id(), AN3.id(), AN4.id());
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newQuad(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4)
{
	Face f = m_faces_container->addQuad(AN1, AN2, AN3, AN4);
	m_face_variable_manager->addEntry(f.id());
	return f;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newPolygon(const std::vector<Node> &ANodes)
{
	std::vector<TCellID> ids;
	ids.resize(ANodes.size());
	for (unsigned int i = 0; i < ANodes.size(); i++)
		ids[i] = ANodes[i].id();

	return newPolygon(ids);
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newPolygon(const std::vector<TCellID> &ANodes)
{
	Face f = m_faces_container->addPolygon(ANodes);
	m_face_variable_manager->addEntry(f.id());
	return f;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newFace(const std::vector<Node> &ANodes)
{
	if (ANodes.size() == 3) {
		return newTriangle(ANodes[0].id(), ANodes[1].id(), ANodes[2].id());
	}
	else if (ANodes.size() == 4) {
		return newQuad(ANodes[0].id(), ANodes[1].id(), ANodes[2].id(), ANodes[3].id());
	}
	else {
		std::vector<TCellID> ids(ANodes.size());
		for (unsigned int iNode = 0; iNode < ANodes.size(); iNode++) {
			ids[iNode] = ANodes[iNode].id();
		}
		return newPolygon(ids);
	}
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newFace(const std::vector<TCellID> &ANodes)
{
	Face f;
	if (ANodes.size() == 3)
		f = newTriangle(ANodes[0], ANodes[1], ANodes[2]);
	else if (ANodes.size() == 4)
		f = newQuad(ANodes[0], ANodes[1], ANodes[2], ANodes[3]);
	else
		f = newPolygon(ANodes);

	return f;
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newTet(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4)
{
	return newTet(AN1.id(), AN2.id(), AN3.id(), AN4.id());
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newTet(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4)
{
	Region r = m_regions_container->addTet(AN1, AN2, AN3, AN4);
	m_region_variable_manager->addEntry(r.id());
	return r;
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newPyramid(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4, const Node &AN5)
{
	return newPyramid(AN1.id(), AN2.id(), AN3.id(), AN4.id(), AN5.id());
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newPyramid(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AN5)
{
	Region r = m_regions_container->addPyramid(AN1, AN2, AN3, AN4, AN5);
	m_region_variable_manager->addEntry(r.id());
	return r;
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newPrism3(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4, const Node &AN5, const Node &AN6)
{
	return newPrism3(AN1.id(), AN2.id(), AN3.id(), AN4.id(), AN5.id(), AN6.id());
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newPrism3(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AN5, const TCellID &AN6)
{
	Region r = m_regions_container->addPrism3(AN1, AN2, AN3, AN4, AN5, AN6);
	m_region_variable_manager->addEntry(r.id());
	return r;
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newHex(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4, const Node &AN5, const Node &AN6, const Node &AN7, const Node &AN8)
{
	return newHex(AN1.id(), AN2.id(), AN3.id(), AN4.id(), AN5.id(), AN6.id(), AN7.id(), AN8.id());
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newHex(const TCellID &AN1,
             const TCellID &AN2,
             const TCellID &AN3,
             const TCellID &AN4,
             const TCellID &AN5,
             const TCellID &AN6,
             const TCellID &AN7,
             const TCellID &AN8)
{
	Region r = m_regions_container->addHex(AN1, AN2, AN3, AN4, AN5, AN6, AN7, AN8);
	m_region_variable_manager->addEntry(r.id());
	return r;
}
/*----------------------------------------------------------------------------*/
TCellID
Mesh::getMaxLocalID(const TInt &ADim) const
{
	TInt max_id = 0;
	if (ADim == 0)
		max_id = m_nodes_container->getMaxID();
	else if (ADim == 1)
		max_id = m_edges_container->getMaxID();
	else if (ADim == 2)
		max_id = m_faces_container->getMaxID();
	else if (ADim == 3)
		max_id = m_regions_container->getMaxID();

	return max_id;
}
/*----------------------------------------------------------------------------*/
void
Mesh::deleteVariable(ECellType AType, const std::string &AName)
{
	switch (AType) {
	case GMDS_NODE: m_node_variable_manager->deleteVariable(AName); break;
	case GMDS_EDGE: m_edge_variable_manager->deleteVariable(AName); break;
	case GMDS_FACE: m_face_variable_manager->deleteVariable(AName); break;
	case GMDS_REGION: m_region_variable_manager->deleteVariable(AName); break;
	default: throw GMDSException("Unmanaged type of value -> impossible to delete a variable");
	}
}
/*----------------------------------------------------------------------------*/
void
Mesh::deleteVariable(ECellType AType, VariableItf *AVar)
{
	switch (AType) {
	case GMDS_NODE: m_node_variable_manager->deleteVariable(AVar); break;
	case GMDS_EDGE: m_edge_variable_manager->deleteVariable(AVar); break;
	case GMDS_FACE: m_face_variable_manager->deleteVariable(AVar); break;
	case GMDS_REGION: m_region_variable_manager->deleteVariable(AVar); break;
	default: throw GMDSException("Unmanaged type of value -> impossible to delete a variable");
	}
}

/*----------------------------------------------------------------------------*/
std::vector<VariableItf *>
Mesh::getAllVariables(ECellType AType) const
{
	switch (AType) {
	case GMDS_NODE: {
		return m_node_variable_manager->getAllVariables();
	} break;
	case GMDS_EDGE: {
		return m_edge_variable_manager->getAllVariables();
	} break;
	case GMDS_FACE: {
		return m_face_variable_manager->getAllVariables();
	} break;
	case GMDS_REGION: {
		return m_region_variable_manager->getAllVariables();
	} break;
	default: throw GMDSException("Unmanaged type of value -> impossible to access to a variable");
	}
}
/*----------------------------------------------------------------------------*/
bool
Mesh::hasVariable(ECellType AType, const std::string &AName) const
{
	switch (AType) {
	case GMDS_NODE: {
		return m_node_variable_manager->doesVariableExist(AName);
	} break;
	case GMDS_EDGE: {
		return m_edge_variable_manager->doesVariableExist(AName);
	} break;
	case GMDS_FACE: {
		return m_face_variable_manager->doesVariableExist(AName);
	} break;
	case GMDS_REGION: {
		return m_region_variable_manager->doesVariableExist(AName);
	} break;
	default: throw GMDSException("Unmanaged type of value -> impossible to access to a variable");
	}

	return false;
}
/*----------------------------------------------------------------------------*/
void
Mesh::clearAndResizeNodeIDContainer(const TInt AMaxID)
{
	m_nodes_container->clear();
	m_nodes_container->resize(AMaxID);
	m_node_variable_manager->clearVariables();
}
/*----------------------------------------------------------------------------*/
void
Mesh::clearAndResizeEdgeIDContainer(const TInt AMaxID)
{
	m_edges_container->clear();
	m_edges_container->resize(AMaxID);
	m_edge_variable_manager->clearVariables();
}
/*----------------------------------------------------------------------------*/
void
Mesh::clearAndResizeFaceIDContainer(const TInt AMaxID)
{
	m_faces_container->clear();
	m_faces_container->resize(AMaxID);
	m_face_variable_manager->clearVariables();
}
/*----------------------------------------------------------------------------*/
void
Mesh::clearAndResizeRegionIDContainer(const TInt AMaxID)
{
	m_regions_container->clear();
	m_regions_container->resize(AMaxID);
	m_region_variable_manager->clearVariables();
}
/*----------------------------------------------------------------------------*/
void
Mesh::updateIDContainers()
{
	m_nodes_container->update();
	m_edges_container->update();
	m_faces_container->update();
	m_regions_container->update();
}
/*----------------------------------------------------------------------------*/
Node
Mesh::newNodeWithID(const TCoord &AX, const TCoord &AY, const TCoord &AZ, const TCellID &AGID)
{
	Node n = m_nodes_container->add(AX, AY, AZ, AGID);
	m_node_variable_manager->addEntry(AGID);
	return n;
}
/*----------------------------------------------------------------------------*/
Edge
Mesh::newEdgeWithID(const TCellID &AV1, const TCellID &AV2, const TCellID &AGID)
{
	Edge e = m_edges_container->add(AV1, AV2, AGID);
	m_edge_variable_manager->addEntry(AGID);
	return e;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newTriangleWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AGID)
{
	Face f = m_faces_container->addTriangle(AN1, AN2, AN3, AGID);
	m_face_variable_manager->addEntry(AGID);
	return f;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newQuadWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AGID)
{
	Face f = m_faces_container->addQuad(AN1, AN2, AN3, AN4, AGID);
	m_face_variable_manager->addEntry(AGID);
	return f;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newPolygonWithID(const std::vector<TCellID> &ANodes, const TCellID &AGID)
{
	Face f = m_faces_container->addPolygon(ANodes, AGID);
	m_face_variable_manager->addEntry(AGID);
	return f;
}
/*----------------------------------------------------------------------------*/
Face
Mesh::newFaceWithID(const std::vector<TCellID> &ANodes, const TCellID &AGID)
{
	Face f;
	if (ANodes.size() == 3)
		f = newTriangleWithID(ANodes[0], ANodes[1], ANodes[2], AGID);
	else if (ANodes.size() == 4)
		f = newQuadWithID(ANodes[0], ANodes[1], ANodes[2], ANodes[3], AGID);
	else
		f = newPolygonWithID(ANodes, AGID);

	return f;
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newTetWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AGID)
{
	Region r = m_regions_container->addTet(AN1, AN2, AN3, AN4, AGID);
	m_region_variable_manager->addEntry(AGID);
	return r;
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newPyramidWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AN5, const TCellID &AGID)
{
	Region r = m_regions_container->addPyramid(AN1, AN2, AN3, AN4, AN5, AGID);
	m_region_variable_manager->addEntry(AGID);
	return r;
}
/*----------------------------------------------------------------------------*/
Region
Mesh::newHexWithID(const TCellID &AN1,
                   const TCellID &AN2,
                   const TCellID &AN3,
                   const TCellID &AN4,
                   const TCellID &AN5,
                   const TCellID &AN6,
                   const TCellID &AN7,
                   const TCellID &AN8,
                   const TCellID &AGID)
{
	Region r = m_regions_container->addHex(AN1, AN2, AN3, AN4, AN5, AN6, AN7, AN8, AGID);
	m_region_variable_manager->addEntry(AGID);
	return r;
}
/*----------------------------------------------------------------------------*/
void
Mesh::initializeGeometryClassification()
{
	classification[0] = newVariable<cad::GeomEntity *, GMDS_NODE>("NodeClassification");
	classification[1] = newVariable<cad::GeomEntity *, GMDS_EDGE>("EdgeClassification");
	classification[2] = newVariable<cad::GeomEntity *, GMDS_FACE>("FaceClassification");
	classification[3] = newVariable<cad::GeomEntity *, GMDS_REGION>("RegionClassification");

	for (auto &i : classification) {
		i->setValuesTo(0);
	}
}
/*----------------------------------------------------------------------------*/
bool
Mesh::doesGeometricClassificationExist(const int ADim)
{

	switch (ADim) {
	case 0: {
		return m_node_variable_manager->doesVariableExist("NodeClassification");
	} break;
	case 1: {
		return m_edge_variable_manager->doesVariableExist("EdgeClassification");
	} break;
	case 2: {
		return m_face_variable_manager->doesVariableExist("FaceClassification");
	} break;
	case 3: {
		return m_region_variable_manager->doesVariableExist("RegionClassification");
	} break;
	default: throw GMDSException("Mesh::doesGeometricClassificationExist : bad ADim");
	}
}
/*----------------------------------------------------------------------------*/
Variable<cad::GeomEntity *> *
Mesh::getGeometricClassification(const int ADim)
{
	return classification[ADim];
}
/*----------------------------------------------------------------------------*/
cad::GeomEntity *
Mesh::getGeometricClassification(const Node &ACell)
{
	return (*(classification[0]))[ACell.id()];
}
/*----------------------------------------------------------------------------*/
cad::GeomEntity *
Mesh::getGeometricClassification(const Edge &ACell)
{
	return (*(classification[1]))[ACell.id()];
}
/*----------------------------------------------------------------------------*/
cad::GeomEntity *
Mesh::getGeometricClassification(const Face &ACell)
{
	return (*(classification[2]))[ACell.id()];
}
/*----------------------------------------------------------------------------*/
cad::GeomEntity *
Mesh::getGeometricClassification(const Region &ACell)
{
	return (*(classification[3]))[ACell.id()];
}
/*----------------------------------------------------------------------------*/
void
Mesh::setGeometricClassification(const Node &ACell, cad::GeomEntity *AGE)
{
	(classification[0])->set(ACell.id(), AGE);
}
/*----------------------------------------------------------------------------*/
void
Mesh::setGeometricClassification(const Edge &ACell, cad::GeomEntity *AGE)
{
	(classification[1])->set(ACell.id(), AGE);
}
/*----------------------------------------------------------------------------*/
void
Mesh::setGeometricClassification(const Face &ACell, cad::GeomEntity *AGE)
{
	(classification[2])->set(ACell.id(), AGE);
}
/*----------------------------------------------------------------------------*/
void
Mesh::setGeometricClassification(const Region &ACell, cad::GeomEntity *AGE)
{
	(classification[3])->set(ACell.id(), AGE);
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Node> *
Mesh::newGroup<Node>(const std::string &AName)
{
	auto it = m_clouds.begin();
	for (; it != m_clouds.end(); it++)
		if ((*it)->name() == AName) {
			const std::string mess = "A node group named " + AName + " already exists !";
			throw GMDSException(mess);
		}

	auto *c = new CellGroup<Node>(this, AName);
	m_clouds.push_back(c);
	return m_clouds.back();
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::deleteGroup<Node>(CellGroup<Node> *ACloud)
{
	auto it = m_clouds.begin();
	for (; it != m_clouds.end(); it++)
		if (ACloud == *it) {
			m_clouds.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Node> *
Mesh::getGroup<Node>(const unsigned int AIndex)
{
	if (AIndex >= m_clouds.size()) {
		std::stringstream mess;
		mess << "There is no node group indexed " << AIndex;
		throw GMDSException(mess.str());
	}
	unsigned int index = 0;
	auto it = m_clouds.begin();
	for (; it != m_clouds.end(); it++, index++)
		if (index == AIndex) return *it;

	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Node> *
Mesh::getGroup<Node>(const std::string &AName)
{

	auto it = m_clouds.begin();
	for (; it != m_clouds.end(); it++)
		if ((*it)->name() == AName) return *it;

	std::stringstream mess;
	mess << "There is no node group with name " << AName;
	throw GMDSException(mess.str());
	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> unsigned int
Mesh::getNbGroups<Node>() const
{
	return m_clouds.size();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Node>
Mesh::groups_begin<Node>()
{
	return m_clouds.begin();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Node>
Mesh::groups_end<Node>()
{
	return m_clouds.end();
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Edge> *
Mesh::newGroup<Edge>(const std::string &AName)
{
	auto it = m_lines.begin();
	for (; it != m_lines.end(); it++)
		if ((*it)->name() == AName) {
			const std::string mess = "An edge group named " + AName + " already exists !";
			throw GMDSException(mess);
		}

	auto *c = new CellGroup<Edge>(this, AName);
	m_lines.push_back(c);
	return m_lines.back();
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::deleteGroup<Edge>(CellGroup<Edge> *AG)
{
	auto it = m_lines.begin();
	for (; it != m_lines.end(); it++)
		if (AG == *it) {
			m_lines.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Edge> *
Mesh::getGroup<Edge>(const unsigned int AIndex)
{
	if (AIndex >= m_lines.size()) {
		std::stringstream mess;
		mess << "There is no edge group indexed " << AIndex;
		throw GMDSException(mess.str());
	}
	unsigned int index = 0;
	auto it = m_lines.begin();
	for (; it != m_lines.end(); it++, index++)
		if (index == AIndex) return *it;

	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Edge> *
Mesh::getGroup<Edge>(const std::string &AName)
{

	auto it = m_lines.begin();
	for (; it != m_lines.end(); it++)
		if ((*it)->name() == AName) return *it;

	std::stringstream mess;
	mess << "There is no edge group with name " << AName;
	throw GMDSException(mess.str());
	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> unsigned int
Mesh::getNbGroups<Edge>() const
{
	return m_lines.size();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Edge>
Mesh::groups_begin<Edge>()
{
	return m_lines.begin();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Edge>
Mesh::groups_end<Edge>()
{
	return m_lines.end();
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Face> *
Mesh::newGroup<Face>(const std::string &AName)
{
	auto it = m_surfaces.begin();
	for (; it != m_surfaces.end(); it++)
		if ((*it)->name() == AName) {
			const std::string mess = "A face group named " + AName + " already exists !";
			throw GMDSException(mess);
		}

	auto *c = new CellGroup<Face>(this, AName);
	m_surfaces.push_back(c);
	return m_surfaces.back();
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::deleteGroup<Face>(CellGroup<Face> *AG)
{
	auto it = m_surfaces.begin();
	for (; it != m_surfaces.end(); it++)
		if (AG == *it) {
			m_surfaces.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Face> *
Mesh::getGroup<Face>(const unsigned int AIndex)
{
	if (AIndex >= m_surfaces.size()) {
		std::stringstream mess;
		mess << "There is no face group indexed " << AIndex;
		throw GMDSException(mess.str());
	}
	unsigned int index = 0;
	auto it = m_surfaces.begin();
	for (; it != m_surfaces.end(); it++, index++)
		if (index == AIndex) return *it;

	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Face> *
Mesh::getGroup<Face>(const std::string &AName)
{

	auto it = m_surfaces.begin();
	for (; it != m_surfaces.end(); it++)
		if ((*it)->name() == AName) return *it;

	std::stringstream mess;
	mess << "There is no face group with name " << AName;
	throw GMDSException(mess.str());
	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> unsigned int
Mesh::getNbGroups<Face>() const
{
	return m_surfaces.size();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Face>
Mesh::groups_begin<Face>()
{
	return m_surfaces.begin();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Face>
Mesh::groups_end<Face>()
{
	return m_surfaces.end();
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Region> *
Mesh::newGroup<Region>(const std::string &AName)
{
	auto it = m_volumes.begin();
	for (; it != m_volumes.end(); it++)
		if ((*it)->name() == AName) {
			const std::string mess = "A region group named " + AName + " already exists !";
			throw GMDSException(mess);
		}

	auto *c = new CellGroup<Region>(this, AName);
	m_volumes.push_back(c);
	return m_volumes.back();
}
/*----------------------------------------------------------------------------*/
template<> void
Mesh::deleteGroup<Region>(CellGroup<Region> *AG)
{
	auto it = m_volumes.begin();
	for (; it != m_volumes.end(); it++)
		if (AG == *it) {
			m_volumes.erase(it);
			return;
		}
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Region> *
Mesh::getGroup<Region>(const unsigned int AIndex)
{
	if (AIndex >= m_volumes.size()) {
		std::stringstream mess;
		mess << "There is no region group indexed " << AIndex;
		throw GMDSException(mess.str());
	}
	unsigned int index = 0;
	auto it = m_volumes.begin();
	for (; it != m_volumes.end(); it++, index++)
		if (index == AIndex) return *it;

	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> CellGroup<Region> *
Mesh::getGroup<Region>(const std::string &AName)
{

	auto it = m_volumes.begin();
	for (; it != m_volumes.end(); it++)
		if ((*it)->name() == AName) return *it;

	std::stringstream mess;
	mess << "There is no region group with name " << AName;
	throw GMDSException(mess.str());
	// useless but allows to remove compilation warnings
	return *it;
}
/*----------------------------------------------------------------------------*/
template<> unsigned int
Mesh::getNbGroups<Region>() const
{
	return m_volumes.size();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Region>
Mesh::groups_begin<Region>()
{
	return m_volumes.begin();
}
/*----------------------------------------------------------------------------*/
template<> Mesh::group_iterator<Region>
Mesh::groups_end<Region>()
{
	return m_volumes.end();
}
/*------------------------------------------------------------------------*/
void
Mesh::deleteNode(const Node &AN)
{
	m_nodes_container->remove(AN.id());
	m_node_variable_manager->removeEntry(AN.id());
}
/*------------------------------------------------------------------------*/
void
Mesh::deleteNode(TCellID n)
{
	m_nodes_container->remove(n);
	m_node_variable_manager->removeEntry(n);
}
void
Mesh::deleteEdge(const Edge &e)
{
	m_edges_container->remove(e.id());
	m_edge_variable_manager->removeEntry(e.id());
}
void
Mesh::deleteEdge(TCellID e)
{
	m_edges_container->remove(e);
	m_edge_variable_manager->removeEntry(e);
}
void
Mesh::deleteFace(const Face &f)
{
	m_faces_container->remove(f.id());
	m_face_variable_manager->removeEntry(f.id());
}
void
Mesh::deleteFace(TCellID f)
{
	m_faces_container->remove(f);
	m_face_variable_manager->removeEntry(f);
}
void
Mesh::deleteRegion(const Region &AR)
{
	m_regions_container->remove(AR.id());
	m_region_variable_manager->removeEntry(AR.id());
}
void
Mesh::deleteRegion(TCellID AR)
{
	m_regions_container->remove(AR);
	m_region_variable_manager->removeEntry(AR);
}

/*------------------------------------------------------------------------*/
void
Mesh::serialize(std::ostream &AStr)
{
	int model = m_model.getDef();
	AStr.write((char *) &model, sizeof(int));
	m_nodes_container->serialize(AStr);
	m_edges_container->serialize(AStr);
	m_faces_container->serialize(AStr);
	m_regions_container->serialize(AStr);

	m_node_variable_manager->serialize(AStr);
	m_edge_variable_manager->serialize(AStr);
	m_face_variable_manager->serialize(AStr);
	m_region_variable_manager->serialize(AStr);

	/*	for(std::list<cloud>::iterator it= m_clouds.begin(); it!=m_clouds.end();it++)
	      it->serialize(AStr);
	   for(std::list<line>::iterator it= m_lines.begin(); it!=m_lines.end();it++)
	      it->serialize(AStr);
	   for(std::list<surface>::iterator it= m_surfaces.begin(); it!=m_surfaces.end();it++)
	      it->serialize(AStr);
	   for(std::list<volume>::iterator it= m_volumes.begin(); it!=m_volumes.end();it++)
	      it->serialize(AStr);
	*/
	/** classification has not to serialized, it is done during variable managers'
	 *  serializations. Idem for m_marks.
	 */
	AStr.write((char *) &m_usedMarks_nodes, sizeof(int32_t));
	AStr.write((char *) &m_usedMarks_edges, sizeof(int32_t));
	AStr.write((char *) &m_usedMarks_faces, sizeof(int32_t));
	AStr.write((char *) &m_usedMarks_regions, sizeof(int32_t));

	AStr.write((char *) &m_maskMarks_nodes, sizeof(int32_t));
	AStr.write((char *) &m_maskMarks_edges, sizeof(int32_t));
	AStr.write((char *) &m_maskMarks_faces, sizeof(int32_t));
	AStr.write((char *) &m_maskMarks_regions, sizeof(int32_t));

	AStr.write((char *) &m_marks_nodes[0], 32 * sizeof(TInt));
	AStr.write((char *) &m_marks_edges[0], 32 * sizeof(TInt));
	AStr.write((char *) &m_marks_faces[0], 32 * sizeof(TInt));
	AStr.write((char *) &m_marks_regions[0], 32 * sizeof(TInt));

	AStr.write((char *) &m_nbUsedMarks_nodes, sizeof(TInt));
	AStr.write((char *) &m_nbUsedMarks_edges, sizeof(TInt));
	AStr.write((char *) &m_nbUsedMarks_faces, sizeof(TInt));
	AStr.write((char *) &m_nbUsedMarks_regions, sizeof(TInt));
}
/*------------------------------------------------------------------------*/
void
Mesh::unserialize(std::istream &AStr)
{

	int def_model;
	AStr.read((char *) &def_model, sizeof(int));
	m_model = MeshModel(def_model);
	m_nodes_container->unserialize(AStr);
	m_edges_container->unserialize(AStr);
	m_faces_container->unserialize(AStr);
	m_regions_container->unserialize(AStr);

	m_node_variable_manager->unserialize(AStr);
	m_edge_variable_manager->unserialize(AStr);
	m_face_variable_manager->unserialize(AStr);
	m_region_variable_manager->unserialize(AStr);

	/*	for(std::list<cloud>::iterator it= m_clouds.begin(); it!=m_clouds.end();it++)
	      it->unserialize(AStr);
	   for(std::list<line>::iterator it= m_lines.begin(); it!=m_lines.end();it++)
	      it->unserialize(AStr);
	   for(std::list<surface>::iterator it= m_surfaces.begin(); it!=m_surfaces.end();it++)
	      it->unserialize(AStr);
	   for(std::list<volume>::iterator it= m_volumes.begin(); it!=m_volumes.end();it++)
	      it->unserialize(AStr);
	*/
	/** classification has not to serialized, it is done during variable managers'
	 *  serializations. Idem for m_marks.
	 */
	AStr.read((char *) &m_usedMarks_nodes, sizeof(int32_t));
	AStr.read((char *) &m_usedMarks_edges, sizeof(int32_t));
	AStr.read((char *) &m_usedMarks_faces, sizeof(int32_t));
	AStr.read((char *) &m_usedMarks_regions, sizeof(int32_t));

	AStr.read((char *) &m_maskMarks_nodes, sizeof(int32_t));
	AStr.read((char *) &m_maskMarks_edges, sizeof(int32_t));
	AStr.read((char *) &m_maskMarks_faces, sizeof(int32_t));
	AStr.read((char *) &m_maskMarks_regions, sizeof(int32_t));

	AStr.read((char *) &m_marks_nodes[0], 32 * sizeof(TInt));
	AStr.read((char *) &m_marks_edges[0], 32 * sizeof(TInt));
	AStr.read((char *) &m_marks_faces[0], 32 * sizeof(TInt));
	AStr.read((char *) &m_marks_regions[0], 32 * sizeof(TInt));

	AStr.read((char *) &m_nbUsedMarks_nodes, sizeof(TInt));
	AStr.read((char *) &m_nbUsedMarks_edges, sizeof(TInt));
	AStr.read((char *) &m_nbUsedMarks_faces, sizeof(TInt));
	AStr.read((char *) &m_nbUsedMarks_regions, sizeof(TInt));
}
/*----------------------------------------------------------------------------*/
void
Mesh::getCommonNodes(const Face &AF1, const Face &AF2, std::vector<TCellID> &ANodes)
{
	ANodes.clear();
	std::vector<TCellID> nodes1 = AF1.getIDs<Node>();
	std::vector<TCellID> nodes2 = AF2.getIDs<Node>();

	bool found;
	for (unsigned int &i1 : nodes1) {
		found = false;
		for (unsigned int i2 = 0; i2 < nodes2.size() && !found; i2++)
			if (i1 == nodes2[i2]) {
				ANodes.push_back(i1);
				found = true;
			}
	}
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>
Mesh::getCommonNodes(const Face &AF1, const Face &AF2)
{
	std::vector<TCellID> nodes;
	getCommonNodes(AF1, AF2, nodes);
	return nodes;
}

/*----------------------------------------------------------------------------*/
void
Mesh::getCommonFaces(const Node &AN1, const Node &AN2, std::vector<TCellID> &AVec)
{
	AVec.clear();
	std::vector<TCellID> set1 = AN1.getIDs<Face>();
	std::vector<TCellID> set2 = AN2.getIDs<Face>();

	bool found;
	for (unsigned int &i1 : set1) {
		found = false;
		for (unsigned int i2 = 0; i2 < set2.size() && !found; i2++)
			if (i1 == set2[i2]) {
				AVec.push_back(i1);
				found = true;
			}
	}
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>
Mesh::getCommonFaces(const Node &AN1, const Node &AN2)
{
	std::vector<TCellID> ids;
	getCommonFaces(AN1, AN2, ids);
	return ids;
}
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TInt
Mesh::getPartID() const
{
	return m_part_id;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#endif     // GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
