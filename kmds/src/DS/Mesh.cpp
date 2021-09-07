/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Mesh.cpp
 *  \author  F. LEDOUX
 *  \date    03/07/2017
 */
/*----------------------------------------------------------------------------*/
#include <algorithm>
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*------------------------------------------------------------------------*/
Mesh::Mesh()
 : m_N(this, DefaultContainerSize)
 , m_E(this, DefaultContainerSize)
 , m_F(this, DefaultContainerSize)
 , m_R(this, DefaultContainerSize)
{
}
/*------------------------------------------------------------------------*/
Mesh::~Mesh()
{
}
/*------------------------------------------------------------------------*/
void
Mesh::updateNodeCapacity(const TSize ASize)
{
        // we resize the node container
        m_N.resize(ASize);
        // and all the associated variables
        m_var_manager[0].resize(ASize);
}
/*------------------------------------------------------------------------*/
void
Mesh::updateEdgeCapacity(const TSize ASize)
{
    // we resize the edge container
    m_E.resize(ASize);
    // and all the associated variables
    m_var_manager[1].resize(ASize);
}
/*------------------------------------------------------------------------*/
void
Mesh::updateFaceCapacity(const TSize ASize)
{
        // we resize the face container
        m_F.resize(ASize);
        // and all the associated variables
        m_var_manager[2].resize(ASize);
}
/*------------------------------------------------------------------------*/
void
Mesh::updateRegionCapacity(const TSize ASize)
{
        // we resize the region container
        m_R.resize(ASize);
        // and all the associated variables
        m_var_manager[3].resize(ASize);
}
/*------------------------------------------------------------------------*/
void
Mesh::doubleNodeCapacity()
{
        updateNodeCapacity(2 * getNodeCapacity());
}
/*------------------------------------------------------------------------*/
void
Mesh::doubleEdgeCapacity()
{
    updateEdgeCapacity(2 * getEdgeCapacity());
}
/*------------------------------------------------------------------------*/
void
Mesh::doubleFaceCapacity()
{
        updateFaceCapacity(2 * getFaceCapacity());
}
/*------------------------------------------------------------------------*/
void
Mesh::doubleRegionCapacity()
{
        updateRegionCapacity(2 * getRegionCapacity());
}
/*------------------------------------------------------------------------*/
void
Mesh::removeAllNodes()
{
        m_N.removeAll();
}
/*------------------------------------------------------------------------*/
void
Mesh::removeAllEdges()
{
        m_E.removeAll();
}
/*------------------------------------------------------------------------*/
void
Mesh::removeAllFaces()
{
        m_F.removeAll();
}
/*------------------------------------------------------------------------*/
void
Mesh::removeAllRegions()
{
        m_R.removeAll();
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addNodes(const int ANb)
{
        TCellID i = m_N.add(ANb);
        m_var_manager[0].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addNode()
{
        TCellID i = m_N.add();
        m_var_manager[0].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newNode(const TCoord AX, const TCoord AY, const TCoord AZ)
{
        TCellID i = m_N.newNode(AX, AY, AZ);
        m_var_manager[0].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newNode(const gmds::math::Point APt)
{
    TCellID i = m_N.newNode(APt.X(), APt.Y(), APt.Z());
    m_var_manager[0].initializeVariables(i);
    return i;
}
/*------------------------------------------------------------------------*/
void
Mesh::removeNode(const TCellID AId)
{
        m_N.remove(AId);
}
/*------------------------------------------------------------------------*/
Node
Mesh::getNode(const TCellID AID) const
{
        return Node(this, AID);
}
/*------------------------------------------------------------------------*/
bool
Mesh::hasNode(const TCellID AID) const
{
        return m_N.has(AID);
}
/*----------------------------------------------------------------------------*/
TSize
Mesh::getNodeCapacity() const
{
        return m_N.capacity();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNodeSupID() const
{
    return m_N.top() - 1;
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbNodes() const
{
        return m_N.nbCells();
}
/*------------------------------------------------------------------------*/
void
Mesh::getNodeLocation(const TCellID AId, TCoord& AX, TCoord& AY, TCoord& AZ) const
{
        m_N.get(AId, AX, AY, AZ);
}
/*------------------------------------------------------------------------*/
gmds::math::Point
Mesh::getNodeLocation(const TCellID AId) const
{
    TCoord xyz[3];
    m_N.get(AId, xyz[0], xyz[1], xyz[2]);
    return gmds::math::Point(xyz[0], xyz[1], xyz[2]);
}
/*------------------------------------------------------------------------*/
void
Mesh::setNodeLocation(const TCellID AId, const TCoord AX, const TCoord AY, const TCoord AZ)
{
        m_N.set(AId, AX, AY, AZ);
}
/*------------------------------------------------------------------------*/
void
Mesh::setNodeLocation(const TCellID AId, const gmds::math::Point APt)
{
    m_N.set(AId, APt.X(), APt.Y(), APt.Z());
}
/*----------------------------------------------------------------------------*/
TCellID
Mesh::addEdges(const int ANb)
{
        TCellID i = m_E.add(ANb);
        m_var_manager[1].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addEdge()
{
        TCellID i = m_E.add();
        m_var_manager[1].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newEdge(const TCellID AN1, const TCellID AN2)
{
        TCellID i = m_E.newEdge(AN1, AN2);
        m_var_manager[1].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newEdge_unsafe(const TCellID AN1, const TCellID AN2)
{
        TCellID i = m_E.newEdge_unsafe(AN1, AN2);
        m_var_manager[1].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
void
Mesh::removeEdge(const TCellID AId)
{
        m_E.remove(AId);
}
/*------------------------------------------------------------------------*/
Edge
Mesh::getEdge(const TCellID AID) const
{
        return Edge(this, AID);
}
/*------------------------------------------------------------------------*/
bool
Mesh::hasEdge(const TCellID AID) const
{
        return m_E.has(AID);
}
/*----------------------------------------------------------------------------*/
TSize
Mesh::getEdgeCapacity() const
{
    return m_E.capacity();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getEdgeSupID() const
{
    return m_E.top() - 1;
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbEdges() const
{
        return m_E.nbCells();
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addTriangles(const int ANb)
{
        TCellID i = m_F.addTriangles(ANb);
        m_var_manager[2].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addQuads(const int ANb)
{
        TCellID i = m_F.addQuads(ANb);
        m_var_manager[2].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addPentagons(const int ANb)
{
        TCellID i = m_F.addPentagons(ANb);
        m_var_manager[2].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addTriangle()
{
        TCellID i = m_F.addTriangle();
        m_var_manager[2].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addQuad()
{
        TCellID i = m_F.addQuad();
        m_var_manager[2].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addPentagon()
{
        TCellID i = m_F.addPentagon();
        m_var_manager[2].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newTriangle(const TCellID AN1, const TCellID AN2, const TCellID AN3)
{
        TCellID i = m_F.newTriangle(AN1, AN2, AN3);
        m_var_manager[2].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newFace(const TCellID AN1, const TCellID AN2, const TCellID AN3)
{
        return newTriangle(AN1, AN2, AN3);
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newQuad(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4)
{
        TCellID i = m_F.newQuad(AN1, AN2, AN3, AN4);
        m_var_manager[2].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newFace(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4)
{
        return newQuad(AN1, AN2, AN3, AN4);
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newPentagon(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5)
{
        TCellID i = m_F.newPentagon(AN1, AN2, AN3, AN4, AN5);
        m_var_manager[2].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newFace(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5)
{
        return newPentagon(AN1, AN2, AN3, AN4, AN5);
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newFace(const FakeFace AFF)
{
        std::vector<TCellID> ids = AFF.node_ids();

        switch(ids.size())
        {
                case 3 :
                        return newTriangle(ids[0], ids[1], ids[2]);
                case 4 :
                        return newQuad(ids[0], ids[1], ids[2], ids[3]);
                case 5 :
                        return newPentagon(ids[0], ids[1], ids[2], ids[3], ids[4]);
                default:
                        std::cerr<<"Mesh::newFace ids.size() "<<ids.size()<<std::endl;
                        throw KException("Mesh::newFace This type of fake face is not handled");
        }
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newFace(const TCellID* AIDs, const TSize ASize)
{
        switch(ASize)
        {
                case 3 :
                        return newTriangle(AIDs[0], AIDs[1], AIDs[2]);
                case 4 :
                        return newQuad(AIDs[0], AIDs[1], AIDs[2], AIDs[3]);
                case 5 :
                        return newPentagon(AIDs[0], AIDs[1], AIDs[2], AIDs[3], AIDs[4]);
                default:
                    std::cerr<<"Mesh::newFace ASize "<<ASize<<std::endl;
                        throw KException("Mesh::newFace This type of fake face is not handled");
        }
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newFace(const Kokkos::View<TCellID*>& AN)
{
        switch(AN.extent(0))
        {
                case 3 :
                        return newTriangle(AN(0), AN(1), AN(2));
                case 4 :
                        return newQuad(AN(0), AN(1), AN(2), AN(3));
                case 5 :
                        return newPentagon(AN(0), AN(1), AN(2), AN(3), AN(4));
                default:
                        std::cerr<<"Mesh::newFace AN.extent(0) "<<AN.extent(0)<<std::endl;
                        throw KException("Mesh::newFace This type of fake face is not handled");
        }
}
/*------------------------------------------------------------------------*/
void
Mesh::removeFace(const TCellID AId)
{
        m_F.remove(AId);
}
/*------------------------------------------------------------------------*/
Face
Mesh::getFace(const TCellID AId) const
{
        return Face(this, AId);
}
/*------------------------------------------------------------------------*/
bool
Mesh::hasFace(const TCellID AId) const
{
        return m_F.has(AId);
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getFaceCapacity() const
{
        return m_F.capacity();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getFaceSupID() const
{
    return m_F.top() - 1;
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbFaces() const
{
        return m_F.nbCells();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbTriangles() const
{
        return m_F.nbTriangles();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbQuads() const
{
        return m_F.nbQuads();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbPentagons() const
{
        return m_F.nbPentagons();
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addTetrahedra(const int ANb)
{
        TCellID i = m_R.addTetrahedra(ANb);
        m_var_manager[3].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addHexahedra(const int ANb)
{
        TCellID i = m_R.addHexahedra(ANb);
        m_var_manager[3].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addPyramids(const int ANb)
{
        TCellID i = m_R.addPyramids(ANb);
        m_var_manager[3].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addPrism3s(const int ANb)
{
        TCellID i = m_R.addPrism3s(ANb);
        m_var_manager[3].initializeVariables(i, ANb);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addTetrahedron()
{
        TCellID i = m_R.addTetrahedron();
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addHexahedron()
{
        TCellID i = m_R.addHexahedron();
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addPyramid()
{
        TCellID i = m_R.addPyramid();
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::addPrism3()
{
        TCellID i = m_R.addPrism3();
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newTetrahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4)
{
        TCellID i = m_R.newTetrahedron(AN1, AN2, AN3, AN4);
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4)
{
        return newTetrahedron(AN1, AN2, AN3, AN4);
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newHexahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                const TCellID AN5, const TCellID AN6, const TCellID AN7, const TCellID AN8)
{
        TCellID i = m_R.newHexahedron(AN1, AN2, AN3, AN4, AN5, AN6, AN7, AN8);
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                const TCellID AN5, const TCellID AN6, const TCellID AN7, const TCellID AN8)
{
        return newHexahedron(AN1, AN2, AN3, AN4, AN5, AN6, AN7, AN8);
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newPyramid(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5)
{
        TCellID i = m_R.newPyramid(AN1, AN2, AN3, AN4, AN5);
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5)
{
        return newPyramid(AN1, AN2, AN3, AN4, AN5);
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newPrism3(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                    const TCellID AN5, const TCellID AN6)
{
        TCellID i = m_R.newPrism3(AN1, AN2, AN3, AN4, AN5, AN6);
        m_var_manager[3].initializeVariables(i);
        return i;
}
/*------------------------------------------------------------------------*/
TCellID
Mesh::newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                const TCellID AN5, const TCellID AN6)
{
        return newPrism3(AN1, AN2, AN3, AN4, AN5, AN6);
}
/*------------------------------------------------------------------------*/
void
Mesh::removeRegion(const TCellID AId)
{
        m_R.remove(AId);
}
/*------------------------------------------------------------------------*/
Region
Mesh::getRegion(const TCellID AId) const
{
        return Region(this, AId);
}
/*------------------------------------------------------------------------*/
bool
Mesh::hasRegion(const TCellID AId) const
{
        return m_R.has(AId);
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getRegionCapacity() const
{
        return m_R.capacity();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getRegionSupID() const
{
        return m_R.top() - 1;
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbRegions() const
{
        return m_R.nbCells();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbTetrahedra() const
{
        return m_R.nbTetrahedra();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbHexahedra() const
{
        return m_R.nbHexahedra();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbPyramids() const
{
        return m_R.nbPyramids();
}
/*------------------------------------------------------------------------*/
TSize
Mesh::getNbPrism3s() const
{
        return m_R.nbPrism3s();
}
/*------------------------------------------------------------------------*/
void
Mesh::getNodeIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        struct SelectNode {
            const kmds::Mesh *m;

            kmds::GrowingView<kmds::TCellID> *result;

            // Constructor takes View by "value"; this does a shallow copy.
            SelectNode(const kmds::Mesh *m_, kmds::GrowingView<kmds::TCellID> *r_)
                    : m(m_), result(r_) {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                    if (m->hasNode(i)) {
                            result->push_back(i);
                    }
            }
        };

        Kokkos::parallel_for(this->getNodeCapacity(), SelectNode(this, ASelection));
}
/*------------------------------------------------------------------------*/
void
Mesh::getNodeIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        int nbCells = this->getNbNodes();
        Kokkos::parallel_for(nbCells, KOKKOS_LAMBDA(const int i){ASelection->set(i,i);});
        ASelection->setTop(nbCells);
}
/*------------------------------------------------------------------------*/
void
Mesh::getEdgeIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        struct SelectEdge {
            const kmds::Mesh *m;

            kmds::GrowingView<kmds::TCellID> *result;

            // Constructor takes View by "value"; this does a shallow copy.
            SelectEdge(const kmds::Mesh *m_, kmds::GrowingView<kmds::TCellID> *r_)
                    : m(m_), result(r_) {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                    if (m->hasEdge(i)) {
                            result->push_back(i);
                    }
            }
        };

        Kokkos::parallel_for(this->getEdgeCapacity(), SelectEdge(this, ASelection));
}
/*------------------------------------------------------------------------*/
void
Mesh::getEdgeIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        int nbCells = this->getNbEdges();
        Kokkos::parallel_for(nbCells, KOKKOS_LAMBDA(const int i){ASelection->set(i,i);});
        ASelection->setTop(nbCells);
}
/*------------------------------------------------------------------------*/
void
Mesh::getFaceIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        struct SelectFace {
            const kmds::Mesh *m;

            kmds::GrowingView<kmds::TCellID> *result;

            // Constructor takes View by "value"; this does a shallow copy.
            SelectFace(const kmds::Mesh *m_, kmds::GrowingView<kmds::TCellID> *r_)
                    : m(m_), result(r_) {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const {
                    if (m->hasFace(i)) {
                            result->push_back(i);
                    }
            }
        };

        Kokkos::parallel_for(this->getFaceCapacity(), SelectFace(this, ASelection));
}
/*------------------------------------------------------------------------*/
void
Mesh::getFaceIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        int nbCells = this->getNbFaces();
        Kokkos::parallel_for(nbCells, KOKKOS_LAMBDA(const int i){ASelection->set(i,i);});
        ASelection->setTop(nbCells);
}
/*------------------------------------------------------------------------*/
void
Mesh::getRegionIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        struct SelectRegion
        {
            const kmds::Mesh* m;

            kmds::GrowingView<kmds::TCellID>* result;
            // Constructor takes View by "value"; this does a shallow copy.
            SelectRegion(const kmds::Mesh* m_, kmds::GrowingView<kmds::TCellID>* r_)
                    : m(m_)
                    , result(r_)
            {
            }

            KOKKOS_INLINE_FUNCTION
            void
            operator()(int i) const
            {
                    if (m->hasRegion(i)) {
                            result->push_back(i);
                    }
            }
        };

        Kokkos::parallel_for(this->getRegionCapacity(), SelectRegion(this, ASelection));
}
/*------------------------------------------------------------------------*/
void
Mesh::getRegionIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        int nbCells = this->getNbRegions();
        Kokkos::parallel_for(nbCells, KOKKOS_LAMBDA(const int i){ASelection->set(i,i);});
        ASelection->setTop(nbCells);
}
/*------------------------------------------------------------------------*/
void
Mesh::getNodeIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        for(kmds::TCellID id=0; id<this->getNodeCapacity(); id++) {
                if(this->hasNode(id)) {
                        ASelection->push_back(id);
                }
        }
}
/*------------------------------------------------------------------------*/
void
Mesh::getEdgeIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        for(kmds::TCellID id=0; id<this->getEdgeCapacity(); id++) {
                if(this->hasEdge(id)) {
                        ASelection->push_back(id);
                }
        }
}
/*------------------------------------------------------------------------*/
void
Mesh::getFaceIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        for(kmds::TCellID id=0; id<this->getFaceCapacity(); id++) {
                if(this->hasFace(id)) {
                        ASelection->push_back(id);
                }
        }
}
/*------------------------------------------------------------------------*/
void
Mesh::getRegionIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const
{
        for(kmds::TCellID id=0; id<this->getRegionCapacity(); id++) {
                if(this->hasRegion(id)) {
                        ASelection->push_back(id);
                }
        }
}
/*------------------------------------------------------------------------*/
Connectivity*
Mesh::createConnectivity(const EMeshDefinition AD)
{
        if (m_connectivity.find(AD) != m_connectivity.end())
                throw KException("Connectivity " + Connectivity::getName(AD) + " already exists");

        // when we create a connectivity, we initialize its size depending
        // on the type of elements it refers to
        int c_size = 0;
        if (AD == N2N || AD == N2E || AD == N2F || AD == N2R) {
                c_size = getNodeCapacity();
        } else if (AD == E2F || AD == E2R || AD == E2E) {
                c_size = getEdgeCapacity();
        } else if (AD == F2E || AD == F2R || AD == F2F || AD == F2F_byN) {
                c_size = getFaceCapacity();
        } else if (AD == R2E || AD == R2F || AD == R2R || AD == R2R_byN) {
                c_size = getRegionCapacity();
        } else {
                throw KException("Mesh::createConnectivity This type of connectivity is not handled");
        }
        Connectivity* c = new Connectivity(this, AD, c_size);

        m_connectivity[AD] = c;
        return m_connectivity[AD];
}
/*------------------------------------------------------------------------*/
Connectivity*
Mesh::getConnectivity(const EMeshDefinition AD)
{
        // TODO: must be done in debug only
        if (m_connectivity.find(AD) == m_connectivity.end())
                createConnectivity(AD);

        return m_connectivity[AD];
}
/*------------------------------------------------------------------------*/
Connectivity*
Mesh::deleteConnectivity(const EMeshDefinition AD)
{
        if (m_connectivity.find(AD) == m_connectivity.end()) {
                throw KException("Mesh::deleteConnectivity connectivity not found.");
        }

        Connectivity* c = this->getConnectivity(AD);
        m_connectivity.erase(AD);
        delete c;
}
/*------------------------------------------------------------------------*/
void
Mesh::deleteVariable(const ECellType AType, std::string AName)
{
        switch (AType) {
                case KMDS_NODE: {
                        m_var_manager[0].deleteVariable(AName);
                } break;
                case KMDS_EDGE: {
                        m_var_manager[1].deleteVariable(AName);
                } break;
                case KMDS_FACE: {
                        m_var_manager[2].deleteVariable(AName);
                } break;
                case KMDS_REGION: {
                        m_var_manager[3].deleteVariable(AName);
                } break;
                default:
                        throw KException("Mesh::deleteVariable Unmanaged type of cell -> impossible to access to a variable");
        }
}
/*------------------------------------------------------------------------*/
VariableManager*
Mesh::getVariableManager(const ECellType AType)
{
        switch (AType) {
                case KMDS_NODE: {
                        return &m_var_manager[0];
                } break;
                case KMDS_EDGE: {
                        return &m_var_manager[1];
                } break;
                case KMDS_FACE: {
                        return &m_var_manager[2];
                } break;
                case KMDS_REGION: {
                        return &m_var_manager[3];
                } break;
                default:
                        throw KException("Mesh::getVariableManager Unmanaged type of cell -> impossible to access to a variable");
        }
}
/*------------------------------------------------------------------------*/