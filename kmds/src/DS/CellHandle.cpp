/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/** \file    CellHandle.cpp
 *  \author  F. LEDOUX
 *  \date    03/08/2017
 */
/*----------------------------------------------------------------------------*/
#include "KM/DS/CellHandle.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Hexahedron.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Prism3.h>
#include <gmds/math/Pyramid.h>
#include <gmds/math/Quadrilateral.h>
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
// KMDS File Headers
/*----------------------------------------------------------------------------*/
#include "KM/DS/Mesh.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*------------------------------------------------------------------------*/
Cell::Cell()
 : owner(nullptr)
 , id(NullID)
{
}
/*------------------------------------------------------------------------*/
Cell::Cell(const Mesh* AM, const TCellID AID)
 : owner(const_cast<Mesh*>(AM))
 , id(AID)
{
}
/*------------------------------------------------------------------------*/
Node::Node()
 : Cell()
{
}
/*------------------------------------------------------------------------*/
Node::Node(const Mesh* AM, const TCellID AID)
 : Cell(AM, AID)
{
}
/*------------------------------------------------------------------------*/
void
Node::setLocation(const TCoord AX, const TCoord AY, const TCoord AZ)
{
        owner->m_N.set(id, AX, AY, AZ);
}
/*------------------------------------------------------------------------*/
void
Node::getLocation(TCoord& AX, TCoord& AY, TCoord& AZ) const
{
        owner->m_N.get(id, AX, AY, AZ);
}
/*------------------------------------------------------------------------*/
void
Node::setEdges(const Kokkos::View<TCellID*>& AV)
{
        owner->getConnectivity(N2E)->set(id, AV);
}
/*------------------------------------------------------------------------*/
void
Node::setFaces(const Kokkos::View<TCellID*>& AV)
{
        owner->getConnectivity(N2F)->set(id, AV);
}
/*------------------------------------------------------------------------*/
void
Node::faceIds(Kokkos::View<TCellID*>& AV) const
{
        owner->getConnectivity(N2F)->get(id, AV);
}
/*------------------------------------------------------------------------*/
void
Node::setRegions(const Kokkos::View<TCellID*>& AV)
{
        owner->getConnectivity(N2R)->set(id, AV);
}
/*------------------------------------------------------------------------*/
void
Node::setRegions(const TCellID* AV, const int ASize)
{
        owner->getConnectivity(N2R)->set(id, AV, ASize);
}
/*------------------------------------------------------------------------*/
void
Node::regionIds(Kokkos::View<TCellID*>& AV) const
{
        owner->getConnectivity(N2R)->get(id, AV);
}
/*------------------------------------------------------------------------*/
ECellType
Node::computeType() const
{
        return KMDS_NODE;
}
/*------------------------------------------------------------------------*/
gmds::math::Point
Node::getPoint() const
{
        TCoord coords[3];
        this->getLocation(coords[0], coords[1], coords[2]);

        return gmds::math::Point (coords[0], coords[1], coords[2]);
}
/*------------------------------------------------------------------------*/
Edge::Edge()
 : Cell()
{
}
/*------------------------------------------------------------------------*/
Edge::Edge(const Mesh* AM, const TCellID AID)
 : Cell(AM, AID)
{
}
/*------------------------------------------------------------------------*/
Face::Face()
 : Cell()
{
}
/*------------------------------------------------------------------------*/
Face::Face(const Mesh* AM, const TCellID AID)
 : Cell(AM, AID)
{
}
/*------------------------------------------------------------------------*/
ECellType
Edge::computeType() const
{
        return KMDS_EDGE;
}
/*------------------------------------------------------------------------*/
void
Edge::nodeIds(Kokkos::View<TCellID*>& AV) const
{
        owner->m_E.get(id, AV);
}
/*------------------------------------------------------------------------*/
void
Edge::nodeIds(TCellID AIds[2]) const
{
        owner->m_E.get(id, AIds[0], AIds[1]);
}
/*------------------------------------------------------------------------*/
void
Edge::setNodes(TCellID AN0, TCellID AN1)
{
        owner->m_E.set(id, AN0, AN1);
}
/*------------------------------------------------------------------------*/
TCoord
Edge::surfvol() const
{
    Kokkos::View<TCellID*> nodes;
    this->nodeIds(nodes);

    return (owner->getNodeLocation(nodes[0])).distance(owner->getNodeLocation(nodes[1]));
}
/*------------------------------------------------------------------------*/
ECellType
Face::computeType() const
{
        ECellType type;
        TInt32 nb_nodes = getNbNodes();
        if (nb_nodes == 3) {
                type = KMDS_TRIANGLE;
        } else if (nb_nodes == 4) {
                type = KMDS_QUAD;
        } else if (nb_nodes == 5) {
                type = KMDS_PENTAGON;
        } else {
                type = KMDS_UNDEF;
        }
        return type;
}
/*------------------------------------------------------------------------*/
void
Face::setNodes(const Kokkos::View<TCellID*>& AN)
{
        owner->m_F.set(id, AN);
}
/*------------------------------------------------------------------------*/
void
Face::setNodes(TCellID* AIds, int ASize)
{
        owner->m_F.setUnsafe(id, AIds, ASize);
}
/*------------------------------------------------------------------------*/
void
Face::setUnsafeNodes(const Kokkos::View<TCellID*>& AN)
{
        owner->m_F.setUnsafe(id, AN);
}
/*------------------------------------------------------------------------*/
void
Face::nodeIds(Kokkos::View<TCellID*>& AV) const
{
        owner->m_F.get(id, AV);
}
/*------------------------------------------------------------------------*/
void
Face::nodeIds(TCellID AIds[5], int* ASize) const
{
        owner->m_F.get(id, AIds, ASize);
}
/*------------------------------------------------------------------------*/
TSize
Face::getNbNodes() const
{
        return owner->m_F.getNbNodes(id);
}
/*------------------------------------------------------------------------*/
void
Face::nodes(Kokkos::View<Node*>& AV) const
{
        Kokkos::View<TCellID*> ids;
        nodeIds(ids);
        int nb_nodes = ids.size();
        for (auto i = 0; i < nb_nodes + 1; i++) {
                AV(i).owner = owner;
                AV(i).id = ids(i);
        }
}
/*------------------------------------------------------------------------*/
void
Face::nodes(Node ANodes[5], int* ASize) const
{
        TCellID ids[5];
        nodeIds(ids,ASize);
        int nb_nodes = *ASize;
        for (auto i = 0; i < nb_nodes + 1; i++) {

                ANodes[i].owner = owner;
                ANodes[i].id = ids[i];
        }
}
/*------------------------------------------------------------------------*/
Kokkos::View<Node*>
Face::nodes() const
{
        Kokkos::View<TCellID*> ids;
        nodeIds(ids);
        int nb_nodes = ids.size();
        Kokkos::View<Node*> AV("N", nb_nodes);
        for (auto i = 0; i < nb_nodes + 1; i++) {
                AV(i).owner = owner;
                AV(i).id = ids(i);
        }
        return AV;
}
/*----------------------------------------------------------------------------*/
TCoord
Face::surfvol() const
{
        TCoord surfvol=0.0;

        Node n[5];
        int nbNodes;
        this->nodes(n, &nbNodes);

        switch (this->computeType()) {
                case KMDS_TRIANGLE:
                        surfvol = gmds::math::Triangle(n[0].getPoint(),
                                                       n[1].getPoint(),
                                                       n[2].getPoint()).area();
                break;
                case KMDS_QUAD:
                        surfvol = gmds::math::Quadrilateral(n[0].getPoint(),
                                                            n[1].getPoint(),
                                                            n[2].getPoint(),
                                                            n[3].getPoint()).area();
                break;
                default:
                        throw KException("Face::surfvol can not be computed for this cell type.");
        }

        return surfvol;
}
/*----------------------------------------------------------------------------*/
TCoord
Face::scaledJacobian() const
{
        TCoord sj=0.0;

        Node n[5];
        int nbNodes;
        this->nodes(n, &nbNodes);

        switch (this->computeType()) {
                case KMDS_TRIANGLE:
                        sj = gmds::math::Triangle(n[0].getPoint(),
                                                       n[1].getPoint(),
                                                       n[2].getPoint()).computeScaledJacobian2D();
                break;
                case KMDS_QUAD:
                        sj = gmds::math::Quadrilateral(n[0].getPoint(),
                                                            n[1].getPoint(),
                                                            n[2].getPoint(),
                                                            n[3].getPoint()).computeScaledJacobian2D();
                break;
                default:
                        throw KException("Face::scaledJacobian can not be computed for this cell type.");
        }

        return sj;
}
/*----------------------------------------------------------------------------*/
TCoord
Face::scaledJacobianAt(kmds::TCellID ANId, gmds::math::Point APt) const
{
    TCoord sj=0.0;

    Node n[5];
    int nbNodes;
    this->nodes(n, &nbNodes);

    // identify the local index of the node
    int nid_l = -1;
    for(int i_n=0; i_n<nbNodes; i_n++) {
        if(ANId == n[i_n].id) {
            nid_l = i_n;
            break;
        }
    }

    switch (this->computeType()) {
        case KMDS_QUAD:

            switch (nid_l) {

                case 0:
                sj = gmds::math::Quadrilateral(APt,
                                               n[1].getPoint(),
                                               n[2].getPoint(),
                                               n[3].getPoint()).computeScaledJacobian2DAt(nid_l);
                break;
                case 1:
                    sj = gmds::math::Quadrilateral(n[0].getPoint(),
                                                   APt,
                                                   n[2].getPoint(),
                                                   n[3].getPoint()).computeScaledJacobian2DAt(nid_l);
                    break;
                case 2:
                    sj = gmds::math::Quadrilateral(n[0].getPoint(),
                                                   n[1].getPoint(),
                                                   APt,
                                                   n[3].getPoint()).computeScaledJacobian2DAt(nid_l);
                    break;
                case 3:
                    sj = gmds::math::Quadrilateral(n[0].getPoint(),
                                                   n[1].getPoint(),
                                                   n[2].getPoint(),
                                                   APt).computeScaledJacobian2DAt(nid_l);
                    break;
            }
            break;
        default:
            throw KException("Face::scaledJacobian can not be computed for this cell type.");
    }

    return sj;
}
/*----------------------------------------------------------------------------*/
gmds::math::Point
Face::midpoint() const
{

        gmds::math::Point pt = 0.0;

        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);
        int nbNodes = nodes.size();

        for(int i=0; i<nodes.size(); i++) {

                pt = pt + owner->getNode(nodes[i]).getPoint();
        }

        pt = 1./nbNodes * pt;

        return pt;
}
/*----------------------------------------------------------------------------*/
bool
Face::isEdgeOrientedOutward(Kokkos::View<kmds::TCellID *> &AIDs) const
{
        assert(AIDs.extent(0) == 2);

        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);

        int nbNodes = nodes.size();

        for(int i=0; i<nodes.size(); i++) {
                // find the first node
                if(nodes(i) == AIDs(0)) {
                        if(nodes((i + 1) % nodes.extent(0)) == AIDs(1)) {
                                return true;
                        } else {
                                if(nodes((i - 1 + nodes.extent(0)) % nodes.extent(0)) == AIDs(1)) {
                                        return false;
                                } else {
                                        throw KException("Face::isEdgeOrientedOutward edge not found.");
                                }
                        }
                }
        }

}
/*----------------------------------------------------------------------------*/
std::vector<FakeEdge>
Face::getFakeEdges() const
{
        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);

        TSize nbNodes = nodes.size();
        std::vector<FakeEdge> fakeEdges(nbNodes);

        for(unsigned int i=0; i<nbNodes; i++) {
                fakeEdges[i] = FakeEdge(nodes[i], nodes[(i+1)%nodes.size()]);
        }

        return fakeEdges;
}
/*----------------------------------------------------------------------------*/
void
Face::getFacetInfo(Kokkos::View<kmds::TCellID *> &AIDs, int &numFacet, bool &isOutward, int &offset) const
{
        // offset makes no sense with edges
        offset = -1;

        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);

        int nbNodes = nodes.size();

        for(int i=0; i<nodes.size(); i++) {
                // find the first node
                if(nodes(i) == AIDs(0)) {
                        if(nodes((i + 1) % nodes.extent(0)) == AIDs(1)) {
                                isOutward = true;
                                numFacet = i;
                                return;
                        } else {
                                if(nodes((i - 1 + nodes.extent(0)) % nodes.extent(0)) == AIDs(1)) {
                                        isOutward = false;
                                        numFacet = (i - 1 + nbNodes) % nodes.extent(0);
                                        return;
                                } else {
                                        throw KException("Face::getFacetInfo facet not found.");
                                }
                        }
                }
        }
}

/*----------------------------------------------------------------------------*/
void
Face::getEdgeInfo(Kokkos::View<TCellID *> &AIDs, int &numFacet, bool &isIJK) const {
    Kokkos::View<TCellID *> nodes;
    this->nodeIds(nodes);

    if ((AIDs[0] == nodes[0]) && (AIDs[1] == nodes[1])) {
        numFacet = 0;
        isIJK = true;
        return;
    }
    if ((AIDs[1] == nodes[0]) && (AIDs[0] == nodes[1])) {
        numFacet = 0;
        isIJK = false;
        return;
    }
    if ((AIDs[0] == nodes[1]) && (AIDs[1] == nodes[2])) {
        numFacet = 1;
        isIJK = true;
        return;
    }
    if ((AIDs[1] == nodes[1]) && (AIDs[0] == nodes[2])) {
        numFacet = 1;
        isIJK = false;
        return;
    }
    if ((AIDs[0] == nodes[2]) && (AIDs[1] == nodes[3])) {
        numFacet = 2;
        isIJK = false;
        return;
    }
    if ((AIDs[1] == nodes[2]) && (AIDs[0] == nodes[3])) {
        numFacet = 2;
        isIJK = true;
        return;
    }
    if ((AIDs[0] == nodes[3]) && (AIDs[1] == nodes[0])) {
        numFacet = 3;
        isIJK = false;
        return;
    }
    if ((AIDs[1] == nodes[3]) && (AIDs[0] == nodes[0])) {
        numFacet = 3;
        isIJK = true;
        return;
    }
}

/*----------------------------------------------------------------------------*/
void
Face::computeBoundingBox(kmds::TCoord minXYZ[3], kmds::TCoord maxXYZ[3]) const
{
    Kokkos::View<kmds::TCellID*> nids;
    this->nodeIds(nids);

    minXYZ[0] = HUGE_VALF;
    minXYZ[1] = HUGE_VALF;
    minXYZ[2] = HUGE_VALF;

    maxXYZ[0] = -HUGE_VALF;
    maxXYZ[1] = -HUGE_VALF;
    maxXYZ[2] = -HUGE_VALF;

    for(int i=0; i<nids.size(); i++) {
        TCoord xyz[3];

        this->owner->getNodeLocation(nids[i], xyz[0], xyz[1], xyz[2]);

        minXYZ[0] = std::min(minXYZ[0], xyz[0]);
        minXYZ[1] = std::min(minXYZ[1], xyz[1]);
        minXYZ[2] = std::min(minXYZ[2], xyz[2]);

        maxXYZ[0] = std::max(maxXYZ[0], xyz[0]);
        maxXYZ[1] = std::max(maxXYZ[1], xyz[1]);
        maxXYZ[2] = std::max(maxXYZ[2], xyz[2]);
    }
}
/*----------------------------------------------------------------------------*/
gmds::math::Triangle
Face::getTriangle(kmds::TCellID AID) const
{
    gmds::math::Point pts[3];
    pts[0] = this->owner->getNodeLocation(AID);

    Kokkos::View<kmds::TCellID*> nids;
    this->nodeIds(nids);

    for(int i=0; i<nids.size(); i++) {
        if(AID == nids[i]) {
            pts[1] = this->owner->getNodeLocation(nids[(i+1)%nids.size()]);
            pts[2] = this->owner->getNodeLocation(nids[(i-1 + nids.size())%nids.size()]);

            return gmds::math::Triangle(pts[0], pts[1], pts[2]);
        }
    }

    std::cout<<"Face::getTriangle "<<std::endl;
    std::cout<<"looking for node AID "<<AID<<std::endl;
    for(int i=0; i<nids.size(); i++) {
        std::cout<<nids[i]<<" ";
    }
    std::cout<<std::endl;

    throw gmds::GMDSException("Face::getTriangle node not found");
}
/*------------------------------------------------------------------------*/
Region::Region()
 : Cell()
{
}
/*------------------------------------------------------------------------*/
Region::Region(const Mesh* AM, const TCellID AID)
 : Cell(AM, AID)
{
}
/*------------------------------------------------------------------------*/
ECellType
Region::computeType() const
{
        ECellType type;
        TInt32 nb_nodes = getNbNodes();
        if (nb_nodes == 4) {
                type = KMDS_TETRA;
        } else if (nb_nodes == 8) {
                type = KMDS_HEX;
        } else if (nb_nodes == 5) {
                type = KMDS_PYRAMID;
        } else if (nb_nodes == 6) {
                type = KMDS_PRISM3;
        } else {
                type = KMDS_UNDEF;
        }
        return type;
}
/*------------------------------------------------------------------------*/
void
Region::setNodes(const Kokkos::View<TCellID*>& AN)
{
    owner->m_R.set(id, AN);
}
/*------------------------------------------------------------------------*/
void
Region::setNodes(TCellID* AIds, int ASize)
{
        owner->m_R.setUnsafe(id, AIds, ASize);
}
/*------------------------------------------------------------------------*/
void
Region::setUnsafeNodes(const Kokkos::View<TCellID*>& AN)
{
        owner->m_R.setUnsafe(id, AN);
}
/*------------------------------------------------------------------------*/
void
Region::nodeIds(Kokkos::View<TCellID*>& AV) const
{
        owner->m_R.get(id, AV);
}
/*------------------------------------------------------------------------*/
void
Region::nodeIds(TCellID AIds[12], int* ASize) const
{
        owner->m_R.get(id, AIds, ASize);
}
/*------------------------------------------------------------------------*/
TSize
Region::getNbNodes() const
{
        return owner->m_R.getNbNodes(id);
}
/*------------------------------------------------------------------------*/
void
Region::nodes(Kokkos::View<Node*>& AV) const
{
        Kokkos::View<TCellID*> ids;
        nodeIds(ids);
        int nb_nodes = ids.size();
        for (auto i = 0; i < nb_nodes + 1; i++) {
                AV(i).owner = owner;
                AV(i).id = ids(i);
        }
}
/*------------------------------------------------------------------------*/
void
Region::nodes(Node ANodes[12], int* ASize) const
{
        TCellID ids[12];
        nodeIds(ids,ASize);
        int nb_nodes = *ASize;
        for (auto i = 0; i < nb_nodes + 1; i++) {

                ANodes[i].owner = owner;
                ANodes[i].id = ids[i];
        }
}
/*------------------------------------------------------------------------*/
Kokkos::View<Node*>
Region::nodes() const
{
        Kokkos::View<TCellID*> ids;
        nodeIds(ids);
        int nb_nodes = ids.size();
        Kokkos::View<Node*> AV("N", nb_nodes);
        for (auto i = 0; i < nb_nodes + 1; i++) {
                AV(i).owner = owner;
                AV(i).id = ids(i);
        }
        return AV;
}
/*----------------------------------------------------------------------------*/
TCoord
Region::surfvol() const
{
        TCoord vol=0.0;

        Node n[12];
        int nbNodes;
        this->nodes(n, &nbNodes);

        switch (this->computeType()) {
            case KMDS_TETRA:
                vol = gmds::math::Tetrahedron(n[0].getPoint(),
                                        n[1].getPoint(),
                                        n[2].getPoint(),
                                        n[3].getPoint()).getVolume();
                break;
            case KMDS_HEX:
                vol = gmds::math::Hexahedron(n[0].getPoint(),
                                             n[1].getPoint(),
                                             n[2].getPoint(),
                                             n[3].getPoint(),
                                             n[4].getPoint(),
                                             n[5].getPoint(),
                                             n[6].getPoint(),
                                             n[7].getPoint()).getVolume();
                break;
            case KMDS_PYRAMID:
                vol = gmds::math::Pyramid(n[0].getPoint(),
                                          n[1].getPoint(),
                                          n[2].getPoint(),
                                          n[3].getPoint(),
                                          n[4].getPoint()).getVolume();

            case KMDS_PRISM3:
                vol = gmds::math::Prism3(n[0].getPoint(),
                                         n[1].getPoint(),
                                         n[2].getPoint(),
                                         n[3].getPoint(),
                                         n[4].getPoint(),
                                         n[5].getPoint()).getVolume();
                break;
            default:
                throw KException("Region::volume can not be computed for this cell type.");
        }

        return vol;
}
/*----------------------------------------------------------------------------*/
TCoord Region::scaledJacobian() const
{
        TCoord vol=0.0;

        Node n[12];
        int nbNodes;
        this->nodes(n, &nbNodes);

        switch (this->computeType()) {
                case KMDS_TETRA:
                        vol = gmds::math::Tetrahedron(n[0].getPoint(),
                                                      n[1].getPoint(),
                                                      n[2].getPoint(),
                                                      n[3].getPoint()).computeScaledJacobian();
                break;
                case KMDS_HEX:
                        vol = gmds::math::Hexahedron(n[0].getPoint(),
                                                     n[1].getPoint(),
                                                     n[2].getPoint(),
                                                     n[3].getPoint(),
                                                     n[4].getPoint(),
                                                     n[5].getPoint(),
                                                     n[6].getPoint(),
                                                     n[7].getPoint()).computeScaledJacobian();
                break;
                case KMDS_PYRAMID:
                        vol = gmds::math::Pyramid(n[0].getPoint(),
                                                  n[1].getPoint(),
                                                  n[2].getPoint(),
                                                  n[3].getPoint(),
                                                  n[4].getPoint()).computeScaledJacobian();

                case KMDS_PRISM3:
                        vol = gmds::math::Prism3(n[0].getPoint(),
                                                 n[1].getPoint(),
                                                 n[2].getPoint(),
                                                 n[3].getPoint(),
                                                 n[4].getPoint(),
                                                 n[5].getPoint()).computeScaledJacobian();
                break;
                default:
                        throw KException("Region::scaledJacobian can not be computed for this cell type.");
        }

        return vol;
}
/*----------------------------------------------------------------------------*/
TCoord Region::scaledJacobianAt(kmds::TCellID ANId, gmds::math::Point APt) const
{
    TCoord qual=0.0;

    Node n[12];
    int nbNodes;
    this->nodes(n, &nbNodes);

    // identify the local index of the node
    int nid_l = -1;
    for(int i_n=0; i_n<nbNodes; i_n++) {
        if(ANId == n[i_n].id) {
            nid_l = i_n;
            break;
        }
    }

    switch (this->computeType()) {
        case KMDS_HEX:

            switch (nid_l) {
                case 0:
                    qual = gmds::math::Hexahedron(APt,
                                                  n[1].getPoint(),
                                                  n[2].getPoint(),
                                                  n[3].getPoint(),
                                                  n[4].getPoint(),
                                                  n[5].getPoint(),
                                                  n[6].getPoint(),
                                                  n[7].getPoint()).computeScaledJacobianAt(nid_l);
                    break;
                case 1:
                    qual = gmds::math::Hexahedron(n[0].getPoint(),
                                                  APt,
                                                  n[2].getPoint(),
                                                  n[3].getPoint(),
                                                  n[4].getPoint(),
                                                  n[5].getPoint(),
                                                  n[6].getPoint(),
                                                  n[7].getPoint()).computeScaledJacobianAt(nid_l);
                    break;
                case 2:
                    qual = gmds::math::Hexahedron(n[0].getPoint(),
                                                  n[1].getPoint(),
                                                  APt,
                                                  n[3].getPoint(),
                                                  n[4].getPoint(),
                                                  n[5].getPoint(),
                                                  n[6].getPoint(),
                                                  n[7].getPoint()).computeScaledJacobianAt(nid_l);
                    break;
                case 3:
                    qual = gmds::math::Hexahedron(n[0].getPoint(),
                                                  n[1].getPoint(),
                                                  n[2].getPoint(),
                                                  APt,
                                                  n[4].getPoint(),
                                                  n[5].getPoint(),
                                                  n[6].getPoint(),
                                                  n[7].getPoint()).computeScaledJacobianAt(nid_l);
                    break;
                case 4:
                    qual = gmds::math::Hexahedron(n[0].getPoint(),
                                                  n[1].getPoint(),
                                                  n[2].getPoint(),
                                                  n[3].getPoint(),
                                                  APt,
                                                  n[5].getPoint(),
                                                  n[6].getPoint(),
                                                  n[7].getPoint()).computeScaledJacobianAt(nid_l);
                    break;
                case 5:
                    qual = gmds::math::Hexahedron(n[0].getPoint(),
                                                  n[1].getPoint(),
                                                  n[2].getPoint(),
                                                  n[3].getPoint(),
                                                  n[4].getPoint(),
                                                  APt,
                                                  n[6].getPoint(),
                                                  n[7].getPoint()).computeScaledJacobianAt(nid_l);
                    break;
                case 6:
                    qual = gmds::math::Hexahedron(n[0].getPoint(),
                                                  n[1].getPoint(),
                                                  n[2].getPoint(),
                                                  n[3].getPoint(),
                                                  n[4].getPoint(),
                                                  n[5].getPoint(),
                                                  APt,
                                                  n[7].getPoint()).computeScaledJacobianAt(nid_l);
                    break;
                case 7:
                    qual = gmds::math::Hexahedron(n[0].getPoint(),
                                                  n[1].getPoint(),
                                                  n[2].getPoint(),
                                                  n[3].getPoint(),
                                                  n[4].getPoint(),
                                                  n[5].getPoint(),
                                                  n[6].getPoint(),
                                                  APt).computeScaledJacobianAt(nid_l);
                    break;

            }
            break;
        default:
            throw KException("Region::scaledJacobianAt can not be computed for this cell type.");
    }

    return qual;
}
/*----------------------------------------------------------------------------*/
gmds::math::Point
Region::midpoint() const
{

        gmds::math::Point pt = 0.0;

        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);
        const int nbNodes = nodes.size();

        for(int i=0; i<nbNodes; i++) {

                pt = pt + owner->getNodeLocation(nodes[i]);
        }

        pt = 1./nbNodes * pt;

        return pt;
}
/*----------------------------------------------------------------------------*/
std::vector<std::vector<TCellID> >
Region::getOrderedNodesFacesIDs() const
{
        std::vector<std::vector<TCellID> > orderedNodesFaces;

        //std::vector<TCellID> nodes = this->getIDs<Node>();
        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);

        switch(this->computeType()) {
                case KMDS_HEX:
                        orderedNodesFaces.resize(6);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[7];
                orderedNodesFaces[1][1] = nodes[4];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[6];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[5];
                orderedNodesFaces[2][1] = nodes[4];
                orderedNodesFaces[2][2] = nodes[0];
                orderedNodesFaces[2][3] = nodes[1];
                orderedNodesFaces[3].resize(4);
                orderedNodesFaces[3][0] = nodes[7];
                orderedNodesFaces[3][1] = nodes[6];
                orderedNodesFaces[3][2] = nodes[2];
                orderedNodesFaces[3][3] = nodes[3];
                orderedNodesFaces[4].resize(4);
                orderedNodesFaces[4][0] = nodes[4];
                orderedNodesFaces[4][1] = nodes[7];
                orderedNodesFaces[4][2] = nodes[3];
                orderedNodesFaces[4][3] = nodes[0];
                orderedNodesFaces[5].resize(4);
                orderedNodesFaces[5][0] = nodes[6];
                orderedNodesFaces[5][1] = nodes[5];
                orderedNodesFaces[5][2] = nodes[1];
                orderedNodesFaces[5][3] = nodes[2];
                break;
                case KMDS_TETRA:
                        orderedNodesFaces.resize(4);
                orderedNodesFaces[0].resize(3);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[2];
                orderedNodesFaces[0][2] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[3];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[0];
                orderedNodesFaces[3][2] = nodes[3];
                break;
                case KMDS_PYRAMID:
                        orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[4];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[4];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[3];
                orderedNodesFaces[3][2] = nodes[4];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[0];
                orderedNodesFaces[4][2] = nodes[4];
                break;
                case KMDS_PRISM3:
                        orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[1];
                orderedNodesFaces[0][2] = nodes[4];
                orderedNodesFaces[0][3] = nodes[3];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[1];
                orderedNodesFaces[1][1] = nodes[2];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[4];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[2];
                orderedNodesFaces[2][1] = nodes[0];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[2][3] = nodes[5];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[0];
                orderedNodesFaces[3][1] = nodes[2];
                orderedNodesFaces[3][2] = nodes[1];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[4];
                orderedNodesFaces[4][2] = nodes[5];
                break;
                default:
                        throw kmds::KException("Region::getOrderedNodesFacesIDs not implemented for this region type");
        }

        return orderedNodesFaces;
}
/*----------------------------------------------------------------------------*/
std::vector<std::vector<TCellID> >
Region::getIJKNodesFacesIDs() const
{
        std::vector<std::vector<TCellID> > orderedNodesFaces;

        //std::vector<TCellID> nodes = this->getIDs<Node>();
        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);

        switch(this->computeType()) {
                case KMDS_HEX:
                        orderedNodesFaces.resize(6);
                // bottom
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[1];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[3];
                // top
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[4];
                orderedNodesFaces[1][1] = nodes[5];
                orderedNodesFaces[1][2] = nodes[6];
                orderedNodesFaces[1][3] = nodes[7];
                // left
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[0];
                orderedNodesFaces[2][1] = nodes[3];
                orderedNodesFaces[2][2] = nodes[7];
                orderedNodesFaces[2][3] = nodes[4];
                // right
                orderedNodesFaces[3].resize(4);
                orderedNodesFaces[3][0] = nodes[1];
                orderedNodesFaces[3][1] = nodes[2];
                orderedNodesFaces[3][2] = nodes[6];
                orderedNodesFaces[3][3] = nodes[5];
                // front
                orderedNodesFaces[4].resize(4);
                orderedNodesFaces[4][0] = nodes[0];
                orderedNodesFaces[4][1] = nodes[1];
                orderedNodesFaces[4][2] = nodes[5];
                orderedNodesFaces[4][3] = nodes[4];
                // back
                orderedNodesFaces[5].resize(4);
                orderedNodesFaces[5][0] = nodes[3];
                orderedNodesFaces[5][1] = nodes[2];
                orderedNodesFaces[5][2] = nodes[6];
                orderedNodesFaces[5][3] = nodes[7];
                break;
                default:
                        throw kmds::KException("Region::getIJKNodesFacesIDs not implemented for this region type");
        }

        return orderedNodesFaces;
}
/*----------------------------------------------------------------------------*/
bool
Region::isFaceOrientedOutward(Kokkos::View<kmds::TCellID *> &AIDs) const
{
        std::vector<std::vector<TCellID> > orderedNodesFaces = this->getOrderedNodesFacesIDs();

        // first find the face
        for(int iFace=0; iFace<orderedNodesFaces.size(); iFace++) {

                if(AIDs.size() == orderedNodesFaces[iFace].size()) {

                        int nbNodesMatched = 0;

                        for(int iNode1=0; iNode1<orderedNodesFaces[iFace].size(); iNode1++) {
                                for(int iNode2=0; iNode2<AIDs.size(); iNode2++) {
                                        if(AIDs[iNode2] == orderedNodesFaces[iFace][iNode1]) {
                                                nbNodesMatched++;
                                        }
                                }
                        }
                        if(nbNodesMatched == AIDs.size()) {
                                // face is found, now check the orientation
                                if(orderedNodesFaces[iFace].size() < 3) {
                                        throw KException("Region::isFaceOrientedOutward face with less than 3 nodes.");
                                }
                                TCellID firstNode = orderedNodesFaces[iFace][0];
                                TCellID secondNode = orderedNodesFaces[iFace][1];

                                for(int iNode2=0; iNode2<AIDs.size(); iNode2++) {
                                        if(AIDs[iNode2] == firstNode) {
                                                if(AIDs[(iNode2+1)%AIDs.size()] == secondNode) {
                                                        return true;
                                                } else  if (AIDs[(AIDs.size()+iNode2-1)%AIDs.size()] == secondNode) {
                                                        return false;
                                                } else {
                                                        throw KException("Region::isFaceOrientedOutward jumbled face.");
                                                }
                                        }
                                }
                        }
                }
        }

        throw KException("Region::isFaceOrientedOutward face not found.");
}
/*----------------------------------------------------------------------------*/
std::vector<FakeFace>
Region::getFakeFaces() const
{
        const std::vector<std::vector<TCellID> > orderedNodesFacesIDs = this->getOrderedNodesFacesIDs();
        std::vector<FakeFace> fakeFaces;

        for(size_t iFace=0; iFace<orderedNodesFacesIDs.size(); iFace++) {

                std::vector<TCellID> ids(orderedNodesFacesIDs[iFace].size());
                for(size_t iNode=0; iNode<ids.size(); iNode++) {
                        ids[iNode] = orderedNodesFacesIDs[iFace][iNode];
                }

                fakeFaces.push_back(FakeFace(ids));
        }

        return fakeFaces;
}
/*----------------------------------------------------------------------------*/
std::vector<bool>
Region::getFakeFacesOrientation() const
{
        const std::vector<std::vector<TCellID> > orderedNodesFacesIDs = this->getOrderedNodesFacesIDs();

        std::vector<bool> isOutward(orderedNodesFacesIDs.size());

        for(size_t iFace=0; iFace<orderedNodesFacesIDs.size(); iFace++) {

                std::vector<TCellID> ids(orderedNodesFacesIDs[iFace].size());
                for(size_t iNode=0; iNode<ids.size(); iNode++) {
                        ids[iNode] = orderedNodesFacesIDs[iFace][iNode];
                }

                FakeFace ff(ids);

                std::vector<TCellID> ffids = ff.node_ids();

                for(size_t iNode=0; iNode<ids.size(); iNode++) {
                        if(ids[iNode] == ffids[0]) {

                                if(ids[(iNode + 1) % ids.size()] == ffids[1]) {
                                        isOutward[iFace] = true;
                                } else {
                                        isOutward[iFace] = false;
                                }

                                break;
                        }
                }
        }



        return isOutward;
}
/*----------------------------------------------------------------------------*/
std::vector<FakeEdge>
Region::getFakeEdges() const
{
        std::vector<FakeEdge> fakeEdges;

        Kokkos::View<TCellID*> nodes;
        this->nodeIds(nodes);

        switch(this->computeType()) {
                case KMDS_HEX:
                        fakeEdges.resize(12);
                fakeEdges[ 0] = FakeEdge(nodes[0],nodes[1]);
                fakeEdges[ 1] = FakeEdge(nodes[1],nodes[2]);
                fakeEdges[ 2] = FakeEdge(nodes[2],nodes[3]);
                fakeEdges[ 3] = FakeEdge(nodes[3],nodes[0]);
                fakeEdges[ 4] = FakeEdge(nodes[4],nodes[5]);
                fakeEdges[ 5] = FakeEdge(nodes[5],nodes[6]);
                fakeEdges[ 6] = FakeEdge(nodes[6],nodes[7]);
                fakeEdges[ 7] = FakeEdge(nodes[7],nodes[4]);
                fakeEdges[ 8] = FakeEdge(nodes[0],nodes[4]);
                fakeEdges[ 9] = FakeEdge(nodes[1],nodes[5]);
                fakeEdges[10] = FakeEdge(nodes[2],nodes[6]);
                fakeEdges[11] = FakeEdge(nodes[3],nodes[7]);
                break;
                case KMDS_TETRA:
                        fakeEdges.resize(6);
                fakeEdges[ 0] = FakeEdge(nodes[0],nodes[1]);
                fakeEdges[ 1] = FakeEdge(nodes[1],nodes[2]);
                fakeEdges[ 2] = FakeEdge(nodes[2],nodes[0]);
                fakeEdges[ 3] = FakeEdge(nodes[0],nodes[3]);
                fakeEdges[ 4] = FakeEdge(nodes[1],nodes[3]);
                fakeEdges[ 5] = FakeEdge(nodes[2],nodes[3]);
                break;
                case KMDS_PYRAMID:
                        fakeEdges.resize(8);
                fakeEdges[ 0] = FakeEdge(nodes[0],nodes[1]);
                fakeEdges[ 1] = FakeEdge(nodes[1],nodes[2]);
                fakeEdges[ 2] = FakeEdge(nodes[2],nodes[3]);
                fakeEdges[ 3] = FakeEdge(nodes[3],nodes[0]);
                fakeEdges[ 4] = FakeEdge(nodes[0],nodes[4]);
                fakeEdges[ 5] = FakeEdge(nodes[1],nodes[4]);
                fakeEdges[ 6] = FakeEdge(nodes[2],nodes[4]);
                fakeEdges[ 7] = FakeEdge(nodes[3],nodes[4]);
                break;
                case KMDS_PRISM3:
                        fakeEdges.resize(9);
                fakeEdges[ 0] = FakeEdge(nodes[0],nodes[1]);
                fakeEdges[ 1] = FakeEdge(nodes[1],nodes[2]);
                fakeEdges[ 2] = FakeEdge(nodes[2],nodes[0]);
                fakeEdges[ 3] = FakeEdge(nodes[3],nodes[4]);
                fakeEdges[ 4] = FakeEdge(nodes[4],nodes[5]);
                fakeEdges[ 5] = FakeEdge(nodes[5],nodes[3]);
                fakeEdges[ 6] = FakeEdge(nodes[0],nodes[3]);
                fakeEdges[ 7] = FakeEdge(nodes[1],nodes[4]);
                fakeEdges[ 8] = FakeEdge(nodes[2],nodes[5]);
                break;
                default:
                        throw kmds::KException("Region::getFakeEdges not implemented for this region type");

        }

        return fakeEdges;
}

/*----------------------------------------------------------------------------*/
void
Region::getEdgeInfo(Kokkos::View<TCellID *> &AIDs, int &numFacet, bool &isIJK) const
{
    Kokkos::View<TCellID *> nodes;
    this->nodeIds(nodes);

    if((AIDs[0] == nodes[0]) && (AIDs[1] == nodes[1])) {
        numFacet = 0;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[0]) && (AIDs[0] == nodes[1])) {
        numFacet = 0;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[1]) && (AIDs[1] == nodes[2])) {
        numFacet = 1;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[1]) && (AIDs[0] == nodes[2])) {
        numFacet = 1;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[2]) && (AIDs[1] == nodes[3])) {
        numFacet = 2;
        isIJK = false;
        return;
    }
    if((AIDs[1] == nodes[2]) && (AIDs[0] == nodes[3])) {
        numFacet = 2;
        isIJK = true;
        return;
    }
    if((AIDs[0] == nodes[3]) && (AIDs[1] == nodes[0])) {
        numFacet = 3;
        isIJK = false;
        return;
    }
    if((AIDs[1] == nodes[3]) && (AIDs[0] == nodes[0])) {
        numFacet = 3;
        isIJK = true;
        return;
    }
    if((AIDs[0] == nodes[0]) && (AIDs[1] == nodes[4])) {
        numFacet = 4;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[0]) && (AIDs[0] == nodes[4])) {
        numFacet = 4;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[1]) && (AIDs[1] == nodes[5])) {
        numFacet = 5;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[1]) && (AIDs[0] == nodes[5])) {
        numFacet = 5;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[2]) && (AIDs[1] == nodes[6])) {
        numFacet = 6;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[2]) && (AIDs[0] == nodes[6])) {
        numFacet = 6;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[3]) && (AIDs[1] == nodes[7])) {
        numFacet = 7;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[3]) && (AIDs[0] == nodes[7])) {
        numFacet = 7;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[4]) && (AIDs[1] == nodes[5])) {
        numFacet = 8;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[4]) && (AIDs[0] == nodes[5])) {
        numFacet = 8;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[5]) && (AIDs[1] == nodes[6])) {
        numFacet = 9;
        isIJK = true;
        return;
    }
    if((AIDs[1] == nodes[5]) && (AIDs[0] == nodes[6])) {
        numFacet = 9;
        isIJK = false;
        return;
    }
    if((AIDs[0] == nodes[6]) && (AIDs[1] == nodes[7])) {
        numFacet = 10;
        isIJK = false;
        return;
    }
    if((AIDs[1] == nodes[6]) && (AIDs[0] == nodes[7])) {
        numFacet = 10;
        isIJK = true;
        return;
    }
    if((AIDs[0] == nodes[7]) && (AIDs[1] == nodes[4])) {
        numFacet = 11;
        isIJK = false;
        return;
    }
    if((AIDs[1] == nodes[7]) && (AIDs[0] == nodes[4])) {
        numFacet = 11;
        isIJK = true;
        return;
    }
}

/*----------------------------------------------------------------------------*/
void
Region::getFaceInfo(const kmds::TCellID AID1,
                    const kmds::TCellID AID2,
                    const kmds::TCellID AID3,
                    const kmds::TCellID AID4,
                    int &ANumFacet,
                    int &AOffset,
                    bool &AIsIJK) const
{
    Kokkos::View<TCellID *> nids;
    this->nodeIds(nids);

    // find which facet it is

    if ((AID1 == nids[0] || AID1 == nids[1] || AID1 == nids[2] || AID1 == nids[3])
        && (AID2 == nids[0] || AID2 == nids[1] || AID2 == nids[2] || AID2 == nids[3])
        && (AID3 == nids[0] || AID3 == nids[1] || AID3 == nids[2] || AID3 == nids[3])
        && (AID4 == nids[0] || AID4 == nids[1] || AID4 == nids[2] || AID4 == nids[3])) {
        // bottom
        ANumFacet = 0;

    } else if ((AID1 == nids[4] || AID1 == nids[5] || AID1 == nids[6] || AID1 == nids[7])
               && (AID2 == nids[4] || AID2 == nids[5] || AID2 == nids[6] || AID2 == nids[7])
               && (AID3 == nids[4] || AID3 == nids[5] || AID3 == nids[6] || AID3 == nids[7])
               && (AID4 == nids[4] || AID4 == nids[5] || AID4 == nids[6] || AID4 == nids[7])) {
        // top
        ANumFacet = 1;

    } else if ((AID1 == nids[0] || AID1 == nids[3] || AID1 == nids[7] || AID1 == nids[4])
               && (AID2 == nids[0] || AID2 == nids[3] || AID2 == nids[7] || AID2 == nids[4])
               && (AID3 == nids[0] || AID3 == nids[3] || AID3 == nids[7] || AID3 == nids[4])
               && (AID4 == nids[0] || AID4 == nids[3] || AID4 == nids[7] || AID4 == nids[4])) {
        // left
        ANumFacet = 2;

    } else if ((AID1 == nids[1] || AID1 == nids[2] || AID1 == nids[6] || AID1 == nids[5])
               && (AID2 == nids[1] || AID2 == nids[2] || AID2 == nids[6] || AID2 == nids[5])
               && (AID3 == nids[1] || AID3 == nids[2] || AID3 == nids[6] || AID3 == nids[5])
               && (AID4 == nids[1] || AID4 == nids[2] || AID4 == nids[6] || AID4 == nids[5])) {
        // right
        ANumFacet = 3;

    } else if ((AID1 == nids[0] || AID1 == nids[1] || AID1 == nids[5] || AID1 == nids[4])
               && (AID2 == nids[0] || AID2 == nids[1] || AID2 == nids[5] || AID2 == nids[4])
               && (AID3 == nids[0] || AID3 == nids[1] || AID3 == nids[5] || AID3 == nids[4])
               && (AID4 == nids[0] || AID4 == nids[1] || AID4 == nids[5] || AID4 == nids[4])) {
        // front
        ANumFacet = 4;

    } else if ((AID1 == nids[3] || AID1 == nids[2] || AID1 == nids[6] || AID1 == nids[7])
               && (AID2 == nids[3] || AID2 == nids[2] || AID2 == nids[6] || AID2 == nids[7])
               && (AID3 == nids[3] || AID3 == nids[2] || AID3 == nids[6] || AID3 == nids[7])
               && (AID4 == nids[3] || AID4 == nids[2] || AID4 == nids[6] || AID4 == nids[7])) {
        // back
        ANumFacet = 5;
    }

    // find the relation between the face and its corresponding region facet
    const std::vector<std::vector<kmds::TCellID> > ijkNodes = this->getIJKNodesFacesIDs();

    for (int i=0; i<ijkNodes[ANumFacet].size(); i++) {
        if(ijkNodes[ANumFacet][i] == AID1) {
            AOffset = i;

            if(ijkNodes[ANumFacet][(i + 1) % ijkNodes[ANumFacet].size()] != AID2) {
                AIsIJK = false;
                break;
            } else {
                AIsIJK = true;
                break;
            }
        }
    }

}
/*----------------------------------------------------------------------------*/
void
Region::computeBoundingBox(kmds::TCoord minXYZ[3], kmds::TCoord maxXYZ[3]) const
{
    Kokkos::View<kmds::TCellID*> nids;
    this->nodeIds(nids);

    minXYZ[0] = HUGE_VALF;
    minXYZ[1] = HUGE_VALF;
    minXYZ[2] = HUGE_VALF;

    maxXYZ[0] = -HUGE_VALF;
    maxXYZ[1] = -HUGE_VALF;
    maxXYZ[2] = -HUGE_VALF;

    for(int i=0; i<nids.size(); i++) {
        TCoord xyz[3];

        this->owner->getNodeLocation(nids[i], xyz[0], xyz[1], xyz[2]);

        minXYZ[0] = std::min(minXYZ[0], xyz[0]);
        minXYZ[1] = std::min(minXYZ[1], xyz[1]);
        minXYZ[2] = std::min(minXYZ[2], xyz[2]);

        maxXYZ[0] = std::max(maxXYZ[0], xyz[0]);
        maxXYZ[1] = std::max(maxXYZ[1], xyz[1]);
        maxXYZ[2] = std::max(maxXYZ[2], xyz[2]);
    }
}
/*------------------------------------------------------------------------*/
