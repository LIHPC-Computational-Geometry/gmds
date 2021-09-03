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
/** \file    CellHandle.h
 *  \author  F. LEDOUX
 *  \date    03/07/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_CELL_HANDLE_H_
#define KMDS_CELL_HANDLE_H_
/*----------------------------------------------------------------------------*/
// STL headers
#include <vector>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>

#include <gmds/math/Point.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <KM/Utils/FakeTypes.h>
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
class Mesh;
struct Cell
{
        /*------------------------------------------------------------------------*/
        /** \brief  Default constructor. Useful for containers initialization
         */
        Cell();

        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        Cell(const Mesh* AM, const TCellID AID);

        /** mesh the cell belongs to */
        Mesh* owner;
        /** cell id in the owner mesh*/
        TCellID id;
};

/*----------------------------------------------------------------------------*/
/* \struct Node
 *
 * \brief Lightweight structure to access nodes through an object-programming
 * style
 */
/*----------------------------------------------------------------------------*/
struct Face;
struct Region;
struct Node : public Cell
{
        /*------------------------------------------------------------------------*/
        /** \brief  Default constructor. Useful for containers initialization
         */
        Node();

        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        Node(const Mesh* AM, const TCellID AID);

        void setLocation(const TCoord AX, const TCoord AY, const TCoord AZ);
        void getLocation(TCoord& AX, TCoord& AY, TCoord& AZ) const;

        void setEdges(const Kokkos::View<TCellID*>& AV);

        void faceIds(Kokkos::View<TCellID*>& AV) const;
        void setFaces(const Kokkos::View<TCellID*>& AV);
        TSize getNbFaces() const;
        Kokkos::View<Node*> faces() const;
        void faces(Kokkos::View<Face*>& AV) const;

        void regionIds(Kokkos::View<TCellID*>& AV) const;
        void setRegions(const Kokkos::View<TCellID*>& AV);
        void setRegions(const TCellID* AV, const int ASize);

        TSize getNbRegions() const;
        Kokkos::View<Node*> regions() const;
        void regions(Kokkos::View<Region*>& AV) const;

        ECellType computeType() const;

        gmds::math::Point getPoint() const;
};
/*----------------------------------------------------------------------------*/
/* \struct Face
 *
 * \brief Lightweight structure to access faces through an object-programming
 * style
 */
/*----------------------------------------------------------------------------*/
struct Edge : public Cell
{
        /*------------------------------------------------------------------------*/
        /** \brief  Default constructor. Useful for containers initialization
         */
        Edge();
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        Edge(const Mesh* AM, const TCellID AID);

        void nodeIds(Kokkos::View<TCellID*>& AV) const;
        void nodeIds(TCellID AIds[2]) const;
        TSize getNbNodes() const;


        void nodes(Kokkos::View<Node*>& AV) const;
        void nodes(Node ANodes[2]) const;


//        void setNodes(TCellID AIds[2]);
        void setNodes(TCellID AN0, TCellID AN1);


    TCoord surfvol() const;


        ECellType computeType() const;
};
/*----------------------------------------------------------------------------*/
/* \struct Face
 *
 * \brief Lightweight structure to access faces through an object-programming
 * style
 */
/*----------------------------------------------------------------------------*/
struct Face : public Cell
{
        /*------------------------------------------------------------------------*/
        /** \brief  Default constructor. Useful for containers initialization
         */
        Face();
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        Face(const Mesh* AM, const TCellID AID);

        void nodeIds(Kokkos::View<TCellID*>& AV) const;
        void nodeIds(TCellID AIds[5], int* ASize) const;
        TSize getNbNodes() const;

        Kokkos::View<Node*> nodes() const;

        void nodes(Kokkos::View<Node*>& AV) const;
        void nodes(Node ANodes[5], int* ASize) const;
        void setNodes(const Kokkos::View<TCellID*>& AN);
        void setNodes(TCellID* AIds, int ASize);
        void setUnsafeNodes(const Kokkos::View<TCellID*>& AN);


    /*------------------------------------------------------------------------*/
    /** \brief  Compute the surface of the face. Only works for triangles and
     *         quads. The computation is simply based on the
     *         decomposition of the face into triangle elements.
     *
     * \return the surface of the face
     */
    TCoord surfvol() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the scaled jacobian of the face.
     *
     * \return the scaled jacobian of the face
     */
    TCoord scaledJacobian() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the scaled jacobian of the face at one node.
     *
     * \return the scaled jacobian of the face
     */
    TCoord scaledJacobianAt(kmds::TCellID ANId, gmds::math::Point APt) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the midpoint of the region, computed as the barycenter of the
         *          nodes of the cell.
         *
         * \return the center of the cell
         */
        gmds::math::Point midpoint() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return whether the edge is outward oriented.
         *
         * \param the ids of the nodes of the edge
         *
         * \return whether the edge is outward oriented
         */
        bool isEdgeOrientedOutward(Kokkos::View<TCellID*>& AV) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the fake edges of the face.
         *
         * \return the fake edges
         */
        std::vector<FakeEdge> getFakeEdges() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return information on a facet (edge) of the face.
     *
     * \param the ids of the nodes of the edge
     * \param the ids of the nodes of the edge
     * \param the ids of the nodes of the edge
     * \param the ids of the nodes of the edge
     *
     */
    void getFacetInfo(Kokkos::View<TCellID*>& AIDs, int& numFacet, bool& isOutward, int& offset) const;


    /*------------------------------------------------------------------------*/
    /** \brief  Return information on an edge of the region.
     *
     * \param AIDs the ids of the nodes of the edge
     *
     */
    void getEdgeInfo(Kokkos::View<TCellID*>& AIDs, int& numFacet, bool& isIJK) const;


    void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const;


    gmds::math::Triangle getTriangle(kmds::TCellID AID)const ;


        ECellType computeType() const;
};
/*----------------------------------------------------------------------------*/
/* \struct Region
 *
 * \brief Lightweight structure to access regions through an object-programming
 * style
 */
/*----------------------------------------------------------------------------*/
struct Region : public Cell
{
        /*------------------------------------------------------------------------*/
        /** \brief  Default constructor. Useful for containers initialization
         */
        Region();
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        Region(const Mesh* AM, const TCellID AID);

        void nodeIds(Kokkos::View<TCellID*>& AV) const;
        void nodeIds(TCellID AIds[12], int* ASize) const;
        TSize getNbNodes() const;

        Kokkos::View<Node*> nodes() const;
        void nodes(Node ANodes[12], int* ASize) const;

        void nodes(Kokkos::View<Node*>& AV) const;
        void setNodes(const Kokkos::View<TCellID*>& AN);
        void setNodes(TCellID* AIds, int ASize);
        void setUnsafeNodes(const Kokkos::View<TCellID*>& AN);

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the volume of the region. Only works for tet, hex,
     *         prism3 and pyramid. The computation is simply based on the
     *         decomposition of the region into tetrahedral elements.
     *
     * \return the volume of the region
     */
    TCoord surfvol() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the scaled jacobian of the region.
     *
     * \return the scaled jacobian of the region
     */
    TCoord scaledJacobian() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the scaled jacobian of the region at one node.
     *
     * \return the scaled jacobian of the region
     */
    TCoord scaledJacobianAt(kmds::TCellID ANId, gmds::math::Point APt) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the midpoint of the region, computed as the barycenter of the
     *          nodes of the cell.
     *
     * \return the center of the cell
     */
    gmds::math::Point midpoint() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return the faces composing the region in the form of
     *      vectors of ordered nodes oriented outwards
     *
     * \return the ordered nodes ids of the faces
     */
    std::vector<std::vector<TCellID> > getOrderedNodesFacesIDs() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return the faces composing the region in the form of
     *      vectors of ordered nodes IJK-oriented
     *
     * \return the IJK-oriented nodes ids of the faces
     */
    std::vector<std::vector<TCellID> > getIJKNodesFacesIDs() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return whether the face is outward oriented.
     *
     * \param the ids of the nodes of the face
     *
     * \return whether the face is outward oriented
     */
    bool isFaceOrientedOutward(Kokkos::View<TCellID*>& AV) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return the fake faces of the region.
     *
     * \return the fake faces
     */
    std::vector<FakeFace> getFakeFaces() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return the fake faces orientation of the region.
     *
     * \return the fake faces orientation compared to outward
     */
    std::vector<bool> getFakeFacesOrientation() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return the fake edges of the region.
     *
     * \return the fake edges
     */
    std::vector<FakeEdge> getFakeEdges() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return information on an edge of the region.
     *
     * \param AIDs the ids of the nodes of the edge
     *
     */
    void getEdgeInfo(Kokkos::View<TCellID*>& AIDs, int& numFacet, bool& isIJK) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Return information on a face (quad) of the region (hex).
     *
     * \param AIDs the ids of the nodes of the face
     *
     */
    void getFaceInfo(const kmds::TCellID AID1,
                     const kmds::TCellID AID2,
                     const kmds::TCellID AID3,
                     const kmds::TCellID AID4,
                     int &numFacet,
                     int &offset,
                     bool &isIJK) const;

    void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const;


        ECellType computeType() const;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* KMDS_CELL_HANDLE_H_ */
/*----------------------------------------------------------------------------*/
