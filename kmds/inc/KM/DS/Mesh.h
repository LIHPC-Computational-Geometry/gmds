/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Mesh.h
 *  \author  F. LEDOUX
 *  \date    03/07/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_MESH_H_
#define KMDS_MESH_H_
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>

/*----------------------------------------------------------------------------*/
// STL headers
#include <cstdio>
#include <map>
#include <string>
#include <typeinfo>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <KM/DS/CellHandle.h>
#include <KM/DS/Connectivity.h>
#include <KM/DS/EContainer.h>
#include <KM/DS/FContainer.h>
#include <KM/DS/NContainer.h>
#include <KM/DS/RContainer.h>
#include <KM/Utils/FakeTypes.h>
#include <KM/Utils/GrowingView.h>
#include <KM/Utils/Variable.h>
#include <KM/Utils/VariableManager.h>
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
/* \class Mesh
 *
 * \brief Mesh data structure
 */
/*----------------------------------------------------------------------------*/
class Mesh
{
 public:
        friend struct Node;
        friend struct Edge;
        friend struct Face;
        friend struct Region;
        /*------------------------------------------------------------------------*/
        /** \brief  Default Constructor
         */
        Mesh();

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~Mesh();

        /*------------------------------------------------------------------------*/
        /** \brief  Resize the node/face container and relative variables
         *
         * \param[in] ASize the new capacity
         */
        void updateNodeCapacity(const TSize ASize);
        void updateEdgeCapacity(const TSize ASize);
        void updateFaceCapacity(const TSize ASize);
        void updateRegionCapacity(const TSize ASize);
        /*------------------------------------------------------------------------*/
        /** \brief  Double the size of node/face container and relative variables
         */
        void doubleNodeCapacity();
        void doubleEdgeCapacity();
        void doubleFaceCapacity();
        void doubleRegionCapacity();

    /*------------------------------------------------------------------------*/
    /** \brief  Delete all mesh entities
     */
    void removeAllNodes();
    void removeAllEdges();
    void removeAllFaces();
    void removeAllRegions();

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate \p ANb new nodes.
         *
         * \param[in] ANb the number of nodes to be created
         * \return    The index of the first created node. Others are contiguous right behind.
         */
        TCellID addNodes(const int ANb);
        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate one new node.
         *
         * \return    The index of the created node
         */
        TCellID addNode();
        /*------------------------------------------------------------------------*/
        /** \brief  Create a node from point (\p AX, \p AY, \p AZ).
         *
         * \param[in] AX X coordinate
         * \param[in] AY Y coordinate
         * \param[in] AZ Z coordinate
         * \return    The index of the created node
         */
        TCellID newNode(const TCoord AX, const TCoord AY, const TCoord AZ);
        /*------------------------------------------------------------------------*/
        /** \brief  Create a node from point (\p AX, \p AY, \p AZ).
         *
         * \param[in] AX X coordinate
         * \param[in] AY Y coordinate
         * \param[in] AZ Z coordinate
         * \return    The index of the created node
         */
        TCellID newNode(const gmds::math::Point APt);
        /*------------------------------------------------------------------------*/
        /** \brief  Remove/Desallocate the node of id \p AId
         *
         * \param[in] The index of the node to be removed.
         */
        void removeNode(const TCellID AId);
        /*------------------------------------------------------------------------*/
        /** \brief  Get the node \p AId
         *
         * \return  the node of id \p AId
         */
        Node getNode(const TCellID AId) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Check if the node \p AId exists in the mesh
         *
         * \return  true if it exists, false otherwise
         */
        bool hasNode(const TCellID AId) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the capacity in terms of nodes
         */
        TSize getNodeCapacity() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns more (superior or equal) than the max id of regions
         */
        TSize getNodeSupID() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of nodes
         */
        TSize getNbNodes() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Get the coordinates of node \p AId.
         *
         * \param[in] AId The index of the node
         * \param[out] AX  X coordinate
         * \param[out] AY  Y coordinate
         * \param[out] AZ  Z coordinate
         */
        void getNodeLocation(const TCellID AId, TCoord& AX, TCoord& AY, TCoord& AZ) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Get the coordinates of node \p AId.
     *
     * \param[in] AId The index of the node
     *
     * \return the node location (as a point)
     */
    gmds::math::Point getNodeLocation(const TCellID AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Set the coordinates of node \p AId.
         *
         * \param[in] AId The index of the node
         * \param[in] AX  X coordinate
         * \param[in] AY  Y coordinate
         * \param[in] AZ  Z coordinate
         */
        void setNodeLocation(const TCellID AId, const TCoord AX, const TCoord AY, const TCoord AZ);

    /*------------------------------------------------------------------------*/
    /** \brief  Set the coordinates of node \p AId.
     *
     * \param[in] AId The index of the node
     * \param[in] APt The new location (as a point)
     *
     */
    void setNodeLocation(const TCellID AId, const gmds::math::Point APt);

    /*------------------------------------------------------------------------*/
    /** \brief  Add/allocate \p ANb new edges.
     *
     * \param[in] ANb the number of edges to be created
     * \return    The index of the first created edge. Others are contiguous right behind.
     */
    TCellID addEdges(const int ANb);
    /*------------------------------------------------------------------------*/
    /** \brief  Add/allocate one new edge.
     *
     * \return    The index of the created edge
     */
    TCellID addEdge();
    /*------------------------------------------------------------------------*/
    /** \brief  Create an edge from nodes (\p AN1, \p AN2).
     *
     * \param[in] AN1 first node id
     * \param[in] AN2 second node id
     * \return    The index of the created edge
     */
    TCellID newEdge(const TCellID AN1, const TCellID AN2);
    TCellID newEdge_unsafe(const TCellID AN1, const TCellID AN2);
    /*------------------------------------------------------------------------*/
    /** \brief  Remove/Desallocate the edge of id \p AId
     *
     * \param[in] The index of the edge to be removed.
     */
    void removeEdge(const TCellID AId);
    /*------------------------------------------------------------------------*/
    /** \brief  Get the edge \p AId
     *
     * \return  the edge of id \p AId
     */
    Edge getEdge(const TCellID AId) const;
    /*------------------------------------------------------------------------*/
    /** \brief  Check if the edge \p AId exists in the mesh
     *
     * \return  true if it exists, false otherwise
     */
    bool hasEdge(const TCellID AId) const;
    /*------------------------------------------------------------------------*/
    /** \brief  Returns the capacity in terms of edges
     */
    TSize getEdgeCapacity() const;
    /*------------------------------------------------------------------------*/
    /** \brief  Returns more (superior or equal) than the max id of regions
     */
    TSize getEdgeSupID() const;
    /*------------------------------------------------------------------------*/
    /** \brief  Returns the number of edges
     */
    TSize getNbEdges() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate \p ANb new triangles/quads/pentagons
         *
         * \param[in] ANb the number of cells to be created
         * \return    The index of the first created cell. Others are contiguous right behind.
         */
        TCellID addTriangles(const int ANb);
        TCellID addQuads(const int ANb);
        TCellID addPentagons(const int ANb);
        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate one new triangle/quad/pentagon.
         *
         * \return    The index of the created cell
         */
        TCellID addTriangle();
        TCellID addQuad();
        TCellID addPentagon();
        /*------------------------------------------------------------------------*/
        /** \brief  Create a cell from node ids
         *
         * \return The index of the created cell
         */
        TCellID newTriangle(const TCellID AN1, const TCellID AN2, const TCellID AN3);
        TCellID newFace(const TCellID AN1, const TCellID AN2, const TCellID AN3);
        TCellID newQuad(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);
        TCellID newFace(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);
        TCellID newPentagon(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                            const TCellID AN5);
        TCellID newFace(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5);

        TCellID newFace(const FakeFace AFF);
        TCellID newFace(const TCellID* AIDs, const TSize ASize);
        TCellID newFace(const Kokkos::View<TCellID*>& AN);

        /*------------------------------------------------------------------------*/
        /** \brief  Remove/Desallocate the face of id \p AId.
         *
         * \param[in] The index of the face to be removed.
         */
        void removeFace(const TCellID AId);
        /*------------------------------------------------------------------------*/
        /** \brief  Get the face \p AId
         *
         * \return  the face of id \p AId
         */
        Face getFace(const TCellID AId) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Check if the face \p AId exists in the mesh
         *
         * \return  true if it exists, fale otherwise
         */
        bool hasFace(const TCellID AId) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the capacity in terms of faces
         */
        TSize getFaceCapacity() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns more (superior or equal) than the max id of faces
         */
        TSize getFaceSupID() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of faces
         */
        TSize getNbFaces() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of triangles
         */
        TSize getNbTriangles() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of quads
         */
        TSize getNbQuads() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of pentagons
         */
        TSize getNbPentagons() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate \p ANb new tetrahedra/hexahedra/pyramids/prism3
         *
         * \param[in] ANb the number of cells to be created
         * \return    The index of the first created cell. Others are contiguous righ behind.
         */
        TCellID addTetrahedra(const int ANb);
        TCellID addHexahedra(const int ANb);
        TCellID addPyramids(const int ANb);
        TCellID addPrism3s(const int ANb);
        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate one new tetrahedra/hexahedra/pyramids/prism3
         *
         * \return    The index of the created cell
         */
        TCellID addTetrahedron();
        TCellID addHexahedron();
        TCellID addPyramid();
        TCellID addPrism3();
        /*------------------------------------------------------------------------*/
        /** \brief  Create a cell from node ids
         *
         * \return The index of the created cell
         */
        TCellID newTetrahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);
        TCellID newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);
        TCellID newHexahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                        const TCellID AN5, const TCellID AN6, const TCellID AN7, const TCellID AN8);
        TCellID newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                        const TCellID AN5, const TCellID AN6, const TCellID AN7, const TCellID AN8);
        TCellID newPyramid(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                            const TCellID AN5);
        TCellID newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5);
        TCellID newPrism3(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                       const TCellID AN5, const TCellID AN6);
        TCellID newRegion(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5, const TCellID AN6);
        /*------------------------------------------------------------------------*/
        /** \brief  Remove/Desallocate the region of id \p AId.
         *
         * \param[in] The index of the region to be removed.
         */
        void removeRegion(const TCellID AId);
        /*------------------------------------------------------------------------*/
        /** \brief  Get the region \p AId
         *
         * \return  the region of id \p AId
         */
        Region getRegion(const TCellID AId) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Check if the region \p AId exists in the mesh
         *
         * \return  true if it exists, fale otherwise
         */
        bool hasRegion(const TCellID AId) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the capacity in terms of regions
         */
        TSize getRegionCapacity() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns more (superior or equal) than the max id of regions
         */
        TSize getRegionSupID() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of regions
         */
        TSize getNbRegions() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of tetrahedra
         */
        TSize getNbTetrahedra() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of hexahedra
         */
        TSize getNbHexahedra() const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of pyramids
         */
        TSize getNbPyramids() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of prism3s
         */
        TSize getNbPrism3s() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns a container with the nodes IDs
         */
        void getNodeIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const;
        void getNodeIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns a container with the edges IDs
         */
        void getEdgeIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const;
        void getEdgeIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns a container with the faces IDs
         */
        void getFaceIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const;
        void getFaceIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns a container with the regions IDs
         */
        void getRegionIDs(kmds::GrowingView<kmds::TCellID>* ASelection) const;
        void getRegionIDs_dummy(kmds::GrowingView<kmds::TCellID>* ASelection) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Returns a container with the nodes IDs
     */
    void getNodeIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Returns a container with the edges IDs
     */
    void getEdgeIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Returns a container with the faces IDs
     */
    void getFaceIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Returns a container with the regions IDs
     */
    void getRegionIDsSeq(kmds::GrowingView<kmds::TCellID>* ASelection) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Create a new m_connectivity, which is associated to the cells of
         *          this mesh. The connectivity can be:
         *           - Downward: R2F, R2E, F2E,
         *           - Upward: N2R,N2F,N2E,E2R,E2F, F2R,
         *           - Dim equal: R2R, F2F, E2E, N2N
         *          Raise an exception if the connectivity cannot be created
         *
         *  \param AD  the type of connectivity to create.
         *
         *  \return A pointer on the connectivity
         */
        Connectivity* createConnectivity(const EMeshDefinition);

        /*------------------------------------------------------------------------*/
        /** \brief  Deletes a new m_connectivity, which is associated to the cells of
         *          this mesh. The connectivity can be:
         *           - Downward: R2F, R2E, F2E,
         *           - Upward: N2R,N2F,N2E,E2R,E2F, F2R,
         *           - Dim equal: R2R, F2F, E2E, N2N
         *          Raise an exception if the connectivity cannot be created
         *
         *  \param AD  the type of connectivity to create.
         *
         *  \return A pointer on the connectivity
         */
        Connectivity* deleteConnectivity(const EMeshDefinition);

        /*------------------------------------------------------------------------*/
        /** \brief  Give access to a connectivity, which is associated to the cells
         *          of this mesh. The connectivity can be:
         *           - Downward: R2F, R2E, F2E,
         *           - Upward: N2R,N2F,N2E,E2R,E2F, F2R,
         *           - Dim equal: R2R, F2F, E2E, N2N
         *          Raise an exception if the connectivity does not exist
         *
         *  \param AD  the type of connectivity to create.
         *
         *  \return A pointer on the connectivity
         */
        Connectivity* getConnectivity(const EMeshDefinition);

        /*------------------------------------------------------------------------*/
        /** \brief  Create a new variable attached to a generic cell type
         * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION)
         *
         *  \param AType the cell type
         *  \param AName the name of the variable. If this name is already used, the
         *  	   variable is not created and an exception is thrown
         *
         *  \return A pointer on the variable
         */
        template <typename T>
        Variable<T>* createVariable(const T ADefaultValue, const ECellType AType, const std::string& AName);

        /*------------------------------------------------------------------------*/
        /** \brief  Access to variable attached to a generic cell type
         * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION)
         *
         *  \param AType the cell type
         *  \param AName the name of the variable. If this name does not exist, the
         *  	   an exception is thrown
         *
         *  \return A pointer on the variable. This pointer can be null if the
         *  		specified type is wrong
         */
        template <typename T>
        Variable<T>* getVariable(const ECellType AType, const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Delete a variable attached to a cell type
     * 			(KMDS_NODE, KMDS_EDGE, KMDS_FACE or KMDS_REGION)
     *
     *  \param AType the cell type
     *  \param AName the name of the variable. If this name does not exist,
     *  	   an exception is thrown
     *
     */
    void deleteVariable(const ECellType AType, std::string AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Returns the variable manager attached to a cell type
     * 			(KMDS_NODE, KMDS_EDGE, KMDS_FACE or KMDS_REGION)
     *
     *  \param AType the cell type
     *
     *  \return the corresponding variable manager
     *
     */
    VariableManager* getVariableManager(const ECellType AType);

        // /*------------------------------------------------------------------------*/
        // /** \brief  Delete a variable attached to a cell type (GMDS_NODE, GMDS_EDGE,
        //  * 			GMDS_FACE or GNDS_REGION)
        //  *
        //  *  \param AType the cell type
        //  *  \param AName the name of the variable to be deleted
        //  */
        // template <typename C>
        // void deleteVariable(const std::string& AName);

        // /** \brief  Delete a variable attached to a cell type (GMDS_NODE, GMDS_EDGE,
        //  * 			GMDS_FACE or GNDS_REGION)
        //  *
        //  *  \param AType the cell type
        //  *  \param AVar a pointer on the variable to be deleted
        //  */
        // template <typename C>
        // void deleteVariable(VariableItf* AVar);

 private:
        /** node container */
        NContainer m_N;
        /** edge container */
        EContainer m_E;
        /** face container */
        FContainer m_F;
        /** region container */
        RContainer m_R;

        /** variable manager, 1 per dimension (node, edge, face, region)*/
        VariableManager m_var_manager[4];
        /** Optional connectivity is stored/and accessible by type*/
        std::map<EMeshDefinition, Connectivity*> m_connectivity;
};

/*----------------------------------------------------------------------------*/
template <typename T>
Variable<T>*
Mesh::createVariable(const T ADefaultValue, const ECellType AType, const std::string& AName)
{
        Variable<T>* v;
        switch (AType) {
        case KMDS_NODE: {
                v = m_var_manager[0].createVariable<T>(ADefaultValue, AName, m_N.capacity());
        } break;
        case KMDS_EDGE: {
                v = m_var_manager[1].createVariable<T>(ADefaultValue, AName, m_E.capacity());
        } break;
        case KMDS_FACE: {
                v = m_var_manager[2].createVariable<T>(ADefaultValue, AName, m_F.capacity());
        } break;
        case KMDS_REGION: {
                v = m_var_manager[3].createVariable<T>(ADefaultValue, AName, m_R.capacity());
        } break;
        default:
                throw KException("Unmanaged type of cell -> impossible to create a variable");
        }

        return v;
}

/*----------------------------------------------------------------------------*/
template <typename T>
Variable<T>*
Mesh::getVariable(ECellType AType, const std::string& AName)
{
        Variable<T>* v;
        switch (AType) {
        case KMDS_NODE: {
                v = m_var_manager[0].getVariable<T>(AName);
        } break;
        case KMDS_EDGE: {
                v = m_var_manager[1].getVariable<T>(AName);
        } break;
        case KMDS_FACE: {
                v = m_var_manager[2].getVariable<T>(AName);
        } break;
        case KMDS_REGION: {
                v = m_var_manager[3].getVariable<T>(AName);
        } break;
        default:
                throw KException("Unmanaged type of cell -> impossible to access to a variable");
        }
        return v;
}
/*----------------------------------------------------------------------------*/
}  // end namespace kmds
/*----------------------------------------------------------------------------*/
#endif /* KMDS_MESH_H_ */
/*----------------------------------------------------------------------------*/
