/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/27/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MESHDOCTOR_H
#define GMDS_MESHDOCTOR_H
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIg_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class IGMeshDoctor
 *
 *  \brief this class provides algorithm to clean, check and modify a mesh.
 *
 */
    class GMDSIg_API MeshDoctor
    {
    public:

        /*------------------------------------------------------------------------*/
        /** \brief Constructor
         * \param AModel a mesh model defining available cells and  connectivities
         */
        explicit MeshDoctor(Mesh* AMesh);

        /*------------------------------------------------------------------------*/
        /** \brief Destructor
         */
        virtual ~MeshDoctor();

        /*------------------------------------------------------------------------*/
        /** \brief orient faces in the 2D case
         */
        int  orient2DFaces();
        /*------------------------------------------------------------------------*/
        /** \brief orient a 2D face
         */
        bool  orient2DFace(Face& AF);

        /*------------------------------------------------------------------------*/
        /** \brief create faces and R2F adjacency
         */
        void  buildFacesAndR2F() const;

        /*------------------------------------------------------------------------*/
        /** @brief create all the edges and X2E adjacency (where X=F or R)
         */
        void  buildEdgesAndX2E() const;

        /*------------------------------------------------------------------------*/
        /** @brief create only boundary edges in 2D and boundary edges and faces in
         *         3D. Warning, in 2D, requires to have N2F connection filled in a
         *         right manner.
         */
        void  buildBoundaryCells() const;

        void  updateUpwardConnectivity() const;

        /*------------------------------------------------------------------------*/
        /** \brief create faces
         */
        void  buildF() const;

        /*------------------------------------------------------------------------*/
        /** \brief create faces
         */
        void  buildE() const;

        /*------------------------------------------------------------------------*/
        /** \brief Fill up the X2Y connectivity for m_mesh
         *  \param ARefModel model that indicates the adjacency relationship that
         *  	   can be used.
         */
        void buildN2E(const MeshModel &ARefModel) const;
        void buildN2F(const MeshModel &ARefModel) const;
        void buildN2R(const MeshModel &ARefModel) const;
        void buildF2R(const MeshModel &ARefModel) const;
        void buildR2E(const MeshModel &ARefModel) const;
        void buildF2E(const MeshModel &ARefModel) const;

        /*------------------------------------------------------------------------*/
        /** \brief Change the mesh algorithms will be applied on
         */
        void setMesh(Mesh* AMesh);

    protected:

        /*------------------------------------------------------------------------*/
        /** \brief create faces and R2F
         */
        void  buildFAndR2F() const;


        static TCoord isLeft(Node& AN1, Node& AN2, Node& AN3);
        /*------------------------------------------------------------------------*/
        /** \brief Utilitary method to add a face if it exists
         *
         * \returns the id of the new or old face
         */
        static TCellID  addFace(Face& AFace,
                         std::map<VirtualFace::FaceID, TCellID>& AFakeFaceMap) ;

        /*------------------------------------------------------------------------*/
        /** \brief Utilitary method to add a face if it does not exist

         * \returns the id of the new or old face
         */
        TCellID  addFace(std::vector<TCellID>& ANodeIDs,
                         std::map<VirtualFace::FaceID, TCellID>& AFakeFaceMap) const;
        /*------------------------------------------------------------------------*/
        /** \brief Utilitary method to add an edge if it does not exist
         */
        void  addEdge(TCellID AN1, TCellID AN2,
                      std::map<VirtualEdge::EdgeID, TCellID>& AFakeEdgeMap) const;
    private:
        Mesh* m_mesh;


    };
/*----------------------------------------------------------------------------*/
} // namespace gmds
#endif //GMDS_MESHDOCTOR_H
