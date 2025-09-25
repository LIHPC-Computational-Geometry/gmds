/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKMESHER_H
#define GMDS_BLOCKMESHER_H
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "gmds/cad/GeomMeshLinker.h"
#include "GMDSBlockMesher_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
    class GMDSBlockMesher_API BlockMesher{

    public:
        /*--------------------------------------------------------------------*/
        /** @enum  Status code for executing algorithms
         */
        typedef enum {
            FAIL,
            SUCCESS,
            NOT_YET_IMPLEMENTED,
            BLOCK_MODEL_ERROR,
            CLASSIFICATION_ERROR,
            BLOCK_VERTEX_MESHING_ERROR,
            BLOCK_EDGE_MESHING_ERROR,
            BLOCK_EDGE_PROJECTION_ERROR,
            BLOCK_FACE_MESHING_ERROR,
            BLOCK_FACE_PROJECTION_ERROR,
            BLOCK_CELL_MESHING_ERROR
        } STATUS;

        /** @brief  Default Constructor
         *  @param ABlocks the block structure we want to mesh
         *  @param ALinker Connection between ABlocks and a geometric model a geometry manager
         */
        BlockMesher(Mesh* ABlocks, cad::GeomMeshLinker* ALinker= NULL);

        /** @brief  Destructor
         */
        virtual ~BlockMesher();

        /*--------------------------------------------------------------------*/
        /** @brief  Mesh the block structure with splitting each block edge
         *          in ANb mesh edges
         *  @param ANb block edge discretization
         */
        STATUS execute(const int ANb = 10);
        /*--------------------------------------------------------------------*/
        /** @brief  Mesh the block structure with splitting each block edge
         *          to get mesh edges having size AEdgeTargetSize
         *  @param AEdgeTargetSize mesh edge target size
         */
        STATUS execute(const double AEdgeTargetSize);

        Mesh* mesh(){return m_mesh;}
        std::map<TCellID,Cell::Data> mesh_node_classification(){return m_mesh_node_classification;}
        std::map<TCellID,Cell::Data> mesh_edge_classification(){return m_mesh_edge_classification;}
        std::map<TCellID,Cell::Data> mesh_face_classification(){return m_mesh_face_classification;}
    private:
        bool meshVertices();
        bool meshEdges();
        bool meshFaces();
        bool meshBlocks();

        /*--------------------------------------------------------------------*/
        /** @brief  Return the edge which goes from @param AN0 to @param AN1
         *          among all the edges available in @param AEdges . The value
         *          of @param AInverted is true if the provided edge is oriented
         *          from @param AN1 to @param AN2, true otherwise.
         *  @param[in]  AN0         first node id
         *  @param[in]  AN1         second node id
         *  @param[in]  AEdges      A collection of edges, one being [AN0, AN1] or the inverse
         *  @param[out] AInverted   true if we get [AN1, AN0], false for [AN0, AN1]
         *  @param[out] AResult     the expected edge
         *  @return true if the edge is found, false otherwise
         */
        bool getEdgeFrom(const TCellID AN0, const TCellID AN1,
                         const std::vector<Edge>& AEdges,
                         Edge  &AResult,
                         bool& AInverted);
    private:
        /** the block structure*/
        Mesh* m_blocks;
        /** the geometry model*/
        cad::GeomMeshLinker* m_linker;
        /** a boolean value indicating that a geometric model is provided*/
        bool m_with_geometry;

        /** the generated mesh*/
        Mesh* m_mesh;

        //mapping
        std::map<TCellID ,TCellID> m_blockV_to_meshN;
        std::map<TCellID , std::vector<TCellID> > m_blockE_to_meshN;
        std::map<TCellID , std::vector<std::vector<TCellID> > > m_blockF_to_meshN;

        std::map<TCellID, Cell::Data> m_mesh_node_classification;
        std::map<TCellID, Cell::Data> m_mesh_edge_classification;
        std::map<TCellID, Cell::Data> m_mesh_face_classification;
        int m_nb_edges_in_discretization;
    };
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_BLOCKMESHER_H
/*----------------------------------------------------------------------------*/


