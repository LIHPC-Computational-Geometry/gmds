/*----------------------------------------------------------------------------*/
#ifndef GMDS_BOUNDARYEXTRACTOR_2D_H
#define GMDS_BOUNDARYEXTRACTOR_2D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  BoundaryExtractor2D
     *  @brief  Class that extracts the boundary of a 2D mesh.
     */
    class GMDSIgAlgo_API BoundaryExtractor2D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         *  @param[in] AFromMesh the 2D surfacic mesh we start from
         *  @param[in] AToMesh   the 2D linear mesh that is built by the algorithm
         */
        BoundaryExtractor2D(Mesh* AFromMesh, Mesh* AToMesh);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~BoundaryExtractor2D();

        /*------------------------------------------------------------------------*/
        /** @brief  check if the input mesh is valid for applying the 2D extraction
         *          algorithm
         */
        bool isValid() const;

        /*------------------------------------------------------------------------*/
        /** @brief execute the extraction algorithm
         */
        void execute();

        /*------------------------------------------------------------------------*/
        /** @brief  When we extract a mesh skin, we can expect to keep a mapping
         *          from the boundary cells of the surface mesh (nodes, edges)
         *          to the cells of the skin mesh. This mapping will be given into
         *          the map given in parameter here.
         *
         *  @param ANodeMap map from boundary nodes of the input mesh to the
         *                  skin nodes.
         *  @param AEdgeMap map from boundary edges of the input mesh to the
         *                  skin edges.
         *  @param ANodeMapInv  map from skin nodes of the input mesh to the
         *                      boundary nodes.
         *  @param AEdgeMapInv  map from skin edges of the input mesh to the
         *                      boundary edges.
         */
        void setMappings(std::map<TCellID,TCellID>* ANodeMap,
                         std::map<TCellID,TCellID>* AEdgeMap,
                         std::map<TCellID,TCellID>* ANodeMapInv,
                         std::map<TCellID,TCellID>* AEdgeMapInv);
        /*------------------------------------------------------------------------*/
        /** @brief  Gives the ability to pass to the algorithm variable to update
         *          indicating classification in AToMesh.
         *
         *          Warning, variables must be living on AToMesh. A node is always
         *          on the lower geometric entity first. For instance, a node
         *          classified onto a curve will not be known as being on the
         *          adjacent surfaces too.
         *
         *  @param ANodeColor node color indicating which point a node is on
         *  @param AEdgeColor edge color indicating which curve an edge is on
         *  @param ANodeColorCurv node color indicating which curve a node is on
         */
        void setColorOption(Variable<int>* ANodeColor,
                            Variable<int>* AEdgeColor,
                            Variable<int>* ANodeColorCurv);


    protected:

        /** a 3D volume mesh */
        Mesh* m_from_mesh;
        /** a 3D surface mesh */
        Mesh* m_to_mesh;
        /** color option is activated (true) or not (false)*/
        bool m_with_color;
        /** indicates to the algorithm to fill in some mapping info*/
        bool m_with_mapping;

        /** node, edge and face color */
        Variable<int>* m_color_node_on_pnt;
        Variable<int>* m_color_edge_on_curv;
        Variable<int>* m_color_node_on_curv;

        /** mappings for each boundary cell type from volume boundary to skin*/
        std::map<TCellID,TCellID>* m_node_map;
        std::map<TCellID,TCellID>* m_edge_map;
        /** mappings for each boundary cell type from skin to volume boundary*/
        std::map<TCellID,TCellID>* m_node_map_inv;
        std::map<TCellID,TCellID>* m_edge_map_inv;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_BOUNDARYEXTRACTOR3D_H
/*----------------------------------------------------------------------------*/

