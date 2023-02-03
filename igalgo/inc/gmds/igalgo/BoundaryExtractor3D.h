/*----------------------------------------------------------------------------*/
#ifndef GMDS_BOUNDARYEXTRACTOR3D_H
#define GMDS_BOUNDARYEXTRACTOR3D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  BoundaryExtractor3D
     *  @brief  Class gathering operations used to get and mark the cells that
     belongs to the boundary of a 2D, 3D mesh
     */
    class GMDSIgAlgo_API BoundaryExtractor3D
            {
            public:

                /*------------------------------------------------------------------------*/
                /** @brief Constructor.
                 *
                 *  @param AFromMesh the 3D volumetric mesh we start from
                 *  @param AToMesh   the 3D surfacic mesh that is built by the algorithm
                 */
                BoundaryExtractor3D(Mesh* AFromMesh, Mesh* AToMesh);

                /*------------------------------------------------------------------------*/
                /** @brief  Destructor.	*/
                virtual ~BoundaryExtractor3D();


                /*------------------------------------------------------------------------*/
                /** @brief  check if the input mesh is valid for applying the 3D extraction
                 *          algorithm
                 */
                bool isValid() const;

                /*------------------------------------------------------------------------*/
                /** @brief execute the extraction algorithm
                 *
                 * @param AAngle surface dot angle
                 */
                void execute(double AAngle=0.72);


                /*------------------------------------------------------------------------*/
                /** @brief  When we extract a mesh skin, we can expect to keep a mapping
                 *          from the boundary cells of the volume mesh (nodes, edges, faces)
                 *          to the cells of the skin mesh. This mapping will be given into
                 *          the map given in parameter here.
                 *
                 *  @param ANodeMap map from boundary nodes of the input mesh to the
                 *                  skin nodes.
                 *  @param AEdgeMap map from boundary edges of the input mesh to the
                 *                  skin edges.
                 *  @param AFaceMap map from boundary faces of the input mesh to the
                 *                  skin faces.*
                 *  @param ANodeMapInv  map from skin nodes of the input mesh to the
                 *                      boundary nodes.
                 *  @param AEdgeMapInv  map from skin edges of the input mesh to the
                 *                      boundary edges.
                 *  @param AFaceMapInv  map from skin faces of the input mesh to the
                 *                      boundary faces.
                 */
                void setMappings(std::map<TCellID,TCellID>* ANodeMap,
                                 std::map<TCellID,TCellID>* AEdgeMap,
                                 std::map<TCellID,TCellID>* AFaceMap,
                                 std::map<TCellID,TCellID>* ANodeMapInv,
                                 std::map<TCellID,TCellID>* AEdgeMapInv,
                                 std::map<TCellID,TCellID>* AFaceMapInv);
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
                 *  @param AFaceColor face color indicating which surface a face is on
                 *  @param ANodeColorSurf node color indicating which surface a node is on
                 *  @param ANodeColorCurv node color indicating which curve a node is on
                 */
                void setColorOption(Variable<int>* ANodeColor,
                                    Variable<int>* AEdgeColor,
                                    Variable<int>* AFaceColor,
                                    Variable<int>* ANodeColorSurf,
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
                Variable<int>* m_color_face_on_surf;
                Variable<int>* m_color_node_on_surf;
                Variable<int>* m_color_node_on_curv;

                /** mappings for each boundary cell type from volume boundary to skin*/
                std::map<TCellID,TCellID>* m_node_map;
                std::map<TCellID,TCellID>* m_edge_map;
                std::map<TCellID,TCellID>* m_face_map;
                /** mappings for each boundary cell type from skin to volume boundary*/
                std::map<TCellID,TCellID>* m_node_map_inv;
                std::map<TCellID,TCellID>* m_edge_map_inv;
                std::map<TCellID,TCellID>* m_face_map_inv;
            };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_BOUNDARYEXTRACTOR3D_H
/*----------------------------------------------------------------------------*/
