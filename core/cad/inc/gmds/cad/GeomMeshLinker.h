/*----------------------------------------------------------------------------*/
/*
 * GeomMeshLinker
 *
 *  Created on:  04/13/2019
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOMMESHLINKER_H
#define GMDS_GEOMMESHLINKER_H
/*----------------------------------------------------------------------------*/
// gmds File headers
/*----------------------------------------------------------------------------*/
#include "gmds/utils/CommonTypes.h"
#include "gmds/utils/Exception.h"

#include "gmds/math/Vector.h"

#include "gmds/cad/GeomCurve.h"
#include "gmds/cad/GeomManager.h"
#include "GMDSCad_export.h"

#include "gmds/ig/Edge.h"
#include "gmds/ig/Node.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
/** @class  GeomMeshLinker
 *  @brief  This class keep the classification link between a mesh of type
 *          Mesh and a geom model GeomManager.
 *
 *          Node, edge and face classification are handled since they are
 *          mandatory for many algorithms.
 *
 */
/*----------------------------------------------------------------------------*/
        class  GMDSCad_API GeomMeshLinker
        {

        public:

            enum eLink { NoLink =0, LinkPoint =1, LinkCurve =2, LinkSurface =3, LinkVolume =4
            };
            GeomMeshLinker();
            /*------------------------------------------------------------------------*/
             /**@brief Constructor, build a link between a mesh and an geometry
              *
              * @param AMesh
              * @param AGeometry
              */
            GeomMeshLinker(Mesh* AMesh, GeomManager* AGeometry);

            /*------------------------------------------------------------------------*/
            /** \brief  Destructor
             */
            virtual ~GeomMeshLinker();

            void clear();
            Mesh* mesh() {return m_mesh;}
            GeomManager* geometry() {return m_geometry;}
            void setMesh(Mesh* AMesh);

            /*------------------------------------------------------------------------*/
            /**@brief For debug purpose, we write an output mesh with only the
             *        classified cells (nodes, edge, faces) that are linked to points,
             *        curves or surfaces. Each cell type will have two exported data
             *        variable, which are classif_entity_dim and classify_entity_id.
             *
             *        As the exported mesh is an extraction of the initial mesh, we also
             *        keep the reference id of each cell
             *        Warning, the file must be a vtk file (.vtk)
             *
             * @param AFileName filename to store the mesh in vtk format
             */
            void writeVTKDebugMesh(const std::string& AFileName);

            void setGeometry(GeomManager* AGeometry);

            /*------------------------------------------------------------------------*/
            /**@brief classify node AN onto point AGeomID
             *
             * @param AN        A node of m_mesh
             * @param AGeomId   point id in the geom manager
             */
            void linkToPoint(const Node& AN, int AGeomId);
            /*------------------------------------------------------------------------*/
            /**@brief classify node AN onto point AGeomID
             *
             * @param AN        A node id of m_mesh
             * @param AGeomId   point id in the geom manager
             */
            void linkNodeToPoint(const TCellID & AN, int AGeomId);
            /*------------------------------------------------------------------------*/
            /**@brief classify node AN onto curve AGeomID
             *
             * @param AN        A node of m_mesh
             * @param AGeomId   curve id in the geom manager
             */
            void linkToCurve(const Node& AN, int AGeomId);
            /*------------------------------------------------------------------------*/
            /**@brief classify node AN onto curve AGeomID
             *
             * @param AN        A node id of m_mesh
             * @param AGeomId   curve id in the geom manager
             */
            void linkNodeToCurve(const TCellID & AN, int AGeomId);

            /*------------------------------------------------------------------------*/
            /**@brief classify edge AE onto curve AGeomID
             *
             * @param AE        An edge of m_mesh
             * @param AGeomId   curve id in the geom manager
             */
            void linkToCurve(const Edge& AE, int AGeomId);
            /*------------------------------------------------------------------------*/
            /**@brief classify edge AE onto curve AGeomID
             *
             * @param AE        An edge id of m_mesh
             * @param AGeomId   curve id in the geom manager
             */
            void linkEdgeToCurve(const TCellID & AE, int AGeomId);
            /*------------------------------------------------------------------------*/
            /**@brief classify node AN onto surface AGeomID
             *
             * @param AN        A node of m_mesh
             * @param AGeomId   surface id in the geom manager
             */
            void linkToSurface(const Node& AN, int AGeomId);

            /*------------------------------------------------------------------------*/
            /**@brief classify node AN onto surface AGeomID
             *
             * @param AN        A node id of m_mesh
             * @param AGeomId   surface id in the geom manager
             */
            void linkNodeToSurface(const TCellID & AN, int AGeomId);
            /*------------------------------------------------------------------------*/
            /**@brief classify edge AE onto surface AGeomID
             *
             * @param AE        An edge of m_mesh
             * @param AGeomId   surface id in the geom manager
             */
            void linkToSurface(const Edge& AE, int AGeomId);

            /*------------------------------------------------------------------------*/
            /**@brief classify edge AE onto surface AGeomID
             *
             * @param AE        An edge id of m_mesh
             * @param AGeomId   surface id in the geom manager
             */
            void linkEdgeToSurface(const TCellID & AE, int AGeomId);
            /*------------------------------------------------------------------------*/
            /**@brief classify face AF onto surface AGeomID
             *
             * @param AF        A face of m_mesh
             * @param AGeomId   surface id in the geom manager
             */
            void linkToSurface(const Face& AF, int AGeomId);

            /*------------------------------------------------------------------------*/
            /**@brief classify face AF onto surface AGeomID
             *
             * @param AF        A face id of m_mesh
             * @param AGeomId   surface id in the geom manager
             */
            void linkFaceToSurface(const TCellID & AF, int AGeomId);

            /*------------------------------------------------------------------------*/
            /**@ brief accessor on the dimension of the geom entity AN is classified on.
             *
             * @param AN a node
             * @return the dimension of the geom entity AN is classified on
             */
	         eLink getGeomDim(const Node& AN);
	         eLink getGeomDim(const Edge& AE);
	         eLink getGeomDim(const Face& AF);
            /*------------------------------------------------------------------------*/
            /**@ brief accessor on the id of the geom entity AN is classified on.
             *
             * @param AN a node
             * @return the id of the geom entity AN is classified on
             */
            int getGeomId(const Node& AN);
            int getGeomId(const Edge& AN);
            int getGeomId(const Face& AN);
            /*------------------------------------------------------------------------*/
            /**@ brief accessor on the dimension of the geom entity AN is classified on.
             *
             * @param AN a node id
             * @return the dimension of the geom entity AN is classified on
             */
            template<typename T> GMDSCad_API eLink getGeomDim(const TCellID& AN);
            /*------------------------------------------------------------------------*/
            /**@ brief accessor on the id of the geom entity AN is classified on.
             *
             * @param AN a node id
             * @return the dimension of the geom entity AN is classified on
             */
            template<typename T> GMDSCad_API int getGeomId(const TCellID& AN);

            /*------------------------------------------------------------------------*/
            /**@ brief accessor on the dimension and id of the geom entity AN is
             *   classified on.
             *
             * @param AN a node
             * @return the (dim, id) of the geom entity AN is classified on
             */
            std::pair<eLink,int> getGeomInfo(const Node& AN);
            std::pair<eLink,int> getGeomInfo(const Edge& AE);
            std::pair<eLink,int> getGeomInfo(const Face& AF);

            /*------------------------------------------------------------------------*/
            /**@ brief accessor on the dimension and id of the geom entity AN is
             *   classified on.
             *
             * @param AN a node id
             * @return the (dim, id) of the geom entity AN is classified on
             */
            template<typename T> GMDSCad_API  std::pair<eLink,int> getGeomInfo(const TCellID& AN);

        private:

            /** global id to define in a robust way each linker instance*/
            static int m_global_link_id;
            /** id to define which linker it is*/
            int m_link_id;
            /** mesh that is cliassifed on the geometry*/
            Mesh*  m_mesh;
            /** geometric model we work on*/
            GeomManager* m_geometry;
            /** variable storing the dim. of the geom entity each node is linked to*/
            Variable<eLink>* m_node_classification_dim{};
            /** variable storing the id of the geom entity each node is linked to*/
            Variable<int>* m_node_classification_id{};
            /** variable storing the dim. of the geom entity each edge is linked to*/
            Variable<eLink>* m_edge_classification_dim{};
            /** variable storing the id of the geom entity each edge is linked to*/
            Variable<int>* m_edge_classification_id{};
            /** variable storing the dim. of the geom entity each face is linked to*/
            Variable<eLink>* m_face_classification_dim{};
            /** variable storing the id of the geom entity each face is linked to*/
            Variable<int>* m_face_classification_id{};
        };

        template <> GMDSCad_API GeomMeshLinker::eLink GeomMeshLinker::getGeomDim<Node>(const TCellID &AN);
        template <> GMDSCad_API GeomMeshLinker::eLink GeomMeshLinker::getGeomDim<Edge>(const TCellID &AN);
        template <> GMDSCad_API GeomMeshLinker::eLink GeomMeshLinker::getGeomDim<Face>(const TCellID &AN);

        template <> GMDSCad_API int GeomMeshLinker::getGeomId<Edge>(const TCellID &AN);
        template <> GMDSCad_API int GeomMeshLinker::getGeomId<Node>(const TCellID &AN);
        template <> GMDSCad_API int GeomMeshLinker::getGeomId<Face>(const TCellID &AN);

        template <> GMDSCad_API std::pair<GeomMeshLinker::eLink,int> GeomMeshLinker::getGeomInfo<Node>(const TCellID& AN);
        template <> GMDSCad_API std::pair<GeomMeshLinker::eLink,int> GeomMeshLinker::getGeomInfo<Edge>(const TCellID& AN);
        template <> GMDSCad_API std::pair<GeomMeshLinker::eLink,int> GeomMeshLinker::getGeomInfo<Face>(const TCellID& AN);

/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_GEOMMESHLINKER_H
/*----------------------------------------------------------------------------*/

