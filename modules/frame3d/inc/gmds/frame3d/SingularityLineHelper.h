/*----------------------------------------------------------------------------*/
#ifndef SH_SINGULARITY_LINE_HELPER_H_
#define SH_SINGULARITY_LINE_HELPER_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <gmds/frame3d/Params.h>
#include <gmds/math/AxisAngleRotation.h>
#include <gmds/math/Triangle.h>
#include "GMDSFrame3d_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
    namespace frame3d {
/*----------------------------------------------------------------------------*/
/** \class  SingularityLineHelper
 *  \brief  Starts from a tet mesh with a frame field, this algorithm extracts
 *          3D singularity lines of the frame field, build hexahedral regions
 *          around them, and finally modify the initial mesh to be the same
 *          but without new hexahedral regions.
 */
        class GMDSFrame3d_API SingularityLineHelper {
        public:

            struct PointVolumetricData{
                gmds::math::Point pnt;      // Point we are located at
                gmds::math::Vector3d dir;   // Direction we propage to
                gmds::Region tet;           // tetra where pnt is located in

                PointVolumetricData(const gmds::math::Point& AP,
                                    const gmds::math::Vector3d& AV,
                                    const gmds::Region& AR)
                        :pnt(AP),dir(AV),tet(AR){}
            };

            /*------------------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] AMesh background tetrahedral mesh we start from
             */
            SingularityLineHelper(Mesh *AMesh, ParamsGlobal &AGlobal, ParamsMark &AMarks);
            /*------------------------------------------------------------------------*/
            /** \brief Destructor.
             */
            virtual ~SingularityLineHelper();

            /*------------------------------------------------------------------------*/
            /** \brief Execute the algorithm
             */

            void execute();


            void moveBoundarySingularity(TCellID AFromFace, TCellID AToFace);

            void updateSmoothingData(std::set<TCellID>& ANodesToUpdate,
                                     Variable<math::Chart>* ANodeChartInfo,
                                     Variable<int>* ANodeVarConstraint,
                                     Variable<int>* ARegionsToUpdate);

            void initFrameSmoothing(std::set<TCellID>& ANodesToUpdate,
                                    Variable<int>* ANodeVarConstraint,
                                    Variable<int>* ARegionsToUpdate);

            bool isSingularTet(const TCellID AId);
        protected:

            double closestBndFace(const TCellID ANodeID, math::Vector3d AV, TCellID& ABndFaceID);
            std::vector<math::Vector3d> defineBoundarySlotsViaAngles(const Face& ABndFace,
                                                                     const math::Point& ASingLoc);

            std::vector<math::Vector3d> defineBoundarySlotsViaVectors(const Face& ABndFace,
                                                                      const math::Point& ASingLoc);

            bool isIn(const math::Point& AP, const Face& ATri,
                      bool& AOnEdge0, bool& AOnEdge1, bool& AOnEdge2);


            math::Chart computeChartIn(const math::Point& APnt,
                                       const Face&        AFace);

            void computeFuzzyHeuns(const math::Point&                AFromPnt,
                                   const math::Vector3d&             AFromDir,
                                   const std::vector<Face>&          AFaces,
                                   const std::vector<math::Triangle>&ATri,
                                   math::Point&                      AToPnt,
                                   math::Vector3d&                   AToDir,
                                   int&                              AToFaceId);
            bool followFlow(const Face AFromFace,
                            const PointVolumetricData& AData,
                            math::Point&               APnt,
                            std::vector<TCellID>& ACrossedTet,
                            std::map<TCellID,math::Chart>& ASingVal);
            /*------------------------------------------------------------------------*/
            /** \brief Detects the singular cells of the frame field
             */
            void detectSingularCells();
            /*------------------------------------------------------------------------*/
            /** \brief TransportBndSingularities
             */
            void transportBndSingularities();
            void transportBndSingularity(TCellID &AFID);

            int computeSingularityIndex(Face& AF);
            /*------------------------------------------------------------------------*/
            /** \brief initialize the used boolean marks specific to this algorithm
             */
            void initializeBooleanMarks();

            /*------------------------------------------------------------------------*/
            /** \brief clean the used boolean marks specific to this algorithm
             */
            void finalizeBooleanMarks();


            void createOneVolumeSingularityPoint(std::vector<gmds::Region> &ACluster);

            void createBoundarySingularityPoints();


            math::Vector3d getOutputNormal(const gmds::Face &AF, const gmds::Region &AR);

            math::Vector3d getInputNormal(const gmds::Face &AF, const gmds::Region &AR);


        protected:
            /** the mesh we work on*/
            Mesh *m_mesh;
            /** fhedo global parameters*/
            ParamsGlobal m_params_gl;


            /** rotation field defined on the vertices of m_mesh*/
            gmds::Variable<math::Chart> *m_chart_field;

            /** after line correction, we get constraints to ensure*/
            gmds::Variable<math::Chart> *m_new_constraint;
            gmds::Variable<int> *m_has_constraint;

            /** set of marks for boundary cells (nodes, edges, faces) */
            ParamsMark m_bm;

            int m_mark_line_sing;
            int m_mark_volume_pnt_sing;
            int m_mark_face_sing;
            int m_mark_edge_sing;
            int m_mark_node_sing;

            /** store the ids of singular region*/
            std::vector<gmds::TCellID> m_sing_tet;
            /** store for each  singular tet, the number of singular faces it has*/
            std::map<gmds::TCellID, int> m_sing_tet_type;
            /** store the ids of singular boundary faces*/
            std::vector<gmds::TCellID> m_sing_bnd_tri;
            std::map<gmds::TCellID, std::vector<math::Vector3d> > m_sing_bnd_seps;
        };
/*----------------------------------------------------------------------------*/
    }
}
/*----------------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------------*/