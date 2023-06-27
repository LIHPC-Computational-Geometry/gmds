/*----------------------------------------------------------------------------*/
/*
 * FACPoint.h
 *
 *  Created on: 27 juin 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDPOINT_H_
#define GMDS_GEOM_FACETEDPOINT_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "gmds/cad/GeomPoint.h"
#include "gmds/math/Point.h"
#include "gmds/ig/Mesh.h"
#include "gmds/ig/Node.h"
#include "GMDSCadFac_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
        class FACCurve;
/*----------------------------------------------------------------------------*/
/** \class FACPoint
 *  \brief This class implements the point services that are required by the
 *  	   mesh to the geometrical model.
 *  \param TBase the data type used to store geometrical data
 */
/*----------------------------------------------------------------------------*/
        class GMDSCadFac_API FACPoint : public GeomPoint {

        public:

            /*------------------------------------------------------------------------*/
            /** \brief  Default Constructor
             *  @param AMesh mesh suport
             */
            FACPoint(Mesh* AMesh);

            /*------------------------------------------------------------------------*/
            /** \brief  Constructor
             *  @param AMesh mesh suport
             *  \param ANode a mesh node id provided by the FacetedModel instance
             */
            FACPoint(Mesh* AMesh, TCellID ANode, const std::string& AName="Unknown point");

            /*------------------------------------------------------------------------*/
            /** \brief  Destructor
             */
            virtual ~FACPoint();


            /*------------------------------------------------------------------------*/
            /** \brief  Access to X coordinate
             *
             *  \return value of the X coordinate
             */
            virtual TCoord X() const;

            /*------------------------------------------------------------------------*/
            /** \brief  Access to Y coordinate
             *
             *  \return value of the Y coordinate
             */
            virtual TCoord Y() const;

            /*------------------------------------------------------------------------*/
            /** \brief  Access to Z coordinate
             *
             *  \return value of the Z coordinate
             */
            virtual TCoord Z() const;

            /*------------------------------------------------------------------------*/
            /** \brief  Access to X, Y and Z coordinates
             *
             *  \param  ACoordinates will receive the value of the X, Y and Z coordinates
             */
            virtual void XYZ(TCoord ACoordinates[3]) const;

            /*------------------------------------------------------------------------*/
            /** \brief  computes the area of the entity.
             */
            virtual TCoord computeArea() const;

            /*------------------------------------------------------------------------*/
            /** \brief  computes the bounding box
             *
             *	\param minXYZ The minimum coordinate of the bounding box.
             *	\param maxXYZ The maximum coordinate of the bounding box.
             */
            virtual void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const;
	         virtual std::tuple<TCoord,TCoord,TCoord,TCoord,TCoord,TCoord>  BBox() const;

            /*------------------------------------------------------------------------*/
            /** \brief  Access to the point as a NumericPoint
             *
             *  \return a numeric point
             */
            gmds::math::Point point() const;

            /*------------------------------------------------------------------------*/
            /** \brief  Set the corresponding mesh node
             *  \param ANode a mesh node provided by the FacetedModel instance
             */
            void set(Node ANode);

            /*------------------------------------------------------------------------*/
            /** \brief  Access to the point as a Mesh Node
             *
             *  \return a pointer onto a node
             */
            Node getNode() const;

            int id() const;

            std::vector<int> getCurveIds() const;


            void activateMultiNodes() {
                isMultiNodes_ = true;
            }

            void setMeshNodes(std::vector<Node> ANodes) {
                pnts_.clear();
                pnts_.resize(ANodes.size());
                for(int i=0; i<ANodes.size(); i++) {
                    pnts_[i] = ANodes[i].id();
                }
            }

            void project(gmds::math::Point& AP) const{
                if(isMultiNodes_) {

                    double dist2_max = HUGE_VALF;
                    gmds::math::Point pt_res(AP);

                    for(int i_p=0; i_p<pnts_.size(); i_p++) {
                        gmds::math::Point pt_tmp = m_support->get<gmds::Node>(pnts_[i_p]).point();

                        double dist2 = AP.distance2(pt_tmp);
                        if(dist2 < dist2_max) {
                            dist2_max = dist2;
                            pt_res = pt_tmp;
                        }
                    }

                    AP = pt_res;

                } else {
                    AP.setXYZ(X(),Y(),Z());
                }
            }

            /**@brief Accessor to the adjacent curves. Warning, there is no
             *  assumption about the ordering
             * @return curves that are adjacent to this point
             */
            virtual std::vector<GeomCurve*>& curves();
            /**@brief Accessor to the adjacent surfaces. Warning, there is no
             *  assumption about the ordering
             * @return surfaces that are adjacent to this point
             */
            virtual std::vector<GeomSurface*>& surfaces();
            /**@brief Accessor to the adjacent volumes. Warning, there is no
             *  assumption about the ordering
             * @return volumes that are adjacent to this point
             */
            virtual std::vector<GeomVolume*>& volumes();

            /**@brief Reset the global id counter to 1.
             */
            static void resetIdCounter();
        private:

            int m_id;
            static int m_next_id;

            Mesh* m_support;
            TCellID m_node_id;

            // for when a geom vertex can have several positions,
            // that happens with deduced geom models
            bool isMultiNodes_;
            std::vector<TCellID> pnts_;


            /** adjacent geometric curves*/
            std::vector<GeomCurve*> m_adjacent_curves;
            /** adjacent geometric surfaces*/
            std::vector<GeomSurface*> m_adjacent_surfaces;
            /** adjacent geometric volumes*/
            std::vector<GeomVolume*> m_adjacent_volumes;


        };
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDPOINT_H_ */
/*----------------------------------------------------------------------------*/
