/*----------------------------------------------------------------------------*/
/** \file    FacetedPoint.t.h
 *  \author  N. LE GOFF
 *  \date    03/04/2012
 */
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/cad/FACPoint.h>
#include <gmds/cad/FACCurve.h>
#include <gmds/ig/Mesh.h>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
	namespace cad{
/*----------------------------------------------------------------------------*/
		int FACPoint::m_next_id=0;
/*----------------------------------------------------------------------------*/
        void FACPoint::resetIdCounter() {m_next_id=1;}
/*----------------------------------------------------------------------------*/
		FACPoint::FACPoint(Mesh* AMesh)
				: m_id(m_next_id++), m_support(AMesh), isMultiNodes_(false)
		{

		}
/*----------------------------------------------------------------------------*/
		FACPoint::FACPoint(Mesh* AMesh,
						   TCellID ANode, const std::string& AName)
				: GeomPoint(AName),m_id(m_next_id++), m_support(AMesh),m_node_id(ANode), isMultiNodes_(false)
		{

		}
/*----------------------------------------------------------------------------*/
		FACPoint::~FACPoint()
		{

		}
/*----------------------------------------------------------------------------*/
		TCoord
		FACPoint::X() const
		{
			return m_support->get<Node>(m_node_id).X();
		}
/*----------------------------------------------------------------------------*/
		TCoord
		FACPoint::Y() const
		{
			return m_support->get<Node>(m_node_id).Y();
		}
/*----------------------------------------------------------------------------*/
		TCoord
		FACPoint::Z() const
		{
			return m_support->get<Node>(m_node_id).Z();
		}
/*----------------------------------------------------------------------------*/
		void
		FACPoint::XYZ(TCoord ACoordinates[3]) const
		{
			math::Point p = m_support->get<Node>(m_node_id).point();
			ACoordinates[0] = p.X();
			ACoordinates[1] = p.Y();
			ACoordinates[2] = p.Z();
		}
/*----------------------------------------------------------------------------*/
		TCoord
		FACPoint::computeArea() const
		{
			return 0;
		}
/*----------------------------------------------------------------------------*/
		void
		FACPoint::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
		{
			math::Point p = m_support->get<Node>(m_node_id).point();

			minXYZ[0]=p.X(); maxXYZ[0]=p.X();
			minXYZ[1]=p.Y(); maxXYZ[1]=p.Y();
			minXYZ[2]=p.Z(); maxXYZ[2]=p.Z();
		}
/*----------------------------------------------------------------------------*/
		gmds::math::Point
		FACPoint::point() const
		{
			return m_support->get<Node>(m_node_id).point();
		}
/*----------------------------------------------------------------------------*/
		void
		FACPoint::set(Node ANode)
		{
			m_node_id = ANode.id();
		}
/*----------------------------------------------------------------------------*/
		Node
		FACPoint::getNode() const
		{
			return m_support->get<Node>(m_node_id);
		}
/*----------------------------------------------------------------------------*/
		int
		FACPoint::id() const
		{
			return m_id;
		}
/*----------------------------------------------------------------------------*/
		std::vector<int> FACPoint::getCurveIds() const{
			Node n = getNode();
			std::vector<int> curve_ids;
			std::vector<Edge> edges = n.get<Edge>();
			Variable<int>* edge_on_crv = m_support->getVariable<int, GMDS_EDGE>("on_curve");
			for(auto e:edges){
				if((*edge_on_crv)[e.id()]!=0){
					curve_ids.push_back((*edge_on_crv)[e.id()]);
				}
			}
            return curve_ids;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomCurve*>& FACPoint::curves() {
            return m_adjacent_curves;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomSurface*>& FACPoint::surfaces() {
            return m_adjacent_surfaces;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomVolume*>& FACPoint::volumes() {
            return m_adjacent_volumes;
		}
/*----------------------------------------------------------------------------*/
	} // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
