/*----------------------------------------------------------------------------*/
/** \file    FacetedCurve.cpp
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#include <gmds/cad/FACCurve.h>
#include <gmds/cad/FACSurface.h>
#include <gmds/cad/FACPoint.h>
#include <gmds/cad/FACVolume.h>
#include <gmds/math/Constants.h>
#include <gmds/math/Plane.h>
#include <gmds/ig/Mesh.h>
#include <set>
#include <map>
#include <gmds/math/Line.h>
#include <gmds/math/Orientation.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
        int FACCurve::m_next_id=1;
/*----------------------------------------------------------------------------*/
        void FACCurve::resetIdCounter() {m_next_id=1;}
/*----------------------------------------------------------------------------*/
        FACCurve::FACCurve(Mesh* AMeshSupport)
                : m_support(AMeshSupport), m_id(m_next_id)
        {
            m_next_id++;
        }
/*----------------------------------------------------------------------------*/
        FACCurve::FACCurve(Mesh* AMeshSupport,
                           std::vector<TCellID>& APoints,
                           std::vector<TCellID>& AEdges,
                           const std::string& AName)
                :GeomCurve(AName),m_support(AMeshSupport),m_id(m_next_id++),
                 m_mesh_nodes(APoints),m_mesh_edges(AEdges)
        {

            for(unsigned int i=0;i<APoints.size()-1;i++){
                Node n1 = m_support->get<Node>(APoints[i]);
                Node n2 = m_support->get<Node>(APoints[i+1]);
                math::Point p1(n1.X(),n1.Y(),n1.Z());
                math::Point p2(n2.X(),n2.Y(),n2.Z());
            }

        }
/*----------------------------------------------------------------------------*/
        FACCurve::FACCurve(Mesh* AMeshSupport,
                           FACPoint* AP1,
                           FACPoint* AP2,
                           std::vector<TCellID>& APoints,
                           std::vector<TCellID>& AEdges,
                           const std::string& AName)
                :GeomCurve(AName),m_support(AMeshSupport),m_id(m_next_id++),
                 m_mesh_nodes(APoints),m_mesh_edges(AEdges)
        {

            for(unsigned int i=0;i<APoints.size()-1;i++){
                Node n1 = m_support->get<Node>(APoints[i]);
                Node n2 = m_support->get<Node>(APoints[i+1]);
                math::Point p1(n1.X(),n1.Y(),n1.Z());
                math::Point p2(n2.X(),n2.Y(),n2.Z());
            }

        }
/*----------------------------------------------------------------------------*/
        FACCurve::
        ~FACCurve()
        {}
/*----------------------------------------------------------------------------*/
        double FACCurve::length() const
        {
            double length = 0.;
            for(auto e_id : m_mesh_edges){
                Edge e = m_support->get<Edge>(e_id);
                length+=e.length();
            }

            return length;
        }
/*----------------------------------------------------------------------------*/
        void FACCurve::project(math::Point& AP) const
        {

            AP = closestPoint(AP);
        }
/*----------------------------------------------------------------------------*/
        GeomCurve::CurvatureInfo FACCurve::getCurvatureInfo() const {
            //We get through all the edges assigned to this curve
            int nb_flat =0;
            int nb_convex=0;
            int nb_concave=0;
            for(auto e_id:m_mesh_edges){
                Edge ei = m_support->get<Edge>(e_id);
                std::vector<Face> adj_faces = ei.get<Face>();
                if(adj_faces.size()==2) {
                    std::vector<Node> ns0 = adj_faces[0].get<Node>();
                    math::Vector3d v0 = adj_faces[0].normal();
                    math::Vector3d v1 = adj_faces[1].normal();
                    math::Point c1  = adj_faces[1].center();
                    bool convex=true;
                    if(math::Orientation::orient3d(c1, ns0[0].point(),
                                                   ns0[1].point(),
                                                   ns0[2].point()) == math::Orientation::NEGATIVE){
                        //its non-convex
                        convex=false;
                    }
                    double a = math::Constants::INVPIDIV180*v0.angle(v1);
                    if(!convex)
                        a+=180;
                    if(a<30)
                        nb_flat++;
                    else if(60<a && a<120)
                        nb_convex++;
                    else if( a> 220)
                        nb_concave++;
                }
            }

            if(nb_convex>nb_concave && nb_convex>nb_flat)
                return GeomCurve::CONVEX;
            if(nb_concave>nb_convex && nb_concave>nb_flat)
                return GeomCurve::CONCAVE;
            return GeomCurve::FLAT;
        }
/*----------------------------------------------------------------------------*/
        math::Point FACCurve::closestPoint(const math::Point &AP) const
        {
            double min_dist = HUGE_VAL;
            TCellID  min_id = NullID;
            for(auto e_id: m_mesh_edges){
                Edge e = m_support->get<Edge>(e_id);
                math:math::Point c= e.center();
                if(AP.distance2(c)<min_dist){
                    min_id = e_id;
                    min_dist = AP.distance2(c);
                }
            }

            //Warning, we could be unable to project
            if(min_id==NullID){
                return AP;
            }

            Edge e = m_support->get<Edge>(min_id);
            std::vector<Node> e_nodes = e.get<Node>();
            math::Segment seg(e_nodes[0].point(),e_nodes[1].point());


            return seg.project(AP);

/*            math::Point closest = e.center();
            double d= AP.distance2(closest);
            if(AP.distance2(e_nodes[0].point()) < d){
                closest = e_nodes[0].point();
                d = AP.distance2(e_nodes[0].point());
            }
            if (AP.distance2(e_nodes[1].point()) < d){
                closest = e_nodes[1].point();
            }
            return closest;*/
        }
/*----------------------------------------------------------------------------*/
        void FACCurve::project(math::Point& AP, math::Vector3d& AV) const
        {
            double min_dist = HUGE_VAL;
            TCellID min_id = NullID;

            for (auto e_id : m_mesh_edges) {
                Edge e = m_support->get<Edge>(e_id);
                std::vector<Node> e_nodes = e.get<Node>();
                math::Segment s(e_nodes[0].point(), e_nodes[1].point());
                double dist = s.distance(AP);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_id = e_id;
                }
            }
            Edge e = m_support->get<Edge>(min_id);
            std::vector<Node> e_nodes = e.get<Node>();
            math::Segment s(e_nodes[0].point(), e_nodes[1].point());
            AP = s.project(AP);
            AV = math::Vector3d(e_nodes[0].point(), e_nodes[1].point());
        }
/*----------------------------------------------------------------------------*/
        TCoord FACCurve::computeArea() const
        {
            throw GMDSException("Not yet implemented!");
        }
/*----------------------------------------------------------------------------*/
        void FACCurve::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
        {
            std::vector<math::Point > pnts;
            for(auto n_id:m_mesh_nodes){
                pnts.push_back(m_support->get<Node>(n_id).point());
            }
            math::Point pi = pnts[0];
            minXYZ[0]=pi.X();
            maxXYZ[0]=pi.X();
            minXYZ[1]=pi.Y();
            maxXYZ[1]=pi.Y();
            minXYZ[2]=pi.Z();
            maxXYZ[2]=pi.Z();
            for(unsigned int i=1;i<pnts.size();i++){
                pi = pnts[i];
                if (pi.X()<minXYZ[0])
                    minXYZ[0]=pi.X();
                else if (pi.X()>maxXYZ[0])
                    maxXYZ[0]=pi.X();

                if (pi.Y()<minXYZ[1])
                    minXYZ[1]=pi.Y();
                else if (pi.Y()>maxXYZ[1])
                    maxXYZ[1]=pi.Y();

                if (pi.Z()<minXYZ[2])
                    minXYZ[2]=pi.Z();
                else if (pi.Z()>maxXYZ[2])
                    maxXYZ[2]=pi.Z();
            }
        }
/*----------------------------------------------------------------------------*/
        bool FACCurve::isALoop()const
        {
            std::map<TCellID , int> nb_occ;
            for(auto e_id:m_mesh_edges){
                Edge e = m_support->get<Edge>(e_id);
                std::vector<TCellID > e_nodes = e.getIDs<Node>();
                nb_occ[e_nodes[0]]++;
                nb_occ[e_nodes[1]]++;
            }

            bool found_not_two = false;
            for(auto o:nb_occ){
                if(o.second==1){
                    found_not_two =true;
                }
            }

            // in case of a loop, all nodes appear twice
            return !found_not_two;
        }
/*----------------------------------------------------------------------------*/
        TCoord FACCurve::computeDihedralAngle() const
        {
            TCoord anglemax = 0.;

            // we go only to size-1 because we are interested in edges (segments)
            for(unsigned int iNode=0; iNode<m_mesh_nodes.size()-1; iNode++)
            {
                std::vector<Edge> edges_possibilities;
                m_support->get<Node>(m_mesh_nodes[iNode]).get<Edge>(edges_possibilities);

                // we discard edges that have only one point on the curve
                Edge current_edge;

                for(auto ei : edges_possibilities)
                {
                    std::vector<TCellID > nodes;
                    ei.getIDs<Node>(nodes);

                    if(nodes[0] == m_mesh_nodes[iNode])
                    {
                        if(nodes[1] == m_mesh_nodes[iNode+1])
                        {
                            current_edge = ei;
                            break;
                        }
                    }
                    else
                    {
                        if(nodes[0] == m_mesh_nodes[iNode+1])
                        {
                            current_edge = ei;
                            break;
                        }
                    }
                }
                std::vector<Face> faces;
                current_edge.get<Face>(faces);
                if(faces.size()!=2)
                    throw GMDSException("Dihedral angles can only be computed for edges having 2 adjacent faces");

                math::Vector3d n0 = faces[0].normal();
                math::Vector3d n1 = faces[1].normal();
                //norm are equals to 1

                // first get angle between [0,PI]
                TCoord angle = acos(-(n0.dot(n1)))*math::Constants::INVPIDIV180;

                // then determine the sign of the sin
                // it works because the degenerate case here would be around Pi, and if the angle is Pi the curve is
                // not a sharp curve.
                math::Point p0 = faces[0].center();
                math::Point p1 = faces[1].center();
                math::Point p = (math::Plane (p1,n1)).project(p0);

                TCoord sign = (math::Vector3d (p0,p)).dot(n1);

                if(sign < 0.)
                    angle = 360. - angle;

                if(angle > anglemax)
                    anglemax  = angle;
            }

            return anglemax;
        }
/*----------------------------------------------------------------------------*/
        void FACCurve::getMeshNodes(std::vector<Node>& ANodes) const {
            ANodes.clear();
            ANodes.reserve(m_mesh_nodes.size());
            for (auto i:m_mesh_nodes)
                ANodes.push_back(m_support->get<Node>(i));
        }
/*----------------------------------------------------------------------------*/
        void
        FACCurve::setMeshEdges(const std::vector<gmds::Edge> &AEdges)
        {
            m_mesh_edges.clear();
            m_mesh_edges.resize(AEdges.size());
            for(int i=0; i<AEdges.size(); i++) {
                m_mesh_edges[i] = AEdges[i].id();
            }
        }
/*----------------------------------------------------------------------------*/
        void FACCurve::getMeshEdges(std::vector<Edge>& AEdges) const
        {
            AEdges.clear();
            AEdges.reserve(m_mesh_edges.size());
            for(auto e_i:m_mesh_edges)
                AEdges.push_back(m_support->get<Edge>(e_i));
        }
/*----------------------------------------------------------------------------*/
        void FACCurve::getTriangulation(std::vector<math::Segment>& ASeg) const
        {
            ASeg.clear();
            ASeg.resize(m_mesh_edges.size());

            for(unsigned int i=0;i<m_mesh_edges.size();i++)
            {
                ASeg[i]=m_support->get<Edge>(m_mesh_edges[i]).segment();
            }
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomPoint*>& FACCurve::points() {
            return m_adjacent_points;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomSurface*>& FACCurve::surfaces() {
            return m_adjacent_surfaces;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomVolume*>& FACCurve::volumes() {
            return m_adjacent_volumes;
        }
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
