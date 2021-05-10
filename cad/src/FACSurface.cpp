/*----------------------------------------------------------------------------*/
/** \file    FacetedSurface.t.h
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
#include "gmds/cad/FACSurface.h"
/*----------------------------------------------------------------------------*/
#include <gmds/cad/FACManager.h>
#include "gmds/cad/FACVolume.h"
#include "gmds/cad/FACCurve.h"
#include "gmds/cad/FACPoint.h"
#include "gmds/math/Ray.h"
#include "gmds/math/Triangle.h"
#include "gmds/utils/RandomGenerator.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
        int FACSurface::m_next_id=1;
/*----------------------------------------------------------------------------*/
        FACSurface::FACSurface(Mesh* AMesh)
                : m_support(AMesh)
        {}
/*----------------------------------------------------------------------------*/

        FACSurface::
        FACSurface(Mesh* AMesh,
                   std::vector<TCellID >& ADiscret,
                   const std::string& AName)
                :GeomSurface(AName),m_support(AMesh),m_id(m_next_id++), m_mesh_faces(ADiscret)
        {
            for(auto f_id : m_mesh_faces){
                std::vector<Node> nodes_fi;
                m_support->get<Face>(f_id).get<Node>(nodes_fi);

            }
        }
/*----------------------------------------------------------------------------*/

        FACSurface::
        ~FACSurface()
        {}
/*----------------------------------------------------------------------------*/
        TCoord FACSurface::computeArea() const
        {
            TCoord totalArea = 0.;

            for(unsigned int i=0;i<m_mesh_faces.size();i++){

                totalArea +=  m_support->get<Face>(m_mesh_faces[i]).area();
            }
            return totalArea;
        }
/*----------------------------------------------------------------------------*/

        void FACSurface::
        computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
        {
            std::set<TCellID> node_ids;
            for(unsigned int i=0;i<m_mesh_faces.size();i++){

                std::vector<TCellID> f_ids = m_support->get<Face>(m_mesh_faces[i]).getIDs<Node>();
                node_ids.insert(f_ids.begin(),f_ids.end());
            }

            // too many comparisons (3 factor)
            math::Point pi = m_support->get<Node>(*node_ids.begin()).getPoint();
            minXYZ[0]=pi.X();
            maxXYZ[0]=pi.X();
            minXYZ[1]=pi.Y();
            maxXYZ[1]=pi.Y();
            minXYZ[2]=pi.Z();
            maxXYZ[2]=pi.Z();
            for(auto i:node_ids){
                math::Point pi = m_support->get<Node>(i).getPoint();
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

        void FACSurface::
        computeNormal( const math::Point& AP,
                       math::Vector3d& AV) const
        {
            Face f0 = m_support->get<Face>(m_mesh_faces[0]);
            TCoord min_dist = f0.distance(AP);
            TCellID index = 0;
            for(unsigned int i=1;i<m_mesh_faces.size();i++){
                TCoord dist = m_support->get<Face>(m_mesh_faces[i]).distance(AP);
                if(dist<min_dist){
                    min_dist=dist;
                    index = i;
                }
            }
            AV = m_support->get<Face>(m_mesh_faces[index]).normal();
        }
/*----------------------------------------------------------------------------*/
        math::Point FACSurface::
        closestPoint(const math::Point& AP) const
        {
            Face f0 = m_support->get<Face>(m_mesh_faces[0]);
            TCoord min_dist = f0.distance(AP);
            TCellID index = 0;
            for(unsigned int i=1;i<m_mesh_faces.size();i++){
                TCoord dist = m_support->get<Face>(m_mesh_faces[i]).distance(AP);
                if(dist<min_dist){
                    min_dist=dist;
                    index = i;
                }
            }
            return m_support->get<Face>(m_mesh_faces[index]).project(AP);

        }
/*----------------------------------------------------------------------------*/
        math::Point FACSurface::center() const
        {
            throw GMDSException("Not yet implemented");

        }
/*----------------------------------------------------------------------------*/
        void FACSurface::
        setMeshFaces(const std::vector<Face>& AFaces)
        {
            m_mesh_faces.clear();
            m_mesh_faces.resize(AFaces.size());
            for(unsigned int i=0; i<AFaces.size(); i++) {
                m_mesh_faces[i] = AFaces[i].id();
            }
        }
/*----------------------------------------------------------------------------*/
        void FACSurface::
        getMeshFaces(std::vector<Face>& AFaces) const
        {
            AFaces.clear();
            AFaces.reserve(m_mesh_faces.size());
            for(auto f_id : m_mesh_faces){
                AFaces.push_back(m_support->get<Face>(f_id));
            }
        }
/*----------------------------------------------------------------------------*/
        void FACSurface::getTriangulation(std::vector<math::Triangle >& ATri) const
        {
            ATri.clear();
            ATri.resize(m_mesh_faces.size());
            for(unsigned int i=0;i<m_mesh_faces.size();i++)
            {
                gmds::Face current = m_support->get<gmds::Face>(m_mesh_faces[i]);
                std::vector<Node> nodes = current.get<Node>();
                ATri[i]=math::Triangle(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint());
            }
        }
/*----------------------------------------------------------------------------*/
        void
        FACSurface::propagateOrient(Face AFace, int AMarkTreatedFaces, int AMarkFacesInverted, Mesh* AMesh)
        {
            std::vector<Face> faces = AFace.get<Face>();

            for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
                if(!AMesh->isMarked(faces[iFace],AMarkTreatedFaces)) {
                    bool faceIsSameOrientation = checkSameOrientFace(AFace,faces[iFace]);

                    if(!faceIsSameOrientation) {
                        invertFace(faces[iFace]);
                        AMesh->mark(faces[iFace],AMarkFacesInverted);
                    }

                    AMesh->mark(faces[iFace],AMarkTreatedFaces);
                    propagateOrient(faces[iFace],AMarkTreatedFaces,AMarkFacesInverted,AMesh);
                }
            }
        }
/*----------------------------------------------------------------------------*/
        bool
        FACSurface::checkSameOrientFace(Face AFaceRef, Face AFaceCheck)
        {
            std::vector<Node> nodes = AFaceRef.get<Node>();
            std::vector<Node> nodesBis = AFaceCheck.get<Node>();

            for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {
                for(unsigned int iNodeBis=0; iNodeBis<nodesBis.size(); iNodeBis++) {
                    if(nodesBis[iNodeBis] == nodes[iNode]) {
                        if(nodesBis[(iNodeBis+1)%nodesBis.size()] == nodes[(iNode+1)%nodes.size()]) {
                            return true;
                        } else {
                            return false;
                        }
                    }
                }
            }

            throw GMDSException("FACSurface::checkOrientFace we should not be in this part of the code.");
        }
/*----------------------------------------------------------------------------*/
        void
        FACSurface::invertAllFaces()
        {
            std::vector<Face> faces;
            getMeshFaces(faces);
            for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
                invertFace(faces[iFace]);
            }
        }
/*----------------------------------------------------------------------------*/
        void
        FACSurface::invertFace(Face AFace)
        {
            std::vector<Node> nodes = AFace.get<Node>();
            std::vector<Node> nodestmp(nodes.rbegin(),nodes.rend());

            AFace.set<Node>(nodestmp);
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomPoint*>& FACSurface::points() {
            return m_adjacent_points;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomCurve*>& FACSurface::curves() {
            return m_adjacent_curves;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomVolume*>& FACSurface::volumes() {
            return m_adjacent_volumes;
        }
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
