/*----------------------------------------------------------------------------*/
/** \file    FacetedSurface.t.h
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
// ANN File Headers
#include "ANN/ANN.h"
/*----------------------------------------------------------------------------*/
#include "gmds/cadfac/FACSurface.h"
/*----------------------------------------------------------------------------*/
#include <gmds/cadfac/FACManager.h>
#include "gmds/cadfac/FACVolume.h"
#include "gmds/cadfac/FACCurve.h"
#include "gmds/cadfac/FACPoint.h"
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
        void FACSurface::resetIdCounter() {m_next_id=1;}
/*----------------------------------------------------------------------------*/
        FACSurface::FACSurface(Mesh* AMesh)
                : m_support(AMesh), m_kd_tree(nullptr)
        {}
/*----------------------------------------------------------------------------*/

        FACSurface::
        FACSurface(Mesh* AMesh,
                   std::vector<TCellID >& ADiscret,
                   const std::string& AName)
                :GeomSurface(AName),m_support(AMesh),m_id(m_next_id++),
                m_mesh_faces(ADiscret),m_kd_tree(nullptr)
        {
            buildANNTree();
        }
/*----------------------------------------------------------------------------*/

        FACSurface::
        ~FACSurface()
        {
            //got to clean some technical attributes used by ANN
            delete m_kd_tree;
            annClose();									// done with ANN

        }
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
            math::Point pi = m_support->get<Node>(*node_ids.begin()).point();
            minXYZ[0]=pi.X();
            maxXYZ[0]=pi.X();
            minXYZ[1]=pi.Y();
            maxXYZ[1]=pi.Y();
            minXYZ[2]=pi.Z();
            maxXYZ[2]=pi.Z();
            for(auto i:node_ids){
                math::Point pi = m_support->get<Node>(i).point();
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
            Variable<int>* var_surf = m_support->getVariable<int, GMDS_FACE>("on_surface");

            Face seed = m_support->get<Face>(getANNClosestTriangle(AP));
            std::set<TCellID> face_ids;
            std::vector<Node> ns = seed.get<Node>();
            for(auto n:ns){
                std::vector<TCellID> f_ids = n.getIDs<Face>();
                for(auto i:f_ids)
                    if(var_surf->value(i)==this->id())
                        face_ids.insert(i);
            }

            Face f0 = m_support->get<Face>(*face_ids.begin());
            TCoord min_dist = f0.distance(AP);
            TCellID closest_face_id = f0.id();
            for(auto f_id:face_ids){
                TCoord dist = m_support->get<Face>(f_id).distance(AP);
                if(dist<min_dist){
                    min_dist=dist;
                    closest_face_id = f_id;
                }
            }
            return m_support->get<Face>(closest_face_id).project(AP);

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

            if(m_kd_tree!=NULL)
                delete m_kd_tree;

            buildANNTree();
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
                ATri[i]=math::Triangle(nodes[0].point(), nodes[1].point(), nodes[2].point());
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
        /*---------------------------------------------------------------------------*/
        void FACSurface::buildANNTree() {

            int dim = 3;                         // dimension
            int maxPts = m_mesh_faces.size();    // maximum number of data points
            ANNpointArray dataPts;               // data points
            dataPts = annAllocPts(maxPts, dim);  // allocate data points

            //========================================================
            // (1) Fill in the  ANN structure for storing points
            // Important: Points in APnts and dataPnts are stored with
            // same index.
            //========================================================
            int nPts=0;  // actual number of data points

            while (nPts < maxPts) {
                math::Point p = m_support->get<Face>(m_mesh_faces[nPts]).center();
                dataPts[nPts][0] = p.X();
                dataPts[nPts][1] = p.Y();
                dataPts[nPts][2] = p.Z();
                nPts++;
            };
            //========================================================
            // (2) Build the search structure
            //========================================================
            if (m_kd_tree != NULL)
                throw GMDSException("FACSurface Issue: the kd tree structure was previously initialized");
            m_kd_tree = new ANNkd_tree(dataPts,// the data points
                                       nPts,   // number of points
                                       dim);   // dimension of space
        }
        /*---------------------------------------------------------------------------*/
        TCellID FACSurface::getANNClosestTriangle(const math::Point& AP) const {
            int k = 1;      // max number of nearest neighbors
            int dim = 3;       // dimension

            ANNpoint queryPt;           // query point
            ANNidxArray nnIdx;          // near neighbor indices
            ANNdistArray dists;         // near neighbor distances
            queryPt = annAllocPt(dim);  // allocate 1 query point
            nnIdx = new ANNidx[k];      // allocate near neigh indices
            dists = new ANNdist[k];     // allocate near neighbor dists

            queryPt[0] = AP.X();
            queryPt[1] = AP.Y();
            queryPt[2] = AP.Z();

            m_kd_tree->annkSearch(	// search
                    queryPt,// query point
                    k,
                    nnIdx,
                    dists,
                    0.01);
            int idx = nnIdx[0];
            delete [] nnIdx;
            delete [] dists;

            return m_mesh_faces[idx];
        }
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
