/*----------------------------------------------------------------------------*/
/** \file    FacetedSurface.t.h
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#include <algorithm>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
// ANN File Headers
#include "ANN/ANN.h"
/*----------------------------------------------------------------------------*/
#include "gmds/cadfac/FACSurface.h"
/*----------------------------------------------------------------------------*/
#include "gmds/cadfac/FACCurve.h"
#include "gmds/cadfac/FACPoint.h"
#include "gmds/cadfac/FACVolume.h"
#include "gmds/math/Ray.h"
#include "gmds/math/Triangle.h"
#include "gmds/utils/RandomGenerator.h"
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace cad {
/*----------------------------------------------------------------------------*/
int FACSurface::m_next_id = 1;
/*----------------------------------------------------------------------------*/
void
FACSurface::resetIdCounter()
{
	m_next_id = 1;
}
/*----------------------------------------------------------------------------*/
FACSurface::FACSurface(Mesh *AMesh) : m_support(AMesh), m_id(m_next_id++), m_kd_tree(nullptr) {}
/*----------------------------------------------------------------------------*/

FACSurface::FACSurface(Mesh *AMesh, std::vector<TCellID> &ADiscret, const std::string &AName) :
  GeomSurface(AName), m_support(AMesh), m_id(m_next_id++), m_mesh_faces(ADiscret), m_kd_tree(nullptr)
{
	buildANNTree();
}
/*----------------------------------------------------------------------------*/

FACSurface::~FACSurface()
{
	// got to clean some technical attributes used by ANN
	delete m_kd_tree;
	annDeallocPts(m_dataPts);
	annClose();     // done with ANN
}
/*----------------------------------------------------------------------------*/
TCoord
FACSurface::computeArea() const
{
	TCoord totalArea = 0.;

	for (unsigned int i = 0; i < m_mesh_faces.size(); i++) {

		totalArea += m_support->get<Face>(m_mesh_faces[i]).area();
	}
	return totalArea;
}

/*----------------------------------------------------------------------------*/
std::tuple<TCoord, TCoord, TCoord, TCoord, TCoord, TCoord>
FACSurface::BBox() const
{
	std::set<TCellID> node_ids;
	for (unsigned int i = 0; i < m_mesh_faces.size(); i++) {

		std::vector<TCellID> f_ids = m_support->get<Face>(m_mesh_faces[i]).getIDs<Node>();
		node_ids.insert(f_ids.begin(), f_ids.end());
	}

	// too many comparisons (3 factor)
	math::Point pi = m_support->get<Node>(*node_ids.begin()).point();
	TCoord  minXYZ[3], maxXYZ[3];
	minXYZ[0] = pi.X();
	maxXYZ[0] = pi.X();
	minXYZ[1] = pi.Y();
	maxXYZ[1] = pi.Y();
	minXYZ[2] = pi.Z();
	maxXYZ[2] = pi.Z();
	for (auto i : node_ids) {
		math::Point pi = m_support->get<Node>(i).point();
		minXYZ[0] = std::min(minXYZ[0], pi.X());
		minXYZ[1] = std::min(minXYZ[1], pi.Y());
		minXYZ[2] = std::min(minXYZ[2], pi.Z());
		maxXYZ[0] = std::max(maxXYZ[0], pi.X());
		maxXYZ[1] = std::max(maxXYZ[1], pi.Y());
		maxXYZ[2] = std::max(maxXYZ[2], pi.Z());
	}
	return {minXYZ[0], minXYZ[1], minXYZ[2], maxXYZ[0], maxXYZ[1], maxXYZ[2]};
}
/*----------------------------------------------------------------------------*/
void
FACSurface::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	auto [min_x, min_y, min_z,max_x,max_y,max_z] = BBox();
	minXYZ[0] = min_x;
	maxXYZ[0] = max_x;
	minXYZ[1] = min_y;
	maxXYZ[1] = max_y;
	minXYZ[2] = min_z;
	maxXYZ[2] = max_z;
}
/*----------------------------------------------------------------------------*/

void
FACSurface::computeNormal(const math::Point &AP, math::Vector3d &AV) const
{
	Face f0 = m_support->get<Face>(m_mesh_faces[0]);
	TCoord min_dist = f0.distance(AP);
	TCellID index = 0;
	for (unsigned int i = 1; i < m_mesh_faces.size(); i++) {
		TCoord dist = m_support->get<Face>(m_mesh_faces[i]).distance(AP);
		if (dist < min_dist) {
			min_dist = dist;
			index = i;
		}
	}
	AV = m_support->get<Face>(m_mesh_faces[index]).normal();
}
/*----------------------------------------------------------------------------*/
math::Point
FACSurface::closestPoint(const math::Point &AP) const
{
	Variable<int> *var_surf = m_support->getVariable<int, GMDS_FACE>("on_surface");

	Face seed = m_support->get<Face>(getANNClosestTriangle(AP));
	std::set<TCellID> face_ids;
	std::vector<Node> ns = seed.get<Node>();
	for (auto n : ns) {
		std::vector<TCellID> f_ids = n.getIDs<Face>();
		for (auto i : f_ids)
			if (var_surf->value(i) == this->id()) face_ids.insert(i);
	}

	Face f0 = m_support->get<Face>(*face_ids.begin());
	TCoord min_dist = f0.distance(AP);
	TCellID closest_face_id = f0.id();
	for (auto f_id : face_ids) {
		TCoord dist = m_support->get<Face>(f_id).distance(AP);
		if (dist < min_dist) {
			min_dist = dist;
			closest_face_id = f_id;
		}
	}
	return m_support->get<Face>(closest_face_id).project(AP);
}
/*----------------------------------------------------------------------------*/
math::Point
FACSurface::center() const
{
	throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void
FACSurface::setMeshFaces(const std::vector<Face> &AFaces)
{
	m_mesh_faces.clear();
	m_mesh_faces.resize(AFaces.size());
	for (unsigned int i = 0; i < AFaces.size(); i++) {
		m_mesh_faces[i] = AFaces[i].id();
	}

	if (m_kd_tree != nullptr) {
		delete m_kd_tree;
		annDeallocPts(m_dataPts);
	}

	buildANNTree();
}
/*----------------------------------------------------------------------------*/
void
FACSurface::getMeshFaces(std::vector<Face> &AFaces) const
{
	AFaces.clear();
	AFaces.reserve(m_mesh_faces.size());
	for (auto f_id : m_mesh_faces) {
		AFaces.push_back(m_support->get<Face>(f_id));
	}
}
/*----------------------------------------------------------------------------*/
void
FACSurface::getTriangulation(std::vector<math::Triangle> &ATri) const
{
	ATri.clear();
	ATri.resize(m_mesh_faces.size());
	for (unsigned int i = 0; i < m_mesh_faces.size(); i++) {
		gmds::Face current = m_support->get<gmds::Face>(m_mesh_faces[i]);
		std::vector<Node> nodes = current.get<Node>();
		ATri[i] = math::Triangle(nodes[0].point(), nodes[1].point(), nodes[2].point());
	}
}
/*----------------------------------------------------------------------------*/
void
FACSurface::propagateOrient(Face AFace, int AMarkTreatedFaces, int AMarkFacesInverted, Mesh *AMesh)
{
	std::vector<Face> faces = AFace.get<Face>();

	for (unsigned int iFace = 0; iFace < faces.size(); iFace++) {
		if (!AMesh->isMarked(faces[iFace], AMarkTreatedFaces)) {
			bool faceIsSameOrientation = checkSameOrientFace(AFace, faces[iFace]);

			if (!faceIsSameOrientation) {
				invertFace(faces[iFace]);
				AMesh->mark(faces[iFace], AMarkFacesInverted);
			}

			AMesh->mark(faces[iFace], AMarkTreatedFaces);
			propagateOrient(faces[iFace], AMarkTreatedFaces, AMarkFacesInverted, AMesh);
		}
	}
}
/*----------------------------------------------------------------------------*/
bool
FACSurface::checkSameOrientFace(Face AFaceRef, Face AFaceCheck)
{
	std::vector<Node> nodes = AFaceRef.get<Node>();
	std::vector<Node> nodesBis = AFaceCheck.get<Node>();

	for (unsigned int iNode = 0; iNode < nodes.size(); iNode++) {
		for (unsigned int iNodeBis = 0; iNodeBis < nodesBis.size(); iNodeBis++) {
			if (nodesBis[iNodeBis] == nodes[iNode]) {
				if (nodesBis[(iNodeBis + 1) % nodesBis.size()] == nodes[(iNode + 1) % nodes.size()]) {
					return true;
				}
				else {
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
	for (unsigned int iFace = 0; iFace < faces.size(); iFace++) {
		invertFace(faces[iFace]);
	}
}
/*----------------------------------------------------------------------------*/
void
FACSurface::invertFace(Face AFace)
{
	std::vector<Node> nodes = AFace.get<Node>();
	std::vector<Node> nodestmp(nodes.rbegin(), nodes.rend());

	AFace.set<Node>(nodestmp);
}
/*----------------------------------------------------------------------------*/
std::vector<GeomPoint *> &
FACSurface::points()
{
	return m_adjacent_points;
}
/*----------------------------------------------------------------------------*/
std::vector<GeomCurve *> &
FACSurface::curves()
{
	return m_adjacent_curves;
}
/*----------------------------------------------------------------------------*/
std::vector<std::vector<GeomCurve *> >
FACSurface::loops(){
  std::vector<std::vector<GeomCurve *> > loops;
  // first, we put a boolean flag on each curve to indicate we do not have
  // traversed it
  std::map<int,bool> traversed_curves;
  for(auto c:m_adjacent_curves)
		traversed_curves[c->id()]=false;

  bool remain_curves = true;
  while(remain_curves){
		remain_curves=false;
		//we pick a non-traversed curve
		int curve_id=-1;
		for (auto i=0;i<m_adjacent_curves.size() & remain_curves==false;i++){
			if(traversed_curves[m_adjacent_curves[i]->id()]==false){
				remain_curves=true;
				curve_id=i;
			}
		}
		if(remain_curves){
			//means, we start a new loop
			std::vector<GeomCurve*> cur_loop;
			GeomCurve* seed = m_adjacent_curves[curve_id];
			cur_loop.push_back(seed);
			traversed_curves[seed->id()]=true;
			std::vector<GeomPoint*> front;
			front.insert(front.end(),seed->points().begin(),seed->points().end());
			if(front.size()>1){
				// the current loop has two end points and we are going to access to adjacent curves
				// in the boundary surface
                while(!front.empty()){
                    auto pnt = front.back();front.pop_back();
                    std::vector<GeomCurve*> curves_of_pnt = pnt->curves();
                    for(auto c: curves_of_pnt){
                        //we get curves that are adjacent to the current surface only
                        std::vector<GeomSurface*> surfaces_of_c = c->surfaces();
                        bool find_current_surface = false;
                        for(auto s:surfaces_of_c){
                            if (s->id()== this->id())
                                find_current_surface=true;
                        }
                        if(find_current_surface){
                            //this curve bounds the current surface and share a point with the previous curve
                            if(!traversed_curves[c->id()]) {
                                cur_loop.push_back(c);
                                traversed_curves[c->id()]= true;
                                std::vector<GeomPoint *> c_pnts = c->points();
                                for (auto c_p: c_pnts)
                                    if (c_p->id() != pnt->id())
                                        front.push_back(c_p);
                            }
                        }
                    }
                }
			}
            //means the current curve form a loop itself
            loops.push_back(cur_loop);
		}
  }
  return loops;
}
/*----------------------------------------------------------------------------*/
std::vector<GeomVolume *> &
FACSurface::volumes()
{
	return m_adjacent_volumes;
}
/*---------------------------------------------------------------------------*/
void
FACSurface::buildANNTree()
{
	if (m_kd_tree != nullptr) throw GMDSException("FACSurface Issue: the kd tree structure was previously initialized");

	int dim = 3;                            // dimension
	int maxPts = m_mesh_faces.size();       // maximum number of data points
	m_dataPts = annAllocPts(maxPts, dim);     // allocate data points

	//========================================================
	// (1) Fill in the  ANN structure for storing points
	// Important: Points in APnts and dataPnts are stored with
	// same index.
	//========================================================
	int nPts = 0;     // actual number of data points

	while (nPts < maxPts) {
		math::Point p = m_support->get<Face>(m_mesh_faces[nPts]).center();
		m_dataPts[nPts][0] = p.X();
		m_dataPts[nPts][1] = p.Y();
		m_dataPts[nPts][2] = p.Z();
		nPts++;
	};
	//========================================================
	// (2) Build the search structure
	//========================================================
	m_kd_tree = new ANNkd_tree(m_dataPts,     // the data points
	                           nPts,        // number of points
	                           dim);        // dimension of space
}
/*---------------------------------------------------------------------------*/
TCellID
FACSurface::getANNClosestTriangle(const math::Point &AP) const
{
	int k = 1;                     // max number of nearest neighbors
	int dim = 3;                   // dimension

	ANNpoint queryPt;              // query point
	ANNidxArray nnIdx;             // near neighbor indices
	ANNdistArray dists;            // near neighbor distances
	queryPt = annAllocPt(dim);     // allocate 1 query point
	nnIdx = new ANNidx[k];         // allocate near neigh indices
	dists = new ANNdist[k];        // allocate near neighbor dists

	queryPt[0] = AP.X();
	queryPt[1] = AP.Y();
	queryPt[2] = AP.Z();

	m_kd_tree->annkSearch(     // search
	   queryPt,                // query point
	   k, nnIdx, dists, 0.01);
	int idx = nnIdx[0];

	annDeallocPt(queryPt);
	delete[] nnIdx;
	delete[] dists;

	return m_mesh_faces[idx];
}
/*----------------------------------------------------------------------------*/
}     // namespace cad
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
