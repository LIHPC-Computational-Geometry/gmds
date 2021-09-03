/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*
 * SingularityGraphBuilder.h
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#ifndef SINGULARITYGRAPHBUILDER_H_
#define SINGULARITYGRAPHBUILDER_H_
/*----------------------------------------------------------------------------*/
#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Quaternion.h>
/*----------------------------------------------------------------------------*/
#include <gmds/singGraphBuild/SingularityGraph.h>
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
/*----------------------------------------------------------------------------*/
class LIB_GMDS_SINGGRAPHBUILD_API SingularityGraphBuilder
{
 public:
  //	static const int mask=DIM3|R|N|F|R2N|F2N|N2F|N2R|R2F|F2R;
  //static const int mask = DIM3 | N | E | F | R | F2N | N2F | E2N | N2E | F2E | E2F | R2N | N2R | R2E | E2R | R2F | F2R;

  /*------------------------------------------------------------------------*/
  /** \brief Constructor
   */
 SingularityGraphBuilder(gmds::Mesh* AMesh)
   :m_mesh(AMesh), m_output_directory_name(""), m_graph(AMesh)
  {
    distanceRaccord_ = 0.005;
  }

  void execute();

  void createVolumeSingularityPoints();

  void createVolumeSingularityLines();
	
  void createSurfaceSingularityLines();


  /*------------------------------------------------------------------------*/
  /** \brief Function to give the directory where we want to put the output
   *		   files (vtk files)
   */
  void setDebugDirectory(const std::string& ADirName)
  {
    m_output_directory_name = ADirName;
  }

private:
	
	void initConfusionBallOfASurfaceSingularityPoint(
		SurfaceSingularityPoint* ASingPoint);
	void extractQuaternions();

	void initMarks();
	void cleanMarks();

	void createOneSingularityLineFrom(
		SingularityPoint* ASingPoint,
		SingularityPoint::Slot* ASlot,
		const int AMarkFacesUsedForSep);

	void colorFaces(const int markF, const int markE);

	void computeVecInSurfacePoint(
		const gmds::Node& ANode, 
		const gmds::math::Vector3d& ANormal, 
		 gmds::math::Vector3d& AFromVec, 
		gmds::math::Vector3d& AToVec);

	void computeBoundSingInfo(
		const int AFaceMark,
		gmds::Face& AFace,
		std::vector<gmds::math::Vector3d >& AVecRep,
		gmds::math::Vector3d& ADirNormal,
		gmds::math::Point& APosSing,
		std::vector<gmds::math::Point  > & posSep,
		std::vector<gmds::math::Vector3d > & dirSep,
		std::vector<int> & sepInFaceID);

	void computeSingInfoAlongEdgeOfTriangle(
		const int AFaceMark,
		gmds::Face& AFace,
		gmds::Node& ANode1, 
		gmds::Node& ANode2,
		gmds::math::Vector3d& AV1,
		gmds::math::Vector3d& AV2,
		gmds::math::Vector3d& ADirNormal,
		const gmds::math::Point& AFaceSingularityPoint,
		std::vector<gmds::math::Point  > & posSep,
		std::vector<gmds::math::Vector3d > & dirSep,
		std::vector<int> & sepInFaceID);

	void computeSingInfoOnEdge(
		const int AFaceMark,
		gmds::Node& ANode1, 
		gmds::Node& ANode2,
		gmds::Face& AFace,
		double cross1[2][3],
		double cross2[2][3],
		gmds::math::Vector3d& ANormalToFace,
		double pos1[3],
		double pos2[3],
		double posSing[3],
		std::vector<gmds::math::Point  >& posSep,
		std::vector<gmds::math::Vector3d >& dirSep,
		std::vector<int> & sepInFaceID);

	void buildInfoAlongEdge(
		const int &markBF,
		gmds::Face& AFace,
		gmds::Node& ANode1,
		gmds::Node& ANode2,
		const double soleq,
		const double x1,
		const double y1,
		const double z1,
		const double x2,
		const double y2,
		const double v1,
		const double w1,
		const double v2,
		const double w2,
		const double x2MINx1,
		const double y2MINy1,
		const double z2MINz1,
		const double x1MINxsRotated,
		const double x2MINx1Rotated,
		const double y1MINysRotated,
		const double y2MINy1Rotated,
		std::vector<gmds::math::Point  >& posSep,
		std::vector<double>& dirSepXRotated,
		std::vector<double>& dirSepYRotated,
		std::vector<int> & sepInFaceID);

	void solveSecondOrderEq(
		const double a,
		const double b,
		const double c,
		const int &markBF,
		gmds::Face& AFace,
		gmds::Node& ANode1,
		gmds::Node& ANode2,
		const double x1,
		const double y1,
		const double z1,
		const double x2,
		const double y2,
		const double v1,
		const double w1,
		const double v2,
		const double w2,
		const double x2MINx1,
		const double y2MINy1,
		const double z2MINz1,
		const double x1MINxsRotated,
		const double x2MINx1Rotated,
		const double y1MINysRotated,
		const double y2MINy1Rotated,
		std::vector<gmds::math::Point  >& posSep,
		std::vector<double>& dirSepXRotated,
		std::vector<double>& dirSepYRotated,
		std::vector<int> & sepInFaceID);
	/*
	 * @param AFace a face
	 * @param ARegion a region adjacent to AFace
	 * @return the vector normal to AFace and going out from ARegion
	 */
	gmds::math::Vector3d getOutputNormal(gmds::Face& AFace, gmds::Region& ARegion);
	gmds::math::Vector3d getOutputNormalOfABoundaryFace(gmds::Face& AFace);
	gmds::math::Vector3d getInputNormal(gmds::Face& AFace, gmds::Region& ARegion);

	/*returns true if it can connect to an existing sing*/
	bool tryToconnectToExistingCurveSingularity(
		int markBordNode,
		int markBordEdge,
		SurfaceSingularityLine* ALine,
		gmds::Edge& AEdge,
		gmds::math::Point& APnt,
		gmds::math::Vector3d& AVec);
	/* Returs if APnt is on the boundary of AFace
	 */
	bool isOnBoundary(const gmds::math::Point& APnt, gmds::Face& AFace);
	bool isIn(const gmds::math::Point& APnt, gmds::Edge& AEdge);
	/* Returs if (APnt,AVec) is a vector pointing into  AFace
	*/
	bool isInside(const gmds::math::Point& APnt,
		const gmds::math::Vector3d& AVec, gmds::Face& AFace);

	int testTriangleSingularity(int ID1, int ID2, int ID3);

	
	// this methods is called to create a surface singularity point during the 
	// volume skeleton creation
	// \param ACreatedFromLine indicates this surface point is created at the
	//			end of volume sing. line
	SingularityPoint* createSurfaceSingularityPointFromVolumeInfo(
		gmds::Face& FCurr,
		const bool ACreatedFromVolumeLine =true);
	void readFile();

	void vec2Triad(gmds::math::Quaternion& AQ,
		double& u1, double& u2, double& u3,
		double& v1, double& v2, double& v3,
		double& w1, double& w2, double& w3);


	void addGeometryToSingularitGraph();

	void writeOutput(const std::string& AFileName);
	void writeOutputSingle(const std::string& AFileName);

	bool intersecDirectionWithSegment(
		double P1[3], double P2[3],
		double PStart[3],
		double dirStart[3],
		double& coeff,
		double& coeffBis,
		double& PEnd0, double& PEnd1, double& PEnd2,
		double epsilon);
	void InvertMatrix(double matrixInit[3][3], double& matrixInv00, double& matrixInv01, double& matrixInv02,
		double& matrixInv10, double& matrixInv11, double& matrixInv12,
		double& matrixInv20, double& matrixInv21, double& matrixInv22);

	void createBoundarySingularityGraphFromField();

	//compute a sep line on the surface
	void computeSurfaceSingularityLine(SingularityPoint* AFromPoint, SingularityPoint::Slot* AFromSlot);


	void computeShortSepLine(
		const gmds::math::Point&  startPnt,
		const gmds::math::Vector3d& startDir,
		gmds::Face& startFace, //next face
		const gmds::math::Point&   singPnt,
		gmds::Face& singFace,
		std::vector<gmds::math::Point >& points,
		int sepNumber);

	void createOneVolumeSingularityPoint(std::vector<gmds::Region>& ACluster);

	gmds::math::Quaternion computeLinearQuaternion(
		gmds::math::Point& APnt,
		gmds::Edge& AEdge);

	gmds::math::Quaternion buildLinearQuaternion(
		gmds::Edge& AEdge, const double AParam);

	void computeOutTriangleOnTheSurface(
		gmds::Face& ACurrentFace,
		gmds::math::Point&  AINPnt,
		gmds::math::Vector3d& AInVec,
		gmds::Face& AOUTFace,
		gmds::math::Point&  AOUTPnt,
		gmds::math::Vector3d& AOUTVec);

	/* \brief compute the intersection of a ray (AINPnt, AINVec) with a node 
	 *		  and two edges. The out point and direction are returned as the
	 *		  the index of the intersected element (0 node, 1 first edge, 2
			  second edge). The out Quaternion is computed along the underlying
			  frame field.
	 */
	int computeOutDataOnTheSurface(
		const gmds::math::Vector3d& ANormal,
		gmds::math::Point&  AINPnt,
		gmds::math::Vector3d& AINVec,
		gmds::Node& AINNode1,
		gmds::Node& AINNode2,
		gmds::Node& AINNode3,
		gmds::Edge& AINEdge1,
		gmds::Edge& AINEdge2,
		gmds::math::Point&  AOUTPnt,
		gmds::math::Vector3d&  AOUTVec);


	/* \brief compute the intersection of a ray (AINPnt, AINVec) with an edge
	 *			AEdge. The intersection is returned in AOUTPnt. IF AOUTParam is
	 *		 included in ]0,1[, the edge is cut. Otherwise, no
	 *	\return true if the edge is cut by (AINPnt, AINVec)
	 */
	bool computeIntersection(
		gmds::math::Point&  APnt,
		gmds::math::Vector3d& AVec,
		gmds::Edge AEdge,
		gmds::math::Point&  AOUTPnt,
		double& AOUTParam,
		double& AOUTParamRay);

	void detectSingularTetrahedra();

	void initSingularityPointOnSurf();

	/* Check if a face is traversed by a singularity line
	*/
	bool isSingularFace(gmds::Face& AF);

	/*	void createSkeleton();
		void createSkeletonNew();
		void writeSkeleton();

		BoundSepInfo getBoundSepInfo(Separatrix3D& sep, int markBordFace);
		void matchCurves( Separatrix3D& s1, Separatrix3D& s2,
		std::vector<Singularity3D>& singularities,
		std::vector<Separatrix3D>& separatrices,
		int markFace, int markEdge, int markNode, int markbordSing,
		int markEdgeEnd);

		double computeBoundingBoxDiag();

		Separatrix3D spreadCurve(Separatrix3D& sep,
		std::vector<Singularity3D>& singularities,
		std::vector<Separatrix3D>& separatrices,
		int markFace, int markEdge, int markNode,
		int markbordSing, int markEndEdge);

		bool getMatchingBoundaryFreeSep(
		Separatrix3D& sep, std::vector<Separatrix3D>& otherSeps,
		int& matchingSepNumber,
		int markBordFace, int markBordEdge, int markBordNode,
		int markEdgeSep, std::vector<Singularity3D>& singularities);

		Edge getEdge(Node n0, Node n1);

		void computeCurveDistances( Edge e, std::vector<Edge>& otherEdges,
		math::Point<double,3> pnt, std::vector<math::Point<double,3> >& otherPoints,
		std::vector<int>& topoDistance, std::vector<double>& geoDistance,
		int markBordFace, int markBordEdge, int markBordNode);

		void getOrderedEdgesOnACurve(Edge e, Node startingNode,
		std::vector<Edge>& curveEdges,
		int markBordFace, int markBordEdge, int markBordNode);


		int testTriangleSingularity(int ID1, int ID2, int ID3, int ID4);
		bool IntersecDirectionWithSegment(double P1[3], double P2[3], double PStart[3], double dirStart[3], double& coeff, double& coeffBis, double& PEnd0, double& PEnd1, double& PEnd2, double epsilon);
		void ComputeSepLine(int markBordEdge, int markbordFace, int markbordSing, double xStart, double yStart, double zStart,
		double dirXStart, double dirYStart, double dirZStart, int FaceStartID,
		int& nodeEndID1, int& nodeEndID2, std::vector<double>& listePointX,
		std::vector<double>& listePointY, std::vector<double>& listePointZ,
		std::vector<Singularity3D>& singularities,
		std::vector<Separatrix3D>& separatrices,
		int sepNumber, double xSing, double ySing,double zSing,
		bool isInFieldLinePropagation=true,
		bool isFromANode = false);


		void ComputeShortSepLine(
		int markBordEdge, int markbordFace, int markbordSing,
		double xStart, double yStart, double zStart,
		//point de depart de la sepline
		double dirXStart, double dirYStart, double dirZStart,
		double xSing, double ySing,double zSing,
		// xSing, ySing et zSing sont les coordonnées de la singularité d'ou l'on part,
		int FaceStartID, int FaceNextID,
		std::vector<double>& listePointX,
		std::vector<double>& listePointY,
		std::vector<double>& listePointZ,
		std::vector<Singularity3D>& singularities,std::vector<Separatrix3D>& separatrices,
		int sepNumber);
		void computePrimal(int markBordFace, int markBordEdge, int markBordNode, int markAloneNodes);
		bool TestDirectionInTriangle(double dirInit[3], int IDFace, int IDNode);

		void InvertMatrix(double matrixInit[3][3], double& matrixInv00, double& matrixInv01, double& matrixInv02,
		double& matrixInv10, double& matrixInv11, double& matrixInv12,
		double& matrixInv20, double& matrixInv21, double& matrixInv22);

		void triad2Vec(double u1, double u2, double u3,
		double v1, double v2, double v3,
		double w1, double w2, double w3,
		double& q1, double& q2, double& q3, double&q4);

		void triadAssociation(int ID1, int ID2, int& assoc1, int& assoc2, int& assoc3);

		void computeSingLines(int markBordFace, int markBordEdge, int markBordNode, int markAloneNodes);



		double computeDistance(double X, double Y, double Z);

		void buildMeshSkeleton(std::vector<Separatrix3D>& separatrices);
		*/


	std::vector<gmds::Face> getAdjacentFacesByNodes(gmds::Face& AFace, const int AMark);
private: 
	gmds::Mesh* m_mesh;
	gmds::Mesh* meshSkel;
	double distanceRaccord_;

	//directory where debug files will be written
	std::string m_output_directory_name;

	gmds::Variable<gmds::math::Quaternion>* m_var_quatern;
	gmds::Variable<double>* m_var_distance;
	gmds::Variable<int>* m_var_color;
	SingularityGraph m_graph;

	std::list<gmds::TCellID> m_3SingTetIDs;
	std::list<gmds::TCellID> m_2SingTetIDsAlongSurf;
	std::map<gmds::TCellID, int> m_region_singularity_type;
	std::map<gmds::TCellID, SurfaceSingularityPoint*> m_face_to_singularity_on_surf;

	//boolean marks for executing the algorithm

	int m_markNodesOnVert;
	int m_markNodesOnCurv;
	int m_markNodesOnSurf;
	int m_markAloneNodes;

	int m_markEdgesOnSurf;
	int m_markEdgesOnCurv;

	int m_markFacesOnSurf;

	int m_markClusterSingDone;
	int m_markBordVolSingForFaces;
	int m_markBordVolSingForEdges;
	int m_markBordVolSingForNodes;

	int m_markBordSurfSingForFaces;//faces containing sing. are marked with it
	int m_markBordSurfSingForEdges;//edges containing sing. are marked with it
	int m_markBordSurfSingForNodes;//nodes containing sing. are marked with it

};
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYGRAPHBUILDER_H_ */
/*----------------------------------------------------------------------------*/
