/*----------------------------------------------------------------------------*/
/*
 * SingularityGraph.h
 *
 *  Created on: 13 juil. 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef SINGULARITYGRAPH_H_
#define SINGULARITYGRAPH_H_
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/singGraphBuild/SingularityPatch.h>
#include <gmds/singGraphBuild/SingularityPoint.h>
#include <gmds/singGraphBuild/SingularityLine.h>
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
/*----------------------------------------------------------------------------*/
class LIB_GMDS_SINGGRAPHBUILD_API SingularityGraph
{
public:
	SingularityGraph(gmds::Mesh* AMesh);
	virtual ~SingularityGraph();

	SingularityPatch* 		newSurfacePatch();

	CurveSingularityLine* 		newCurveLine();
	SurfaceSingularityLine* 	newSurfaceLine();
	VolumeSingularityLine* 		newVolumeLine();

	VertexSingularityPoint*		newVertexPoint();
	CurveSingularityPoint*		newCurvePoint();
	SurfaceSingularityPoint* 	newSurfacePoint();
	VolumeSingularityPoint* 	newVolumePoint();

	std::vector<SingularityPoint*>& getPoints();
	std::vector<SingularityLine* >& getLines();

	std::vector<VertexSingularityPoint* >  getVertexPoints();
	std::vector<CurveSingularityPoint* >   getCurvePoints();
	std::vector<SurfaceSingularityPoint* > getSurfacePoints();
	std::vector<VolumeSingularityPoint* >  getVolumePoints();

	std::vector<CurveSingularityLine* >   getCurveLines();
	std::vector<SurfaceSingularityLine* > getSurfaceLines();
	std::vector<VolumeSingularityLine* >  getVolumeLines();

	std::vector<SingularityPatch* >  getSurfacePatchs();

	int getNbPoints() const;
	int getNbVolumePoints() const;
	int getNbSurfacePoints() const;
	int getNbCurvePoints() const;
	int getNbVertexPoints() const;

	int getNbLines() const;
	int getNbVolumeLines() const;
	int getNbSurfaceLines() const;
	int getNbCurveLines() const;

	int getNbSurfacePatches() const;
	/*------------------------------------------------------------------------*/
	/* \brief This operation tries to put a singularity point onto existing 
	 *        geometrical curves.
	 *
	 * \param APnt the point the singularity is located
	 * \param AVec the direction the splitting comes from
	 * \param AEdge the edge the split must occur (it contains APnt)
	 * \param [OUT] APnt the created singularity point
	 * \param [OUT] ASlot the associated slot directed towards AInVec
	 */
	void splitCurveLine(const gmds::math::Point&  APnt, 
			    const gmds::math::Vector3d& AVec, 
			    const gmds::Edge&         AEdge,
			    SingularityPoint*&         ASing,
			    SingularityPoint::Slot*&   ASlot);


	/*------------------------------------------------------------------------*/
	/* \brief This operation tries to put a singularity point onto existing
	 *        geometrical curves.
	 *
	 * \param APnt the point the singularity is located
	 * \param AVec the direction the splitting comes from
	 * \param ANode the node the split must occur (it must be located at APnt)
	 * \param [OUT] APnt the created singularity point
	 * \param [OUT] ASlot the associated slot directed towards AInVec
	 */
	void splitCurveLine(const gmds::math::Point&  APnt, 
			    const gmds::math::Vector3d& AVec, 
			    const gmds::Node&         ANode,
			    SingularityPoint*&         ASing,
			    SingularityPoint::Slot*&   ASlot);
	
	/*------------------------------------------------------------------------*/
	/* \brief This operation splits a surface sing. line by inserting a sing
	 *        point
	 *
	 * \param APnt the singularity point
	 * \param ALine the singularity line
	 */
	void splitSurfaceLine(SurfaceSingularityPoint* APnt,
			      SurfaceSingularityLine* ALine);
	
	/*------------------------------------------------------------------------*/
	/* \brief This operation splits a surface sing. line by inserting a sing
	 *        point. Unlike splitSurfaceLine, it dooes not recompute traversed tri.       
	 *
	 * \param APnt the singularity point
	 * \param ALine the singularity line
	 * \param APrevPointID the ID of the point within the line discretization located just before the singPoint
	 */
	void splitSurfaceLineSimple(SurfaceSingularityPoint *APnt, SurfaceSingularityLine *ALine, unsigned int APrevPointID);
	/*------------------------------------------------------------------------*/
	/* \brief Removes the line ALine from the singularity graph. Connected singularity
	 *        points are updated too.
	 *
	 * \param ALine a singularity line of the graph
	 */
	void removeLine(SingularityLine* ALine);
	
	
	/*------------------------------------------------------------------------*/
	/* \brief Removes the singularity point APoint from the singularity graph. Connected singularity
	 *        points are updated too. it also removed the singularity lines connected to this point.
	 *
	 * \param APoint a singularity point of the graph
	 */
	void removePoint(SingularityPoint* APoint);
	
	/*------------------------------------------------------------------------*/
	/* \brief Removes the curve line ALine from the singularity graph. Updated singularity graph Line Numbers
	 *
	 * \param ALine a singularity curve line of the graph
	 */
	void removeCurveLine(SingularityLine* ALine);
	
	/*------------------------------------------------------------------------*/
	/* \brief Build the surface patchs from the set of lines and
	 *        points. Warning, this function assumes that lines and points
	 *        are coherent. No check is done.
	 */
	void buildSurfacePatchs();
        
        /*------------------------------------------------------------------------*/
	/* \brief Build the surface patchs from the set of lines and
	 *        points. Warning, this function assumes that lines and points
	 *        are coherent. No check is done. The surface patches will be curve.
	 */
	void buildCurveSurfacePatchs(unsigned int& number_of_control_points);
              
     double BernsteinBasis(unsigned int j, unsigned int k, double n);
     
	double binomialCoeff(unsigned int& n, unsigned int& k);
	
	gmds::math::Point interpolateBtwPoints(gmds::math::Point& point1, gmds::math::Point& point2, double& interpFact);
            
	void createVTKOutputFile(const std::string& AFileName) const;
        
     void createVTKOutputFile(const std::string& AFileName, bool &curvePatches) const;

	 // compute group of linked edges (singularity lines)
	 void computeLinkedEdges();

	 size_t getNLinkedEdgeGroup() const { return m_nGroupLinkedEdges; };
private:

	gmds::Mesh* m_mesh;
	std::vector<SingularityPatch*> m_patchs;
	std::vector<SingularityLine*> m_lines;
	std::vector<SingularityPoint*> m_points;
	size_t m_nGroupLinkedEdges = 0;

};
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYGRAPH_H_ */
/*----------------------------------------------------------------------------*/
