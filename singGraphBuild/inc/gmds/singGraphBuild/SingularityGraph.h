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
 * SingularityGraph.h
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
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
	void splitSurfaceLine(SurfaceSingularityPoint*& APnt,
			      SurfaceSingularityLine*& ALine);
	

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

private:

	gmds::Mesh* m_mesh;
	std::vector<SingularityPatch*> m_patchs;
	std::vector<SingularityLine*> m_lines;
	std::vector<SingularityPoint*> m_points;

};
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYGRAPH_H_ */
/*----------------------------------------------------------------------------*/
