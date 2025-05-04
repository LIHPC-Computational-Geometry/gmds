/*----------------------------------------------------------------------------*/
/*
 * SingGraphBuilder2DOriginal.h
 *
 *  Created on: April 10 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#pragma once
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_SING_GRAPH_BUILD_export.h"
#include <gmds/singGraphBuild/SingularityGraph.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
#include <gmds/singGraphBuild/Tools.h>
/*----------------------------------------------------------------------------*/

class LIB_GMDS_SING_GRAPH_BUILD_API SingGraphBuilder2DOriginal : public SingularityGraphBuilder2D
{
 private:
	/** technical container, which is used to store all the free slots during the
	 *  algorithm */
	std::list<SingularityPoint::Slot *> m_free_slots;
 public:
	SingGraphBuilder2DOriginal(gmds::Mesh *AMesh, gmds::Variable<gmds::math::Cross2D> *AField,
	                               const bool ABuildGeomSing = true);

	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creation of singularity lines 
	*   formely SingularityGraphBuilder2D::createSingularityLines
	*/
	void createSingularityLines();
	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Compute the singularity line starting from point ASingPnt in the
	 *         direction of ASingSlot
	 *
	 * \param ASingPnt  the singularity point we start from
	 * \param ASingSlot the associated slot we really consider
	 * \param ARemovedSlot a slot that has been removed (if any)
	 */
	void computeSingularityLine(SingularityPoint *AFromPoint, SingularityPoint::Slot *AFromSlot, unsigned int &cont,
	                            SingularityPoint::Slot *&ARemovedSlot);

		/*----------------------------------------------------------------------------------------------------*/
	/** \brief Compute the singularity line starting from point \param AFromPnt in the
	 *         direction of \param AFromSlot
	 *
	 * \param[in] AFromPnt      the singularity point we start from
	 * \param[in] AFromSlot     the associated slot we consider
	 * \param[out] AToSingPnt    the singularity point we arrive at, if it exists
	 * \param[out] AToSlot       the associated slot
	 * \param[out] AToPnt        the location we arrive at the end
	 * \param[out] AToDir        the last direction to arrive at the end
	 * \param[out] APoints       a discretization of the stream line
	 * \param[out] ATriangles    the list of traversed triangles
	 * \param[out] AToCellDim    the dimension of the cell we finish on (0,1,2)
	 * \param[out] AToCellID     the id of the cell we finish on
	 * \param[out] streamlineDeviation     the deviation of the computed streamline
	 * \param[out] AEndOnBnd     indicates if we finish on the boundary (true)
	 * \param[out] AToSlotIsFree indicates if we finish onto a free slot (true)
	 * \param[out] APntToCreate  indicates if we must create the end point (true)
	 */
	/*----------------------------------------------------------------------------------------------------*/
	void computeStreamLine(SingularityPoint *AFromPnt,
	                       SingularityPoint::Slot *AFromSlot,
	                       SingularityPoint *&AToSingPnt,
	                       SingularityPoint::Slot *&AToSlot,
	                       gmds::math::Point &AToPnt,
	                       gmds::math::Vector3d &AToDir,
	                       std::vector<gmds::math::Point> &APoints,
	                       std::vector<gmds::TCellID> &ATriangles,
	                       int &AToCellDim,
	                       gmds::TCellID &AToCellID,
	                       double &streamlineDeviation,
	                       bool &AEndOnBnd,
	                       bool &AToSlotIsFree,
	                       bool &APntToCreate);

	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Improve a singularity line to avoid side effects due to the way
	 *         singularity points are considered.
	 *
	 * \param ALine     the singularity line we work on
	 * \param AFromPnt  the singularity point we start from
	 * \param AFromSlot the associated slot we really consider
	 * \param AToPnt    the singularity point we go to
	 * \param AToSlot   the associated slot we really consider
	 * \param to_dir    the direction of the singularity line
	 * \param to_dir    the deviation of the newly computed singularity line
	 * \param foundBackTrackPath   boolean value indicating wheather we have found a path
	 */
	void backtrackSingularityLine(SurfaceSingularityLine *ALine,
	                              SingularityPoint *AFromPnt,
	                              SingularityPoint::Slot *AFromSlot,
	                              SingularityPoint *AToPnt,
	                              SingularityPoint::Slot *AToSlot,
	                              gmds::math::Vector3d &to_dir,
	                              bool &foundBackTrackPath);
};