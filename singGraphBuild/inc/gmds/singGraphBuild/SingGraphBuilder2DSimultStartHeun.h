/*----------------------------------------------------------------------------*/
/*
 * SingGraphBuilder2DSimultStartHeun.h
 *
 *  Created on: April 10 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#pragma once
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_SINGGRAPHBUILD_export.h"
#include <gmds/singGraphBuild/SingularityGraph.h>
#include <gmds/singGraphBuild/SingularityGraphBuilder2D.h>
#include <gmds/singGraphBuild/Tools.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

class LIB_GMDS_SINGGRAPHBUILD_API SingGraphBuilder2DSimultStartHeun : public SingularityGraphBuilder2D
{
 public:
	SingGraphBuilder2DSimultStartHeun(Mesh *AMesh, Variable<math::Cross2D> *AField, const bool ABuildGeomSing = true);

	/*----------------------------------------------------------------------------------------------------*/
	/** \brief Creation of singularity lines starting simultaneously from all singularities
	 *   (using Heun's computation)
	 *   formely SingularityGraphBuilder2D::createSingularityLinesSimultaneousStart
	 */
	void createSingularityLines();
};
}
