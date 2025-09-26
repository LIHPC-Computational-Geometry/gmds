//
// Created by aero on 26/11/2021.
//

#ifndef GMDS_GRID_SMOOTH2D_H
#define GMDS_GRID_SMOOTH2D_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Blocking2D.h>
/*----------------------------------------------------------------------------*/
#include <string>
#include <map>
/*----------------------------------------------------------------------------*/
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSAero_API Grid_Smooth2D {
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param AMesh the mesh where we work on
         *  @param AVarBnd node variable (value 1 means constrained node)
         *  @param ANbIterations nb max iterations
	 */
	Grid_Smooth2D(Blocking2D *AMesh,
	         const int ANbIterations = 100);

	/*-------------------------------------------------------------------*/
	/** @brief Set the max number of iterations
         *  @param[in] ANbIterations
	 */
	void setNbIterations(const int ANbIterations);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
 private:
	/*-------------------------------------------------------------------*/
	/** @brief Find the middle of a branch between 3 points
	 */
	math::Point FindMidBranche(const math::Point A, const math::Point B, const math::Point C);
 private:
	/** mesh we work on */
	Blocking2D *m_mesh;
	/** nb max iterations */
	int m_nb_max_iterations;

	/** free nodes in the mesh */
	std::vector<TCellID> m_free_nodes;
	typedef struct {
		unsigned int val[3][3];
	} stencil;
    /** stencil for each node id */
	std::map<TCellID, stencil> m_stencil;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_GRID_SMOOTH2D_H
