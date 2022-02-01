//
// Created by rochec on 01/02/2022.
//

#ifndef GMDS_POINTFOLLOWINGVECTORFIELDONNODES_H
#define GMDS_POINTFOLLOWINGVECTORFIELDONNODES_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API PointFollowingVectorFieldOnNodes
{
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
	 */
	PointFollowingVectorFieldOnNodes(Mesh *AMesh, math::Point A_Pstart, double A_distance, Variable<math::Vector3d>* A_gradient);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Return the closest node of M is
	 */
	TCellID closestNode(math::Point M);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the minimal edge's lenght
	 */
	double minEdgeLenght();
	/*-------------------------------------------------------------------*/
	/** @brief Write the discrete path in a vtk field
	 */
	void writeDiscretePathInVTK();
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** starting point */
	math::Point m_Pstart ;
	/** ending point */
	math::Point m_Pend ;
	/** the distance */
	double m_distance ;
	/** the vector field to follow */
	Variable<math::Vector3d>* m_gradient ;
	/** liste des points intermédiaires */
	std::vector<math::Vector3d> m_discrete_path;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_POINTFOLLOWINGVECTORFIELDONNODES_H
