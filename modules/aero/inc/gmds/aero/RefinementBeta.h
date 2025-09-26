//
// Created by rochec on 27/07/2022.
//

#ifndef GMDS_REFINEMENTBETA_H
#define GMDS_REFINEMENTBETA_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include "gmds/ig/Mesh.h"
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSAero_API RefinementBeta
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
         *  @param[in] APoints point vector
         *  @param[in] Asize_first_edge size of the first edge
         *
	 */
	RefinementBeta(std::vector<math::Point>& APoints, double Asize_first_edge);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief Get the new vector of positions
		*
		* \return
	 */
	std::vector<math::Point> GetNewPositions();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the Beta parameter
	 	* \param[in] first_edge_size the size imposed on the first edge
	 	* \param[in] sum_edge_sizes the sum of the edges size
	 	* \param[in] Nbr_Points final number of points
		*
		* \return  a double Beta computed by a Newton algorithm
	 */
	static double ComputeBeta(double first_edge_size, double sum_edge_sizes, int Nbr_Points);
	/*-------------------------------------------------------------------*/


 private:
	/** Initial point vector */
	std::vector<math::Point> m_Points;
	/** Refined vector of points */
	std::vector<math::Point> m_PointsRefined;
	/** Size first edge */
	double m_size_first_edge;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_REFINEMENTBETA_H
