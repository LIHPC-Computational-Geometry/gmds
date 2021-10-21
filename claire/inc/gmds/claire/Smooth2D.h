#ifndef GMDS_SMOOTH2D_H
#define GMDS_SMOOTH2D_H

#include <gmds/ig/Mesh.h>
#include <string>
#include <map>
namespace gmds {
class Smooth2D
{
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
	 *  @param AMesh the mesh where we work on
	 *  @param AVarBnd node variable (value 1 means constrained node)
	 *  @param ANbIterations nb max iterations
	 */
	Smooth2D(Mesh *AMesh,
	         const Variable<int>* AVarBnd,
	         const int ANbIterations = 100);

	/*------------------------------------------------------------------------*/
	/** @brief Set the max number of iterations
	 *  @param[in] ANbIterations
	 */
	void setNbIterations(const int ANbIterations);

	/*------------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** Variable to store which nodes are constrained*/
	const Variable<int>* m_node_constrained;
	/** nb max iterations */
	int m_nb_max_iterations;

	typedef struct{
		unsigned int val[3][3];
	}stencil;
	std::map<TCellID,stencil> m_node_to_stencil;
};

}     // namespace gmds

#endif     // GMDS_SMOOTH2D_H
