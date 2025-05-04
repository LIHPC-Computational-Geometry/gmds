/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSSmoothy_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ABSTRACT_SMOOTHER_H
#	define GMDS_ABSTRACT_SMOOTHER_H
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace smoothy {
/*----------------------------------------------------------------------------*/
/** \class AbstractSmoother
 *  \brief This class provides attributes and methods that are shared by all
 *         smoothing algorithms as child classes. It is abstract and so cannot
 *         be instanciated.
 */
/*----------------------------------------------------------------------------*/
class GMDSSmoothy_API AbstractSmoother
{
 public:
	/**@brief constructor
	 * @param ALinker the linker that has the knowledge and connection
	 *                between geometry and mesh
	 */
	AbstractSmoother(Mesh *AMesh);

	/**@brief validation of the mesh model
	 * @return true if the model fits the algorithms requirement. False otherwise
	 */
	virtual bool isValid() const = 0;
	/** Executes the smoothing algorithm
	 * @return 1 if it succeeds, 0 if it fails. Other values could be returned by
	 * specific implementations.
	 */
	virtual int smooth() = 0;
	/**@brief specify the nodes to smooth
	 * @param ANodeIds the list of node ids to smooth
	 */
	void setNodes(std::vector<TCellID> &ANodeIds);
	/**@brief to smooth the whole mesh. In this case, the algorithm could be written on
	 * the full mesh, but each algorithm can decide what to do with constrained nodes on the
	 * boundary for instance.
	 */
	void setNodesAll();
	/**@brief Set the number of smoothing iterations
	 * @param ANbIterations the expected number of iterations
	 */
	void setNbIterations(const int ANbIterations);

 protected:
	/** the mesh we work on*/
	Mesh *m_mesh;
	/** the number of smoothing iterations*/
	int m_nb_iterations;
	/** the nodes we wants to smooth */
	std::vector<TCellID> m_nodes;
};
/*----------------------------------------------------------------------------*/
}     // namespace smoothy
/*----------------------------------------------------------------------------*/
}     // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_ABSTRACT_SMOOTHER_H
/*----------------------------------------------------------------------------*/
