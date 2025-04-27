/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSALGORITHM_H
#define GMDS_MCTSALGORITHM_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSTree.h>
#include <gmds/rlBlocking/MCTSAgent.h>
#include <gmds/rlBlocking/MCTSStatePolycube.h>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSAlgorithm
 *  @brief  Class that provides ....
 */
class LIB_GMDS_RLBLOCKING_API MCTSAlgorithm
{
	MCTSTree *tree;
	int max_iter, max_seconds,max_same_quality;
 public:

	/*------------------------------------------------------------------------*/
	/** @brief Constructor.
         * @param
	 */
	MCTSAlgorithm(gmds::cad::GeomManager *AGeom,gmds::blocking::CurvedBlocking *ABlocking,int max_iter = 100000, int max_seconds = 30,int max_same_quality = 10);

	/*------------------------------------------------------------------------*/
	/** @brief  Destructor.	*/
	~MCTSAlgorithm();

	/*------------------------------------------------------------------------*/
	/** @brief  Performs the MCTS algorithm
	 */
	void execute(std::string ANameGeo);

 private:
	/** a geom */
	gmds::cad::GeomManager *m_geom;
	/** a blocking */
	gmds::blocking::CurvedBlocking *m_blocking;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSALGORITHM_H
/*----------------------------------------------------------------------------------------*/
