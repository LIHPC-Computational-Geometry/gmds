/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_MCTSSTATE_POLYCUBE_H
#define GMDS_MCTSSTATE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include <gmds/rlBlocking/MCTSState.h>
#include <gmds/rlBlocking/MCTSMovePolycube.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------------------*/
/** @class  MCTSState
 *  @brief  Class that provides the interface to be implemented for performing the
 *  MCST algorithm
 */
class LIB_GMDS_RLBLOCKING_API MCTSStatePolycube: public MCTSState{

 public:
	/*------------------------------------------------------------------------*/
	/** @brief  Constructore
	 */
	MCTSStatePolycube(gmds::cad::GeomManager *Ageom, gmds::blocking::CurvedBlocking *ABlocking,
	                  std::vector<double> AHist,std::string ANameGeom);
	   /*------------------------------------------------------------------------*/
	/** @brief  Destructor
	 */
	~MCTSStatePolycube();
	/*------------------------------------------------------------------------*/
	/** @brief  Gives the set of actions that can be tried from the current state
	 */
	std::deque<MCTSMove *> *actions_to_try() const ;
	/*------------------------------------------------------------------------*/
	/** @brief  Performs the @p AMove to change of states
	 * @param[in] AMove the movement to apply to get to a new state
	 */
	MCTSState *next_state(const MCTSMove *AMove) const;
	/*------------------------------------------------------------------------*/
	/** @brief Rollout from this state (random simulation)
	 *  @return the rollout status
	 */
	double state_rollout() const;

	/** @brief check the history of qualities
	 *  @return nb of same quality from the history
	 */
	int check_nb_same_quality() const;
	/** @brief  check the result of a terminal state
	 * @return Win = all elements are capt, Lose: parent_quality &lt enfant_quality,
	 * Draw : same quality for a long time
	 */
	ROLLOUT_STATUS result_terminal() const;
	/*------------------------------------------------------------------------*/
	/** @brief  Indicate if we have a terminal state (win=true, fail=false)
	 * @return true if we have a leaf (in the sense of a traditional tree)
	 */
	bool is_terminal() const;
	/** @brief return the blocking quality
	 * */
	double get_quality() const;
	/** @brief return the geom */
	gmds::cad::GeomManager *get_geom();
	/** @brief return the current blocking */
	gmds::blocking::CurvedBlocking *get_blocking();
	/** @brief return the current classifier */
	gmds::blocking::CurvedBlockingClassifier *get_class();
	/** @brief return the current classification */
	gmds::blocking::ClassificationErrors get_errors();
	/** @brief return the history of the parents quality */
	std::vector<double> get_history() const;

	/** @brief update the classification of a state */
	void update_class();

 private :
	/** @brief the curved blocking of the current state */
	gmds::blocking::CurvedBlocking* m_blocking;
	gmds::cad::GeomManager* m_geom;
	gmds::blocking::CurvedBlockingClassifier* m_class_blocking;
	gmds::blocking::ClassificationErrors m_class_errors;
	std::vector<double> m_history;
	std::string m_name_geom;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------------------*/
#endif     // GMDS_MCTSSTATE_POLYCUBE_H
/*----------------------------------------------------------------------------------------*/
