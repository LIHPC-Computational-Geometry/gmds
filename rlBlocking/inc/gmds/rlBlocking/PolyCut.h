#ifndef GMDS_POLYCUT_H
#define GMDS_POLYCUT_H
/*---------------------------------------------------------------------------*/
#include "LIB_GMDS_RLBLOCKING_export.h"
#include "gmds/rlBlocking/IStateSPAM.h"
#include "gmds/rlBlocking/IActionSPAM.h"
#include "gmds/rlBlocking/IRewardFunctionSPAM.h"
#include <gmds/utils/CommonTypes.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <iostream>
#include <memory>
/*---------------------------------------------------------------------------*/

class LIB_GMDS_RLBLOCKING_API PolyCutAction : public IAction{
 public:
	enum ActionType{
		None, Cut, Delete, Classification
	};
	gmds::TCellID m_AIdEdge;
	gmds::TCellID m_AIdBlock;
	double m_AParamCut;
	ActionType m_typeAction;

	/** @brief  Overloaded ==
	 */
	PolyCutAction(gmds::TCellID AIdEdge = -1,gmds::TCellID AIdBlock = -1 , double AParamCut = 0,ActionType ATypeAction = None );
	bool operator==(const IAction& AOther) const override;
	friend std::ostream& operator<<(std::ostream& os, const PolyCutAction& PCA);

	void print();
};
/*---------------------------------------------------------------------------*/
class LIB_GMDS_RLBLOCKING_API PolyCutState : public IState{
 public:
	/** @brief the curved blocking of the current state */
	gmds::blocking::CurvedBlocking* m_blocking;
	gmds::cad::GeomManager* m_geom;
	/** @brief the classification of the blocking */
	gmds::blocking::CurvedBlockingClassifier* m_class_blocking;
	/** @brief points, curves, surfaces not yet captured */
	gmds::blocking::ClassificationErrors m_class_errors;
	std::vector<double> m_history;
	std::vector<std::shared_ptr<IAction>> m_actions_list;


	PolyCutState(gmds::cad::GeomManager *AGeom, gmds::blocking::CurvedBlocking *ABlocking, std::vector<double> AHist);
	PolyCutState(const PolyCutState& AState); //copy
	/** @brief Destructor
	 */
	 ~PolyCutState();


	std::vector<std::shared_ptr<IAction>> creat_actions_list();
	std::vector<std::shared_ptr<IAction>> get_actions() const override;
	/**@brief Computes the state reach from the current one when applying @p AAction
     * @param AAction the action to apply
     * @return the state that is built from this one when applying @p AAction
	 */
	std::shared_ptr<IState> apply(std::shared_ptr<IAction>  AAction) const override;

	bool is_terminal() const override;
	bool  lost() const;
	bool  win() const;
	bool  draw(int nbSameQuality) const;

	int get_nb_blocks();

	void write(const std::string& AFileName,
	                   const int AStageIndex,
	                   const int ANodeId,
	                   const int ADepth) const;


	std::vector<gmds::TCellID> get_blocks_id() const;
	std::vector<gmds::TCellID> get_faces_id() const;
	std::vector<gmds::TCellID> get_edges_id() const;
	std::vector<gmds::TCellID> get_nodes_id() const;

	void get_info_blocking() const;
	friend std::ostream& operator<<(std::ostream& os, const PolyCutState& PCS);


	/** @brief check the history of qualities
	 *  @return nb of same quality from the history
	 */
	int check_nb_same_quality() const;
	/** @brief return the blocking quality
	 * */
	double get_quality() const;
	/** @brief update the classification of a state */
	void update_class();
};
/*---------------------------------------------------------------------------*/
struct PolyCutRewardFunction: public IRewardFunction{
	double evaluate(std::shared_ptr<IState> AState) const override;
};

/*---------------------------------------------------------------------------*/
#endif     // GMDS_POLYCUT_H
/*---------------------------------------------------------------------------*/
