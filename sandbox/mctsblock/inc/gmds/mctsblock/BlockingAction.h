/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_BLOCKING_ACTION_H
#define GMDS_MCTS_BLOCKING_ACTION_H
/*----------------------------------------------------------------------------*/
#include <GMDSMctsBlock_export.h>
#include "mcts/IAction.h"
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace mctsblock {
/*----------------------------------------------------------------------------*/
/**@brief this class encapsulate the action of cutting a block edge
 */
class GMDSMctsBlock_API EdgeCutAction : public IAction {
 public:
	/**@brief Computes the state reach from @p AState by applying the current
     * action
     * @param[in] AState the state we start from
     * @return the state that is built from this @p AState when applying the
     * action
	 */
	std::shared_ptr<IState> apply_on(std::shared_ptr<IState> AState) const override;

	EdgeCutAction(const TCellID AEdgeId, const double AParam);

	bool operator==(const IAction& other) const override;
	std::string get_description() const override;

 private:
	/** the id of the edge to cut */
	TCellID m_edge_id;
	/** the parameter where the edge is cutted */
	double m_cut_param;
};

/*----------------------------------------------------------------------------*/
/**@brief this class encapsulate the action of cutting a block edge
 */
class GMDSMctsBlock_API BlockRemovalAction : public IAction {
 public:
	/**@brief Computes the state reach from @p AState by applying the current
     * action
     * @param[in] AState the state we start from
     * @return the state that is built from this @p AState when applying the
     * action
	 */
	std::shared_ptr<IState> apply_on(std::shared_ptr<IState> AState) const override;

	BlockRemovalAction(const TCellID ABlockID);

	bool operator==(const IAction& other) const override;
	std::string get_description() const override;

 private:
	/** the id of the block to remove*/
	TCellID m_block_id;
};

/*----------------------------------------------------------------------------*/
/**@brief this class encapsulate the action of trying to classify boundary edges
 * and faces
 */
class GMDSMctsBlock_API CaptureAction : public IAction {
 public:
	/**@brief Computes the state reach from @p AState by applying the current
     * action
     * @param[in] AState the state we start from
     * @return the state that is built from this @p AState when applying the
     * action
	 */
	std::shared_ptr<IState> apply_on(std::shared_ptr<IState> AState) const override;

	CaptureAction();

	bool operator==(const IAction& other) const override;
	std::string get_description() const override;

};
/*----------------------------------------------------------------------------*/
}     // namespace blocking
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MCTS_BLOCKING_ACTION_H
