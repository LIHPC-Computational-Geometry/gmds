/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/BlockingAction.h>
#include <gmds/mctsblock/BlockingState.h>
#include <gmds/mctsblock/BlockingClassifier.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::mctsblock;
/*----------------------------------------------------------------------------*/
/** EDGE CUT ACTION */
/*----------------------------------------------------------------------------*/
EdgeCutAction::EdgeCutAction(const gmds::TCellID AEdgeId, const double AParam, const gmds::math::Point APoint)
   :m_edge_id(AEdgeId), m_cut_param(AParam),m_capt_point(APoint)
{}
/*----------------------------------------------------------------------------*/
bool
EdgeCutAction::operator==(const IAction &AOther) const
{
	const auto &o = (const EdgeCutAction &) AOther;
	return m_edge_id == o.m_edge_id && m_cut_param == o.m_cut_param;
}
/*----------------------------------------------------------------------------*/
std::shared_ptr<IState>
EdgeCutAction::apply_on(std::shared_ptr<IState> AState) const
{
	auto current_state = std::dynamic_pointer_cast<BlockingState>(AState);
	auto new_blocking = std::make_shared<Blocking>(* current_state->get_blocking());
	new_blocking->cut_sheet(m_edge_id, m_cut_param);
	auto next_state =  std::make_shared<BlockingState>(new_blocking, current_state->get_depth()+1, current_state->get_memory());
	return next_state;
}
/*----------------------------------------------------------------------------*/
std::string
EdgeCutAction::get_description() const
{
	return "Cut edge "+std::to_string(m_edge_id)+" with param "+std::to_string(m_cut_param)+", try to capt point XYZ("+std::to_string(m_capt_point.X())+"), "+
	       std::to_string(m_capt_point.Y())+", "+std::to_string(m_capt_point.Z());
}
/*----------------------------------------------------------------------------*/
/** BLOCK REMOVAL ACTION */
/*----------------------------------------------------------------------------*/
BlockRemovalAction::BlockRemovalAction(const gmds::TCellID ABlockID)
: m_block_id(ABlockID)
{
}
/*----------------------------------------------------------------------------*/
bool
BlockRemovalAction::operator==(const IAction &AOther) const
{
	const auto &o = (const BlockRemovalAction &) AOther;
	return m_block_id == o.m_block_id;
}
/*----------------------------------------------------------------------------*/
std::shared_ptr<IState>
BlockRemovalAction::apply_on(std::shared_ptr<IState> AState) const
{
	   auto current_state = std::dynamic_pointer_cast<BlockingState>(AState);
	   auto new_blocking = std::make_shared<Blocking>(* current_state->get_blocking());
	   new_blocking->remove_block(m_block_id);
	   auto next_state =  std::make_shared<BlockingState>(new_blocking, current_state->get_depth()+1, current_state->get_memory());
	   return next_state;
}
/*----------------------------------------------------------------------------*/
std::string
BlockRemovalAction::get_description() const
{
	   return "Remove block "+std::to_string(m_block_id);
}
/*----------------------------------------------------------------------------*/
/** CLASSIFICATION ACTION */
/*----------------------------------------------------------------------------*/
CaptureAction::CaptureAction()
{}
/*----------------------------------------------------------------------------*/
bool
CaptureAction::operator==(const IAction &AOther) const
{
	return true;
}
/*----------------------------------------------------------------------------*/
std::shared_ptr<IState>
CaptureAction::apply_on(std::shared_ptr<IState> AState) const
{


	auto current_state = std::dynamic_pointer_cast<BlockingState>(AState);
	auto new_blocking = std::make_shared<Blocking>(* current_state->get_blocking());
	std::set<TCellID> nids, eids,fids;
	new_blocking->extract_boundary(nids,eids,fids);
	//	BlockingClassifier(new_blocking.get()).try_and_classify_nodes(nids);
	//	BlockingClassifier(b).try_and_capture_curves(eids);

	auto next_state =  std::make_shared<BlockingState>(new_blocking, current_state->get_depth()+1, current_state->get_memory());
	return next_state;
}



/*----------------------------------------------------------------------------*/
std::string
CaptureAction::get_description() const
{
	return "Capture";
}
/*----------------------------------------------------------------------------*/
