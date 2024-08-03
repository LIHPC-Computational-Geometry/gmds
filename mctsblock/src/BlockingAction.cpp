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
EdgeCutAction::EdgeCutAction(const gmds::TCellID AEdgeId, const double AParam)
: m_edge_id(AEdgeId), m_cut_param(AParam)
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
	auto s = std::dynamic_pointer_cast<BlockingState>(AState);
	auto b = s->get_blocking();
	b->cut_sheet(m_edge_id, m_cut_param);
	return std::make_shared<BlockingState>(b, s->get_depth()+1, s->get_memory());
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
		auto s = std::dynamic_pointer_cast<BlockingState>(AState);
		auto b = s->get_blocking();
		b->remove_block(m_block_id);
		return std::make_shared<BlockingState>(b, s->get_depth()+1, s->get_memory());
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
	auto s = std::dynamic_pointer_cast<BlockingState>(AState);
	auto b = s->get_blocking();
	std::set<TCellID> nids, eids,fids;
	b->extract_boundary(nids,eids,fids);
	BlockingClassifier(b).try_and_classify_nodes(nids);
//	BlockingClassifier(b).try_and_capture_curves(eids);
	return std::make_shared<BlockingState>(b, s->get_depth()+1, s->get_memory());
}
/*----------------------------------------------------------------------------*/
