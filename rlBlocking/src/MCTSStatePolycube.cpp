/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSStatePolycube.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
MCTSStatePolycube::~MCTSStatePolycube() noexcept
{}
/*----------------------------------------------------------------------------*/
std::queue<MCTSMove *> *
MCTSStatePolycube::actions_to_try() const
{}
/*----------------------------------------------------------------------------*/
MCTSState *
MCTSStatePolycube::next_state(const gmds::MCTSMove *AMove) const
{}
/*----------------------------------------------------------------------------*/
MCTSState::ROLLOUT_STATUS
MCTSStatePolycube::rollout() const
{
	return MCTSState::WIN;
}
/*----------------------------------------------------------------------------*/
bool
MCTSStatePolycube::is_terminal() const
{}
/*----------------------------------------------------------------------------*/
