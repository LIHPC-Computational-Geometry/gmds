/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSAlgorithm.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
MCTSAlgorithm::MCTSAlgorithm(gmds::cad::GeomManager* AGeom,gmds::blocking::CurvedBlocking* ABlocking,int max_iter, int max_seconds, int max_same_quality)
	: m_geom(AGeom),m_blocking(ABlocking),max_iter(max_iter), max_seconds(max_seconds),max_same_quality(max_same_quality)
{	std::vector<double>  hist_empty;
	MCTSStatePolycube *init_state = new MCTSStatePolycube(m_geom,m_blocking,hist_empty);
	tree = new MCTSTree(init_state);}
/*----------------------------------------------------------------------------*/
MCTSAlgorithm::~MCTSAlgorithm(){;}
/*----------------------------------------------------------------------------*/
void MCTSAlgorithm::execute()
{
	std::cout<<"==========================================================="<<std::endl;
	std::cout<<"==================== BEGIN EXECUTE ALGO===================="<<std::endl;
	std::cout<<"==========================================================="<<std::endl;

	// This is just to test the example, not ideally what one should be doing... TODO
	bool done;

	MCTSState *state = new MCTSStatePolycube(this->m_geom, this->m_blocking, std::vector<double> ());
	//state->print();                           // IMPORTANT: state will be garbage after advance_tree()
	MCTSAgent agent(state, 100);
	do {
		agent.feedback();
		agent.genmove();
		// TODO: This way we don't check if the enemy move ends the game but it's our responsibility to check that, not the tree's...
		const MCTSState *new_state = agent.get_current_state();
		new_state->print();
//		if (new_state->is_terminal()) {
//			winner = ((const TicTacToe_state *) new_state)->get_winner();
//			break;
//		}
		done = new_state->is_terminal();
	} while (!done);


	std::cout<<"==========================================================="<<std::endl;
	std::cout<<"===================== END EXECUTE ALGO====================="<<std::endl;
	std::cout<<"==========================================================="<<std::endl;
}
/*----------------------------------------------------------------------------*/
