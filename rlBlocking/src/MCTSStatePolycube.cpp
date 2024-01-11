/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSStatePolycube.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;

/*----------------------------------------------------------------------------*/
MCTSStatePolycube::MCTSStatePolycube(gmds::cad::GeomManager* AGeom, gmds::blocking::CurvedBlocking* ABlocking,
                                     std::vector<double> hist )
	:m_geom(AGeom),m_blocking(ABlocking),m_history(hist)
{
	gmds::blocking::CurvedBlockingClassifier classifier(m_blocking);
	m_class_blocking = new blocking::CurvedBlockingClassifier(classifier);
	m_class_errors = m_class_blocking->classify();
	;}
/*----------------------------------------------------------------------------*/
MCTSStatePolycube::~MCTSStatePolycube() noexcept
{ delete m_class_blocking;}
/*----------------------------------------------------------------------------*/
std::queue<MCTSMove *> *
MCTSStatePolycube::actions_to_try() const
{
	std::queue<MCTSMove *> *Q = new std::queue<MCTSMove *>();
	if (m_class_errors.non_captured_points.size()== 0){
		if(m_class_errors.non_captured_curves.size()==0){
			auto blocks = m_blocking->get_all_id_blocks();
			for(auto b : blocks){
				Q->push(new MCTSMovePolycube(NullID,b,0,0));
			}
		}
		else{
			auto listPossibleCuts = m_class_blocking->list_Possible_Cuts();
			for(auto c : listPossibleCuts){
				Q->push(new MCTSMovePolycube(c.first,NullID,c.second,1));
			}
		}
	}
	else{
		auto listPossibleCuts = m_class_blocking->list_Possible_Cuts();
		for(auto c : listPossibleCuts){
			Q->push(new MCTSMovePolycube(c.first,NullID,c.second,1));
		}
	}
	return Q;
}
/*----------------------------------------------------------------------------*/
MCTSState
   *MCTSStatePolycube::next_state(const gmds::MCTSMove *AMove) const
{
	MCTSMovePolycube *m = (MCTSMovePolycube *) AMove;
	std::vector<double> hist_update = get_history();
	hist_update.push_back(get_quality());
	MCTSStatePolycube *new_state = new MCTSStatePolycube(this->m_geom,this->m_blocking,hist_update);
	if(m->m_typeMove == 0){
		new_state->m_blocking->remove_block(m->m_AIdBlock);
		return new_state;
	}
	else if(m->m_typeMove ==1) {
		new_state->m_blocking->cut_sheet(m->m_AIdEdge,m->m_AParamCut);
		return new_state;
	}
	else{
		std::cerr << "Warning: Bad type move !" << std::endl;
		return new_state;
	}

}
/*----------------------------------------------------------------------------*/
double
MCTSStatePolycube::state_rollout() const
{
	std::cout<<"STATE ROLLOUT"<<std::endl;
	double res;
	if(is_terminal()){
		if(MCTSStatePolycube::result_terminal() == WIN){
			res=1;
			return res;
		}
		else if (MCTSStatePolycube::result_terminal() == LOSE) {
			res=-1;
			return res;
		}
		else{
			//Draw
			res=0;
			return res;
		}
	}
	std::queue<MCTSMove *> *list_action = actions_to_try();
	long long r;
	int a;
	MCTSStatePolycube *curstate = (MCTSStatePolycube *) this;   // TODO: ignore const...
	srand(time(NULL));
	bool first = true;
	do {
		if (list_action->empty()) {
			std::cerr << "Warning: Ran out of available moves and state is not terminal?";
			return 0.0;
		}
		//Get first move/action
		//But, maybe, better to take rand move if its a delete move...
		MCTSMove *firstMove = list_action->front(); //TODO: implement random move when only delete moves is possible
		list_action->pop();

		MCTSStatePolycube *old = curstate;
		curstate = (MCTSStatePolycube *) curstate->next_state(firstMove);
		if (!first) {
			delete old;
		}
		first = false;
	} while (!curstate->is_terminal());

	if(MCTSStatePolycube::result_terminal() == WIN){
		res=1;
	}
	else if (MCTSStatePolycube::result_terminal() == LOSE) {
		res=-1;
	}
	else{
		//Draw
		res=0;
	}
	delete curstate;
	return res;
}
/*----------------------------------------------------------------------------*/
MCTSStatePolycube::ROLLOUT_STATUS
MCTSStatePolycube::result_terminal() const
{
	int max_nb_same = 3;
	if (get_quality() == 0) {
		return WIN;
	}
	else if (check_nb_same_quality() >= max_nb_same){
		return DRAW;
	}
	else if (m_history.back() < get_quality()){
		return LOSE;
	}
	std::cerr << "ERROR: NOT terminal state !" << std::endl;
	return DRAW;
}
/*----------------------------------------------------------------------------*/
int MCTSStatePolycube::check_nb_same_quality() const
{
	int nb_same = 0;
	double state_quality = get_quality();
	for (int i = m_history.size() - 1; i >= 0; --i) {
		if(m_history[i] == state_quality){
			nb_same++;
		}
		else{
			break;
		}
	}
	return nb_same;
}
/*----------------------------------------------------------------------------*/
bool
MCTSStatePolycube::is_terminal() const
{
	if (m_class_errors.non_captured_points.empty() && m_class_errors.non_captured_curves.empty() && m_class_errors.non_captured_surfaces.empty()) {
		return true;
	}
	else if(check_nb_same_quality() >= 3){
		return true;
	}
	else if(!m_history.empty() && m_history.back() < get_quality()){
		return true;
	}
	else {
		std::cout<<"NOT TERMINAL STATE"<<std::endl;
		return false;
	}
}
/*----------------------------------------------------------------------------*/
double
   MCTSStatePolycube::get_quality() const
{
	return m_class_errors.non_captured_points.size() * 0.8 + m_class_errors.non_captured_curves.size() * 0.6 +
	       m_class_errors.non_captured_surfaces.size() * 0.4;
}
/*----------------------------------------------------------------------------*/
gmds::cad::GeomManager* MCTSStatePolycube::get_geom(){
	return m_blocking->geom_model();
}

/*----------------------------------------------------------------------------*/
gmds::blocking::CurvedBlocking *MCTSStatePolycube::get_blocking()
{
	return m_blocking;
}

/*----------------------------------------------------------------------------*/
gmds::blocking::CurvedBlockingClassifier *MCTSStatePolycube::get_class()
{
	return m_class_blocking;
}
/*----------------------------------------------------------------------------*/
gmds::blocking::ClassificationErrors MCTSStatePolycube::get_errors()
{
	return m_class_errors;
}
/*----------------------------------------------------------------------------*/
std::vector<double> MCTSStatePolycube::get_history() const
{
	return m_history;
}
/*----------------------------------------------------------------------------*/
