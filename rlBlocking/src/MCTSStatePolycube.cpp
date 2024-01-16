/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSStatePolycube.h>
#include <queue>
/*----------------------------------------------------------------------------*/
using namespace gmds;


/*----------------------------------------------------------------------------*/
MCTSStatePolycube::MCTSStatePolycube(gmds::cad::GeomManager* AGeom, gmds::blocking::CurvedBlocking* ABlocking,
                                     std::vector<double> hist )
	:m_geom(AGeom),m_history(hist)
{
	m_blocking = new blocking::CurvedBlocking(*ABlocking);
	gmds::blocking::CurvedBlockingClassifier classifier(m_blocking);
	m_class_blocking = new blocking::CurvedBlockingClassifier(classifier);
	m_class_errors = m_class_blocking->classify(0.2);
	;}
/*----------------------------------------------------------------------------*/
MCTSStatePolycube::~MCTSStatePolycube() noexcept
{
	delete m_class_blocking;
	delete m_blocking;
}
/*----------------------------------------------------------------------------*/
std::deque<MCTSMove *> *
MCTSStatePolycube::actions_to_try() const
{
	std::deque<MCTSMove *> *Q = new std::deque<MCTSMove *>();
	if (m_class_errors.non_captured_points.size()== 0){
		std::cout<<"POINTS CAPT :"<<std::endl;
		if(m_class_errors.non_captured_curves.size()==0){
			std::cout<<"CURVES CAPT :"<<std::endl;
			auto blocks = m_blocking->get_all_id_blocks();
			for(auto b : blocks){
				Q->push_back(new MCTSMovePolycube(NullID,b,0,2));
			}
		}
		else{
			std::cout<<"NB CURVES CAPT :"<< m_class_errors.non_captured_curves.size()<<std::endl;
			auto listPossibleCuts = m_class_blocking->list_Possible_Cuts();
			for(auto c : listPossibleCuts){
				Q->push_back(new MCTSMovePolycube(c.first,NullID,c.second,1));
			}
		}
	}
	else{
		std::cout<<"POINTS NO CAPT :"<<std::endl;
		for(auto p : m_class_errors.non_captured_points){
			std::cout<<"p :"<<p<<std::endl;
		}
		auto listPossibleCuts = m_class_blocking->list_Possible_Cuts();
		for(auto c : listPossibleCuts){
			Q->push_back(new MCTSMovePolycube(c.first,NullID,c.second,1));
		}
	}
	std::cout<<"LIST ACTIONS :"<<std::endl;
	for(auto &m : *Q){
		m->print();
	}
	return Q;
}
/*----------------------------------------------------------------------------*/
MCTSState
   *MCTSStatePolycube::next_state(const gmds::MCTSMove *AMove) const
{
	std::cout<<"==================== EXECUTE ACTION ! ===================="<<std::endl;
	MCTSMovePolycube *m = (MCTSMovePolycube *) AMove;
	std::vector<double> hist_update = get_history();
	hist_update.push_back(get_quality());
	gmds::blocking::CurvedBlocking* new_b = new gmds::blocking::CurvedBlocking(*m_blocking);
	MCTSStatePolycube *new_state = new MCTSStatePolycube(this->m_geom,new_b,hist_update);
	if(m->m_typeMove == 2){
		//TODO ERROR, sometimes, block select not in the current blocks list...Check why !!!
		std::cout<<"LIST BLOCK BLOCKING : "<<std::endl;
		bool b_in_list = false;
		for(auto b : m_blocking->get_all_id_blocks()){
			std::cout<<b<<std::endl;
			if(b == m->m_AIdBlock){
				b_in_list = true;
				break;
			}
		}
		if(b_in_list){
			std::cout<<"BLOCK A DELETE :"<<m->m_AIdBlock<<std::endl;
			new_state->m_blocking->remove_block(m->m_AIdBlock);
		}
		else{
			std::cout<<"BLOCK A DELETE :"<<m_blocking->get_all_id_blocks().front()<<std::endl;
			new_state->m_blocking->remove_block(m_blocking->get_all_id_blocks().front());
		}

		new_state->update_class();
		//SAVE Blocking vtk
		std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/cb2/";
		std::string id_act = std::to_string(m->m_AIdEdge);
		std::string name_file = "cb2_action"+ id_act +".vtk";
		new_state->m_blocking->save_vtk_blocking(name_save_folder+name_file);
		return new_state;
	}
	else if(m->m_typeMove ==1) {
		new_state->m_blocking->cut_sheet(m->m_AIdEdge,m->m_AParamCut);
		new_state->update_class();
		//SAVE Blocking vtk
		std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/cb2/";
		std::string id_act = std::to_string(m->m_AIdEdge);
		std::string name_file = "cb2_action"+ id_act +".vtk";
		new_state->m_blocking->save_vtk_blocking(name_save_folder+name_file);
		return new_state;
	}
	else{
		std::cerr << "Warning: Bad type move ! \n Type move :" << m->m_typeMove << " & ID " << m->m_AIdEdge<< std::endl;
		return new_state;
	}
	std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/";
	std::string id_act = std::to_string(m->m_AIdEdge);
	std::string name_file = "M1_action"+ id_act +".vtk";
	m_blocking->save_vtk_blocking(name_save_folder+name_file);

}
/*----------------------------------------------------------------------------*/
double
MCTSStatePolycube::state_rollout() const
{
	std::cout<<"STATE ROLLOUT"<<std::endl;
	std::cout<<"VAL TERMINAL begin "<<is_terminal()<<std::endl;
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

	long long r;
	int a;
	MCTSStatePolycube *curstate = (MCTSStatePolycube *) this;   // TODO: ignore const...
	srand(time(NULL));
	bool first = true;
	do {
		std::deque<MCTSMove *> *list_action = actions_to_try();
		//Get first move/action
		//But, maybe, better to take rand move if its a delete move...
		MCTSMove *firstMove = list_action->front(); //TODO: implement random move when only delete moves is possible
		list_action->pop_front();

		MCTSStatePolycube *old = curstate;
		std::cout<<"===== SIZE UNTRIED ACTIONS : "<<list_action->size()+1<<" ====="<<std::endl;
		curstate = (MCTSStatePolycube *) curstate->next_state(firstMove);
		if (!first) {
			delete old;
		}
		first = false;
	} while (!curstate->is_terminal());

	if(curstate->result_terminal() == WIN){
		res=1;
	}
	else if (curstate->result_terminal() == LOSE) {
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
	if (m_class_errors.non_captured_points.empty() && m_class_errors.non_captured_curves.empty() && m_class_errors.non_captured_surfaces.empty()) {
		return WIN;
	}
	else if (check_nb_same_quality() >= 3){
		return DRAW;
	}
	else if (!m_history.empty() && m_history.back() < this->get_quality()){
		return LOSE;
	}
	else if (this->actions_to_try()->empty()){
		return LOSE;
	}
	std::cerr << "ERROR: NOT terminal state ..." << std::endl;
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
	else if(!m_history.empty() && m_history.back() < this->get_quality()){
		return true;
	}
	else if(this->actions_to_try()->empty()){
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
void MCTSStatePolycube::update_class()
{
	gmds::blocking::CurvedBlockingClassifier classifier(m_blocking);
	m_class_blocking = new blocking::CurvedBlockingClassifier(classifier);
	m_class_errors = m_class_blocking->classify(0.2);
}
/*----------------------------------------------------------------------------*/
