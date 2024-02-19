/*----------------------------------------------------------------------------*/
#include <gmds/rlBlocking/MCTSStatePolycube.h>
#include <queue>
/*----------------------------------------------------------------------------*/
using namespace gmds;


/*----------------------------------------------------------------------------*/
MCTSStatePolycube::MCTSStatePolycube(gmds::cad::GeomManager* AGeom, gmds::blocking::CurvedBlocking* ABlocking,
                                     std::vector<double> hist,std::string ANameGeom )
	:m_geom(AGeom),m_history(hist),m_name_geom(ANameGeom)
{
	m_blocking = new blocking::CurvedBlocking(*ABlocking);
	gmds::blocking::CurvedBlockingClassifier classifier(m_blocking);
	m_class_blocking = new blocking::CurvedBlockingClassifier(classifier);
	m_class_blocking->detect_classification_errors(m_class_errors);
	//m_class_errors = m_class_blocking->classify(0.2); /TODO add detected_errors()
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
	std::deque<MCTSMove *> *listActions = new std::deque<MCTSMove *>();
		//Add cuts in the list of actions
	   auto listPossibleCuts = m_class_blocking->list_Possible_Cuts();
	   for(auto c : listPossibleCuts){
		   listActions->push_back(new MCTSMovePolycube(c.first,NullID,c.second,1));
	   }
	   //Add deletes block in the list of actions
	   auto blocks = m_blocking->get_all_id_blocks();
	   for(auto b : blocks){
		   listActions->push_back(new MCTSMovePolycube(NullID,b,0,2));
	   }

	   listActions->push_back(new MCTSMovePolycube(NullID,NullID,0,3));

	std::cout<<"LIST ACTIONS :"<<std::endl;
	for(auto &m : *listActions){
		m->print();
	}
	return listActions;
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
	MCTSStatePolycube *new_state = new MCTSStatePolycube(this->m_geom,new_b,hist_update,this->m_name_geom);
	if(m->m_typeMove == 2){
		//Delete block action
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
			new_state->m_class_errors = this->m_class_errors;
		}
		else{
			std::cout<<"BLOCK A DELETE not in :"<<m_blocking->get_all_id_blocks().front()<<std::endl;
			new_state->m_blocking->remove_block(m_blocking->get_all_id_blocks().front());
		}

		//new_state->update_class();
		//SAVE Blocking vtk
		std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/"+m_name_geom+"/";
		std::string id_act = std::to_string(m->m_AIdBlock);
		std::string name_file = m_name_geom+"_action"+ id_act;
		//new_state->m_blocking->save_vtk_blocking(name_save_folder+name_file);
		return new_state;
	}
	else if(m->m_typeMove ==1) {
		//Cut sheet action
		new_state->m_blocking->cut_sheet(m->m_AIdEdge,m->m_AParamCut);
		new_state->update_class();
		//SAVE Blocking vtk
		std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/"+m_name_geom+"/";
		std::string id_act = std::to_string(m->m_AIdEdge);
		std::string name_file = m_name_geom+"_action"+ id_act;
		//new_state->m_blocking->save_vtk_blocking(name_save_folder+name_file);
		return new_state;
	}
	else if(m->m_typeMove == 3){
		//Update Classification
		new_state->update_class();
		return new_state;
	}
	else{
		std::cerr << "Warning: Bad type move ! \n Type move :" << m->m_typeMove << " & ID " << m->m_AIdEdge<< std::endl;
		return new_state;
	}
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
			//SAVE Blocking vtk
			std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/"+m_name_geom+"/";
			//std::string id_act = std::to_string(m->m_AIdEdge);
			std::string name_file = m_name_geom+"_action_win";
			this->m_blocking->save_vtk_blocking(name_save_folder+name_file);
			res=1;
			return res;
		}
		else if (MCTSStatePolycube::result_terminal() == LOSE) {
			//SAVE Blocking vtk
			std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_Correction_Class_Dev/saveResults/"+m_name_geom+"/";
			//std::string id_act = std::to_string(m->m_AIdEdge);
			std::string name_file = m_name_geom+"_action_lose.vtk";
			this->m_blocking->save_vtk_blocking(name_save_folder+name_file);
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
		std::deque<MCTSMove *> *list_action = curstate->actions_to_try();
		if(list_action->size()==1){
			std::cout<<list_action->front()<<std::endl;
		}
		//Get random move/action
		MCTSMove *randMove = randomMove(*list_action) ; //TODO: implement random move when only delete moves is possible
		std::cout<<"ACTION SELECT : "<<std::endl;
		randMove->print();

		MCTSStatePolycube *old = curstate;
		std::cout<<"===== SIZE UNTRIED ACTIONS : "<<list_action->size()+1<<" ====="<<std::endl;
		curstate = (MCTSStatePolycube *) curstate->next_state(randMove);
		if (!first) {
			delete old;
		}
		first = false;
	} while (!curstate->is_terminal());

	if(curstate->result_terminal() == WIN){
		//SAVE Blocking vtk
		std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_fix_cut_sheet/saveResults/"+m_name_geom+"/";
		//std::string id_act = std::to_string(m->m_AIdEdge);
		std::string name_file = m_name_geom+"_action_win.vtk";
		curstate->m_blocking->save_vtk_blocking(name_save_folder+name_file);
		res=1;
	}
	else if (curstate->result_terminal() == LOSE) {
		//SAVE Blocking vtk
		std::string name_save_folder = "/home/bourmaudp/Documents/PROJETS/gmds/gmds_fix_cut_sheet/saveResults/"+m_name_geom+"/";
		//std::string id_act = std::to_string(m->m_AIdEdge);
		std::string name_file = m_name_geom+"_action_lose.vtk";
		curstate->m_blocking->save_vtk_blocking(name_save_folder+name_file);
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
void MCTSStatePolycube::update_class()
{
	gmds::blocking::CurvedBlockingClassifier classifier(m_blocking);
	m_class_blocking = new blocking::CurvedBlockingClassifier(classifier);
	m_class_blocking->classify(0.2);
	m_class_errors.non_captured_points.clear();
	m_class_errors.non_captured_curves.clear();
	m_class_errors.non_captured_surfaces.clear();
	m_class_blocking->detect_classification_errors(m_class_errors);
}
/*----------------------------------------------------------------------------*/
bool
MCTSStatePolycube::is_terminal() const
{
	if (this->get_quality() == 0) {
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
MCTSMove* MCTSStatePolycube::randomMove(std::deque<MCTSMove *> AListActions) const
{
	if (AListActions.empty()) {
		return nullptr; //The deque is empty
	}
	// Générateur de nombres aléatoires
	std::random_device rd;
	std::mt19937 gen(rd());

	// Sélection d'un index aléatoire dans la deque
	std::uniform_int_distribution<size_t> distribution(0, AListActions.size() - 1);
	size_t randomIndex = distribution(gen);

	// Retourne l'élément correspondant à l'index sélectionné
	return AListActions[randomIndex];

}
/*----------------------------------------------------------------------------*/
