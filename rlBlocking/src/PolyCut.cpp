/*---------------------------------------------------------------------------*/
#include <iostream>
#include "gmds/rlBlocking/PolyCut.h"
/*---------------------------------------------------------------------------*/
PolyCutAction::PolyCutAction(gmds::TCellID AIdEdge, gmds::TCellID AIdBlock, double AParamCut, PolyCutAction::ActionType ATypeAction)
  :m_AIdEdge(AIdEdge),m_AIdBlock(AIdBlock), m_AParamCut(AParamCut), m_typeAction(ATypeAction)
{}
/*---------------------------------------------------------------------------*/
bool PolyCutAction::operator==(const IAction &AOther) const
{
	const auto &o = (const PolyCutAction &) AOther;
	return m_AIdEdge == o.m_AIdEdge && m_AIdBlock == o.m_AIdBlock
	       && m_AParamCut == o.m_AParamCut && m_typeAction== o.m_typeAction;
}
/*---------------------------------------------------------------------------*/
PolyCutState::PolyCutState(gmds::cad::GeomManager *AGeom, gmds::blocking::CurvedBlocking *ABlocking, std::vector<double> AHist)
	:m_geom(AGeom),m_history(AHist)
{
	m_blocking = new gmds::blocking::CurvedBlocking(*ABlocking);
	gmds::blocking::CurvedBlockingClassifier classifier(m_blocking);
	m_class_blocking = new gmds::blocking::CurvedBlockingClassifier(classifier);
	m_class_blocking->detect_classification_errors(m_class_errors);
}
/*---------------------------------------------------------------------------*/
PolyCutState::PolyCutState(const PolyCutState &AState)
{
	m_blocking = AState.m_blocking;
	m_geom = AState.m_geom;
	m_class_blocking = AState.m_class_blocking;
	m_class_errors = AState.m_class_errors;
	m_history = AState.m_history;
}
/*---------------------------------------------------------------------------*/
PolyCutState::~PolyCutState() noexcept
{
	delete m_class_blocking;
	delete m_blocking;
}
/*---------------------------------------------------------------------------*/
std::vector<std::shared_ptr<IAction>> PolyCutState::get_actions() const
{
	std::vector<std::shared_ptr<IAction>> actions;
	//Add cuts in the actions list
	auto listCuts = m_class_blocking->list_Possible_Cuts();
	for(auto c : listCuts){
		actions.push_back(std::make_shared<PolyCutAction>(c.first,gmds::NullID,c.second,PolyCutAction::Cut));
	}
	//Add deletes block in the list of actions
	auto blocks = m_blocking->get_all_id_blocks();
	for(auto b : blocks){
		actions.push_back(std::make_shared<PolyCutAction>(gmds::NullID,b,0,PolyCutAction::Delete));
	}
	actions.push_back(std::make_shared<PolyCutAction>(gmds::NullID,gmds::NullID,0,PolyCutAction::Classification));

	return actions;
}
/*----------------------------------------------------------------------------*/
double
PolyCutState::get_quality() const
{
	return m_class_errors.non_captured_points.size() * 0.8 + m_class_errors.non_captured_curves.size() * 0.6 +
	       m_class_errors.non_captured_surfaces.size() * 0.4;
}
/*----------------------------------------------------------------------------*/
void PolyCutState::update_class()
{
	gmds::blocking::CurvedBlockingClassifier classifier(m_blocking);
	m_class_blocking = new gmds::blocking::CurvedBlockingClassifier(classifier);
	m_class_blocking->classify(0.2);
	m_class_errors.non_captured_points.clear();
	m_class_errors.non_captured_curves.clear();
	m_class_errors.non_captured_surfaces.clear();
	m_class_blocking->detect_classification_errors(m_class_errors);
}
/*---------------------------------------------------------------------------*/
std::shared_ptr<IState> PolyCutState::apply(const std::shared_ptr<IAction> AAction) const {
	auto next_state = std::make_shared<PolyCutState>(*this);
	auto a = std::dynamic_pointer_cast<PolyCutAction>(AAction);
	if(a->m_typeAction == PolyCutAction::Cut){

	}
	else if(a->m_typeAction==PolyCutAction::Delete){
		//Delete block action
		//TODO ERROR, sometimes, block select not in the current blocks list...Check why !!!
		std::cout<<"LIST BLOCK BLOCKING : "<<std::endl;
		bool b_in_list = false;
		for(auto b : m_blocking->get_all_id_blocks()) {
			std::cout << b << std::endl;
			if (b == a->m_AIdBlock) {
				b_in_list = true;
				break;
			}
		}
		if(b_in_list){
			next_state->m_blocking->remove_block(a->m_AIdBlock);
			next_state->m_class_errors = this->m_class_errors;
		}
		else{
			std::cerr<<"BLOCK A DELETE not in :"<<m_blocking->get_all_id_blocks().front()<<std::endl;
			next_state->m_blocking->remove_block(m_blocking->get_all_id_blocks().front());
		}
	}
	else if(a->m_typeAction==PolyCutAction::Classification){
		//Update Classification
		next_state->update_class();
	}
	else{
		std::cerr << "Warning: Bad type action ! \n Type action :"<<a->m_typeAction<<std::endl;
	}
	return next_state;
}
/*---------------------------------------------------------------------------*/
bool PolyCutState::lost() const
{
	if (!m_history.empty() && m_history.back() < this->get_quality()){
		return true;
	}
	else if (get_actions().empty()){
		return true;
	}
	return false;
}
/*---------------------------------------------------------------------------*/
bool PolyCutState::win() const
{
	if(m_class_errors.non_captured_points.empty() && m_class_errors.non_captured_curves.empty() && m_class_errors.non_captured_surfaces.empty()){
		return true;
	}
	return false;
}
/*----------------------------------------------------------------------------*/
int PolyCutState::check_nb_same_quality() const
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
/*---------------------------------------------------------------------------*/
bool PolyCutState::draw(int nbSameQuality) const
{
	if(check_nb_same_quality() >= nbSameQuality){
		return true;
	}
	return false;
}
/*---------------------------------------------------------------------------*/
bool PolyCutState::is_terminal() const
{
	//the state is terminal if
	// 1. The quality is egal to 0, so, all the points, the curves and the surfaces are captured
	// 2.
	return lost() || win() || draw(3);
}
/*---------------------------------------------------------------------------*/
double PolyCutRewardFunction::evaluate(std::shared_ptr<IState> AState) const {
	if(std::dynamic_pointer_cast<PolyCutState>(AState)->win())
		return 1;
	else if(std::dynamic_pointer_cast<PolyCutState>(AState)->lost())
		return -1;
	return 0;
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const PolyCutAction& PCA){
	if(PCA.m_typeAction== PolyCutAction::Cut){
		os<<"$$$$ Edge id : "<<PCA.m_AIdEdge<<" - Cut Action with param : "<<PCA.m_AParamCut;
	}
	else if(PCA.m_typeAction== PolyCutAction::Delete){
		os<<"$$$$ Block id : "<<PCA.m_AIdBlock<<" - Delete Action";
	}
	else if(PCA.m_typeAction== PolyCutAction::Classification){
		os<<"$$$$ Classification Action";
	}

	return os;
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const PolyCutState& PCS){

	//TODO
}
/*---------------------------------------------------------------------------*/
