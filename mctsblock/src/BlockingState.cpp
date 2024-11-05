/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/BlockingClassifier.h>
#include <gmds/mctsblock/BlockingState.h>
#include <gmds/mctsblock/BlockingAction.h>
#include <gmds/utils/Exception.h>
#include <gmds/cad/GeomManager.h>
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::mctsblock;
/*----------------------------------------------------------------------------*/
const int BlockingState::m_memory_depth = 4;
double BlockingState::weight_nodes = 10000;
double BlockingState::weight_edges = 100;
double BlockingState::weight_faces = 1;
/*----------------------------------------------------------------------------*/
BlockingState::BlockingState(std::shared_ptr<Blocking> AB, int ADepth, std::deque<double> APrevScores)
  : m_blocking(AB), m_depth(ADepth), m_memory_scores(APrevScores)
{
//	m_blocking->reset_classification();
	m_blocking->extract_boundary(m_boundary_node_ids, m_boundary_edge_ids, m_boundary_face_ids);
	BlockingClassifier(m_blocking.get()).try_and_capture(m_boundary_node_ids,m_boundary_edge_ids, m_boundary_face_ids);
	updateMemory(computeScore());

	m_expected_optimal_score = weight_nodes * m_blocking->geom_model()->getNbPoints()
	                           + weight_edges * m_blocking->geom_model()->getNbCurves()
	                           + weight_faces * m_blocking->geom_model()->getNbSurfaces();
}
/*----------------------------------------------------------------------------*/
BlockingState::BlockingState(const gmds::mctsblock::BlockingState &AState) :
  m_blocking(AState.m_blocking),
  m_depth(AState.m_depth),
  m_memory_scores(AState.m_memory_scores),
  m_boundary_node_ids(AState.m_boundary_node_ids),
  m_boundary_edge_ids(AState.m_boundary_edge_ids),
  m_boundary_face_ids(AState.m_boundary_face_ids)
{
}
/*----------------------------------------------------------------------------*/
std::vector<std::shared_ptr<IAction>>
BlockingState::get_actions() const
{
	auto actions = get_possible_cuts();
	auto block_removals = get_possible_block_removals();
	actions.insert(actions.end(),block_removals.begin(),block_removals.end());

	/*TODO Here we could filter equivalent actions!!!
	 * For edge cut, it means for instance, when a cut on an edge will be
	 * equivalent to the cut on a parallel edge by propagation*/
	return actions;
}
/*----------------------------------------------------------------------------*/
bool
BlockingState::win() const
{
	/* we win if we don't have anymore classification errors. It means that the
	   state score, which is the last element of the memory scores, is equal to
	   0.*/
	return (m_memory_scores.back() == m_expected_optimal_score);
} /*----------------------------------------------------------------------------*/
bool
BlockingState::lost() const
{
	bool val = false;
	//we lost if we don't have actions to perform
	if (!win() 	&& this->get_actions().empty())
		return true;

	return false;
	// We lost if the new score have a worst quality than the previous one
	/*if (m_memory_scores.size() > 1 && m_memory_scores[m_memory_scores.size()] < m_memory_scores[m_memory_scores.size()-1]) {
		std::cout << "memo last : " << m_memory_scores[4] << " & memo last -1 :" << m_memory_scores[3] << std::endl;
		return true;
	}


	else
		return false;*/
/*
	// We lost if our score doesn't improve during the last steps.
	// if the memory stack is not full we keep working, so we do not lost
	if (m_memory_scores.size() < BlockingState::m_memory_depth) return false;     // not lost

	// now we check the value of the current score in comparisons to the score history
	return (m_memory_scores[4] >= m_memory_scores[3] ||
	        m_memory_scores[4] >= m_memory_scores[2] ||
	        m_memory_scores[4] >= m_memory_scores[1] ||
	        m_memory_scores[4] >= m_memory_scores[0]);*/
}
/*----------------------------------------------------------------------------*/
bool
BlockingState::draw() const
{return !win() && !lost();}
/*----------------------------------------------------------------------------*/
bool
BlockingState::is_terminal() const
{
	return lost() || win();
}
/*----------------------------------------------------------------------------*/
std::string
BlockingState::write(const std::string &AFileName, const int AStageIndex, const int ANodeId, const int ADepth) const
{
	return "";
}
/*----------------------------------------------------------------------------*/
double
BlockingState::computeScore()
{
   auto errors = BlockingClassifier(m_blocking.get()).detect_classification_errors();

	return weight_nodes*(m_blocking->geom_model()->getNbPoints()-errors.non_captured_points.size())+
	       weight_edges*(m_blocking->geom_model()->getNbCurves()-errors.non_captured_curves.size())+
	       weight_faces*(m_blocking->geom_model()->getNbSurfaces()-errors.non_captured_surfaces.size());
}
/*----------------------------------------------------------------------------*/
void
BlockingState::updateMemory(double AScore)
{
	// We add a new score in the memory stack
	// the oldest score is in 0, the newest is in 3
	if (m_memory_scores.size() == 4) {
		// the stack is full we remove the first element
		for (int i = 0 ; i <3;i++){
			m_memory_scores[i] = m_memory_scores[i+1];
		}
		m_memory_scores[4]= AScore;
	}
	// we add the new score at the back
	m_memory_scores.push_back(AScore);
}
/*----------------------------------------------------------------------------*/
std::vector<std::shared_ptr<IAction>>
BlockingState::get_possible_block_removals() const
{
	auto nodes = m_blocking->get_all_nodes();
	std::set<TCellID> blocks_to_keep;
	for (auto n : nodes) {
		if (n->info().geom_dim == 0) {
			// n is classified onto a geom point
			// if it belongs to a single block, the corresponding block cannot
			// be removed
			auto n_blocks = m_blocking->get_blocks_of_node(n);
			if (n_blocks.size() == 1) {
				blocks_to_keep.insert(n_blocks[0]->info().topo_id);
			}
		}
	}

	// identify blocks with their centroid inside the geometry
	cad::GeomManager *geom = m_blocking->geom_model();

	auto blocks = m_blocking->get_all_blocks();
	for (auto b : blocks) {
		gmds::math::Point pt = m_blocking->get_center_of_block(b);
		bool is_inside = geom->is_in(pt);
		if(is_inside) {
			blocks_to_keep.insert(m_blocking->get_block_id(b));
		}
	}

	// all the other blocks can be removed
	auto all_blocks = m_blocking->get_all_id_blocks();
	std::vector<std::shared_ptr<IAction> > actions;
	for (auto b : all_blocks) {
		//We check if b is a block to keep?
		if (blocks_to_keep.find(b) == blocks_to_keep.end()){
			// b is not to keep, we remove it
			actions.push_back(std::make_shared<BlockRemovalAction>(b));
		}
	}
	return  actions;
}
/*----------------------------------------------------------------------------*/
std::vector<std::shared_ptr<IAction>>
BlockingState::get_possible_cuts() const
{
	double epsilon = 1e-2;
	std::set<std::pair<TCellID, double>> list_cuts;
	// We look which geometric points and curves are not yet captured
	auto non_captured_entities = BlockingClassifier(m_blocking.get()).detect_classification_errors();
	auto non_captured_points = non_captured_entities.non_captured_points;
	auto non_captured_curves = non_captured_entities.non_captured_curves;

	auto all_block_edges = m_blocking->get_all_edges();
	for (auto p : non_captured_points) {
		// the point p is not captured. We look for a cut along a block edge to capture it
		auto action = m_blocking->get_cut_info(p,all_block_edges);
		auto param_cut = std::get<1>(action);
		// We only consider cut of edge that occur inside the edge and not on one of
		// its end points (so for param O or 1)
		if (epsilon < param_cut  && param_cut < 1.0-epsilon) {
			auto e2cut = m_blocking->get_edge(std::get<0>(action)->info().topo_id);
			list_cuts.insert(std::make_pair(e2cut->info().topo_id, std::get<1>(action)));
			// TODO: verify if we always have a cut that is possible for a point???
			// As we cut the edge std::get<0>(action)->info().topo_id, we remove it to avoid multiple actions on the same edge
			std::vector<Blocking::Edge> edges_to_remove;
			m_blocking->get_all_sheet_edges(e2cut,edges_to_remove);

			// Loop through each element to delete
			for (auto e : edges_to_remove) {
				// Find the element in the vector
				auto it = find(all_block_edges.begin(), all_block_edges.end(), e);
				// If the element is found, erase it
				if (it != all_block_edges.end()) {
					all_block_edges.erase(it);
				}
			}
		}
	}

	for (auto c_id : non_captured_curves) {
		auto c = m_blocking->geom_model()->getCurve(c_id);
		gmds::TCoord minXYX[3];
		gmds::TCoord maxXYX[3];
		c->computeBoundingBox(minXYX, maxXYX);

		math::Point bb_corners[8]={
		   math::Point(minXYX[0], minXYX[1], minXYX[2]),
		   math::Point(maxXYX[0], maxXYX[1], maxXYX[2]),
		   math::Point(minXYX[0], minXYX[1], maxXYX[2]),
		   math::Point(maxXYX[0], minXYX[1], maxXYX[2]),
		   math::Point(maxXYX[0], minXYX[1], minXYX[2]),
		   math::Point(minXYX[0], maxXYX[1], maxXYX[2]),
		   math::Point(minXYX[0], maxXYX[1], minXYX[2]),
		   math::Point(maxXYX[0], maxXYX[1], minXYX[2])
		};
		for (auto p:bb_corners) {
			auto cut_info_p = m_blocking->get_cut_info(p, all_block_edges);
			auto param_cut = std::get<1>(cut_info_p);
			// We only consider cut of edge that occur inside the edge and not on one of
			// its end points (so for param O or 1)
			if (epsilon < param_cut  && param_cut < 1.0-epsilon) {
				auto e2cut = m_blocking->get_edge(std::get<0>(cut_info_p)->info().topo_id);
				list_cuts.insert(std::make_pair(e2cut->info().topo_id, std::get<1>(cut_info_p)));
				// TODO: verify if we always have a cut that is possible for a point???
				// As we cut the edge std::get<0>(action)->info().topo_id, we remove it to avoid multiple actions on the same edge
				std::vector<Blocking::Edge> edges_to_remove;
				m_blocking->get_all_sheet_edges(e2cut,edges_to_remove);

				// Loop through each element to delete
				for (auto e : edges_to_remove) {
					// Find the element in the vector
					auto it = find(all_block_edges.begin(), all_block_edges.end(), e);
					// If the element is found, erase it
					if (it != all_block_edges.end()) {
						all_block_edges.erase(it);
					}
				}
			}

		}
	}

	std::vector<std::shared_ptr<IAction> > actions;
	for(auto cut_info:list_cuts)
		actions.push_back(std::make_shared<EdgeCutAction>(cut_info.first,cut_info.second));

	return actions;
}