/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_BLOCKING_CLASSIFIER_H
#define GMDS_MCTS_BLOCKING_CLASSIFIER_H
/*----------------------------------------------------------------------------*/
#include <GMDSMctsBlock_export.h>
#include <gmds/mctsblock/Blocking.h>
/*----------------------------------------------------------------------------*/
#include <string>
#include <tuple>
#include <type_traits>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace mctsblock {
/*----------------------------------------------------------------------------*/
/**@struct ClassificationErrors
 * @brief This structure lists the classification errors when a curved blocking
 * 		 structure is classified onto a geometrical model.
 */
struct ClassificationErrors
{
	/*** ids of geometric points that are not captured by the blocking */
	std::vector<int> non_captured_points;
	/*** ids of geometric curves that are not captured by the blocking */
	std::vector<int> non_captured_curves;
	/*** ids of geometric surfaces that are not captured by the blocking */
	std::vector<int> non_captured_surfaces;
	/*** ids of block nodes thar are not classified  */
	std::vector<int> non_classified_nodes;
	/*** ids of block edges thar are not classified  */
	std::vector<int> non_classified_edges;
	/*** ids of block faces thar are not classified  */
	std::vector<int> non_classified_faces;
};
/*----------------------------------------------------------------------------*/
/**@class BlockingClassifier
 * @brief Provide functions to check and update the classification of a curved
 * 		 blocking structure onto a geometrical model.
 *
 */
class GMDSMctsBlock_API BlockingClassifier
{
 public:
	/** @brief Constructor
	 * @param[in] ABlocking a blocking structure
	 */
	BlockingClassifier(Blocking *ABlocking);

	/** @brief  Destructor
	 */
	virtual ~BlockingClassifier();

	/**@brief This operation detects error in the classification.
	 * Errors are geometrical cells that are not catch by the
	 * blocking cells.
	 */
	ClassificationErrors detect_classification_errors();

	/**@brief this method clears the classification of the blocking
	 * structure. All the cells are then unclassified (geom_dim=4
	 * and geom_id=NullID)
	 *
	 */
	void clear_classification();

	/**@brief This method tries and classifies the unclassified nodes of the
	 * blocking onto geometrical model entities.
	 * @param[in] ANodeIds the ids of the node we want to classify
	 * @param[in] ATolerance maximal distance to allow projections.
	 * @return the number of unclassified nodes that remain.
	 */
	int try_and_classify_nodes(std::set<TCellID> &ANodeIds, const double ATolerance = 0.1);

	/**@brief This methods classify edge eges onto the geometric model. It is designed to be called
	 * on a boundary edges only
	 * @param[out] AErrors 		list of errors done during the classification
	 */
	int try_and_capture(std::set<TCellID> &ANodeIds, std::set<TCellID> &AEdgeIds, std::set<TCellID> &AFaceIds);

	/**\brief return if the capt is possible
	 * * @param[in] AId an id of a not captured element, we want to split something to captured it
	 * @param[in] ADim the dim of the element not capt
	 * @return return true if the capt is possible, else, return false
	 */
	//bool check_capt_element(const int AId, const int ADim);

	/**@brief Split the sheet defined by edge @p AE
	 * @param[in] AnIdElement an id of a not captured element, we want to split something to captured it
	 * @param[in] ADim the dim of the element not captured
	 */
	//void capt_element(const int AnIdElement, const int ADim);

	/**\brief return if a cut is possible
	 * @param[in] pointId 		A point id
	 * @param[in] AllEdges 	all the edges of the blocking
	 * @return return true if a cut is possible, else, return false
	 */
	//bool check_cut_possible(int pointId, std::vector<std::vector<Blocking::Edge>> &AllEdges);

	/**@brief This methods classifies all the nodes onto the geometric model with snapping
	 * @param[out] AErrors 		list of errors done during the classification
	 * @param[in] AMaxDistance maximal distance to allow projections.
	 * @param[in] APointSnapDistance under this distance we collapse to the point
	 */
	void classify_nodes(ClassificationErrors &AErrors, const double AMaxDistance, const double APointSnapDistance);

	/**@brief This methods classify the cells of the blocking structure
	 * on the geometrical model. This classification process implies to
	 * project block nodes, edges, faces when mandatory.
	 * The @p AMaxDistance parameter is used to accept the projection. If the
	 * distance to the potential geometric cell to project on is greater than
	 * @p AMaxDistance, the projection is not done, and the cell can remain
	 * unclassified (geom_dim = 4). After projecting the node, a correction
	 * stage is done to check the distance to point. If this distance is lower
	 * than @p APointSnapDistance, the node is snapped onto the corresponding
	 * point
	 *
	 * @param[in] AMaxDistance maximal distance to allow projections
	 * @param[in] APointSnapDistance under this distance we collapse to the point
	 * @return list of errors. If empty, it means that the classification is
	 * a complete success.
	 */
	ClassificationErrors classify(const double AMaxDistance = 0.01, const double APointSnapDistance = 0.1);

	/**@brief This method colored all the faces with the same color by surfaces
	 * @return a map with the boundary faces colored
	 */
	std::map<Blocking::Face, int> blocking_color_faces();

	/**@brief This method check if the classification between the model and the blocking is valid (all the elements of the model are captured)
	 * @return If valid, return True, else, return False;
	 */
	bool checkValidity(ClassificationErrors &AErrors);

 private:
	/**@brief This method check if a 0-cell of the blocking structure is classified
	 * on the geometrical point @p AP.
	 * @param AP a geometrical point
	 * @return true and the node classified on @p AP, false otherwise and a random node
	 */
	std::pair<bool, Blocking::Node> find_node_classified_on(cad::GeomPoint *AP);

	/**@brief This method check if a 2-cell of the blocking structure is classified
	 * on the geometrical curve @s AS.
	 * @param AS a geometrical surface
	 * @return true and the faces classified on @s AS, false otherwise
	 */
	std::pair<bool, std::vector<Blocking::Face>> find_faces_classified_on(cad::GeomSurface *AS);

	/**@brief This methods try to capture all points of the geometric model. It is called internally
	 * vy the method *classify*.
	 * @param[out] AErrors 		list of errors done during the classification
	 * @param[in] AMaxDistance maximal distance to allow projections.
	 * @param[in] APointSnapDistance under this distance we collapse the node to the point
	 */
	void captured_points(ClassificationErrors &AErrors, const double AMaxDistance, const double APointSnapDistance);
	/**@brief This methods classify all eges onto the geometric model. It is called internally
	 * by the method *classify*.
	 * @param[out] AErrors 		list of errors done during the classification
	 */
	void classify_edges(ClassificationErrors &AErrors);

	/**@brief This methods classify all faces onto the geometric model. It is called internally
	 * by the method *classify*. We do this method after all points and curves are captured
	 * @param[out] AErrors 		list of errors done during the classification
	 */
	void classify_faces(ClassificationErrors &AErrors);
	/**@brief Generic method that gives among a collection of geometrical entities of same
	 * dimension, the cloest entity to point @p AP.
	 * @param[in] AP the point we consider
	 * @param[in] AGeomCells the list of geometrical cells we want to project @p AP on
	 * @return a tuple containing the distance to the closest cell in @p AGeomCells, the id of
	 * this cell and the projection of @p AP on @p AGeomCells
	 */
	std::tuple<double, int, math::Point> get_closest_cell(const math::Point &AP, const std::vector<cad::GeomEntity *> &AGeomCells);

	/**@brief This methods check if all boundary elements of a surface are captured.
	 * @return True if the boundary is captured
	 */
	bool boundary_surface_captured(cad::GeomSurface *AS);

 private:
	/**@brief This method check if there exists a series of 1-cells of the blocking
	 * that is classified on the geometrical curve @c AC, and capture it totally. We make the
	 * assumption that the curve is a simple connected 1-manifold. This method is a sub-case of
	 * method find_edges_classified_on_curve() where the curve has two end points.
	 * @param[in] AC a geometrical curve
	 * @param[in] ANodesOnAC nodes known as classified on AC
	 * @param[in] AEdgesOnAC edges known as classified on AC
	 * @return true and the series of edges classified on @c AC, false otherwise
	 */
	std::pair<bool, std::vector<Blocking::Edge>> find_edges_classified_on_curve_with_2_end_points(cad::GeomCurve *AC,
	                                                                                              std::vector<Blocking::Node> &ANodesOnAC,
	                                                                                              std::vector<Blocking::Edge> &AEdgesOnAC);

	/**@brief This method check if there exists a series of 1-cells of the blocking
	 * that is classified on the geometrical curve @c AC, and capture it totally. We make the
	 * assumption that the curve is a simple connected 1-manifold. This method is a sub-case of
	 * method find_edges_classified_on_curve() where the curve has only one end point1.
	 * @param[in] AC a geometrical curve
	 * @param[in] ANodesOnAC nodes known as classified on AC
	 * @param[in] AEdgesOnAC edges known as classified on AC
	 * @return true and the series of edges classified on @c AC, false otherwise
	 */
	std::pair<bool, std::vector<Blocking::Edge>> find_edges_classified_on_curve_with_1_end_point(cad::GeomCurve *AC,
	                                                                                             std::vector<Blocking::Node> &ANodesOnAC,
	                                                                                             std::vector<Blocking::Edge> &AEdgesOnAC);

	/**@brief This method check if there exists a series of 1-cells of the blocking
	 * that is classified on the geometrical curve @c AC, and capture it totally. We make the
	 * assumption that the curve is a simple connected 1-manifold. This method is a sub-case of
	 * method find_edges_classified_on_curve() where the curve has no end points.
	 * @param[in] AC a geometrical curve
	 * @param[in] ANodesOnAC nodes known as classified on AC
	 * @param[in] AEdgesOnAC edges known as classified on AC
	 * @return true and the series of edges classified on @c AC, false otherwise
	 */
	std::pair<bool, std::vector<Blocking::Edge>> find_edges_classified_on_curve_without_end_points(cad::GeomCurve *AC,
	                                                                                               std::vector<Blocking::Node> &ANodesOnAC,
	                                                                                               std::vector<Blocking::Edge> &AEdgesOnAC);
 bool build_edge_path_for_curve(cad::GeomCurve* AC, cad::GeomPoint* AP1, cad::GeomPoint* AP2,
	                               Blocking::Node& AN1, Blocking::Node& AN2,
	                               std::vector<TCellID>& ANodeIDs,
	                               std::vector<TCellID>& AEdgeIDs);

	 std::pair<bool, Blocking::Edge> find_aligned_edge(cad::GeomPoint* APoint,
											  const math::Vector3d& ATangent,
											  std::set<TCellID> &AEdgeIds);

 private:

	/*** the associated geometric model*/
	Blocking *m_blocking;
	/*** the associated geometric model*/
	cad::GeomManager *m_geom_model;
};
/*----------------------------------------------------------------------------*/
}     // namespace mctsblock
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MCTS_BLOCKING_CLASSIFIER_H
