/*----------------------------------------------------------------------------*/
#ifndef GMDS_CURVED_BLOCKING_CLASSIFIER_H
#	define GMDS_CURVED_BLOCKING_CLASSIFIER_H
/*----------------------------------------------------------------------------*/
#	include <LIB_GMDS_BLOCKING_export.h>
#	include <gmds/blocking/CurvedBlocking.h>
/*----------------------------------------------------------------------------*/
#	include <string>
#	include <type_traits>
#  include <tuple>

/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
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
/**@class CurvedBlockingClassifier
 * @brief Provide functions to check and update the classification of a curved
 * 		 blocking structure onto a geometrical model.
 *
 */
class LIB_GMDS_BLOCKING_API CurvedBlockingClassifier
{
 public:
	/** @brief Constructor that a curve blocking
	 * @param[in] ABlocking a blocking structure
	 */
	CurvedBlockingClassifier(CurvedBlocking* ABlocking);

	/** @brief  Destructor
	 */
	virtual ~CurvedBlockingClassifier();

	/**@brief This operation detects error in the classification.
	 * Errors are geometrical cells that are not catch by the
	 * blocking cells.
	 */
	ClassificationErrors detect_classification_errors(ClassificationErrors& AErrors);

	/**@brief this method clears the classification of the blocking
	 * structure. All the cells are then unclassified (geom_dim=4
	 * and geom_id=NullID)
	 *
	 */
	void clear_classification();
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
	ClassificationErrors classify(const double AMaxDistance=0.01, const double APointSnapDistance=0.1);

	/**@brief This method colored all the faces with the same color by surfaces
	 * @return a map with the boundary faces colored
	 */
	std::map<CurvedBlocking::Face,int> blocking_color_faces();

    /**@brief This method return all the possible cut
	 * @return a vector with only pair in, the first (pair.first) is the edge, and the second (pair.second) is the param to cut
	 */

    std::vector<std::pair<TCellID ,double>> list_Possible_Cuts();





 private:
	/**@brief This method check if a 0-cell of the blocking structure is classified
	 * on the geometrical point @p AP.
	 * @param AP a geometrical point
	 * @return true and the node classified on @p AP, false otherwise and a random node
	 */
	std::pair<bool, CurvedBlocking::Node> find_node_classified_on(cad::GeomPoint* AP);


	/**@brief This method check if a 1-cell of the blocking structure is classified
	 * on the geometrical curve @c AC.
	 * @param AC a geometrical curve
	 * @return true and the edges classified on @c AC, false otherwise
	 */
	std::pair<bool, std::vector<CurvedBlocking::Edge>> find_edge_classified_on(cad::GeomCurve* AC);

	/**@brief This method check if a 2-cell of the blocking structure is classified
	 * on the geometrical curve @s AS.
	 * @param AS a geometrical surface
	 * @return true and the faces classified on @s AS, false otherwise
	 */
	std::pair<bool, std::vector<CurvedBlocking::Face>> find_face_classified_on(cad::GeomSurface* AS);


	/**@brief This methods classify all nodes onto the geometric model. It is called internally
	 * vy the method *classify*.
	 * @param[out] AErrors 		list of errors done during the classification
	 * @param[in] AMaxDistance maximal distance to allow projections.
	 * @param[in] APointSnapDistance under this distance we collapse to the point
	 */
	void classify_nodes(ClassificationErrors& AErrors, const double AMaxDistance,const double APointSnapDistance);
	/**@brief This methods classify all eges onto the geometric model. It is called internally
	 * by the method *classify*.
	 * @param[out] AErrors 		list of errors done during the classification
	 */
	void classify_edges(ClassificationErrors& AErrors);

	/**@brief This methods classify all faces onto the geometric model. It is called internally
	 * by the method *classify*. We do this method after all points and curves are captured
	 * @param[out] AErrors 		list of errors done during the classification
	 */
	void classify_faces(ClassificationErrors& AErrors);
	/**@brief Generic method that gives among a collection of geometrical entities of same
	 * dimension, the cloest entity to point @p AP.
	 * @param[in] AP the point we consider
	 * @param[in] AGeomCells the list of geometrical cells we want to project @p AP on
	 * @return a tuple containing the distance to the closest cell in @p AGeomCells, the id of
	 * this cell and the projection of @p AP on @p AGeomCells
	 */
	std::tuple<double, int, math::Point> get_closest_cell(const math::Point& AP,
	                                         const std::vector<cad::GeomEntity*>& AGeomCells);

	/**@brief This methods check if all boundary elements of a surface are captured.
	 * @return True if the boundary is captured
	 */
	bool boundary_surface_captured(cad::GeomSurface* AS);




 private:
	/*** the associated geometric model*/
	CurvedBlocking* m_blocking;
	/*** the associated geometric model*/
	cad::GeomManager* m_geom_model;
};
/*----------------------------------------------------------------------------*/
}     // namespace blocking
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif// GMDS_BLOCKING_H
