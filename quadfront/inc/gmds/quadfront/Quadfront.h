#ifndef GMDS_QUADFRONT_QUADFRONT_H
#define GMDS_QUADFRONT_QUADFRONT_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_QUADFRONT_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace quadfront{
/*----------------------------------------------------------------------------*/
/** \class  dummy
 *  \brief  dummy class.
 */
class LIB_GMDS_QUADFRONT_API Quadfront{

 public:
	/*-------------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;
	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit Quadfront(const std::string& AFilename);
	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~Quadfront() =default;
	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	Quadfront::STATUS execute();

	/*-------------------------------------------------------------------------*/
	/** \brief  compute angle on a node between two boundary edge
	 */
	double angle(Node Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief set the number of the  Boundary of a edge
	 */
	void initialise_Boundary();
	/*-------------------------------------------------------------------------*/
	/** \brief set the number of the  Boundary of a edge
	 */
	void update_Boundary2(Edge& Aedge_base, Edge& Aedge_up);
	/*-------------------------------------------------------------------------*/
	/** \brief set in m_nodeBoundary the adjacent edge on the Boundary and the status of the node
	 */
	void set_nodeBoundary(Node& Anode, int& i, Edge& Aedge0, Edge& Aedge1);
	/*-------------------------------------------------------------------------*/
	/** \brief get from m_nodeBoundary the adjacent edge on the Boundary and the status of the node
	 */
	std::tuple<int, Edge, Edge> get_nodeBoundary(Node& Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief Test if a node is on a boundary
	 */
	bool find_nodeBoundary(Node& Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief get the map of adjacents edges of a node inside a face
	 */
	std::map<Node, std::pair<Edge,Edge>> get_edgeInFace(Face Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief get the opposit edge of a node inside a face
	 */

	Edge get_opposit2Node(Node& Anode, Face& Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief get the opposit node of an edge inside a face
	 */
	Node get_opposit2Edge(Edge& Aedge, Face& Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief convert Edge to vector
	 */
	math::Vector3d node2Vector(Edge& Aedge, Node& Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief  get the bissectrice of two vectors
	 */
	math::Vector3d get_bissectrice(Node& Anode, Edge& Aedge0, Edge& Aedge1);
	/*-------------------------------------------------------------------------*/
	/** \brief  get the intersection bet
	 */
	math::Vector3d  intersectionVec2Edge(math::Vector3d& vec, Node& Anode, Edge& Aedge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge  sideEdge(Node& Anode, Edge& Aedge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void removeFace(Face& Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void removeEdge(Edge& Aedge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge swap(Node &Anode, Face& Aface, Face& Aface_opposit);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge split(Node& AnewNode, Node& Anode, Face& Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge split_operation(Node &AnewNode, Node &Anode, Face& Aface, Face& Aface_opposit);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	std::vector<Edge>  interestionS(Node& AnodeRecovery0, Node& AnodeRecovery1);
	/*-------------------------------------------------------------------------*/
	/** \brief  create Quad
	 */
	Edge recoveryEdge(std::vector<Edge>  listInside, Node& nodeRecovery0, Node& nodeRecovery1);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void testInside(Face& faceRef, Face& faceInside, std::vector<Face>& listInside);
	/*-------------------------------------------------------------------------*/
	/** \brief  create Quad
	 */
	Face createQuad(Node& Anode, Edge& Aedge, Edge& sideEdge0, Edge& sideEdge1, Edge& recoveryEdge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge closingSeam(Edge& Aedge0, Edge& Aedge1, Edge& Aedge_triangle0, Edge& Aedge_triangle1);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge transitionSeam(Node& Anode, Edge& Aedge_max, Edge& Aedge_min);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge transitionSplit(Node& Anode, Edge& Aedge_max, Edge& Aedge_min);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void seam(Node& Anode, std::vector<Node>& listNode);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Face qmorphAlgo(Edge& Aedge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Node get_oppositNode_Quad(Node &Anode, Face &Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge get_oppositEdge_Quad(Edge &Aedge, Face &Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	std::vector<double> statFront(std::vector<Edge> listEdge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void blSmoothing(Node& Anode, double average, double tr);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void qmorphSmoothing(Node& Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void interiorSmoothing(Node& Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void triangleSmoothing(Edge Aedge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void initCorners(Node& Anode, std::vector<Face>& listQuad);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void petitAngle(Face& Aface, std::vector<Face>& listFace, std::vector<Edge>& listEdge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void initEnds(Node& Anode, std::vector<Face>& listQuad);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void initRevers(Node& Anode, std::vector<Face>& listQuad);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void initialiseQuad(Node& Anode, std::vector<Face>& listQuad);
	/*-------------------------------------------------------------------------*/

	Mesh* m_mesh;
	//List of gmds Node on the boundary
	std::map<Node, std::pair<Edge, Edge>> m_nodeBoundary;
	//List of gmds edges on the front
	std::vector<std::vector<Edge>> m_listFront;
	//List of gmds nodes on the front
	std::map<Node, std::tuple<int, Edge, Edge>> m_nodeFront;

 protected:



};

/*----------------------------------------------------------------------------*/
}  // end namespace quadfront
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // GMDS_QUADFRONT_QUADFRONT_H