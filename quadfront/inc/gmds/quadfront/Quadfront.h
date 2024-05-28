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
	/** \brief set the number of the  Boundary of a edge
	 */
	void set_edgeBoundary(Edge& Aedge, int i );
	/*-------------------------------------------------------------------------*/
	/** \brief get the number of the Boundary of a edge
	 */
	int get_edgeBoundary(Edge& edge);
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
	/** \brief Test if a edge is on a boundary
	 */
	bool find_edgeBoundary(Edge& Aedge);
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
	static Node get_opposit2Edge(Edge& Aedge, Face& Aface);
	/*-------------------------------------------------------------------------*/
	/** \brief convert Edge to vector
	 */
	math::Vector3d node2Vector(Edge& Aedge, Node& Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief  compute angle on a node between two boundary edge
	 */
	double angle(Node Anode);
	/*-------------------------------------------------------------------------*/
	/** \brief  get the bissectrice of two vectors
	 */
	math::Vector3d get_bissectrice(math::Vector3d Avec1, math::Vector3d Avec2);
	/*-------------------------------------------------------------------------*/
	/** \brief set status of a node in m_nodeBoundary
	 */
	void initialise_nodeBoundary();
	/*-------------------------------------------------------------------------*/
	/** \brief  get the intersection bet
	 */
	 std::tuple<double, double, double> intersectionVec2Edge(math::Vector3d& vec, Node& Anode, Edge& Aedge);
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
	Edge swap(Edge& Aedge_inter,
				 Node& Anode,
				 Node& Anode_opposit,
				 std::map<Node, std::pair<Edge, Edge>>& node2edge_opposit,
				 Edge& Aedge_min1,
				 Edge& Aedge_min2);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	Edge split(Node& AnewNode,
	           Edge& Aedge_inter,
	           Node& Anode,
	           Node& Anode_opposit,
	           std::map<Node, std::pair<Edge, Edge>>& node2edge_opposit,
	           Edge& Aedge_min1,
	           Edge& Aedge_min2);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	std::vector<Face> get_superposed(Node& Anode, Edge& Aedge, Edge& sideEdge0, Edge& sideEdge1, Edge& upEdge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	std::vector<Edge>  interestionS(Node& AnodeRecovery0, Node& AnodeRecovery1);
	/*-------------------------------------------------------------------------*/
	/** \brief  create Quad
	 */
	void testInside(Face& faceRef, Face& faceInside, std::vector<Face>& listInside);
	/*-------------------------------------------------------------------------*/
	/** \brief  create Quad
	 */
	std::vector<Face> createQuad(Node& Anode, Edge& Aedge, Edge& sideEdge0, Edge& sideEdge1, Edge& recoveryEdge);
	/*-------------------------------------------------------------------------*/
	/** \brief  Status code for executing algorithms
	 */
	void  edgeRecovery(Edge Aedge);
	/*-------------------------------------------------------------------------*/

	Mesh* m_mesh;
	//List of gmds edges on the boundary
	std::map<Edge, int> m_edgeBoundary;
	//List of gmds nodes on the boundary
	std::map<Node, std::tuple<int, Edge, Edge>> m_nodeBoundary;

 protected:



};

/*----------------------------------------------------------------------------*/
}  // end namespace quadfront
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // GMDS_QUADFRONT_QUADFRONT_H