/*----------------------------------------------------------------------------*/
#ifndef GMDS_CURVED_BLOCKING_H
#	define GMDS_CURVED_BLOCKING_H
/*----------------------------------------------------------------------------*/
#	include <CGAL/Cell_attribute.h>
#	include <CGAL/Generalized_map.h>
#	include <LIB_GMDS_BLOCKING_export.h>
#	include <gmds/math/Point.h>
#	include <gmds/utils/CommonTypes.h>
#	include <gmds/utils/Exception.h>
# include <gmds/cad/GeomManager.h>
#	include <string>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace blocking {
/*----------------------------------------------------------------------------*/
/**@struct CellInfo
 * @brief This structure gather the pieces of data that are shared by any
 * 		 blocking cell. Each cell is defined by:
* 		 	- its dimension @p topo_dim, which is 0 for a node, 1 for an edge, 2
* 		     for a face, and 3 for a block
* 		   - its id @p topo_id, which is unique. Each time a cell is created, the
 * 		  blocking structure assigns an id to it.
* 		   - The data relative to the geometric cell it is classified on. A
 * 		  geometrical cell is defined by the couple (@p geom_dim, @p geom_id).
 * 		  A value of 4 for @p geom_dim means that the block cell is not classified.
 *
 */
struct CellInfo {
	/*** dimension of the topological cell */
	int topo_dim;
	/*** unique id of the topological cell */
	int topo_id;
	/*** dimension of the geometrical cell we are classifid on */
	int geom_dim;
	/*** unique id of the geomtrical cell */
	int geom_id;
	/** @brief Constructor
	 * @param ATopoDim Cell dimension
	 * @param ATopoId  Cell unique id
	 * @param AGeomDim on-classify geometric cell dimension (4 if not classified)
	 * @param AGeomId on-classify geometric cell unique id
	 */
	CellInfo(const int ATopoDim = 4, const int ATopoId = NullID, const int AGeomDim = 4, const int AGeomId = NullID) :
	  topo_dim(ATopoDim), topo_id(ATopoId), geom_dim(AGeomDim), geom_id(AGeomId)
	{}
};
/*----------------------------------------------------------------------------*/
/**@struct NodeInfo
 * @brief Specific structure for describing a node. It extends the "CellInfo"
 *        structure by giving a location (x,y,z).
 */
struct NodeInfo:CellInfo{
	/*** node location in space, i.e. a single point */
	math::Point point;
	/** @brief Constructor
	 * @param ATopoId  Node unique id
	 * @param AGeomDim on-classify geometric cell dimension (4 if not classified)
	 * @param AGeomId  on-classify geometric cell unique id
	 */
	NodeInfo(const int ATopoId = NullID,
	         const int AGeomDim = 4,
	         const int AGeomId = NullID,
	         const math::Point &APoint = math::Point(0, 0, 0)) :
	  CellInfo(0,ATopoId,AGeomDim,AGeomId), point(APoint)
	{}
};
/*----------------------------------------------------------------------------*/
/** When we merge two cells*/
struct MergeFunctor
{
	template<class Cell_attribute> void operator()(Cell_attribute &ca1, Cell_attribute &ca2)
	{
		if (ca1.info().geom_dim == ca2.info().geom_dim) {
			// the cells are classifed on the same dim geom entity
			if (ca1.info().geom_id == ca2.info().geom_id) {
				ca1.info().geom_dim = ca2.info().geom_dim;
			}
			else
				throw GMDSException("Classification error!!!");
		}
		else if (ca1.info().geom_dim < ca2.info().geom_dim) {
			// the cells are classifed on the same dim geom entity
			ca1.info().geom_dim = ca1.info().geom_dim;
			ca1.info().geom_id = ca1.info().geom_id;
		}
		else {     // third case: ca1.info().dim>ca2.info().dim
			ca1.info().geom_dim = ca2.info().geom_dim;
			ca1.info().geom_id = ca2.info().geom_id;
		}
	}
};
/*----------------------------------------------------------------------------*/
struct MergeFunctorNode
{
	template<class Cell_attribute> void operator()(Cell_attribute &ca1, Cell_attribute &ca2)
	{
		if (ca1.info().geom_dim == ca2.info().geom_dim) {
			// the cells are classifed on the same dim geom entity
			if (ca1.info().geom_id == ca2.info().geom_id) {
				ca1.info().point = 0.5 * (ca1.info().point + ca2.info().point);
				// TODO: add the projection stage on the geometric entity
			}
			else
				throw GMDSException("Classification error!!!");
		}
		else if (ca1.info().geom_dim < ca2.info().geom_dim) {
			// the cells are classifed on the same dim geom entity
			ca1.info().geom_dim = ca1.info().geom_dim;
			ca1.info().geom_id = ca1.info().geom_id;
			ca1.info().point = ca1.info().point;
		}
		else {     // third case: ca1.info().dim>ca2.info().dim
			ca1.info().geom_dim = ca2.info().geom_dim;
			ca1.info().geom_id = ca2.info().geom_id;
			ca1.info().point = ca2.info().point;
		}
	}
};
/*----------------------------------------------------------------------------*/
struct SplitFunctor
{
	template<class Cell_attribute> void operator()(Cell_attribute &ca1, Cell_attribute &ca2)
	{
		ca1.info().geom_dim = ca1.info().geom_dim;
		ca1.info().geom_id = ca1.info().geom_id;
		ca2.info() = ca1.info();
	}
};
/*----------------------------------------------------------------------------*/
struct SplitFunctorNode
{
	template<class Cell_attribute> void operator()(Cell_attribute &ca1, Cell_attribute &ca2)
	{
		ca2.info().geom_dim = ca1.info().geom_dim;
		ca2.info().geom_id = ca1.info().geom_id;
		ca2.info().point = ca1.info().point;
	}
};
/*----------------------------------------------------------------------------*/
struct CellData
{
	template<class GMap> struct Dart_wrapper
	{
		typedef CGAL::Cell_attribute<GMap, NodeInfo, CGAL::Tag_true, MergeFunctorNode, SplitFunctorNode> Node_attribute;     // A weight
		typedef CGAL::Cell_attribute<GMap, CellInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor> Edge_attribute;             // A weight
		typedef CGAL::Cell_attribute<GMap, CellInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor> Face_attribute;             // A weight
		typedef CGAL::Cell_attribute<GMap, CellInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor> Block_attribute;            // A weight
		typedef std::tuple<Node_attribute, Edge_attribute, Face_attribute, Block_attribute> Attributes;
	};
};
// Definition of my generalized map.
typedef CGAL::Generalized_map<3, CellData> GMap3;
typedef GMap3::Dart_descriptor DartDescriptor;

/*----------------------------------------------------------------------------*/
/**@class CurvedBlocking
 * @brief Provide a curved blocking data structure using the 3-G-Map model
 * 		 as described and provided by CGAL.
 * 		 (see https://doc.cgal.org/latest/Generalized_map/index.html)
 */
class LIB_GMDS_BLOCKING_API CurvedBlocking
{
 public:
	typedef GMap3::Attribute_handle<3>::type Block;
	typedef GMap3::Attribute_handle<2>::type Face;
	typedef GMap3::Attribute_handle<1>::type Edge;
	typedef GMap3::Attribute_handle<0>::type Node;
	/** @brief  Default Constructor
	 */
	CurvedBlocking();

	/** @brief  Destructor
	 */
	virtual ~CurvedBlocking();

	/** Create a single hexahedral block in the blocking structure
	 * @return The created block
	 */
	Block createBlock(math::Point &AP1, math::Point &AP2,
	                  math::Point &AP3, math::Point &AP4,
	                  math::Point &AP5, math::Point &AP6,
	                  math::Point &AP7, math::Point &AP8);

	/** Get all the faces of a block. If it is a hexahedral block,
	 * we have 6 faces, the first and the second are opposite, idem
	 * for the third and fourth, and the fifth and sixth.
	 * @param AB a block
	 * @return the set of faces of the block.
	 */
	std::vector<Face> get_faces_of_block(const Block AB);
	/** Get all the nodes of a block. If it is a hexahedral block,
	 * we have 8 nodes, given as usual in gmds
	 * @param AB a block
	 * @return the set of nodes of the block.
	 * */
	std::vector<Node> get_nodes_of_block(const Block AB);
	/** Get all the nodes of a face (ordered).
	 * @param AF a face
	 * @return the set of nodes of the face.
	 * */
	std::vector<Node> get_nodes_of_face(const Face AF);
	/** Return the face center
	 * @param AF a face
	 * @return the center point of AF
	 * */
	math::Point get_center_of_face(const Face AF);
	/** Return the face center
	 * @param AF a face
	 * @return a point which is the center of AF
	 * */
	math::Point get_center_of_block(const Block AB);

	void sew(DartDescriptor d1, DartDescriptor d2)
	{
		m_gmap.sew<3>(d1, d2);
	}

	std::string info() const;

	bool isValidTopology() const
	{
		return m_gmap.is_valid();
	}

 private:
	Node createNode(const int AGeomDim, const int AGeomId, math::Point &APoint);
	Edge createEdge(const int AGeomDim, const int AGeomId);
	Face createFace(const int AGeomDim, const int AGeomId);
	Block createBlock(const int AGeomDim, const int AGeomId);
	GMap3 m_gmap;

	static int m_counter_nodes;
	static int m_counter_edges;
	static int m_counter_faces;
	static int m_counter_blocks;
};
/*----------------------------------------------------------------------------*/
}     // namespace blocking
/*----------------------------------------------------------------------------*/
}     // namespace gmds

/*----------------------------------------------------------------------------*/
#endif     // GMDS_BLOCKING_H
/*----------------------------------------------------------------------------*/