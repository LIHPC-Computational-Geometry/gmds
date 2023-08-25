/*----------------------------------------------------------------------------*/
#ifndef GMDS_CURVED_BLOCKING_H
#define GMDS_CURVED_BLOCKING_H
/*----------------------------------------------------------------------------*/
#include <CGAL/Cell_attribute.h>
#include <CGAL/Generalized_map.h>
#include <LIB_GMDS_BLOCKING_export.h>
#include <gmds/cad/GeomManager.h>
#include <gmds/ig/Mesh.h>
#include <gmds/math/Point.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#include <string>
#include <tuple>
#include <type_traits>

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
 */
struct CellInfo
{
	/*** dimension of the topological cell */
	int topo_dim;
	/*** unique id of the topological cell */
	int topo_id;
	/*** dimension of the geometrical cell we are classifid on */
	int geom_dim;
	/*** unique id of the geomtrical cell */
	int geom_id;
	/*** global counter used to assign an unique id to each block */
	static int m_counter_global_id;

	/** @brief Constructor
	 * @param ATopoDim Cell dimension
	 * @param AGeomDim on-classify geometric cell dimension (4 if not classified)
	 * @param AGeomId on-classify geometric cell unique id
	 */
	CellInfo(const int ATopoDim = 4, const int AGeomDim = 4, const int AGeomId = NullID) :
	  topo_dim(ATopoDim), topo_id(m_counter_global_id++), geom_dim(AGeomDim), geom_id(AGeomId)
	{
	}
};
/*----------------------------------------------------------------------------*/
/**@struct NodeInfo
 * @brief Specific structure for describing a node. It extends the "CellInfo"
 *        structure by giving a location (x,y,z).
 */
struct NodeInfo : CellInfo
{
	/*** node location in space, i.e. a single point */
	math::Point point;
	/** @brief Constructor
	 * @param AGeomDim on-classify geometric cell dimension (4 if not classified)
	 * @param AGeomId  on-classify geometric cell unique id
	 */
	NodeInfo(const int AGeomDim = 4, const int AGeomId = NullID, const math::Point &APoint = math::Point(0, 0, 0)) :
	  CellInfo(0, AGeomDim, AGeomId), point(APoint)
	{
	}
};
/*----------------------------------------------------------------------------*/
/** @struct MergeFunctor
 * @brief This structure provides a function for the Gmap class in order to
 * handle attributes when i-cells are merged. For instance, when two blocks B1
 * and B2 are glued along the face F1 for B1 and F2 for B2, those faces are merged
 * into a single one. Attributes of faces F1 and F2 are so merged too (it will
 * also be the case of the nodes and edges that are incident to faces F1 and F2).
 * common face.
 *
 * This functor defines the behavior for the attributes of Edges, Faces and Blocks.
 * More specifically, when two i-cells are merged, we consider the first one as the
 * master and the second one as the slave. The resulting i-cell keep the topological
 * data of the first one.
 *
 * For the geometrical data, we look at the lowest dimension of the associated geometrical
 * cell. We are always constraint to preserving classification of lowest dimensional cell.
 * For instance, if an edge E1, classified on a Surface, is merged with an edge E2, classified
 * on a curve, we preserve the curve classification.
 */
struct MergeFunctor
{
	/** @brief Merging function for general cells (not a node).
	 * @tparam Cell_attribute
	 * @param[in] ACA1 the first attribute
	 * @param[in] ACA2 the second attribute
	 */
	template<class Cell_attribute> void operator()(Cell_attribute &ACA1, Cell_attribute &ACA2)
	{
		if (ACA1.info().geom_dim == ACA2.info().geom_dim) {
			// the cells are classifed on the same dim geom entity
			if (ACA1.info().geom_id == ACA2.info().geom_id) {
				ACA1.info().geom_dim = ACA2.info().geom_dim;
			}
			else
				throw GMDSException("Classification error!!!");
		}
		else if (ACA1.info().geom_dim < ACA2.info().geom_dim) {
			// the cells are classifed on the same dim geom entity
			ACA1.info().geom_dim = ACA1.info().geom_dim;
			ACA1.info().geom_id = ACA1.info().geom_id;
		}
		else {     // third case: ca1.info().dim>ca2.info().dim
			ACA1.info().geom_dim = ACA2.info().geom_dim;
			ACA1.info().geom_id = ACA2.info().geom_id;
		}
		// nothing to do for the topological information of
	}
};
/*----------------------------------------------------------------------------*/
/** @struct MergeFunctorNode
 * @brief This structure provides a function for the Gmap class in order to
 * handle attributes when O-cells are merged. For instance, when two blocks B1
 * and B2 are glued along the face F1 for B1 and F2 for B2, those faces are merged
 * into a single one. Incident nodes of F1 and F2 are then merged too.
 *
 * This functor defines the behavior when merging the attributes of nodes.
 *  The resulting attribute keep the topological data of the first one.
 *
 * For the geometrical data, we look at the lowest dimension of the associated geometrical
 * cell. We are always constraint to preserving classification of lowest dimensional cell.
 * For instance, if a node N1, classified on a point, is merged with a node N2, classified
 * on a surface, we preserve the point classification.
 *
 * The node location is also impacted. If the merged nodes are on the same geometrical
 * cell, we take the average location on the geometrical entity. Otherwise, we keep the
 * most constrained location
 */
struct MergeFunctorNode
{

	/** @brief Merging function for nodes
	 * @tparam Cell_attribute
	 * @param[in] ACA1 the first node attribute
	 * @param[in] ACA2 the second node attribute
	 */
	template<class Cell_attribute> void operator()(Cell_attribute &ACA1, Cell_attribute &ACA2)
	{
		if (ACA1.info().geom_dim == ACA2.info().geom_dim) {
			if (ACA1.info().geom_id == ACA2.info().geom_id) {
				// nodes are clasified on the same geometrical entity!
				ACA1.info().point = 0.5 * (ACA1.info().point + ACA2.info().point);
				// TODO: add the projection stage on the geometric entity
			}
			else
				throw GMDSException("Classification error!!!");
		}
		else if (ACA1.info().geom_dim < ACA2.info().geom_dim) {
			// the cells are classifed on the same dim geom entity
			ACA1.info().geom_dim = ACA1.info().geom_dim;
			ACA1.info().geom_id = ACA1.info().geom_id;
			ACA1.info().point = ACA1.info().point;
		}
		else {     // third case: ca1.info().dim>ca2.info().dim
			ACA1.info().geom_dim = ACA2.info().geom_dim;
			ACA1.info().geom_id = ACA2.info().geom_id;
			ACA1.info().point = ACA2.info().point;
		}
	}
};
/*----------------------------------------------------------------------------*/
/** @struct SplitFunctor
 * @brief This structure provides a function for the Gmap class in order to
 * handle attributes when i-cells are split. It simply consists in duplicating
 * the attributes.
 */
struct SplitFunctor
{

	/** @brief Splitting function for general cells (including nodes).
	 * @tparam Cell_attribute
	 * @param[in] ACA1 the first cell attribute
	 * @param[in] ACA2 the second cell attribute
	 */
	template<class Cell_attribute> void operator()(Cell_attribute &ca1, Cell_attribute &ca2)
	{
		ca1.info().geom_dim = ca1.info().geom_dim;
		ca1.info().geom_id = ca1.info().geom_id;

		ca2.info().geom_dim = ca1.info().geom_dim;
		ca2.info().geom_id = ca1.info().geom_id;
		ca2.info().topo_dim = ca1.info().topo_dim;
		ca2.info().topo_id = CellInfo::m_counter_global_id++;

	}
};

/*----------------------------------------------------------------------------*/
/** @struct SplitFunctorNode
 * @brief This structure provides a function for the Gmap class in order to
 * handle attributes when nodes are split. It simply consists in duplicating
 * the attributes.
 */
struct SplitFunctorNode
{
	template<class Cell_attribute> void operator()(Cell_attribute &ca1, Cell_attribute &ca2)
	{
		ca2.info().geom_dim = ca1.info().geom_dim;
		ca2.info().geom_id = ca1.info().geom_id;
		ca2.info().point = ca1.info().point;
		ca2.info().topo_dim = ca1.info().topo_dim;
		ca2.info().topo_id = CellInfo::m_counter_global_id++;
	}
};
/*----------------------------------------------------------------------------*/
/**@struct CellData
 * @brief This structure provides the expected pieces of information that are
 * require by the CGAL N-G-Map implementation in order to automatically
 * handle attributes during n-G-map modification operations.
 */
struct CellData
{
	template<class GMap> struct Dart_wrapper
	{
		using Node_attribute = CGAL::Cell_attribute<GMap, NodeInfo, CGAL::Tag_true, MergeFunctorNode, SplitFunctorNode>;
		using Edge_attribute = CGAL::Cell_attribute<GMap, CellInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor>;
		using Face_attribute = CGAL::Cell_attribute<GMap, CellInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor>;
		using Block_attribute = CGAL::Cell_attribute<GMap, CellInfo, CGAL::Tag_true, MergeFunctor, SplitFunctor>;
		/** ordered tuple of attributes to indicate to the GMap what attribute corresponds to each cell.*/
		using Attributes = std::tuple<Node_attribute, Edge_attribute, Face_attribute, Block_attribute>;
	};
};
/*----------------------------------------------------------------------------*/
/** Definition of my generalized map.*/
using GMap3 = CGAL::Generalized_map<3, CellData>;
/** Definition of the dart type.*/
using Dart3 = GMap3::Dart_handle;
/*----------------------------------------------------------------------------*/
/**@class CurvedBlocking
 * @brief Provide a curved blocking data structure using the 3-G-Map model
 * 		 as described and provided by CGAL.
 * 		 (see https://doc.cgal.org/latest/Generalized_map/index.html)
 */
class LIB_GMDS_BLOCKING_API CurvedBlocking
{
 public:
	/** Inner type to see i-attributes as the respective block cells whe want
	 * to handle in practice (considering a cellular view of the block structure
	 * and not its underlying topological maps).
	 */
	using Block = GMap3::Attribute_handle<3>::type;
	using Face = GMap3::Attribute_handle<2>::type;
	using Edge = GMap3::Attribute_handle<1>::type;
	using Node = GMap3::Attribute_handle<0>::type;
	/** @brief Constructor that takes a geom model as an input. A
	 * blocking is always used for partitioning a geometric domain.
	 * @param[in] AGeomModel the geometric model we want to block
	 * @param[in] AInitAsBoundingBox indicates that the block structure must remain
	 * empty (false) or be initialized as the bounding box of @p AGeomModel
	 */
	CurvedBlocking(cad::GeomManager *AGeomModel, bool AInitAsBoundingBox = false);

	/** @brief  Destructor
	 */
	virtual ~CurvedBlocking();
	/**@brief gives access to the underlying gmap structure
	 * @return the internal 3-G-map
	 */
	GMap3 *gmap();
	/**@brief gives access to the associated geom model
	 * @return the internal geom model
	 */
	cad::GeomManager *geom_model();
	/**@brief Gives the number of @p TDim-cells in the blocking structure.
	 * 		 @p TDIM must be comprised in [0,3].
	 *
	 * @tparam TDim the dimension of cells we want the number
	 * @return the number of @p TDim-cells in the block structure
	 */
	template<int TDim> int get_nb_cells() const
	{
		return m_gmap.number_of_attributes<TDim>();
	}
	/** Create a single hexahedral block in the blocking structure
	 * @return The created block
	 */
	Block create_block(
	   const math::Point &AP1,
	   const math::Point &AP2,
	   const math::Point &AP3,
	   const math::Point &AP4,
	   const math::Point &AP5,
	   const math::Point &AP6,
	   const math::Point &AP7,
	   const math::Point &AP8);
	/** Removes the block @AB from the structure
	 * @param[in] AB the block to remove
	 */
	void remove_block(Block AB);

	/** Return the info for the node of id @p ANodeId
	 * @param[in] ANodeId topological node id
	 * @return a tuple where the first parameter is the geom_dim, the second its geom_id, and the third is location
	 */
	std::tuple<int, int, math::Point> get_node_info(const int ANodeId);
	/** Return the info for the edge of id @p AEdgeId
	 * @param[in] AEdgeId topological edge id
	 * @return a tuple where the first parameter is the geom_dim, the second its geom_id
	 */
	std::tuple<int, int> get_edge_info(const int AEdgeId);
	/** Return the info for the face of id @p AFaceId
	 * @param[in] AFaceId topological face id
	 * @return a tuple where the first parameter is the geom_dim, the second its geom_id
	 */
	std::tuple<int, int> get_face_info(const int AEdgeId);
	/** Return the info for the block of id @p ABlockId
	 * @param[in] ABlockId topological block id
	 * @return a tuple where the first parameter is the geom_dim, the second its geom_id
	 */
	std::tuple<int, int> get_block_info(const int ABlockId);
	/**@brief Non-optimal method to get all the blocks of the structure. The
	 * best option is to traverse the block structure through the gmap
	 * structure (iterators on attributes)
	 * @return a vector of blocks
	 */
	std::vector<Block> get_all_blocks();
	std::vector<TCellID> get_all_id_blocks();
	/**@brief Non-optimal method to get all the faces of the structure. The
	 * best option is to traverse the block structure through the gmap
	 * structure (iterators on attributes)
	 * @return a vector of faces
	 */
	std::vector<Face> get_all_faces();
	std::vector<TCellID> get_all_id_faces();
	/**@brief Non-optimal method to get all the edges of the structure. The
	 * best option is to traverse the block structure through the gmap
	 * structure (iterators on attributes)
	 * @return a vector of edges
	 */
	std::vector<Edge> get_all_edges();
	std::vector<TCellID> get_all_id_edges();
	/**@brief Non-optimal method to get all the nodes of the structure. The
	 * best option is to traverse the block structure through the gmap
	 * structure (iterators on attributes)
	 * @return a vector of nodes
	 */
	std::vector<Node> get_all_nodes();
	std::vector<TCellID> get_all_id_nodes();
	/**@brief moves node @p AN towards the expected new location @p ALoc.
	 * If @p AN is classified onto a geometrical cell, the node @p AN
	 * is first moved to @p ALoc, then it is projected onto the
	 * geometrical cell it is classified on.
	 *
	 * @param AN the node to move
	 * @param ALoc the new node location
	 */
	void move_node(Node AN, math::Point &ALoc);
	/** Get all the edges adjacent to a node
	 * @param[in] AN a node
	 * @return the set of edges adjacent to the node.
	 */
	std::vector<Edge> get_edges_of_node(const Node AN);
	/** Get all the faces adjacent to a node
	 * @param[in] AN a node
	 * @return the set of faces adjacent to the node.
	 */
	std::vector<Face> get_faces_of_node(const Node AN);
	/** Get all the blocks adjacent to a node
	 * @param[in] AN a node
	 * @return the set of blocks adjacent to the node.
	 */
	std::vector<Block> get_blocks_of_node(const Node AN);
	/** Get all the faces adjacent to an edge
	 * @param[in] AE an edge
	 * @return the set of faces adjacent to the edge.
	 */
	std::vector<Face> get_faces_of_edge(const Edge AE);
	/** Get all the blocks adjacent to an edge
	 * @param[in] AE an edge
	 * @return the set of blocks adjacent to the edge.
	 */
	std::vector<Block> get_blocks_of_edge(const Edge AE);
	/** Get all the edges adjacent to a face
	 * @param[in] AF a face
	 * @return the set of edges adjacent to the face.
	 */
	std::vector<Edge> get_edges_of_face(const Face AF);
	/** Get all the blocks adjacent to a face
	 * @param[in] AF a face
	 * @return the set of blocks adjacent to the face.
	 */
	std::vector<Block> get_blocks_of_face(const Face AF);

	/** Get all the faces of a block. If it is a hexahedral block,
	 * we have 6 faces, the first and the second are opposite, idem
	 * for the third and fourth, and the fifth and sixth.
	 * @param[in] AB a block
	 * @return the set of faces of the block.
	 */
	std::vector<Face> get_faces_of_block(const Block AB);
	/** Get all the edges of a block. If it is a hexahedral block,
	 * we have 12 faces.
	 * @param[in] AB a block
	 * @return the set of edges of the block.
	 */
	std::vector<Edge> get_edges_of_block(const Block AB);
	/** Get all the nodes of a block. If it is a hexahedral block,
	 * we have 8 nodes, given as usual in gmds
	 * @param[in] AB a block
	 * @return the set of nodes of the block.
	 */
	std::vector<Node> get_nodes_of_block(const Block AB);
	/** Get all the nodes of a face (ordered).
	 * @param[in] AF a face
	 * @return the set of nodes of the face.
	 */
	std::vector<Node> get_nodes_of_face(const Face AF);
	/** Get the ending nodes of an edge.
	 * @param[in] AE an edge
	 * @return the ending nodes of the edge
	 */
	std::vector<Node> get_nodes_of_edge(const Edge AE);
	/** Return the face center
	 * @param AF a face
	 * @return the center point of @p AF
	 */
	math::Point get_center_of_face(const Face AF);
	/** Return the face normal. No orientation is provided the
	 *  normal can be in any direction.
	 * @param AF a face
	 * @return the normal to @p AF.
	 */
	math::Vector3d get_normal_of_face(const Face AF);
	/** Return the edge center
	 * @param AE an edge
	 * @return the center point of @p AE
	 */
	math::Point get_center_of_edge(const Edge AE);
	/** Return the block center
	 * @param[in] AB a block
	 * @return a point which is the center of @p AB
	 */
	math::Point get_center_of_block(const Block AB);
	/**@brief Get all the parallel edges composing the sheet defined from edge @p AE
	 * @param[in]  AE			the edge we start from
	 * @param[out] AEdges	all the edges of the sheet defined by @p AE
	 */
	void get_all_sheet_edges(const Edge AE, std::vector<Edge> &AEdges);

	/**@brief Get all the parallel edges composing sheets. Each set of parallel
	 * 		 edges is an item pf @p ASheetEdges
	 * return all the edges gathered by sheet
	 */
	std::vector<std::vector<Edge>> get_all_sheet_edge_sets();
	/**@brief Get one dart per  parallel edges composing the sheet defined from edge
	 * @p AE. All the returned darts are on the same side of each edge.
	 * @param[in]  AE			the edge we start from
	 * @param[out] ADarts	one dart per edge of the sheet defined by @p AE
	 */
	void get_all_sheet_darts(const Edge AE, std::vector<Dart3> &ADarts);
	/**@brief Get one dart per parallel face composing the chord defined from face
	 * @p AF. All the returned darts belong to different blocks.
	 * @param[in]  AF			the face we start from
	 * @param[out] ADarts	one dart per face of the sheet defined by @p AE
	 */
	void get_all_chord_darts(const Face AF, std::vector<Dart3> &ADarts);
	/**@brief Split the sheet defined by edge @p AE at the closest position of @p AP.
	 * 		 The method project @p AP on @p AE. It gives us a cut ratio, we use on
	 * 		 all the parallel edges.
	 * @param[in] AE an edge we want to split in two edges
	 * @param[in] AP a point we use to define where AE must be cut.
	 */
	void cut_sheet(const Edge AE, const math::Point &AP);
	/**@brief Split the sheet defined by edge @p AE at the parameter @p AParam, which is included
	 * 		 in ]0,1[. The first end point of @p AE is at parameter 0, the second one at parameter 1.
	 *
	 * @param[in] AE an edge we want to split in two edges
	 * @param[in] AParam a parameter included in ]0,1[
	 */
	void cut_sheet(const Edge AE, const double AParam);
	/**@brief Split the sheet defined by edge @p AE
	 * @param[in] AE an edge we want to split in two edges
	 */
	void cut_sheet(const Edge AE);

	/**@brief Low level operation that @p TDim-sew two darts
	 * @tparam TDim sewing dimension
	 * @param[in] AD1 First dart
	 * @param[in] AD2 Second dart
	 */
	template<int TDim> void sew(Dart3 AD1, Dart3 AD2)
	{
		static_assert(TDim >= 0 && TDim < 4, "The parameter must be included in [0,3]");
		m_gmap.sew<TDim>(AD1, AD2);
	}
	/**@brief intiialize the block structure from a gmds cellular mesh. The provided
	 * 		 mesh @p ACellMesh must have the following characteristics:
	 * 		 - DIM3, N, E, F, R, R2N, F2N, E2N
	 *
	 * @param[in,out] ACellMesh A cellular mesh
	 */
	void init_from_mesh(Mesh &ACellMesh);

	/**@brief Convert the block structure into a gmds cellular mesh. The provided
	 * 		 mesh @p ACellMesh must have the following characteristics:
	 * 		 - DIM3, N, E, F, R, R2N, F2N, E2N
	 * 		 Attributes are exported into @p ACellMesh as variable
	 *
	 * @param[in,out] ACellMesh A cellular mesh
	 */
	void convert_to_mesh(Mesh &ACellMesh);
	/**@brief Provides a list of information about the blocking structure
	 * @return a string containing the expected pieces of information
	 */
	std::string info() const;

	/**@brief Check the topological validity of the block structur
	 *
	 * @return true if valid, false otherwise
	 */
	bool is_valid_topology() const;
	/**\brief Order the edges of @p AEdges accordingly to their distance to @p AP.
	 * 		 More specifically, we orthogonnaly project @p AP on each edge of @p AEdges.
	 * 		 For each edge, we store the distance between @p AP and the edge and the coordinate
	 * 		 the projected point inside the edge (between 0 and 1).
	 * @param[in] AP 		A Point to project on each edge
	 * @param[in] AEdges 	A set of edges
	 * @return for each edge, we get the distance (first) and the coordinate (second)
	 */
	std::vector<std::pair<double, double>> 	get_projection_info(math::Point &AP, std::vector<CurvedBlocking::Edge> &AEdges);

 private:
	/**@brief Create a node attribute in the n-gmap
	 *
	 * @param[in] AGeomDim dimension of the associated geometric cell
	 * @param[in] AGeomId id of the associated geometric cell
	 * @param[in] APoint  spatial location of the node
	 * @return the created node attribute
	 */
	Node create_node(const int AGeomDim, const int AGeomId, const math::Point &APoint);
	/**@brief Create an edge attribute in the n-gmap
	 *
	 * @param AGeomDim dimension of the associated geometric cell
	 * @param AGeomId id of the associated geometric cell
	 * @return the created edge attribute
	 */
	Edge create_edge(const int AGeomDim, const int AGeomId);
	/**@brief Create a face attribute in the n-gmap
	 *
	 * @param AGeomDim dimension of the associated geometric cell
	 * @param AGeomId id of the associated geometric cell
	 * @return the created face attribute
	 */
	Face create_face(const int AGeomDim, const int AGeomId);
	/**@brief Create a block attribute in the n-gmap
	 *
	 * @param AGeomDim dimension of the associated geometric cell
	 * @param AGeomId id of the associated geometric cell
	 * @return the created block attribute
	 */
	Block create_block(const int AGeomDim, const int AGeomId);

 private:
	/*** the associated geometric model*/
	cad::GeomManager *m_geom_model;
	/*** the underlying n-g-map model*/
	GMap3 m_gmap;
};
/*----------------------------------------------------------------------------*/
}     // namespace blocking
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_BLOCKING_H
