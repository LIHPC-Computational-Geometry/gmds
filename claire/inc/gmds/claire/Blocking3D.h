//
// Created by rochec on 29/11/23.
//

#ifndef GMDS_BLOCKING3D_H
#define GMDS_BLOCKING3D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIg_export.h"
#include <map>
#include <gmds/math/DiscretizationScheme1D.h>
#include <gmds/utils/Array.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** @class Blocking3D
 *  @brief This structure is built from Blocking2D.
 *
 *  		  The blocking3D structure is a full unstructured mesh made of hex
 *         blocks, quad faces, edges and nodes. Block, faces and edges only
 *         exist at the macro level, or block level. Nodes are block corners
 *         but also inner-block nodes.
 *
 p)*         Each node knows it is a corner-block (block level) or an inner-node
 *         using the node on-node variable called "embedding". An embedding is
 *         a Cell::Data that corresponds to the block entity the node is
 *         embedded in. It is made of a dim and an id.
 *         It can be on a corner-block(0), face-block(2), edge-block (1) and
 *         inner-block(3).
 *
 *         Each block as a grid structure. It must be so possible to use bracket
 *         notation [i+,j-1] to acces to a neighbor nodes. This type of traversal
 *         should be possible almost everywhere and allows the user to traverse
 *         several blocks in a row. The only issue is when you meet singular
 *         corners (ie with valence not equal to 4).
 *
 *         Add a gridview that is linked to a blocking2D objet and starting
 *         from a node encapsulates the grid view traversal.
 */
/*----------------------------------------------------------------------------*/
class GMDSIg_API Blocking3D : public Mesh
{
 public:
	class GMDSIg_API Block{
	 public:
		friend class Blocking3D;
		/**
             *
             * @return the id of the origin corner node, which the first face node
		 */
		TCellID origin();

		Node getNode(const int& AIndex);
		Node operator()(const int AI, const int AJ, const int AK);
		Edge getEdgeI();
		Edge getEdgeJ();
		Edge getEdgeK();
		math::Vector3d getUnitVectorI();
		math::Vector3d getUnitVectorJ();
		math::Vector3d getUnitVectorK();

		TCellID id() {return m_region.id();}
		/**@brief Only uniform discretization is supported right now. It will be possible
             * to extend it easily in a short future. Default is 10.
		 */
		void setNbDiscretizationI(const int AN);
		void setNbDiscretizationJ(const int AN);
		void setNbDiscretizationK(const int AN);
		int getNbDiscretizationI() const;
		int getNbDiscretizationJ() const;
		int getNbDiscretizationK() const;

		/** @brief Get the indices I and J of node @p AID in the block.
		       *
		       * @param AID a node id
		       * @return A pair containing the indice I and J
		 */
		std::tuple<int,int,int> getIndices(const TCellID AID);

		/** @brief Use to know if the edge @p AID is one of the two edges on the I axe.
		       *
		       * @param AID an edge id
		       * @return true if the edge is on I, false otherwise
		 */
		bool isEdgeOnI(TCellID AID);
		/** @brief Use to know if the edge @p AID is one of the two edges on the J axe.
		       *
		       * @param AID an edge id
		       * @return true if the edge is on J, false otherwise
		 */
		bool isEdgeOnJ(TCellID AID);
		/** @brief Use to know if the edge @p AID is one of the two edges on the K axe.
		       *
		       * @param AID an edge id
		       * @return true if the edge is on K, false otherwise
		 */
		bool isEdgeOnK(TCellID AID);

	 private:
		/** Access to the edge with local index @p AI and @p AJ in the
             *  current face
             *
             *  An exception is throww if @p AI==@p AJ and if @p AI or @p AJ is not included in [0,3]
             *
             * @param AI local index of a node in the face
             * @param AJ local index of a node in the face
             * @return the edge connecting the two input nodes
		 */
		Edge getEdge(const int AI, const int AJ);
		/** Access to the face with local index @p AI, @p AJ, @p AK and @p AL in the
             *  current hex
             *
             *  An exception is throww if @p AI==@p AJ and if @p AI, @p AJ, @p AK or @p AL is not included in [0,7]
             *
             * @param AI local index of a node in the face
             * @param AJ local index of a node in the face
             * @param AK local index of a node in the face
             * @param AL local index of a node in the face
             * @return the edge connecting the two input nodes
		 */
		Face getFace(const int AI, const int AJ, const int AK, const int AL);
	 private:
		/** @brief private constructor, only a Blocking2D instance can
             *         construct such an object.
		 */
		Block(const Region&   ARegion, Blocking3D* ASupport);
		Block(const TCellID ARegionID, Blocking3D* ASupport);

		/** Face at the block level. We keep it instead of the id because
             * it is more rich and in particular stores the reference to the
             * mesh it belongs to.
		 */
		Region m_region;
		Blocking3D* m_support;
		Array3D<TCellID>* m_grid_view;

	};

	/** @brief Constructor
	 */
	Blocking3D();

	/** @brief Constructor that take a mesh in input to create the blocks.
	      *
	      * @param AMesh a mesh with MeshModel (DIM3|F|E|N|F2N|E2N|F2E|E2F|N2E|N2F)
	 */
	explicit Blocking3D(const Mesh &AMesh);
	/** @brief Destructor
	 */
	virtual ~Blocking3D();

	Node newBlockCorner(const math::Point& APnt);
	Node newBlockCorner(const TCoord AX,const TCoord AY, const TCoord AZ=0);

	/** Block creation from 8 nodes. The region is oriented according to the nodes.
         *  the origin node is the first one and the I direction is given by @p AN1 -> @p AN2,
         *  the J direction by @p AN1 -> @p AN3 and the K direction by @p AN1 -> @p AN5.
         *  Once created, it can not be changed.
         *
         *  Block edges are created on the fly if it wasn't yet available in the mesh
         *
         * @param AN1 a first node
         * @param AN2 a second node
         * @param AN3 a third node
         * @param AN4 a fourth node
         * @param AN5 a fifth node
         * @param AN6 a sixth node
         * @param AN7 a seventh node
         * @param AN8 a eighth node
         * @return A 3D block, which is a hex region with extra attributes
	 */
	Block newBlock(const Node& AN1, const Node& AN2, const Node& AN3, const Node& AN4,
	               const Node& AN5, const Node& AN6, const Node& AN7, const Node& AN8);
	Block newBlock(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
	               const TCellID AN5, const TCellID AN6, const TCellID AN7, const TCellID AN8);

	/** Return a handle block for block of id @p AId
         * @param AId block id
         * @return the expected block
	 */
	Block block(const TCellID AId);
	/** Returns all the blocks
         * @return a collection of block handlers
	 */
	std::vector<Block> allBlocks();
	/** Create all the edge nodes. Once this method call the block structure
         * can not be modified anymore. We should check the conditions to allow it in the future.
	 */
	void initializeEdgesPoints();
	/** Create all the inner faces nodes, requires initializeEdgesPoints. Once this method call the block structure
         * can not be modified anymore. We should check the conditions to allow it in the future.
	 */
	void initializeFacesPoints();
	/** Create all the inner blocks nodes, requires initializeFacesPoints. Once this method call the block structure
         * can not be modified anymore. We should check the conditions to allow it in the future.
	 */
	void initializeBlocksPoints();

 private:
	/** Give the number of subdivision for edge @p AEdge in region @p ARegion
          *
          * @param ARegion the region we look into
          * @param AEdge the edge that must be adjacent to @p ARegion
          * @return the discretisation of @p AEdge
	 */
	int getNbDiscretization(const Region& ARegion, const Edge& AEdge) const;
	/** Check if the discretization of an edge is set the same in each adjacent face.
         *
	 */
	bool checkDiscretizationValidity() const;
	/** Check if the edge [@p AN1, @p AN2] exists in the mesh. If not, it is created.
         * We use and update the variable m_n2e to do that.
         *
         * @param AN1 first node
         * @param AN2 second node
         * @return the id of the edge connecting @p AN1 and @p AN2
	 */
	TCellID getEdge(const TCellID AN1, const TCellID AN2);
	/** Check if the face [@p AN1, @p AN2, @p AN3, @p AN4] exists in the mesh. If not, it is created.
         * We use and update the variable m_n2e to do that.
         *
         * @param AN1 first node
         * @param AN2 second node
         * @param AN3 third node
         * @param AN4 fourth node
         * @return the id of the face connecting @p AN1, @p AN2, @p AN3 and @p AN4
	 */
	TCellID getFace(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);
 private:
	/** the embedding variable, which allows us to know where each mesh node is*/
	Variable<int> *m_embedding_dim;
	Variable<TCellID> *m_embedding_id;
	/** variable associated to regions*/
	Variable<int> *m_discretization_I;
	Variable<int> *m_discretization_J;
	Variable<int> *m_discretization_K;
	/** variables that store the grids nodes for blocks entities*/
	Variable<Array3D<TCellID> *> *m_region_grids;
	Variable<Array2D<TCellID> *> *m_face_grids;
	Variable<std::vector<TCellID> *> *m_edge_grids;
};
/*----------------------------------------------------------------------------*/
}
#endif     // GMDS_BLOCKING3D_H