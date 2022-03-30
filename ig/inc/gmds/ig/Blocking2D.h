/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_H
#define GMDS_BLOCKING_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIg_export.h"
#include <map>
#include <gmds/math/DiscretizationScheme1D.h>
#include <gmds/utils/Array.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*----------------------------------------------------------------------------*/
    /** @class Blocking2D
     *  @brief The blocking2D structure is a full unstructured mesh made of hex
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
    class GMDSIg_API Blocking2D: public Mesh{
    public:

        class GMDSIg_API Block{
        public:
            friend class Blocking2D;
            /**
             *
             * @return the id of the origin corner node, which the first face node
             */
            TCellID origin();

            Node getNode(const int& AIndex);
            Node operator()(const int AI, const int AJ);
            Edge getEdgeI();
            Edge getEdgeJ();
            math::Vector3d getUnitVectorI();
            math::Vector3d getUnitVectorJ();

            TCellID id() {return m_face.id();}
            /**@brief Only uniform discretization is supported right now. It will be possible
             * to extend it easily in a short future. Default is 10.
             */
            void setNbDiscretizationI(const int AN);
            void setNbDiscretizationJ(const int AN);
            int getNbDiscretizationI() const;
            int getNbDiscretizationJ() const;



         private:
            /** Access to the edge with local index @p AI nad @p AJ in the
             *  current face
             *
             *  An exception is throww if @p AI==@p AJ and if @p AI or @p AJ is not included in [0,3]
             *
             * @param AI local index of a node in the face
             * @param AJ local index of a node in the face
             * @return the edge connecting the two input nodes
             */
            Edge getEdge(const int AI, const int AJ);
        private:
            /** @brief private constructor, only a Blocking2D instance can
             *         construct such an object.
             */
            Block(const Face&   AFace  , Blocking2D* ASupport);
            Block(const TCellID AFaceID, Blocking2D* ASupport);

            /** Face at the block level. We keep it instead of the id because
             * it is more rich and in particular stores the reference to the
             * mesh it belongs to.
             */
            Face m_face;
            Blocking2D* m_support;
            Array2D<TCellID>* m_grid_view;

        };
        /** @brief Constructor
         */
        Blocking2D();

	     /** @brief Constructor that take a mesh in input to create the blocks and set the I/J discretization value for the blocks to @p AN.
	      *
	      * @param AMesh a mesh with MeshModel (DIM3|F|E|N|F2N|E2N|F2E|E2F|N2E|N2F)
	      * @param AN value of the I/J discretization, if not informed set by default to 10
	      */
	     Blocking2D(const Mesh& AMesh, int AN = 10);
        /** @brief Destructor
         */
        virtual ~Blocking2D();

        Node newBlockCorner(const math::Point& APnt);
        Node newBlockCorner(const TCoord AX,const TCoord AY, const TCoord AZ=0);
        /** Block creation from 4 nodes. The face is oriented according to the nodes.
         *  the origin node is the first one and the I direction is given by @p AN1 -> @p AN2,
         *  and the J direction by @p AN1 -> @p AN3. Once created, it can not be changed.
         *
         *  Block edges are created on the fly if it wasn't yet available in the mesh
         *
         * @param AN1 a first node
         * @param AN2 a second node
         * @param AN3 a third node
         * @param AN4 a fourth node
         * @return A 2D block, which is a quad face with extra attributes
         */
        Block newBlock(const Node& AN1, const Node& AN2, const Node& AN3, const Node& AN4);
        Block newBlock(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);

        /** Return a handle block for block of id @p AId
         * @param AId block id
         * @return the expected block
         */
        Block block(const TCellID AId);

        /** Returns all the blocks
         * @return a collection of block handlers
         */
        std::vector<Block> allBlocks();
        /** create all the edge and block inner nodes. Once this method call the block structure
         * can not be modified anymore. We should check the conditions to allow it in the future.
         */
        void initializeGridPoints();

	     /** Get the neighboring nodes id of node of id @p AId in I+/-1, J+/-1 if they exist. Does not support singular node.
	      *
	      * @param AId node id
	      * @return a collection of nodes ids
	      */
	     std::vector<TCellID> getNodeNeighbors(TCellID AId);

	     /** Build the blocks from the cells of mesh's object then discretize the created blocks in I/J with the value of @p AN.
	      *
	      * @param AN the value of the discretization for the created blocks in I and J
	      */
	     void buildBlocks(const int AN);

	     /** Get the id of the blocking entity (node, edge, face) on which the grid node @p AId is classified.
	      *
	      * @param AId id of the grid node
	      * @return Id of the blocking entity
	      */
	     TCellID getBlockingId(TCellID AId);

	     /** Get the dimension of the blocking entity (node, edge, face) on which the grid node @p AId is classified
	      *
	      * @param AId id of the grid node
	      * @return Dimension of the blocking entity
	      */
	     int getBlockingDim(TCellID AId);

    private:

        /**
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

        /** Give the number of subdivision for edge @p AEdge in face @p AFace
          *
          * @param AFace the face we look into
          * @param AEdge the edge that must be adjacent to @p AFace
          * @return the discretisation of @p AEdge
          */
        int getNbDiscretization(const Face& AFace, const Edge& AEdge) const;

    private:
        /** the embedding variable, which allows us to know where each mesh node is*/
        Variable<int>* m_embedding_dim;
        Variable<TCellID>* m_embedding_id;
        /** variable associated to faces*/
        Variable<int>* m_discretization_I;
        Variable<int>* m_discretization_J;
	     /** variables that store the grids nodes for blocks entities*/
        Variable<Array2D<TCellID>* >* m_face_grids;
        Variable<std::vector<TCellID>* >* m_edge_grids;

    };
}
#endif //GMDS_BLOCKING_H
