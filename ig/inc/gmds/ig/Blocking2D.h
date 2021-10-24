/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKING_H
#define GMDS_BLOCKING_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIg_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*----------------------------------------------------------------------------*/
    /** @class Blocking2D
     *  @brief The blocking2D structure is a full unstructured mesh made of hex
     *         blocks, quad faces, edges and nodes. Block, faces and edges only
     *         exist at the macro level, or block level. Nodes are block corners
     *         but also inner-block nodes.
     *
     *         Each node knows it is a corner-block (block level) or an inner-node
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
    class GMDSIg_API Blocking2D{
            public:

            class Block2D{
                public:
                friend class Blocking2D;
                private:
                /** @brief private constructor, only a Blocking2D instance can
                 *         construct such an object.
                 */
                Block2D();
                /** Face at the block level. We keep it instead of the id because
                 * it is more rich and in particular stores the mesh pointer to.
                 */
                Face m_block;

            };
            /*------------------------------------------------------------------------*/
            /** @brief Constructor
             */
            Blocking2D();

            Node newBlockCorner(const math::Point& APnt);
            Node newBlockCorner(const TCoord AX,const TCoord AY, const TCoord AZ=0);
            Region newBlock(const Node& AN1, const Node& AN2, const Node& AN3, const Node& AN4);

            private:
            Mesh m_blocks;
            Variable<Cell::Data>* m_embedding;
    }
}
#endif //GMDS_BLOCKING_H
