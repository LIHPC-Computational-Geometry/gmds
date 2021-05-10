/*----------------------------------------------------------------------------*/
#ifndef GMDS_SHEET_OPERATOR_2D_H
#define GMDS_SHEET_OPERATOR_2D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  Operator2D
     *  @brief  Abstract class that gathers methods and behaviors that are common
     *          to all the sheet operators in 2D (select, pillow, collapse)
     *
     *          The mesh we work on  must have F and N and the F2N fields and not R
     */
    class EXPORT_GMDS Operator2D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Check if the mesh fits algorithm requirements, which are:
         *         - to be a full quad mesh
         */
        bool isValid();

    protected:

        /*------------------------------------------------------------------------*/
        /** @brief  Constructor.
         * @param AMesh the mesh to apply the operator on
         */
        Operator2D(Mesh* AMesh);
        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.
         */
        virtual ~Operator2D();


        /*------------------------------------------------------------------------*/
        /** @brief  Check if there exist an edge connecting AN1 and AN2
         * @param[in] AN1 id of the first edge extremity
         * @param[in] AN2 id of the second edge extremity
         */
        bool isAnEdge(const TCellID AN1, const TCellID AN2);

        /*------------------------------------------------------------------------*/
        /** @brief Get the faces sharing edge \p AE
         * @param[in] AE a virtual edge
         * @return the adjacent faces (1 or 2)
         */
        std::vector<TCellID> getAdjacentFaces(const VirtualEdge& AE);
        /*------------------------------------------------------------------------*/
        /** @brief  get the set of quads that are both adjacent to AN1 and AN2
         * @param[in] AN1 A first node id
         * @param[in] AN2 A second node id
         */
        std::vector<TCellID> getAdjacentFaces(const TCellID AN1,const TCellID AN2);



    protected:
        /** a mesh */
        Mesh* m_mesh;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_SHEET_COLLAPSE_H
