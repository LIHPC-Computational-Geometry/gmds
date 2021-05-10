/*----------------------------------------------------------------------------*/
#ifndef GMDS_THEXBUILDER_H
#define GMDS_THEXBUILDER_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  THexBuilder
     *  @brief  Class that provides a way to transform a tet mesh into a hex one
     *          by splitting each tet into 4 hexes.
     *
     *          The mesh that is transformed must only have R and N and the R2N
     *          fields
     *          Geometric classification is not handled yet.
     */
    class EXPORT_GMDS THexBuilder
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         */
        THexBuilder(Mesh* AMesh);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~THexBuilder();

        /*------------------------------------------------------------------------*/
        /** @brief Check if the mesh fits algorithm requirements
         */

        bool isValid() const;
        /*------------------------------------------------------------------------*/
        /** @brief  Performs the mesh transformation
         */

        void execute();
    protected:
        /* a mesh */
        Mesh* m_mesh;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_THEXBUILDER_H
