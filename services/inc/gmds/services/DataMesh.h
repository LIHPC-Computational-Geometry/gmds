/*----------------------------------------------------------------------------*/
// Created by ledouxf on 2/6/19.
/*----------------------------------------------------------------------------*/
#ifndef GMDS_DATAMESH_H
#define GMDS_DATAMESH_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <gmds/services/AbstractData.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

    /*------------------------------------------------------------------------*/
    /**@class DataMesh
     * @brief Defines a service data template for a mesh object
     */
    class DataMesh: public AbstractData {
    public:

        /** @brief Constructor
         * @param AI    Scalar value
         * @param AName Data name
         */
        DataMesh(Mesh* AMesh, const std::string AName = "Mesh data");

        /**@brief Value getter
         *
         * @return data value
         */
        Mesh* mesh() const;

        /**@brief Value setter
         *
         * @param AV new value
         */
        void set( Mesh* AM);
    private:
        /** Stored value **/
        Mesh* m_mesh;
    };

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_DATAMESH_H
/*----------------------------------------------------------------------------*/

