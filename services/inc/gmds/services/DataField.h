/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 2/6/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_DATAFIELD_H
#define GMDS_DATAFIELD_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <gmds/services/AbstractData.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    /**@class DataFieldOnMesh
     * @brief Defines a service data template for an object attached to a mesh
     */

    class DataFieldOnMesh: public AbstractData{
    public:

        /** @brief Constructor
        * @param AI    Scalar value
        * @param AName Data name
        */
        DataFieldOnMesh(Mesh* AMesh, const std::string AName = "Mesh data")
        :AbstractData(AName), m_mesh(AMesh){}

        /**@brief Value getter
         *
         * @return data value
         */
        Mesh* mesh() const { return m_mesh;}

        /**@brief Value setter
         *
         * @param AV new value
         */
        void set( Mesh* AM) {m_mesh = AM;}

    private:
        /** Stored value **/
        Mesh* m_mesh;
    };
    /*------------------------------------------------------------------------*/
    /**@class DataField
     * @brief Defines a service data template for a field object
     */
    template<typename T>
    class DataField: public DataFieldOnMesh {
    public:

        /** @brief Constructor
         * @param AI    Scalar value
         * @param AName Data name
         */
        DataField(Mesh* AMesh, const std::string AName = "Mesh data")
        :DataFieldOnMesh(AMesh, AName){}

  };

using DataFieldInt       = DataField<int>;
using DataFieldDouble    = DataField<double>;
using DataFieldFloat     = DataField<float>;

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_DATAFIELD_H
/*----------------------------------------------------------------------------*/

