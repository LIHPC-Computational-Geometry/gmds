/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/9/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_DATA_SCALAR_H
#define GMDS_DATA_SCALAR_H
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#include <gmds/services/AbstractData.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

    /*------------------------------------------------------------------------*/
    /**@class DataScalar
     * @brief Defines a service data template for scalar type
     */
    template<typename T> class DataScalar: public AbstractData {
    public:

        /** @brief Constructor
         * @param AI    Scalar value
         * @param AName Data name
         */
        DataScalar(const T AI, const std::string AName = "int data")
                :AbstractData(AName), m_value(AI){}

        /**@brief Value getter
         *
         * @return data value
         */
        T value() const { return m_value;}

        /**@brief Value setter
         *
         * @param AV new value
         */
        void set(const T AV){m_value=AV;}
    private:
        /** Stored value **/
        T m_value;
    };

    using DataInt       = DataScalar<int>;
    using DataDouble    = DataScalar<double>;
    using DataFloat     = DataScalar<float>;
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_DATA_SCALAR_H
/*----------------------------------------------------------------------------*/
