/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/9/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ABSTRACTDATA_H
#define GMDS_ABSTRACTDATA_H
/*----------------------------------------------------------------------------*/
#include <vector>
#include <string>
/*----------------------------------------------------------------------------*/
#include "gmds/services/AbstractService.h"
/*----------------------------------------------------------------------------*/
namespace gmds {

    /*------------------------------------------------------------------------*/
    /** \class AbstractData
     *
     * \brief Any data that we provide to a service must extend this abstract
     *        class
     *
     */
    class AbstractData {


    public:
        /** @brief provides the name of the data
         * @return a name
         */
        std::string name() const;

    protected:

        /** Constructor method
         *
         * @param AName name given to the data
         */
        AbstractData(const std::string AName="default_name");

        virtual ~AbstractData(){}
    private:

        /** @brief data name **/
        std::string m_name;

    };
}

/*----------------------------------------------------------------------------*/
#endif //GMDS_ABSTRACTDATA_H
/*----------------------------------------------------------------------------*/
