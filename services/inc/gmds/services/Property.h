/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/9/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_PROPERTY_H
#define GMDS_PROPERTY_H

#include "GMDSServices_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{

    class AbstractData;
    /*------------------------------------------------------------------------*/
    /** \class Property
     *
     * \brief Interface that must be extended by any property
     *
     */class GMDSServices_API Property {
    public:

        /*------------------------------------------------------------------------*/
        /** \brief  Check if the property is verified for data @AData
         *
         * @param AData the data we check the property on
         * @return true if the property is okay, false otherwise
         */
        virtual bool isValid(const AbstractData* AData) const =0;

        virtual ~Property(){}
    protected:
         /** @brief Default constructor
          */
        Property();
    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_PROPERTY_H
/*----------------------------------------------------------------------------*/
