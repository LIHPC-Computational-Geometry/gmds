/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/9/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_PROPERTYSCALAR_H
#define GMDS_PROPERTYSCALAR_H
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/services/Property.h>
#include <gmds/services/DataScalar.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

    //forward declaration
    template<typename TScalar> class PropertyScalarInRange;
    /*------------------------------------------------------------------------*/
    /**@class PropertyScalarBuilder
     * @brief Builder of scalar properties. Each builder stores and has the
     *        responsability of destroying properties
     */
    class PropertyScalarBuilder {
    public:
        typedef enum {
            POSITIVE,
            NEGATIVE,
            POSITIVE_STRICTLY,
            NEGATIVE_STRICTLY
        } type;


        virtual ~PropertyScalarBuilder(){
            for(auto p:m_builded_props)
                delete p;
        }

        template<typename TScalar> Property* build(const type AType);
        template<typename TScalar> PropertyScalarInRange<TScalar>* buildInRange();


    private:

        std::vector<Property*> m_builded_props;

    };


    /*------------------------------------------------------------------------*/
    /**@class PropertyScalarPositive
     * @brief Properties for a scalar data to be positive
     */
    template<typename T>
    class PropertyScalarPositive : public Property {
        friend class PropertyScalarBuilder;
    private:
        PropertyScalarPositive() : Property(){;}
    public:
        bool isValid(const AbstractData* AData) const {
            const  DataScalar<T>* scalar = dynamic_cast<const DataScalar<T>*>(AData);
            if(scalar==0)
                throw GMDSException("Incompatible data and property");

            return scalar->value()>=0;
        }
    };
    /*------------------------------------------------------------------------*/
    /**@class PropertyScalarNegative
     * @brief Properties for a data to be negative
     */
    template<typename T>
    class PropertyScalarNegative: public Property {
        friend class PropertyScalarBuilder;
    private:
        PropertyScalarNegative() : Property(){;}
    public:
        bool isValid(const AbstractData* AData) const {
            const DataScalar<T>* scalar = dynamic_cast<const DataScalar<T>*>(AData);
            if(scalar==0)
                throw GMDSException("Incompatible data and property");

            return scalar->value()<=0;
        }

    };
    /*------------------------------------------------------------------------*/
    /**@class PropertyScalarPositiveStrictly
     * @brief Properties for a data to be strictly positive
     */
    template<typename T>
    class PropertyScalarPositiveStrictly : public Property {
        friend class PropertyScalarBuilder;
    private:
        PropertyScalarPositiveStrictly() : Property(){;}
    public:
        bool isValid(const AbstractData* AData) const {
            const DataScalar<T>* scalar = dynamic_cast<const DataScalar<T>*>(AData);
            if(scalar==0)
                throw GMDSException("Incompatible data and property");

            return scalar->value()>0;
        }
    };

    /*------------------------------------------------------------------------*/
    /**@class PropertyScalarNegativeStrictly
     * @brief Properties for a data to be strictly negative
     */
    template<typename T>
    class PropertyScalarNegativeStrictly : public Property {
        friend class PropertyScalarBuilder;
    private:
        PropertyScalarNegativeStrictly() : Property(){;}
    public:
        bool isValid(const AbstractData* AData) const {
            const  DataScalar<T>* scalar = dynamic_cast<const DataScalar<T>*>(AData);
            if(scalar==0)
                throw GMDSException("Incompatible data and property");

            return scalar->value()<0;
        }
    };

    /*------------------------------------------------------------------------*/
    /**@class PropertyScalarInRange
     * @brief Properties for a data to be in a specified range of values
     */
    template<typename T>
    class PropertyScalarInRange : public Property {
        friend class PropertyScalarBuilder;
    private:
        PropertyScalarInRange() : Property(), m_init(false), m_min(0), m_max(0){;}
    public:
        void setRange(const T AMin, const T AMax){
            m_init=true;
            m_min = AMin;
            m_max = AMax;
        }

        bool isValid(const AbstractData * AData) const  {
          const  DataScalar<T>* scalar = dynamic_cast<const DataScalar<T>*>(AData);

            if(scalar==0)
                throw GMDSException("Incompatible data and property");

            if(m_init==false)
                throw GMDSException("Min and Max values of a scalar InRange Property are not set!");

            return (m_min<=scalar->value()) && (scalar->value()<=m_max);
        }


    private:
        bool m_init;
        T m_min;
        T m_max;
    };


    template<typename TScalar> Property*
    PropertyScalarBuilder::build(const PropertyScalarBuilder::type AType){
        Property *prop;
        switch (AType) {
            case PropertyScalarBuilder::POSITIVE:
                prop = new  PropertyScalarPositive<TScalar>(); break;
            case PropertyScalarBuilder::NEGATIVE:
                prop = new PropertyScalarNegative<TScalar>(); break;
            case PropertyScalarBuilder::POSITIVE_STRICTLY:
                prop = new PropertyScalarPositiveStrictly<TScalar>(); break;
            case PropertyScalarBuilder::NEGATIVE_STRICTLY:
                prop = new PropertyScalarNegativeStrictly<TScalar>(); break;
            default:
                throw GMDSException("Unknown mesh property type");
        };

        m_builded_props.push_back(prop);
        return prop;
    }

    template<typename TScalar> PropertyScalarInRange<TScalar>*
    PropertyScalarBuilder::buildInRange(){
        PropertyScalarInRange<TScalar>* prop = new PropertyScalarInRange<TScalar>();
        m_builded_props.push_back(prop);
        return prop;

    }

}


#endif //GMDS_PROPERTYSCALAR_H
