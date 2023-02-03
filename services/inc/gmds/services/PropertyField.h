/*----------------------------------------------------------------------------*/
// Created by ledouxf on 2/6/19.
/*----------------------------------------------------------------------------*/
#ifndef GMDS_PROPERTYMESHVARIABLE_H
#define GMDS_PROPERTYMESHVARIABLE_H
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/services/Property.h>
#include <gmds/services/DataField.h>
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

    /*------------------------------------------------------------------------*/
    /**@class PropertyFieldBuilder
     * @brief Builder of mesh variable properties. Each builder stores and has the
     *        responsability of destroying properties
     */
    class PropertyFieldBuilder {
    public:
        typedef enum {
            ON_NODE,
        } type;


        virtual ~PropertyFieldBuilder(){
            for(auto p:m_builded_props)
                delete p;
        }

        template<typename T> Property* build(Mesh* AMesh, const type AType);

    private:

        std::vector<Property*> m_builded_props;

    };

    /*------------------------------------------------------------------------*/
    /**@class PropertyFieldOnMesh
     * @brief Property that a field lived on a specified mesh
     */
    class PropertyFieldOnMesh: public Property {
    protected:
        PropertyFieldOnMesh(Mesh* AMesh) : Property(), m_mesh(AMesh){}

    public:

        virtual bool isValidOnMesh(const AbstractData* AData) const=0;
        bool isValid(const AbstractData* AData) const {
            const DataFieldOnMesh *field = dynamic_cast<const DataFieldOnMesh*>(AData);
            if (field == 0)
                throw GMDSException("Incompatible data and property");

            return (field->mesh()== m_mesh) && isValidOnMesh(AData);
        }
    private:
        Mesh* m_mesh;
    };

    /*------------------------------------------------------------------------*/
    /**@class PropertyOnNode
     * @brief Property attached to a node
     */
    template<typename T>
    class PropertyFieldOnNode: public PropertyFieldOnMesh {
        friend class PropertyFieldBuilder;
    private:
        PropertyFieldOnNode(Mesh* AMesh) : PropertyFieldOnMesh(AMesh){}

    public:
        bool isValidOnMesh(const AbstractData* AData) const {
            const  DataField<T>* field = dynamic_cast<const DataField<T>*>(AData);
            if(field==0)
                throw GMDSException("Incompatible data and property");

            Mesh* m = field->mesh();

            try{
                m->getVariable<T, GMDS_NODE>(field->name());
            }
            catch(GMDSException& e) {
                return false;
            }
            return true;
        }
    };


    template<typename T> Property*
    PropertyFieldBuilder::build(Mesh* AMesh,
                                const PropertyFieldBuilder::type AType){
        Property *prop;
        switch (AType) {
            case PropertyFieldBuilder::ON_NODE:
                prop = new  PropertyFieldOnNode<T>(AMesh); break;
            default:
                throw GMDSException("Unknown mesh property type");
        };

        m_builded_props.push_back(prop);
        return prop;
    }



}


#endif //GMDS_PROPERTYMESHVARIABLE_H
