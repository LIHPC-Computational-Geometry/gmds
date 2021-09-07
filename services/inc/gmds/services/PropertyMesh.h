/*----------------------------------------------------------------------------*///
// Created by ledouxf on 2/6/19.
/*----------------------------------------------------------------------------*/
#ifndef GMDS_PROPERTYMESH_H
#define GMDS_PROPERTYMESH_H
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/services/Property.h>
#include <gmds/services/DataMesh.h>

/*----------------------------------------------------------------------------*/
namespace gmds {

    /*------------------------------------------------------------------------*/
    /**@class PropertyMeshBuilder
     * @brief Builder of mesh properties. Each builder stora and has the
     *        responsability of
     */
    class PropertyMeshBuilder {
    public:
        typedef enum {
            IS_EMPTY,
            IS_FULL_HEX,
            IS_FULL_TET,
            IS_FULL_QUAD,
            IS_FULL_TRI,
        } type;

        virtual ~PropertyMeshBuilder();

        Property* build(const type AType);

    private:

        std::vector<Property*> m_builded_props;

    };

    /*------------------------------------------------------------------------*/
    class PropertyMeshEmpty : public Property{
        friend class PropertyMeshBuilder;
    private:
        PropertyMeshEmpty(){;}
    public:
        bool isValid(const AbstractData* AData) const;
    };
    /*------------------------------------------------------------------------*/
    class PropertyMeshFullHex: public Property{
        friend class PropertyMeshBuilder;
    private:
        PropertyMeshFullHex(){;}
    public:
        bool isValid(const AbstractData* AData) const;
    };
    /*------------------------------------------------------------------------*/
    class PropertyMeshFullQuad: public Property{
        friend class PropertyMeshBuilder;
    private:
        PropertyMeshFullQuad(){;}
    public:
        bool isValid(const AbstractData* AData) const;
    };
    /*------------------------------------------------------------------------*/
    class PropertyMeshFullTet: public Property{
        friend class PropertyMeshBuilder;
    private:
        PropertyMeshFullTet(){;}
    public:
        bool isValid(const AbstractData* AData) const;
    };
    /*------------------------------------------------------------------------*/
    class PropertyMeshFullTri: public Property{
        friend class PropertyMeshBuilder;
    private:
        PropertyMeshFullTri(){;}
    public:
        bool isValid(const AbstractData* AData) const;
    };

}
#endif //GMDS_PROPERTYMESH_H
