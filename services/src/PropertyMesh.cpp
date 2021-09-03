/*----------------------------------------------------------------------------*/
// Created by ledouxf on 2/6/19.
/*----------------------------------------------------------------------------*/
#include <gmds/services/PropertyMesh.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
PropertyMeshBuilder::~PropertyMeshBuilder() {
    for(auto p:m_builded_props)
        delete p;
}
/*----------------------------------------------------------------------------*/
Property* PropertyMeshBuilder::
build(const gmds::PropertyMeshBuilder::type AType) {
    Property* prop = 0;
    switch (AType) {
        case IS_EMPTY:
            prop = new PropertyMeshEmpty(); break;
        case IS_FULL_HEX:
            prop = new PropertyMeshFullHex(); break;
        case IS_FULL_QUAD:
            prop = new PropertyMeshFullQuad(); break;
        case IS_FULL_TET:
            prop = new PropertyMeshFullTet(); break;
        case IS_FULL_TRI:
            prop = new PropertyMeshFullTri(); break;
        default:
            throw GMDSException("Unknown mesh property type");
    };

    m_builded_props.push_back(prop);
    return prop;
}
/*----------------------------------------------------------------------------*/
bool PropertyMeshEmpty::isValid(const AbstractData* AData) const {
    const  DataMesh* data = dynamic_cast<const DataMesh*>(AData);
    if(data==0)
        throw GMDSException("Incompatible data and property");

    Mesh* mesh = data->mesh();

    return (mesh->getNbNodes()==0 &&
            mesh->getNbEdges()==0 &&
            mesh->getNbFaces()==0 &&
            mesh->getNbRegions()==0 );

}
/*----------------------------------------------------------------------------*/
bool PropertyMeshFullHex::isValid(const AbstractData* AData) const {
    const  DataMesh* data = dynamic_cast<const DataMesh*>(AData);
    if(data==0)
        throw GMDSException("Incompatible data and property");
    return (data->mesh()->getNbHexahedra()==data->mesh()->getNbRegions());
}
/*----------------------------------------------------------------------------*/
bool PropertyMeshFullTet::isValid(const AbstractData* AData) const {
    const  DataMesh* data = dynamic_cast<const DataMesh*>(AData);
    if(data==0)
        throw GMDSException("Incompatible data and property");
    return (data->mesh()->getNbTetrahedra()==data->mesh()->getNbRegions());


}
/*----------------------------------------------------------------------------*/
bool PropertyMeshFullQuad::isValid(const AbstractData* AData) const {
    const  DataMesh* data = dynamic_cast<const DataMesh*>(AData);
    if(data==0)
        throw GMDSException("Incompatible data and property");

    return (data->mesh()->getNbQuadrilaterals()==data->mesh()->getNbFaces());

}
/*----------------------------------------------------------------------------*/
bool PropertyMeshFullTri::isValid(const AbstractData* AData) const {
    const  DataMesh* data = dynamic_cast<const DataMesh*>(AData);
    if(data==0)
        throw GMDSException("Incompatible data and property");

    return (data->mesh()->getNbTriangles()==data->mesh()->getNbFaces());

}
/*----------------------------------------------------------------------------*/
