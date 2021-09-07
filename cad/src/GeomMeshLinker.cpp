/*----------------------------------------------------------------------------*/
//
// Created by totoro on 2019-04-13.
//
/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomMeshLinker.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
int GeomMeshLinker::m_global_link_id=0;
/*----------------------------------------------------------------------------*/
GeomMeshLinker::GeomMeshLinker()
        :m_link_id(m_global_link_id++),m_mesh(NULL),m_geometry(NULL)
{}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::GeomMeshLinker(Mesh* AMesh, GeomManager* AGeometry)
:m_link_id(m_global_link_id++)
{
setMesh(AMesh);
setGeometry(AGeometry);
}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::~GeomMeshLinker(){

    clear();
}

/*----------------------------------------------------------------------------*/
void GeomMeshLinker::clear(){
    if(m_mesh!=NULL) {
        m_mesh->deleteVariable(GMDS_NODE, m_node_classification_dim);
        m_mesh->deleteVariable(GMDS_NODE, m_node_classification_id);
        m_mesh->deleteVariable(GMDS_EDGE, m_node_classification_dim);
        m_mesh->deleteVariable(GMDS_EDGE, m_edge_classification_id);
        m_mesh->deleteVariable(GMDS_FACE, m_face_classification_dim);
        m_mesh->deleteVariable(GMDS_FACE, m_face_classification_id);
    }
    m_mesh = NULL;
    m_geometry=NULL;
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::setMesh(Mesh* AMesh){
    if(m_mesh!=NULL){
        clear();
    }

    m_mesh = AMesh;

    m_node_classification_dim = m_mesh->newVariable<ELink, GMDS_NODE>("geom_link_dim_"+std::to_string(m_link_id));
    m_node_classification_id  = m_mesh->newVariable<int, GMDS_NODE>("geom_link_id_"+std::to_string(m_link_id));

    m_edge_classification_dim = m_mesh->newVariable<ELink, GMDS_EDGE>("geom_link_dim_"+std::to_string(m_link_id));
    m_edge_classification_id  = m_mesh->newVariable<int, GMDS_EDGE>("geom_link_id_"+std::to_string(m_link_id));

    m_face_classification_dim = m_mesh->newVariable<ELink, GMDS_FACE>("geom_link_dim_"+std::to_string(m_link_id));
    m_face_classification_id  = m_mesh->newVariable<int, GMDS_FACE>("geom_link_id_"+std::to_string(m_link_id));

}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::setGeometry(GeomManager* AGeometry){
    m_geometry=AGeometry;
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkToPoint(const Node& AN, const int AGeomId){
    linkNodeToPoint(AN.id(),AGeomId);
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkNodeToPoint(const TCellID & AN, const int AGeomId){
    (*m_node_classification_dim)[AN]=LINK_POINT;
    (*m_node_classification_id)[AN]=AGeomId;

}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkToCurve(const Edge& AE, const int AGeomId){
    linkEdgeToCurve(AE.id(),AGeomId);
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkEdgeToCurve(const TCellID & AE, const int AGeomId){
    (*m_edge_classification_dim)[AE]=LINK_CURVE;
    (*m_edge_classification_id )[AE]=AGeomId;

}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkToCurve(const Node& AN, const int AGeomId){
    linkNodeToCurve(AN.id(),AGeomId);
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkNodeToCurve(const TCellID & AN, const int AGeomId){
    (*m_node_classification_dim)[AN]=LINK_CURVE;
    (*m_node_classification_id)[AN]=AGeomId;

}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkToSurface(const Node& AN, const int AGeomId){
    linkNodeToSurface(AN.id(),AGeomId);
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkNodeToSurface(const TCellID & AN, const int AGeomId){
    (*m_node_classification_dim)[AN]=LINK_SURFACE;
    (*m_node_classification_id)[AN]=AGeomId;

}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkToSurface(const Edge& AE, const int AGeomId){
    linkEdgeToSurface(AE.id(),AGeomId);
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkEdgeToSurface(const TCellID & AE, const int AGeomId){
    (*m_edge_classification_dim)[AE]=LINK_SURFACE;
    (*m_edge_classification_id )[AE]=AGeomId;

}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkToSurface(const Face& AF, const int AGeomId){
    linkFaceToSurface(AF.id(),AGeomId);
}
/*----------------------------------------------------------------------------*/
void GeomMeshLinker::linkFaceToSurface(const TCellID &AF, const int AGeomId) {
    (*m_face_classification_dim)[AF]=LINK_SURFACE;
    (*m_face_classification_id )[AF]=AGeomId;

}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::ELink GeomMeshLinker::getGeomDim(const Node& AN){return getGeomDim<Node>(AN.id());}
/*----------------------------------------------------------------------------*/
int GeomMeshLinker::getGeomId(const Node& AN){return getGeomId<Node>(AN.id());}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::ELink GeomMeshLinker::getGeomDim(const Edge& AE){return getGeomDim<Edge>(AE.id());}
/*----------------------------------------------------------------------------*/
int GeomMeshLinker::getGeomId(const Edge& AE){return getGeomId<Edge>(AE.id());}
/*----------------------------------------------------------------------------*/
GeomMeshLinker::ELink GeomMeshLinker::getGeomDim(const Face& AF){return getGeomDim<Face>(AF.id());}
/*----------------------------------------------------------------------------*/
int GeomMeshLinker::getGeomId(const Face& AF){return getGeomId<Face>(AF.id());}
/*----------------------------------------------------------------------------*/
std::pair<GeomMeshLinker::ELink,int>  GeomMeshLinker::getGeomInfo(const Node& AN){
    return getGeomInfo<Node>(AN.id());
}
/*----------------------------------------------------------------------------*/
std::pair<GeomMeshLinker::ELink,int>  GeomMeshLinker::getGeomInfo(const Edge& AN){
    return getGeomInfo<Edge>(AN.id());
}
/*----------------------------------------------------------------------------*/
std::pair<GeomMeshLinker::ELink,int>  GeomMeshLinker::getGeomInfo(const Face& AN){
    return getGeomInfo<Face>(AN.id());
}
/*----------------------------------------------------------------------------*/
template <> GMDSCad_API GeomMeshLinker::ELink GeomMeshLinker::getGeomDim<Node>(const TCellID &AN){
    return (*m_node_classification_dim)[AN];
}
/*----------------------------------------------------------------------------*/
template <> GMDSCad_API int GeomMeshLinker::getGeomId<Node>(const TCellID &AN){
    return (*m_node_classification_id)[AN];
}
/*----------------------------------------------------------------------------*/
template <> GMDSCad_API GeomMeshLinker::ELink GeomMeshLinker::getGeomDim<Edge>(const TCellID &AN){
    return (*m_edge_classification_dim)[AN];
}
/*----------------------------------------------------------------------------*/
template <>  GMDSCad_API int GeomMeshLinker::getGeomId<Edge>(const TCellID &AN){
    return (*m_edge_classification_id)[AN];
}
/*----------------------------------------------------------------------------*/
template <> GMDSCad_API GeomMeshLinker::ELink GeomMeshLinker::getGeomDim<Face>(const TCellID &AN){
    return (*m_face_classification_dim)[AN];
}
/*----------------------------------------------------------------------------*/
template <> int GeomMeshLinker::getGeomId<Face>(const TCellID &AN){
    return (*m_face_classification_id)[AN];
}
/*----------------------------------------------------------------------------*/
template <> GMDSCad_API
std::pair<GeomMeshLinker::ELink,int> GeomMeshLinker::getGeomInfo<Node>(const TCellID &AN){
    return std::make_pair((*m_node_classification_dim)[AN],
                          (*m_node_classification_id)[AN]);
}
/*----------------------------------------------------------------------------*/
template <> GMDSCad_API
std::pair<GeomMeshLinker::ELink,int> GeomMeshLinker::getGeomInfo<Edge>(const TCellID &AE){
    return std::make_pair((*m_edge_classification_dim)[AE],
                          (*m_edge_classification_id )[AE]);
}
/*----------------------------------------------------------------------------*/
template <> GMDSCad_API
std::pair<GeomMeshLinker::ELink,int> GeomMeshLinker::getGeomInfo<Face>(const TCellID &AF){
    return std::make_pair((*m_edge_classification_dim)[AF],
                          (*m_edge_classification_id )[AF]);
}
/*----------------------------------------------------------------------------*/
