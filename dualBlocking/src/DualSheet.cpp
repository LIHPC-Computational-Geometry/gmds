//
// Created by calderans on 16/12/19.
//


#include "gmds/dualBlocking/DualSheet.h"


using namespace gmds;
using namespace db;

DualSheet::DualSheet(std::vector<TCellID> ASurface,int AID, Mesh* AMesh):m_mesh_surface(AMesh){
  m_surface = ASurface;
  m_sheet_ID = AID;
}
/*------------------------------------------------------------------------*/
DualSheet::DualSheet(int AID, Mesh* AMesh):m_mesh_surface(AMesh){
    m_sheet_ID = AID;
}
/*------------------------------------------------------------------------*/
DualSheet::DualSheet(const DualSheet &ASheet){
  this->m_surface = ASheet.m_surface;
  this->m_sheet_ID = ASheet.m_sheet_ID;
  this->boundary = ASheet.boundary;
  this->m_mesh_surface = ASheet.m_mesh_surface;
}
/*------------------------------------------------------------------------*/
void DualSheet::setSurface(std::vector<gmds::TCellID> ASurface){
  for(auto id : ASurface){
    m_surface.push_back(id);
  }
}
/*------------------------------------------------------------------------*/
DualSheet::~DualSheet(){

}
/*------------------------------------------------------------------------*/
void DualSheet::setSurface(std::map<gmds::TCellID,gmds::math::Vector3d> ASurface){
  for(auto id : ASurface){
    m_surface.push_back(id.first);
  }
}
/*------------------------------------------------------------------------*/
std::vector<gmds::TCellID> DualSheet::getSurface(){
  return m_surface;
}
/*------------------------------------------------------------------------*/
int DualSheet::getID(){
  return m_sheet_ID;
}
/*------------------------------------------------------------------------*/
void DualSheet::boundarySheet() {
    boundary = true;
}
/*------------------------------------------------------------------------*/
bool DualSheet::isBoundary() {
    return boundary;
}
/*------------------------------------------------------------------------*/
void DualSheet::buildSurface(gmds::Mesh* AMesh) {

    std::cout<<"BUILD Nb points on surf "<<AMesh->getNbNodes()<<std::endl;
    std::cout<<"BUILD Nb faces on surf "<<AMesh->getNbFaces()<<std::endl;

    std::map<TCellID,TCellID> amesh_to_surface;

    for(auto n : AMesh->nodes()){
        Node node = AMesh->get<Node>(n);
        Node surface_node = m_mesh_surface->newNode(node.getPoint());
        amesh_to_surface.emplace(n,surface_node.id());
    }


    for(auto f : AMesh->faces()){
        Face face = AMesh->get<Face>(f);
        std::vector<TCellID > nodes;
        for(auto n : face.getIDs<Node>()){
            nodes.push_back(amesh_to_surface[n]);
        }
        m_mesh_surface->newFace(nodes);
    }
    std::cout<<"surface Nb points on surf "<<m_mesh_surface->getNbNodes()<<std::endl;
    std::cout<<"surface Nb faces on surf "<<m_mesh_surface->getNbFaces()<<std::endl;
}
/*------------------------------------------------------------------------*/
Mesh* DualSheet::getSurfaceMesh(){
    return m_mesh_surface;
}
