/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Operator2D.h>
#include <gmds/ig/MeshDoctor.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Operator2D::Operator2D(Mesh* AMesh): m_mesh(AMesh)
{}
/*----------------------------------------------------------------------------*/
Operator2D::~Operator2D()
{}
/*----------------------------------------------------------------------------*/
bool Operator2D::isValid()
{
    for(auto f_id:m_mesh->faces()){
        Face f=m_mesh->get<Face>(f_id);
        if (f.type()!=GMDS_QUAD){
            return false;
        }
    }

    return  (!m_mesh->getModel().has(R) &&
             m_mesh->getModel().has(F)&&
             m_mesh->getModel().has(N)&&
             m_mesh->getModel().has(F2N)&&
             m_mesh->getModel().has(N2F) &&
             m_mesh->getModel().has(E2N));
}
/*----------------------------------------------------------------------------*/
bool Operator2D::isAnEdge(const TCellID AN1, const TCellID AN2) {
    //we pick the first adjacent hex of AN1
    std::vector<TCellID> f1_ids = m_mesh->get<Node>(AN1).getIDs<Face>();
    for(auto i1:f1_ids) {
        Face f = m_mesh->get<Face>(i1);
        TCellID adj_na=NullID;
        TCellID adj_nb=NullID;
        f.getAdjacentNodes(AN1,adj_na,adj_nb);
        if(adj_na==AN2 || adj_nb==AN2) {
            return true;
        }
    }
    return false;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Operator2D::getAdjacentFaces(const TCellID AN1,
                                                       const TCellID AN2)
{
    std::vector<TCellID> q1 = m_mesh->get<Node>(AN1).getIDs<Face>();
    std::vector<TCellID> q2 = m_mesh->get<Node>(AN2).getIDs<Face>();
    std::vector<TCellID> res;
    for(auto i1:q1){
        for(auto i2:q2){
            if(i1==i2){
                res.push_back(i1);
            }
        }
    }
    return res;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Operator2D::getAdjacentFaces(const VirtualEdge &AE) {
    return getAdjacentFaces(AE.first(),AE.second());
}
