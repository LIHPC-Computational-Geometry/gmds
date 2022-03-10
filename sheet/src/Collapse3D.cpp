/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Collapse3D.h>
#include <gmds/utils/LocalCellTopology.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Collapse3D::Collapse3D(Mesh* AMesh): m_mesh(AMesh)
{
    buildLocalN2R();
}
/*----------------------------------------------------------------------------*/
Collapse3D::~Collapse3D()
{}
/*----------------------------------------------------------------------------*/
bool Collapse3D::isValid() const
{
    for(auto r_id:m_mesh->regions()){
        Region r=m_mesh->get<Region>(r_id);
        if (r.type()!=GMDS_HEX){
            return false;
        }
    }
    return (m_mesh->getModel().has(DIM3) &&
            m_mesh->getModel().has(R)&&
            m_mesh->getModel().has(N)&&
            m_mesh->getModel().has(R2N));
}
/*----------------------------------------------------------------------------*/
void Collapse3D::execute(const gmds::TCellID AN1, const gmds::TCellID AN2)
{
}
/*----------------------------------------------------------------------------*/
void Collapse3D::buildLocalN2R() {
    m_N2R.clear();
    std::vector<TCellID> n_ids;
    for(auto r_id:m_mesh->regions()){
        Region r=m_mesh->get<Region>(r_id);
        r.getIDs<Node>(n_ids);
        for(auto n_id:n_ids){
            m_N2R[n_id].push_back(r_id);
        }
    }
}
