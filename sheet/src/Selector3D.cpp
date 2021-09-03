/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Selector3D.h>
#include <gmds/utils/LocalCellTopology.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Selector3D::Selector3D(Mesh* AMesh)
        :m_mesh(AMesh)
{
    buildLocalN2R();
}
/*----------------------------------------------------------------------------*/
Selector3D::~Selector3D()
{}
/*----------------------------------------------------------------------------*/
bool Selector3D::isValid() const
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
void Selector3D::execute(const gmds::TCellID AN1, const gmds::TCellID AN2)
{
    //We start from a valid edge
    if(!isAnEdge(AN1,AN2))
        return;

    //We clean the sheet set
    m_sheet_hexes.clear();
    //We use a mark for knowing sheet hexes
    int mark_sheet = m_mesh->newMark<Region>();
    std::vector<TCellID> initial_hexes = getHexesSharingEdge(AN1, AN2);
    std::vector<HexSheetInfo> advancing_set;
    //We mark all the regions as being in the sheet
    for(auto h_id:initial_hexes){
        m_mesh->mark<Region>(h_id,mark_sheet);
        HexSheetInfo info;
        info.hex_id=h_id;
        info.e1=AN1;
        info.e2=AN2;
        advancing_set.push_back(info);
    }

    while(!advancing_set.empty()){
        HexSheetInfo current_info = advancing_set.back();
        advancing_set.pop_back();
        Region current_hex = m_mesh->get<Region>(current_info.hex_id);
        m_sheet_hexes.push_back(current_hex.id());

        std::vector<TCellID> current_nodes = current_hex.getIDs<Node>();
        //We get the local index of the input edge
        int local_from_edge_id = getLocalEdgeIndex(current_info);
        //We get opposite edge in the hex to propagate
        for(auto i=0;i<3;i++){
            size_t opp_edge_i= LocalHexTopology::OppositeEdges[local_from_edge_id][i];
            //We get the global node id for the opp edge
            TCellID opp_n0 = current_nodes[LocalHexTopology::E2N[opp_edge_i][0]];
            TCellID opp_n1 = current_nodes[LocalHexTopology::E2N[opp_edge_i][1]];
            std::vector<TCellID> adj_hex_01 = getHexesSharingEdge(opp_n0, opp_n1);
            for(auto opp_h:adj_hex_01){
                if(!m_mesh->isMarked<Region>(opp_h,mark_sheet)) {
                    HexSheetInfo opp_info;
                    opp_info.hex_id=opp_h;
                    opp_info.e1=opp_n0;
                    opp_info.e2=opp_n1;
                    advancing_set.push_back(opp_info);
                    m_mesh->mark<Region>(opp_h,mark_sheet);

                }
            }
        }
    }
    //We clean the used mark
    for(auto h_id:m_sheet_hexes){
        m_mesh->unmark<Region>(h_id,mark_sheet);
    }
    //and we free the mark
    m_mesh->freeMark<Region>(mark_sheet);
}
/*----------------------------------------------------------------------------*/
void Selector3D::buildLocalN2R() {
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
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Selector3D::getSheetCells() const {
    return m_sheet_hexes;
}
/*----------------------------------------------------------------------------*/
void Selector3D::getSheetCells(std::vector<TCellID>& AHexes) const {
    AHexes=m_sheet_hexes;
}
/*----------------------------------------------------------------------------*/
bool Selector3D::isAnEdge(const TCellID AN1, const TCellID AN2) {

    std::vector<TCellID> adj_reg_1 = m_N2R[AN1];
    for(auto  r_id:adj_reg_1) {
        Region h = m_mesh->get<Region>(r_id);
        std::vector<TCellID> hn_ids = h.getIDs<Node>();
        for (auto i = 0; i < 12; i++) {
            TCellID glob_n0 = hn_ids[LocalHexTopology::E2N[i][0]];
            TCellID glob_n1 = hn_ids[LocalHexTopology::E2N[i][1]];
            if ((glob_n0 == AN1 && glob_n1 == AN2) ||
                (glob_n0 == AN2 && glob_n1 == AN1)) {
                return true;
            }
        }
    }
    return false;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Selector3D::getHexesSharingEdge(const TCellID AN1,
                                                          const TCellID AN2) {
    std::vector<TCellID> h1 = m_N2R[AN1];
    std::vector<TCellID> h2 = m_N2R[AN2];
    std::set<TCellID> h12(h1.begin(), h1.end());
    h12.insert(h2.begin(), h2.end());

    std::vector<TCellID> res;
    //Now we check if we [AN1,AN2] is an edge in every hex of h12.
    for(auto h:h12){
        std::vector<TCellID > hn = m_mesh->get<Region>(h).getIDs<Node>();
        bool found_edge = false;
        for(auto i=0;i<12;i++){
            TCellID glob_n0 = hn[LocalHexTopology::E2N[i][0]];
            TCellID glob_n1 = hn[LocalHexTopology::E2N[i][1]];
            if( (glob_n0==AN1 && glob_n1==AN2) ||
                (glob_n0==AN2 && glob_n1==AN1)){
                found_edge=true;
            }
        }
        if(found_edge){
            res.push_back(h);
        }
    }
    return res;
}
/*----------------------------------------------------------------------------*/
int Selector3D::getLocalEdgeIndex(const HexSheetInfo &AInfo)
{
    Region r = m_mesh->get<Region>(AInfo.hex_id);
    std::vector<TCellID> nids = r.getIDs<Node>();
    TCellID e1 = AInfo.e1;
    TCellID e2 = AInfo.e2;
    for(auto i_edge=0; i_edge<12; i_edge++){
        TCellID ei_n1 = nids[LocalHexTopology::E2N[i_edge][0]];
        TCellID ei_n2 = nids[LocalHexTopology::E2N[i_edge][1]];
        if( (ei_n1==e1 && ei_n2==e2) ||
            (ei_n1==e2 && ei_n2==e1)) {
            return i_edge;
        }
    }

    throw GMDSException("Error in Selector3D::getLocalEdgeIndex");
}