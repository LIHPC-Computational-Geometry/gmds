/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Selector2D.h>
#include <gmds/cad/GeomMeshLinker.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Selector2D::Selector2D(gmds::Mesh *AMesh)
        :Operator2D(AMesh)
{}
/*----------------------------------------------------------------------------*/
Selector2D::~Selector2D()
{}
/*----------------------------------------------------------------------------*/
void Selector2D::execute(const gmds::TCellID AN1, const gmds::TCellID AN2)
{
    //We start from a valid edge
    if(!isAnEdge(AN1,AN2))
        return;

    //We clean the sheet set
    m_sheet_quad_infos.clear();
    //We use a mark for knowing sheet hexes
    int mark_sheet = m_mesh->newMark<Face>();
    std::vector<TCellID> initial_quads = getAdjacentFaces(AN1, AN2);
    std::vector<QuadSheetInfo> advancing_set;
    //We mark all the regions as being in the sheet
    for(auto q_id:initial_quads){
        m_mesh->mark<Face>(q_id,mark_sheet);
        QuadSheetInfo info;
        info.quad_id=q_id;
        info.e1=AN1;
        info.e2=AN2;
        advancing_set.push_back(info);
    }

    while(!advancing_set.empty()){
        QuadSheetInfo current_info = advancing_set.back();
        advancing_set.pop_back();
        m_sheet_quad_infos.push_back(current_info);
        Face current_quad = m_mesh->get<Face>(current_info.quad_id);

        std::vector<TCellID> current_nodes = current_quad.getIDs<Node>();
        //We get the local index of the input edge
        int local_from_edge_id = getLocalEdgeIndex(current_info);
        //We get opposite edge in the quad to propagate
        size_t opp_edge_i= (local_from_edge_id+2)%4;
        //We get the global node id for the opp edge
        TCellID opp_n0 = current_nodes[opp_edge_i];
        TCellID opp_n1 = current_nodes[(opp_edge_i+1)%4];
        std::vector<TCellID> adj_quad_01 = getAdjacentFaces(opp_n0, opp_n1);
        for(auto opp_q:adj_quad_01){
            if(!m_mesh->isMarked<Face>(opp_q,mark_sheet)) {
                QuadSheetInfo opp_info;
                opp_info.quad_id=opp_q;
                opp_info.e1=opp_n0;
                opp_info.e2=opp_n1;
                advancing_set.push_back(opp_info);
                m_mesh->mark<Face>(opp_q,mark_sheet);

            }
        }

    }
    //We clean the used mark
    for(auto info:m_sheet_quad_infos){
        m_mesh->unmark<Face>(info.quad_id,mark_sheet);
    }
    //and we free the mark
    m_mesh->freeMark<Face>(mark_sheet);
}

/*----------------------------------------------------------------------------*/
std::vector<TCellID> Selector2D::getSheetCells()  const {
    std::vector<TCellID> cell_ids;
    cell_ids.reserve(m_sheet_quad_infos.size());
    for(auto info:m_sheet_quad_infos){
        cell_ids.push_back(info.quad_id);
    }
    return cell_ids;
}
/*----------------------------------------------------------------------------*/
void Selector2D::getSheetCells(std::vector<TCellID>& AQuads) const {
    AQuads.clear();
    AQuads.reserve(m_sheet_quad_infos.size());
    for(auto info:m_sheet_quad_infos){
        AQuads.push_back(info.quad_id);
    }
}
/*----------------------------------------------------------------------------*/

std::vector<VirtualEdge> Selector2D::getSheetTraversedEdges() const {
    std::set<VirtualEdge> set_edges;
    for(auto info:m_sheet_quad_infos){
        VirtualEdge e1(info.e1,info.e2);
        Face f = m_mesh->get<Face>(info.quad_id);
        std::vector<TCellID > fns = f.getIDs<Node>();
        int edge_local_idx = -1;
        for(auto i=0; i<4 && edge_local_idx==-1; i++){
            TCellID id1 = fns[i];
            TCellID id2 = fns[(i+1)%4];
            if((id1 == e1.first() && id2==e1.second()) ||
               (id2 == e1.first() && id1==e1.second())){
                edge_local_idx=i;
            }
        }
        if(edge_local_idx==-1){
            std::cout<<"ERRROR"<<std::endl;
        }
        VirtualEdge opp_e1(fns[(edge_local_idx+2)%4],
                           fns[(edge_local_idx+3)%4]);

        set_edges.insert(opp_e1);
        set_edges.insert(e1);
    }

    std::vector<VirtualEdge> edges(set_edges.begin(),set_edges.end());

    return edges;
}
/*----------------------------------------------------------------------------*/
int Selector2D::getLocalEdgeIndex(const QuadSheetInfo &AInfo)
{
    Face f = m_mesh->get<Face>(AInfo.quad_id);
    std::vector<TCellID> nids = f.getIDs<Node>();
    TCellID e1 = AInfo.e1;
    TCellID e2 = AInfo.e2;
    for(auto i_edge=0; i_edge<4; i_edge++){
        TCellID ei_n1 = nids[i_edge];
        TCellID ei_n2 = nids[(i_edge+1)%4];
        if( (ei_n1==e1 && ei_n2==e2) ||
            (ei_n1==e2 && ei_n2==e1)) {
            return i_edge;
        }
    }

    throw GMDSException("Error in Selector2D::getLocalEdgeIndex");
}