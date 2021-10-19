/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Collapse2D.h>
#include <gmds/sheet/Selector2D.h>
#include <gmds/cad/GeomMeshLinker.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Collapse2D::Collapse2D(Mesh* AMesh, gmds::cad::GeomMeshLinker *ALinker)
: Operator2D(AMesh), m_geom_linker(ALinker)
{}
/*----------------------------------------------------------------------------*/
Collapse2D::~Collapse2D()
{}
/*----------------------------------------------------------------------------*/
void Collapse2D::execute(const gmds::TCellID AN1, const gmds::TCellID AN2)
{
    Selector2D selector(m_mesh);
    selector.execute(AN1,AN2);
    std::vector<VirtualEdge> traversed_edges = selector.getSheetTraversedEdges();
    std::vector<TCellID> traversed_quads = selector.getSheetCells();

    //We need to split node edges in two sets of each side of the collapse to
    //sheet
    int mark_quad_sheets = m_mesh->newMark<Face>();

    for(auto q_id:traversed_quads){
        m_mesh->mark<Face>(q_id,mark_quad_sheets);
    }

    //we are going to color each extremity of the traversed edges. While both
    //extremities of all the edges do not have colors, we keep going
    std::map<TCellID, int> color;
    //we store the opposite node(s) for each edge extremity. We get a vector
    //to handle self-intersecting and self-touching sheets.
    std::map<TCellID , std::vector<TCellID> > sheet_n2n;
    for(auto e:traversed_edges){
        color[e.first() ]=0;
        color[e.second()]=0;
        sheet_n2n[e.first()].push_back(e.second());
        sheet_n2n[e.second()].push_back(e.first());
    }

    bool keep_coloring = true;
    int current_color=0;

    while (keep_coloring){
        keep_coloring=false;
        TCellID seed_node_id = NullID;
        for(auto i=0; i<traversed_edges.size() &&!keep_coloring; i++){
            VirtualEdge e = traversed_edges[i];
            if (color[e.first()]==0){
                keep_coloring=true;
                seed_node_id=e.first();
            }
            else if (color[e.second()]==0){
                keep_coloring=true;
                seed_node_id=e.second();
            }
        }
        if(keep_coloring){
            //we pick a new color
            current_color++;

            //we color starting from the seed
            std::vector<TCellID> adv_front;
            adv_front.push_back(seed_node_id);
            while (!adv_front.empty()){
                TCellID cur_node_id = adv_front.back();
                adv_front.pop_back();
                color[cur_node_id]=current_color;

                //we pick the adjacent nodes on this side of the sheet
                Node cur_node = m_mesh->get<Node>(cur_node_id);
                std::vector<TCellID> all_current_face_ids = cur_node.getIDs<Face>();
                std::vector<TCellID> sheet_face_ids;
                for(auto f_id:all_current_face_ids){
                    if(m_mesh->isMarked<Face>(f_id,mark_quad_sheets)){
                        //means that it is a face on this side
                        sheet_face_ids.push_back(f_id);
                    }
                }
                std::vector<TCellID> adj_candidate_nodes;
                std::vector<TCellID> cur_opp_ids = sheet_n2n[cur_node_id];
                for(auto sheet_fid : sheet_face_ids){
                    Face shf = m_mesh->get<Face>(sheet_fid);
                    TCellID sh_n1, sh_n2;
                    shf.getAdjacentNodes(cur_node_id, sh_n1, sh_n2);
                    if(color[sh_n1]==0 &&
                       std::find(cur_opp_ids.begin(),cur_opp_ids.end(),sh_n1)==cur_opp_ids.end()){
                        //means sh_n1 is on the boudary
                        adv_front.push_back(sh_n1);
                    }
                    else if(color[sh_n2]==0 &&
                            std::find(cur_opp_ids.begin(),cur_opp_ids.end(),sh_n2)==cur_opp_ids.end()){
                        //means sh_n2 is on the boudary
                        adv_front.push_back(sh_n2);
                    }
                }
            }

        }
    }

    //now all the nodes of each boundary must have a color
    if (current_color!=2){
        std::string mess = "Color error for sheet collapsing: "+std::to_string(current_color);
        mess+="\n Only 2 colors handled, i.e. no self-intersecting and self-touching sheet";
        throw GMDSException(mess);
    }

    // Now we verified that each edge has its end point with different colors
    for(auto e:traversed_edges){
        if(color[e.first()]==color[e.second()]){
            std::string mess = "Color error: a traversed edge must have ";
            mess +="different color for its end points\n";
            mess +="Edge ["+std::to_string(e.first())+
                   ", "+std::to_string(e.second())+"]";
            throw GMDSException(mess);
        }
    }

    if(!checkGeometricClassification(traversed_edges, color)){
        std::cout<<"Impossible to collapse because of geometric classification";
        //We unmark the sheet faces that are still in the mesh
        for(auto q_id:traversed_quads){
            m_mesh->unmark<Face>(q_id,mark_quad_sheets);
        }
    }
    else{
        //We can collapse
        collapseAndReconnect(sheet_n2n, traversed_quads, mark_quad_sheets);
        //We don't have to unmark the sheet quads since they have been removed
    }

    //we free the used marks
    m_mesh->freeMark<Face>(mark_quad_sheets);

}
/*----------------------------------------------------------------------------*/
void Collapse2D::
collapseAndReconnect(const std::map<TCellID , std::vector<TCellID> >& AToCollapse,
                     const std::vector<TCellID>& ASheetQuads,
                     const int AMarkSheetQuads)
{
    //we store, which node has been treated, none at the beginning
    std::map<TCellID,bool> done;
    for(auto i:AToCollapse){
        done[i.first]=false;
    }

    //We remove all the quads in the sheet
    for(auto q_id: ASheetQuads){
        Face q = m_mesh->get<Face>(q_id);
        std::vector<Node> ns = q.get<Node>();
        for(auto n:ns){
            n.remove(q);
        }
        m_mesh->deleteFace(q);
    }
    //Now quads are removed, leaving some holes in the mesh
    //We connect them now
    for(auto item:AToCollapse){
        if(!done[item.first]){
            //means node n1 not collapsed
            TCellID n1_id = item.first;
            std::vector<TCellID> n2_ids = item.second;
            if(n2_ids.size()!=1){
                throw GMDSException("Collapse 2D only work for simple loop");
            }
            else{
                //simple case
                TCellID n2_id = n2_ids.front();
                Node n1 = m_mesh->get<Node>(n1_id);
                Node n2 = m_mesh->get<Node>(n2_id);

                done[n1_id]=true;
                done[n2_id]=true;
                //We check the geometric classification of each node
                int geo_dim1 = m_geom_linker->getGeomDim(n1);
                int geo_dim2 = m_geom_linker->getGeomDim(n2);

                TCellID to_keep;
                TCellID to_erase;
                if(geo_dim1==geo_dim2){
                    //Considering tests done before getting in this operation
                    //only two options are possible
                    // 1. dim = 1 on the same curve
                    // 2. dim 2 same surface
                    to_keep = n1_id;
                    to_erase = n2_id;
                    math::Point p1= n1.point();
                    math::Point p2= n2.point();

                    if(geo_dim1==cad::GeomMeshLinker::LINK_CURVE){
                        //On the same curve, same dim
                        int curve_id = m_geom_linker->getGeomId(n1);
                        cad::GeomManager* geo_model = m_geom_linker->geometry();
                        math::Point new_point = 0.5*(p1+p2);

                        geo_model->getCurve(curve_id)->project(new_point);
                        n1.setPoint(new_point);
                    }
                    else if(geo_dim1!=cad::GeomMeshLinker::LINK_SURFACE){
                        throw GMDSException("collapseAndReconnect: Unexpected configuration");
                    }
                }
                else if(geo_dim1<geo_dim2){
                    to_keep = n1_id;
                    to_erase = n2_id;
                }
                else{ // geo_dim1 > geo_dim2
                    to_keep = n2_id;
                    to_erase = n1_id;
                }

                Node n_keep = m_mesh->get<Node>(to_keep);
                Node n_erase = m_mesh->get<Node>(to_erase);
                //We get all the faces around n_erase and we update their connectivity
                std::vector<Face> faces_erase = n_erase.get<Face>();
                for(auto f:faces_erase){
                    f.replace(n_erase, n_keep);
                    n_keep.add(f);
                }
                //and we suppress n_erase
                m_mesh->deleteNode(n_erase);

            }
        }
    }

}
/*----------------------------------------------------------------------------*/
bool Collapse2D::
checkGeometricClassification(const std::vector<VirtualEdge>& AEdges,
                             const std::map<TCellID,int>& ASideColor)
{
    /* for each side of the sheet, we store the lowest dimension of geom entity
     * a node is classified on. */
    std::map<int, int> min_geom_entity_dim;
    std::map<int, int> min_geom_entity_id;

    // First we store default value
    for(auto sc:ASideColor){
        min_geom_entity_dim[sc.second]=5;
        min_geom_entity_id [sc.second]=-1;
    }

    for(auto e:AEdges){
        TCellID n1_id = e.first(), n2_id = e.second();
        int geom_dim_1  = m_geom_linker->getGeomDim<Node>(n1_id);
        int geom_dim_2  = m_geom_linker->getGeomDim<Node>(n2_id);

        //if classified on same dim geom entity
        if(geom_dim_1==geom_dim_2){
            if(geom_dim_1==cad::GeomMeshLinker::LINK_POINT){
                return false;
            }
            else if(geom_dim_1==cad::GeomMeshLinker::LINK_CURVE){
                // Means both extremities are on a curve.
                // Is it the same, if yes and that we are on the boudary of the
                // sheet, i.e the edge is adjacent to only one quad sheet, the
                // it's okay.
                int geom_id_1 = m_geom_linker->getGeomId<Node>(n1_id);
                int geom_id_2 = m_geom_linker->getGeomId<Node>(n2_id);
                if(geom_id_1!=geom_id_2) {
                    return false;
                }
                //same curve, we check that we are on the sheet extremity,
                // i.e the edge is on the boundary.
                std::vector<TCellID > f12 = getAdjacentFaces(n1_id,n2_id);
                if(f12.size()==2){
                    return false;
                }
            }
            // we are inside a surface, no specific constraint
        }

        //No constraint local to e, but a global constraint can pe possible
        //We update the colors for each side

        //We check the color of the first node
        int c1 = ASideColor.at(n1_id);
        if(min_geom_entity_dim[c1]>geom_dim_1){
            min_geom_entity_dim[c1]=geom_dim_1;
        }
        //... and the second node
        int c2 = ASideColor.at(n2_id);
        if(min_geom_entity_dim[c2]>geom_dim_2){
            min_geom_entity_dim[c2]=geom_dim_2;
        }
    }

    // At this point in this method, all the local configurations to an edge
    // have been handled.
    //Do we have??

    return true;
}