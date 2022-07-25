/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Pillow2D.h>
#include <gmds/utils/LocalCellTopology.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Pillow2D::Pillow2D(Mesh* AMesh, const bool AEscapeBnd)
        :Operator2D(AMesh), m_escape_bnd(AEscapeBnd)
{}
/*----------------------------------------------------------------------------*/
Pillow2D::~Pillow2D()
{}
/*----------------------------------------------------------------------------*/
void Pillow2D::setBoundaryBehaviour(const bool AEscape) {
    m_escape_bnd=AEscape;
}
/*----------------------------------------------------------------------------*/
bool Pillow2D::checkCurve(std::vector<VirtualEdge> &AEdges) {

    //the curve is valid if:
    // - each node must belong to exactly two edges of AEdge or be on the
    //   boundary
    //  AND
    // - it partitions quads in exactly two sets

    //-------------------------------------------------------------
    // (1) Check node property
    //-------------------------------------------------------------
     std::map<VirtualEdge::EdgeID, bool> done;
    for(auto e:AEdges){
        done[e.getID()]=false;
    }

    std::vector<VirtualEdge> front;
    front.push_back(AEdges[0]);

    while(!front.empty()){
        VirtualEdge seed = front.back();
        front.pop_back();
        std::vector<TCellID> seed_nodes = seed.node_ids();

        done[seed.getID()]=true;

        //we ge the adjacent edges which are not yet handled
        // for each node of the seed, we store the number of edges sharing it
        int found_adj[2]={0,0};
        for(auto i_e=0;i_e<AEdges.size() &&
                found_adj[0]<2 &&
                found_adj[1]<2 ; i_e++){
            //we check that we do not use the seed edge
            if(AEdges[i_e]==seed)
                continue;
            //the current edge must be not already treated
            if(done[AEdges[i_e].getID()])
                continue;
            //ok, cur is not yet handled. We ms
            std::vector<TCellID> cur_nodes = AEdges[i_e].node_ids();
            //we check if cur & seed share some point
            for(auto i=0; i<2;i++){
                TCellID si = seed_nodes[i];
                for(auto j=0;j<2;j++){
                    TCellID cj = cur_nodes[j];
                    if(si==cj){
                        //common point
                        done[AEdges[i_e].getID()]=true;
                        front.push_back(AEdges[i_e]);
                        found_adj[i]++;
                    }
                }
            }
        }
        //we got through all the unmarked faces sharing an edge with seed
        if(found_adj[0]>1 ||found_adj[1]>1){
            //non-manifold curve
            return false;
        }
    }

    //-------------------------------------------------------------
    // (2) Check the volume partitioning
    //-------------------------------------------------------------
    int mark_quad = m_mesh->newMark<Face>();
    //we need to color all the quads with exactly two colors
    int color = 0;
    std::vector<TCellID> adv_front;
    //we put the first quad into it
    TCellID first_quad_id = *(m_mesh->faces().begin());
    adv_front.push_back(first_quad_id);
    m_mesh->mark<Face>(first_quad_id,mark_quad);
    bool mark_all=false;
    do {
        while (!adv_front.empty()) {
            TCellID cur_quad  = adv_front.back();
            adv_front.pop_back();
            std::vector<TCellID> n_ids = m_mesh->get<Face>(cur_quad).getIDs<Node>();
            for (auto i_edge = 0; i_edge < 4; i_edge++) {
                TCellID n0 = n_ids[i_edge];
                TCellID n1 = n_ids[(i_edge+1)%4];

                //check if it is a surface face
                VirtualEdge cur_edge(n0,n1);
                if (std::find(AEdges.begin(), AEdges.end(), cur_edge) == AEdges.end()) {
                    //so it it not a separation edge, we get the adj face
                    std::vector<TCellID> adj_fac = getAdjacentFaces(cur_edge);
                    if (adj_fac.size() == 1) {
                        //means boundary edge
                        //nothing to do
                        ;
                    } else if (adj_fac.size() == 2) {
                        //we get the opposite face
                        TCellID opp_fac = NullID;
                        if (adj_fac[0] == cur_quad)
                            opp_fac = adj_fac[1];
                        else
                            opp_fac = adj_fac[0];

                        if (!m_mesh->isMarked<Face>(opp_fac, mark_quad)) {
                            adv_front.push_back(opp_fac);
                            m_mesh->mark<Face>(opp_fac, mark_quad);
                        }
                    }
                }
            }
        }
        //we have finished to advance in a part
        //Do we have cover all the mesh?
        bool found_unmarked=false;
        for(Mesh::faces_iterator it = m_mesh->faces_begin();
            it != m_mesh->faces_end(); ++it){
            TCellID i = *it;
            if(!m_mesh->isMarked<Face>(i,mark_quad))
                found_unmarked=true;
        }
        if (found_unmarked){
            color++;
            if(color>1){
                //means we have at least 3 parts
                return false;
            }
        }
        else{
            //all marked
            mark_all=true;
        }
    }
    while(!mark_all);
    //clean used marks
    m_mesh->negateMaskMark<Face>(mark_quad);
    m_mesh->freeMark<Face>(mark_quad);

    return true;
}
/*----------------------------------------------------------------------------*/
bool Pillow2D::execute(const std::vector<TCellID>& AFaceIDs, bool AWithCheck)
{
    //We need to get the set of faces that surrounds the set of hexes given
    //as an input

    // we mark all the cells of the pillow set with a mark
    int mark_pillow_set = m_mesh->newMark<Face>();
    for(auto f_id:AFaceIDs){
        m_mesh->mark<Face>(f_id,mark_pillow_set);
    }

    std::set<VirtualEdge> bnd_edges;

    //======================================================================
    // (1) We get the node bounding the set of faces to be shrunk
    //======================================================================
    for(auto f_id:AFaceIDs) {
        Face f = m_mesh->get<Face>(f_id);
        std::vector<TCellID> n_ids = f.getIDs<Node>();
        for(auto i_edge = 0 ; i_edge<4; i_edge++){
            VirtualEdge ei(n_ids[i_edge],n_ids[(i_edge+1)%4]);
            std::vector<TCellID> adj_fac = getAdjacentFaces(ei);
            if(adj_fac.size()==1){
                //means boundary face
                if(m_escape_bnd==false){
                    bnd_edges.insert(ei);
                }
            }
            else if(adj_fac.size()==2){
                //we get the opposite face
                TCellID opp_quad=NullID;
                if(adj_fac[0]==f_id)
                    opp_quad=adj_fac[1];
                else
                    opp_quad=adj_fac[0];

                if(!m_mesh->isMarked<Face>(opp_quad,mark_pillow_set)){
                    //means internal boundary face so
                    bnd_edges.insert(ei);
                }
            }
        }
    }
    std::vector<VirtualEdge> be(bnd_edges.begin(), bnd_edges.end());
    bool success =execute(be, AWithCheck);
    //======================================================================
    //we clean the used marks
    //======================================================================
    for(auto f_id:AFaceIDs){
        m_mesh->unmark<Face>(f_id,mark_pillow_set);
    }
    m_mesh->freeMark<Face>(mark_pillow_set);

    return success;
}
/*----------------------------------------------------------------------------*/
bool Pillow2D::execute(std::vector<VirtualEdge>& AEdges, bool AWithCheck){

    bool check = true;
    if(AWithCheck){
        check= checkCurve(AEdges);
    }
    if(!check){
        return false;
    }
    //TODO We must ensure that the set of edges given as an input defines
    // a simple-loop closed curve

    //Then the aim is to:
    // (1) extract all the nodes that are on the pillow interface
    // (2) define a set of quad to shrink (i.e. a side of the curve to pillow
    //     to.
    //======================================================================
    //We need to build a set of quad that we will use to perform the pillowing
    //======================================================================

    //======================================================================
    // (1) We get the node of the curve
    //======================================================================
    // And we build a local connectivity from  the shrink nodes to the edges
    std::map<TCellID , std::vector<VirtualEdge> > n2e;
    std::set<TCellID> shrink_nodes;
    for(auto e:AEdges){
        std::vector<TCellID> enids = e.node_ids();
        shrink_nodes.insert(enids.begin(),enids.end());
        for(auto ni:enids){
            n2e[ni].push_back(e);
        }
    }

    //we mark all the shrink nodes
    int mark_shrink_nodes = m_mesh->newMark<Node>();
    for(auto i:shrink_nodes){
        m_mesh->mark<Node>(i,mark_shrink_nodes);
    }

    //Edges are now oriented in a common manner
    orient(AEdges);

    //now we check if an edge is on the boundary, if yes the orientation
    //must be done inward. If not, we reorder all the curve
    bool found_bnd_edge=false;

    for(auto i= 0;i<AEdges.size() && !found_bnd_edge;i++) {
        if (isABoundaryEdge(AEdges[i])) {
            found_bnd_edge = true;
            //check this orientation
            //We get the single face adjacent to this edge
            //if the edge is oriented towards the face we're good
            //otherwise, we inverse
            bool must_inverse = false;
            TCellID adj_quad_id = getAdjacentFaces(AEdges[i])[0];

            math::Point p0 = m_mesh->get<Node>(AEdges[i].node_ids()[0]).point();
            math::Point p1 = m_mesh->get<Node>(AEdges[i].node_ids()[1]).point();
            Face adj_quad = m_mesh->get<Face>(adj_quad_id);
            math::Vector3d quad_normal = adj_quad.normal();
            math::Vector3d to_quad=adj_quad.center()-p0;
            math::Vector v01_normal = (p1-p0).cross(quad_normal);
            if(v01_normal.dot(to_quad)<0){
                must_inverse=true;
            }

            if (must_inverse) {
                for (auto j = 0; j < AEdges.size(); j++) {
                    AEdges[j].reverse();
                }
            }
        }
    }

    //Now we get the inward quads, as the curve is oriented in a consistent manner
    std::set<TCellID> boundary_shrink_quads;
    std::set<TCellID> boundary_out_quads;
    int mark_shrink_quads = m_mesh->newMark<Face>();
    int mark_out_quads = m_mesh->newMark<Face>();

    for(auto i= 0;i<AEdges.size();i++) {
        std::vector<TCellID> adj_quad_ids = getAdjacentFaces(AEdges[i]);

        if (adj_quad_ids.size()==1){
            boundary_shrink_quads.insert(adj_quad_ids[0]);
            m_mesh->mark<Face>(adj_quad_ids[0],mark_shrink_quads);
        }
        else if(adj_quad_ids.size()==2){
            TCellID  adj_quad_id = adj_quad_ids[0];
            Face adj_quad = m_mesh->get<Face>(adj_quad_id);
            math::Point p0 = m_mesh->get<Node>(AEdges[i].oriented_node_ids()[0]).point();
            math::Point p1 = m_mesh->get<Node>(AEdges[i].oriented_node_ids()[1]).point();
            math::Vector3d quad_normal = adj_quad.normal();
            math::Vector3d to_quad=adj_quad.center()-p0;
            math::Vector v01_normal =(p1-p0).cross(quad_normal);
            if(v01_normal.dot(to_quad)<0){
                boundary_shrink_quads.insert(adj_quad_ids[1]);
                m_mesh->mark<Face>(adj_quad_ids[1],mark_shrink_quads);
                boundary_out_quads.insert(adj_quad_ids[0]);
                m_mesh->mark<Face>(adj_quad_ids[0],mark_out_quads);
            }
            else{
                boundary_shrink_quads.insert(adj_quad_ids[0]);
                m_mesh->mark<Face>(adj_quad_ids[0],mark_shrink_quads);
                boundary_out_quads.insert(adj_quad_ids[1]);
                m_mesh->mark<Face>(adj_quad_ids[1],mark_out_quads);
            }
        }
        else {
            throw GMDSException ("ERROR: Adjacency E2F is wrong");
        }
    }
    //======================================================================
    // All the faces of the shrink set that have an edge in AEdges are now
    // in the shrink set. Those sharing only a vertex with such
    // edges are missing, we have now to add them

    for(auto n_id:shrink_nodes){
        std::vector<TCellID> q_ids = m_mesh->get<Node>(n_id).getIDs<Face>();
        bool all_done = true;
        do{
            all_done=true;
            for(auto q:q_ids){
                if(!m_mesh->isMarked<Face>(q, mark_out_quads) &&
                   !m_mesh->isMarked<Face>(q, mark_shrink_quads)){
                    all_done=false;
                    //we get faces sharing an edge with q. They could
                    //be unmarked, in or out the shrink set but not both
                    std::vector<TCellID> f_ids = getFacesSharingAnEdge(q);
                    bool found_in=false;
                    bool found_out=false;
                    for(auto f:f_ids){
                        if(m_mesh->isMarked<Face>(f,mark_out_quads)){
                            found_out=true;
                        }
                        if(m_mesh->isMarked<Face>(f,mark_shrink_quads)){
                            found_in=true;
                        }
                    }
                    if(found_in && found_out){
                        throw GMDSException("ERROR in the IN/OUT process");
                    }
                    if(found_in){
                        m_mesh->mark<Face>(q,mark_shrink_quads);
                        boundary_shrink_quads.insert(q);
                    }
                    if(found_out){
                        m_mesh->mark<Face>(q,mark_out_quads);
                        boundary_out_quads.insert(q);
                    }
                }
            }
        }
        while (!all_done);
    }
    //======================================================================
    //We have the edges and the shrink set
    //======================================================================
    // (2) We build the pillow layer of quads
    //======================================================================
    // (2.1) For each boundary node, we build its twin to be connected to
    std::map<TCellID ,TCellID > n2n;
    for(auto n_id:shrink_nodes){
        math::Point p = m_mesh->get<Node>(n_id).point();
        Node sibling = m_mesh->newNode(p);
        // We change node positions now
        //For that, we need to know if the node is on the mesh boundary, and
        //if it is, we have to deal with multiple cases

        // We change node positions now
        std::vector<TCellID> adj_f =  m_mesh->get<Node>(n_id).getIDs<Face>();
        std::vector<TCellID> adj_shrink_f;
        for(auto q_id:adj_f){
            if(m_mesh->isMarked<Face>(q_id,mark_shrink_quads)){
                adj_shrink_f.push_back(q_id);
            }
        }
        if(adj_shrink_f.empty()){
            throw GMDSException("ERROR: Empty Shrink set");
        }

        if(isABoundaryNode(n_id)){
            // we get the boundary edges (in the geometric meaning)
            // so we go through the adjacent faces
            std::vector<TCellID> adj_bnd_edges[2];
            for(auto q_id:adj_shrink_f) {
                Face q = m_mesh->get<Face>(q_id);
                std::vector<TCellID> q_nodes = q.getIDs<Node>();
                for (auto i_edge = 0; i_edge < 4; i_edge++) {
                    //We check if this edge contains n_id first and is a bnd edge
                    if((q_nodes[i_edge]==n_id ||
                            q_nodes[(i_edge+1)%4]==n_id )  &&
                       isABoundaryEdge(q_nodes[i_edge],
                                       q_nodes[(i_edge+1)%4])) {

                        adj_bnd_edges[0].push_back(q_nodes[i_edge]);
                        adj_bnd_edges[1].push_back(q_nodes[(i_edge+1)%4]);
                    }
                }
            }

            math::Point c(0, 0, 0);
            for (auto i = 0; i < adj_bnd_edges[0].size(); i++) {
                math::Point face_center = 0.25 * (m_mesh->get<Node>(adj_bnd_edges[0][i]).point() +
                                                  m_mesh->get<Node>(adj_bnd_edges[1][i]).point());
                c = c + face_center;
            }
            c = (1.0 / adj_bnd_edges[0].size()) * c;
            //we only move the input bnd node towards the shrink set
            m_mesh->get<Node>(n_id).setPoint(0.95 * p + 0.05 * c);


        }
        else
        {

            std::vector<math::Point> shrink_centers;

            for(auto f_id:adj_shrink_f){
                shrink_centers.push_back(m_mesh->get<Face>(f_id).center());
            }

            //in-volume node
            //We get the centers of the adjacent hexes that belongs to
            //the shrink set
            math::Point c(0,0,0);
            for (auto ci:shrink_centers){
                c = c+ci;
            }
            c = (1.0/shrink_centers.size())*c;
            //we only move the input bnd node towards the shrink set
            m_mesh->get<Node>(n_id).setPoint(0.95*p+0.05*c);

        }
        n2n[n_id]=sibling.id();
    }
    //======================================================================
    //======================================================================
    // (2.2) For each edge, we build the corresponding quad
    //======================================================================
    std::set<TCellID > n2update;
    for(auto e:AEdges) {
        std::vector<TCellID> e_node_ids = e.node_ids();

        //it's a boundary edge of the shrink set
        //We create a new quad
        m_mesh->newQuad(e_node_ids[0],
                        e_node_ids[1],
                        n2n[e_node_ids[0]],
                        n2n[e_node_ids[1]]);
        n2update.insert(e_node_ids[0]);
        n2update.insert(e_node_ids[1]);
    }
    //We've got a topologic inconsistency to solve
    for(auto old_id: n2update){
        //shrink nodes have been duplicated but some adjacent quads
        //should now be connected to their twin and not to them.

        //we accees to all those quads
        std::vector<TCellID> old_fids =  m_mesh->get<Node>(old_id).getIDs<Face>();
        for(auto f_id:old_fids) {
            if(m_mesh->isMarked<Face>(f_id,mark_out_quads)) {
                //means we have to change the R2N connectivity for this region
                Face f = m_mesh->get<Face>(f_id);
                f.replace<Node>(old_id, n2n[old_id]);
            }
        }
    }


    //We clean and release boolean marks
    for(auto i:shrink_nodes){
        m_mesh->unmark<Node>(i,mark_shrink_nodes);
    }
    for(auto i:boundary_shrink_quads) {
        m_mesh->unmark<Face>(i, mark_shrink_quads);
    }
    for(auto i:boundary_out_quads) {
        m_mesh->unmark<Face>(i, mark_out_quads);
    }
    m_mesh->freeMark<Node>(mark_shrink_nodes);
    m_mesh->freeMark<Face>(mark_shrink_quads);
    m_mesh->freeMark<Face>(mark_out_quads);
    return true;
}
/*----------------------------------------------------------------------------*/
void Pillow2D::orient(std::vector<VirtualEdge> &AEdges) {
    std::map<VirtualEdge::EdgeID, bool> done;
    for(auto e:AEdges){
        done[e.getID()]=false;
    }

    std::vector<VirtualEdge> front;
    front.push_back(AEdges[0]);

    while(!front.empty()){
        VirtualEdge seed = front.back();
        front.pop_back();
        std::vector<TCellID> seed_nodes = seed.oriented_node_ids();

        done[seed.getID()]=true;

        //we ge the adjacent edges which are not yet handled
        int found_bnd=0;
        for(auto i_e=0;i_e<AEdges.size() && found_bnd<2; i_e++){
            //we check that we do not use the seed edge
            if(AEdges[i_e]==seed)
                continue;
            //the current edge must be not already treated
            if(done[AEdges[i_e].getID()])
                continue;
            //ok, cur is not yet handled. We ms
            std::vector<TCellID> cur_nodes = AEdges[i_e].oriented_node_ids();
            //we check if cur & seed share some point
            if(seed_nodes[0]==cur_nodes[0] || seed_nodes[0]==cur_nodes[1] ||
               seed_nodes[1]==cur_nodes[0] || seed_nodes[1]==cur_nodes[1] )
            {
                //means the current and seed edges share a point
                found_bnd++;
                done[AEdges[i_e].getID()]=true;
                front.push_back(AEdges[i_e]);
                //we check now if the orientation of current must be
                //reversed
                if(seed_nodes[0]==cur_nodes[0] ||  seed_nodes[1]==cur_nodes[1] )
                    AEdges[i_e].reverse();
            }
        }
    }

    //done, reoriented
}
/*----------------------------------------------------------------------------*/
bool Pillow2D::isABoundaryNode(const TCellID ANodeID) {
    //we get adjacent faces
    std::vector<TCellID> adj_faces =  m_mesh->get<Node>(ANodeID).getIDs<Face>();
    // For each ffcee, we get the edges containing the ANodeID
    for(auto f_id:adj_faces){
        Face f =m_mesh->get<Face>(f_id);
        TCellID  adj_n0, adj_n1;
        f.getAdjacentNodes(ANodeID, adj_n0, adj_n1);
        if(isABoundaryEdge(ANodeID,adj_n0) || isABoundaryEdge(ANodeID,adj_n0) )
            return true;
    }
    return false;
}
/*----------------------------------------------------------------------------*/
bool Pillow2D::isABoundaryEdge(const gmds::VirtualEdge &AE){
    return isABoundaryEdge(AE.first(), AE.second());
}
/*----------------------------------------------------------------------------*/
bool Pillow2D::isABoundaryEdge(const TCellID AN1, const TCellID AN2)
{
    return getAdjacentFaces(AN1,AN2).size()==1;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Pillow2D::getFacesSharingAnEdge(const TCellID AFaceID)
{

    std::vector<TCellID> adj_faces;
    Face f = m_mesh->get<Face>(AFaceID);
    std::vector<TCellID> n_ids = f.getIDs<Node>();
    for(auto i = 0 ; i< 4; i++){
        std::vector<TCellID> adj_f = getAdjacentFaces( n_ids[i],
                                                       n_ids[(i+1)%4]);
        if(adj_f.size()==1){
            //means boundary face
            //nothing to do
            ;
        }
        else if(adj_f.size()==2){
            //we get the opposite regions
            TCellID opp_face=NullID;
            if(adj_f[0]==AFaceID)
                opp_face=adj_f[1];
            else
                opp_face=adj_f[0];

            adj_faces.push_back(opp_face);
        }
    }
    return adj_faces;
}
/*----------------------------------------------------------------------------*/
