/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Pillow3D.h>
#include <gmds/utils/LocalCellTopology.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Pillow3D::Pillow3D(Mesh* AMesh, const bool AEscapeBnd)
        :m_mesh(AMesh), m_escape_bnd(AEscapeBnd)
{}
/*----------------------------------------------------------------------------*/
Pillow3D::~Pillow3D()
{}
/*----------------------------------------------------------------------------*/
void Pillow3D::setBoundaryBehaviour(const bool AEscape) {
    m_escape_bnd=AEscape;
}
/*----------------------------------------------------------------------------*/
bool Pillow3D::isValid() const
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
bool Pillow3D::checkSurface(std::vector<VirtualFace> &AFaces)
{
    //the surface is valid if:
    // - each edge must belong to exactly two faces of AFaces or be on the
    //   boundary
    //  AND
    // - it partitions hexes in exactly two sets

    //-------------------------------------------------------------
    // (1) Check edge property
    //-------------------------------------------------------------
     std::map<VirtualFace::FaceID, bool> done;
    for(auto f:AFaces){
        done[f.getID()]=false;
    }

    std::vector<VirtualFace> front;
    front.push_back(AFaces[0]);

    while(!front.empty()){
        VirtualFace seed = front.back();
        front.pop_back();
        std::vector<TCellID> seed_nodes = seed.node_ids();

        done[seed.getID()]=true;

        //we ge the adjacent faces which are not yet handled
        // for each edge of the seed, we store the number of faces sharing the edge
        int found_adj[4]={0,0,0,0};
        for(auto i_f=0;i_f<AFaces.size() &&
                found_adj[0]<2 &&
                found_adj[1]<2 &&
                found_adj[2]<2 &&
                found_adj[3]<2; i_f++){
            //we check that we do not use the seed face
            if(AFaces[i_f]==seed)
                continue;
            //the current face must be not already treated
            if(done[AFaces[i_f].getID()])
                continue;
            //ok, cur is not yet handled. We ms
            std::vector<TCellID> cur_nodes = AFaces[i_f].node_ids();
            //we check if cur & seed share some edge
            for(auto i=0; i<4;i++){
                TCellID si = seed_nodes[i];
                TCellID sj = seed_nodes[(i+1)%4];
                for(auto j=0;j<4;j++){
                    TCellID ci = cur_nodes[j];
                    TCellID cj = cur_nodes[(j+1)%4];
                    if((si==ci && sj==cj) ||(si==cj && sj==ci)){
                        //common edge
                        done[AFaces[i_f].getID()]=true;
                        front.push_back(AFaces[i_f]);
                        found_adj[i]++;
                    }
                }
            }
        }
        //we got through all the unmarked faces sharing an edge with seed
        if(found_adj[0]>1 ||found_adj[1]>1 ||found_adj[2]>1 ||found_adj[3]>1 ){
            //non-manifold curve
            return false;
        }
    }

    //-------------------------------------------------------------
    // (2) Check the volume partitioning
    //-------------------------------------------------------------
    int mark_hex = m_mesh->newMark<Region>();
    //we need to color all the hexes with exactly two colors
    int color = 0;
    std::vector<TCellID> adv_front;
    //we put the first hex into it
    TCellID first_hex_id = *(m_mesh->regions().begin());
    adv_front.push_back(first_hex_id);
    m_mesh->mark<Region>(first_hex_id,mark_hex);
    bool mark_all=false;
    do {
        while (!adv_front.empty()) {
            TCellID cur_hex = adv_front.back();
            adv_front.pop_back();
            std::vector<TCellID> n_ids = m_mesh->get<Region>(cur_hex).getIDs<Node>();
            std::vector<TCellID> face_nids;
            face_nids.resize(4);
            for (auto i_face = 0; i_face < 6; i_face++) {
                face_nids[0] = n_ids[LocalHexTopology::F2N[i_face][0]];
                face_nids[1] = n_ids[LocalHexTopology::F2N[i_face][1]];
                face_nids[2] = n_ids[LocalHexTopology::F2N[i_face][2]];
                face_nids[3] = n_ids[LocalHexTopology::F2N[i_face][3]];

                //check if it is a surface face
                VirtualFace cur_face(face_nids);
                if (std::find(AFaces.begin(), AFaces.end(), cur_face) == AFaces.end()) {
                    //so it it not a separation face, we get the adj region
                    std::vector<TCellID> adj_reg = getAdjacentRegions(face_nids[0],
                                                                      face_nids[1],
                                                                      face_nids[2],
                                                                      face_nids[3]);
                    if (adj_reg.size() == 1) {
                        //means boundary face
                        //nothing to do
                        ;
                    } else if (adj_reg.size() == 2) {
                        //we get the opposite regions
                        TCellID opp_reg = NullID;
                        if (adj_reg[0] == cur_hex)
                            opp_reg = adj_reg[1];
                        else
                            opp_reg = adj_reg[0];

                        if (!m_mesh->isMarked<Region>(opp_reg, mark_hex)) {
                            adv_front.push_back(opp_reg);
                            m_mesh->mark<Region>(opp_reg, mark_hex);
                        }
                    }
                }
            }
        }
        //we have finished to advance in a part
        //Do we have cover all the mesh?
        bool found_unmarked=false;
        for(Mesh::regions_iterator it = m_mesh->regions_begin(); it != m_mesh->regions_end(); ++it){
            TCellID i = *it;
            if(!m_mesh->isMarked<Region>(i,mark_hex))
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
    m_mesh->negateMaskMark<Region>(mark_hex);
    m_mesh->freeMark<Region>(mark_hex);

    return true;
}
/*----------------------------------------------------------------------------*/
bool Pillow3D::execute(const std::vector<TCellID>& ARegionIDs, bool AWithCheck)
{
    //We need to get the set of faces that surrounds the set of hexes given
    //as an input

    buildLocalN2R();

    // we mark all the cells of the pillow set with a mark
    int mark_pillow_set = m_mesh->newMark<Region>();
    for(auto r_id:ARegionIDs){
        m_mesh->mark<Region>(r_id,mark_pillow_set);
    }

    std::set<VirtualFace> bnd_faces;

    //======================================================================
    // (1) We get the node bounding the set of regions to be shrunk
    //======================================================================
    std::vector<TCellID> face_nids;
    face_nids.resize(4);
    for(auto r_id:ARegionIDs) {

        Region r = m_mesh->get<Region>(r_id);
        std::vector<TCellID> n_ids = r.getIDs<Node>();
        for(auto i_face = 0 ; i_face<6; i_face++){
            face_nids[0] = n_ids[LocalHexTopology::F2N[i_face][0]];
            face_nids[1] = n_ids[LocalHexTopology::F2N[i_face][1]];
            face_nids[2] = n_ids[LocalHexTopology::F2N[i_face][2]];
            face_nids[3] = n_ids[LocalHexTopology::F2N[i_face][3]];
            std::vector<TCellID> adj_reg = getAdjacentRegions(face_nids[0],
                                                              face_nids[1],
                                                              face_nids[2],
                                                              face_nids[3]);
            if(adj_reg.size()==1){
                //means boundary face
                if(m_escape_bnd==false){
                    VirtualFace f(face_nids);
                    bnd_faces.insert(f);
                }
            }
            else if(adj_reg.size()==2){
                //we get the opposite regions
                TCellID opp_reg=NullID;
                if(adj_reg[0]==r_id)
                    opp_reg=adj_reg[1];
                else
                    opp_reg=adj_reg[0];

                if(!m_mesh->isMarked<Region>(opp_reg,mark_pillow_set)){
                    //means internal boundary face so
                    VirtualFace f(face_nids);
                    bnd_faces.insert(f);
                }
            }
        }
    }
    std::vector<VirtualFace> bf(bnd_faces.begin(), bnd_faces.end());
    bool success =execute(bf, AWithCheck);
    //======================================================================
    //we clean the used marks
    //======================================================================
    for(auto r_id:ARegionIDs){
        m_mesh->unmark<Region>(r_id,mark_pillow_set);
    }
    m_mesh->freeMark<Region>(mark_pillow_set);

    return success;
}
/*----------------------------------------------------------------------------*/
void Pillow3D::orient(std::vector<VirtualFace> &AFaces) {

    std::map<VirtualFace::FaceID, bool> done;
    for(auto f:AFaces){
        done[f.getID()]=false;
    }

    std::vector<VirtualFace> front;
    front.push_back(AFaces[0]);

    while(!front.empty()){
        VirtualFace seed = front.back();
        front.pop_back();
        std::vector<TCellID> seed_nodes = seed.oriented_node_ids();

        done[seed.getID()]=true;

        //we ge the adjacent faces which are not yet handled
        int found_bnd=0;
        for(auto i_f=0;i_f<AFaces.size() && found_bnd<4; i_f++){
            //we check that we do not use the seed face
            if(AFaces[i_f]==seed)
                continue;
            //the current face must be not already treated
            if(done[AFaces[i_f].getID()])
                continue;
            //ok, cur is not yet handled. We ms
            std::vector<TCellID> cur_nodes = AFaces[i_f].oriented_node_ids();
            //we check if cur & seed share some edge
            for(auto i=0; i<4;i++){
                TCellID si = seed_nodes[i];
                TCellID sj = seed_nodes[(i+1)%4];
                for(auto j=0;j<4;j++){
                    TCellID ci = cur_nodes[j];
                    TCellID cj = cur_nodes[(j+1)%4];
                    if(si==ci && sj==cj){
                        //common edge, same order
                        //we reverse, bad orientation
                        AFaces[i_f].reverse();
                        done[AFaces[i_f].getID()]=true;
                        front.push_back(AFaces[i_f]);
                        found_bnd++;
                    }
                    else if(si==cj && sj==ci){
                        //common edge, reverse order, nothing to do
                        done[AFaces[i_f].getID()]=true;
                        front.push_back(AFaces[i_f]);
                        found_bnd++;
                    }
                }
            }
        }
    }

    //done, reoriented
}
/*----------------------------------------------------------------------------*/
bool Pillow3D::execute(std::vector<VirtualFace>& AFaces, bool AWithCheck){

    bool check = true;
    if(AWithCheck){
        check= checkSurface(AFaces);
    }
    if(!check){
        return false;
    }
    //TODO We must ensure that the set of faces given as an input defines
    // a surface

    //Then the aim is to:
    // (1) extract all the nodes that are on the pillow interface
    // (2) define a set of hex to shrink (i.e. a side of the surface to pillow
    //     to.
    //======================================================================
    //We need to build a set of hexes that we will use to perform the pillowing
    //======================================================================
    buildLocalN2R();
    //======================================================================
    // (1) We get the node of the surface
    //======================================================================
    // And we build a local connectivity from  the shrink nodes to the faces
    std::map<TCellID , std::vector<VirtualFace> > n2f;
    std::set<TCellID> shrink_nodes;
    for(auto f:AFaces){
        std::vector<TCellID> fnids = f.node_ids();
        shrink_nodes.insert(fnids.begin(),fnids.end());
        for(auto ni:fnids){
            n2f[ni].push_back(f);
        }
    }

    //we mark all the shrink nodes
    int mark_shrink_nodes = m_mesh->newMark<Node>();
    for(auto i:shrink_nodes){
        m_mesh->mark<Node>(i,mark_shrink_nodes);
    }

    //Faces are now oriented in a common manner
    orient(AFaces);

    //now we check if a face is on the boundary, if yes the orientation
    //must be done inward. If not, we reorder all the surface
    bool found_bnd_face=false;

    for(auto i= 0;i<AFaces.size() && !found_bnd_face;i++) {
        if (isABoundaryFace(AFaces[i])) {
            found_bnd_face = true;
            //check this orientation
            //We get the single region adjacent to this face
            //if the face is oriented towards the region we're good
            //otherwise, we inverse
            bool must_inverse = false;
            TCellID adj_hex_id = getAdjacentRegions(AFaces[i])[0];

            math::Point p0 = m_mesh->get<Node>(AFaces[i].node_ids()[0]).point();
            math::Point p1 = m_mesh->get<Node>(AFaces[i].node_ids()[1]).point();
            math::Point p2 = m_mesh->get<Node>(AFaces[i].node_ids()[2]).point();
            math::Triangle p012(p0,p1,p2);
            math::Vector3d to_hex=m_mesh->get<Region>(adj_hex_id).center()-p0;
            if(p012.getNormal().dot(to_hex)<0){
                must_inverse=true;
            }

            if (must_inverse) {
                for (auto j = 0; j < AFaces.size(); j++) {
                    AFaces[j].reverse();
                }
            }
        }
    }

    //Now we get the inward hexes, as the surface is oriented in a consistent manner
    std::set<TCellID> boundary_shrink_hexes;
    std::set<TCellID> boundary_out_hexes;
    int mark_shrink_hexes = m_mesh->newMark<Region>();
    int mark_out_hexes = m_mesh->newMark<Region>();

    for(auto i= 0;i<AFaces.size();i++) {
        std::vector<TCellID>adj_hex_ids = getAdjacentRegions(AFaces[i]);

        if (adj_hex_ids.size()==1){
            boundary_shrink_hexes.insert(adj_hex_ids[0]);
            m_mesh->mark<Region>(adj_hex_ids[0],mark_shrink_hexes);
        }
        else if(adj_hex_ids.size()==2){
            TCellID  adj_hex_id = adj_hex_ids[0];
            math::Point p0 = m_mesh->get<Node>(AFaces[i].oriented_node_ids()[0]).point();
            math::Point p1 = m_mesh->get<Node>(AFaces[i].oriented_node_ids()[1]).point();
            math::Point p2 = m_mesh->get<Node>(AFaces[i].oriented_node_ids()[2]).point();
            math::Triangle p012(p0,p1,p2);
            math::Vector3d to_hex=m_mesh->get<Region>(adj_hex_id).center()-p0;
            if(p012.getNormal().dot(to_hex)<0){
                boundary_shrink_hexes.insert(adj_hex_ids[1]);
                m_mesh->mark<Region>(adj_hex_ids[1],mark_shrink_hexes);
                boundary_out_hexes.insert(adj_hex_ids[0]);
                m_mesh->mark<Region>(adj_hex_ids[0],mark_out_hexes);
            }
            else{
                boundary_shrink_hexes.insert(adj_hex_ids[0]);
                m_mesh->mark<Region>(adj_hex_ids[0],mark_shrink_hexes);
                boundary_out_hexes.insert(adj_hex_ids[1]);
                m_mesh->mark<Region>(adj_hex_ids[1],mark_out_hexes);
            }
        }
        else {
            throw GMDSException ("ERROR: Adjacency F2R is wrong");
        }
    }
    //======================================================================
    // All the regions of the shrink set that have a face in AFaces are now
    // in the shrink set. Those sharing only a vertex or an edge with such
    // faces are missing, we have now to add them

    for(auto n_id:shrink_nodes){
        std::vector<TCellID> h_ids = m_N2R[n_id];
        bool all_done = true;
        do{
            all_done=true;
            for(auto h:h_ids){
                if(!m_mesh->isMarked<Region>(h,mark_out_hexes) &&
                   !m_mesh->isMarked<Region>(h,mark_shrink_hexes)){
                    all_done=false;
                    //we get regions sharing a face with h. They could
                    //be unmarked, in or out the shrink set but not both
                    std::vector<TCellID> r_ids = getRegionsSharingAFace(h);
                    bool found_in=false;
                    bool found_out=false;
                    for(auto r:r_ids){
                        if(m_mesh->isMarked<Region>(r,mark_out_hexes)){
                            found_out=true;
                        }
                        if(m_mesh->isMarked<Region>(r,mark_shrink_hexes)){
                            found_in=true;
                        }
                    }
                    if(found_in && found_out){
                        throw GMDSException("ERROR in the IN/OUT process");
                    }
                    if(found_in){
                        m_mesh->mark<Region>(h,mark_shrink_hexes);
                        boundary_shrink_hexes.insert(h);
                    }
                    if(found_out){
                        m_mesh->mark<Region>(h,mark_out_hexes);
                        boundary_out_hexes.insert(h);
                    }
                }
            }
        }
        while (!all_done);
    }
    //======================================================================
    //We have the faces and the shrink set


    //======================================================================
    // (2) We build the pillow layer of hexes
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
        std::vector<TCellID> adj_r = m_N2R[n_id];
        std::vector<TCellID> adj_shrink_r;
        for(auto h_id:adj_r){
            if(m_mesh->isMarked<Region>(h_id,mark_shrink_hexes)){
                adj_shrink_r.push_back(h_id);
            }
        }
        if(adj_shrink_r.empty()){
            throw GMDSException("ERROR: Empty Shrink set");
        }


        if(isABoundaryNode(n_id)){
            // we get the boundary faces (in the geometric meaning)
            // so we go through the adjacent regions
            std::vector<TCellID> adj_bnd_faces[4];
            for(auto h_id:adj_shrink_r) {
                Region h = m_mesh->get<Region>(h_id);
                std::vector<TCellID> h_nodes = h.getIDs<Node>();
                for (auto i_face = 0; i_face < 6; i_face++) {
                    //We check if this face contains n_id first and is a bnd face
                    if((h_nodes[LocalHexTopology::F2N[i_face][0]]==n_id ||
                        h_nodes[LocalHexTopology::F2N[i_face][1]]==n_id ||
                        h_nodes[LocalHexTopology::F2N[i_face][2]]==n_id ||
                        h_nodes[LocalHexTopology::F2N[i_face][3]]==n_id)  &&
                       isABoundaryFace(h_nodes[LocalHexTopology::F2N[i_face][0]],
                                       h_nodes[LocalHexTopology::F2N[i_face][1]],
                                       h_nodes[LocalHexTopology::F2N[i_face][2]],
                                       h_nodes[LocalHexTopology::F2N[i_face][3]])) {

                        adj_bnd_faces[0].push_back(h_nodes[LocalHexTopology::F2N[i_face][0]]);
                        adj_bnd_faces[1].push_back(h_nodes[LocalHexTopology::F2N[i_face][1]]);
                        adj_bnd_faces[2].push_back(h_nodes[LocalHexTopology::F2N[i_face][2]]);
                        adj_bnd_faces[3].push_back(h_nodes[LocalHexTopology::F2N[i_face][3]]);

                    }
                }
            }
            //We check if there exist an edge common to all the bnd faces
            //If yes, we shrink along this edge
            bool single_edge =false;
            TCellID single_edge_id1 = NullID;
            TCellID single_edge_id2 = NullID;
            if(adj_bnd_faces[0].size()>1) {
                for (auto i_edge=0; i_edge<4 && !single_edge; i_edge++) {
                    TCellID  edge_ext1 = adj_bnd_faces[i_edge][0];
                    TCellID  edge_ext2 = adj_bnd_faces[(i_edge+1)%4][0];

                    if(edge_ext1!=n_id && edge_ext2!=n_id)
                        continue;
                    bool found_in_all=true;
                    //for each other face, we look for this edge
                    for (auto i = 1; i < adj_bnd_faces[0].size()  && found_in_all; i++) {
                        //face i
                        found_in_all=false;
                        for (auto j_edge = 0; j_edge < 4 ; j_edge++) {

                            TCellID edge_ext1_i = adj_bnd_faces[j_edge][i];
                            TCellID edge_ext2_i = adj_bnd_faces[(j_edge + 1) % 4][i];
                            if (edge_ext1_i != n_id && edge_ext2_i != n_id)
                                continue;
                            if ((edge_ext1 == edge_ext1_i && edge_ext2 == edge_ext2_i) ||
                                (edge_ext1 == edge_ext2_i && edge_ext2 == edge_ext1_i)) {
                                found_in_all=true;
                            }
                        }

                    }
                    if (found_in_all) {
                        single_edge = true;
                        single_edge_id1 = edge_ext1;
                        single_edge_id2 = edge_ext2;
                    }

                }
            }
            //Otherwise, we build the center of mass of all the faces
            //and we shrink toward it.

            math::Point c(0, 0, 0);
            if(single_edge){
                c = 0.5 * (m_mesh->get<Node>(single_edge_id1).point() +
                           m_mesh->get<Node>(single_edge_id2).point());
            }
            else {
                for (auto i = 0; i < adj_bnd_faces[0].size(); i++) {
                    math::Point face_center = 0.25 * (m_mesh->get<Node>(adj_bnd_faces[0][i]).point() +
                                                      m_mesh->get<Node>(adj_bnd_faces[1][i]).point() +
                                                      m_mesh->get<Node>(adj_bnd_faces[2][i]).point() +
                                                      m_mesh->get<Node>(adj_bnd_faces[3][i]).point());
                    c = c + face_center;
                }
                c = (1.0 / adj_bnd_faces[0].size()) * c;
            }
            //we only move the input bnd node towards the shrink set
            m_mesh->get<Node>(n_id).setPoint(0.95 * p + 0.05 * c);


        }
        else
        {

            std::vector<math::Point> shrink_centers;

            for(auto h_id:adj_shrink_r){
                shrink_centers.push_back(m_mesh->get<Region>(h_id).center());
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
    // (2.2) For each face, we build the corresponding hex
    //======================================================================
    std::set<TCellID > n2update;
    for(auto f:AFaces){
        std::vector<TCellID> f_node_ids = f.node_ids();

        //it's a boundary face of the shrink set
        //We create a new hex
        m_mesh->newHex(f_node_ids[0],
                       f_node_ids[1],
                       f_node_ids[2],
                       f_node_ids[3],
                       n2n[f_node_ids[0]],
                       n2n[f_node_ids[1]],
                       n2n[f_node_ids[2]],
                       n2n[f_node_ids[3]]);
        n2update.insert(f_node_ids[0]);
        n2update.insert(f_node_ids[1]);
        n2update.insert(f_node_ids[2]);
        n2update.insert(f_node_ids[3]);



    }
    //We've got a topologic inconsistency to solve
    for(auto old_id: n2update){
        //shrink nodes have been duplicated but some adjacent hexes
        //should now be connected to their twin and not to them.

        //we accees to all those hexes
        std::vector<TCellID> old_rids = m_N2R[old_id];
        for(auto r_id:old_rids) {
            if(m_mesh->isMarked<Region>(r_id,mark_out_hexes)) {
                //means we have to change the R2N connectivity for this region
                Region r = m_mesh->get<Region>(r_id);
                r.replace<Node>(old_id, n2n[old_id]);
            }
        }
    }


    //We clean and release boolean marks
    for(auto i:shrink_nodes){
        m_mesh->unmark<Node>(i,mark_shrink_nodes);
    }
    for(auto i:boundary_shrink_hexes) {
        m_mesh->unmark<Region>(i, mark_shrink_hexes);
    }
    for(auto i:boundary_out_hexes) {
        m_mesh->unmark<Region>(i, mark_out_hexes);
    }
    m_mesh->freeMark<Node>(mark_shrink_nodes);
    m_mesh->freeMark<Region>(mark_shrink_hexes);
    m_mesh->freeMark<Region>(mark_out_hexes);
    return true;
}

/*----------------------------------------------------------------------------*/
void Pillow3D::buildLocalN2R() {
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
bool Pillow3D::isABoundaryNode(const TCellID ANodeID) {
    //we get adjacent region
    std::vector<TCellID> adj_hexes = m_N2R[ANodeID];
    // For each region, we get the faces containing the node
    for(auto h_id:adj_hexes){
        Region h = m_mesh->get<Region>(h_id);
        std::vector<TCellID> h_nodes = h.getIDs<Node>();
        for(auto i_face = 0;  i_face<6; i_face++) {
            if ((h_nodes[LocalHexTopology::F2N[i_face][0]]==ANodeID ||
                 h_nodes[LocalHexTopology::F2N[i_face][1]]==ANodeID ||
                 h_nodes[LocalHexTopology::F2N[i_face][2]]==ANodeID ||
                 h_nodes[LocalHexTopology::F2N[i_face][3]]==ANodeID)  &&
                isABoundaryFace(h_nodes[LocalHexTopology::F2N[i_face][0]],
                                h_nodes[LocalHexTopology::F2N[i_face][1]],
                                h_nodes[LocalHexTopology::F2N[i_face][2]],
                                h_nodes[LocalHexTopology::F2N[i_face][3]]))
                return true;
        }
    }
    // if at least one face is on the boundary, then the node is too
    return false;
}
/*----------------------------------------------------------------------------*/
bool Pillow3D::isABoundaryFace(const VirtualFace& AF){
    std::vector<TCellID> fids = AF.node_ids();
    return isABoundaryFace(fids[0], fids[1], fids[2], fids[3]);
}
/*----------------------------------------------------------------------------*/
bool Pillow3D::isABoundaryFace(const TCellID AN1, const TCellID AN2,
                               const TCellID AN3, const TCellID AN4)
{
    return getAdjacentRegions(AN1,AN2,AN3,AN4).size()==1;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Pillow3D::getAdjacentRegions(const TCellID AN1,
                                                  const TCellID AN2,
                                                  const TCellID AN3,
                                                  const TCellID AN4)
{
    std::vector<TCellID> r0 = m_N2R[AN1];
    std::vector<TCellID> r1 = m_N2R[AN2];
    std::vector<TCellID> r2 = m_N2R[AN3];
    std::vector<TCellID> r3 = m_N2R[AN4];
    std::vector<TCellID > r0123;
    //get the intersection
    for(auto i0:r0){
        for(auto i1:r1){
            if(i0==i1){
                for(auto i2:r2){
                    if(i0==i2){
                        for(auto i3:r3){
                            if(i0==i3)
                                r0123.push_back(i0);
                        }
                    }
                }
            }
        }
    }
    return r0123;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Pillow3D::getAdjacentRegions(const VirtualFace &AF) {
    return getAdjacentRegions(AF.node_ids()[0],
                              AF.node_ids()[1],
                              AF.node_ids()[2],
                              AF.node_ids()[3]);
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Pillow3D::getRegionsSharingAFace(const TCellID ARegionID)
{

    std::vector<TCellID> adj_regions;
    std::vector<TCellID> face_nids;
    face_nids.resize(4);
    Region r = m_mesh->get<Region>(ARegionID);
    std::vector<TCellID> n_ids = r.getIDs<Node>();
    for(auto i_face = 0 ; i_face<6; i_face++){
        face_nids[0] = n_ids[LocalHexTopology::F2N[i_face][0]];
        face_nids[1] = n_ids[LocalHexTopology::F2N[i_face][1]];
        face_nids[2] = n_ids[LocalHexTopology::F2N[i_face][2]];
        face_nids[3] = n_ids[LocalHexTopology::F2N[i_face][3]];
        std::vector<TCellID> adj_reg = getAdjacentRegions(face_nids[0],
                                                          face_nids[1],
                                                          face_nids[2],
                                                          face_nids[3]);
        if(adj_reg.size()==1){
            //means boundary face
            //nothing to do
            ;
        }
        else if(adj_reg.size()==2){
            //we get the opposite regions
            TCellID opp_reg=NullID;
            if(adj_reg[0]==ARegionID)
                opp_reg=adj_reg[1];
            else
                opp_reg=adj_reg[0];

            adj_regions.push_back(opp_reg);
        }
    }
    return adj_regions;
}
/*----------------------------------------------------------------------------*/
