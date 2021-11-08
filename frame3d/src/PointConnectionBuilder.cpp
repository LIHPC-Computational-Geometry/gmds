// GMDS File Headers
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Triangle.h>
#include <gmds/math/Plane.h>
#include <gmds/utils/OrientedGraph.h>
/*---------------------------------------------------------------------------*/
// STL File Headers
#include <set>
#include <algorithm>
/*---------------------------------------------------------------------------*/
// Frame3D File Headers
#include <gmds/frame3d/PointConnectionBuilder.h>
using namespace gmds;
/*---------------------------------------------------------------------------*/
void computeCombinationsRec(std::vector<int>& cmb, int n, int p , int i, int k,
                            std::vector<std::vector<int> >& solutions)
{
    if (k == p) {
        solutions.push_back(cmb);
    }
    else if (i < n) {
        computeCombinationsRec(cmb,n,p,i+1,k,solutions);
        cmb[k] = i;
        computeCombinationsRec(cmb,n,p,i+1,k+1,solutions);
    }
}
/*---------------------------------------------------------------------------*/
// Returns in ASol the Cnp combinations of AP integer into [0;AN]
void computeCombinations(int AN, int AP,
                         std::vector<std::vector<int> >& ASol)
{
    std::vector<int> witness;
    witness.resize(3);
    computeCombinationsRec(witness, AN, AP,0,0, ASol);
}

/*---------------------------------------------------------------------------*/
PointConnectionBuilder::
PointConnectionBuilder(Mesh*                        AMesh,
                const std::vector<math::Point>&     APnts,
                const std::vector<math::Chart>&     ACharts,
                const std::vector<Cell::Data>&      AData,
                const std::vector<int>&             ATypes,
                const std::vector<int>&             AClass,
                const std::vector<int>&             ACurv,
                const std::vector<int>&             ASurf,
                const std::map<int, math::Vector3d>&ANormal):
m_mesh(AMesh),
m_pnt(APnts),
m_chart(ACharts),
m_mesh_data(AData),
m_type(ATypes),
m_classification(AClass),
m_curve(ACurv),
m_surface(ASurf),
m_normal(ANormal),
m_hexes(Mesh(MeshModel(DIM3 | R | F | N | R2N | R2F | F2N | F2R | N2F))),
m_dot_tolerance(0.75),
m_spacing(1)
{
    if(m_pnt.size()!=m_chart.size() ||
       m_pnt.size()!=m_type.size() ||
       m_pnt.size()!=m_classification.size()||
       m_pnt.size()!=m_surface.size()) {
        throw GMDSException("Incompatible size vector as input!");
    }

    m_used.resize(m_pnt.size(), 0);
    m_with_debug_info=false;
    m_output_dir = ".";
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
setDebugInfo(const bool &AWithDebug, const std::string& AOutputDir)
{
    m_with_debug_info = AWithDebug;
    m_output_dir = AOutputDir;
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::execute()
{


    //======================================================================
    // STEP 1 - PREFILTER ON THE POINT DISTANCE
    // SHOULD BE IMPROVED
    //======================================================================
    createDistanceFilter();

    //======================================================================
    // STEP 2 - Compute the exact relation between points and cells of the
    //          tetrahedral background mesh. Indeed, this relation is loose
    //          at the end of the previous algorithm (point generation). We
    //          only have a close tet for each point (the tet in which the
    //          point was generated)
    //======================================================================
    computeMeshAssociation();

    if(m_with_debug_info) {
        writeInput();
    }
        //======================================================================
    // STEP 3 - For each point, build best-fit oriented edges to "adjacent"
    //          points
    //======================================================================
    std::vector<std::vector<OrientedEdge> > oriented_edges_init;
    buildOrientedEdges(oriented_edges_init);


    if(m_with_debug_info) {
        writeEdges(oriented_edges_init, m_output_dir + "/EDGES_LOCAL.vtk");
    }
    //======================================================================
    // STEP 4 - Correction of the oriented-edges to create edges
    //======================================================================
    buildEdges(oriented_edges_init,m_edges);


    if(m_with_debug_info)
        writeEdges(m_edges, m_output_dir+"/EDGES_GLOBAL.vtk");

    //======================================================================
    // STEP 5 - For each stable point compute and store its hex-corners
    //======================================================================
    buildHexCorners(m_edges);

    //======================================================================
    // STEP 6 - Build stable hexahedral elts
    //======================================================================
    buildHexahedral();
    std::cout<<"Nb created hexes: "<<m_hexes.getNbHexahedra()<<std::endl;
    if(m_with_debug_info) {
        writeHexes();
    }
}


/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::createDistanceFilter()
{
    double max_distance = 2*m_spacing;

    for(auto i=0; i<m_pnt.size(); i++){
        math::Point pi = m_pnt[i];
        if(m_type[i]!=FRAME_SING){
            for(unsigned int j=i+1; j<m_pnt.size(); j++){
                math::Point pj = m_pnt[j];
                if(m_type[j]!=FRAME_SING){
                    if(pi.distance(pj)<max_distance){
                        m_filter[i].push_back(j);
                        m_filter[j].push_back(i);
                    }
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::computeMeshAssociation()
{
    for(auto i=0; i<m_pnt.size(); i++){
        computeMeshAssociation(i);
    }
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::computeMeshAssociation(const int AID)
{
    int geom_dim = m_classification[AID];
    TCellID start_tet_id = m_mesh_data[AID].id;
    double distance = m_spacing;

    bool global_surface = false;
    if(geom_dim==3){
        m_mesh_data[AID] = getRegionContaining(m_pnt[AID],
                                               start_tet_id);
    }
    else if (geom_dim==2){
        int surf_id = m_surface[AID];
        if(global_surface){
            m_mesh_data[AID]= Cell::Data(2,closestFace(m_pnt[AID],surf_id).id());
        }
        else{
            m_mesh_data[AID] = getBoundaryFaceContaining(m_pnt[AID],
                                                         start_tet_id,
                                                         distance,
                                                         surf_id);
        }

    }
    else if (geom_dim==1){
        int curv_id = m_curve[AID];
        if(global_surface){
            m_mesh_data[AID]= Cell::Data(1,closestEdge(m_pnt[AID],curv_id).id());
        }
        else{
            m_mesh_data[AID] = getBoundaryEdgeContaining(m_pnt[AID],
                                                         start_tet_id,
                                                         distance,
                                                         curv_id);
        }

    }
    // Otherwise we do nothing since it is useless for our algorithm
    // since oriented edges starting from points classified on curves and
    // vertices are not built on frame interpolation
}
/*---------------------------------------------------------------------------*/
Face PointConnectionBuilder::closestFace(math::Point& AP, const int ASurfID)
{
    Variable<int>* color = m_mesh->getVariable<int,GMDS_FACE>("BND_SURFACE_COLOR");
    Face closest_face;
    double closest_dist = 10000;
    for(auto f_id: m_mesh->faces()){
        Face f = m_mesh->get<Face>(f_id);
        if((*color)[f.id()]==ASurfID){
            std::vector<Node> f_nodes = f.get<Node>();
            math::Triangle t(f_nodes[0].getPoint(),
                             f_nodes[1].getPoint(),
                             f_nodes[2].getPoint());
            math::Point p = t.project(AP);
            if(p.distance2(AP)<closest_dist){
                closest_dist=p.distance2(AP);
                closest_face=f;
            }
        }
    }
    // The point is moved too!!

    std::vector<Node> closest_nodes = closest_face.get<Node>();
    math::Triangle t(closest_nodes[0].getPoint(),
                     closest_nodes[1].getPoint(),
                     closest_nodes[2].getPoint());

    AP = t.project(AP);
    return closest_face;

}

/*---------------------------------------------------------------------------*/
Edge PointConnectionBuilder::closestEdge(math::Point& AP, const int ACurvID)
{
    Variable<int>* color = m_mesh->getVariable<int,GMDS_EDGE>("BND_CURVE_COLOR");
    Edge closest_edge;
    double closest_dist = 10000;
    for(auto e_id: m_mesh->edges()){
        Edge e = m_mesh->get<Edge>(e_id);
        if((*color)[e.id()]==ACurvID){
            std::vector<Node> e_nodes = e.get<Node>();
            math::Segment s(e_nodes[0].getPoint(),
                            e_nodes[1].getPoint());
            math::Point p = s.project(AP);
            if(p.distance2(AP)<closest_dist){
                closest_dist=p.distance2(AP);
                closest_edge=e;
            }
        }
    }
    // The point is moved too!!
    std::vector<Node> closest_nodes = closest_edge.get<Node>();
    math::Segment s(closest_nodes[0].getPoint(),
                    closest_nodes[1].getPoint());

    AP = s.project(AP);
    return closest_edge;

}
/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::getOppositeFace(const TCellID ANodeID,
                                      const Region&  AR,
                                      Face & AOut)
{
    std::vector<Face> fs = AR.get<Face>();
    for(const auto& f:fs){
        std::vector<TCellID> f_nids = f.getIDs<Node>();

        if(f_nids[0]!=ANodeID && f_nids[1]!=ANodeID && f_nids[2]!=ANodeID){
            AOut = f;
            return true;
        }

    }
    return false;
}
/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::getOppositeRegion(const TCellID ANodeID,
                                        const Region&  AR,
                                        Region & AOut)
{
    Face opp_face;
    if(getOppositeFace(ANodeID, AR, opp_face)){
        std::vector<Region> f_r = opp_face.get<Region>();
        if(f_r.size()==1){
            return false;
        }
        if(f_r[0].id()==AR.id()){
            AOut =f_r[1];
            return true;
        }
        else{
            AOut = f_r[0];
            return true;
        }
    }

    return false;
}
/*---------------------------------------------------------------------------*/
Cell::Data PointConnectionBuilder::getRegionContaining(const math::Point& APnt,
                                                const TCellID      ARegionID)
{
    gmds::Region current_r = m_mesh->get<Region>(ARegionID);

    while(true){
        std::vector<Node> n = current_r.get<Node>();
        math::Point p[4] ={
                n[0].getPoint(), n[1].getPoint(),
                n[2].getPoint(), n[3].getPoint()
        };

        double coeff[4]={0, 0, 0, 0};

        math::Point::computeBarycentric(p[0], p[1], p[2], p[3], APnt,
                                        coeff[0],coeff[1],
                                        coeff[2],coeff[3]);

        if(coeff[0]>=0 && coeff[1]>=0 && coeff[2]>=0 &&  coeff[3]>=0)
            return Cell::Data(3,current_r.id());

        Region next_r;
        bool on_bnd;
        if(coeff[0]<0)
            on_bnd=getOppositeRegion(n[0].id(), current_r, next_r);
        else if(coeff[1]<0)
            on_bnd=getOppositeRegion(n[1].id(), current_r, next_r);
        else if(coeff[2]<0)
            on_bnd=getOppositeRegion(n[2].id(), current_r, next_r);
        else //(coeff[3]<0)
            on_bnd=getOppositeRegion(n[3].id(), current_r, next_r);

        if(!on_bnd){
            current_r=next_r;
        }
        else{
            return Cell::Data(3,current_r.id());
        }
    }
}
/*---------------------------------------------------------------------------*/
Cell::Data PointConnectionBuilder::
getBoundaryFaceContaining(math::Point&  AP,
                          const TCellID ATetID,
                          const double  ADistance,
                          const int     ASurfID)
{
    gmds::Region current_r = m_mesh->get<Region>(ATetID);

    std::set<TCellID> close_reg = getCloseRegionsFrom(AP,current_r,ADistance);
    std::set<TCellID> close_bnd_faces;

    Variable<int>* color = m_mesh->getVariable<int,GMDS_FACE>("BND_SURFACE_COLOR");
    for(auto i:close_reg){
        Region r = m_mesh->get<Region>(i);
        std::vector<Face> r_faces = r.get<Face>();
        for(const auto& f:r_faces){
            if((*color)[f.id()]==ASurfID){
                close_bnd_faces.insert(f.id());
            }
        }
    }

    Face closest_face;
    double closest_dist = 10000;
    for(auto i:close_bnd_faces){
        Face f = m_mesh->get<Face>(i);
        std::vector<Node> f_nodes = f.get<Node>();
        math::Triangle t(f_nodes[0].getPoint(),
                         f_nodes[1].getPoint(),
                         f_nodes[2].getPoint());
        math::Point p = t.project(AP);
        if(p.distance2(AP)<closest_dist){
            closest_dist=p.distance2(AP);
            closest_face=f;
        }
    }
    // The point is moved too!!

    std::vector<Node> closest_nodes = closest_face.get<Node>();
    math::Triangle t(closest_nodes[0].getPoint(),
                     closest_nodes[1].getPoint(),
                     closest_nodes[2].getPoint());

    AP = t.project(AP);
    return Cell::Data(2,closest_face.id());
}

/*---------------------------------------------------------------------------*/
Cell::Data PointConnectionBuilder::
getBoundaryEdgeContaining(math::Point&  AP,
                          const TCellID ATetID,
                          const double  ADistance,
                          const int     ACurvID)
{
    gmds::Region current_r = m_mesh->get<Region>(ATetID);

    std::set<TCellID> close_reg = getCloseRegionsFrom(AP,current_r,ADistance);
    std::set<TCellID> close_bnd_edges;

    Variable<int>* color = m_mesh->getVariable<int,GMDS_EDGE>("BND_CURVE_COLOR");
    for(auto i:close_reg){
        Region r = m_mesh->get<Region>(i);
        std::vector<Edge> r_edges = r.get<Edge>();
        for(const auto& e:r_edges){
            if((*color)[e.id()]==ACurvID){
                close_bnd_edges.insert(e.id());
            }
        }
    }

    Edge closest_edge;
    double closest_dist = 10000;
    for(auto i:close_bnd_edges){
        Edge e = m_mesh->get<Edge>(i);
        std::vector<Node> e_nodes = e.get<Node>();
        math::Segment s(e_nodes[0].getPoint(),
                        e_nodes[1].getPoint());
        math::Point p = s.project(AP);
        if(p.distance2(AP)<closest_dist){
            closest_dist=p.distance2(AP);
            closest_edge=e;
        }
    }

    std::vector<Node> closest_nodes = closest_edge.get<Node>();
    math::Segment s(closest_nodes[0].getPoint(),
                    closest_nodes[1].getPoint());

    // The point is moved too!!
    AP = s.project(AP);

    return Cell::Data(1,closest_edge.id());
}
/*---------------------------------------------------------------------------*/
std::set<TCellID>
PointConnectionBuilder::getCloseRegionsFrom(const math::Point& AFromPnt,
                                     const Region& AFromTet,
                                     const double AEpsilon)
{
    std::set<TCellID> region_ids;

    //starting from AFromTet, we look for all tet t such that the distance
    //to one face of t is < to AEpsilon
    std::vector<Region> candidates;
    candidates.push_back(AFromTet);
    std::set<TCellID> done;


    while(!candidates.empty()){
        Region c = candidates.back();
        candidates.pop_back();

        done.insert(c.id());

        std::vector<Face> fs = c.get<Face>();
        for(const auto & f:fs){
            std::vector<Node> n = f.get<Node>();
            math::Plane pl(n[0].getPoint(),
                           n[1].getPoint(),
                           n[2].getPoint());
            if(pl.project(AFromPnt).distance(AFromPnt)<AEpsilon){
                //means close to this face

                //Add the region
                region_ids.insert(c.id());

                //test to add the opposite region
                std::vector<TCellID> f_regs = f.getIDs<Region>();
                if(f_regs.size()>1){
                    //otherwise nothing to do we are on the bnd
                    TCellID opp_reg_id = (c.id()==f_regs[0])?f_regs[1]:f_regs[0];
                    if(region_ids.find(opp_reg_id)==region_ids.end() &&
                       done.find(opp_reg_id)==done.end()){
                        //we add it as a candidate since it is not
                        candidates.push_back(m_mesh->get<Region>(opp_reg_id));
                    }
                }
            }
        }
    }
    return region_ids;
}
/*---------------------------------------------------------------------------*/
std::vector<TCellID> PointConnectionBuilder::
getFaces(const gmds::Node& ANI, const gmds::Node& ANJ)
{
    std::vector<TCellID> fij, fi, fj;
    fi = ANI.getIDs<Face>();
    fj = ANJ.getIDs<Face>();
    for(auto i:fi){
        for(auto j:fj){
            if(i==j){
                fij.push_back(i);
            }
        }
    }
    return fij;
}
/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::isIn(const TCellID AI,
                           const std::vector<TCellID>& AV)
{
    return std::find(AV.begin(),AV.end(),AI)!=AV.end();
}

/*---------------------------------------------------------------------------*/
Face PointConnectionBuilder::getFace(const Face& AFrom,
                              const Region& AR,
                              const TCellID ANI,
                              const TCellID ANJ)
{
    std::vector<Face> faces = AR.get<Face>();
    for(const auto& f:faces){
        std::vector<TCellID> fn = f.getIDs<Node>();
        auto found_ni = false, found_nj=false;
        for(auto n_id:fn){
            if(n_id==ANI)
                found_ni=true;
            else if(n_id==ANJ)
                found_nj=true;
        }
        if(found_ni && found_nj && f.id()!=AFrom.id())
            return f;
    }
    throw GMDSException("PointConnectionBuilder::getFace(..) - No next face found");
}
/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::isIn(const int AFrom,
                           const int ATo,
                           const std::vector<OrientedEdge>& AEdgeSet,
                           OrientedEdge& AOutEdge) const
{
    OrientedEdge witness(AFrom,ATo);
    for(unsigned int i=0; i<AEdgeSet.size(); i++){
        if(witness  == AEdgeSet[i]){
            AOutEdge = AEdgeSet[i];
            return true;
        }
    }
    return false;
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
buildEdges(std::vector<std::vector<OrientedEdge> >& AInEdges,
           std::vector<std::vector<OrientedEdge> >& AOutEdges)
{
    int nb_points = m_pnt.size();

    AOutEdges.clear();
    AOutEdges.resize(nb_points);
    std::vector<int> on_curves;

    //======================================================================
    // STEP 1 - BUILD GLOBAL EDGES
    //======================================================================
    for(int i=0; i<nb_points; i++){

        if(m_classification[i]==ON_CURVE)
            on_curves.push_back(i);

        std::vector<OrientedEdge> orient_edges = AInEdges[i];

        for(auto e:orient_edges){
            int from = e.first;
            int to   = e.second;
            OrientedEdge e_out;
            //===============================================================
            //We check if ei has already been put in the out set of edges
            if(isIn(from,to,AOutEdges[from], e_out)){
                //YES, nothing to do so
                continue;
            }
            if(m_classification[from]>m_classification[to]){
                // Wse build the edge systematically
                // And the edge can be kept as the opposite
                AOutEdges[from].push_back(e);

                OrientedEdge e_inv(to,from);
                AOutEdges[to].push_back(e_inv);
            }
            else if(m_classification[from]==m_classification[to]){
                //We connect only if the connection exists in both sides
                OrientedEdge e_inv;
                if(isIn(to, from,
                        AInEdges[to], e_inv)){
                    //YES, We connect
                    AOutEdges[from].push_back(e);
                    AOutEdges[to  ].push_back(e_inv);
                }
            }
        }//for(unsigned int j=0; j<oriented_edges_i.size(); j++)

    }//for(int i=0; i<nb_points; i++)
    //======================================================================
    // STEP 2 - CLEAN SOME EDGES
    //======================================================================
    for(auto i:on_curves){
        //i= index of a point
        std::vector<OrientedEdge> edges = AOutEdges[i];
        std::vector<OrientedEdge> final_edges;
        for(auto j=0; j<edges.size(); j++){
            OrientedEdge ej = edges[j];

            if(m_classification[ej.second]==ON_CURVE ||
               m_classification[ej.second]==ON_VERTEX ){
                //Another edge connected to a surface can be a better choice
                bool found_best=false;
                for(auto k=0; k<edges.size(); k++){
                    OrientedEdge ek = edges[k];
                    if(k==j ||
                       m_classification[ek.second]==ON_CURVE ||
                       m_classification[ek.second]==ON_VERTEX){
                        continue;
                    }
                    //We compare length and direction of ej and ek
                    math::Vector3d vj(m_pnt[ej.first], m_pnt[ej.second]);
                    math::Vector3d vk(m_pnt[ek.first], m_pnt[ek.second]);
                    if(vk.norm2()<vj.norm2()){
                        vj.normalize();
                        vk.normalize();
                        if(vj.dot(vk)>0.8)
                            found_best=true;
                    }
                }//for(auto k=0; k<edges.size(); k++)
                if(!found_best){
                    final_edges.push_back(ej);
                }
            }
            else{
                final_edges.push_back(ej);
            }
        }//for(auto j=0; j<edges.size(); j++)

        //The cleaned set of edges is assigned to i
        AOutEdges[i]=final_edges;

    }
}
/*---------------------------------------------------------------------------*/
int PointConnectionBuilder::
filterPointsForBuildingOrientedEdge(const std::vector<int>& AID,
                                    const std::vector<double>& ADot,
                                    const std::vector<double>& ADist){
    //======================================================================
    //We get the best aligned candidate
    //======================================================================

    int    best_candidate_id =-1;
    double best_dot          = 0;
    for(auto i_c=0; i_c<AID.size();i_c++){
        if(ADot[i_c]>best_dot){
            best_dot=ADot[i_c];
            best_candidate_id=i_c;
        }
    }
    //Then we check if a closer one is in a reasonable angle tolerance
    int            ref_id     = best_candidate_id;
    double         ref_dist   = ADist[ref_id];

   // std::cout<<"\t best aligned: "<<best_candidate_id<<std::endl;

    int    new_best_id   = ref_id;
    double new_best_dist = ref_dist;

    for(unsigned int i_c=0; i_c<AID.size();i_c++){
        //Only closer point deserve to be checked
        if(i_c==ref_id || ADist[i_c]>=new_best_dist)
            continue;

        double dot_prod = abs(best_dot-ADot[i_c]);
        if((dot_prod<0.045) ||//5 degree of diff, we switch to this one,
           //we are quite close in angle
           (dot_prod<0.33 && ADist[i_c]<0.7*ref_dist) )// only switch if distance is really dimished
        {
            new_best_dist=ADist[i_c];
            new_best_id = i_c;
        }
    }//for(unsigned int i_c=0; i_c<final_candidates.size();i_c++)

    return new_best_id;
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
buildOrientedEdgesInVolume(const int            APntID,
                           const math::Vector3d AVec[][2],
                           math::Point          APnt[][2],
                           int                  AIndex[][2],
                           bool                 AFound[][2])
{
    int i = APntID;
    math::Point pi = m_pnt[i];

    // Means the point is in a stable area. So it is connected to 6 other
    // points if it is far way from a singularity area.
    double tol = 0.8;

    //======================================================================
    // STEP 1 - WE DETECT THE POSSIBLE POINTS
    //======================================================================
    for(int axis=0; axis<3; axis++){

        for(int dir=0; dir<2; dir++){

            math::Vector3d v_axis=AVec[axis][dir];
            v_axis.normalize(); //likely useless

            math::Point next_pi;
            if(!computeVolumePointFrom(i, v_axis,next_pi)){
                continue; //means we reached a FF-singular area
            }
            math::Vector3d v(pi,next_pi);
            std::vector<int> candidates = m_filter[i];
            math::Point pw(pi.X()+v.X(),
                           pi.Y()+v.Y(),
                           pi.Z()+v.Z());

            v.normalize();

            std::vector<int>            final_candidates;
            std::vector<double>         final_candidates_dot;
            std::vector<double>         final_candidates_dist;

            for(unsigned int j=0; j<candidates.size(); j++){

                math::Point pj = m_pnt[candidates[j]];
                math::Vector3d vij(pi,pj);
                vij.normalize();
                if(v.dot(vij)>tol) {
                    //pj is candidate so
                    final_candidates.push_back(candidates[j]);
                    final_candidates_dot.push_back(v.dot(vij));
                    final_candidates_dist.push_back(pi.distance(pj));

                }//if(v.dot(vij)>m_dot_tolerance)

            }//for(unsigned int j=0; j<candidates.size(); j++)

            //=========================================================
            // Now, we have to parse candidates to find the best one.
            // We start from the best aligned one, and we decrease while satisfying
            // the m_spacing criterion.
            // But if two candidates are only 10° different, we select via
            // distance.
            if(final_candidates.empty()){
                //We don't have found a candidate point
                continue;
            }

            int j = filterPointsForBuildingOrientedEdge(final_candidates,
                                                        final_candidates_dot,
                                                        final_candidates_dist);


            APnt  [axis][dir] = m_pnt[final_candidates[j]];
            AFound[axis][dir] = true;
            AIndex[axis][dir] = final_candidates[j];

        }//for(int dir=0;dir<2;dir++)

    }//for(int axis=0; axis<3; axis++)
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
buildOrientedEdgesOnSurface(const int            APntID,
                            const math::Vector3d AVec[][2],
                            math::Point          APnt[][2],
                            int                  AIndex[][2],
                            bool                 AFound[][2])
{
    int i = APntID;
    math::Point pi = m_pnt[i];

    // Means the point is in a stable area. So it is connected to 6 other
    // points if it is far way from a singularity area.
    double tol = 0.8;

    //We are on a surface, we must take care of points outside of the domain
    math::Vector3d n=m_normal[i];
    n.normalize();
    int surf_id = m_surface[i];

    // as the point is on the boundary, the classification process assigned
    // it to a starting triangle
    if(m_mesh_data[i].dim!=2)
        throw GMDSException("No face in buildOrientedEdgesOnSurface");
    gmds::Face t = m_mesh->get<Face>(m_mesh_data[i].id);

    //======================================================================
    // STEP 1 - WE DETECT THE POSSIBLE POINTS, WHICH ARE RESTRICTED TO BE ON
    //          THE BOUNDARY
    //======================================================================
    for(auto axis=0; axis<3; axis++){

        for(auto dir=0; dir<2; dir++){

            math::Vector3d v_axis=AVec[axis][dir];
            v_axis.normalize(); //likely useless
            math::Vector3d v;
            if((v_axis.dot(n))>tol){
                //We go out of the domain
                continue;
            }
            else if((v_axis.dot(n))<-tol){
                //we go inside the domain
                math::Point next_pi;
                if(computeVolumePointFrom(i, v_axis,next_pi))
                    v=math::Vector3d(pi,next_pi);
                else
                    continue; //it means we reached a FF-singular area
            }
            else{
                //we are on the surface
                math::Point next_pi;
                if(computeSurfacePointFrom(i, v_axis,next_pi))
                    v=math::Vector3d(pi,next_pi);
                else
                    continue; //it means we reached a FF-singular area

            }
            math::Point pw(pi.X()+v.X(),
                           pi.Y()+v.Y(),
                           pi.Z()+v.Z());
            std::vector<int> candidates = m_filter[i];

            v.normalize();

            std::vector<int>            final_candidates;
            std::vector<double>         final_candidates_dot;
            std::vector<double>         final_candidates_dist;


            for(unsigned int j=0; j<candidates.size(); j++){
                math::Point pj = m_pnt[candidates[j]];
                math::Vector3d vij(pi,pj);
                vij.normalize();

                //only volume and boundary well-aligned points are added
                if(v.dot(vij)>tol){

                    if(m_classification[candidates[j]]==ON_SURFACE){
                        //if the point is on a surface it must be
                        // on the same surface! No it can be a thin layer

                        if ( m_surface[candidates[j]]==surf_id){
                            //If they are on the same surface, we must discard
                            //configuration with cylinder shapes.
                            math::Vector3d nj=m_normal[candidates[j]];
                            if(nj.dot(n)<=0.0)
                                continue;
                        }
                    }

                    //pj is candidate so
                    final_candidates.push_back(candidates[j]);
                    final_candidates_dot.push_back(v.dot(vij));
                    final_candidates_dist.push_back(pi.distance(pj));
                }//if(v.dot(vij)>m_dot_tolerance)

            }//for(unsigned int j=0; j<candidates.size(); j++)
            //=========================================================
            // Now, we have to parse candidates to find the best one.
            // We start from the best aligned one, and we decrease while satisfying
            // the m_spacing criterion.
            // But if two candidates are only 10° different, we select via
            // distance.
            if(final_candidates.empty()){
                //We don't have found a candidate point
                continue;
            }

            int j = filterPointsForBuildingOrientedEdge(final_candidates,
                                                        final_candidates_dot,
                                                        final_candidates_dist);

            APnt  [axis][dir] = m_pnt[final_candidates[j]];
            AFound[axis][dir] = true;
            AIndex[axis][dir] = final_candidates[j];
        }//for(int dir=0;dir<2;dir++)

    }//for(int axis=0; axis<3; axis++)

}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
buildOrientedEdgesOnCurve(const int            APntID,
                          const math::Vector3d AVec[][2],
                          math::Point          APnt[][2],
                          int                  AIndex[][2],
                          bool                 AFound[][2])
{
    int i = APntID;
    math::Point pi = m_pnt[i];

//    std::cout<<"ORIENTED EDGES FROM "<<APntID<<" ON CURVE"<<std::endl;
    // Means the point is in a stable area. So it is connected to 6 other
    // points if it is far way from a singularity area.
    double tol = 0.9;

    //======================================================================
    // STEP 1 - WE DETECT THE POSSIBLE POINTS, WHICH ARE RESTRICTED TO BE ON
    //          THE BOUNDARY
    //======================================================================
    for(int axis=0; axis<3; axis++){

        for(int dir=0; dir<2; dir++){

            math::Vector3d v=AVec[axis][dir];
            v.normalize(); //likely useless

            std::vector<int>            candidates = m_filter[i];
            std::vector<int>            final_candidates;
            std::vector<double>         final_candidates_dot;
            std::vector<double>         final_candidates_dist;

            for(unsigned int j=0; j<candidates.size(); j++){
                math::Point pj = m_pnt[candidates[j]];
                math::Vector3d vij(pi,pj);
                vij.normalize();

                //only point well-aligned and on the boundary are taken
                //into account
                if(v.dot(vij)>tol &&
                   (m_classification[candidates[j]]==ON_CURVE  ||
                    m_classification[candidates[j]]==ON_VERTEX)){


                       //pj is candidate so
                       final_candidates.push_back(candidates[j]);
                       final_candidates_dot.push_back(v.dot(vij));
                       final_candidates_dist.push_back(pi.distance(pj));
                   }//if(v.dot(vij)>m_dot_tolerance)

            }//for(unsigned int j=0; j<candidates.size(); j++)
            //=========================================================
            // Now, we have to parse candidates to find the best one.
            // We start from the best aligned one, and we decrease while satisfying
            // the m_spacing criterion.
            // But if two candidates are only 10° different, we select via
            // distance.
            if(final_candidates.empty()){
                //We don't have found a candidate point
                continue;
            }

            int j = filterPointsForBuildingOrientedEdge(final_candidates,
                                                        final_candidates_dot,
                                                        final_candidates_dist);
            APnt  [axis][dir] = m_pnt[final_candidates[j]];
            AFound[axis][dir] = true;
            AIndex[axis][dir] = final_candidates[j];

        }//for(int dir=0;dir<2;dir++)

    }//for(int axis=0; axis<3; axis++)
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
buildOrientedEdges(std::vector<std::vector<OrientedEdge> >& AEdges)
{
    //One vector of edges per point. This vector can be empty so.
    AEdges.clear();
    AEdges.resize(m_pnt.size());

    //We go through all the points and we store for each of them the list of
    // oriented edges
    for(auto i=0; i<m_pnt.size(); i++){

        math::Point pi = m_pnt[i];

        if(m_type[i]==FRAME_SING){
            //only regular and param sing points are considered for
            //hex formation
            continue;
        }
        //==================================================================
        //We get the 6 vectors starting from pi following its chart
        math::Chart ci = m_chart[i];
        math::Vector3d ci_vectors[3][2] = {
            {ci.X(), -ci.X()},
            {ci.Y(), -ci.Y()},
            {ci.Z(), -ci.Z()}
        };

        //==================================================================
        // 6 corresponding points in directions x [0], y[1] and z[2]
        // the point in [0] is in positive direction
        // the point in [1] is in negative direction
        // For instance p[1][1] is the pnt starting from pi following -ci[1]
        math::Point p[3][2];
        int         p_index[3][2];
        bool        p_found[3][2]={{false,false},{false,false},{false,false}};


        if(m_classification[i]==IN_VOLUME){
            buildOrientedEdgesInVolume(i, ci_vectors, p, p_index, p_found);
        }
        else if(m_classification[i]==ON_SURFACE){
            buildOrientedEdgesOnSurface(i, ci_vectors, p, p_index, p_found);;
        }
        else if(m_classification[i]==ON_CURVE){
            buildOrientedEdgesOnCurve(i, ci_vectors, p, p_index, p_found);
        }

        //======================================================================
        // STEP 2 - ORIENTED EDGES ARE NOW CREATED
        //======================================================================
        std::vector<OrientedEdge> edges_i;
        for(auto axis=0; axis<3; axis++){
            for(auto dir=0; dir<2; dir++){
                if(p_found[axis][dir]==false)
                    continue;

                //We have an edge to create
                OrientedEdge e(i,                 // first point index
                               p_index[axis][dir],// second point index
                               axis,              // used chart axis in i
                               dir                // used chart direction in i
                               );

                edges_i.push_back(e);

            } //for(int dir=0;dir<2;dir++)

        } //for(auto axis=0; axis<3; axis++)
        AEdges[i]=edges_i;

    }
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
getEdges(std::vector<std::pair<int,int > >& AEdges) {
    for(auto edge_set:m_edges){
        for(auto e : edge_set){
            AEdges.push_back(std::make_pair(e.first,e.second));
        }
    }
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::getHexes(std::vector<std::vector<int> > &AHexes)
{
    Variable<int>* var_id = m_hexes.getVariable<int,GMDS_NODE>("HD_ID");

    for(auto h_id:m_hexes.regions()){
        std::vector<TCellID> node_ids = m_hexes.get<Region>(h_id).getIDs<Node>();
        std::vector<int> ids;
        ids.resize(8);
        for(auto i=0;i<8;i++) {
            ids[i]= var_id->value(node_ids[i]);
        }
        AHexes.push_back(ids);
    }
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
buildHexCorners(std::vector<std::vector<OrientedEdge> >& AEdges)
{
    for(unsigned int i=0; i<m_pnt.size(); i++){
        math::Point pi = m_pnt[i];
        std::vector<OrientedEdge> edges = AEdges[i];
        //We get the 6 vectors starting from pi following its chart
        math::Chart ci = m_chart[i];
        math::Vector3d ci_vectors[3][2] = {
            {ci.X(), -ci.X()},
            {ci.Y(), -ci.Y()},
            {ci.Z(), -ci.Z()}
        };

        int  p_index[3][2];
        bool p_found[3][2]= {
            {false,false},
            {false,false},
            {false,false}
        };

        bool found_undef_axis =false;
        int nb_undef_axi=0;
        for(unsigned int j=0; j<edges.size();j++){
            OrientedEdge ej = edges[j];
            if(ej.axis==-1 || ej.dir==-1){
                found_undef_axis=true;
                nb_undef_axi++;
            }
            else{
                p_found[ej.axis][ej.dir]=true;
                p_index[ej.axis][ej.dir]=ej.second;
            }
        }
        if(found_undef_axis){
            buildCornersAsSolidAngles(i, edges);
        }
        else{
            //now hex corner structure are created from 1 to 8 for pi
            //CORNER 1 - (X,Y,Z)
            if(p_found[0][0] && p_found[1][0] && p_found[2][0]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 2 - (X,Y,-Z)
            if(p_found[0][0] && p_found[1][0] && p_found[2][1]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][1], ci_vectors[2][1]);
            }
            //CORNER 3 - (X,Y,-Z)
            if(p_found[0][0] && p_found[1][1] && p_found[2][0]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 4 - (X,-Y,-Z)
            if(p_found[0][0] && p_found[1][1] && p_found[2][1]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][1], ci_vectors[2][1]);
            }
            //CORNER 5 - (-X,Y,Z)
            if(p_found[0][1] && p_found[1][0] && p_found[2][0]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 6 - (-X,Y,-Z)
            if(p_found[0][1] && p_found[1][0] && p_found[2][1]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][1], ci_vectors[2][1]);
            }
            //CORNER 7 - (-X,Y,-Z)
            if(p_found[0][1] && p_found[1][1] && p_found[2][0]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 8 - (-X,-Y,-Z)
            if(p_found[0][1] && p_found[1][1] && p_found[2][1]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][1], ci_vectors[2][1]);
            }
        }
    }//for(unsigned int i=0; i<m_pnt.size(); i++)...
}

/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
buildCornersAsSolidAngles(const int AOrigin,
                          std::vector<OrientedEdge> & AEdges)
{

    int origin_id = AOrigin;

    //=====================================================================
    // STEP 1 - Compute combination of 3 integers among [0;AEdges.size()]
    //=====================================================================
    std::vector<std::vector<int> > triplets;
    computeCombinations(AEdges.size(),3, triplets);

    //All edges come from the same points
    math::Point origin_pnt(m_pnt[origin_id]);
    //=====================================================================
    // STEP 2 - Filter valid combination, that is vector triplets such that
    // no other vector is inside their positive space quadrant
    //=====================================================================
    for(unsigned int i=0; i<triplets.size(); i++){
        std::vector<int> ti = triplets[i];
        math::Point p[3] = {
            m_pnt[AEdges[ti[0]].second],
            m_pnt[AEdges[ti[1]].second],
            m_pnt[AEdges[ti[2]].second]
        };

        // Points in p were only selected on a combinatorial way. If they
        // are geometrically lying on the same plan, they do not define a
        // valid corner

        math::Vector3d v12(origin_pnt,p[0]);
        math::Vector3d v13(origin_pnt,p[1]);
        math::Vector3d v14(origin_pnt,p[2]);


//      if (origin_pnt.areCoplanar(p[0], p[1], p[2])){
      if (fabs(v12.dot(v13.cross(v14)))<1e-11){
            continue;
        }

        bool valid = true;

        for(unsigned int j=0; j<AEdges.size() && valid; j++){

            int j_index = AEdges[j].second;

            //We don't check the vectors of the current basis
            if(j_index==AEdges[ti[0]].second ||
               j_index==AEdges[ti[1]].second ||
               j_index==AEdges[ti[2]].second)
                continue;

            math::Point pj =m_pnt[j_index];

            double coord[4];
            try {
                math::Point::computeBarycentric(origin_pnt,
                                                p[0], p[1], p[2], pj,
                                                coord[0], coord[1], coord[2],
                                                coord[3]);
            } catch (GMDSException& e) {
                valid=false;
            }
            if(coord[1]>0 && coord[2]>0 && coord[3]>0)
                valid=false;

            //Find where is vj comparing to v[0], v[1] and v[2]
        }
        if(valid){

            //So we can build a hex corner from this set of edges
            math::Vector3d v[3] = {
                math::Vector3d(origin_pnt,p[0]),
                math::Vector3d(origin_pnt,p[1]),
                math::Vector3d(origin_pnt,p[2])
            };
            addCorner(origin_id,
                      AEdges[ti[0]].second, v[0],
                      AEdges[ti[1]].second, v[1],
                      AEdges[ti[2]].second, v[2]);
        }
    }



}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
addCorner(const int AIndex,
          const int AIndex1, const math::Vector3d& AV1,
          const int AIndex2, const math::Vector3d& AV2,
          const int AIndex3, const math::Vector3d& AV3)
{
    HexCorner corner;
    corner.p = AIndex;
    corner.index =m_hc_mapping[AIndex].size();
    corner.free=true;

    corner.adj[0] = AIndex1;
    corner.vec[0] = AV1;

    corner.adj[1] = AIndex2;
    corner.vec[1] = AV2;

    corner.adj[2] = AIndex3;
    corner.vec[2] = AV3;

    m_hc_mapping[AIndex].push_back(corner);
}
/*---------------------------------------------------------------------------*/
int PointConnectionBuilder::findFreeCorners(const std::vector<HexCorner>& AIn,
                                     const int AFrom,
                                     std::vector<HexCorner>& AOut)
{
    AOut.clear();
    for(auto const& c:AIn){
        if(!c.free) {
            //this corner is already used for a building another hex
            continue;
        }
        //c must point to the point AFrom
        bool found = false;
        for(int i=0; i<3; i++){
            if(c.adj[i]==AFrom){
                found=true;
            }
        }
        if(found)
            AOut.push_back(c);
    }
    return AOut.size();
}
/*---------------------------------------------------------------------------*/
std::vector<int> PointConnectionBuilder::findCommonPoints(const int AFrom,
                                                   const int AI,
                                                   const int AJ)
{
    std::vector<HexCorner> corners_i, corners_j;
    findFreeCorners(m_hc_mapping[AI], AFrom, corners_i);
    findFreeCorners(m_hc_mapping[AJ], AFrom, corners_j);

    std::vector<int> common;

    for(auto ci:corners_i){
        for(auto cj:corners_j){

            for(int i=0; i<3; i++) {
                if(ci.adj[i]==AFrom)
                    continue;
                for(int j=0; j<3; j++) {
                    if(cj.adj[j]==AFrom)
                        continue;

                    if(ci.adj[i]==cj.adj[j])
                        common.push_back(ci.adj[i]);
                }
            }
        }
    }
    return common;
}

/*---------------------------------------------------------------------------*/
std::vector<PointConnectionBuilder::HexCorner>
PointConnectionBuilder::findFreeCorners(const std::vector<int>& AOrigin,
                                 const int AI,
                                 const int AJ)
{
    std::vector<HexCorner>   out;

    for(auto& org: AOrigin){
        std::vector<HexCorner> corners_i;
        findFreeCorners(m_hc_mapping[org], AI, corners_i);
        for(auto ci:corners_i) {

            bool found = false;

            for(int i=0; i<3; i++) {
                if(ci.adj[i]==AJ)
                    found =true;
            }
            if(found)
                out.push_back(ci);
        }
    }
    return out;
}

/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::getCorner(const int AOrigin,
                                const int AI,
                                const int AJ,
                                const int AK,
                                HexCorner& AOut)
{
    for(auto& c: m_hc_mapping[AOrigin]){
        if((c.adj[0]==AI && c.adj[1]==AJ && c.adj[2]==AK) ||
           (c.adj[0]==AI && c.adj[2]==AJ && c.adj[1]==AK) ||
           (c.adj[1]==AI && c.adj[0]==AJ && c.adj[2]==AK) ||
           (c.adj[1]==AI && c.adj[2]==AJ && c.adj[0]==AK) ||
           (c.adj[2]==AI && c.adj[0]==AJ && c.adj[1]==AK) ||
           (c.adj[2]==AI && c.adj[1]==AJ && c.adj[0]==AK) ){
            AOut = c;
            return true;
        }
    }
    return false;
}

/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::isCorner(const HexCorner& AC,
                               const int AOrigin,
                               const int AI,
                               const int AJ,
                               const int AK )
{
    if(AC.p!=AOrigin)
        return false;

    if((AC.adj[0]==AI && AC.adj[1]==AJ && AC.adj[2]==AK) ||
       (AC.adj[0]==AI && AC.adj[2]==AJ && AC.adj[1]==AK) ||
       (AC.adj[1]==AI && AC.adj[0]==AJ && AC.adj[2]==AK) ||
       (AC.adj[1]==AI && AC.adj[2]==AJ && AC.adj[0]==AK) ||
       (AC.adj[2]==AI && AC.adj[0]==AJ && AC.adj[1]==AK) ||
       (AC.adj[2]==AI && AC.adj[1]==AJ && AC.adj[0]==AK) ){
        return true;
    }
    return false;
}
/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::findCommmonLastCorner(const HexCorner& ACorner1,
                                                   const HexCorner& ACorner2,
                                                   const HexCorner& ACorner3,
                                                   HexCorner&       ACornerOut)
{
    //We look for the point shared by ACorner1, ACorner2 and ACorner3
    int common_pnt_12[2]={-1,-1};
    int i_12=0;
    for(int i1=0; i1<3; i1++){
        int adj1= ACorner1.adj[i1];
        for(int i2=0; i2<3; i2++){
            if(adj1 == ACorner2.adj[i2]){
                common_pnt_12[i_12++]=adj1;
            }
        }
    }//for(int i1=0; i1<3; i1++)

    int common_pnt=-1;
    for(int i3=0; i3<3; i3++){
        int adj3= ACorner3.adj[i3];
        for(int i12=0; i12<2; i12++){
            if(adj3 == common_pnt_12[i12]){
                common_pnt =adj3;
            }
        }
    }//for(int i1=0; i1<3; i1++){

    if(common_pnt==-1){
        return false;
    }

    // Among all corner, we need the one pointing to ACorner1, ACorner2 and
    // ACorner3
    std::vector<HexCorner> candidates = m_hc_mapping[common_pnt];
    for(unsigned int i = 0; i<candidates.size(); i++){
        HexCorner ci= candidates[i];
        if(!ci.free)
            continue;
        int nb_same=0;
        for(int j=0; j<3; j++){
            if(ci.adj[j]==ACorner1.p ||
               ci.adj[j]==ACorner2.p ||
               ci.adj[j]==ACorner3.p )
                nb_same++;

        }
        if(nb_same==3){
            ACornerOut = ci;
            return true;
        }
    }
    return false;
}
/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::buildHexahedral()
{
    Variable<int>* var_id = m_hexes.newVariable<int,GMDS_NODE>("HD_ID");
    Variable<int>* var_type = m_hexes.newVariable<int,GMDS_REGION>("TYPE");
    Variable<int>* var_cl = m_hexes.newVariable<int,GMDS_NODE>("classification");
    Variable<int>* var_surf = m_hexes.newVariable<int,GMDS_NODE>("surface_id");
    Variable<int>* var_curv = m_hexes.newVariable<int,GMDS_NODE>("curve_id");

    /* We go trought all the points and we build hexahedral elements
     * for each associated free corner */
    for(unsigned int i=0; i<m_pnt.size(); i++){

        math::Point pi = m_pnt[i];
        std::vector<HexCorner> corners_i = m_hc_mapping[i];
        //Now we go throught the associated corners

        for(unsigned int j=0; j<corners_i.size();j++){
            HexCorner cj = corners_i[j];

            if(!cj.free) //this corner is already used for a hex
                continue;
            //=====================================================
            // STEP 1 - We look for 3 adjacent compatible corners
            // to form a hexahedral element
            //=====================================================
            std::vector<int> common_01 =findCommonPoints(cj.p,cj.adj[0],cj.adj[1]);
            std::vector<int> common_02 =findCommonPoints(cj.p,cj.adj[0],cj.adj[2]);
            std::vector<int> common_12 =findCommonPoints(cj.p,cj.adj[1],cj.adj[2]);

            if(common_01.empty() || common_02.empty() || common_12.empty())
                continue;

            std::vector<HexCorner> corners_01 = findFreeCorners(common_01, cj.adj[0],cj.adj[1]);
            std::vector<HexCorner> corners_02 = findFreeCorners(common_02, cj.adj[0],cj.adj[2]);
            std::vector<HexCorner> corners_12 = findFreeCorners(common_12, cj.adj[1],cj.adj[2]);

            if(corners_01.empty()|| corners_02.empty()|| corners_12.empty())
                continue;

            bool found_last_corner = false;

            HexCorner final_01, final_02, final_12, final_012;
            for(auto c_01:corners_01){
                for(auto c_02:corners_02){
                    for(auto c_12:corners_12){
                        if(c_01.p==c_02.p || c_01.p==c_12.p ||c_02.p==c_12.p)
                            continue;

                        HexCorner c_012;
                        if(findCommmonLastCorner(c_01,c_02,c_12,c_012)){
                            found_last_corner=true;
                            final_01=c_01;
                            final_02=c_02;
                            final_12=c_12;
                            final_012=c_012;
                        }
                    }
                }
            }
            if(!found_last_corner)
                continue;

            HexCorner final_0, final_1, final_2;

            if(!getCorner(cj.adj[0], cj.p,final_01.p,final_02.p, final_0) ||
               !getCorner(cj.adj[1], cj.p,final_01.p,final_12.p, final_1) ||
               !getCorner(cj.adj[2], cj.p,final_02.p,final_12.p, final_2))
                continue;

            //Now we have the 8 corners, we can build the hexahedron
            int idx[8] = {static_cast<int>(i),
                final_0.p,
                final_01.p,
                final_1.p,
                final_2.p,
                final_02.p,
                final_012.p,
                final_12.p
            };

            Node n[8];
            for(int i_n=0; i_n<8; i_n++){
                int index=idx[i_n];
                std::map<int, Node>::iterator it=m_node_mapping.find(index);
                if(it==m_node_mapping.end()){
                    //New node to add
                    n[i_n]=m_hexes.newNode(m_pnt[index]);
                    (*var_id)[n[i_n].id()] = index;
                    (*var_cl)[n[i_n].id()] = m_classification[index];
                    (*var_surf)[n[i_n].id()] = m_surface[index];
                    (*var_curv)[n[i_n].id()] = m_curve[index];
                    m_node_mapping[index]=n[i_n];
                    m_used[index] = true;
                }
                else {
                    //existing node
                    n[i_n] = it->second;
                }
            }

            Region r = m_hexes.newHex(n[0], n[1], n[2], n[3],
                                       n[4], n[5], n[6], n[7]);
            (*var_type)[r.id()]=GMDS_HEX;
            //Eventually, used corners are no more free!
            m_hc_mapping[cj.p][cj.index].free=false;
            m_hc_mapping[final_0.p][final_0.index].free=false;
            m_hc_mapping[final_1.p][final_1.index].free=false;
            m_hc_mapping[final_2.p][final_2.index].free=false;
            m_hc_mapping[final_12.p][final_12.index].free=false;
            m_hc_mapping[final_01.p][final_01.index].free=false;
            m_hc_mapping[final_02.p][final_02.index].free=false;
            m_hc_mapping[final_012.p][final_012.index].free=false;

        }//for(unsigned j=0; j<corners_i.size();j++)

    }//for(unsigned int i=0; i<m_pnt.size(); i++
}
/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::computeVolumePointFrom(const int APntIndex,
                                             const math::Vector3d& AV,
                                             math::Point& APnt)
{
    math::Point p = m_pnt[APntIndex];
    double max_dist = m_spacing;
    Region t;
    APnt.X() = p.X()+AV.X();
    APnt.Y() = p.Y()+AV.Y();
    APnt.Z() = p.Z()+AV.Z();
    return true;

}
/*---------------------------------------------------------------------------*/
bool PointConnectionBuilder::computeSurfacePointFrom(const int APntIndex,
                                              const math::Vector3d& AV,
                                              math::Point& APnt)
{
    math::Point p = m_pnt[APntIndex];
    APnt.X() = p.X()+AV.X();
    APnt.Y() = p.Y()+AV.Y();
    APnt.Z() = p.Z()+AV.Z();
    return true;
}

/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::writeInput()
{
    MeshModel model(DIM3 | R  | N | R2N );
    Mesh mesh_pnts  (model);
    Mesh mesh_charts(model);
    Variable<int>* var_pnts   = mesh_pnts.newVariable<int,GMDS_REGION>("type");
    Variable<int>* var_charts_type = mesh_charts.newVariable<int,GMDS_REGION>("type");
    Variable<int>* var_charts_idx = mesh_charts.newVariable<int,GMDS_REGION>("index");
    Variable<int>* var_charts_cl = mesh_charts.newVariable<int,GMDS_REGION>("classification");
    Variable<int>* var_curv = mesh_charts.newVariable<int,GMDS_REGION>("curv_number");
    Variable<int>* var_surf = mesh_charts.newVariable<int,GMDS_REGION>("surf_number");
    Variable<int>* var_used = mesh_charts.newVariable<int,GMDS_REGION>("used");


    double cube_size = m_spacing/5.0;

    for(unsigned int i=0; i<m_pnt.size();i++){
        math::Point pi = m_pnt[i];
        math::Chart ci = m_chart[i];
        int ti = m_type[i];

        //POINT OUTPUT
        Node   ni  = mesh_pnts.newNode(pi);
        Region ri = mesh_pnts.newTet(ni,ni,ni,ni);
        (*var_pnts)[ri.id()]= ti;

        //CHART OUTPUT
        math::Vector vx = ci.X();
        math::Vector vy = ci.Y();
        math::Vector vz = ci.Z();
        math::Point center = pi;
        math::Point p1 = center + (vx + vy - vz)*cube_size;
        Node n1 = mesh_charts.newNode(p1);
        math::Point p2 = center + (vx - vy - vz)*cube_size;
        Node n2 = mesh_charts.newNode(p2);
        math::Point p3 = center + (vx + vy + vz).opp()*cube_size;
        Node n3 = mesh_charts.newNode(p3);
        math::Point p4 = center + (vy - vx - vz)*cube_size;
        Node n4 = mesh_charts.newNode(p4);

        math::Point p5 = center + (vx + vy + vz)*cube_size;
        Node n5 = mesh_charts.newNode(p5);
        math::Point p6 = center + (vx - vy + vz)*cube_size;
        Node n6 = mesh_charts.newNode(p6);
        math::Point p7 = center + (vx + vy - vz).opp()*cube_size;
        Node n7 = mesh_charts.newNode(p7);
        math::Point p8 = center + (vy - vx + vz)*cube_size;
        Node n8 = mesh_charts.newNode(p8);
        Region r = mesh_charts.newHex(n1, n2, n3, n4, n5, n6, n7, n8);
        (*var_charts_type)[r.id()]= ti;
        (*var_charts_idx)[r.id()] = i;
        (*var_charts_cl)[r.id()] = m_classification[i];
        (*var_curv)[r.id()] = m_curve[i];
        (*var_surf)[r.id()] = m_surface[i];
        (*var_used)[r.id()] = m_used[i];
    }
    static int nb_file=0;
    std::string file_pnts, file_charts;
    file_pnts=m_output_dir+"/PCB_INPUT_PNTS_"+to_string(nb_file)+".vtk";
    file_charts=m_output_dir+"/PCB_INPUT_CHARTS_"+to_string(nb_file)+".vtk";
    nb_file++;

    IGMeshIOService ioService(&mesh_pnts);
    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::R);
    vtkWriter.write(file_pnts);

    IGMeshIOService ioService2(&mesh_charts);
    VTKWriter vtkWriter2(&ioService2);
    vtkWriter2.setCellOptions(gmds::N|gmds::R);
    vtkWriter2.setDataOptions(gmds::N|gmds::R);
    vtkWriter2.write(file_charts);

}

/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::
writeEdges(std::vector<std::vector<OrientedEdge> >& AEdges,
           const std::string& AFileName)
{
    MeshModel model(DIM3 | F  | N | F2N );
    Mesh mesh_edges  (model);
    for(unsigned int i=0; i<m_pnt.size();i++){

        std::vector<OrientedEdge> edges_i = AEdges[i];

        for(unsigned int j=0; j<edges_i.size(); j++){
            math::Point p0 = m_pnt[edges_i[j].first];
            math::Point p1 = m_pnt[edges_i[j].second];
            Node   n0 = mesh_edges.newNode(p0);
            Node   n1 = mesh_edges.newNode(p1);
            mesh_edges.newTriangle(n0,n1,n1);

        }
    }
    IGMeshIOService ioService(&mesh_edges);
    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write(AFileName);
}


/*---------------------------------------------------------------------------*/
void PointConnectionBuilder::writeHexes() {
    IGMeshIOService ioService(&m_hexes);
    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(N | F | R);
    vtkWriter.setDataOptions(N | F | R);
    std::string file_name = m_output_dir + "/EXTRACTED_HEX.vtk";
    vtkWriter.write(file_name);
}
/*----------------------------------------------------------------------------*/
math::Vector3d PointConnectionBuilder::
getOutputNormal(Face& AFace, Region& ARegion)
{
    std::vector<Node> region_nodes = ARegion.get<Node>();
    std::vector<Node> face_nodes = AFace.get<Node>();

    math::Point  reg_center  = ARegion.center();
    math::Point  face_center = AFace.center();
    math::Vector face_normal = AFace.normal();
    math::Vector v_ref(face_center, reg_center);

    if(face_normal.dot(v_ref)>0)
        return math::Vector3d(-face_normal.X(), -face_normal.Y(), -face_normal.Z());

    return math::Vector3d(face_normal.X(), face_normal.Y(), face_normal.Z());

}
