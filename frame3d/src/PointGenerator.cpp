/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Triangle.h>
#include <gmds/utils/Log.h>

/*---------------------------------------------------------------------------*/
// FRAME File Headers
#include <gmds/frame3d/PointGenerator.h>
/*----------------------------------------------------------------------------*/
// OpenNL File Headers
#include "OpenNL_psm.h"
/*----------------------------------------------------------------------------*/
// ANN File Headers
#include "ANN/ANN.h"
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <set>
#include <algorithm>
#include <queue>
#include <gmds/math/Tetrahedron.h>
#include <gmds/io/IGMeshIOService.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
struct SPointData{
    math::Point point;
    int classification;
    int type;
    int index;
} ;
/*---------------------------------------------------------------------------*/
// \brief Structure to sort points stored in a std::vector
struct SPointDataStrictWeakOrder{
    bool operator()(SPointData AP1, SPointData AP2) {
        return AP1.classification < AP2.classification;
    }
} pointDataStrictWeakOrder;
/*---------------------------------------------------------------------------*/
// \brief Structure to sort points stored in a std::vector
struct SPointStrictWeakOrder{
    bool operator()(math::Point AP1, math::Point AP2) {
        return AP1 < AP2;
    }
};
/*---------------------------------------------------------------------------*/
PointGenerator::
PointGenerator(Mesh* AMesh,
               const ParamsGlobal& AParamGl,
               std::map<gmds::TCellID, gmds::math::Vector3d>& ANormal,
               const ParamsMark& AMarks,
               const double ASpacing,
               const double ACurl):
m_mesh(AMesh),
m_param_gl(AParamGl),
m_normal(ANormal),
m_bm(AMarks),
m_spacing(ASpacing),
m_curl(0.35),
m_nb_unknowns(0)
{
    m_rotation_field = m_mesh->getVariable<math::AxisAngleRotation, GMDS_NODE>("rotation_field");

    m_surface_face_color = m_mesh->getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
    m_curve_edge_color = m_mesh->getVariable<int,GMDS_EDGE>("BND_CURVE_COLOR");

    m_surface_node_color= m_mesh->getVariable<int,GMDS_NODE>("BND_SURFACE_COLOR");
    m_curve_node_color= m_mesh->getVariable<int, GMDS_NODE>("BND_CURVE_COLOR");

}
/*---------------------------------------------------------------------------*/
math::Chart::Mapping PointGenerator::getRij(const TCellID AFrom,
                                            const TCellID ATo) const
{
    math::AxisAngleRotation rotation_from = (*m_rotation_field)[AFrom];
    math::AxisAngleRotation rotation_to   = (*m_rotation_field)[ATo];

    math::Chart chart_from  = rotation_from.toChart();
    math::Chart chart_to    = rotation_to.toChart();

    return math::Chart::Mapping(chart_from,chart_to);

}
/*---------------------------------------------------------------------------*/
math::Chart::Mapping PointGenerator::getRij(const Node& AFrom,
                                            const Node& ATo) const
{
    return getRij(AFrom.id(),ATo.id());
}
/*---------------------------------------------------------------------------*/
math::Vector3d PointGenerator::computeGij(OrientedEdge& AE)  {
    Node ni = AE.first;
    Node nj = AE.second;
    
    TCellID i = ni.id();
    TCellID j = nj.id();
    
    math::Vector3d ref_XYZ[3] = {
        math::Vector3d(1., 0., 0.),
        math::Vector3d(0., 1., 0.),
        math::Vector3d(0., 0., 1.)
    };
    math::AxisAngleRotation rot_i = (*m_rotation_field)[i];
    math::AxisAngleRotation rot_j = (*m_rotation_field)[j];
    
    math::Chart::Mapping r_ij = getRij(i,j);
    
    //vector from pi to pj
    math::Vector3d xij(ni.point(), nj.point());
    
    math::Vector3d gij_component[3];
    for (int d = 0; d<3; d++) {
        math::Vector3d tmp = (rot_j*(r_ij.inverse()*ref_XYZ[d]) +
                              rot_i*ref_XYZ[d]);
        gij_component[d] = tmp;
    }
    
    math::Vector3d gij;
    for (int i_dim = 0; i_dim<3; i_dim++) {
        gij[i_dim] = 0.5*xij.dot(gij_component[i_dim]);
    }
    
    return 2*M_PI*gij/m_spacing;
    //2PI since used for cos and sin
    //1/m_spacing since Bi and Bj should have a m_spacing length and so Bi-1
    //and Bj-1 have the inverse
}
/*---------------------------------------------------------------------------*/
math::Vector3d PointGenerator::computeCij(OrientedEdge& AE)
{
    if(m_curl == 0.0) {
        return math::Vector3d(0, 0, 0);
    }
    
    Edge e = AE.edge;
    TCellID e_id = e.id();
    
    math::Chart::Mapping r_ij = getRij(AE.first,AE.second);
    
    math::Vector3d correction(0, 0, 0);
    math::Vector3d edge_correction = m_corr[e_id];
    for (int d = 0;d<3;d++){
        correction[d] = edge_correction[d];
    }
    if (!AE.isWellOriented()){
        correction = -(r_ij*correction);
    }
    
    return correction;
}

/*---------------------------------------------------------------------------*/
std::vector<PointGenerator::OrientedEdge> PointGenerator::
getOrientedEdge(const gmds::Face& AF){
    std::vector<Node> n = AF.get<Node>();
    std::vector<Edge> e = AF.get<Edge>();
    
    std::vector<OrientedEdge> oriented_edges;
    oriented_edges.resize(3);
    
    for(int i=0; i<3; i++){
        int i0 = n[i].id();
        int i1 = n[(i+1)%3].id();
        
        bool found_edge = false;
        Edge ei;
        for(int j=0; j<3 && !found_edge; j++){
            Edge ej = e[j];
            
            std::vector<TCellID> ej_nodes = ej.getIDs<Node>();
            
            if((ej_nodes[0]==i0 && ej_nodes[1]==i1) ||
               (ej_nodes[0]==i1 && ej_nodes[1]==i0)){
                ei=ej;
                found_edge=true;
            }
        }
        OrientedEdge oei(ei,n[i],n[(i+1)%3]);
        oriented_edges[i]=oei;
    }
    return oriented_edges;
}

/*---------------------------------------------------------------------------*/
void PointGenerator::computeCurlCorrection()
{
    int nb_variables = 3*m_mesh->getNbEdges();

    if(m_curl == 0.0) {
        Log::mng()<< "No curl-correction (provided target is null)\n";
        Log::mng().flush();
        return;
    }
    
    //======================================================================
    // 1 - The correction term of each edge is set to (0,0,0)
    //======================================================================
    for(auto e_id : m_mesh->edges()){
        m_corr[e_id] = math::Vector3d(0, 0, 0);
    }

    //======================================================================
    // 2 - System building
    //======================================================================
    Log::mng()<< "Curl correction - Matrix Assembly\n";
    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, NLint(nb_variables));
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    //======================================================================
    // 2-a - Constraint on boundary edge are set
    //======================================================================
    Log::mng()<<"Curl correction  - boundary constraints on edges\n";
    Log::mng().flush();
    
    nlBegin(NL_SYSTEM);
    for(auto e_id : m_mesh->edges()){

        Edge e = m_mesh->get<Edge>(e_id);
        std::vector<TCellID> e_node_ids = e.getIDs<Node>();
        TCellID i0 =e_node_ids[0];
        TCellID i1 =e_node_ids[1];
        
        //We get the mapping transforming the chart in i0 to those in i1
        math::Chart::Mapping r_01  = getRij(i0,i1);
        
        //We get the constraint assigned to each end nodes of e
        math::Vector3d bnd_constr[2] = {
            m_bnd_constraint[i0],
            m_bnd_constraint[i1]
        };
        //Constraint in i1 are "seen" from the basis in i0
        bnd_constr[1] = r_01*bnd_constr[1];
        for (int i = 0; i < 3; i++) {
            bnd_constr[1][i] = abs(bnd_constr[1][i]);
        }
        
        for (int i = 0; i < 3; i++) {
            if (bnd_constr[0][i] > 0 && bnd_constr[1][i] > 0 ) {
                nlSetVariable (3 * m_edge_id[e.id()] + i, 0);
                nlLockVariable(3 * m_edge_id[e.id()] + i);
            }
        }

    }//for(it_e=m_mesh->edges_begin(); !it_e.isDone() ; it_e.next())
    
    nlBegin(NL_MATRIX);
    
    //======================================================================
    // 2-b - 1-form closure of the faces (Disc. Ext. Calculus)
    //======================================================================
    Log::mng() << "Curl correction - 1-form closure for faces\n";
    for (auto f_id:m_mesh->faces()){
        Face f = m_mesh->get<Face>(f_id);

        // All the computation are performed in thefirst vertex's basis
        std::vector<TCellID> n_id = f.getIDs<Node>();
        
        std::vector<OrientedEdge> orient_e = getOrientedEdge(f);
        
         // init chg basis for each edge following the pattern
        math::Chart::Mapping r_ij[3];
        for (int ie = 0; ie<3; ie++)  {
            if (orient_e[ie].isWellOriented()){
                r_ij[ie] = getRij(n_id[0], n_id[ie]);
            }
            else{
                r_ij[ie] = getRij(n_id[0], n_id[(ie + 1) % 3]);
                r_ij[ie].m_dir *= -1;
            }
        }
        
        // Righ-hand side term in the AX=b system we build in
        math::Vector3d b(0, 0, 0);
        // add each edge's contribution
        for (int i = 0; i < 3; i++) {
            OrientedEdge ei = orient_e[i];
            if(!ei.isWellOriented()) {
                ei =OrientedEdge(ei.edge,ei.second,ei.first);
            }
            b -= r_ij[i] * computeGij(ei);
        }
        
        
        //Add a line into the system
        for (int i_dim = 0; i_dim<3; i_dim++){
            nlBegin(NL_ROW);
            for (int i_edge = 0; i_edge < 3; i_edge++) {
                nlCoefficient(3*m_edge_id[orient_e[i_edge].edge.id()] +
                              r_ij[i_edge].m_map[i_dim],
                              r_ij[i_edge].m_dir[i_dim]) ;
            }
            nlRightHandSide(b[i_dim]);
            nlEnd(NL_ROW);
        }
        
    }//for (;!it_f.isDone();it_f.next())

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    
    
    //======================================================================
    // 3 - System resolution
    //======================================================================
    solveSystem();
    
    //======================================================================
    // 4 - Solution is put in the m_corr attributed
    //======================================================================
    std::map<TCellID, gmds::math::Vector3d >::iterator corr_it;
    for(corr_it=m_corr.begin(); corr_it!=m_corr.end();corr_it++){
        int lid = m_edge_id[corr_it->first];
        math::Vector3d cor_vec(nlGetVariable(3*lid  ),
                               nlGetVariable(3*lid+1),
                               nlGetVariable(3*lid+2));
        
        corr_it->second=cor_vec;
    }
    
    nlDeleteContext(nlGetCurrent());

    //======================================================================
    // 5 - Obtained solution is scale in order to avoid too large
    //     distortions
    //======================================================================
    Log::mng()<< "Curl correction - Scaling result\n";
    Log::mng().flush();
    for(auto e_id : m_mesh->edges()){
        Edge e = m_mesh->get<Edge>(e_id);
        OrientedEdge oe(e);
        
        math::Vector3d g_ij = computeGij(oe);
        math::Vector3d c_ij = computeCij(oe);
        double scale =1.;
        if (c_ij.norm()>m_curl * g_ij.norm())  {
            scale = m_curl * g_ij.norm() / c_ij.norm();
        }
        m_corr[e.id()] = scale * m_corr[e.id()];
    }
}

/*---------------------------------------------------------------------------*/
void PointGenerator::execute()
{

    buildLocalIDs();

    //======================================================================
    // PART 0 - Initialization (boundary constraints + curl correction)
    //======================================================================
    setBoundaryConstraint();

    computeCurlCorrection();

    //======================================================================
    // PART 1 - PARAMETRIZATION BUILDING
    //======================================================================
    // NUMBER OF UNKNOWNS
    // Unknowns of the system are ui=Mi(Xi), Rij and tij
    // Rij is solved independly considering optimal basis change along each
    // mesh edges (as done with quaternion comparison in previous works)
    // tij is a 3-dim vector corresponding to a translation
    // ui is a "mapping function" in Xi, so only the result of this mapping
    // functin in Xi, that is a 3D point in the parametric space
    // We have so ui to find in each point in a (cos, sin) representation
    // so 3*2 = 6 variables per node
    //======================================================================
    m_nb_unknowns = 6*m_mesh->getNbNodes();

    initSystem();

    buildSystem();

    solveSystem();

    getUiSolution();

    //======================================================================
    // PART 2 - POINT GENERATION
    //======================================================================
    std::vector<math::Point> points;
    extractPoints();


    cleanSystem();
    if(m_param_gl.with_debug_files){
        writeOutput();
    }
}
/*---------------------------------------------------------------------------*/
void PointGenerator::solveSystem(){
    
    NLint NNZ;
    nlGetIntegerv(NL_NNZ, &NNZ);
    
    Log::mng()<< "Solve with nb unknowns =" << m_nb_unknowns
    << " NNZ = " << NNZ
    << "...\n";
    Log::mng().flush();
    
    nlSolve();
    
    NLint nb_iter;
    NLdouble Gflops;
    NLdouble elapsed;
    nlGetIntegerv(NL_USED_ITERATIONS, &nb_iter);
    nlGetDoublev(NL_GFLOPS, &Gflops);
    nlGetDoublev(NL_ELAPSED_TIME, &elapsed);
    
    Log::mng()<<"Elapsed time:" << elapsed <<" s";
    Log::mng()<<" (in " << nb_iter << " iterations)\n";
    Log::mng().flush();
    
}
/*---------------------------------------------------------------------------*/
void PointGenerator::getUiSolution(){
  
    for(auto g_id:m_mesh->nodes()) {
         int l_id = m_id[g_id];
        math::Vector3d ui(0, 0, 0);
        double c = (.5 / M_PI);
        
        for(int k=0; k<3; k++){
            ui[k] = c* atan2(nlGetVariable(6 * l_id + 2 * k + 1),
                             nlGetVariable(6 * l_id + 2 * k ));
        }
        m_Ui[g_id]=ui;
    }
}
/*---------------------------------------------------------------------------*/
math::Vector3d PointGenerator::computeGijWithCurl(Edge& AEdge,
                                                  Node& AFrom,
                                                  Node& ATo)
{
    OrientedEdge oe(AEdge,AFrom,ATo);
    return computeGij(oe) + computeCij(oe);
}
/*---------------------------------------------------------------------------*/
void PointGenerator::buildLocalIDs()
{
    //======================================================================
    // For this first version all the mesh nodes are considered,
    // it should change in the future
    //======================================================================
    // Node numbering
    int local_index=0;
    for(auto n_id:m_mesh->nodes()) {
        m_id[n_id]=local_index++;
    }
    //======================================================================
    // Edge numbering
    local_index=0;
    for(auto e_id:m_mesh->edges()) {
        m_edge_id[e_id]=local_index++;
    }
}

/*---------------------------------------------------------------------------*/
void PointGenerator::setBoundaryConstraint()
{
    //======================================================================
    // Only boundary nodes are constrained.
    // We initialize their constraint
    //======================================================================
    for(auto n_id : m_mesh->nodes()) {
        if(m_mesh->isMarked<Node>(n_id,m_bm.mark_node_on_pnt)){
            m_bnd_constraint[n_id]=math::Vector3d(1,1,1);
        }
        else{
            m_bnd_constraint[n_id] = math::Vector3d(0,0,0);
        }
    }
    //======================================================================
    // We go through each boundary face to constrained its
    // incicent nodes (so nodes on surf, curves and points)
    //======================================================================
    for(auto f_id : m_mesh->faces()) {
        
        Face f = m_mesh->get<Face>(f_id);
        //if f is not on the surface, we pass it
        if(!m_mesh->isMarked(f,m_bm.mark_face_on_surf)){
            continue;
        }
        
        std::vector<Node> f_nodes=f.get<Node>();
        
        math::Vector nf = f.normal();
        math::Vector3d f_normal(nf.X(), nf.Y(), nf.Z());
        
        //for each node of a boudary face, we look for its chart alignment
        //with the normal to the surface defined by f
        for(auto ni : f_nodes){
            TCellID ni_id = ni.id();
            
            math::AxisAngleRotation ri = (*m_rotation_field)[ni_id];
            math::Chart ci = ri.toChart();
            
            double best_fit_value = 0;
            int    best_fid_id    = -1;
            for(int i_axis=0; i_axis<3; i_axis++){
                double  value = fabs(ci[i_axis].dot(f_normal));
                if(value>0.9){
                    best_fit_value = value;
                    best_fid_id    = i_axis;
                }
            }
            if(best_fid_id!=-1){
                m_bnd_constraint[ni_id][best_fid_id]=1.;
            }
            
        }
    }
    
#ifdef DEBUG_GMDS
    //=================================
    // FOR INFORMATION PURPOSE ONLY
    //=================================
    int one_lock=0, two_locks=0, three_locks=0, zero_lock=0;
    for(IGMesh::node_iterator it_n = m_mesh->nodes_begin();
        !it_n.isDone(); it_n.next()) {
        Node n = it_n.value();
        math::Vector3d v = m_bnd_constraint[n.id()];
        int nb_zero=0;
        if(v.X()==0)
            nb_zero++;
        if(v.Y()==0)
            nb_zero++;
        if(v.Z()==0)
            nb_zero++;
        if(nb_zero==0)
            three_locks++;
        else if(nb_zero==1)
            two_locks++;
        else if(nb_zero==2)
            one_lock++;
        else if(nb_zero==3)
            zero_lock++;
    }
    Log::mng()<<"---------------------\n";
    Log::mng()<<"locked frame vertex   : "<<two_locks+three_locks<<"\n";
    Log::mng()<<"\t 2-axis lock  "<<two_locks<<"\n";
    Log::mng()<<"\t 3-axis lock  "<<three_locks<<"\n";
    Log::mng()<<"one locked axis vertex: "<<one_lock<<"\n";
    Log::mng()<<"no locked axis vertex : "<<zero_lock<<"\n";
    Log::mng()<<"total number vertices : "<<m_mesh->getNbNodes()<<"\n";
    Log::mng()<<"---------------------\n";;
    Log::mng().flush();
#endif //DEBUG GMDS
}
/*---------------------------------------------------------------------------*/
void PointGenerator::initSystem()
{
    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, NLint(m_nb_unknowns));
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
}
/*---------------------------------------------------------------------------*/
void PointGenerator::buildSystem()
{
    nlBegin(NL_SYSTEM);
    //======================================================================
    // First, boundary constraint are pushed into the system
    //======================================================================
    for(auto n_id:m_mesh->nodes()) {
        
        int  i = m_id[n_id];
        
        if(m_mesh->isMarked<Node>(n_id,m_bm.mark_node_on_surf)||
           m_mesh->isMarked<Node>(n_id,m_bm.mark_node_on_curv)||
           m_mesh->isMarked<Node>(n_id,m_bm.mark_node_on_pnt)){
            math::Vector3d ci= m_bnd_constraint[n_id];
            for (int d = 0; d < 3; d++) {
                if (ci[d] == 1.0) {
                    nlSetVariable (6 * i + 2 * d    , 1.); //cos=1
                    nlSetVariable (6 * i + 2 * d + 1, 0.); //sin=0
                    nlLockVariable(6 * i + 2 * d );
                    nlLockVariable(6 * i + 2 * d + 1);
                }
            }//for (int d = 0; d < 3; d++)
            
        }//if(m_mesh->isMarked(n,m_bm.mark_node_on_surf))
    }

    nlBegin(NL_MATRIX);
    //======================================================================
    // Second parametrization constraint to be solved are added
    //======================================================================
    // The system is based on edges
    for(auto e_id:m_mesh->edges()) {
        
        Edge e = m_mesh->get<Edge>(e_id);
        std::vector<Node> e_nodes = e.get<Node>();
        int i = m_id[e_nodes[0].id()];
        int j = m_id[e_nodes[1].id()];
        
        Node ni = e_nodes[0];
        Node nj = e_nodes[1];
        math::Chart::Mapping R_ij = getRij(ni,nj);
        math::Vector3d g_ij = computeGijWithCurl(e,ni,nj);

        for (int d = 0; d<3; d++) {
            double c = cos(g_ij[d]);
            double s = sin(g_ij[d]);
            
            int off0 = 6 * i + 2 * d;
            int off1 = 6 * j + 2 * R_ij.getPermutations()[d];

            nlRowScaling(1.);
            nlBegin(NL_ROW); // First line of Eq. (5)
            nlCoefficient(off0    ,  c  );
            nlCoefficient(off0 + 1,  s  );
            nlCoefficient(off1    , -1.0);
            nlEnd(NL_ROW);
            
            nlRowScaling(1.);
            nlBegin(NL_ROW); // Second line of Eq. (5)
            nlCoefficient(off0    , -s);
            nlCoefficient(off0 + 1,  c);
            nlCoefficient(off1 + 1, -R_ij.getDirections()[d]);
            nlEnd(NL_ROW);
        }
    }
    
    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
}
/*---------------------------------------------------------------------------*/
void PointGenerator::cleanSystem(){
    nlDeleteContext(nlGetCurrent());
}
/*---------------------------------------------------------------------------*/
void PointGenerator::extractPoints(){
    
    Log::mng()<<"> Point extraction\n";
    Log::mng().flush();
    
    //======================================================================
    // STEP 1. For each region we extract a set of points
    //  EASY TO BE PARALLELIZED
    //======================================================================
    std::vector<math::Point>     extracted_points;
    std::vector<math::Chart>     extracted_charts;
    std::vector<Cell::Data>      extracted_data;
    std::vector<int>             extracted_types;
    std::vector<int>             extracted_classification;
    std::vector<double>          extracted_spacing;
    std::vector<int>             extracted_surfaces;
    std::vector<int>             extracted_curves;
    std::map<int,math::Vector3d> extracted_normal;

    std::vector<Region> wrong_tets;
    std::vector<int>    wrong_fail;
    
    
    std::map<TCellID, std::pair<int,int> > tet2pnts;
    
    for(auto r_id: m_mesh->regions()) {
        
        Region r = m_mesh->get<Region>(r_id);
        
        int status =extractPoints(r,
                                  extracted_points,
                                  extracted_charts,
                                  extracted_data,
                                  extracted_types,
                                  extracted_classification,
                                  extracted_spacing,
                                  extracted_surfaces,
                                  extracted_curves,
                                  extracted_normal);
        if(status!=0){
            wrong_tets.push_back(r);
            wrong_fail.push_back(status);
        }
    }
    if(!wrong_tets.empty()){
        Log::mng()<<"Nb wrong param. tets= "<<wrong_tets.size()<<"\n";
        Log::mng().flush();
        
        Variable<int>* var_wrong =
        m_mesh->newVariable<int,GMDS_REGION>("WrongParametrization");
        
        for(auto ti : wrong_tets){
            (*var_wrong)[ti.id()]=1;
            
        }
    }
    Log::mng()<<"\t "<<extracted_points.size()<<" extracted points\n";
    Log::mng().flush();

    //======================================================================
    // STEP 2. Extracted points are filtered to merge those, which are
    //         too close to each other.
    //======================================================================
    // Merging parameter in function of the expected spacing
    // This factor (0.2) could be a parameter too???
    double eps = 0.25*m_spacing;
    double eps2 =eps*eps;
    


    //======================================================================
    //containers are resized to the size of extracted points size
    // only a few number will be removed. It is a better strategy than
    // the one of ST vectors which can be more expensive.
    // Moreover, it prevents the algorithm to perform to many memory
    // reallocation
    
    int nb_extracted_points = extracted_points.size();
    
    m_points.clear();
    m_charts.clear();
    m_point_data.clear();
    m_point_types.clear();
    m_point_classification.clear();
    m_point_curves.clear();
    m_point_surfaces.clear();
    m_point_surface_normals.clear();
    m_point_spacing.clear();
    m_point_vertex_to_curves.clear();
    
    m_points.reserve(nb_extracted_points);
    m_charts.reserve(nb_extracted_points);
    m_point_data.reserve(nb_extracted_points);
    m_point_types.reserve(nb_extracted_points);
    m_point_classification.reserve(nb_extracted_points);
    m_point_curves.reserve(nb_extracted_points);
    m_point_surfaces.reserve(nb_extracted_points);
    m_point_spacing.reserve(nb_extracted_points);
    

    Log::mng()<<"> Point filtering to merge close points (eps= "<<eps<<")\n";
    //for each created point, we store the indices of the points
    // that are close to them
    for (int i=0; i<nb_extracted_points; i++){
        
        bool already_in=false;
        math::Point pi = extracted_points[i];
        for (int j=i+1;j<extracted_points.size() && !already_in; j++) {
            math::Vector3d vij(pi,extracted_points[j]);
            if (vij.norm2()< eps2)  {
                //we only keep element with a smaller (or equal classification)
                if(extracted_classification[j]<=
                   extracted_classification[i]){
                    already_in=true;
                }
            }
        }
        
        if (!already_in){
            m_points.push_back(pi);
            
            m_point_data.push_back(extracted_data[i]);
            m_charts.push_back(extracted_charts[i]);
            m_point_types.push_back(extracted_types[i]);
            m_point_classification.push_back(extracted_classification[i]);
            m_point_curves.push_back(extracted_curves[i]);
            m_point_surfaces.push_back(extracted_surfaces[i]);
            if(extracted_classification[i]==2)
                m_point_surface_normals[m_points.size()-1]=extracted_normal[i];
            m_point_spacing.push_back(extracted_spacing[i]);
            
            if(extracted_classification[i]==0){
                //We get the curves the current point belongs to
            }
            
        }
    }

    Log::mng()<<"\t "<<m_points.size()
              <<" points after filtering (eps="<<eps<<")"
              <<" with spacing: "<<m_spacing<<"and nb charts:"<<m_charts.size()<<" done\n";
    Log::mng().flush();
    std::cout<<"\t "<<m_points.size()
              <<" points after filtering (eps="<<eps<<")"
              <<" with spacing: "<<m_spacing<<"and nb charts:"<<m_charts.size()<<" done\n";
    std::cout.flush();

}
/*---------------------------------------------------------------------------*/
std::vector<Region>
PointGenerator::getTetSharingOneNode(const Region& ATet)
{
    //===================================================================
    // (1) We get all the regions sharing a node with ATet by filling in
    // a set of regions
    //===================================================================
    std::set<Region> candidates;
    std::vector<Node> adj_nodes = ATet.get<Node>();
    for(unsigned int i=0; i<adj_nodes.size();i++){
        std::vector<Region> current_reg = adj_nodes[i].get<Region>();
        candidates.insert(current_reg.begin(),current_reg.end());
    }
    //===================================================================
    // (2) We extract a vector of regions that do not contain ATeet
    //===================================================================
    std::vector<Region> adj;
    int index=0;
    adj.resize(candidates.size()-1); //-1 as ATet must be removed
    for(std::set<Region>::iterator itr = candidates.begin();
        itr!=candidates.end(); itr++){
        if(itr->id()!=ATet.id())
            adj[index++]=*itr;
    }
    return adj;
}

/*---------------------------------------------------------------------------*/
bool PointGenerator::hasNodes(const gmds::Face& AF,
                              const gmds::TCellID AN0,
                              const gmds::TCellID AN1,
                              const gmds::TCellID AN2)
{
    std::vector<TCellID> n = AF.getIDs<Node>();
    return ((n[0]==AN0 && n[1]==AN1 && n[2]==AN2) ||
            (n[0]==AN0 && n[2]==AN1 && n[1]==AN2) ||
            (n[1]==AN0 && n[0]==AN1 && n[2]==AN2) ||
            (n[1]==AN0 && n[2]==AN1 && n[0]==AN2) ||
            (n[2]==AN0 && n[1]==AN1 && n[0]==AN2) ||
            (n[2]==AN0 && n[0]==AN1 && n[1]==AN2) );

}
/*---------------------------------------------------------------------------*/
int PointGenerator::
extractPoints(const Region&                         ATet,
              std::vector<math::Point>&             APnts,
              std::vector<math::Chart>&             ACharts,
              std::vector<Cell::Data>&              AData,
              std::vector<int>&                     ATypes,
              std::vector<int>&                     AClass,
              std::vector<double>&                  ASpacing,
              std::vector<int>&                     ASurf,
              std::vector<int>&                     ACurv,
              std::map<int,gmds::math::Vector3d>&   ASurfNormal)
{
    //===============================================================
    // We compute the type of tet we are in for output purpose
    //===============================================================
    int region_type = getSingularityType(ATet);
    int pnt_type=0;
    if(0<region_type && region_type<5)
        pnt_type=1;
    else if(region_type>=5)
        pnt_type=2;
    
   
    std::vector<Node> n = ATet.get<Node>();
    std::vector<Edge> e = ATet.get<Edge>();
    
    
    //===============================================================
    // If our tet has at least one node on the boundary, we need to
    // build local edge and face structure for computing point
    // classification
    //===============================================================
    bool need_classification =false;
    if(m_mesh->isMarked(n[0], m_bm.mark_node_on_surf) ||
       m_mesh->isMarked(n[1], m_bm.mark_node_on_surf) ||
       m_mesh->isMarked(n[2], m_bm.mark_node_on_surf) ||
       m_mesh->isMarked(n[3], m_bm.mark_node_on_surf) ){
        need_classification=true;
    }
    
    std::vector<Edge> local_edges;
    std::vector<Face> local_faces;
    local_edges.resize(6);
    local_faces.resize(4);
    
    if(need_classification){
        std::vector<Face> f = ATet.get<Face>();
        
        //We build edge O1, O2, O3
        TCellID ref_id = n[0].id();
        for(int i_n=0; i_n<3; i_n++){
            int cur_id = n[i_n+1].id();
            bool found_edge=false;
            for(int i_e=0; i_e<6 && !found_edge; i_e++){
                Edge ei = e[i_e];
                
                std::vector<TCellID> ei_nodes = ei.getIDs<Node>();
                
                if((ei_nodes[0]==ref_id && ei_nodes[1]==cur_id) ||
                   (ei_nodes[0]==cur_id && ei_nodes[1]==ref_id)){
                    found_edge=true;
                    local_edges[i_n]=ei;
                }
            }
        }
        //now edges 12 and 23
        ref_id = n[1].id();
        for(int i_n=0; i_n<2; i_n++){
            int cur_id = n[i_n+2].id();
            bool found_edge=false;
            for(int i_e=0; i_e<6 && !found_edge; i_e++){
                Edge ei = e[i_e];
                
                std::vector<TCellID> ei_nodes = ei.getIDs<Node>();
                
                if((ei_nodes[0]==ref_id && ei_nodes[1]==cur_id) ||
                   (ei_nodes[0]==cur_id && ei_nodes[1]==ref_id)){
                    found_edge=true;
                    local_edges[3+i_n]=ei;
                }
            }
        }
        
        // and edge 23
        ref_id = n[2].id();
        int cur_id = n[3].id();
        bool found_edge=false;
        for(int i_e=0; i_e<6 && !found_edge; i_e++){
            Edge ei = e[i_e];
            
            std::vector<TCellID> ei_nodes = ei.getIDs<Node>();
            
            if((ei_nodes[0]==ref_id && ei_nodes[1]==cur_id) ||
               (ei_nodes[0]==cur_id && ei_nodes[1]==ref_id)){
                found_edge=true;
                local_edges[5]=ei;
            }
        }
        // and finally local faces
        for(int i_f=0; i_f<4; i_f++){
            Face fi = f[i_f];
            if(hasNodes(fi, n[1].id(),n[2].id(),n[3].id()))
                local_faces[0]=fi;
            else if(hasNodes(fi, n[0].id(),n[2].id(),n[3].id()))
                local_faces[1]=fi;
            else if(hasNodes(fi, n[0].id(),n[1].id(),n[3].id()))
                local_faces[2]=fi;
            else if(hasNodes(fi, n[0].id(),n[1].id(),n[2].id()))
                local_faces[3]=fi;
        }
    }

    //===============================================================
    // Build local variables to work on. All the computation will
    // be done with node n[0] as the reference. So we need to get
    // data relatively to n[0] (like edges coming from this node)
    //===============================================================
    math::Vector3d lU    [4]; // parametric location of each tet corner
    math::Chart    lChart[4]; //chart corresponding to each tet corner
    math::Point    lX    [4]; // geometric location of each tet corner
    TCellID        n_id  [4]; // ids of each tet corner
    
    // init from stored
    for (int i=0; i<4; i++){
        n_id  [i] = n[i].id();
        lX    [i] = n[i].point();
        lU    [i] = m_Ui[n_id[i]];
        lChart[i] = (*m_rotation_field)[n[i].id()].toChart();
    }
    //===============================================================
    //We need to build oriented edges from the reference point 0 used
    // for computations
    //===============================================================
    std::vector<OrientedEdge> oe;
    oe.resize(3);
    for(int i=0;i<3;i++){
        int i0 = n_id[0];
        int i1 = n_id[i+1];
        
        bool found_edge = false;
        Edge ei;
        for(int j=0; j<e.size() && !found_edge; j++){
            Edge ej = e[j];
            
            std::vector<TCellID> ej_nodes = ej.getIDs<Node>();
            
            if((ej_nodes[0]==i0 && ej_nodes[1]==i1) ||
               (ej_nodes[0]==i1 && ej_nodes[1]==i0)){
                ei=ej;
                found_edge=true;
            }
        }
        OrientedEdge oei(ei,n[0],n[i+1]);
        oe[i]=oei;
    }
    
    //===============================================================
    // transports all the value in nodes 1, 2, 3 to the reference
    // chart of vertex 0.
    //===============================================================
    for (int i=0; i<3; i++) {
        //We look for the edge from vO to v[i+1]
        OrientedEdge oei = oe[i];
        
        math::Vector3d g_0i = computeGijWithCurl(oei.edge,n[0],n[i+1])/(2*M_PI);
        
        lU[i+1]= getRij(n[0],n[i+1]) * lU[i+1];
        
        alignTij(lU[i+1], m_Ui[n_id[0]], g_0i);
    }
    //===============================================================
    // We store the four parametric points corresponding to our tet
    // nodes, all seen from the n[0] point of view
    const math::Point PU[4] ={
        math::Point(lU[0].X(),lU[0].Y(),lU[0].Z()),
        math::Point(lU[1].X(),lU[1].Y(),lU[1].Z()),
        math::Point(lU[2].X(),lU[2].Y(),lU[2].Z()),
        math::Point(lU[3].X(),lU[3].Y(),lU[3].Z())
    };

    //===============================================================
    // computation of an integer-value bounding box in the parametric
    // domain U
    //===============================================================
    // X, Y, Z min coordinate of the bounding box
    int bb_min[3] = {INT_MAX, INT_MAX, INT_MAX};
    // X, Y, Z max coordinate of the bounding box
    int bb_max[3] = {INT_MIN, INT_MIN, INT_MIN};
    
    for (int i_n=0;i_n<4;i_n++) {
        for (int i_coord=0; i_coord<3; i_coord++){
            int floor_val = std::floor(lU[i_n][i_coord]);
            int ceil_val  = std::ceil (lU[i_n][i_coord]);
            
            if(floor_val-1<bb_min[i_coord])
                bb_min[i_coord] = floor_val-1;
            
            if(ceil_val+1>bb_max[i_coord])
                bb_max[i_coord] = ceil_val+1;
        }
    }
    
    
    //===============================================================
    // Considering all the integer point located in the bounding box,
    // we look for those that belongs to the parametric tet and we
    // project them bach to the geometric space.
    //
    // At the end of the projection, a corrective projection is
    // applied for the points that could be generated outside of the
    // domain
    //===============================================================
    int nb_failures = 0;
    for (int x=bb_min[0]-1; x<=bb_max[0]+1; x++){
        for (int y=bb_min[1]-1; y<=bb_max[1]+1; y++){
            for (int z=bb_min[2]-1; z<=bb_max[2]+1; z++){
                
                math::Point    p_param(x,y,z);
                math::Point    p_physic;
                math::Vector4d coord; //init to (0,0,0,0) by default
                try{
                    if(computeBarycenter(p_param, PU, coord)) {
                        
                        //=======================================
                        //We compute the physical point
                        p_physic = (coord[0]*lX[0] +
                                    coord[1]*lX[1] +
                                    coord[2]*lX[2] +
                                    coord[3]*lX[3]);
                        
                        //======================================
                        //We compute the corresponding chart
                        std::vector<math::Quaternion> qs;
                        qs.resize(4);
                        qs[0]= math::Quaternion(lChart[0]);
                        qs[1]= math::Quaternion(lChart[1]);
                        qs[2]= math::Quaternion(lChart[2]);
                        qs[3]= math::Quaternion(lChart[3]);
                        
                        std::vector<TCoord> ws;
                        ws.resize(4);
                        ws[0]=coord[0];
                        ws[1]=coord[1];
                        ws[2]=coord[2];
                        ws[3]=coord[3];
                        
                        math::Quaternion q = math::Quaternion::mean(qs,
                                                                    ws);
                        
                        
                        //======================================
                        //We deduce the geometric classification
                        int cl=3;
                        int surf_id=-1;
                        int curv_id =-1;
                        TCellID node_id=-1;
                        math::Vector3d nv;
                        if(need_classification){
                            //projection is performed here too!!
                            cl = computeClassification(coord,
                                                       n,
                                                       local_edges,
                                                       local_faces,
                                                       surf_id,
                                                       curv_id,
                                                       node_id,
                                                       nv);
                        }
                        //======================================
                        //We update the containers
                        APnts.push_back  (p_physic);
                        ATypes.push_back (pnt_type);
                        AData.push_back(Cell::Data(3,ATet.id()));
                        ACharts.push_back(math::Chart(q));
                        AClass.push_back(cl);
                        ACurv.push_back(curv_id);
                        ASurf.push_back(surf_id);
                        ASpacing.push_back(m_spacing);
                        if(cl==2){
                            ASurfNormal[APnts.size()-1]=nv;
                        }
                        if(cl==0){
                            //classification on a geometric point
                            m_point_vertex_to_node[APnts.size()-1]=node_id;
                        }

                    }
                }
                catch(GMDSException& e){
                    nb_failures++;
                }
            }
        }
    }
    return nb_failures;
}
/*---------------------------------------------------------------------------*/
bool PointGenerator::computeBarycenter(const math::Point& APntParam,
                                       const math::Point* ATetParam,
                                       math::Vector4d& ACoord) const
{
    //==============================================================
    // First we compute the location of APntParam ito ATetParam
    //==============================================================
    double coeff[4]={0, 0, 0, 0};
    
    math::Point::computeBarycentric(ATetParam[0],ATetParam[1],
                                    ATetParam[2],ATetParam[3],
                                    APntParam,
                                    coeff[0],coeff[1],
                                    coeff[2],coeff[3]);
    
    //==============================================================
    // If APntParam is "almost" into ATetParam, then we provide its
    // Barycentric coordinates
    //==============================================================
    double tolerance = -0.01;
    if(coeff[0]<tolerance || coeff[1]<tolerance ||
       coeff[2]<tolerance || coeff[3]<tolerance)
        return false;
    
    ACoord[0] = coeff[0];
    ACoord[1] = coeff[1];
    ACoord[2] = coeff[2];
    ACoord[3] = coeff[3];
    
    return true;
}

/*---------------------------------------------------------------------------*/
int PointGenerator::getSingularityType(const Face& AF)
{
    if(isFFSingular(AF)){
       return 2;
    }
    else if (isPGPSingular((AF))){
        return 1;
    }
    return 0;
}
/*---------------------------------------------------------------------------*/
bool PointGenerator::isFFSingular(const Face& AF)
{
    std::vector<TCellID> n = AF.getIDs<Node>();
    math::Chart::Mapping m01 = getRij(n[0], n[1]);
    math::Chart::Mapping m12 = getRij(n[1], n[2]);
    math::Chart::Mapping m20 = getRij(n[2], n[0]);
    
    return !(m20*m12*m01).isIdentity();
}


/*---------------------------------------------------------------------------*/
bool PointGenerator::isOnCurve(const Node& AN)
{
    return (m_mesh->isMarked(AN, m_bm.mark_node_on_curv) ||
            m_mesh->isMarked(AN, m_bm.mark_node_on_pnt));
}

/*---------------------------------------------------------------------------*/
int PointGenerator::computeClassification(const math::Vector4d&     ACoord,
                                          const std::vector<Node>&  AN,
                                          const std::vector<Edge>&  AE,
                                          const std::vector<Face>&  AF,
                                          int &                     ASID,
                                          int &                     ACID,
                                          TCellID&                  ANID,
                                          math::Vector3d&           ANormal)
{
    math::Vector4d coord = ACoord;
    double tol = 0.1;
    
    //this 3 boolean marks indicate if the point is located on a node, edge
    // or face of the current tet.
    bool on_node = false;
    bool on_edge = false;
    bool on_face = false;

    //index of the found cell
    int index = -1;
    //ON POINTS
    if(abs(coord[0])<tol && abs(coord[1])<tol && abs(coord[2])<tol){
        on_node = true;
        index   = 3;
    }
    else if(abs(coord[0])<tol && abs(coord[1])<tol && abs(coord[3])<tol){
        on_node = true;
        index   = 2;
    }
    else if(abs(coord[0])<tol && abs(coord[2])<tol && abs(coord[3])<tol){
        on_node = true;
        index   = 1;
    }
    else if(abs(coord[1])<tol && abs(coord[2])<tol && abs(coord[3])<tol){
        on_node = true;
        index   = 0;
    } // THEN ON EDGES
    else if(abs(coord[0])<tol && abs(coord[1])<tol){
         // point generated on edge [n[2],n[3]]
        on_edge = true;
        index   = 5;
    }
    else if(abs(coord[0])<tol && abs(coord[2])<tol){
        // point generated on edge [n[1],n[3]]
        on_edge = true;
        index   = 4;
        
    }
    else if(abs(coord[0])<tol && abs(coord[3])<tol){
        // point generated on edge [n[1],n[2]]
        on_edge = true;
        index   = 3;
    }
    else if(abs(coord[1])<tol && abs(coord[2])<tol){
        // point generated on edge [n[0],n[3]]
        on_edge = true;
        index   = 2;
    }
    else if(abs(coord[1])<tol && abs(coord[3])<tol){
        // point generated on edge [n[0],n[2]]
        on_edge = true;
        index   = 1;
    }
    else if(abs(coord[2])<tol && abs(coord[3])<tol){
        // point generated on edge [n[0],n[1]]
        on_edge = true;
        index   = 0;
    }
    // THEN ON FACE
    else if(abs(coord[0])<tol){
        // point on face [n[1],n[2],n[3]]
        on_face = true;
        index   = 0;
    }
    else if(abs(coord[1])<tol){
        // point on face [n[0],n[2],n[3]]
        on_face = true;
        index   = 1;
    }
    else if(abs(coord[2])<tol){
        // point on face [n[0],n[1],n[3]]
        on_face = true;
        index   = 2;
    }
    else if(abs(coord[3])<tol){
        // point on face [n[0],n[1],n[2]]
        on_face = true;
        index   = 3;
    }
    if(on_node){
        Node n = AN[index];
        
        if(m_mesh->isMarked(n,m_bm.mark_node_on_pnt)){
            ANID = n.id();
            return 0;
        }
        else if(m_mesh->isMarked(n,m_bm.mark_node_on_curv)){
            ACID = (*m_curve_node_color)[n.id()];
            return 1;
        }
        else if(m_mesh->isMarked(n,m_bm.mark_node_on_surf)){
            ASID = (*m_surface_node_color)[n.id()];
            ANormal = m_normal[n.id()];
            ANormal.normalize();
            return 2;
        }
    }
    else if (on_edge){
        Edge e = AE[index];
        std::vector<Node> e_nodes = e.get<Node>();

        if(m_mesh->isMarked(e,m_bm.mark_edge_on_curv)){
            ACID = (*m_curve_edge_color)[e.id()];

            return 1;
        }
        else if(m_mesh->isMarked(e,m_bm.mark_edge_on_surf)){

            int c0 =(*m_surface_node_color)[e_nodes[0].id()];
            int c1 =(*m_surface_node_color)[e_nodes[1].id()];
            if(c0!=0)
                ASID=c0;
            else if(c1!=0)
                ASID=c1;
            else{
                //it means that the edge connect to nodes lying
                //on curves, so we use adjacent bnd faces
                std::vector<Face> fes = e.get<Face>();
                //we look for the first boundary face
                Face bnd_f;
                auto found_bnd = true;
                for(auto i_f=0; i_f<fes.size() && found_bnd; i_f++){
                    Face fi = fes[i_f];
                    if(m_mesh->isMarked(fi,m_bm.mark_face_on_surf)){
                        bnd_f=fi;
                        found_bnd=true;
                    }
                }
                if(!found_bnd) {
                    throw GMDSException("Surface classification issue");
                }
                ASID = (*m_surface_face_color)[bnd_f.id()];
                
            }
            if(isOnCurve(e_nodes[0]) && !isOnCurve(e_nodes[1])){
                ANormal =m_normal[e_nodes[1].id()];
            }
            else if(isOnCurve(e_nodes[1]) && !isOnCurve(e_nodes[0])){
                ANormal =m_normal[e_nodes[0].id()];
            }
            else
                ANormal = 0.5*(m_normal[e_nodes[0].id()]+
                               m_normal[e_nodes[1].id()]);
            
            ANormal.normalize();
         
            return 2;
        }
    }
    else if(on_face){
        Face f = AF[index];
        if(m_mesh->isMarked(f,m_bm.mark_face_on_surf)){
            ASID = (*m_surface_face_color)[f.id()];
            ANormal = getOutputNormalOfABoundaryFace(f);
            ANormal.normalize();
         
            return 2;
        }
    }
    return 3;
}

/*----------------------------------------------------------------------------*/
math::Vector3d PointGenerator::
getOutputNormal(Face& AFace, Region& ARegion)
{
    std::vector<Node> region_nodes = ARegion.get<Node>();
    std::vector<Node> face_nodes = AFace.get<Node>();
    
    if (region_nodes.size() != 4)
        throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on tetrahedral regions");
    if (face_nodes.size() != 3)
        throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on triangular faces");
    
    //we go through all the nodes of ARegion to find the one that do not belong
    //to AFAce
    for (auto n: region_nodes)   {
        if (n != face_nodes[0] && n != face_nodes[1] && n != face_nodes[2]) {
            //n is the node opposite to the face AFace
            Node n0 = face_nodes[0];
            Node n1 = face_nodes[1];
            Node n2 = face_nodes[2];
            math::Vector normal_to_face = AFace.normal();
            math::Vector in_vector(n0.point(), n.point());
            if (normal_to_face.dot(in_vector)>0.0)
            {
                return math::Vector3d(-normal_to_face.get(0),
                                      -normal_to_face.get(1),
                                      -normal_to_face.get(2));
            }
            else
            {
                return normal_to_face;
            }
            
        } //if (n != face_nodes[0] && n != face_nodes[1] && n != face_nodes[2])
        
    }//for (unsigned int i = 0; i<region_nodes.size(); i++)
    return math::Vector3d(0, 0, 0);
}
/*----------------------------------------------------------------------------*/
math::Vector3d PointGenerator::
getOutputNormalOfABoundaryFace(Face& AFace)
{
    std::vector<Region> adj_regions = AFace.get<Region>();
    if (adj_regions.size() != 1)
        throw GMDSException("A boundary face must be adjacent to only 1 region!!!");
    
    return getOutputNormal(AFace, adj_regions[0]);
}

/*---------------------------------------------------------------------------*/
bool PointGenerator::isPGPSingular(const Face& AF)
{
    //We get the nodes of the faces
    std::vector<Node> n = AF.get<Node>();
    //We build corresponding oriented edges
    std::vector<OrientedEdge> oe = getOrientedEdge(AF);
    
    //===========================================================
    // STEP 1 - We compute the difference between U[2] and U[1]
    // in the basis of 1
    //===========================================================
    //We get the parameteriation value in each point
    math::Vector3d lU[3] ={
        m_Ui[n[0].id()],
        m_Ui[n[1].id()],
        m_Ui[n[2].id()]
    };
    
    Node n1 = n[1], n2 = n[2];
    
    //We get the Gij vector along edge [1,2]
    math::Vector3d g_12 = computeGijWithCurl(oe[1].edge,n1,n2)/(2.*M_PI);

    //We mode U[2] in the reference chart of U[1]
    lU[2] = getRij(n1,n2) * lU[2];
    
    //We compute the Tij translation
    alignTij(lU[2], lU[1], g_12);
    
    double diff_12 = (lU[1] - lU[2]).norm();

    //===========================================================
    // STEP 1 - We compute the difference between U[2] and U[1]
    // in the basis of 0 now
    //===========================================================
    //lU has been modified, we reinitalize it
    lU[0] = m_Ui[n[0].id()];
    lU[1] = m_Ui[n[1].id()];
    lU[2] = m_Ui[n[2].id()];

    //We compute lU[1] and lU[2] from the point of view of n[0]
    for (int i = 1; i < 3; i++) {
        Node n0 = n[0], ni = n[i];
        Edge ei =(i==1)?oe[0].edge:oe[2].edge;
        
        OrientedEdge oei(ei,n0,ni);
        math::Vector3d g_0i = computeGijWithCurl(ei,n0,ni)/(2.*M_PI);
 
        lU[i] = getRij(n0,ni) * lU[i];
        alignTij(lU[i],lU[0],g_0i);
    }
    
    double diff_12_from_0 = (lU[1] - lU[2]).norm();
    
    return(abs(diff_12_from_0 - diff_12) > .001);
    
}
/*---------------------------------------------------------------------------*/
void PointGenerator::alignTij(math::Vector3d& AToBeAligned,
                              const math::Vector3d& ARef,
                              const math::Vector3d& AGeomDeviation)
{
    math::Vector3d ref = ARef-AGeomDeviation;
    for (int i = 0; i < 3; i++){
        while (AToBeAligned[i] - ref[i] > .5) {
            AToBeAligned[i] -= 1.0;
        }
        while (AToBeAligned[i] - ref[i] < -.5) {
            AToBeAligned[i] += 1.0;
        }
    }
}
/*---------------------------------------------------------------------------*/
int PointGenerator::getSingularityType(const Region& AR)
{
    std::vector<Face> f = AR.get<Face>();
    
    int sing_pgp=0;
    int sing_ff =0;
    for(int i=0; i<4; i++){
        int si = getSingularityType(f[i]);
        if(si==1)
            sing_pgp++;
        else if(si==2)
            sing_ff++;
    }
    if(sing_ff!=0){
        return sing_ff+4; //offset of 4 since 1 to 4 are for PGP sing
    }
    else{
        return sing_pgp;
    }
}
/*---------------------------------------------------------------------------*/
void PointGenerator::buildANNTree(const std::vector<math::Point>&     APnts,
                                  std::vector<math::Chart>&           ACharts,
                                  std::vector<Cell::Data>&            AMeshData,
                                  std::vector<int>&                   ATypes,
                                  std::vector<int>&                   AClass,
                                  std::vector<int>&                   ACurv,
                                  std::vector<int>&                   ASurf,
                                  std::map<int,gmds::math::Vector3d>& ANormal,
                                  std::vector<double>&                ASpacing,
                                  std::vector<std::vector<int> >&     AOutPnts)
{
    AOutPnts.clear();
    
    int nb_pnts = APnts.size();
    int	k		= 20;      // max number of nearest neighbors
    int	dim		= 3;       // dimension
    int	maxPts	= nb_pnts; // maximum number of data points
    
    int					nPts;					// actual number of data points
    ANNpointArray		dataPts;				// data points
    ANNpoint			queryPt;				// query point
    ANNidxArray			nnIdx;					// near neighbor indices
    ANNdistArray		dists;					// near neighbor distances
    ANNkd_tree*			kdTree;					// search structure
    
    
    queryPt = annAllocPt(dim);					// allocate 1 query point
    dataPts = annAllocPts(maxPts, dim);			// allocate data points
    nnIdx = new ANNidx[k];						// allocate near neigh indices
    dists = new ANNdist[k];						// allocate near neighbor dists

    //========================================================
    // (1) Fill in the  ANN structure for storing points
    //
    // Important: Points in APnts and dataPnts are stored with
    // same index.
    //========================================================
    nPts = 0;
    std::cout << "Data Points:\n";
    while (nPts < maxPts) {
        math::Point p = APnts[nPts];
        dataPts[nPts][0] = p.X();
        dataPts[nPts][1] = p.Y();
        dataPts[nPts][2] = p.Z();
        nPts++;
    };
    //========================================================
    // (2) Build the search structure
    //========================================================
    kdTree = new ANNkd_tree(dataPts,	// the data points
                            nPts,		// number of points
                            dim);		// dimension of space

    //========================================================
    // (2) Search
    //========================================================
    double eps = 0.25*m_spacing;
    double eps2 =eps*eps;
    
    std::vector<bool> flag(nb_pnts,false);

    // For each point, we look for
    for(int i=0; i<nb_pnts;i++){
        if(flag[i]==true){
            //already in the vicinity of a point
            continue;
        }
        math::Point pi = APnts[i];
        queryPt = dataPts[i];
        // annkFRSearch returns the number of points at a
        // squared distance of queryPt minus k. Putting
        // k=0, we get the exact number of pnts in the
        // eps2-radius ball and we get them in a second
        // step
//        int nb_pnts =  kdTree->annkFRSearch(queryPt,// query point
//                                            eps2, //  max. squared radius
//                                            0,  //k
//                                            nnIdx,
//                                            dists);
        
        kdTree->annkSearch(		// search
                           queryPt,// query point
                           k,
                           nnIdx,
                           dists,
                           0.1);
        //            kdTree->annkFRSearch(queryPt,// query point
        //                                 eps2, //  max. squared radius
        //                                 nb_pnts,
        //                                 nnIdx, //point index
        //                                 dists //distance
        //                                 );
        //some other points in the vicinity!
        std::vector<int> close_pi;
        bool stop = false;
        for(int i_c=0; i_c<k && !stop; i_c++){
            if(dists[i_c]>eps2)
                stop= true;
            else
            {
                if (dists[i_c]>1e-8){
                    close_pi.push_back(nnIdx[i_c]);
                    flag[nnIdx[i_c]]=true;
                }
            }
        }
        if(close_pi.empty()) {//nothing is close enough
            //no other point in the vicinity ball
            m_points.push_back(pi);
            m_charts.push_back(ACharts[i]);
            m_point_data.push_back(AMeshData[i]);
            m_point_types.push_back(ATypes[i]);
            m_point_classification.push_back(AClass[i]);
            m_point_curves.push_back(ACurv[i]);
            m_point_surfaces.push_back(ASurf[i]);
            if(AClass[i]==2){
                m_point_surface_normals[m_points.size()-1]=ANormal[i];
            }
            m_point_spacing.push_back(ASpacing[i]);
            
        }
        else{
            close_pi.push_back(i);
            
            AOutPnts.push_back(close_pi);
        }
    }
    delete [] nnIdx;							// clean things up
    delete [] dists;
    delete kdTree;
    annClose();									// done with ANN
    

}
/*---------------------------------------------------------------------------*/
void PointGenerator::writeOutput(){
    static int nb_file=0;

    Variable<int>* var_sing = 0;
    Variable<int>* var_cl = 0;
    Variable<int>* var_cl_face = 0;

    try{
        var_sing = m_mesh->getVariable<int,GMDS_REGION>("PG_sing");
        var_cl = m_mesh->getVariable<int,GMDS_NODE>("classification");
        var_cl_face = m_mesh->getVariable<int,GMDS_FACE>("classification");
    }
    catch (GMDSException& e){
        var_sing = m_mesh->newVariable<int,GMDS_REGION>("PG_sing");
        var_cl = m_mesh->newVariable<int,GMDS_NODE>("classification");
        var_cl_face = m_mesh->newVariable<int,GMDS_FACE>("classification");

    }
    
    //=========================================================================
    //  LOOP ON NODES TO GET CLASSIFICATION
    //=========================================================================
    for (auto n_id:m_mesh->nodes()){
        if(m_mesh->isMarked<Node>(n_id, m_bm.mark_node_on_pnt))
            (*var_cl)[n_id] = 0;
        else if(m_mesh->isMarked<Node>(n_id, m_bm.mark_node_on_curv))
            (*var_cl)[n_id] = 1;
        else if(m_mesh->isMarked<Node>(n_id, m_bm.mark_node_on_surf))
            (*var_cl)[n_id] = 2;
        else
            (*var_cl)[n_id] = 3;
    }
    
    //=========================================================================
    //  LOOP ON FACES TO GET CLASSIFICATION
    //=========================================================================
    for (auto f_id:m_mesh->faces()){
        if(m_mesh->isMarked<Face>(f_id, m_bm.mark_face_on_surf))
            (*var_cl_face)[f_id] = 1;
        else
            (*var_cl_face)[f_id] = 0;
    }
    //=========================================================================
    //  LOOP ON REGIONS TO GET ALL SING. TETS
    //=========================================================================
    int nbColoredTet = 0;
    for (auto r_id:m_mesh->regions()){
        Region r = m_mesh->get<Region>(r_id);
        int color = getSingularityType(r);
        if (color!=0 )
            nbColoredTet++;
        (*var_sing)[r_id] = color;
        
        
    }

    gmds::IGMeshIOService ioService(m_mesh);
    gmds::VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N| gmds::F| gmds::R);
    writer.setDataOptions(gmds::N| gmds::F| gmds::R);
    std::string file_name=m_param_gl.output_dir+ "/PG_DEBUG_"+ to_string(nb_file);
    std::cout<<"WRITE "<<file_name<<std::endl;
    writer.write(file_name);
    std::cout<<"done"<<std::endl;
    nb_file++;

}
/*---------------------------------------------------------------------------*/
PointGenerator::OrientedEdge::OrientedEdge()
{;}
/*---------------------------------------------------------------------------*/
PointGenerator::OrientedEdge::OrientedEdge(gmds::Edge& AE,
                                           gmds::Node& AF,
                                           gmds::Node& AS):
edge(AE), first(AF),second(AS){;}
/*---------------------------------------------------------------------------*/
PointGenerator::OrientedEdge::OrientedEdge(gmds::Edge& AE)
{
    edge=AE;
    std::vector<Node> n = edge.get<Node>();
    first  = n[0];
    second = n[1];
}
/*---------------------------------------------------------------------------*/
bool PointGenerator::OrientedEdge::isWellOriented(){
    return (edge.getIDs<gmds::Node>()[0]==first.id());
}

/*---------------------------------------------------------------------------*/


