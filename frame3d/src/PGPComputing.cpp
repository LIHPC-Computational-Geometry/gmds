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
#include <gmds/frame3d/PGPComputing.h>
/*----------------------------------------------------------------------------*/
// OpenNL File Headers
#include "OpenNL_psm.h"
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
PGPComputing::
PGPComputing(gmds::Mesh *AMesh, const ParamsGlobal &AParamGl,
             const ParamsMark &AMarks, const double ASpacing,
             const double ACurl)
: m_mesh(AMesh),
m_param_gl(AParamGl),
m_bm(AMarks),
m_spacing(ASpacing),
m_curl(0.35),
m_nb_unknowns(0)

{
    m_rotation_field = m_mesh->getVariable<math::AxisAngleRotation, GMDS_NODE>("rotation_field");
}
/*---------------------------------------------------------------------------*/
void PGPComputing::execute() {

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

    cleanSystem();
}
/*---------------------------------------------------------------------------*/
void PGPComputing::buildLocalIDs()
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
void PGPComputing::computeCurlCorrection()
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
        m_corr[e_id] = math::Vector3d({0, 0, 0});
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
        math::Vector3d b({0, 0, 0});
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
        math::Vector3d cor_vec({nlGetVariable(3*lid  ),
                               nlGetVariable(3*lid+1),
                               nlGetVariable(3*lid+2)});

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
void PGPComputing::setBoundaryConstraint()
{
    //======================================================================
    // Only boundary nodes are constrained.
    // We initialize their constraint
    //======================================================================
    for(auto n_id : m_mesh->nodes()) {
        if(m_mesh->isMarked<Node>(n_id,m_bm.mark_node_on_pnt)){
            m_bnd_constraint[n_id]={1,1,1};
        }
        else{
            m_bnd_constraint[n_id] = {0,0,0};
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
        math::Vector3d f_normal({nf.X(), nf.Y(), nf.Z()});

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
void PGPComputing::initSystem()
{
    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, NLint(m_nb_unknowns));
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
}
/*---------------------------------------------------------------------------*/
void PGPComputing::buildSystem()
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
void PGPComputing::cleanSystem(){
    nlDeleteContext(nlGetCurrent());
}
/*---------------------------------------------------------------------------*/
void PGPComputing::solveSystem(){

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
void PGPComputing::getUiSolution(){

    for(auto g_id:m_mesh->nodes()) {
        int l_id = m_id[g_id];
        math::Vector3d ui({0, 0, 0});
        double c = (.5 / M_PI);

        for(int k=0; k<3; k++){
            ui[k] = c* atan2(nlGetVariable(6 * l_id + 2 * k + 1),
                             nlGetVariable(6 * l_id + 2 * k ));
        }
        m_Ui[g_id]=ui;
    }
}
/*---------------------------------------------------------------------------*/
math::Vector3d PGPComputing::
computeGijWithCurl(Edge& AEdge,
                   Node& AFrom,
                   Node& ATo)
{
    OrientedEdge oe(AEdge,AFrom,ATo);
    return computeGij(oe) + computeCij(oe);
}

/*---------------------------------------------------------------------------*/
math::Chart::Mapping PGPComputing::getRij(const Node& AFrom,
                                          const Node& ATo) const
{
    return getRij(AFrom.id(),ATo.id());
}
/*---------------------------------------------------------------------------*/
math::Vector3d PGPComputing::
computeGij(PGPComputing::OrientedEdge& AE)  {
    Node ni = AE.first;
    Node nj = AE.second;

    TCellID i = ni.id();
    TCellID j = nj.id();

    math::Vector3d ref_XYZ[3] = {
            math::Vector3d({1., 0., 0.}),
            math::Vector3d({0., 1., 0.}),
            math::Vector3d({0., 0., 1.})
    };
    math::AxisAngleRotation rot_i = (*m_rotation_field)[i];
    math::AxisAngleRotation rot_j = (*m_rotation_field)[j];

    math::Chart::Mapping r_ij = getRij(i,j);

    //vector from pi to pj
    math::Vector3d xij= nj.point()-ni.point();

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
math::Vector3d PGPComputing::
computeCij(PGPComputing::OrientedEdge& AE)
{
    if(m_curl == 0.0) {
        return {0, 0, 0};
    }

    Edge e = AE.edge;
    TCellID e_id = e.id();

    math::Chart::Mapping r_ij = getRij(AE.first,AE.second);

    math::Vector3d correction={0, 0, 0};
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
void PGPComputing::alignTij(math::Vector3d& AToBeAligned,
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
math::Chart::Mapping PGPComputing::getRij(const TCellID AFrom,
                                          const TCellID ATo) const
{
    math::AxisAngleRotation rotation_from = (*m_rotation_field)[AFrom];
    math::AxisAngleRotation rotation_to   = (*m_rotation_field)[ATo];

    math::Chart chart_from  = rotation_from.toChart();
    math::Chart chart_to    = rotation_to.toChart();

    return math::Chart::Mapping(chart_from,chart_to);

}
/*---------------------------------------------------------------------------*/
std::vector<PGPComputing::OrientedEdge>
PGPComputing::getOrientedEdge(const gmds::Face& AF){
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
PGPComputing::OrientedEdge::OrientedEdge()
{;}
/*---------------------------------------------------------------------------*/
PGPComputing::OrientedEdge::OrientedEdge(gmds::Edge& AE,
                                         gmds::Node& AF,
                                         gmds::Node& AS):
        edge(AE), first(AF),second(AS){;}
/*---------------------------------------------------------------------------*/
PGPComputing::OrientedEdge::OrientedEdge(gmds::Edge& AE)
{
    edge=AE;
    std::vector<Node> n = edge.get<Node>();
    first  = n[0];
    second = n[1];
}
/*---------------------------------------------------------------------------*/
bool PGPComputing::OrientedEdge::isWellOriented(){
    return (edge.getIDs<gmds::Node>()[0]==first.id());
}

/*---------------------------------------------------------------------------*/


