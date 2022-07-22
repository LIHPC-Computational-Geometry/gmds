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
#include <gmds/frame3d/PGPStructPointExtraction.h>
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <set>
#include <queue>
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
PGPStructPointExtraction::
PGPStructPointExtraction(Mesh* AMesh,
                         const ParamsGlobal& AParamGl,
                         std::map<gmds::TCellID, gmds::math::Vector3d>& ANormal,
                         const ParamsMark& AMarks,
                         std::map<gmds::TCellID, gmds::math::Vector3d >& AUI,
                         const double ASpacing,
                         const double ACurl):
        m_mesh(AMesh),
        m_param_gl(AParamGl),
        m_bm(AMarks),
        m_Ui(AUI),
        m_spacing(ASpacing),
        m_curl(0.35)
{
    m_rotation_field = m_mesh->getVariable<math::AxisAngleRotation, GMDS_NODE>("rotation_field");

}
/*---------------------------------------------------------------------------*/
math::Chart::Mapping PGPStructPointExtraction::getRij(const TCellID AFrom,
                                            const TCellID ATo) const
{
    math::AxisAngleRotation rotation_from = (*m_rotation_field)[AFrom];
    math::AxisAngleRotation rotation_to   = (*m_rotation_field)[ATo];

    math::Chart chart_from  = rotation_from.toChart();
    math::Chart chart_to    = rotation_to.toChart();

    return math::Chart::Mapping(chart_from,chart_to);
}
/*---------------------------------------------------------------------------*/
math::Chart::Mapping PGPStructPointExtraction::getRij(const Node& AFrom,
                                            const Node& ATo) const
{
    return getRij(AFrom.id(),ATo.id());
}
/*---------------------------------------------------------------------------*/
math::Vector3d PGPStructPointExtraction::computeGij(OrientedEdge& AE)  {
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
math::Vector3d PGPStructPointExtraction::computeCij(OrientedEdge& AE)
{
    if(m_curl == 0.0) {
        return math::Vector3d({0, 0, 0});
    }
    
    Edge e = AE.edge;
    TCellID e_id = e.id();
    
    math::Chart::Mapping r_ij = getRij(AE.first,AE.second);
    
    math::Vector3d correction({0, 0, 0});
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
std::vector<PGPStructPointExtraction::OrientedEdge> PGPStructPointExtraction::
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
void PGPStructPointExtraction::execute()
{
    Log::mng()<<"> Point extraction\n";
    Log::mng().flush();

    //======================================================================
    // STEP 1. For each region we extract a set of points
    //  EASY TO BE PARALLELIZED
    //======================================================================
    std::vector<math::Point>     extracted_points;

    std::vector<Region> wrong_tets;
    std::vector<int>    wrong_fail;

    std::map<TCellID, std::pair<int,int> > tet2pnts;

    for(auto r_id: m_mesh->regions()) {

        Region r = m_mesh->get<Region>(r_id);

        int status =extractPoints(r, extracted_points);
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

    m_points = extracted_points;
}
/*---------------------------------------------------------------------------*/
math::Vector3d PGPStructPointExtraction::computeGijWithCurl(Edge& AEdge,
                                                  Node& AFrom,
                                                  Node& ATo)
{
    OrientedEdge oe(AEdge,AFrom,ATo);
    return computeGij(oe) + computeCij(oe);
}
/*---------------------------------------------------------------------------*/
int PGPStructPointExtraction::
extractPoints(const Region&                         ATet,
              std::vector<math::Point>&             APnts)
{
    std::vector<Node> n = ATet.get<Node>();
    std::vector<Edge> e = ATet.get<Edge>();

    std::vector<Edge> local_edges;
    std::vector<Face> local_faces;
    local_edges.resize(6);
    local_faces.resize(4);
    

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
    
    for (auto & u : lU) {
        for (int i=0; i < 3; i++){
            int floor_val = std::floor(u[i]);
            int ceil_val  = std::ceil (u[i]);
            
            if(floor_val-1<bb_min[i])
                bb_min[i] = floor_val - 1;
            
            if(ceil_val+1>bb_max[i])
                bb_max[i] = ceil_val + 1;
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
                        
                        //We update the containers
                        APnts.push_back  (p_physic);


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
bool PGPStructPointExtraction::computeBarycenter(const math::Point& APntParam,
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
void PGPStructPointExtraction::alignTij(math::Vector3d& AToBeAligned,
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
PGPStructPointExtraction::OrientedEdge::OrientedEdge()
{;}
/*---------------------------------------------------------------------------*/
PGPStructPointExtraction::OrientedEdge::OrientedEdge(gmds::Edge& AE,
                                           gmds::Node& AF,
                                           gmds::Node& AS):
edge(AE), first(AF),second(AS){;}
/*---------------------------------------------------------------------------*/
PGPStructPointExtraction::OrientedEdge::OrientedEdge(gmds::Edge& AE)
{
    edge=AE;
    std::vector<Node> n = edge.get<Node>();
    first  = n[0];
    second = n[1];
}
/*---------------------------------------------------------------------------*/
bool PGPStructPointExtraction::OrientedEdge::isWellOriented(){
    return (edge.getIDs<gmds::Node>()[0]==first.id());
}

/*---------------------------------------------------------------------------*/


