/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
#include <gmds/polyblock/PolyBlock2D.h>
#include <gmds/utils/CommonTypes.h>
/*---------------------------------------------------------------------------*/
// Usage of the Eigen Template library
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/math/Numerics.h>
/*----------------------------------------------------------------------------*/
// OpenNL File Headers
#include "OpenNL_psm.h"
#include<Eigen/IterativeLinearSolvers>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
const math::Vector3d PolyBlock2D::REF[]= {
        {1, 0, 0},
        {-1, 0, 0},
        {0, 1, 0},
        {0, -1, 0},
        {0, 0, 1},
        {0, 0, -1}};
/*----------------------------------------------------------------------------*/
PolyBlock2D::
PolyBlock2D(Mesh* AMesh, const int ABndNodeMark, const int ABndEdgeMark)
: m_mesh(AMesh), m_bnd_node_mark(ABndNodeMark), m_bnd_edge_mark(ABndEdgeMark)
{
    m_bnd_face_mark = m_mesh->newMark<Face>();

    for(auto f_id:m_mesh->faces()) {
        std::vector<TCellID> edge_ids = m_mesh->get<Face>(f_id).getIDs<Edge>();
        for(auto e_id:edge_ids){
            if(m_mesh->isMarked<Edge>(e_id,m_bnd_edge_mark)){
               m_mesh->mark<Face>(f_id, m_bnd_face_mark);
            }
        }
    }

    m_bnd_normal_node = m_mesh->newVariable<math::Vector3d, GMDS_NODE>("bnd_normal");
    m_bnd_normal_edge = m_mesh->newVariable<math::Vector3d, GMDS_EDGE>("bnd_normal");
    m_target_XYZ_node = m_mesh->newVariable<math::Vector3d, GMDS_NODE>("target");
    m_target_XYZ_edge = m_mesh->newVariable<math::Vector3d, GMDS_EDGE>("target");
    m_cross_face = m_mesh->newVariable<math::Cross2D, GMDS_FACE>("cross_field");
    m_v1_face = m_mesh->newVariable<math::Vector3d, GMDS_FACE>("V1");
    m_v2_face = m_mesh->newVariable<math::Vector3d, GMDS_FACE>("V2");
    m_x_value = m_mesh->newVariable<double, GMDS_FACE>("X");
    m_y_value = m_mesh->newVariable<double, GMDS_FACE>("Y");
    m_U = m_mesh->newVariable<double, GMDS_NODE>("U");
    m_V = m_mesh->newVariable<double, GMDS_NODE>("V");
    m_patch_number = m_mesh->newVariable<int, GMDS_NODE>("BndPatch");

    computeBndNormals();
    computeTargetXYZ();
}
/*----------------------------------------------------------------------------*/
PolyBlock2D::~PolyBlock2D() {
    m_mesh->unmarkAll<Face>(m_bnd_face_mark);
    m_mesh->freeMark<Face>(m_bnd_face_mark);
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::execute()  {
    std::cout<<"PolyBlock2D - start point"<<std::endl;
    int nb_iterations=1;
    for(auto i=0;i<nb_iterations;i++) {


        //cross field calculation is done first
        std::cout << "==== Stage 1: Polycube-like cross field (face-centered least-square)" << std::endl;
        computeCrossField();
        //UV Parameterization
//    std::cout<<"==== Stage 2: compoute boundary patch"<<std::endl;
        //  computeBoundaryPatchs();
        std::cout << "==== Stage 3: cross-field aligned UV parameterization (node-based least-square)" << std::endl;
        computeUVParam();
        writeDebugFile();
        deformOmega();
        computeBndNormals();

    }
    std::cout<<"PolyBlock2D - end point "<<std::endl;

}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::computeBndNormals() {
    //first we work on edges
    for(auto e_id:m_mesh->edges()){
        Edge e = m_mesh->get<Edge>(e_id);
        if(m_mesh->isMarked(e,m_bnd_edge_mark)){
            //boundary edge
            std::vector<Node> end_nodes = e.get<Node>();
            math::Point p0 = end_nodes[0].point();
            math::Point p1 = end_nodes[1].point();
            math::Vector3d v01=p1-p0;
            v01.normalize();
            math::Vector normal({v01.Y(), -v01.X(), 0});
            //check if we execute inside or outside
            Face adj_face = e.get<Face>()[0];
            math::Vector to_inside=adj_face.center()-e.center();
            if(to_inside.dot(normal)>0){
                normal = normal.opp();
            }
            m_bnd_normal_edge->set(e_id,normal);
        }
        else{
            //internal edge
            m_bnd_normal_edge->set(e_id,math::Vector3d({0, 0, 0}));
        }
    }

    //and we assign to nodes
    for(auto n_id:m_mesh->nodes()){
        Node n = m_mesh->get<Node>(n_id);
        if(m_mesh->isMarked(n,m_bnd_node_mark)){
            //boundary node, we get the adjacent edges that are on the boundary
            // must be 2 in principle
            std::vector<TCellID> n_edges = n.getIDs<Edge>();
            std::vector<TCellID> bnd_edges;
            for(auto e_id:n_edges){
                if(m_mesh->isMarked<Edge>(e_id,m_bnd_edge_mark)){
                    bnd_edges.push_back(e_id);
                }
            }
            math::Vector3d v0 = m_bnd_normal_edge->value(bnd_edges[0]);
            math::Vector3d v1 = m_bnd_normal_edge->value(bnd_edges[1]);
            m_bnd_normal_node->set(n_id,0.5*(v0+v1));
        }
        else{
            //internal edge
            m_bnd_normal_node->set(n_id,math::Vector3d({0, 0, 0}));

        }
    }
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::computeTargetXYZ() {
    for(auto e_id:m_mesh->edges()){
        if(m_mesh->isMarked<Edge>(e_id,m_bnd_edge_mark)) {
            math::Vector3d v = m_bnd_normal_edge->value(e_id);
            auto closest_ref = 0;
            auto closest_dot = REF[0].dot(v);
            for(auto i=1;i<6;i++){
                auto di = REF[i].dot(v);
                if(di>closest_dot){
                    closest_dot=di;
                    closest_ref=i;
                }
            }
            m_target_XYZ_edge->set(e_id, REF[closest_ref]);
        }
       else{
            m_target_XYZ_edge->set(e_id, math::Vector3d({0, 0, 0}));
       }
    }

    //and we assign to nodes
    for(auto n_id:m_mesh->nodes()) {
        Node n = m_mesh->get<Node>(n_id);
        if (m_mesh->isMarked(n, m_bnd_node_mark)) {
            //boundary node, we get the adjacent edges that are on the boundary
            // must be 2 in principle
            std::vector<TCellID> n_edges = n.getIDs<Edge>();
            std::vector<TCellID> bnd_edges;
            for (auto e_id: n_edges) {
                if (m_mesh->isMarked<Edge>(e_id, m_bnd_edge_mark)) {
                    bnd_edges.push_back(e_id);
                }
            }
            math::Vector3d v0 = m_target_XYZ_edge->value(bnd_edges[0]);
            math::Vector3d v1 = m_target_XYZ_edge->value(bnd_edges[1]);
            m_target_XYZ_node->set(n_id, v0);
        } else {
            //internal edge
            m_target_XYZ_node->set(n_id, math::Vector3d({0, 0, 0}));

        }
    }
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::computeCrossField() {
    int nb_iter_max=6;
    for(auto i=1;i<=nb_iter_max;i++){
        writeDebugFile();
        //we recompute all the cross that are define along a boundary
        updateCrossFieldBoundary();
        double alpha = 1-((double)i/((double)nb_iter_max+1));
        solveCrossField(alpha);
    }
    writeDebugFile();
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::updateCrossFieldBoundary() {
    math::Vector3d z_vec({0, 0, 1});
    for(auto f_id:m_mesh->faces()) {
        std::vector<TCellID> edge_ids = m_mesh->get<Face>(f_id).getIDs<Edge>();
        for(auto e_id:edge_ids){
            if(m_mesh->isMarked<Edge>(e_id,m_bnd_edge_mark)){
                math::Vector3d v1 = m_bnd_normal_edge->value(e_id);
                math::Vector3d v2 = v1.cross(z_vec);
                math::Cross2D c(v1,v2);
                m_cross_face->set(f_id,c);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::solveCrossField(const double AAlpha) {
    TInt nb_unknowns = m_mesh->getNbFaces();
    //initialize the NL solver
    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, nb_unknowns);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);

    for (auto f_id: m_mesh->faces()) {
        //we get the value that was computed at the previous iteration
        math::Cross2D ci = m_cross_face->value(f_id);
        nlSetVariable(f_id,ci.referenceAngle());
    }
    nlBegin(NL_MATRIX);

    double alpha = AAlpha;
    double beta = 1-alpha;

    //we initialize the connectivity relations
    for (auto e_id:m_mesh->edges()) {
        if(!m_mesh->isMarked<Edge>(e_id,m_bnd_edge_mark)) {
            Edge e = m_mesh->get<Edge>(e_id);
            std::vector<TCellID> e_faces = e.getIDs<Face>();
            int i = e_faces[0];
            int j = e_faces[1];
            nlBegin(NL_ROW);
            nlCoefficient(i, 1);
            nlCoefficient(j, -1);
            nlRightHandSide(0);
            nlEnd(NL_ROW);
        }
    }
    //we minimize the alignment with the prescribed boundary normal
    for(auto f_id:m_mesh->faces()){
        if(m_mesh->isMarked<Face>(f_id,m_bnd_face_mark)){
            double angle = m_cross_face->value(f_id).referenceAngle();
            nlBegin(NL_ROW);
            nlCoefficient(f_id,  beta);
            nlRightHandSide(beta*angle);
            nlEnd(NL_ROW);
        }
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);

    nlSolve();
    NLint    nb_iter;
    NLdouble elapsed_time;

    nlGetIntegerv(NL_USED_ITERATIONS, &nb_iter);
    nlGetDoublev (NL_ELAPSED_TIME   , &elapsed_time);

    std::cout<<"cross-field computation in "<<elapsed_time
             <<"s and "<<nb_iter<<" iterations with (alpha, beta)=("<<alpha<<", "<<beta<<")"<<std::endl;

    //get the solution now!!!
    for (auto f_id: m_mesh->faces()) {
        double ai = nlGetVariable(f_id);
        m_x_value->set(f_id,cos(ai));
        m_y_value->set(f_id,sin(ai));
        math::Cross2D c(ai);
        m_cross_face->set(f_id,c);
        m_v1_face->set(f_id,c.componentVectors()[0]);
        m_v2_face->set(f_id,c.componentVectors()[1]);
    }

}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::computeUVParam() {
    for(auto i=0; i<1; i++)
        solveUV();
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::solveUV() {
    //computation is performed on the nodes
    TInt  nb_nodes = m_mesh->getNbNodes();
    TInt nb_unknowns = 2*nb_nodes;
    TCellID lock1 = 0;
    TCellID lock2 = 2;

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, nb_unknowns);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);

    nlSetVariable(2*lock1,0);
    nlSetVariable(2*lock1+1,0);
    nlSetVariable(2*lock2,10);
    nlSetVariable(2*lock2+1,10);
    nlLockVariable(2*lock1);
    nlLockVariable(2*lock1+1);
    nlLockVariable(2*lock2);
    nlLockVariable(2*lock2+1);

    nlBegin(NL_MATRIX);
    //we initialize the connectivity relations
    for (auto f_id:m_mesh->faces()) {
        Face f = m_mesh->get<Face>(f_id);
        std::vector<TCellID> nodes = f.getIDs<Node>();
        math::Point zi, zj, zk;
        projectInUV(f_id,zi,zj,zk);
        TCellID i = nodes[0];
        TCellID j = nodes[1];
        TCellID k = nodes[2];
        math::Vector3d ejk=zk-zj;
        math::Vector3d eki=zi-zk;
        math::Vector3d eij=zj-zi;

        nlBegin(NL_ROW);
        nlCoefficient(i, ejk.X());
        nlCoefficient(j, eki.X());
        nlCoefficient(k, eij.X());
        nlRightHandSide(0);
        nlEnd(NL_ROW);

        nlBegin(NL_ROW);
        nlCoefficient(2*i, ejk.Y());
        nlCoefficient(2*j, eki.Y());
        nlCoefficient(2*k, eij.Y());

        nlCoefficient(2*i+1, -ejk.X());
        nlCoefficient(2*j+1, -eki.X());
        nlCoefficient(2*k+1, -eij.X());
        nlRightHandSide(0);
        nlEnd(NL_ROW);

    }
    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);

    nlSolve();
    NLint    nb_iter;
    NLdouble elapsed_time;

    nlGetIntegerv(NL_USED_ITERATIONS, &nb_iter);
    nlGetDoublev (NL_ELAPSED_TIME   , &elapsed_time);

    std::cout<<"... elapsed time "<<elapsed_time<<" in "<<nb_iter<<" iterations "<<std::endl;

    //get the solution now!!!
    for (auto n_id: m_mesh->nodes()) {
        double ui = nlGetVariable(2*n_id);
        double vi = nlGetVariable(2*n_id+1);
        m_U->set(n_id, ui);
        m_V->set(n_id,vi);
        std::cout<<n_id<<" (u,v) = ("<<ui<<", "<<vi<<")"<<std::endl;
    }

}
/*----------------------------------------------------------------------------*/
int PolyBlock2D::computeUVBndEdge(const TCellID AEdgeID){
    if(!m_mesh->isMarked<Edge>(AEdgeID,m_bnd_edge_mark))
        return 0;

    Edge e = m_mesh->get<Edge>(AEdgeID);
    //we get is directional info
    math::Cross2D c = m_cross_face->value(e.getIDs<Face>()[0]);
    std::vector<Node> e_nodes = e.get<Node>();
    math::Vector3d dir=e_nodes[1].point()-e_nodes[0].point();
    dir.normalize();
    if(dir.dot(c.componentVectors()[0])<dir.dot(c.componentVectors()[1])){
        //means aligned with V direction, so
        return 1;
    }
    //so a V edge
    return 2;
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::computeBoundaryPatchs() {
    //TODO Not finished to be implemented
    int mark_done = m_mesh->newMark<Edge>();
    int patch_number=0;
    for(auto e_id:m_mesh->edges()){
        if(m_mesh->isMarked<Node>(e_id,m_bnd_edge_mark)) {
            if(!m_mesh->isMarked<Edge>(e_id, mark_done)){
                //new unvisited bnd edge, we mark it first
                m_mesh->mark<Edge>(e_id,mark_done);
                patch_number+=1;
                Edge current = m_mesh->get<Edge>(e_id);
                std::vector<TCellID> cn_ids = current.getIDs<Node>();
                m_patch_number->set(cn_ids[0],patch_number);
                m_patch_number->set(cn_ids[1],patch_number);
                //we get is directional info
                int bnd_val = computeUVBndEdge(e_id);
                std::vector<Edge> to_visit = getAdjBndEdges(e_id);
                while (!to_visit.empty()){
                    Edge current = to_visit.back();
                    to_visit.pop_back();
                    if(computeUVBndEdge(current.id())==bnd_val){
                        //same patch!!
                        m_mesh->mark(current,mark_done);
                    }
                }
            }
        }
        else{
            m_mesh->mark<Edge>(e_id,mark_done);
        }
    }
    //we unmarl all the nodes
    m_mesh->negateMaskMark<Edge>(mark_done);
    m_mesh->unmarkAll<Edge>(mark_done);
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::writeDebugFile(){
    static int nb_file=1;
    std::string vtk_file ="polyblock2D_"+ to_string(nb_file++)+".vtk";
    IGMeshIOService ioService(m_mesh);
    VTKWriter w(&ioService);
    w.setCellOptions(gmds::N|gmds::F);
    w.setDataOptions(gmds::N|gmds::F);
    w.write(vtk_file);

}

/*----------------------------------------------------------------------------*/
std::vector<Edge> PolyBlock2D::getAdjBndEdges(const TCellID AEdgeId) {
    if(!m_mesh->isMarked<Edge>(AEdgeId, m_bnd_edge_mark))
        throw GMDSException("Error, not a boundary edge");

    Edge e = m_mesh->get<Edge>(AEdgeId);
    std::vector<Node> e_nodes = e.get<Node>();
    std::vector<Edge> bnd_edges;
    for(auto n:e_nodes){
        std::vector<TCellID> n_edge = n.getIDs<Edge>();
        for(auto e_id:n_edge){
            if(e_id!=AEdgeId && m_mesh->isMarked<Edge>(e_id,m_bnd_edge_mark))
                bnd_edges.push_back(m_mesh->get<Edge>(e_id));
        }
    }
    if(bnd_edges.size()!=2)
        throw GMDSException("Error, a boundary edge must be adjacent to exactly 2 boundary edges");

    return bnd_edges;
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::deformOmega() {
    for(auto n_id:m_mesh->nodes()){
        Node n = m_mesh->get<Node>(n_id);
        n.setX(m_U->value(n_id));
        n.setY(m_V->value(n_id));
    }
}
/*----------------------------------------------------------------------------*/
void PolyBlock2D::projectInUV(const TCellID& AFaceID, math::Point& AZ1, math::Point& AZ2, math::Point& AZ3){
    Face f = m_mesh->get<Face>(AFaceID);
    math::Cross2D uv = m_cross_face->value(AFaceID);
    math::Vector3d u({1, 0, 0});// = uv.componentVectors()[0];
    math::Vector3d v({0, 1, 0});// = uv.componentVectors()[1];

    AZ1.setX(0);
    AZ1.setY(0);

    std::vector<Node> ns = f.get<Node>();
    math::Point p0 = ns[0].point();
    math::Point p1 = ns[1].point();
    math::Point p2 = ns[2].point();

    math::Vector3d v01=p1-p0;
    AZ2.setX(u.dot(v01));
    AZ2.setY(v.dot(v01));

    math::Vector3d v02=p2-p0;
    AZ3.setX(u.dot(v02));
    AZ3.setY(v.dot(v02));

}