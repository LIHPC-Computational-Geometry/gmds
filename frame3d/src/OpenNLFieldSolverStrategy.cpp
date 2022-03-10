/*----------------------------------------------------------------------------*/
// STL File Headers
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
// OpenNL File Headers
#include "OpenNL_psm.h"
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include <gmds/frame3d/OpenNLFieldSolverStrategy.h>
#include "gmds/math/SHarmonicL4.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
OpenNLFieldSolverStrategy::
OpenNLFieldSolverStrategy(const bool AWithCurvature)
:m_with_curvature(AWithCurvature)
{}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::solve(){
    if(m_surface_mode)
        std::cout<<"OpenNL solving within surface mode ("<<m_mesh->getNbNodes()<<", "<<m_surface_nodes->size()<<")"<<std::endl;
    /* Solve and get solution */
    nlSolve();
    NLint    nb_iter;
    NLdouble elapsed_time;
    
    nlGetIntegerv(NL_USED_ITERATIONS, &nb_iter);
    nlGetDoublev (NL_ELAPSED_TIME   , &elapsed_time);
    
    std::cout<<"... elapsed time "<<elapsed_time
    <<" in "<<nb_iter<<" iterations"<<std::endl;

}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::initializeAssembly(){
    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, m_nb_unknowns);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::finalizeAssembly(){
    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::clean()
{
    nlDeleteContext(nlGetCurrent());
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::setX()
{
    /*if(m_surface_mode){
        for (auto n:*m_surface_nodes) {
            int i = (*m_ordering)[n.id()];
            math::SHarmonicL4 shi = (*m_harmonic_field)[n.id()];
            for(int k=0;k<9;k++)
                nlSetVariable(9*i+k, shi[k]);
        }

    }
    else*/
        for (auto n_id: m_mesh->nodes()) {
            int i = (*m_ordering)[n_id];
            math::SHarmonicL4 shi = (*m_harmonic_field)[n_id];
            for (int k = 0; k < 9; k++)
                nlSetVariable(9 * i + k, shi[k]);
        }

}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::getFeasibleSolution()
{
    /*if(m_surface_mode){
        for (auto n:*m_surface_nodes) {
            if (!m_mesh->isMarked(n, m_markNodeLocked)) {

                int local_id = 9 * (*m_ordering)[n.id()];
                double t[9] = {
                        nlGetVariable(local_id),
                        nlGetVariable(local_id + 1),
                        nlGetVariable(local_id + 2),
                        nlGetVariable(local_id + 3),
                        nlGetVariable(local_id + 4),
                        nlGetVariable(local_id + 5),
                        nlGetVariable(local_id + 6),
                        nlGetVariable(local_id + 7),
                        nlGetVariable(local_id + 8)};

                math::Vector9d v(t);

                math::AxisAngleRotation r;

                math::SHarmonicL4 sh = math::SHarmonicL4::closest(v, r);

                math::Vector3d vx(1, 0, 0);
                math::Vector3d vy(0, 1, 0);
                math::Vector3d vz(0, 0, 1);
                vx = r * vx;
                vy = r * vy;
                vz = r * vz;

                (*m_harmonic_field)[n.id()] = sh;
                (*m_chart_field)[n.id()] = math::Chart(math::Vector3d(vx[0], vx[1], vx[2]),
                                                       math::Vector3d(vy[0], vy[1], vy[2]),
                                                       math::Vector3d(vz[0], vz[1], vz[2]));

            }
        }
    }
    else {*/
        for (auto n_id: m_mesh->nodes()) {
            if (!m_mesh->isMarked<Node>(n_id, m_markNodeLocked)) {
                int local_id = 9 * (*m_ordering)[n_id];
                double t[9] = {
                        nlGetVariable(local_id),
                        nlGetVariable(local_id + 1),
                        nlGetVariable(local_id + 2),
                        nlGetVariable(local_id + 3),
                        nlGetVariable(local_id + 4),
                        nlGetVariable(local_id + 5),
                        nlGetVariable(local_id + 6),
                        nlGetVariable(local_id + 7),
                        nlGetVariable(local_id + 8)};

                math::Vector9d v(t);

                math::AxisAngleRotation r;

                math::SHarmonicL4 sh = math::SHarmonicL4::closest(v, r);

                math::Vector3d vx(1, 0, 0);
                math::Vector3d vy(0, 1, 0);
                math::Vector3d vz(0, 0, 1);
                vx = r * vx;
                vy = r * vy;
                vz = r * vz;

                (*m_harmonic_field)[n_id] = sh;
                (*m_chart_field)[n_id] = math::Chart(vx,vy,vz);

            }
        }
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addLocalOptimConstraints(){
    int step = 9*m_nb_free_nodes+2*m_surface_nodes->size();
   /* if(m_surface_mode) {
        for (auto n:*m_surface_nodes) {
            //      See what to do with fixed nodes
            //        if(m_mesh->isMarked(n,m_markNodeFixed))
            //            continue;

            TCellID n_id = n.id();
            int i = (*m_ordering)[n_id];

            math::SHarmonicL4 prevSH = (*m_harmonic_field)[n_id];
            math::SHarmonicL4 ai = prevSH;
            math::SHarmonicL4 cx = math::EX() * ai;
            math::SHarmonicL4 cy = math::EY() * ai;
            math::SHarmonicL4 cz = math::EZ() * ai;
            //We create the smoothing constraint on the 9D rep vectors of i and j
            for (int k = 0; k < 9; k++) {
                nlRowScaling(m_penalty);
                nlBegin(NL_ROW);
                nlCoefficient(9 * i + k, 1.0);
                nlCoefficient(step + 3 * i, -cx[k]);
                nlCoefficient(step + 3 * i + 1, -cy[k]);
                nlCoefficient(step + 3 * i + 2, -cz[k]);
                nlRightHandSide(ai[k]);
                nlEnd(NL_ROW);
            }
        }
    }
    else*/

        for (auto n_id:m_mesh->nodes()){

            int i = (*m_ordering)[n_id];

            math::SHarmonicL4 prevSH = (*m_harmonic_field)[n_id];
            math::SHarmonicL4 ai = prevSH;
            math::SHarmonicL4 cx = math::EX() * ai;
            math::SHarmonicL4 cy = math::EY() * ai;
            math::SHarmonicL4 cz = math::EZ() * ai;
            //We create the smoothing constraint on the 9D rep vectors of i and j
            for (int k = 0; k < 9; k++) {
                nlRowScaling(m_penalty);
                nlBegin(NL_ROW);
                nlCoefficient(9 * i + k, 1.0);
                nlCoefficient(step + (3 * i)    , -cx[k]);
                nlCoefficient(step + (3 * i) + 1, -cy[k]);
                nlCoefficient(step + (3 * i) + 2, -cz[k]);
                nlRightHandSide(ai[k]);
                nlEnd(NL_ROW);
            }
        }

}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addLockedTerms() {
   /* if (m_surface_mode) {

        for (auto ni:*m_surface_nodes) {

            if (m_mesh->isMarked(ni, m_markNodeLocked)) {
                //Variable to be locked
                TCellID i = ni.id();
                int local_i = (*m_ordering)[i];

                math::SHarmonicL4 sh = (*m_harmonic_field)[ni.id()];
                for (int k = 0; k < 9; k++) {
                    nlSetVariable(9 * local_i + k, sh[k]);
                    nlLockVariable(9 * local_i + k);
                }m_surface_mode

            }//if(m_mesh->isMarked(ni,m_markNodeFixed ))

        }//for (IGMesh::node_iterator it_n = m_mesh->nodes_begin()....
    }
    else {*/
        for (auto n_id : m_mesh->nodes()) {
            if (m_mesh->isMarked<Node>(n_id, m_markNodeLocked)) {
                //Variable to be locked

                int local_i = (*m_ordering)[n_id];

                math::SHarmonicL4 sh = (*m_harmonic_field)[n_id];
                for (int k = 0; k < 9; k++) {
                    nlSetVariable(9 * local_i + k, sh[k]);
                    nlLockVariable(9 * local_i + k);
                }

            }//if(m_mesh->isMarked(ni,m_markNodeFixed ))
        }//for (IGMesh::node_iterator it_n = m_mesh->nodes_begin()....


//TODO LOOK WHY???
    nlBegin(NL_MATRIX);
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addBoundaryConstraints(){
    
    Variable<double>* cot_w = m_mesh->getVariable<double,GMDS_EDGE>("cot_weight");

    if(m_with_curvature){
        //curvature must be also taken into account for boundary constraint
        //We use the same penalty term than for boundary constraint
        //TODO check if another weighted term could not be added.
    }
    for(auto ni: *m_surface_nodes){

      //  std::vector<TCellID> e_ni = ni.getIDs<Edge>();
      /*  double wi=0;
        for(auto ei:e_ni){
            wi+=(*cot_w)[ei]/2.0;
        }
        wi=1;*/
        TCellID i   = ni.id();
        int local_i = (*m_ordering)[i];
        
        math::SHarmonicL4 h0 = (*m_H0)[i];
        math::SHarmonicL4 h4 = (*m_H4)[i];
        math::SHarmonicL4 h8 = (*m_H8)[i];
        //We create the smoothing constraint on the 9D rep vectors of i and j
        
        for(int k=0;k<9;k++){
            
            nlBegin(NL_ROW);
            nlCoefficient(9*local_i+k                  ,  m_penalty);
            nlCoefficient(9*m_nb_free_nodes+2*local_i  , -m_penalty*h8[k]);
            nlCoefficient(9*m_nb_free_nodes+2*local_i+1, -m_penalty*h0[k]);
            nlRightHandSide(m_penalty*sqrt(7.0/12.0)*h4[k]);
            nlEnd(NL_ROW);
            
        }
    }


//    std::cout<<"Alignment constraint!!!"<<std::endl;
//    for(IGMesh::node_iterator it_n = m_mesh->nodes_begin();
//        !it_n.isDone(); it_n.next())
//    {
//        Node ni = it_n.value();
//        int local_i = (*m_ordering)[ni.id()];
//        
//        math::Point pi = ni.point();
//        if((!m_mesh->isMarked(ni, m_markNodeLocked)) &&
//           (!m_mesh->isMarked(ni, m_markNodeOnSurf)))
//        {
//            
//            //inner node so
//            Node closest_bnd_node =(*m_surface_nodes)[0];
//            double closest_dist = pi.distance(closest_bnd_node.point());
//            
//            for(unsigned int i=1; i<m_surface_nodes->size(); i++){
//                Node bndi = (*m_surface_nodes)[i];
//                double di = pi.distance(bndi.point());
//                if(di<closest_dist){
//                    closest_dist=di;
//                    closest_bnd_node=bndi;
//                }
//                
//            }//for(unsigned int i=1; i<m_surface_nodes->size(); i++)
//            
//            //We have the closest boundary node, we constraint the local
//            //frame onto the surface normal defined in that bnd node.
//            math::SHarmonicL4 h0 = (*m_H0)[closest_bnd_node.id()];
//            math::SHarmonicL4 h4 = (*m_H4)[closest_bnd_node.id()];
//            math::SHarmonicL4 h8 = (*m_H8)[closest_bnd_node.id()];
//            //We create the smoothing constraint on the 9D rep vectors of i and j
//            
//            for(int k=0;k<9;k++){
//                
//                nlBegin(NL_ROW);
//                nlCoefficient(9*local_i+k                  ,  1);
//                nlCoefficient(9*m_nb_free_nodes+2*local_i  , -h8[k]);
//                nlCoefficient(9*m_nb_free_nodes+2*local_i+1, -h0[k]);
//                nlRightHandSide(sqrt(7.0/12.0)*h4[k]);
//                nlEnd(NL_ROW);
//                
//            }
//        }
//        
//    }//for(IGMesh::node_iterator it_n = m_mesh....
//    std::cout<<"\t > Initialized...OK"<<std::endl;
}
/*---------------------------------------------------------------------------*/
void OpenNLFieldSolverStrategy::addSmoothingTerms() {

    Variable<double> *cot_w = m_mesh->getVariable<double,GMDS_EDGE>("cot_weight");

   /* if (m_surface_mode) {
        for (auto e:*m_surface_edges) {
            std::vector<Node> e_nodes = e.get<Node>();
            Node ni = e_nodes[0];
            Node nj = e_nodes[1];
            int i = (*m_ordering)[ni.id()];
            int j = (*m_ordering)[nj.id()];

            double c = (*cot_w)[e.id()];

            for (int k = 0; k < 9; k++) {
                nlBegin(NL_ROW);
                nlCoefficient(9 * i + k, c);
                nlCoefficient(9 * j + k, -c);
                nlRightHandSide(0);
                nlEnd(NL_ROW);
            }
        }
    } else */
        for (auto e_id:m_mesh->edges()) {
            Edge e = m_mesh->get<Edge>(e_id);
            std::vector<Node> e_nodes = e.get<Node>();
            Node ni = e_nodes[0];
            Node nj = e_nodes[1];
            int i = (*m_ordering)[ni.id()];
            int j = (*m_ordering)[nj.id()];

            double c = 1;//(*cot_w)[e.id()];

            for (int k = 0; k < 9; k++) {
                nlBegin(NL_ROW);
                nlCoefficient(9 * i + k, c);
                nlCoefficient(9 * j + k, -c);
                nlRightHandSide(0);
                nlEnd(NL_ROW);
            }
        }



}
/*----------------------------------------------------------------------------*/
