/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Numerics.h>
#include <gmds/math/FE.h>
#include <gmds/math/Point.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Log.h>
#include <gmds/math/SHarmonicL4.h>
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include <gmds/frame3d/FEMSolverStrategy.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
FEMSolverStrategy::FEMSolverStrategy(){}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::solve() {

    //We build the system AX=b in this function. It gathers 9 independant
    // sytem to solve. One Laplace equation per SH4 component
    int nb_unknowns[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (auto i = 0; i < 9; i++) {
        m_dirichlet_mark[i] = m_mesh->newMark<Node>();
    }

    for (auto i = 0; i < 9; i++) {
        //we get the fixed and free nodes, i.e. nodes holding dirichlet
        // conditions and other nodes
        for (auto n_id : m_mesh->nodes()) {
            Node n = m_mesh->get<Node>(n_id);
            math::SHarmonicL4 shi = (*m_harmonic_field)[n.id()];
            if (m_mesh->isMarked(n, m_markNodeOnSurf)) {
                //fixed frame, so all the 9 components are boundary condition
                m_mesh->mark(n, m_dirichlet_mark[i]);
            }
//            else if (m_mesh->isMarked(n, m_markNodeOnSurf)) {
//                //only 2 components are fixed, the others are unknowns
//                bool fixed = false;
//                //We check if this coordinate is known or not
//
//
//
//                TCellID global_id = n.id();
//
//                math::SHarmonicL4 h4 = (*m_H4)[global_id];
//                math::SHarmonicL4 h8 = (*m_H8)[global_id];
//                //only one component of h4 and one of h8 is non-null
//                // it gives a boundary conditions to fullfill
//                int non_null_4 = -1;
//                int non_null_8 = -1;
//                for(auto k=0;k<9;k++){
//                    if(h4[k]!=0)
//                        non_null_4=k;
//                    if(h8[k]!=0)
//                        non_null_8=k;
//                }
//                if(non_null_4==i ||non_null_8==i)
//                    fixed = true;
//fixed = false;
//                if (fixed) {
//                    m_mesh->mark(n, m_dirichlet_mark[i]);
//                } else {
//                    (*m_ordering)[global_id]=nb_unknowns[i];
//                    nb_unknowns[i]++;
//                }
//            }
            else {
                (*m_ordering)[n.id()]=nb_unknowns[i];
                nb_unknowns[i]++;
            }

        }

    }
    // We begin by initialize the Sparse Matrix m_A
    TInt max_adj = 0;
    for (auto n_id : m_mesh->nodes()) {
        Node n = m_mesh->get<Node>(n_id);
        std::vector<Edge> edges_i = n.get<Edge>();
        if (edges_i.size() > max_adj)
            max_adj = edges_i.size();
    }

    //We build the 9 systems AX=b for components 0 to 9.
    for (auto i = 0; i < 9; i++) {
        Eigen::VectorXd Xi(nb_unknowns[i]);
        Eigen::VectorXd bi(nb_unknowns[i]);
        Eigen::SparseMatrix<double> Ai(nb_unknowns[i], nb_unknowns[i]);

        Ai.reserve(Eigen::VectorXi::Constant(nb_unknowns[i], max_adj));
        buildStiffnessMatrix(Ai, i);
        Ai.makeCompressed();

        buildRightHandDirichlet(bi, i);

        std::cout << "System solving for component " << i << std::endl;
        Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > cholcos(Ai);
        Xi = cholcos.solve(bi);
        m_X[i] = Xi;
        std::cout<<Xi<<std::endl;
    }


    for (auto i = 0; i < 9; i++) {
        m_mesh->unmarkAll<Node>(m_dirichlet_mark[i]);
        m_mesh->freeMark<Node>(m_dirichlet_mark[i]);
    }

}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::initializeAssembly() {

}
/*---------------------------------------------------------------------------*/
void
FEMSolverStrategy::buildStiffnessMatrix(Eigen::SparseMatrix<double>& AM,
                                        const int AI)
{
    //==================================================================
    // We compute the global stiffness
    //==================================================================
    // First loop on the tets
    for (auto r_id: m_mesh->regions()) {
        Region r = m_mesh->get<Region>(r_id);
        std::vector<Node> current_nodes;
        r.get<Node>(current_nodes);
        math::Point p1 = current_nodes[0].point();
        math::Point p2 = current_nodes[1].point();
        math::Point p3 = current_nodes[2].point();
        math::Point p4 = current_nodes[3].point();
        math::Matrix<4, 4, double> s = math::TetrahedronP1::stiffnessMatrix(p1,p2,p3,p4);
        // Second loop on the internal  nodes
        // with i,j their index locally to the face, and I,J the index global
        // to the mesh
        for (auto i = 0; i < current_nodes.size(); i++) {
            Node n_i = current_nodes[i];
            if (m_mesh->isMarked(n_i, m_dirichlet_mark[AI]))
                continue;
            TCellID I =(*m_ordering)[n_i.id()];

            for (unsigned int j = 0; j < current_nodes.size(); j++) {
                Node n_j = current_nodes[j];
                if (m_mesh->isMarked(n_j, m_dirichlet_mark[AI]))
                    continue;
                TCellID J = (*m_ordering)[n_j.id()];

                AM.coeffRef(I, J) += s.get(i, j);
            }
        }
    }
}
/*---------------------------------------------------------------------------*/
void
FEMSolverStrategy::buildRightHandDirichlet(Eigen::VectorXd& AB, const int AI)
{
    for (auto n_id : m_mesh->nodes()) {
        Node n = m_mesh->get<Node>(n_id);
        math::SHarmonicL4 shi = (*m_harmonic_field)[n.id()];
        if (!m_mesh->isMarked(n, m_dirichlet_mark[AI])) {
            //se we have an inner node
            TCellID I = (*m_ordering)[n.id()];
            AB[I]=0;
            std::vector<Region> vicinity = n.get<Region>();

            double null_vec[4]={0,0,0,0};
            for (auto r:vicinity) {
                std::vector<Node> nodes_of_r = r.get<Node>();
                bool is_boundary_region = false;

                math::Vector4d Ud(null_vec);
                int n_index_in_r = -1;
                // we store if the i^th node of r is on the boundary or not

                for (auto k = 0; k < 4; k++) {
                    Node n_k = nodes_of_r[k];
                    if (m_mesh->isMarked(n_k, m_dirichlet_mark[AI])) {
                        is_boundary_region = true;
                        math::SHarmonicL4 shi = (*m_harmonic_field)[n.id()];
                        Ud[k] = shi[AI];
                    }
                    if (n_k.id() == n.id()) {
                        n_index_in_r = k;
                    }
                }
                if (is_boundary_region) {
                    math::Point pi = nodes_of_r[0].point();
                    math::Point pj = nodes_of_r[1].point();
                    math::Point pk = nodes_of_r[2].point();
                    math::Point pl = nodes_of_r[3].point();
                    math::Matrix<4, 4, double> s = math::TetrahedronP1::stiffnessMatrix(pi, pj, pk,pl);

                    // Ud is build
                    double val_tab[4]={s.get(0, n_index_in_r),
                                       s.get(1, n_index_in_r),
                                       s.get(2, n_index_in_r),
                                       s.get(3, n_index_in_r)};

                    math::Vector4d v(val_tab);
                    AB[I] += -v.dot(Ud);
                }
            }

        }
    }
}

/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::finalizeAssembly(){

}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::clean()
{
}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::setX()
{
}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::getFeasibleSolution()
{
    for (auto n_id:m_mesh->nodes()) {

        if (!m_mesh->isMarked<Node>(n_id, m_markNodeLocked)) {

            int local_id = (*m_ordering)[n_id];
            double t[9] = {
                    m_X[0][local_id],
                    m_X[1][local_id],
                    m_X[2][local_id],
                    m_X[3][local_id],
                    m_X[4][local_id],
                    m_X[5][local_id],
                    m_X[6][local_id],
                    m_X[7][local_id],
                    m_X[8][local_id]};

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
            (*m_chart_field)[n_id] = math::Chart(math::Vector3d(vx[0], vx[1], vx[2]),
                                                   math::Vector3d(vy[0], vy[1], vy[2]),
                                                   math::Vector3d(vz[0], vz[1], vz[2]));

        }
    }

}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::addLocalOptimConstraints(){

}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::addLockedTerms() {
   }
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::addBoundaryConstraints(){
}
/*---------------------------------------------------------------------------*/
void FEMSolverStrategy::addSmoothingTerms() {
}
/*----------------------------------------------------------------------------*/