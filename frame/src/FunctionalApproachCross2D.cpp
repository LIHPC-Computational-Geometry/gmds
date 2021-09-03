/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Point.h>
#include <gmds/math/Numerics.h>
/*---------------------------------------------------------------------------*/
#include "gmds/frame/FunctionalApproachCross2D.h"
#include "gmds/frame/CrossFieldGeneration2D.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
FunctionalApproachCross2D::
FunctionalApproachCross2D(CrossFieldGeneration2D* AParent,
                          gmds::Mesh* AMesh,
                          gmds::Variable<gmds::math::Cross2D>* AField,
                          std::vector<gmds::Node>& AFromNodes,
                          std::vector<gmds::Node>& AToNodes,
                          std::vector<gmds::Face>& AFaces)
: m_parent(AParent),
m_mesh(AMesh), m_field(AField),
m_from_nodes(AFromNodes),
m_to_nodes(AToNodes)
{
    m_nodes.resize(AFromNodes.size()+AToNodes.size());
    int index=0;
    int index_to_I=0;
    
    for(unsigned int i=0;i<AFromNodes.size();i++) {
        m_nodes[index] = AFromNodes[i];
        m_id_to_local_index[AFromNodes[i].getID()]=index;
        m_id_to_I[AFromNodes[i].getID()]=index_to_I;
        index++;
        index_to_I++;
    }
    for(unsigned int i=0;i<AToNodes.size();i++) {
        m_nodes[index] = AToNodes[i];
        m_id_to_I[AToNodes[i].getID()]=index_to_I;
        index++;
        index_to_I++;
    }
    
    /* Now we build the set of edges we will work on */
    std::set<TCellID> set_edge_ids;
    for(unsigned int i=0;i<AFaces.size();i++){
        Face fi = AFaces[i];
        std::vector<TCellID> fi_edge_ids = fi.getIDs<Edge>();
        set_edge_ids.insert(fi_edge_ids.begin(),fi_edge_ids.end());
    }

    std::set<TCellID>::iterator it_e = set_edge_ids.begin();
    m_edges.resize(set_edge_ids.size());
    int i_edge=0;
    for(;it_e!=set_edge_ids.end();it_e++){
        m_edges[i_edge++] = m_mesh->get<Edge>(*it_e);
    }
}
/*---------------------------------------------------------------------------*/
void FunctionalApproachCross2D::execute()
{
    m_mark_dirichlet = m_mesh->getNewMark<Node>();
    
    for(unsigned int i=0; i<m_from_nodes.size(); i++)
        m_mesh->mark(m_from_nodes[i],m_mark_dirichlet);
    
    
    TInt nb_columns = 2*m_nodes.size();
    Eigen::VectorXd X(nb_columns);
    std::cout<<"===================="<<std::endl;
    std::cout<<" INITIALIZATION"<<std::endl;

    initialization(X);
    //==================================================================
    // STEP 3 / Assign the value onto the mesh nodes
    //==================================================================
    std::cout<<"===================="<<std::endl;
    std::cout<<"Value assignment to mesh node"<<std::endl;
    for(unsigned int i=0; i<m_nodes.size();i++) {
        Node    n_i = m_nodes[i];
        TCellID n_gmds_id = n_i.getID();
        int     n_local_id =  m_id_to_I[n_gmds_id];
        
        double  cosA = X[2*n_local_id];
        double  sinA = X[2*n_local_id+1];
        
        double  a = atan2(sinA,cosA);
        
        (*m_field)[n_gmds_id] = math::Cross2D(a);
    }
    std::cout<<"\t --> E value = "<<getFieldCurvature()<<std::endl;

    std::cout<<"===================="<<std::endl;
    std::cout<<" SMOOTHING"<<std::endl;
    smoothing(X);
    //==================================================================
    // STEP 3 / Assign the value onto the mesh nodes
    //==================================================================
    std::cout<<"===================="<<std::endl;
    std::cout<<"Value assignment to mesh node"<<std::endl;
    for(unsigned int i=0; i<m_nodes.size();i++) {
        Node    n_i = m_nodes[i];
        TCellID n_gmds_id = n_i.getID();
        int     n_local_id =  m_id_to_I[n_gmds_id];
        
        double  cosA = X[2*n_local_id];
        double  sinA = X[2*n_local_id+1];
        
        double  a = atan2(sinA,cosA);
        
        (*m_field)[n_gmds_id] = math::Cross2D(a);
    }
    std::cout<<"\t --> E value = "<<getFieldCurvature()<<std::endl;
    m_parent->writeForDebug();
    
    m_mesh->unmarkAll<Node>(m_mark_dirichlet);
    m_mesh->freeMark<Node> (m_mark_dirichlet);
}
/*---------------------------------------------------------------------------*/
void FunctionalApproachCross2D::initialization(Eigen::VectorXd& AX,
                                               const int ANbMaxIterations)
{
    // X vector is built from all the nodes with two component per nodes
    TInt nb_columns = AX.size();
    TInt nb_rows   = 2*(m_mesh->getNbEdges()+m_from_nodes.size());
    
    //==================================================================
    // STEP 1 / System assembly
    //==================================================================
    std::cout<<"System assembly"<<std::endl;
    //nb_lines*nb_columns
    Eigen::SparseMatrix<double> A(nb_rows,nb_columns);
    Eigen::VectorXd             b(nb_rows);
    Eigen::VectorXd             X0(nb_columns);
    
    Eigen::SparseMatrix<double> AT(nb_columns,nb_rows);
    
    std::cout<<"NB EDGES: "<<m_edges.size()<<std::endl;
    //we reserve a room for 2 non_zeros per columns
//    A.reserve(Eigen::VectorXi::Constant(nb_columns,2));
//    AT.reserve(Eigen::VectorXi::Constant(nb_rows,2));
//    
    double sqrtPI = sqrt(math::Constants::PI);
    double lambda = 100;
    //=================================
    //First loop on the edges
    //=================================
    for(unsigned int k=0; k < m_edges.size(); k++){
        //we fill in rows 2k and 2k+1
        
        std::vector<Node> current_nodes;
        m_edges[k].get<Node>(current_nodes);
        Node n0 = current_nodes[0];
        Node n1 = current_nodes[1];
        
        TCellID I =  m_id_to_I[n0.getID()];
        TCellID J =  m_id_to_I[n1.getID()];
        
        A.coeffRef(2*k  , 2*I  ) =  sqrtPI;
        A.coeffRef(2*k  , 2*J  ) = -sqrtPI;
        A.coeffRef(2*k+1, 2*I+1) =  sqrtPI;
        A.coeffRef(2*k+1, 2*J+1) = -sqrtPI;
        b[2*k] = 0;
        b[2*k+1] = 0;
    }
    //=================================
    //Second loop for boundary conditions
    //=================================
    int nb_lines_E = m_edges.size();
    for(unsigned int i=0; i <  m_from_nodes.size(); i++){
        int line_index = nb_lines_E+i;
        //n_i is a boundary node of the mesh
        Node n_i = m_from_nodes[i];
        int n_id = m_id_to_I[n_i.getID()];
        
        A.coeffRef(2*line_index  , 2*n_id  ) = lambda;
        A.coeffRef(2*line_index+1, 2*n_id+1) = lambda;
        
        math::Cross2D cross = (*m_field)[n_i.getID()];
        double ref_angle = cross.referenceAngle();
        b[2*line_index  ] = lambda*cos(ref_angle);
        b[2*line_index+1] = lambda*sin(ref_angle);
        
        
    }
    
    // Initial guess
    for(unsigned int i=0; i<m_nodes.size();i++) {
        Node    n_i = m_nodes[i];
        TCellID n_gmds_id = n_i.getID();
        int     n_local_id =  m_id_to_I[n_gmds_id];
        if(m_mesh->isMarked(n_i,m_mark_dirichlet)) {
            math::Cross2D cross = (*m_field)[n_gmds_id];
            double ref_angle = cross.referenceAngle();
            X0[2*n_local_id]  =cos(ref_angle);
            X0[2*n_local_id+1]=sin(ref_angle);
            
        }
        else {
            X0[2*n_local_id]=1;
            X0[2*n_local_id+1]=0;
        }
    }
    
    //    std::cout<<A<<std::endl;
    AT= A.transpose();
    A=AT*A;
    //We finalize the Sparse Matrix DS
    A.makeCompressed();
    AT.makeCompressed();
    
    //==================================================================
    // STEP 2 / System Solving
    //==================================================================
    std::cout<<"System solving for initialization"<<std::endl;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > conj_grad;
    
    conj_grad.compute(A);
  // conj_grad.setMaxIterations(ANbMaxIterations);
    AX = conj_grad.solveWithGuess(AT*b,X0);
    std::cout<<" Solved in nb. iterations= "<<conj_grad.iterations()<<", with error "<<conj_grad.error()<<std::endl;
    
    std::cout<<"\t --> E value = "<<getFieldCurvature()<<std::endl;
    
    //==================================================================
    // STEP 3 / Solution reprojection
    //==================================================================
    for(unsigned int i=0; i<m_nodes.size();i++) {
        Node    n_i = m_nodes[i];
        TCellID n_gmds_id = n_i.getID();
        int     n_id =  m_id_to_I[n_gmds_id];
        double norm = sqrt(AX[2*n_id]*AX[2*n_id] + AX[2*n_id+1]*AX[2*n_id+1]);
        AX[2*n_id]  =AX[2*n_id]/norm;
        AX[2*n_id+1]=AX[2*n_id+1]/norm;
    }
}
/*---------------------------------------------------------------------------*/
void FunctionalApproachCross2D::smoothing(Eigen::VectorXd& AX,
                                          const int AMaxNbIterations,
                                          const int AMaxInnerLoopIterations)
{
    TInt nb_columns = AX.size();
    TInt nb_rows   = 2*m_mesh->getNbEdges()+ 2*m_from_nodes.size() + m_nodes.size();
    
    double sqrtPI = sqrt(math::Constants::PI);
    double lambda = 100;
    
    
    //==================================================================
    // STEP 1 / System assembly
    //==================================================================
    std::cout<<"System assembly"<<std::endl;
    //nb_lines*nb_columns
    Eigen::SparseMatrix<double> A(nb_rows,nb_columns);
    Eigen::VectorXd             b(nb_rows);
    Eigen::VectorXd             X0(nb_columns);
    
    Eigen::SparseMatrix<double> AT(nb_columns,nb_rows);
    
    //we reserve a room for 2 non_zeros per columns
    A.reserve(Eigen::VectorXi::Constant(nb_columns,2));
    AT.reserve(Eigen::VectorXi::Constant(nb_rows,2));
    
    //=================================
    //First loop on the edges
    //=================================
    for(unsigned int k=0; k < m_edges.size(); k++){
        //we fill in rows 2k and 2k+1
        
        std::vector<Node> current_nodes;
        m_edges[k].get<Node>(current_nodes);
        Node n0 = current_nodes[0];
        Node n1 = current_nodes[1];
        
        TCellID I =  m_id_to_I[n0.getID()];
        TCellID J =  m_id_to_I[n1.getID()];
        
        A.coeffRef(2*k  , 2*I  ) =  sqrtPI;
        A.coeffRef(2*k  , 2*J  ) = -sqrtPI;
        A.coeffRef(2*k+1, 2*I+1) =  sqrtPI;
        A.coeffRef(2*k+1, 2*J+1) = -sqrtPI;
        b[2*k] = 0;
        b[2*k+1] = 0;
    }
    //=================================
    //Second loop for boundary conditions
    //=================================
    int nb_lines_E = m_edges.size();
    for(unsigned int i=0; i <  m_from_nodes.size(); i++){
        int line_index = nb_lines_E+i;
        //n_i is a boundary node of the mesh
        Node n_i = m_from_nodes[i];
        int n_id = m_id_to_I[n_i.getID()];
        
        A.coeffRef(2*line_index  , 2*n_id  ) = lambda;
        A.coeffRef(2*line_index+1, 2*n_id+1) = lambda;
        
        math::Cross2D cross = (*m_field)[n_i.getID()];
        double ref_angle = cross.referenceAngle();
        b[2*line_index  ] = lambda*cos(ref_angle);
        b[2*line_index+1] = lambda*sin(ref_angle);
        
        
    }
    std::cout<<"Start iteration"<<std::endl;
    for (int i=0; i<AMaxNbIterations; i++) {
        X0 = AX;
        
        //=================================
        //Third loop for feasibility constraint
        //=================================
        nb_lines_E  = 2*(m_edges.size() + m_from_nodes.size());
        for(unsigned int i=0; i <  m_nodes.size(); i++){
            int line_index = nb_lines_E+i;
            //n_i is a  node of the mesh
            Node n_i = m_nodes[i];
            int n_id = m_id_to_I[n_i.getID()];
            
            A.coeffRef(line_index, 2*n_id  ) = lambda*X0[2*n_id];
            A.coeffRef(line_index, 2*n_id+1) = lambda*X0[2*n_id+1];
            
            b[line_index] = lambda;
            
        }
        
        //==================================================================
        // STEP 2 / System Solving
        //==================================================================
     //   std::cout<<"System solving for smoothing"<<std::endl;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > conj_grad;
        
        AT= A.transpose();
        Eigen::SparseMatrix<double> ATA(nb_columns,nb_columns);

        ATA=AT*A;
        //We finalize the Sparse Matrix DS
        ATA.makeCompressed();
        
        conj_grad.compute(ATA);
        conj_grad.setMaxIterations(AMaxInnerLoopIterations);
        AX = conj_grad.solveWithGuess(AT*b,X0);
        
        //==================================================================
        // STEP 3 / Solution reprojection
        //==================================================================
        for(unsigned int i=0; i<m_nodes.size();i++) {
            Node    n_i = m_nodes[i];
            TCellID n_gmds_id = n_i.getID();
            int     n_id =  m_id_to_I[n_gmds_id];
            double norm = sqrt(AX[2*n_id]*AX[2*n_id] + AX[2*n_id+1]*AX[2*n_id+1]);
            AX[2*n_id]  =AX[2*n_id]/norm;
            AX[2*n_id+1]=AX[2*n_id+1]/norm;
        }


    }
}
/*---------------------------------------------------------------------------*/
double FunctionalApproachCross2D::getFieldCurvature(){
    double curvature=0;
    for(unsigned int k=0; k < m_edges.size(); k++){
        //we fill in rows 2k and 2k+1
        
        std::vector<Node> current_nodes;
        m_edges[k].get<Node>(current_nodes);
        Node n0 = current_nodes[0];
        Node n1 = current_nodes[1];

        math::Cross2D cross0 = (*m_field)[n0.getID()];
        math::Cross2D cross1 = (*m_field)[n1.getID()];
        double rot_angle = cross0.angle(cross1);
        curvature += rot_angle*rot_angle;
        
    }

    return curvature;
}
/*---------------------------------------------------------------------------*/