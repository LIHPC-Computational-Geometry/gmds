/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Plane.h>
/*---------------------------------------------------------------------------*/
// HLBFGS Header files
#include "HLBFGS.h"
/*---------------------------------------------------------------------------*/
// FRAME File Headers
#include "FieldSmoother.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
IGMesh* FieldSmoother::m_mesh = 0;
int FieldSmoother::m_mark_node_on_point = 0;
int FieldSmoother::m_mark_node_on_curve = 0;
int FieldSmoother::m_mark_node_on_surface= 0;
std::vector<gmds::Node> FieldSmoother::m_candidates;
std::vector<gmds::Edge> FieldSmoother::m_edges;
std::map<gmds::TCellID, int> FieldSmoother::m_candidates_index;
std::map<gmds::TCellID, gmds::math::Chart> FieldSmoother::m_boundary_triad;
/*---------------------------------------------------------------------------*/
FieldSmoother::
FieldSmoother(IGMesh* AMesh,
               gmds::Variable<gmds::math::AxisAngleRotation>* ARotField,
               std::map<gmds::TCellID, gmds::math::Vector3d>& ANormal,
               const int AMarkNodeOnPoint,
               const int AMarkNodeOnCurve,
               const int AMarkNodeOnSurface)
:m_rotation_field(ARotField),
m_normal(ANormal), m_mark_candidates(0),
m_candidate_mark_initialized(false)
{
    m_mesh = AMesh;
    m_mark_node_on_point  = AMarkNodeOnPoint;
    m_mark_node_on_curve  = AMarkNodeOnCurve;
    m_mark_node_on_surface= AMarkNodeOnSurface;
}
/*---------------------------------------------------------------------------*/
void FieldSmoother::selectNodes(const int AMark)
{
    m_mark_candidates = AMark;
    m_candidate_mark_initialized = true;
}
/*----------------------------------------------------------------------------*/
void FieldSmoother::evalF(int AN,
                          double* AX, double *APrevX,
                          double* AFunc, double* AGrad)
{
    //===========================================================
    // Initialization of local array and output values
    //===========================================================
    int nb_vars = AN;
    
    
    // sin and cosinus value of Euler angles are locally stored to avoid
    // further computations.
    std::vector<double> cos_val, sin_val;
    cos_val.resize(nb_vars);
    sin_val.resize(nb_vars);
    
    for (int i = 0; i < nb_vars; i++){
        cos_val[i] = cos(AX[i]);
        sin_val[i] = sin(AX[i]);
    }
    
    *AFunc = 0.0;
    
    // AGrad is a double N-size array already initialized with the right
    // size in HLBFGS. We initialize all its entries to null
    for (unsigned int i = 0;i < nb_vars;i++)
        AGrad[i] = 0.0;
    
    
    //===========================================================
    // Matrix initalization (system build on edges)
    //===========================================================
    for (unsigned int i = 0; i < m_edges.size(); i++) {
        std::vector<Node> edge_nodes = m_edges[i].get<Node>();
        Node n0 = edge_nodes[0];
        Node n1 = edge_nodes[1];
        
        int i0 = m_candidates_index[n0.getID()];
        int i1 = m_candidates_index[n1.getID()];
        
        gmds::math::Matrix<3, 3, double> Mix, Miy, Miz;
        gmds::math::Matrix<3, 3, double> Mjx, Mjy, Mjz;
        gmds::math::Matrix<3, 3, double> DMix, DMiy, DMiz;
        gmds::math::Matrix<3, 3, double> DMjx, DMjy, DMjz;
        gmds::math::Matrix<3, 3, double> M, Mi, Mj, Mi1, Mi2, Mi3;
        // Matrix constructor must ensure that all the 3x3 matrices
        // previously defined contain only  null values
        
        // Matrix computation for node n0
        initMatrix(n0, cos_val, sin_val,
                   Mix, Miy, Miz,
                   DMix, DMiy, DMiz);
        
        Mi = Mix * Miy * Miz;
        Mi1 = DMix * Miy * Miz;
        Mi2 = Mix * DMiy * Miz;
        Mi3 = Mix * Miy * DMiz;
        
        // Matrix computation for node n1
        initMatrix(n1, cos_val, sin_val,
                   Mjx, Mjy, Mjz,
                   DMjx, DMjy, DMjz);
        
        Mj = Mjx * Mjy * Mjz;
        M = Mi.transpose() * Mj;
        
        double diff0[3][3][3] = {
            { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
            { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
            { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } };
        double diff1[3][3][3] = {
            { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
            { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
            { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } };
        double diff_i[3] = { 0, 0, 0 };
        double diff_j[3] = { 0, 0, 0 };
        double e[18];

        gmds::math::Matrix<3, 3, double> tmp_mat = (Mi.transpose() *
                                                    (DMjx * Mjy * Mjz));
        tmp_mat.getTab(diff1[0]);
        
        tmp_mat = (Mi.transpose() * (Mjx * DMjy * Mjz));
        tmp_mat.getTab(diff1[1]);
        
        tmp_mat = (Mi.transpose() * (Mjx * Mjy * DMjz));
        tmp_mat.getTab(diff1[2]);
        
        (Mi1.transpose() * Mj).getTab(diff0[0]);
        (Mi2.transpose() * Mj).getTab(diff0[1]);
        (Mi3.transpose() * Mj).getTab(diff0[2]);
        
        e[0]  = M(0, 0) * M(0, 1);
        e[1]  = M(0, 0) * M(0, 2);
        e[2]  = M(0, 1) * M(0, 2);
        e[3]  = M(1, 0) * M(1, 1);
        e[4]  = M(1, 0) * M(1, 2);
        e[5]  = M(1, 1) * M(1, 2);
        e[6]  = M(2, 0) * M(2, 1);
        e[7]  = M(2, 0) * M(2, 2);
        e[8]  = M(2, 1) * M(2, 2);
        e[9]  = M(0, 0) * M(1, 0);
        e[10] = M(0, 0) * M(2, 0);
        e[11] = M(1, 0) * M(2, 0);
        e[12] = M(0, 1) * M(1, 1);
        e[13] = M(0, 1) * M(2, 1);
        e[14] = M(1, 1) * M(2, 1);
        e[15] = M(0, 2) * M(1, 2);
        e[16] = M(0, 2) * M(2, 2);
        e[17] = M(1, 2) * M(2, 2);
        
        double m_sigma = 1.0;
        double local_f = 0.0;
        
        static const int eflag[9][2][2] = {
            { { 0, 0 }, { 0, 1 } }, { { 0, 0 }, { 0, 2 } }, { { 0, 1 }, { 0, 2 } },
            { { 1, 0 }, { 1, 1 } }, { { 1, 0 }, { 1, 2 } }, { { 1, 1 }, { 1, 2 } },
            { { 2, 0 }, { 2, 1 } }, { { 2, 0 }, { 2, 2 } }, { { 2, 1 }, { 2, 2 } } };
        
        
        for (int k = 0; k < 9; k++) {
            double v_k = exp(e[k + 9] * e[k + 9] / m_sigma);
            
            local_f += 2 * (v_k - 1);
            
            for (int l = 0; l < 3; l++){
                
                double d_i =((diff0[l][eflag[k][0][1]][eflag[k][0][0]] *
                              M(eflag[k][1][1], eflag[k][1][0])) +
                             (diff0[l][eflag[k][1][1]][eflag[k][1][0]] *
                              M(eflag[k][0][1], eflag[k][0][0])));

                double d_j =((diff1[l][eflag[k][0][1]][eflag[k][0][0]] *
                              M(eflag[k][1][1], eflag[k][1][0])) +
                             (diff1[l][eflag[k][1][1]][eflag[k][1][0]] *
                              M(eflag[k][0][1], eflag[k][0][0])));
                
                diff_i[l] += 2 * (e[k + 9] * d_i * v_k / m_sigma);
                diff_j[l] += 2 * (e[k + 9] * d_j * v_k / m_sigma);
            }
        }
        
        
        *AFunc += (local_f / 2.0);
        for (int k = 0; k < 3; k++)
        {
            AGrad[3 * i0 + k] += diff_i[k];
            AGrad[3 * i1 + k] += diff_j[k];
        }
        
    }//for (unsigned int i = 0; i < m_edges.size(); i++)
}
/*---------------------------------------------------------------------------*/
void FieldSmoother::initMatrix(Node& ANode,
                               std::vector<TCoord>& ACos,
                               std::vector<TCoord>& ASin,
                               gmds::math::Matrix<3, 3, double>& AMix,
                               gmds::math::Matrix<3, 3, double>& AMiy,
                               gmds::math::Matrix<3, 3, double>& AMiz,
                               gmds::math::Matrix<3, 3, double>& ADMix,
                               gmds::math::Matrix<3, 3, double>& ADMiy,
                               gmds::math::Matrix<3, 3, double>& ADMiz)
{
    Node n = ANode;
    int i = m_candidates_index[n.getID()];
    
    if (m_mesh->isMarked(n, m_mark_node_on_curve)) {
        //===========================================================
        //n is classified on a geometric curve
        //===========================================================
        AMix(0, 0) =  1.;
        AMix(1, 1) =  ACos[3*i];
        AMix(2, 2) =  ACos[3*i];
        AMix(1, 2) = -ASin[3*i];
        AMix(2, 1) =  ASin[3*i];
        
        AMiy(1, 1) =  1.;
        AMiy(0, 0) =  ACos[3*i + 1];
        AMiy(2, 2) =  ACos[3*i + 1];
        AMiy(0, 2) =  ASin[3*i + 1];
        AMiy(2, 0) = -ASin[3*i + 1];
        
        AMiz(2, 2) =  1.;
        AMiz(0, 0) =  ACos[3*i + 2];
        AMiz(1, 1) =  ACos[3*i + 2];
        AMiz(0, 1) = -ASin[3*i + 2];
        AMiz(1, 0) =  ASin[3*i + 2];
        
        // Stationnary values for the scheme so null derivatices in this node
    }
    else if (m_mesh->isMarked(n, m_mark_node_on_surface)){
        //===========================================================
        // n is classifed on a geometric surface
        //===========================================================
        AMix(0, 0) = 1.;
        AMix(1, 1) = 1.;
        AMix(2, 2) = 1.;
        
        AMiy(0, 0) = 1.;
        AMiy(1, 1) = 1.;
        AMiy(2, 2) = 1.;
        
        
        math::Chart t = m_boundary_triad[n.getID()];
        math::Vector tX = t.X();
        math::Vector tY = t.Y();
        math::Vector tZ = t.Z();
        
        AMiz(0, 0) =  ACos[3*i + 2] * tY[0] + ASin[3*i + 2] * tZ[0];
        AMiz(1, 0) =  ACos[3*i + 2] * tY[1] + ASin[3*i + 2] * tZ[1];
        AMiz(2, 0) =  ACos[3*i + 2] * tY[2] + ASin[3*i + 2] * tZ[2];
        AMiz(0, 1) = -ASin[3*i + 2] * tY[0] + ACos[3*i + 2] * tZ[0];
        AMiz(1, 1) = -ASin[3*i + 2] * tY[1] + ACos[3*i + 2] * tZ[1];
        AMiz(2, 1) = -ASin[3*i + 2] * tY[2] + ACos[3*i + 2] * tZ[2];
        AMiz(0, 2) = tX[0];
        AMiz(1, 2) = tX[1];
        AMiz(2, 2) = tX[2];
        
        ADMiz(0, 0) = -ASin[3*i + 2] * tY[0] + ACos[3*i + 2] * tZ[0];
        ADMiz(1, 0) = -ASin[3*i + 2] * tY[1] + ACos[3*i + 2] * tZ[1];
        ADMiz(2, 0) = -ASin[3*i + 2] * tY[2] + ACos[3*i + 2] * tZ[2];
        ADMiz(0, 1) = -ACos[3*i + 2] * tY[0] - ASin[3*i + 2] * tZ[0];
        ADMiz(1, 1) = -ACos[3*i + 2] * tY[1] - ASin[3*i + 2] * tZ[1];
        ADMiz(2, 1) = -ACos[3*i + 2] * tY[2] - ASin[3*i + 2] * tZ[2];
        ADMiz(0, 2) = 0.;
        ADMiz(1, 2) = 0.;
        ADMiz(2, 2) = 0.;
        
    }
    else{
        //===========================================================
        // General case inside the volume
        //===========================================================
        //point pas sur le bord, remplissage classique des matrices
        //			std::cout<<"vi inside"<<std::endl;
        AMix(0, 0) =  1.;
        AMix(1, 1) =  ACos[3*i];
        AMix(2, 2) =  ACos[3*i];
        AMix(1, 2) = -ASin[3*i];
        AMix(2, 1) =  ASin[3*i];
        
        AMiy(1, 1) =  1.;
        AMiy(0, 0) =  ACos[3*i + 1];
        AMiy(2, 2) =  ACos[3*i + 1];
        AMiy(0, 2) =  ASin[3*i + 1];
        AMiy(2, 0) = -ASin[3*i + 1];
        
        AMiz(2, 2) =  1.;
        AMiz(0, 0) =  ACos[3*i + 2];
        AMiz(1, 1) =  ACos[3*i + 2];
        AMiz(0, 1) = -ASin[3*i + 2];
        AMiz(1, 0) =  ASin[3*i + 2];
        
        ADMix(1, 1) = -ASin[3*i];
        ADMix(2, 2) = -ASin[3*i];
        ADMix(1, 2) = -ACos[3*i];
        ADMix(2, 1) =  ACos[3*i];
        
        ADMiy(0, 0) = -ASin[3*i + 1];
        ADMiy(2, 2) = -ASin[3*i + 1];
        ADMiy(0, 2) =  ACos[3*i + 1];
        ADMiy(2, 0) = -ACos[3*i + 1];
        
        
        ADMiz(0, 0) = -ASin[3*i + 2];
        ADMiz(1, 1) = -ASin[3*i + 2];
        ADMiz(0, 1) = -ACos[3*i + 2];
        ADMiz(1, 0) =  ACos[3*i + 2];
    }
}
/*----------------------------------------------------------------------------*/
void FieldSmoother::
newiteration(int ANbIter, int AIter, double* AX,
             double* AFunc, double* AGrad, double* AGradNorm)
{
    std::cout << ANbIter << ": " << AIter << " " << *AFunc << " "
    << *AGradNorm << std::endl;
}
/*----------------------------------------------------------------------------*/
void FieldSmoother:: rebuildAxisAngleRotations(double*& AEuler)
{
    
    for (unsigned int i = 0; i < m_candidates.size(); i++){
        
        Node n = m_candidates[i];
        TCoord angleX = AEuler[3 * i];
        TCoord angleY = AEuler[3 * i + 1];
        TCoord angleZ = AEuler[3 * i + 2];
        
        math::Quaternion q;
        if (m_mesh->isMarked(n,m_mark_node_on_curve)) {
            //==============================================
            // CASE 1 - Classified on a curve
            //==============================================
            q.setFromEulerAngle(angleX, angleY, angleZ);
        }
        else if (m_mesh->isMarked(n, m_mark_node_on_surface)) {
            //==============================================
            // CASE 2 - Classified on a surface
            //==============================================
            math::Chart t = m_boundary_triad[n.getID()];
            
            math::Vector new_y = cos(angleZ)*t.Y();
            math::Vector new_z = sin(angleZ)*t.Z();
            math::Vector v1 = new_y+new_z;
            
            math::Vector v2 = t.X().cross(v1);
            
            math::Quaternion q(math::Chart(t.X(), v1, v2));

            // the quaternion must be aligned with the normal
            // this last projection step ensures it
            math::Vector3d normal = m_normal[n.getID()];
            q=q.alignWith(math::Vector(normal.X(),
                                       normal.Y(),
                                       normal.Z()));

        }
        else {
            //==============================================
            // CASE 3 - Classified into a volume
            //==============================================
            q.setFromEulerAngle(angleX, angleY, angleZ);
        }
        

        //Update of the variable at node n
        math::AxisAngleRotation old_v = (*m_rotation_field)[n.getID()];
        (*m_rotation_field)[n.getID()] = math::AxisAngleRotation(q);
        
        std::cout<<old_v.axis()<<" --> "
        <<(*m_rotation_field)[n.getID()].axis()<<std::endl;
        
    }
}

/*---------------------------------------------------------------------------*/
void FieldSmoother::execute()
{
    //=====================================================================
    // Static variables are initialized, i.e. the nodes we will work on
    //=====================================================================
    if (m_candidate_mark_initialized) {
        std::cout<<"HLBGFS smoothing performed on marked mesh nodes nodes!"
        <<std::endl;
    }
    else {
        std::cout<<"HLBGFS smoothing performed on all the mesh nodes!"
        <<std::endl;
    }

    //Now, we collect the nodes we want to smooth the Quaternion on
    initCandidates();
    
    //=====================================================================
    // Now the smoothing can be performed
    //=====================================================================
    double* euler_angles = new double[3 * m_candidates.size()];

    buildEulerAngles(euler_angles);
    
    double parameter[20];
    int info[20];
    // nb variables = 3 Euler angle for each vertex we work on
    int N = 3 * m_candidates.size();
    // Choice of the optimisation algorithm
    int M = 7;
    int T = 0;
    int max_nb_iter = 1000;
    bool with_hessian = false;

    //=====================================================================
    // Initialization of HLBFGS library parameters
    //=====================================================================
    // We use default values provided by HLBFGS tutorial webpage
    INIT_HLBFGS(parameter, info);
    info[0] = 20;
    info[4] = max_nb_iter;
    info[6] = T;
    info[7] = with_hessian ? 1 : 0;
    info[10] = 0;
    info[11] = 1;
    
    //=====================================================================
    // Call to the HLBFGS solver
    //=====================================================================
    HLBFGS(N, M, &euler_angles[0], (this->evalF),
           0, HLBFGS_UPDATE_Hessian,
           this->newiteration, parameter, info);
    //=====================================================================
    // Rebuilt of the solution
    //=====================================================================
    rebuildAxisAngleRotations(euler_angles);
    
    //Local memory cleaning
    delete[] euler_angles;
    
    writeSolution();
}
/*---------------------------------------------------------------------------*/
void FieldSmoother::initCandidates()
{
    
    //=====================================================================
    // NODES INITIALIZATION
    //=====================================================================
    m_candidates.clear();
    IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
    if(m_candidate_mark_initialized) {
        //================================================================
        // Case 1 - Only marked nodes are taken into account
        //================================================================
        for (; !it_nodes.isDone(); it_nodes.next()) {
            Node n = it_nodes.value();
            if (m_mesh->isMarked(n,m_mark_candidates))
                m_candidates.push_back(it_nodes.value());
        }
    }
    else {
        //================================================================
        // Case 2 - All nodes are taken into account
        //================================================================
        m_candidates.resize(m_mesh->getNbNodes());
        int index=0;
        for (; !it_nodes.isDone(); it_nodes.next()) {
            m_candidates[index++]=it_nodes.value();
        }
    }
    
    std::cout << std::endl << "HLBGFS smoothing: Nb vertex candidates ("
    << m_candidates.size() << " / " << m_mesh->getNbNodes() <<")"<< std::endl;
    
    //=====================================================================
    // EDGES INITIALIZATION
    //=====================================================================
    
    //Definition of the support Edges
    m_edges.clear();
    if(m_candidate_mark_initialized) {
        //================================================================
        // Case 1 - Only marked nodes are taken into account
        //================================================================
        for (unsigned int i = 0; i < m_candidates.size(); i++){
            Node n = m_candidates[i];
            std::vector<Edge> adj_edges = n.get<Edge>();
            
            for (unsigned int j = 0; j < adj_edges.size(); j++)
            {
                Edge e_j = adj_edges[j];
                
                std::vector<Node> e_j_nodes = e_j.get<Node>();
                Node other_node =
                (e_j_nodes[0].getID() == n.getID()) ? e_j_nodes[1] : e_j_nodes[0];
                
                if (m_mesh->isMarked(other_node, m_mark_candidates))
                {
                    if (other_node.getID() < n.getID()) //to store the edge only once
                        m_edges.push_back(e_j);
                }
            }
        }
    }
    else {
        //================================================================
        // Case 2 - All nodes are taken into account
        //================================================================
        m_edges.resize(m_mesh->getNbEdges());
        IGMesh::edge_iterator it_edges = m_mesh->edges_begin();

        int index=0;
        for (; !it_edges.isDone(); it_edges.next()) {
            m_edges[index++]=it_edges.value();
        }
    }
    
    
    std::cout << std::endl << "HLBGFS smoothing: Nb edge candidates ("
    << m_edges.size() << " / " << m_mesh->getNbEdges() <<")"<< std::endl;

}
/*---------------------------------------------------------------------------*/
void FieldSmoother::buildEulerAngles(double*& AEuler)
{
    //===============================================================
    // We compute Euler angles for any candidate vertex
    //===============================================================
    for (unsigned int i = 0; i < m_candidates.size(); i++) {
        
        Node n = m_candidates[i];
        TCellID n_id = n.getID();
        math::AxisAngleRotation r= (*m_rotation_field)[n_id];
        math::Quaternion q = r.quaternion();
        
        //we store the candidate index of n_i
        m_candidates_index[n.getID()] = i;
        
        if (m_mesh->isMarked(n, m_mark_node_on_curve)){
            //===================================================
            //Case 1 - node on a curve
            //===================================================
            q.toEulerAngle(AEuler[3 * i],
                           AEuler[3 * i + 1],
                           AEuler[3 * i + 2]);
        } //if (m_mesh->isMarked(n, m_mark_node_on_curve))
        else if (m_mesh->isMarked(n, m_mark_node_on_surface)) {
            //===================================================
            //Case 2 - node on a surface
            //===================================================
            math::Chart t(q);
            math::Vector X[6];
            X[0] = t.X();
            X[1] = t.Y();
            X[2] = t.Z();
            X[3] = X[0].opp();
            X[4] = X[1].opp();
            X[5] = X[2].opp();
            
            math::Vector3d normal_3D = m_normal[n_id];
            math::Vector normal(normal_3D.X(),
                                normal_3D.Y(),
                                normal_3D.Z());
            
            // We get among X, the vector that is the most aligned with
            // the vertex normal
            int index_normal = 0;
            double max_dot = normal.dot(X[0]);
            for (unsigned int j = 1; j<6; j++){
                TCoord current_dot = normal.dot(X[j]);
                if (current_dot>max_dot){
                    index_normal = j;
                    max_dot = current_dot;
                }
            }
            
            //Outward chart is so built
            math::Vector X_out[3];
            X_out[0] = X[index_normal];
            if (index_normal == 0){
                X_out[1] = X[1];
                X_out[2] = X[2];
            }
            else if (index_normal == 1){
                X_out[1] = X[2];
                X_out[2] = X[0];
            }
            else if (index_normal == 2){
                X_out[1] = X[0];
                X_out[2] = X[1];
            }
            else if (index_normal == 3){
                X_out[1] = X[2];
                X_out[2] = X[1];
            }
            else if (index_normal == 4){
                X_out[1] = X[0];
                X_out[2] = X[2];
            }
            else if (index_normal == 5){
                X_out[1] = X[0];
                X_out[2] = X[4];
            }
            
            m_boundary_triad[n_id] = math::Chart(X_out[0], X_out[1], X_out[2]);
            
            math::Quaternion q_out(m_boundary_triad[n_id]);
            
            q_out.toEulerAngle(AEuler[3 * i],
                               AEuler[3 * i + 1],
                               AEuler[3 * i + 2]);
        } //else if (m_mesh->isMarked(n, m_mark_node_on_surface))
        else {
            //===================================================
            //Case 3 - node inside the volume
            //===================================================
            q.toEulerAngle(AEuler[3 * i],
                           AEuler[3 * i + 1],
                           AEuler[3 * i + 2]);
        }
    }
}
/*---------------------------------------------------------------------------*/
void FieldSmoother::writeSolution(){
    static int nb_file = 0;
    
    IGMesh::node_iterator it = m_mesh->nodes_begin();
    double bound = 100000;
    double x_min = bound;
    double y_min = bound;
    double z_min = bound;
    double x_max = -bound;
    double y_max = -bound;
    double z_max = -bound;
    for (; !it.isDone(); it.next())
    {
        Node n = it.value();
        math::Point p = n.getPoint();
        if (p.X() < x_min)
            x_min = p.X();
        if (p.X() > x_max)
            x_max = p.X();
        
        if (p.Y() < y_min)
            y_min = p.Y();
        if (p.Y() > y_max)
            y_max = p.Y();
        
        if (p.Z() < z_min)
            z_min = p.Z();
        if (p.Z() > z_max)
            z_max = p.Z();
    }
    double dist_x = x_max - x_min;
    double dist_y = y_max - y_min;
    double dist_z = z_max - z_min;
    
    double cube_size = 0;
    if (dist_x <= dist_y && dist_x <= dist_z){
        cube_size = dist_x;
    }
    else if (dist_y <= dist_x && dist_y <= dist_z){
        cube_size = dist_y;
    }
    else
        cube_size = dist_z;
    
    //    cube_size /= 20;
    
    MeshModel model_cube(DIM3 | R | N | R2N);
    IGMesh mesh_cube(model_cube);
    
    Variable<int>* v = mesh_cube.newVariable<int>(GMDS_REGION, "classification");
    
    
    for (it = m_mesh->nodes_begin(); !it.isDone(); it.next())
    {
        Node n = it.value();
        math::Point center = n.getPoint();
        
        math::AxisAngleRotation rot =(*m_rotation_field)[n.getID()];
        math::Chart c = rot.toChart();
        
        math::Vector3d evx(c.X()[0],c.X()[1],c.X()[2]);
        math::Vector3d evy(c.Y()[0],c.Y()[1],c.Y()[2]);
        math::Vector3d evz(c.Z()[0],c.Z()[1],c.Z()[2]);
        
        evx.normalize();
        evy.normalize();
        evz.normalize();
        
        if(evx.X()>100){
            std::cout<<"Error Nan for node "<<n.getID()<<std::endl;
            continue;
        }
        math::Vector vx(evx.X(),evx.Y(),evx.Z());
        math::Vector vy(evy.X(),evy.Y(),evy.Z());
        math::Vector vz(evz.X(),evz.Y(),evz.Z());
        
        math::Point p1 = center + (vx + vy - vz)*cube_size;
        
        Node n1 = mesh_cube.newNode(p1);
        math::Point p2 = center + (vx - vy - vz)*cube_size;
        Node n2 = mesh_cube.newNode(p2);
        math::Point p3 = center + (vx + vy + vz).opp()*cube_size;
        Node n3 = mesh_cube.newNode(p3);
        math::Point p4 = center + (vy - vx - vz)*cube_size;
        Node n4 = mesh_cube.newNode(p4);
        
        math::Point p5 = center + (vx + vy + vz)*cube_size;
        Node n5 = mesh_cube.newNode(p5);
        math::Point p6 = center + (vx - vy + vz)*cube_size;
        Node n6 = mesh_cube.newNode(p6);
        math::Point p7 = center + (vx + vy - vz).opp()*cube_size;
        Node n7 = mesh_cube.newNode(p7);
        math::Point p8 = center + (vy - vx + vz)*cube_size;
        Node n8 = mesh_cube.newNode(p8);
        //        std::cout<<"Cube from "<<std::endl;
        //        std::cout<<"\t "<<p1<<std::endl;
        //        std::cout<<"\t "<<p2<<std::endl;
        //        std::cout<<"\t "<<p3<<std::endl;
        //        std::cout<<"\t "<<p4<<std::endl;
        //        std::cout<<"\t "<<p5<<std::endl;
        //        std::cout<<"\t "<<p6<<std::endl;
        //        std::cout<<"\t "<<p7<<std::endl;
        //        std::cout<<"\t "<<p8<<std::endl;
        Region r = mesh_cube.newHex(n1, n2, n3, n4, n5, n6, n7, n8);
        if(m_mesh->isMarked(n, m_mark_node_on_point))
            (*v)[r.getID()]=0;
        else if(m_mesh->isMarked(n, m_mark_node_on_curve))
            (*v)[r.getID()]=1;
        else if(m_mesh->isMarked(n, m_mark_node_on_surface))
            (*v)[r.getID()]=2;
        else
            (*v)[r.getID()]=3;
        
    }
    VTKWriter<IGMesh> writer_cube(mesh_cube);
    
    std::stringstream file_name_cube;
    file_name_cube<<"HLBFGS_DEBUG_CUBE_" << nb_file;
    writer_cube.write(file_name_cube.str(), DIM3 | R | N);
    
    
    
    Variable<int>* var_sing = 0;
    try{
        var_sing = m_mesh->getVariable<int>(GMDS_REGION, "sing_tet");
    }
    catch (GMDSException& e){
        var_sing = m_mesh->newVariable<int>(GMDS_REGION, "sing_tet");
    }
    
    IGMesh::region_iterator itr = m_mesh->regions_begin();
    
    int nbColoredTet = 0;
    //=========================================================================
    // INTERN SKELETON CREATION
    //=========================================================================
    // FIRST LOOP ON REGIONS TO GET ALL THE 3-SING. TETS
    //=========================================================================
    for (; !itr.isDone(); itr.next()){
        Region current_region = itr.value();
        std::vector<Node> nodes = current_region.get<Node>();
        bool onPnt = false;
        for (unsigned int i_node = 0; i_node < nodes.size(); i_node++){
            Node ni = nodes[i_node];
            if (m_mesh->isMarked(ni, m_mark_node_on_point))
                onPnt = true;
        }
        if (onPnt){
            (*var_sing)[current_region.getID()] = 0;
        }
        else
        {
            std::vector<TCellID> node_ids = current_region.getIDs<Node>();
            math::Quaternion q[4];
            for(int i_n=0; i_n<4; i_n++){
                math::AxisAngleRotation rot =(*m_rotation_field)[node_ids[i_n]];
                math::Chart ci = rot.toChart();

                q[i_n]= math::Quaternion(ci);
            }
            
            int sing_type = math::Chart::testSingularity(q[0],q[1],q[2],q[3]);
            if (sing_type != 0)
                nbColoredTet++;
            
            (*var_sing)[current_region.getID()] = sing_type;
        }
    }
    std::cout << "Nb. colored tetrahedra: " << nbColoredTet << std::endl;
    
    
    
    VTKWriter<IGMesh> writer(*m_mesh);
    
    std::stringstream file_name;
    file_name<<"HLBFGS_DEBUG_" << nb_file;
    writer.write(file_name.str(), DIM3 | R| N);
    
    
    nb_file++;
}

/*---------------------------------------------------------------------------*/


