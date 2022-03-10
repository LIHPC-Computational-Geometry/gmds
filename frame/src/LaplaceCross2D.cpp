/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Point.h>
#include <gmds/math/Numerics.h>
/*---------------------------------------------------------------------------*/
#include <gmds/frame/LaplaceCross2D.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
LaplaceCross2D::LaplaceCross2D(gmds::Mesh* AMesh,
			       gmds::Variable<gmds::math::Cross2D>* AField,
			       std::vector<gmds::Node>& AFromNodes,
			       std::vector<gmds::Node>& AToNodes,
			       std::vector<gmds::Face>& AFaces)
  : m_mesh(AMesh), m_field(AField), 
    m_from_nodes(AFromNodes),
    m_to_nodes(AToNodes)
{
  m_nodes.resize(AFromNodes.size()+AToNodes.size());
  int index=0;
  
  for(unsigned int i=0;i<AFromNodes.size();i++) {
    m_nodes[index] = AFromNodes[i];
    m_id_to_local_index[AFromNodes[i].id()]=index;
    index++;
  }
  int index_to_I=0;
  for(unsigned int i=0;i<AToNodes.size();i++) {
    m_nodes[index] = AToNodes[i];
    m_id_to_I[AToNodes[i].id()]=index_to_I;
    index++;
    index_to_I++;
  }

  /* We have the mesh nodes stored in m_nodes. It gives us a compact 
   * numerotation that we will use. 
   */

  m_faces = AFaces;
}
/*---------------------------------------------------------------------------*/
std::vector<Edge> LaplaceCross2D::
getEdgesOnCurve(const Node& ANode) const
{
  std::vector<Edge> edges_on_curve;
  std::vector<Edge> adj_edges = ANode.get<Edge>();
  for (unsigned int i = 0; i < adj_edges.size(); i++)
    {
      Edge ei = adj_edges[i];
      std::vector<Node> nodes = ei.get<Node>();
      std::vector<Face> faces = ei.get<Face>();
      if (m_mesh->isMarked(nodes[0], m_mark_dirichlet) && 
	  m_mesh->isMarked(nodes[1], m_mark_dirichlet) &&
	  faces.size()==1)
	edges_on_curve.push_back(ei);
    }
  return edges_on_curve;
}
/*---------------------------------------------------------------------------*/
Node LaplaceCross2D::
getNeighboorOn(const Node& ANode, const Edge& AEdge) const
{
  std::vector<Node> nodes = AEdge.get<Node>();
  if (nodes[0].id() == ANode.id())
    return nodes[1];

  return nodes[0];
}
/*---------------------------------------------------------------------------*/
void LaplaceCross2D::execute()
{
  m_mark_dirichlet = m_mesh->newMark<Node>();
  
  for(unsigned int i=0; i<m_from_nodes.size(); i++)
    m_mesh->mark(m_from_nodes[i],m_mark_dirichlet);
 

  // Unknown value are on the internal nodes, which are stored in 
  // m_to_nodes.
  TInt nb_unknowns = m_to_nodes.size();

  //==================================================================
  // STEP 1 / Matrix assembly  
  //==================================================================
  //std::cout<<"System assembly"<<std::endl;
  Eigen::SparseMatrix<double> Acos(nb_unknowns,nb_unknowns); 
  Eigen::VectorXd             bcos(nb_unknowns); 
  Eigen::SparseMatrix<double> Asin(nb_unknowns,nb_unknowns); 
  Eigen::VectorXd             bsin(nb_unknowns); 

  //We begin by initialize the Sparse Matrix DS
  int max_adj = getMaxAdjacency();

  //we reserve a room for max_adj non_zeros per columns
  Acos.reserve(Eigen::VectorXi::Constant(nb_unknowns,max_adj));
  Asin.reserve(Eigen::VectorXi::Constant(nb_unknowns,max_adj));

  buildStiffnessMatrix(Acos);
  buildStiffnessMatrix(Asin);

  //We finalize the Sparse Matrix DS
  Acos.makeCompressed();
  Asin.makeCompressed();

  buildRightHandDirichlet(bcos,bsin);
  //==================================================================
  // STEP 2 / System Solving
  //==================================================================
  //std::cout<<"System solving for cosinus"<<std::endl;
  Eigen::VectorXd   xcos(nb_unknowns);
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > cholcos(Acos);
  xcos = cholcos.solve(bcos);

  /*Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg_cos;
  cg_cos.compute(Acos);
  xcos = cg_cos.solve(bcos);
  */

  //std::cout<<"System solving for sinus"<<std::endl;
  Eigen::VectorXd   xsin(nb_unknowns);
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > cholsin(Asin);
  xsin = cholsin.solve(bsin);

  /*Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg_sin;
  cg_sin.compute(Asin);
  xsin = cg_sin.solve(bsin);*/
  //==================================================================
  // STEP 3 / Assign the value onto the mesh nodes
  //==================================================================
  //std::cout<<"Value assignment to mesh node"<<std::endl; 
  for(unsigned int i=0; i<m_to_nodes.size();i++) {
    Node    n_i = m_to_nodes[i];
    TCellID n_gmds_id = n_i.id();
    int     n_local_id =  m_id_to_I[n_gmds_id];
    
    double  cosA = xcos[n_local_id];
    double  sinA = xsin[n_local_id];

    double  a = atan2(sinA,cosA);

     (*m_field)[n_gmds_id] = math::Cross2D(a);
  }

  m_mesh->unmarkAll<Node>(m_mark_dirichlet);
  m_mesh->freeMark<Node> (m_mark_dirichlet);

}
/*---------------------------------------------------------------------------*/
void LaplaceCross2D::buildStiffnessMatrix(Eigen::SparseMatrix<double>& AM)
{ 
  //==================================================================
  // We compute the global stiffness
  //==================================================================
  //First loop on the faces (triangles)
  for(unsigned int k=0; k < m_faces.size(); k++){
      std::vector<Node> current_nodes;
      m_faces[k].get<Node>(current_nodes);
      math::Point p1 = current_nodes[0].point();
      math::Point p2 = current_nodes[1].point();
      math::Point p3 = current_nodes[2].point();
      math::Matrix<3,3,double> s = math::stiffnessMatrix2D(p1,p2,p3);
      
      // Second loop on the internal  nodes
      // with i,j their index locally to the face, and I,J the index global
      // to the mesh
      for(unsigned int i=0; i < current_nodes.size(); i++){
          Node n_i = current_nodes[i];
          if(m_mesh->isMarked(n_i,m_mark_dirichlet))
              continue;
          TCellID I =  m_id_to_I[n_i.id()];
          
          for(unsigned int j=0; j < current_nodes.size(); j++){
              Node n_j = current_nodes[j];
              if(m_mesh->isMarked(n_j,m_mark_dirichlet))
                  continue;
              TCellID J =  m_id_to_I[n_j.id()];
              
              AM.coeffRef(I,J) +=s.get(i,j); 
              
          }
          
      }
  }

}
/*---------------------------------------------------------------------------*/
void LaplaceCross2D::buildRightHandDirichlet(Eigen::VectorXd& ABCos, 
					     Eigen::VectorXd& ABSin )
{
  for(unsigned int i=0; i <  m_to_nodes.size(); i++){
    //n_i is an inner node of the mesh
    Node n_i = m_to_nodes[i]; 
    int n_id = m_id_to_I[n_i.id()];
    ABCos[n_id] = 0.0;
    ABSin[n_id] = 0.0;
    std::vector<Face> faces_of_n_i = n_i.get<Face>();
   
    for(unsigned int j=0; j<faces_of_n_i.size(); j++){

      Face f_j = faces_of_n_i[j];
      std::vector<Node> nodes_of_f_j = f_j.get<Node>();
      bool is_boundary_face=false;

      math::Vector3d Ud_cos(0,0,0);
      math::Vector3d Ud_sin(0,0,0);
      int n_i_index_in_f_j=-1;
      //we store if the i^th node of f_j is on the boundary or not

      for(unsigned int k=0; k < nodes_of_f_j.size(); k++){
	
	Node n_k = nodes_of_f_j[k];
	if(m_mesh->isMarked(n_k, m_mark_dirichlet)) {
	  is_boundary_face = true;
	  math::Cross2D cross = (*m_field)[n_k.id()];
	  double ref_angle = cross.referenceAngle();
	  double a =math::modulo2PI(cross.referenceAngle());
	  Ud_cos[k] = cos(a);
	  Ud_sin[k] = sin(a);
	  }
	if(n_k.id()==n_i.id()){
	  n_i_index_in_f_j=k;
	}
      }
      if(is_boundary_face){
	math::Point pi = nodes_of_f_j[0].point();
	math::Point pj = nodes_of_f_j[1].point();
	math::Point pk = nodes_of_f_j[2].point();
	math::Matrix<3,3,double> s = math::stiffnessMatrix2D(pi,pj,pk);

	//Ud is build
	math::Vector3d v(s.get(0,n_i_index_in_f_j),
		       s.get(1,n_i_index_in_f_j),
		       s.get(2,n_i_index_in_f_j));
	ABCos[n_id] = -v.dot(Ud_cos);
	ABSin[n_id] = -v.dot(Ud_sin);
      }
    
    }

  }
}
/*---------------------------------------------------------------------------*/
int LaplaceCross2D::getMaxAdjacency() const
{
  double m = 0;
  for(unsigned int i=0; i < m_nodes.size(); i++){
    Node n_i = m_nodes[i];
    std::vector<Edge> edges_i = n_i.get<Edge>();
    if(edges_i.size()>m)
      m = edges_i.size();
  }
  return m;
}

/*---------------------------------------------------------------------------*/
