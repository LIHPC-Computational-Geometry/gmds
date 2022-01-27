//
// Created by rochec on 27/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/LeastSquaresGradientComputation.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

LeastSquaresGradientComputation::LeastSquaresGradientComputation(Mesh *AMesh, Variable<double>* Adistance) {
	m_mesh = AMesh;
	m_distance = Adistance;
	m_gradient2D = m_mesh->newVariable<math::Vector3d, GMDS_NODE>("gradient_2D");
}


/*------------------------------------------------------------------------*/
LeastSquaresGradientComputation::STATUS
LeastSquaresGradientComputation::execute()
{

	for(auto node_id:m_mesh->nodes()) {
		std::cout << "------------------------" << std::endl;
		std::cout << "NOEUD TRAITé : " << node_id << std::endl;
		math::Vector3d Gradient = computeGradientOnSimpleVertex(node_id);
		m_gradient2D->set(node_id, Gradient) ;
	}

	return LeastSquaresGradientComputation::SUCCESS;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
math::Vector3d
LeastSquaresGradientComputation::computeGradientOnSimpleVertex(TCellID node_id){

	math::Vector3d Gradient;

	int dim=3;
	// Construction de la matrice et des vecteurs
	Eigen::SparseMatrix<double> A(dim,dim);
	Eigen::VectorXd x(dim);
	Eigen::VectorXd b(dim);
	buildMatrix(node_id, A, b);

	// Résolution du système avec gradient conjugué
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
	cg.compute(A);
	x = cg.solve(b);
	std::cout << "#iterations:     " << cg.iterations() << std::endl;
	std::cout << "estimated error: " << cg.error()      << std::endl;

	Gradient.setXYZ( x[0], x[1], x[2] );

	return Gradient;
}
/*------------------------------------------------------------------------*/



/*------------------------------------------------------------------------*/
void
LeastSquaresGradientComputation::buildMatrix(TCellID n_id, Eigen::SparseMatrix<double> &M, Eigen::VectorXd &b){

	// Dimensions des matrices
	int dim=3;
	int ki;

	Node n = m_mesh->get<Node>(n_id);
	math::Point p_i = n.point() ;
	std::vector<Edge> adjacent_edges = n.get<Edge>() ;
	std::vector<TCellID> adjacent_nodes ;
	std::vector<math::Point> adjacent_points ;
	for(auto edge:adjacent_edges) {
		TCellID n_ik = edge.getOppositeNodeId(n);
		adjacent_nodes.push_back(n_ik);
		Node node = m_mesh->get<Node>(n_ik);
		adjacent_points.push_back(node.point());
	}
	// Normalement, de cette manière, l'arête stockée à l'indice j dans adjacent_edges
	// est celle qui correspond au noeud opposé dans le vecteur adjacent_nodes.
	ki = adjacent_nodes.size();

	// Initialisation des matrices
	Eigen::VectorXd Di (ki);
	Eigen::SparseMatrix<double> Wi(ki,ki);
	Eigen::SparseMatrix<double> Ai(ki, dim);
	Eigen::SparseMatrix<double> AiT(dim, ki);

	for (int j=0; j < ki; j++){
		// Initialisation de la matrice A
		Ai.coeffRef(j, 0) = adjacent_points[j].X() - p_i.X();
		Ai.coeffRef(j, 1) = adjacent_points[j].Y() - p_i.Y();
		Ai.coeffRef(j, 2) = adjacent_points[j].Z() - p_i.Z();
		// Initialisation de la matrice des poids Wi
		Wi.coeffRef(j, j) = 1.0/pow(adjacent_edges[j].length(),2);
		// Initialisation du vecteur colonne Di
		double Dij = (m_distance->value( adjacent_nodes[j]) - m_distance->value(n_id) );
		Di[j]= Dij;
	}

	AiT= Ai.transpose();
	M = AiT*Wi*Ai;
	M.makeCompressed();

	b = AiT*Wi*Di;

}
/*------------------------------------------------------------------------*/