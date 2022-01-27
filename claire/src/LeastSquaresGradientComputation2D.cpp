//
// Created by rochec on 27/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/LeastSquaresGradientComputation2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


LeastSquaresGradientComputation2D::LeastSquaresGradientComputation2D(Mesh *AMesh, Variable<double>* Adistance) {
	m_mesh = AMesh;
	m_distance = Adistance;
	m_gradient2D = m_mesh->newVariable<math::Vector3d ,GMDS_FACE>("gradient_2D");
}


/*------------------------------------------------------------------------*/
LeastSquaresGradientComputation2D::STATUS LeastSquaresGradientComputation2D::execute()
{

	for(auto face_id:m_mesh->faces()) {
		//math::Vector3d Gradient = computeGradientOnSimpleFace(face_id);
		//m_gradient2D->set(face_id, Gradient) ;
	}

	return LeastSquaresGradientComputation2D::SUCCESS;
}
/*------------------------------------------------------------------------*/



/*------------------------------------------------------------------------*/
void LeastSquaresGradientComputation2D::buildMatrix(TCellID n_id, Eigen::SparseMatrix<double> M, Eigen::VectorXd b){

	Eigen::VectorXd Di;
	Eigen::SparseMatrix<double> Wi;
	Eigen::SparseMatrix<double> Ai;
	Eigen::SparseMatrix<double> AiT;

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

	for (int j=0; j < adjacent_nodes.size(); j++){
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