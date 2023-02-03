//
// Created by rochec on 24/10/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/DiffusionEquation2D.h>
#include <Eigen/Sparse>
#include <Eigen/Eigen>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

DiffusionEquation2D::DiffusionEquation2D(Mesh *AMesh, int AmarkFrontNodes_int, int AmarkFrontNodes_out, Variable<double>* Adistance) :
  m_mesh(AMesh),
  m_mass(m_mesh->getNbNodes(), m_mesh->getNbNodes()),
  m_stiffness(m_mesh->getNbNodes(), m_mesh->getNbNodes()),
  m_sol_0(m_mesh->getNbNodes()),
  m_sol_n(m_mesh->getNbNodes())
{
	m_markNodes_int = AmarkFrontNodes_int;
	m_markNodes_out = AmarkFrontNodes_out;
	m_distance = Adistance;

	// m_id_local_index gives us a compact numerotation of the nodes of the mesh
	int index=0;
	for(auto n_id:m_mesh->nodes()) {
		m_id_local_index[n_id]=index;
		index++;
	}

	m_dt = 10.0;
	m_sigma = 1.0;
	m_it_max = 1000;

}


/*------------------------------------------------------------------------*/
DiffusionEquation2D::STATUS
DiffusionEquation2D::execute()
{
	initialisation();
	// Boucle en temps
	std::cout << "Time loop..." << std::endl;
	for (int it = 1; it <= m_it_max; it++)
	{
		oneTimeStep();
	}

	for (auto n_id:m_mesh->nodes())
	{
		m_distance->set(n_id, m_sol_n.coeffRef(m_id_local_index[n_id]));
	}

	return DiffusionEquation2D::SUCCESS;
}
/*------------------------------------------------------------------------*/






//
//		GEOMETRY FUNCTIONS
//

/*------------------------------------------------------------------------*/
// Construction des 3 fonctions de base
// (phi_hat_1, phi_hat_2, phi_hat_3)
double DiffusionEquation2D::phi_hat(int hat_i, math::Point X_hat)
{
	double value;
	switch (hat_i)
	{
	case 1:
		value=1-X_hat.X()-X_hat.Y();
		break;
	case 2:
		value=X_hat.X();
		break;
	case 3:
		value=X_hat.Y();
		break;
	default:
		std::cout << "Please choose a valid number of hat_i for the reference triangle (1, 2 or 3)" << std::endl;
		abort();
	}
	return value;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction des 3 gradients des fonctions de base
// (gradphihat1, gradphihat2, gradphihat3)
// (indépendant du vecteur X_hat)
Eigen::Vector2d DiffusionEquation2D::grad_phi_hat(int hat_i)
{
	Eigen::Vector2d gradphi;
	switch (hat_i)
	{
	case 1:
		gradphi(0) = -1; gradphi(1) = -1;
		break;
	case 2:
		gradphi(0) = 1; gradphi(1) = 0;
		break;
	case 3:
		gradphi(0) = 0; gradphi(1) = 1;
		break;
	default:
		std::cout << "Please choose a valid number of hati for the reference triangle (1, 2 or 3)" << std::endl;
		abort();
	}
	return gradphi;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction de la fonction de transformation
// du triangle de référence K_hat en K
Eigen::Vector2d DiffusionEquation2D::FK(TCellID triK_id, math::Point X_hat)
{
	Face triK = m_mesh->get<Face>(triK_id);
	std::vector<Node> tri_nodes = triK.get<Node>() ;

	Eigen::Vector2d X;
	X(0) = tri_nodes[0].X()* phi_hat(1, X_hat) + tri_nodes[1].X()* phi_hat(2, X_hat) + tri_nodes[2].X()* phi_hat(3, X_hat);
	X(1) = tri_nodes[0].Y()* phi_hat(1, X_hat) + tri_nodes[1].Y()* phi_hat(2, X_hat) + tri_nodes[2].Y()* phi_hat(3, X_hat);

	return X;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction de la jacobienne de la fonction de transformation
// du triangle de référence K_hat en K
// (indépendant du vecteur X_hat)
Eigen::Matrix2d DiffusionEquation2D::JFK(TCellID triK_id)
{
	Face triK = m_mesh->get<Face>(triK_id);
	std::vector<Node> tri_nodes = triK.get<Node>() ;

	Eigen::Matrix2d X;
	X(0,0) = tri_nodes[1].X() - tri_nodes[0].X();
	X(0,1) = tri_nodes[2].X() - tri_nodes[0].X();
	X(1,0) = tri_nodes[1].Y() - tri_nodes[0].Y();
	X(1,1) = tri_nodes[2].Y() - tri_nodes[0].Y();

	return X;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction de la valeur absolue du déterminant de la jacobienne
// de la fonction de transformation
// du triangle de référence hatK en K
// (indépendant du vecteur hatX)
double DiffusionEquation2D::absdetJFK(TCellID triK_id)
{
	return fabs((JFK(triK_id)).determinant());
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction du gradient
// (gradphi_i) avec i qui est le hati-ème sommet de K
// appliqué en F(hatX)
Eigen::Vector2d DiffusionEquation2D::grad_phi(int i_hat, TCellID triK_id)
{
	return (((JFK(triK_id)).inverse()).transpose())*grad_phi_hat(i_hat);
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Points et poids de quadrature de la formule des sommets (intégration 2D)
void DiffusionEquation2D::quadraturePointsAndWeightsSummitpointFormula(std::vector<double> &weights, std::vector<math::Point> &points)
{
	weights.push_back(1./3.);
	weights.push_back(1./3.);
	weights.push_back(1./3.);

	points.push_back(math::Point(0.,0.));
	points.push_back(math::Point(1.,0.));
	points.push_back(math::Point(0.,1.));
}
/*------------------------------------------------------------------------*/



//
// Assemble Matrix
//

/*------------------------------------------------------------------------*/
void DiffusionEquation2D::assembleMassAndStiffnessMatrices()
{
	std::vector<Eigen::Triplet<double>> tripletsMass, tripletsStiffness;
	m_mass.setZero();
	m_stiffness.setZero();

	std::vector<double> weights;
	std::vector<math::Point> quadrature_points;
	quadraturePointsAndWeightsSummitpointFormula(weights, quadrature_points);

	// Loop on the element K in T_h
	for (auto triK_id:m_mesh->faces())
	{
		Face K = m_mesh->get<Face>(triK_id);
		std::vector<Node> K_nodes = K.get<Node>();
		// Loop on the quadrature points
		for (int p_id=0; p_id < quadrature_points.size(); p_id++)
		{
			math::Point p = quadrature_points[p_id];
			double w = weights[p_id];
			// Loop on the nodes of the triangle K
			for (int i_hat=1; i_hat<4; i_hat++)
			{
				// Loop on the nodes of the triangle K
				for (int j_hat=1; j_hat<4; j_hat++)
				{
					// i = global id of the node i_hat of triangle K
					int i( m_id_local_index[K_nodes[i_hat-1].id()] );
					// j = global id of the node j_hat of triangle K
					int j( m_id_local_index[K_nodes[j_hat-1].id()] );

					// M(i,j) = M(i,j) + contribution de la _masse (hati,hatj) en hatX_q
					m_mass.coeffRef(i,j) += 0.5				// Area of the ref triangle K_hat
					                         			*w					// Weight of the quadrature points
					   										* phi_hat(i_hat, p)
					   										* phi_hat(j_hat, p)
					   										* absdetJFK(K.id());

					// K(i,j) = K(i,j) + contribution de la rigidité (hati,hatj) en hatX_q
					m_stiffness.coeffRef(i,j) += 0.5		// Area of the ref triangle K_hat
					   								*w						// Weight of the quadrature points
					   								* grad_phi(i_hat, K.id()).dot(grad_phi(j_hat, K.id()))
					   								* absdetJFK(K.id());

				}
			}
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void DiffusionEquation2D::applyBCToSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> & systemMatrix)
{
	Eigen::SparseVector<double> zeroRow(systemMatrix.cols());
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (m_mesh->isMarked(n, m_markNodes_int) || m_mesh->isMarked(n, m_markNodes_out) )
		{
			int index = m_id_local_index[n_id];
			systemMatrix.row(index) = zeroRow;
			systemMatrix.coeffRef(index, index) = 1.0;
		}
	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void DiffusionEquation2D::applyBCToSecondMember(Eigen::SparseVector<double> & secondMember)
{
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (m_mesh->isMarked(n, m_markNodes_int) )
		{
			int index = m_id_local_index[n_id];
			secondMember.coeffRef(index) = 0.0;
		}
		if (m_mesh->isMarked(n, m_markNodes_out))
		{
			int index = m_id_local_index[n_id];
			secondMember.coeffRef(index) = 1.0;
		}
	}
}
/*------------------------------------------------------------------------*/






//
// Time Scheme
//

/*------------------------------------------------------------------------*/
void DiffusionEquation2D::initialisation()
{
	// Construction de la matrice du système

	std::cout << "Assembling" << std::endl;
	assembleMassAndStiffnessMatrices();
	Eigen::SparseMatrix<double,Eigen::RowMajor> systemMatrix(m_mass+m_dt*m_sigma*m_stiffness);

	// Application des conditions aux bords sur la matrice
	applyBCToSystemMatrix(systemMatrix);

	// On envoie la matrice du système au solver
	std::cout << "Preprocessing of the matrix" << std::endl;
	//Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > solver;
	// Compute the ordering permutation vector from the structural pattern of A
	m_solver.analyzePattern(systemMatrix);
	// Compute the numerical factorization
	m_solver.factorize(systemMatrix);

	// On recupère la solution initiale
	m_sol_0.setZero();
	for (auto n_id:m_mesh->nodes())
	{
		int index = m_id_local_index[n_id];
		m_sol_0.coeffRef(index) = 0.0;
		m_sol_n.coeffRef(index) = 0.0;
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void DiffusionEquation2D::oneTimeStep()
{
	// On avance le temps
	m_t += m_dt;
	m_it += 1;
	// Construction du second membre
	Eigen::SparseVector<double> secondMember(m_mass*m_sol_n);

	// Application des conditions sur le terme de droite
	applyBCToSecondMember(secondMember);
	// Résolution du système
	//std::cout << "Solve system at time " << m_t << std::endl;
	m_sol_n = m_solver.solve(secondMember);
}
/*------------------------------------------------------------------------*/
