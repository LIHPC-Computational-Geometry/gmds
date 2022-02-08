//
// Created by rochec on 04/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AdvectedPointRK4_2D.h>
#include <gmds/claire/ResLU.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AdvectedPointRK4_2D::AdvectedPointRK4_2D(Mesh *AMesh, math::Point A_Pstart, double A_d0, Variable<double>* A_distance, Variable<math::Vector3d>* A_gradient2D) {
	m_mesh = AMesh;
	m_Pstart = A_Pstart;
	m_Pend = A_Pstart;
	m_d0 = A_d0;
	m_distance = A_distance;
	m_gradient2D = A_gradient2D;
	m_discrete_path.push_back(m_Pstart);
}


/*------------------------------------------------------------------------*/
AdvectedPointRK4_2D::STATUS AdvectedPointRK4_2D::execute()
{
	double minLenght = minEdgeLenght();
	std::cout << "min : " << minLenght << std::endl;

	double dt = 0.9*minLenght;
	int iterations=0;
	int max_iterations=10000;

	double dist(0);
	math::Vector3d Grad;
	TCellID face_id_P = inWhichTriangle(m_Pstart) ;
	computeInterpolatedDistanceAndGrad(m_Pstart, face_id_P, dist, Grad);
	m_Pend = m_Pstart;
	std::cout << "distance initiale : " << dist << std::endl;
	std::cout << "gradient initial : " << Grad << std::endl;

	while ( (abs(dist-m_d0) > pow(10,-6)) && iterations < max_iterations ) {
		std::cout << "-------- ITERATION " << iterations << " ---------" << std::endl;
		math::Point M = RungeKutta4(m_Pend, Grad.normalize(), dt);
		double dist_M;
		math::Vector3d Grad_M;
		TCellID face_id = inWhichTriangle(M) ;
		if(face_id != NullID) {
			computeInterpolatedDistanceAndGrad(M, face_id, dist_M, Grad_M);
			// Rajouter condition comme ça : inWhichTriangle(M) != std::numeric_limits<int>::max()
			if ((dist_M < m_d0)) {
				m_Pend = M;
				Grad = Grad_M;
				dist = dist_M;
				std::cout << "distance : " << dist << std::endl;
				m_discrete_path.push_back(m_Pend);
			}
			else {
				// Alors on a dépassé la distance souhaitée, m_Pend n'est pas mis à jour, on retranche le pas de temps
				// et on continue la boucle
				dt = 0.95 * dt;
			}
		}
		else{
			dt = 0.95 * dt;
		}
		iterations += 1;
	}

	face_id_P = inWhichTriangle(m_Pend) ;
	computeInterpolatedDistanceAndGrad(m_Pend, face_id_P, dist, Grad);
	std::cout << "----------------------------------" << std::endl;
	std::cout << "Point final : " << m_Pend << std::endl;
	std::cout << "Distance finale : " << dist << std::endl;

	writeDiscretePathInVTK();

	return AdvectedPointRK4_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/





/*------------------------------------------------------------------------*/
bool AdvectedPointRK4_2D::isInTriangle(TCellID face_id, math::Point M){
	bool isInFace(false);
	Face face = m_mesh->get<Face>(face_id);
	std::vector<TCellID> face_nodes_ids = face.getIDs<Node>();
	TCellID vi_id = face_nodes_ids[0];
	TCellID vj_id = face_nodes_ids[2];
	TCellID vk_id = face_nodes_ids[1];

	Node vi = m_mesh->get<Node>(vi_id);
	Node vj = m_mesh->get<Node>(vj_id);
	Node vk = m_mesh->get<Node>(vk_id);

	math::Point vi_coord = vi.point();
	math::Point vj_coord = vj.point();
	math::Point vk_coord = vk.point();

	math::Vector3d vij = vj_coord-vi_coord ;
	math::Vector3d vjk = vk_coord-vj_coord ;
	math::Vector3d vki = vi_coord-vk_coord ;
	math::Vector3d viM = M-vi_coord ;
	math::Vector3d vjM = M-vj_coord ;
	math::Vector3d vkM = M-vk_coord ;

	double d1 = ( vij.cross(viM) ).dot( viM.cross(-vki) ) ;
	double d2 = ( -vij.cross(vjM) ).dot( vjM.cross(vjk) ) ;
	double d3 = ( vki.cross(vkM) ).dot( vkM.cross(-vjk) ) ;

	if (d1 >= 0 && d2 >= 0 && d3 >= 0) {
		isInFace = true;
	}
	return isInFace;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
TCellID AdvectedPointRK4_2D::inWhichTriangle(math::Point M){
	TCellID face_id;
	bool isInFace(false);
	/*
	bool found_item=false;
	for(auto f_it = m_mesh->faces_begin(); f_it!= m_mesh->faces_end() && !found_item;++f_it){
		TCellID f_id = *f_it;
		found_item=true;
	}
	 */

	for (auto f_id:m_mesh->faces()){
		if(!isInFace){
			isInFace = isInTriangle(f_id, M);
			if(isInFace){
				face_id = f_id;
			}
		}
	}
	if (!isInFace){
		face_id = NullID;
		std::cout << "triangle : " << face_id << std::endl;
	}
	return face_id;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AdvectedPointRK4_2D::minEdgeLenght(){
	// Initialisation avec une arête prise au hasard. Pb : si l'arête 0 a été
	// retirée
	Edge edge_0 = m_mesh->get<Edge>(0);
	double minLenght(edge_0.length());
	for (auto edge_id:m_mesh->edges()){
		Edge edge = m_mesh->get<Edge>(edge_id);
		if(edge.length() < minLenght){
			minLenght = edge.length() ;
			//std::cout << "Edge id : " << edge_id << ", Taille : " << minLenght << std::endl;
		}
	}
	return minLenght;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AdvectedPointRK4_2D::buildSystemMatrix(TCellID face_id, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b1,
                                       Eigen::VectorXd &b2, Eigen::VectorXd &b3){

	Face f = m_mesh->get<Face>(face_id);
	std::vector<Node> nodes = f.get<Node>() ;
	math::Point p0 = nodes[0].point() ;
	math::Point p1 = nodes[1].point() ;
	math::Point p2 = nodes[2].point() ;

	TCellID n0_id = nodes[0].id();
	TCellID n1_id = nodes[1].id();
	TCellID n2_id = nodes[2].id();

	// Assemble la matrice composée des coordonnées des points des noeuds de la face
	// x1 y1 1
	// x2 y2 1
	// x3 y3 1
	A.coeffRef(0,0) = p0.X() ;
	A.coeffRef(0,1) = p0.Y() ;
	A.coeffRef(0,2) = 1.0 ;
	A.coeffRef(1,0) = p1.X() ;
	A.coeffRef(1,1) = p1.Y() ;
	A.coeffRef(1,2) = 1.0 ;
	A.coeffRef(2,0) = p2.X() ;
	A.coeffRef(2,1) = p2.Y() ;
	A.coeffRef(2,2) = 1.0 ;

	// Initialisation du vecteur des distances aux noeuds
	b1[0] = m_distance->value(n0_id)  ;
	b1[1] = m_distance->value(n1_id)  ;
	b1[2] = m_distance->value(n2_id)  ;

	//Initialisation du vecteur des composantes du grad dans la direction x
	b2[0] = m_gradient2D->value(n0_id).X()  ;
	b2[1] = m_gradient2D->value(n1_id).X()  ;
	b2[2] = m_gradient2D->value(n2_id).X()  ;

	//Initialisation du vecteur des composantes du grad dans la direction y
	b3[0] = m_gradient2D->value(n0_id).Y()  ;
	b3[1] = m_gradient2D->value(n1_id).Y()  ;
	b3[2] = m_gradient2D->value(n2_id).Y()  ;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AdvectedPointRK4_2D::computeInterpolatedDistanceAndGrad(math::Point M, TCellID face_id, double &distance_M, math::Vector3d &Grad_M){

	// On commence par chercher dans quel élément (face) est le point M.
	//TCellID face_id = inWhichTriangle(M) ;
	std::cout << "Noeud dans le triangle : " << face_id << std::endl;

	// En 2D, résolution de 3 systèmes 3x3. Un pour la distance, et deux pour le gradient (un pour chaque composante)
	// Assemblage de la matrice A commune aux 3 systèmes
	Eigen::SparseMatrix<double> A(3,3);
	Eigen::VectorXd dist(3);	// Vecteur des distances aux trois noeuds
	Eigen::VectorXd grad_x(3);	// Vecteur des gradients aux trois noeuds dans la direction x
	Eigen::VectorXd grad_y(3);	// Vecteur des gradients aux trois noeuds dans la direction y
	buildSystemMatrix(face_id, A, dist, grad_x, grad_y);

	// Définition de la distance exacte au point M
	distance_M = interpolationDistance(A, dist, M);

	// Calcul du gradient exact par interpolation sur la face
	Grad_M = interpolationGradient(A, grad_x, grad_y, M);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AdvectedPointRK4_2D::interpolationDistance(Eigen::SparseMatrix<double> A, Eigen::VectorXd b, math::Point M){
	double dist_interpolee;

	ResLU resolutionlu(A, b);
	resolutionlu.execute();
	Eigen::VectorXd coeff = resolutionlu.getSolution();

	dist_interpolee = coeff[0]*M.X() + coeff[1]*M.Y() + coeff[2] ;

	return dist_interpolee;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d AdvectedPointRK4_2D::interpolationGradient(Eigen::SparseMatrix<double> A, Eigen::VectorXd bx, Eigen::VectorXd by, math::Point M){
	math::Vector3d Gradient_M;

	// Définition du gradient exact dans la direction x au point M
	ResLU resolutionlu(A, bx);
	resolutionlu.execute();
	Eigen::VectorXd coeff_grad_x = resolutionlu.getSolution();
	Gradient_M.setX(coeff_grad_x[0]*M.X() + coeff_grad_x[1]*M.Y() + coeff_grad_x[2]) ;


	// Définition du gradient exact dans la direction y au point M
	resolutionlu = ResLU(A, by);
	resolutionlu.execute();
	Eigen::VectorXd coeff_grad_y = resolutionlu.getSolution();
	Gradient_M.setY(coeff_grad_y[0]*M.X() + coeff_grad_y[1]*M.Y() + coeff_grad_y[2]) ;

	// Définition du gradient exact dans la direction z au point M
	Gradient_M.setZ(0.0);

	return Gradient_M;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Point AdvectedPointRK4_2D::RungeKutta4(math::Point yn, math::Vector3d grad_yn, double dt){
	math::Point ynew;

	math::Vector3d k1 = grad_yn ;
	math::Vector3d k2 = grad_yn+(dt/2.0)*k1 ;
	math::Vector3d k3 = grad_yn+(dt/2.0)*k2 ;
	math::Vector3d k4 = grad_yn + dt*k3 ;

	ynew = yn + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4) ;

	return ynew;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AdvectedPointRK4_2D::writeDiscretePathInVTK(){

	// First, we create the file where we are going to store the info
	std::ofstream stream= std::ofstream("GMDS_discrete_path.vtk", std::ios::out);
	//set the numerical precision (number of digits)
	stream.precision(15);
	//Header indicating which type of file it is
	//stream << "# Claire PHD Debug Version 1.0\n\n";

	stream << "# vtk DataFile Version 2.0\n";
	stream << "Test moche écriture VTK\n";
	stream << "ASCII\n\n";

	stream << "DATASET UNSTRUCTURED_GRID\n\n";

	stream << "POINTS ";
	stream << m_discrete_path.size() ;
	stream << " float\n";
	for (int i=0; i< m_discrete_path.size(); i++){
		stream << m_discrete_path[i].X() << " " << m_discrete_path[i].Y() << " " << m_discrete_path[i].Z() << "\n";
	}

	stream << "\n";

	stream << "POINT_DATA " << m_discrete_path.size() << "\n" ;
	stream << "SCALARS GMDS_discrete_path float 1\n";
	stream << "LOOKUP_TABLE default\n";
	for (int i=0; i< m_discrete_path.size(); i++){
		stream << 1.0 << "\n";
	}

	stream.close();
}
/*------------------------------------------------------------------------*/