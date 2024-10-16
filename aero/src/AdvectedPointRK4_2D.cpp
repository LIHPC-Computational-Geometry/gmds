//
// Created by rochec on 04/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/aero/AdvectedPointRK4_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
AdvectedPointRK4_2D::AdvectedPointRK4_2D(Mesh *AMesh, FastLocalize *A_fl, const math::Point& A_Pstart, double A_d0, Variable<double>* A_distance, Variable<math::Vector3d>* A_gradient2D) :
	m_mesh(AMesh),
  	m_fl(A_fl)
{
	m_Pstart = A_Pstart;
	m_Pend = A_Pstart;
	m_d0 = A_d0;
	m_distance = A_distance;
	m_gradient2D = A_gradient2D;
	m_discrete_path.push_back({m_Pstart.X(),m_Pstart.Y(),m_Pstart.Z()});
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AdvectedPointRK4_2D::STATUS AdvectedPointRK4_2D::execute()
{
	double minLenght = minEdgeLenght();

	double dt = 0.9*minLenght;
	double err = pow(10,-6);
	int max_iterations=10000;

	int iterations=0;
	double dist(0);
	math::Vector3d Grad;
	Eigen::Matrix3d Mat_A_Inv;

	// Initialisation
	TCellID face_id = inWhichTriangle(m_Pstart, NullID) ;		// Dans quel triangle est le point de départ
	//std::cout << "Face id " << face_id << std::endl;
	if (face_id==NullID)
	{
		gmds::Cell::Data data = m_fl->find(m_Pstart);
		TCellID n_closest_id = data.id;
		Node n_closest = m_mesh->get<Node>(n_closest_id);
		m_Pstart = n_closest.point() ;
		m_Pend = n_closest.point() ;
		Face f = n_closest.get<Face>()[0] ;
		face_id = f.id() ;
		std::cout << "Face id " << face_id << std::endl;
	}

	Mat_A_Inv = getInvMatrixA(face_id);
	//std::cout << "Mat inversée" << std::endl;
	dist = interpolationDistance(face_id, Mat_A_Inv, m_Pstart);	// A quelle distance est le point de départ
	//std::cout << "dist " << dist << std::endl;
	Grad = interpolationGradient(face_id, Mat_A_Inv, m_Pstart);	// Quel est le gradient à ce point

	//std::cout << "grad " << Grad << std::endl;

	while ( (abs(dist-m_d0) > err) && iterations < max_iterations ) {
		math::Point M = RungeKutta4(m_Pend, Grad.normalize(), dt);	// Calcule la position du point à l'itération n+1 avec un RK4
		// On vérifie ensuite si cette position est "valide"
		face_id = inWhichTriangle(M, face_id) ;

		// Si le noeud M calculé est bien dans le maillage,
		// alors on regarde la distance et le gradient.
		// Sinon, on retranche le pas de temps et on recommence.
		if(face_id != NullID) {
			Mat_A_Inv = getInvMatrixA(face_id);		// Calcul la matrice inverse à A pour résoudre le système Ax=b
			double dist_M = interpolationDistance(face_id, Mat_A_Inv, M);	// A quelle distance est le point M

			// Si la distance du point M est plus petite que celle souhaitée,
			// alors on met à jour la valeur du gradient et on continue le RK4
			// Sinon, on retranche le pas de temps et on recommence
			if ((dist_M < m_d0)) {
				//std::cout << "Triangle : " << face_id << std::endl;
				m_Pend = M;
				Grad = interpolationGradient(face_id, Mat_A_Inv, M);	// Mise à jour du gradient
				dist = dist_M;
				//std::cout << "distance : " << dist << std::endl;
				m_discrete_path.push_back({m_Pend.X(),m_Pend.Y(),m_Pend.Z()});
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

	// Pour le noeud Pend, calcule de la distance finale.
	// Normalement, il n'y en a pas besoin, mais c'est pour l'affichage.
	face_id = inWhichTriangle(m_Pend, face_id) ;
	Mat_A_Inv = getInvMatrixA(face_id);
	dist = interpolationDistance(face_id, Mat_A_Inv, m_Pend);

	writeDiscretePathInVTK();

	return AdvectedPointRK4_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
bool AdvectedPointRK4_2D::isInTriangle(TCellID face_id, math::Point &M){
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
TCellID AdvectedPointRK4_2D::inWhichTriangle(math::Point &M, TCellID f0_id){
	TCellID face_id;
	bool isInFace(false);

	// Si un f0_id a été donné en entrée, on regarde dans les triangles voisins à
	// celui ci si M y est, avant de regarder sur la totalité du maillage
	if (f0_id != NullID){
		// On regarde quels sont les faces voisines de la face r0_id
		std::vector<TCellID> adjacent_faces;
		Face f0 = m_mesh->get<Face>(f0_id);
		std::vector<Node> face_nodes = f0.get<Node>();
		for (auto const& n:face_nodes){
			std::vector<Face> node_faces = n.get<Face>();
			for (auto const& f1:node_faces){
				bool alreadyinvector(false);
				for (auto const& fv_id:adjacent_faces){
					if(f1.id() == fv_id){
						alreadyinvector = true;
					}
				}
				if (!alreadyinvector) {
					adjacent_faces.push_back(f1.id());
				}
			}
		}
		// On regarde si le point M est dans un des triangles
		for (auto const& f_adj_id:adjacent_faces){
			if(!isInFace){
				isInFace = isInTriangle(f_adj_id, M);
				if(isInFace){
					face_id = f_adj_id;
				}
			}
		}
	}


	// Use FastLocalize to check the tetras around the closest node
	// to the point M
	if (!isInFace)
	{
		gmds::Cell::Data data = m_fl->find(M);
		TCellID n_closest_id = data.id;
		Node n_closest = m_mesh->get<Node>(n_closest_id);
		std::vector<Face> n_closest_tri = n_closest.get<Face>();
		for (auto const& f:n_closest_tri)
		{
			if (!isInFace) {
				isInFace = isInTriangle(f.id(), M);
				if (isInFace) {
					face_id = f.id();
				}
			}
		}
	}

	/*
	for(auto f_it = m_mesh->faces_begin(); f_it!= m_mesh->faces_end() && !isInFace;++f_it){
		TCellID f_id = *f_it;
		isInFace = isInTriangle(f_id, M);
		if(isInFace){
			face_id = f_id;
		}
	}
	 */


	if (!isInFace){
		face_id = NullID;
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
	for (auto const& edge_id:m_mesh->edges()){
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
Eigen::Matrix3d AdvectedPointRK4_2D::getInvMatrixA(TCellID face_id){

	Eigen::Matrix3d Mat_A;

	Face f = m_mesh->get<Face>(face_id);
	std::vector<Node> nodes = f.get<Node>() ;
	math::Point p0 = nodes[0].point() ;
	math::Point p1 = nodes[1].point() ;
	math::Point p2 = nodes[2].point() ;

	// x1 y1 1
	// x2 y2 1
	// x3 y3 1
	// Remplissage de la matrice A du système Ax=b
	Mat_A(0,0) = p0.X() ;
	Mat_A(0,1) = p0.Y() ;
	Mat_A(0,2) = 1.0 ;
	Mat_A(1,0) = p1.X() ;
	Mat_A(1,1) = p1.Y() ;
	Mat_A(1,2) = 1.0 ;
	Mat_A(2,0) = p2.X() ;
	Mat_A(2,1) = p2.Y() ;
	Mat_A(2,2) = 1.0 ;

	return Mat_A.inverse();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AdvectedPointRK4_2D::interpolationDistance(TCellID face_id, Eigen::Matrix3d &Mat_A_Inv, math::Point M){
	double dist_interpolee;

	Face f = m_mesh->get<Face>(face_id);
	std::vector<Node> nodes = f.get<Node>() ;

	TCellID n0_id = nodes[0].id();
	TCellID n1_id = nodes[1].id();
	TCellID n2_id = nodes[2].id();

	// Initialisation du vecteur des distances aux noeuds
	Eigen::Vector3d b;
	b[0] = m_distance->value(n0_id)  ;
	b[1] = m_distance->value(n1_id)  ;
	b[2] = m_distance->value(n2_id)  ;

	Eigen::Vector3d coef = Mat_A_Inv * b;
	dist_interpolee = coef[0]*M.X() + coef[1]*M.Y() + coef[2] ;

	return dist_interpolee;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d AdvectedPointRK4_2D::interpolationGradient(TCellID face_id, Eigen::Matrix3d &Mat_A_Inv, math::Point M){
	math::Vector3d Gradient_M;

	Face f = m_mesh->get<Face>(face_id);
	std::vector<Node> nodes = f.get<Node>() ;

	TCellID n0_id = nodes[0].id();
	TCellID n1_id = nodes[1].id();
	TCellID n2_id = nodes[2].id();

	// Initialisation du vecteur des distances aux noeuds
	Eigen::Vector3d bx, by;
	//Initialisation du vecteur des composantes du grad dans la direction x
	bx[0] = m_gradient2D->value(n0_id).X()  ;
	bx[1] = m_gradient2D->value(n1_id).X()  ;
	bx[2] = m_gradient2D->value(n2_id).X()  ;

	//Initialisation du vecteur des composantes du grad dans la direction y
	by[0] = m_gradient2D->value(n0_id).Y()  ;
	by[1] = m_gradient2D->value(n1_id).Y()  ;
	by[2] = m_gradient2D->value(n2_id).Y()  ;

	// Définition du gradient exact dans la direction x au point M
	Eigen::Vector3d coeff_grad_x = Mat_A_Inv * bx;
	Gradient_M.setX(coeff_grad_x[0]*M.X() + coeff_grad_x[1]*M.Y() + coeff_grad_x[2]) ;

	// Définition du gradient exact dans la direction y au point M
	Eigen::Vector3d coeff_grad_y = Mat_A_Inv * by;
	Gradient_M.setY(coeff_grad_y[0]*M.X() + coeff_grad_y[1]*M.Y() + coeff_grad_y[2]) ;

	// Définition du gradient exact dans la direction z au point M
	Gradient_M.setZ(0.0);

	return Gradient_M;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Point AdvectedPointRK4_2D::RungeKutta4(math::Point &yn, math::Vector3d grad_yn, double dt){
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
	for (auto const& point_local:m_discrete_path){
		stream << point_local.X() << " " << point_local.Y() << " " << point_local.Z() << "\n";
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


/*-------------------------------------------------------------------*/
math::Point AdvectedPointRK4_2D::getPend(){
	return m_Pend;
}
/*-------------------------------------------------------------------*/