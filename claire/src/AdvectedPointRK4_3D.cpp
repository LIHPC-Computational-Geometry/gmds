//
// Created by rochec on 24/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AdvectedPointRK4_3D.h>
#include <gmds/claire/Utils.h>
#include <limits>
#include <chrono>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
AdvectedPointRK4_3D::AdvectedPointRK4_3D(Mesh *AMesh, FastLocalize* A_fl, math::Point A_Pstart, double A_d0, Variable<double>* A_distance, Variable<math::Vector3d>* A_gradient2D) :
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
AdvectedPointRK4_3D::STATUS AdvectedPointRK4_3D::execute()
{
	double minLenght = minEdgeLenght();

	double dt = 0.9*minLenght;
	double err = pow(10,-6);
	int max_iterations=10000;

	int iterations=0;
	double dist(0);
	math::Vector3d Grad;
	Eigen::Matrix4d Mat_A_Inv;

	// Initialisation
	TCellID region_id = inWhichTetra(m_Pstart) ;		// Dans quel triangle est le point de départ

	if (region_id == NullID)
	{
		std::cout << "Tetra pas trouvé, point: " << m_Pstart << std::endl;
		//exit(1);
		// The node considered is not in the domain. Then, we replace it by the closest node of the tetra mesh.
		gmds::Cell::Data data = m_fl->find(m_Pstart);
		TCellID n_closest_id_fl = data.id;
		Node n_closest_fl = m_mesh->get<Node>(n_closest_id_fl);
		m_Pstart = n_closest_fl.point();
		m_Pend = n_closest_fl.point();
		if (region_id==NullID && n_closest_fl.get<Region>().size() == 0 )
		{
			std::cout << "ATTENTION AdvectedPointRK4_3D: Point de départ hors du domaine." << std::endl;
		}
		Region r = n_closest_fl.get<Region>()[0];
		region_id = r.id();
	}

	Mat_A_Inv = getInvMatrixA(region_id);
	dist = interpolationDistance(region_id, Mat_A_Inv, m_Pstart);	// A quelle distance est le point de départ
	Grad = interpolationGradient(region_id, Mat_A_Inv, m_Pstart);	// Quel est le gradient à ce point

	if ( abs(Grad.X()) <= pow(10,-6) && abs(Grad.Y())<=pow(10,-6) && abs(Grad.Z())<=pow(10,-6))
	{
		Region r_test = m_mesh->get<Region>(region_id);
		for (auto n_test:r_test.get<Node>())
		{
			std::cout << "Node " << n_test << ", distance: " << m_distance->value(n_test.id()) << std::endl;
		}

		gmds::Cell::Data data = m_fl->find(m_Pstart);
		Node n = m_mesh->get<Node>(data.id);
		std::vector<Edge> edges = n.get<Edge>();
		for (auto e:edges)
		{
			Node n_opp = e.getOppositeNode(n);
			Grad = Grad + m_gradient2D->value(n_opp.id());
		}
		Grad = Grad/edges.size();
	}

	if ( abs(Grad.X()) <= pow(10,-6) && abs(Grad.Y())<=pow(10,-6) && abs(Grad.Z())<=pow(10,-6))
	{
		std::cout << "ATTENTION AdvectedPointRK4_3D: Starting gradient vector equal to 0." << std::endl;
	}
	if (dist > m_d0)
	{
		std::cout << "ATTENTION AdvectedPointRK4_3D: Distance de départ supérieure à la distance cible." << std::endl;
		std::cout << "Starting point: " << m_Pstart << std::endl;
		std::cout << "Tetra : " << region_id << std::endl;
		std::cout << "Point 1 : " << m_mesh->get<Region>(region_id).get<Node>()[0].point() << std::endl;
		std::cout << "Point 2 : " << m_mesh->get<Region>(region_id).get<Node>()[1].point() << std::endl;
		std::cout << "Point 3 : " << m_mesh->get<Region>(region_id).get<Node>()[2].point() << std::endl;
		std::cout << "Point 4 : " << m_mesh->get<Region>(region_id).get<Node>()[3].point() << std::endl;
		std::cout << "Distance cible: " << m_d0 << ", distance de départ: " << dist << std::endl;
		exit(1);
	}

	while ( (abs(dist-m_d0) > err) && iterations < max_iterations ) {
		math::Point M = RungeKutta4(m_Pend, Grad.normalize(), dt);	// Calcule la position du point à l'itération n+1 avec un RK4
		// On vérifie ensuite si cette position est "valide"
		region_id = inWhichTetra(M, region_id) ;

		// Si le noeud M calculé est bien dans le maillage,
		// alors on regarde la distance et le gradient.
		// Sinon, on retranche le pas de temps et on recommence.
		if(region_id != NullID) {
			Mat_A_Inv = getInvMatrixA(region_id);		// Calcul la matrice inverse à A pour résoudre le système Ax=b
			double dist_M = interpolationDistance(region_id, Mat_A_Inv, M);	// A quelle distance est le point M

			// Si la distance du point M est plus petite que celle souhaitée,
			// alors on met à jour la valeur du gradient et on continue le RK4
			// Sinon, on retranche le pas de temps et on recommence
			if ((dist_M < m_d0)) {
				m_Pend = M;
				Grad = interpolationGradient(region_id, Mat_A_Inv, M);	// Mise à jour du gradient
				dist = dist_M;
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
	region_id = inWhichTetra(m_Pend, region_id) ;
	Mat_A_Inv = getInvMatrixA(region_id);
	dist = interpolationDistance(region_id, Mat_A_Inv, m_Pend);

	//writeDiscretePathInVTK();

	return AdvectedPointRK4_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
bool AdvectedPointRK4_3D::isInTetra(TCellID region_id, math::Point M){

	Region r = m_mesh->get<Region>(region_id);
	std::vector<Node> Tetra_nodes = r.get<Node>();
	return math::Utils::isInTetra(Tetra_nodes[0].point(), Tetra_nodes[1].point(), Tetra_nodes[2].point(), Tetra_nodes[3].point(), M);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
TCellID AdvectedPointRK4_3D::inWhichTetra(math::Point M, TCellID r0_id){
	TCellID region_id;
	/*
	bool isInRegion(false);

	// Si un r0_id a été donné en entrée, on regarde dans les tétras voisins à
	// celui ci si M y est, avant de regarder sur la totalité du maillage
	if (r0_id != NullID){
		// On regarde quelles sont les régions voisines de la région r0_id
		std::vector<TCellID> adjacent_regions;
		Region r0 = m_mesh->get<Region>(r0_id);
		std::vector<Node> region_nodes = r0.get<Node>();
		for (auto n:region_nodes){
			std::vector<Region> node_regions = n.get<Region>();
			for (auto r1:node_regions){
				bool alreadyinvector(false);
				for (auto rv_id:adjacent_regions){
					if(r1.id() == rv_id){
						alreadyinvector = true;
					}
				}
				if (!alreadyinvector) {
					adjacent_regions.push_back(r1.id());
				}
			}
		}
		// On regarde si le point M est dans un des tétra
		for (auto r_adj_id:adjacent_regions){
			if(!isInRegion){
				isInRegion = isInTetra(r_adj_id, M);
	         if(isInRegion){
					region_id = r_adj_id;
				}
			}
		}
	}

	// Use FastLocalize to check the tetras around the closest node
	// to the point M
	if (!isInRegion)
	{
		gmds::Cell::Data data = m_fl->find(M);
		TCellID n_closest_id = data.id;
		Node n_closest = m_mesh->get<Node>(n_closest_id);
		std::vector<Region> n_closest_tetras = n_closest.get<Region>();
		//std::cout << "Dim of the closest Cell: " << data.dim << std::endl;
		//std::cout << "Nbr closest tetras: " << n_closest_tetras.size() << std::endl;
		for (auto r:n_closest_tetras)
		{
			if (!isInRegion) {
				isInRegion = isInTetra(r.id(), M);
				if (isInRegion) {
					region_id = r.id();
				}
			}
		}
	}

	if (!isInRegion){
		region_id = NullID;
		//std::cout << "Attention: APRK4, region not found." << std::endl;
	}


	return region_id;
	*/
	return m_fl->findTetra(M);
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AdvectedPointRK4_3D::minEdgeLenght(){
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
Eigen::Matrix4d AdvectedPointRK4_3D::getInvMatrixA(TCellID region_id){

	Eigen::Matrix4d Mat_A;

	Region r = m_mesh->get<Region>(region_id);
	std::vector<Node> nodes = r.get<Node>() ;
	math::Point p0 = nodes[0].point() ;
	math::Point p1 = nodes[1].point() ;
	math::Point p2 = nodes[2].point() ;
	math::Point p3 = nodes[3].point() ;

	TCellID n0_id = nodes[0].id();
	TCellID n1_id = nodes[1].id();
	TCellID n2_id = nodes[2].id();
	TCellID n3_id = nodes[3].id();

	// x1 y1 z1 1
	// x2 y2 z2 1
	// x3 y3 z3 1
	// x4 y4 z4 1
	// Remplissage de la matrice A du système Ax=b
	Mat_A(0,0) = p0.X() ;
	Mat_A(0,1) = p0.Y() ;
	Mat_A(0,2) = p0.Z() ;
	Mat_A(0,3) = 1.0 ;
	Mat_A(1,0) = p1.X() ;
	Mat_A(1,1) = p1.Y() ;
	Mat_A(1,2) = p1.Z() ;
	Mat_A(1,3) = 1.0 ;
	Mat_A(2,0) = p2.X() ;
	Mat_A(2,1) = p2.Y() ;
	Mat_A(2,2) = p2.Z() ;
	Mat_A(2,3) = 1.0 ;
	Mat_A(3,0) = p3.X() ;
	Mat_A(3,1) = p3.Y() ;
	Mat_A(3,2) = p3.Z() ;
	Mat_A(3,3) = 1.0 ;

	return Mat_A.inverse();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AdvectedPointRK4_3D::interpolationDistance(TCellID region_id, Eigen::Matrix4d Mat_A_Inv, math::Point M){
	double dist_interpolee;

	Region r = m_mesh->get<Region>(region_id);
	std::vector<Node> nodes = r.get<Node>() ;

	TCellID n0_id = nodes[0].id();
	TCellID n1_id = nodes[1].id();
	TCellID n2_id = nodes[2].id();
	TCellID n3_id = nodes[3].id();

	// Initialisation du vecteur des distances aux noeuds
	Eigen::Vector4d b;
	b[0] = m_distance->value(n0_id)  ;
	b[1] = m_distance->value(n1_id)  ;
	b[2] = m_distance->value(n2_id)  ;
	b[3] = m_distance->value(n3_id)  ;

	Eigen::Vector4d coef = Mat_A_Inv * b;
	dist_interpolee = coef[0]*M.X() + coef[1]*M.Y() + coef[2]*M.Z() + coef[3] ;

	return dist_interpolee;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d AdvectedPointRK4_3D::interpolationGradient(TCellID region_id, Eigen::Matrix4d Mat_A_Inv, math::Point M){
	math::Vector3d Gradient_M;

	Region r = m_mesh->get<Region>(region_id);
	std::vector<Node> nodes = r.get<Node>() ;

	TCellID n0_id = nodes[0].id();
	TCellID n1_id = nodes[1].id();
	TCellID n2_id = nodes[2].id();
	TCellID n3_id = nodes[3].id();

	// Initialisation du vecteur des distances aux noeuds
	Eigen::Vector4d bx, by, bz;
	//Initialisation du vecteur des composantes du grad dans la direction x
	bx[0] = m_gradient2D->value(n0_id).X()  ;
	bx[1] = m_gradient2D->value(n1_id).X()  ;
	bx[2] = m_gradient2D->value(n2_id).X()  ;
	bx[3] = m_gradient2D->value(n3_id).X()  ;

	//Initialisation du vecteur des composantes du grad dans la direction y
	by[0] = m_gradient2D->value(n0_id).Y()  ;
	by[1] = m_gradient2D->value(n1_id).Y()  ;
	by[2] = m_gradient2D->value(n2_id).Y()  ;
	by[3] = m_gradient2D->value(n3_id).Y()  ;

	//Initialisation du vecteur des composantes du grad dans la direction z
	bz[0] = m_gradient2D->value(n0_id).Z()  ;
	bz[1] = m_gradient2D->value(n1_id).Z()  ;
	bz[2] = m_gradient2D->value(n2_id).Z()  ;
	bz[3] = m_gradient2D->value(n3_id).Z()  ;

	// Définition du gradient exact dans la direction x au point M
	Eigen::Vector4d coeff_grad_x = Mat_A_Inv * bx;
	Gradient_M.setX(coeff_grad_x[0]*M.X() + coeff_grad_x[1]*M.Y() + coeff_grad_x[2]*M.Z() + coeff_grad_x[3]) ;

	// Définition du gradient exact dans la direction y au point M
	Eigen::Vector4d coeff_grad_y = Mat_A_Inv * by;
	Gradient_M.setY(coeff_grad_y[0]*M.X() + coeff_grad_y[1]*M.Y() + coeff_grad_y[2]*M.Z() + coeff_grad_y[3]) ;

	// Définition du gradient exact dans la direction z au point M
	Eigen::Vector4d coeff_grad_z = Mat_A_Inv * bz;
	Gradient_M.setZ(coeff_grad_z[0]*M.X() + coeff_grad_z[1]*M.Y() + coeff_grad_z[2]*M.Z() + coeff_grad_z[3]) ;

	return Gradient_M;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Point AdvectedPointRK4_3D::RungeKutta4(math::Point yn, math::Vector3d grad_yn, double dt){
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
void AdvectedPointRK4_3D::writeDiscretePathInVTK(){

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


/*-------------------------------------------------------------------*/
math::Point AdvectedPointRK4_3D::getPend(){
	return m_Pend;
}
/*-------------------------------------------------------------------*/