//
// Created by rochec on 01/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/PointFollowingVectorFieldOnNodes.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PointFollowingVectorFieldOnNodes::PointFollowingVectorFieldOnNodes(Mesh *AMesh, math::Point A_Pstart, double A_d0,
                                                                   Variable<double> *A_distance,
                                                                   Variable<math::Vector3d>* A_gradient) {
	m_mesh = AMesh;
	m_Pstart = A_Pstart;
	m_d0 = A_d0;
	m_distance = A_distance;
	m_gradient = A_gradient;
	m_discrete_path.push_back(m_Pstart);
}


/*------------------------------------------------------------------------*/
PointFollowingVectorFieldOnNodes::STATUS PointFollowingVectorFieldOnNodes::execute()
{
	m_Pend = m_Pstart;
	double minLenght = minEdgeLenght();
	//std::cout << "min : " << minLenght << std::endl;

	TCellID closest_node_id = closestNode(m_Pend);
	double d_closest_node = m_distance->value(closest_node_id);

	while (d_closest_node < m_d0){
		math::Vector3d Grad_local = m_gradient->value(closest_node_id) ;
		m_Pend = m_Pend + minLenght*Grad_local.normalize();
		closest_node_id = closestNode(m_Pend);
		d_closest_node = m_distance->value(closest_node_id);

		//std::cout << "Point intermédiaire : " << m_Pend << std::endl;
		//std::cout << "Distance : " << d_closest_node << std::endl;
		m_discrete_path.push_back(m_Pend);
	}

	//std::cout << "Point final : " << m_Pend << std::endl;

	// Ecriture du chemin suivi
	writeDiscretePathInVTK();

	return PointFollowingVectorFieldOnNodes::SUCCESS;
}
/*------------------------------------------------------------------------*/



/*------------------------------------------------------------------------*/
TCellID PointFollowingVectorFieldOnNodes::closestNode(math::Point M){
	TCellID closest_node_id;
	double min_norme(std::numeric_limits<double>::max());
	for (auto n_id:m_mesh->nodes()){
		Node n = m_mesh->get<Node>(n_id);
		math::Vector3d Vec = M-n.point() ;
		double norme = Vec.norm() ;
		if(norme < min_norme){
			min_norme = norme;
			closest_node_id = n_id;
		}
	}
	return closest_node_id;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
double PointFollowingVectorFieldOnNodes::minEdgeLenght(){
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
void PointFollowingVectorFieldOnNodes::writeDiscretePathInVTK(){

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