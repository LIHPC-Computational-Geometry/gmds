//
// Created by rochec on 14/04/2022.
//

/*----------------------------------------------------------------------------*/
#include <gmds/claire/AeroMeshQuality.h>
#include <gmds/claire/Utils.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {


/*------------------------------------------------------------------------*/
double AeroMeshQuality::oppositeedgeslenghtratio(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id){

	//std::cout << "ids : " << n0_id << ", " << n1_id << ", " << n2_id << ", " << n3_id << std::endl;

	double r1 = Utils::distFromNodeIds(AMesh, n0_id, n1_id)/Utils::distFromNodeIds(AMesh, n2_id, n3_id) ;
	//if(r1 >= 1){
	//	r1 = 1.0/r1;
	//}

	/*
	double r2 = Utils::distFromNodeIds(AMesh, n1_id, n2_id)/Utils::distFromNodeIds(AMesh, n3_id, n0_id) ;
	if(r2 >= 1){
		r2 = 1.0/r2;
	}

	std::cout << "R1 et R2 " << r1 << " " << r2 << std::endl;
	std::cout << "--------------------------" << std::endl;
	 */

	//return std::min(r1, r2);
	return r1;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::angleouverture(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id){

	//std::cout << "ids : " << n0_id << ", " << n1_id << ", " << n2_id << ", " << n3_id << std::endl;

	Node n0 = AMesh->get<Node>(n0_id);
	Node n1 = AMesh->get<Node>(n1_id);
	Node n2 = AMesh->get<Node>(n2_id);
	Node n3 = AMesh->get<Node>(n3_id);

	math::Point p0 = n0.point();
	math::Point p1 = n1.point();
	math::Point p2 = n2.point();
	math::Point p3 = n3.point();

	Vector3d v1 = p1-p0;
	Vector3d v2 = p2-p0;
	Vector3d v3 = p3-p0;

	v1 = v1.normalize();
	v2 = v2.normalize();
	v3 = v3.normalize();

	double angle1 = acos(v1.dot(v2)) ;

	double angle2 = acos(v1.dot(v3)) ;

	return angle1+angle2;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::minlenghtedge(Mesh *AMesh, TCellID n0_id, TCellID n1_id, TCellID n2_id, TCellID n3_id){

	//std::cout << "ids : " << n0_id << ", " << n1_id << ", " << n2_id << ", " << n3_id << std::endl;

	Node n0 = AMesh->get<Node>(n0_id);
	Node n1 = AMesh->get<Node>(n1_id);
	Node n2 = AMesh->get<Node>(n2_id);
	Node n3 = AMesh->get<Node>(n3_id);

	math::Point p0 = n0.point();
	math::Point p1 = n1.point();
	math::Point p2 = n2.point();
	math::Point p3 = n3.point();

	Vector3d v1 = p1-p0;
	Vector3d v2 = p1-p2;
	Vector3d v3 = p2-p3;
	Vector3d v4 = p3-p0;

	return std::min(std::min(v1.norm(), v2.norm()), std::min(v3.norm(), v4.norm()));
}
/*------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
      /*----------------------------------------------------------------------------*/