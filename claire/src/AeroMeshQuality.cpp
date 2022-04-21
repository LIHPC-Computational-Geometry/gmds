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

	double r1 = Utils::distFromNodeIds(AMesh, n0_id, n1_id)/Utils::distFromNodeIds(AMesh, n2_id, n3_id) ;
	if(r1 >= 1){
		r1 = 1.0/r1;
	}

	double r2 = Utils::distFromNodeIds(AMesh, n1_id, n2_id)/Utils::distFromNodeIds(AMesh, n3_id, n0_id) ;
	if(r2 >= 1){
		r2 = 1.0/r2;
	}

	return std::min(r1, r2);
}
/*------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
      /*----------------------------------------------------------------------------*/