/*----------------------------------------------------------------------------*/
#include "gmds/elgmorphing/ElgMorphing.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace elgmorphing {
/*----------------------------------------------------------------------------*/
ElgMorphing::ElgMorphing()
{

}
/*----------------------------------------------------------------------------*/
void ElgMorphing::computeBoundingBox(gmds::Mesh *AMesh,
												 std::string AGroupName,
												 gmds::TCoord minXYZ[3],
												 gmds::TCoord maxXYZ[3]) const
{
	minXYZ[0] =  HUGE_VALF;
	minXYZ[1] =  HUGE_VALF;
	minXYZ[2] =  HUGE_VALF;
	maxXYZ[0] = -HUGE_VALF;
	maxXYZ[1] = -HUGE_VALF;
	maxXYZ[2] = -HUGE_VALF;

	gmds::CellGroup<gmds::Node>* grp = AMesh->getGroup<gmds::Node>(AGroupName);
	for(auto nid: grp->cells()) {
		gmds::math::Point pt = AMesh->get<gmds::Node>(nid).point();

		if(minXYZ[0] > pt.X()) {
			minXYZ[0] = pt.X();
		}
		if(minXYZ[1] > pt.Y()) {
			minXYZ[1] = pt.Y();
		}
		if(minXYZ[2] > pt.Z()) {
			minXYZ[2] = pt.Z();
		}
		if(maxXYZ[0] < pt.X()) {
			maxXYZ[0] = pt.X();
		}
		if(maxXYZ[1] < pt.Y()) {
			maxXYZ[1] = pt.Y();
		}
		if(maxXYZ[2] < pt.Z()) {
			maxXYZ[2] = pt.Z();
		}
	}
}
/*----------------------------------------------------------------------------*/
gmds::math::Point ElgMorphing::bbtransform2d(gmds::math::Point APt,
														 gmds::TCoord minXYZ_orig[3], gmds::TCoord maxXYZ_orig[3],
														 gmds::TCoord minXYZ_dest[3], gmds::TCoord maxXYZ_dest[3]) const
{
	gmds::TCoord xt = ((APt.X() - minXYZ_orig[0]) / (maxXYZ_orig[0] - minXYZ_orig[0])) * (maxXYZ_dest[0] - minXYZ_dest[0]) + minXYZ_dest[0];
	gmds::TCoord yt = ((APt.Y() - minXYZ_orig[1]) / (maxXYZ_orig[1] - minXYZ_orig[1])) * (maxXYZ_dest[1] - minXYZ_dest[1]) + minXYZ_dest[1];
//	gmds::TCoord zt = ((APt.Z() - minXYZ_orig[2]) / (maxXYZ_orig[2] - minXYZ_orig[2])) * (maxXYZ_dest[2] - minXYZ_dest[2]) + minXYZ_dest[2];

//	return gmds::math::Point(xt,yt,zt);
	return gmds::math::Point(xt,yt,0.);
}
/*----------------------------------------------------------------------------*/
}  // namespace elgmorphing
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/