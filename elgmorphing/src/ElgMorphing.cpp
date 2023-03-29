/*----------------------------------------------------------------------------*/
#include "gmds/elgmorphing/ElgMorphing.h"

#include <algorithm>

#include <gmds/smoothy/EllipticSmoother2D.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace elgmorphing {
/*----------------------------------------------------------------------------*/
ElgMorphing::ElgMorphing()
{

}
/*----------------------------------------------------------------------------*/
ElgMorphing::STATUS ElgMorphing::execute()
{
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================
//	gmds::smoothy::EllipticSmoother2D smoother2D(mesh_orig_);
//	smoother2D.lock(mark_fixed_nodes_);
//	smoother2D.execute();

	return ElgMorphing::SUCCESS;
}
/*----------------------------------------------------------------------------*/
std::set<TCellID> ElgMorphing::compoundGrpNodes(gmds::Mesh *AMesh,
																std::string AString)
{
	std::set<TCellID> nodeset;

	const int nbgrp = AMesh->getNbGroups<gmds::Face>();
	for(int igrp=0; igrp<nbgrp; igrp++) {
		gmds::CellGroup<gmds::Face>* grp = AMesh->getGroup<gmds::Face>(igrp);

		// check whether AString is a suffix of the group name
		if(0 == grp->name().find(AString)) {

			for (auto cid : grp->cells()) {
				gmds::Face c = AMesh->get<gmds::Face>(cid);
				std::vector<gmds::TCellID> nids = c.getIDs<gmds::Node>();
				for(auto n: nids) {
					nodeset.insert(n);
				}
			}
		}
	}

	return nodeset;
}
/*----------------------------------------------------------------------------*/
std::set<TCellID> ElgMorphing::compoundGrpBnd(gmds::Mesh *AMesh,
															 std::string AString)
{
	std::set<TCellID> nodeset;
	std::set<TCellID> cellset;
	const int nbgrp = AMesh->getNbGroups<gmds::Face>();
	for(int igrp=0; igrp<nbgrp; igrp++) {
		gmds::CellGroup<gmds::Face>* grp = AMesh->getGroup<gmds::Face>(igrp);

		// check whether AString is a suffix of the group name
		if(0 == grp->name().find(AString)) {

			for (auto cid : grp->cells()) {
				cellset.insert(cid);
			}
		}
	}

	for(auto cid: cellset) {
		std::vector<gmds::Edge> edges = AMesh->get<gmds::Face>(cid).get<gmds::Edge>();

		for(auto e: edges) {
			std::vector<gmds::TCellID> cids = e.getIDs<gmds::Face>();
			if (1 == cids.size()) {
				std::vector<gmds::TCellID> nids = e.getIDs<gmds::Node>();
				nodeset.insert(nids[0]);
				nodeset.insert(nids[1]);
			}
			else {
				if ((cellset.find(cids[0]) == cellset.end()) || (cellset.find(cids[1]) == cellset.end())) {
					std::vector<gmds::TCellID> nids = e.getIDs<gmds::Node>();
					nodeset.insert(nids[0]);
					nodeset.insert(nids[1]);
				}
			}
		}
	}

	return nodeset;
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

		minXYZ[0] = std::min(minXYZ[0],pt.X());
		minXYZ[1] = std::min(minXYZ[1],pt.Y());
		minXYZ[2] = std::min(minXYZ[2],pt.Z());
		maxXYZ[0] = std::max(maxXYZ[0],pt.X());
		maxXYZ[1] = std::max(maxXYZ[1],pt.Y());
		maxXYZ[2] = std::max(maxXYZ[2],pt.Z());
	}
}
/*----------------------------------------------------------------------------*/
void ElgMorphing::computeBoundingBox(const gmds::Mesh *AMesh,
												 const std::set<TCellID>& AIDs,
												 gmds::TCoord minXYZ[3],
												 gmds::TCoord maxXYZ[3]) const
{
	minXYZ[0] =  HUGE_VALF;
	minXYZ[1] =  HUGE_VALF;
	minXYZ[2] =  HUGE_VALF;
	maxXYZ[0] = -HUGE_VALF;
	maxXYZ[1] = -HUGE_VALF;
	maxXYZ[2] = -HUGE_VALF;

	for(auto nid: AIDs) {
		gmds::math::Point pt = AMesh->get<gmds::Node>(nid).point();
		minXYZ[0] = std::min(minXYZ[0],pt.X());
		minXYZ[1] = std::min(minXYZ[1],pt.Y());
		minXYZ[2] = std::min(minXYZ[2],pt.Z());
		maxXYZ[0] = std::max(maxXYZ[0],pt.X());
		maxXYZ[1] = std::max(maxXYZ[1],pt.Y());
		maxXYZ[2] = std::max(maxXYZ[2],pt.Z());
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