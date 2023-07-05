/*------------------------------------------------------------------------*/
#include "gmds/elgmorphing/ElgMorphing.h"

#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/LimaReader.h>
#include <gmds/io/MdlReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/smoothy/EllipticSmoother2D.h>
#include <iostream>
#include <sstream>
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	std::cout << "=== ELGMORPHING ALGO ====" << std::endl;

	std::string filename_mli;
	std::string filename_mdl;
	std::string matname_orig;
	std::string matname_dest;
	std::istringstream iss_1(argv[1]);
	iss_1 >> filename_mli;
	std::istringstream iss_2(argv[2]);
	iss_2 >> filename_mdl;
	std::istringstream iss_3(argv[3]);
	iss_3 >> matname_orig;
	std::istringstream iss_4(argv[4]);
	iss_4 >> matname_dest;

	std::cout<<"source mesh "<<filename_mli<<std::endl;
	std::cout<<"target mesh "<<filename_mdl<<std::endl;
	std::cout<<"mat_orig    "<<matname_orig<<std::endl;
	std::cout<<"mat_dest    "<<matname_dest<<std::endl;

	gmds::Mesh m_dest(gmds::MeshModel(gmds::DIM2|gmds::E|gmds::N|gmds::E2N));
	gmds::MdlReader reader_dest(m_dest, matname_dest);
	reader_dest.read(filename_mdl);
	reader_dest.createVariablesFromGroups();

	gmds::elgmorphing::ElgMorphing elg;
	gmds::TCoord minXYZ_dest[3];
	gmds::TCoord maxXYZ_dest[3];
	elg.computeBoundingBox(&m_dest, matname_dest, minXYZ_dest, maxXYZ_dest);

	gmds::Mesh m_orig(gmds::MeshModel(gmds::DIM2 | gmds::F | gmds::E | gmds::N | gmds::F2N | gmds::F2E | gmds::E2F | gmds::E2N | gmds::N2E | gmds::N2F));
	gmds::LimaReader reader(m_orig);
	reader.read(filename_mli,gmds::F|gmds::N);

	gmds::MeshDoctor doc(&m_orig);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();
	for(auto f_id:m_orig.faces()) {
		gmds::Face f=m_orig.get<gmds::Face>(f_id);
		if (f.normal().dot(gmds::math::Vector3d({.0, .0, 1.0})) <= 0) {
			std::vector<gmds::TCellID> ns = f.getIDs<gmds::Node>();
			std::vector<gmds::TCellID> ns2(ns.size());
			for (auto i = 0; i < ns.size(); i++)
				ns2[ns.size() - 1 - i] = ns[i];
			f.set<gmds::Node>(ns2);
		}
	}

	//==================================================================
	// MARK ALL THE BOUNDARY CELL OF THE  MESH
	//==================================================================
	auto nb_locked = 0;
	auto mark_bnd_nodes = m_orig.newMark<gmds::Node>();
	gmds::math::Point origin({0,0,0});
	gmds::Variable<int>* var_bnd = m_orig.newVariable<int,gmds::GMDS_NODE>("bnd");
	std::set<gmds::TCellID> mat_nids = elg.compoundGrpNodes(&m_orig, matname_orig);
	std::set<gmds::TCellID> mat_bnd = elg.compoundGrpBnd(&m_orig, matname_orig);
	std::cout<<"nbnodes    "<<matname_orig<<" "<<mat_nids.size()<<std::endl;
	std::cout<<"nbnodesBnd "<<matname_orig<<" "<<mat_bnd.size()<<std::endl;

	gmds::BoundaryOperator2D bo(&m_orig);
	std::vector<gmds::TCellID> bnd_nids;
	bo.getBoundaryNodes(bnd_nids);
	std::cout<<"bndnodes "<<bnd_nids.size()<<std::endl;

	for(auto nid:bnd_nids) {
		nb_locked += 1;
		m_orig.mark<gmds::Node>(nid, mark_bnd_nodes);
		var_bnd->set(nid, 1);
	}
	for(auto nid:mat_bnd) {
		nb_locked += 1;
		m_orig.mark<gmds::Node>(nid, mark_bnd_nodes);
		var_bnd->set(nid, 1);
	}

	bool debug_mode = true;
	gmds::IGMeshIOService ioService(&m_orig);
	gmds::VTKWriter w(&ioService);
	w.setCellOptions(gmds::N|gmds::F);
	w.setDataOptions(gmds::N|gmds::F);
	if(debug_mode)
		w.write("elgmorphing_marked.vtk");

	gmds::TCoord minXYZ_orig[3];
	gmds::TCoord maxXYZ_orig[3];
	elg.computeBoundingBox(&m_orig, mat_nids, minXYZ_orig, maxXYZ_orig);
	for(auto nid: mat_bnd) {
		gmds::math::Point pt = elg.bbtransform2d(m_orig.get<gmds::Node>(nid).point(), minXYZ_orig, maxXYZ_orig, minXYZ_dest, maxXYZ_dest);
		m_orig.get<gmds::Node>(nid).setPoint(pt);
	}

	if(debug_mode)
		w.write("elgmorphing_moved.vtk");

	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================
	gmds::smoothy::EllipticSmoother2D smoother2D(&m_orig);
	smoother2D.lock(mark_bnd_nodes);
	smoother2D.execute();

	if(debug_mode)
		w.write("elgmorphing_out.vtk");

	return 0;
}