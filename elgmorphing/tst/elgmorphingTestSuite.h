/*----------------------------------------------------------------------------*/
#include <gmds/elgmorphing/ElgMorphing.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/MdlReader.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Point.h>
#include <gmds/smoothy/EllipticSmoother2D.h>

#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, bbtransform2d)
{
	gmds::math::Point pt(1,1,1);
	gmds::TCoord minXYZ_orig[3] = {-1,-2,-3};
	gmds::TCoord maxXYZ_orig[3] = {2,3,4};
	gmds::TCoord minXYZ_dest[3] = {10,12,13};
	gmds::TCoord maxXYZ_dest[3] = {11,15,15};

	gmds::elgmorphing::ElgMorphing elg;
	gmds::math::Point ptnew = elg.bbtransform2d(pt, minXYZ_orig, maxXYZ_orig, minXYZ_dest, maxXYZ_dest);

	ASSERT_DOUBLE_EQ(ptnew.X(), 10.+1.*(2./3.));
	ASSERT_DOUBLE_EQ(ptnew.Y(), 12.+3.*(3./5.));
//	ASSERT_DOUBLE_EQ(ptnew.Z(), 13.+2.*(4./7.));
}
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, boundingbox2d)
{
	gmds::Mesh m_dest(gmds::MeshModel(gmds::DIM2|gmds::E|gmds::N|gmds::E2N));

	std::string dir(TEST_SAMPLES_DIR);

	gmds::MdlReader reader_dest(m_dest);
	reader_dest.read(dir + "/test4Mdl.mdl");
	reader_dest.createVariablesFromGroups();

	gmds::elgmorphing::ElgMorphing elg;
	gmds::TCoord minXYZ_dest[3];
	gmds::TCoord maxXYZ_dest[3];
	elg.computeBoundingBox(&m_dest, std::string("ALU"), minXYZ_dest, maxXYZ_dest);

	gmds::math::Point ptmin(minXYZ_dest[0],minXYZ_dest[1],minXYZ_dest[2]);
	gmds::math::Point ptmax(maxXYZ_dest[0],maxXYZ_dest[1],maxXYZ_dest[2]);

	ASSERT_EQ(gmds::math::Point(0.,-1.,0.), ptmin);
	ASSERT_EQ(gmds::math::Point(4.,2.,0.), ptmax);
}
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, DISABLED_morphing2d)
{
	gmds::Mesh m_dest(gmds::MeshModel(gmds::DIM2|gmds::E|gmds::N|gmds::E2N));

	std::string dir(TEST_SAMPLES_DIR);

	gmds::MdlReader reader_dest(m_dest);
	reader_dest.read(dir + "/test4Mdl.mdl");
	reader_dest.createVariablesFromGroups();

	gmds::elgmorphing::ElgMorphing elg;
	gmds::TCoord minXYZ_dest[3];
	gmds::TCoord maxXYZ_dest[3];
	elg.computeBoundingBox(&m_dest, std::string("ALU"), minXYZ_dest, maxXYZ_dest);

	//==================================================================
	// ORIGIN MESH
	//==================================================================
	auto debug_mode=true;
	gmds::Mesh m(gmds::MeshModel(gmds::DIM2 | gmds::F | gmds::E | gmds::N | gmds::F2N | gmds::F2E | gmds::E2F | gmds::E2N | gmds::N2E | gmds::N2F));

	std::string vtk_file = dir+"/circle_in_square.vtk";
	gmds::IGMeshIOService ioService(&m);

	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();
	for(auto f_id:m.faces()) {
		gmds::Face f=m.get<gmds::Face>(f_id);
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
	// we get all the nodes that are on the mesh boundary (including internal ones)
	auto nb_locked = 0;
	auto mark_bnd_nodes = m.newMark<gmds::Node>();
	gmds::math::Point origin({0,0,0});
	gmds::math::Vector perturb({5,5,0});
	gmds::Variable<int>* var_bnd = m.newVariable<int,gmds::GMDS_NODE>("bnd");
	gmds::CellGroup<gmds::Node>* grp = m.newGroup<gmds::Node>("circle");
	for (auto n_id : m.nodes()) {
		gmds::math::Point p = m.get<gmds::Node>(n_id).point();

		if (p.Y()==-10 ||p.Y()==10 ||p.X()==10 ||p.X()==-10 ||
		    fabs(p.distance(origin)-2)<0.1) {
			nb_locked += 1;
			m.mark<gmds::Node>(n_id, mark_bnd_nodes);
			var_bnd->set(n_id, 1);

			// fill a group with the node on the circle
			if (fabs(p.distance(origin) - 2) < 0.1) {
				grp->add(n_id);
			}
		}
		else{
			var_bnd->set(n_id,0);
		}
	}

	gmds::TCoord minXYZ_orig[3];
	gmds::TCoord maxXYZ_orig[3];
	elg.computeBoundingBox(&m, std::string("circle"), minXYZ_orig, maxXYZ_orig);
	for(auto nid: grp->cells()) {
		gmds::math::Point pt = elg.bbtransform2d(m.get<gmds::Node>(nid).point(), minXYZ_orig, maxXYZ_orig, minXYZ_dest, maxXYZ_dest);
		m.get<gmds::Node>(nid).setPoint(pt);
	}

	gmds::VTKWriter w(&ioService);
	w.setCellOptions(gmds::N|gmds::F);
	w.setDataOptions(gmds::N|gmds::F);
	if(debug_mode)
		w.write("elgmorphing_move_circle_marked.vtk");
	//==================================================================
	// PERFORM THE MESH SMOOTHING NOW
	//==================================================================
	gmds::smoothy::EllipticSmoother2D smoother2D(&m);
	smoother2D.lock(mark_bnd_nodes);
	smoother2D.execute();

	if(debug_mode)
		w.write("elgmorphing_move_circle_out.vtk");

	for(auto f_id:m.faces()){
		gmds::Face f = m.get<gmds::Face>(f_id);
		std::vector<gmds::Node> f_nodes = f.get<gmds::Node>();
		// We can do it because we know the mesh is full quad
//		quality::QuadQuality qe = quality::QuadQuality::build(f_nodes[0].point(),
//		                                                      f_nodes[1].point(),
//		                                                      f_nodes[2].point(),
//		                                                      f_nodes[3].point());
//		ASSERT_NEAR(qe.scaledJacobian(),1,0.9);
	}
}
/*----------------------------------------------------------------------------*/