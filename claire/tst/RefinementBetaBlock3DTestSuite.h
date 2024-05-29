/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/claire/Blocking3D.h>
#include <gmds/claire/RefinementBetaBlock3D.h>
#include <gmds/claire/Utils.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(RefinementBetaBlock3DTestSuite, RefinementBetaBlock3D_1)
{
	Blocking3D b;
	Node n1 = b.newBlockCorner(0,0, 0);
	Node n2 = b.newBlockCorner(1,0, 0);
	Node n3 = b.newBlockCorner(1,1, 0);
	Node n4=  b.newBlockCorner(0,1, 0);
	Node n5=  b.newBlockCorner(0,0, 1);
	Node n6=  b.newBlockCorner(1,0, 1);
	Node n7=  b.newBlockCorner(1,1, 1);
	Node n8=  b.newBlockCorner(0,1, 1);

	Blocking3D::Block b0 = b.newBlock(n1,n2,n3,n4,n5,n6,n7,n8);

	b0.setNbDiscretizationI(4);
	b0.setNbDiscretizationJ(8);
	b0.setNbDiscretizationK(12);
	b.initializeGridPoints();

	b0 = b.block(0);

	TCellID f_i0_id = math::Utils::CommonFace(&b, n1.id(), n4.id(), n8.id(), n5.id());
	TCellID f_imax_id = math::Utils::CommonFace(&b, n2.id(), n3.id(), n7.id(), n6.id());
	TCellID f_j0_id = math::Utils::CommonFace(&b, n1.id(), n2.id(), n6.id(), n5.id());
	TCellID f_jmax_id = math::Utils::CommonFace(&b, n3.id(), n4.id(), n8.id(), n7.id());
	TCellID f_k0_id = math::Utils::CommonFace(&b, n1.id(), n2.id(), n3.id(), n4.id());
	TCellID f_kmax_id = math::Utils::CommonFace(&b, n5.id(), n6.id(), n7.id(), n8.id());

	double eps(pow(10,-6));

	// FACE I=0
	Blocking3D::BlockFace bf = b.blockFace(f_i0_id);
	RefinementBetaBlock3D algo_ref_i0 = RefinementBetaBlock3D(&b0, &bf, 0.01);
	RefinementBetaBlock3D::STATUS algo_res_i0 = algo_ref_i0.execute();

	ASSERT_EQ(algo_res_i0,RefinementBetaBlock3D::SUCCESS);
	ASSERT_NEAR(b0(1,2,2).X(),0.01, eps);

	// FACE I=MAX
	bf = b.blockFace(f_imax_id);
	RefinementBetaBlock3D algo_ref_imax = RefinementBetaBlock3D(&b0, &bf, 0.01);
	RefinementBetaBlock3D::STATUS algo_res_imax = algo_ref_imax.execute();

	ASSERT_EQ(algo_res_imax,RefinementBetaBlock3D::SUCCESS);
	ASSERT_NEAR(b0(b0.getNbDiscretizationI()-2,3,7).X(), 0.99, eps);

	// FACE J=0
	bf = b.blockFace(f_j0_id);
	RefinementBetaBlock3D algo_ref_j0 = RefinementBetaBlock3D(&b0, &bf, 0.01);
	RefinementBetaBlock3D::STATUS algo_res_j0 = algo_ref_j0.execute();

	ASSERT_EQ(algo_res_j0,RefinementBetaBlock3D::SUCCESS);
	ASSERT_NEAR(b0(0,1,8).Y(),0.01, eps);

	// FACE J=MAX
	bf = b.blockFace(f_jmax_id);
	RefinementBetaBlock3D algo_ref_jmax = RefinementBetaBlock3D(&b0, &bf, 0.01);
	RefinementBetaBlock3D::STATUS algo_res_jmax = algo_ref_jmax.execute();

	ASSERT_EQ(algo_res_jmax,RefinementBetaBlock3D::SUCCESS);
	ASSERT_NEAR(b0(3,b0.getNbDiscretizationJ()-2,10).Y(), 0.99, eps);

	// FACE K=0
	bf = b.blockFace(f_k0_id);
	RefinementBetaBlock3D algo_ref_k0 = RefinementBetaBlock3D(&b0, &bf, 0.01);
	RefinementBetaBlock3D::STATUS algo_res_k0 = algo_ref_k0.execute();

	ASSERT_EQ(algo_res_k0,RefinementBetaBlock3D::SUCCESS);
	ASSERT_NEAR(b0(2,7,1).Z(),0.01, eps);

	// FACE K=MAX
	bf = b.blockFace(f_kmax_id);
	RefinementBetaBlock3D algo_ref_kmax = RefinementBetaBlock3D(&b0, &bf, 0.01);
	RefinementBetaBlock3D::STATUS algo_res_kmax = algo_ref_kmax.execute();

	ASSERT_EQ(algo_res_kmax,RefinementBetaBlock3D::SUCCESS);
	ASSERT_NEAR(b0(3,4,b0.getNbDiscretizationK()-2).Z(), 0.99, eps);


	gmds::IGMeshIOService ioService(&b);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("RefinementBetaBlock3DTestSuite.vtk");
}