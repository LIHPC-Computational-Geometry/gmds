/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, init)
{
	gmds::blocking::CurvedBlocking bl;
	gmds::math::Point p000(0,0,0);
	gmds::math::Point p010(0,1,0);
	gmds::math::Point p110(1,1,0);
	gmds::math::Point p100(1,0,0);

	gmds::math::Point p001(0,0,1);
	gmds::math::Point p011(0,1,1);
	gmds::math::Point p111(1,1,1);
	gmds::math::Point p101(1,0,1);

	gmds::math::Point p002(0,0,2);
	gmds::math::Point p012(0,1,2);
	gmds::math::Point p112(1,1,2);
	gmds::math::Point p102(1,0,2);
	auto b1 = bl.createBlock(p000, p010, p110, p100, p001, p011, p111, p101);
	auto b2 = bl.createBlock(p001, p011, p111, p101, p002, p012, p112, p102);
	std::cout<<" ---------- before sewing 2 hexes ------------"<<std::endl;
	std::cout<<bl.info()<<std::endl;

	bl.sew<3>(b1->dart(),b2->dart());
	std::cout<<" ---------- and after now!!  ------------"<<std::endl;
	std::cout<<bl.info()<<std::endl;
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, single_block)
{
	gmds::blocking::CurvedBlocking bl;
	gmds::math::Point p000(0, 0, 0);
	gmds::math::Point p010(0, 1, 0);
	gmds::math::Point p110(1, 1, 0);
	gmds::math::Point p100(1, 0, 0);

	gmds::math::Point p001(0, 0, 1);
	gmds::math::Point p011(0, 1, 1);
	gmds::math::Point p111(1, 1, 1);
	gmds::math::Point p101(1, 0, 1);

	gmds::math::Point p002(0, 0, 2);
	gmds::math::Point p012(0, 1, 2);
	gmds::math::Point p112(1, 1, 2);
	gmds::math::Point p102(1, 0, 2);
	auto b = bl.createBlock(p000, p010, p110, p100, p001, p011, p111, p101);
	ASSERT_EQ(b->info().geom_dim, 4);
	ASSERT_EQ(b->info().geom_id, -1);
	auto fs = bl.get_faces_of_block(b);
	ASSERT_EQ(fs.size(),6);
	auto block_center = bl.get_center_of_block(b);
	for(auto f:fs){
		gmds::math::Point face_center = bl.get_center_of_face(f);
		ASSERT_NEAR(block_center.distance(face_center),0.5, 1e-8);
	}
}