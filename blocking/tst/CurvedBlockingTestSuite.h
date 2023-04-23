/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingTestSuite, init)
{
	gmds::blocking::CurvedBlocking bl;
	gmds::blocking::DartHandler d = bl.createHex();
	gmds::blocking::DartHandler d2 = bl.createHex();
	std::cout<<" ---------- before sewing 2 hexes ------------"<<std::endl;
	std::cout<<bl.info()<<std::endl;

	bl.sew(d,d2);
	std::cout<<" ---------- and after now!!  ------------"<<std::endl;
	std::cout<<bl.info()<<std::endl;
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
