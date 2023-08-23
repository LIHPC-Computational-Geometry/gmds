/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include<KM/Utils/GrowingView.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class GrowingViewTest : public ::testing::Test
{
protected:
    GrowingViewTest()
    {
        ;
    }
    virtual ~GrowingViewTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
		  Kokkos::InitializationSettings kargs;
		  kargs.set_num_threads(3);
//        Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
        Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
        // Kokkos::Serial::finalize();
        // Kokkos::Threads::finalize();
//        Kokkos::OpenMP::finalize();
        Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(GrowingViewTest, creation)
{

    kmds::GrowingView<int> gv("GV_NAME", 11);
    EXPECT_EQ(0, gv.getNbElems());

    gv.push_back(1);
    gv.push_back(3);
    EXPECT_EQ(2, gv.getNbElems());
    EXPECT_EQ(1, gv.get(0));

    gv.resize(17);
    EXPECT_EQ(17, gv.getNbElems());

    gv.push_back(4);
    EXPECT_EQ(18, gv.getNbElems());

    kmds::TInt32 first = gv.addElems(20);
    EXPECT_EQ(18, first);
    EXPECT_EQ(38, gv.getNbElems());
    EXPECT_EQ(4, gv.get(17));

    gv.set(17, 5);
    EXPECT_EQ(5, gv.get(17));
}
/*----------------------------------------------------------------------------*/
TEST_F(GrowingViewTest, management)
{

    kmds::GrowingView<int> gv("GV_NAME", 11);
    EXPECT_GE(11, gv.capacity());

    gv.reserve(10);
    EXPECT_GE(gv.capacity(), 10);

    gv.reserve(320);
    EXPECT_GE(gv.capacity(), 320);

    gv.setTop(217);
    kmds::TInt32 top = gv.push_back(47);
    EXPECT_EQ(217, top);

    gv.clear();
    EXPECT_EQ(0, gv.getNbElems());

    kmds::TInt32 first = gv.push_back(7);
    EXPECT_EQ(1, gv.getNbElems());
    EXPECT_EQ(0, first);
}
/*----------------------------------------------------------------------------*/