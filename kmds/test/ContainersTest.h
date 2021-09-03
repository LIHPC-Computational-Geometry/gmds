/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
//#include <KMDS/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class ContainersTest : public ::testing::Test
{
 protected:
        ContainersTest()
        {
                ;
        }
        virtual ~ContainersTest()
        {
                ;
        }
};
/*----------------------------------------------------------------------------*/
TEST_F(ContainersTest, init)
{
//        int num_threads = 4;
//        int use_numa = 1;
//        int use_core = 1;
//        Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
//        // Mesh m;
//        EXPECT_EQ(true, true);
//
//        Kokkos::OpenMP::finalize();
}
/*----------------------------------------------------------------------------*/
