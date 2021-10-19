/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class MaterialAssignmentTest : public ::testing::Test
{
 protected:
        MaterialAssignmentTest()
        {
                ;
        }
        virtual ~MaterialAssignmentTest()
        {
                ;
        }

    static void
    SetUpTestCase()
    {
            // Kokkos::Serial::initialize();
            // Kokkos::Threads::initialize();
            Kokkos::InitArguments kargs;
            kargs.num_threads = 3;
//            int num_threads = 3;
//            int use_numa = 1;
//            int use_core = 1;
//            Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
            Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
            // Kokkos::Serial::finalize();
            // Kokkos::Threads::finalize();
//            Kokkos::OpenMP::finalize();
            Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(MaterialAssignmentTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialAssignmentTest, createMat)
{
        elg3d::MaterialAssignment ma;

        ma.createMaterial("iron");

        EXPECT_EQ(ma.getNbMaterials(), 1);
        EXPECT_EQ(ma.getMaterialID("iron"), 0);
        EXPECT_EQ(ma.getMaterialName(0), "iron");

        ma.createMaterial("copper");
        EXPECT_EQ(ma.getNbMaterials(), 2);
        EXPECT_EQ(ma.getMaterialID("copper"), 1);
        EXPECT_EQ(ma.getMaterialName(1), "copper");
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialAssignmentTest, values)
{
        elg3d::MaterialAssignment ma(2);

        ma.createMaterial("iron");
        ma.createMaterial("copper");

        ma.setMaterial(0, 0);
        ma.setMaterial(0, 1);
        EXPECT_EQ(ma.getMaterial(0),0);
        //EXPECT_EQ(ma.getCells(0).size(), 2);

        ma.changeMaterial(1,0);
        EXPECT_EQ(ma.getMaterial(0),1);
        //EXPECT_EQ(ma.getCells(0).size(), 1);
        //EXPECT_EQ(ma.getCells(1).size(), 1);
}
/*----------------------------------------------------------------------------*/