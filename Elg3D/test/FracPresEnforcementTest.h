/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/FracPresEnforcement.h"
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class FracPresEnforcementTest : public ::testing::Test
{
 protected:
        FracPresEnforcementTest()
        {
                ;
        }
        virtual ~FracPresEnforcementTest()
        {
                ;
        }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
        Kokkos::InitArguments kargs;
        kargs.num_threads = 1;
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
TEST_F(FracPresEnforcementTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresEnforcementTest, maintain)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);


    elg3d::FracPres fp_new;

    kmds::GrowingView <kmds::TCellID> cellIDs("CELLIDS", 2);
    cellIDs.push_back(0);
    cellIDs.push_back(1);
    std::vector<std::string> names;
    names.push_back(std::string("iron"));

    elg3d::FracPresEnforcement_maintain_xD(&cellIDs, &fp, &fp_new, names);

    EXPECT_EQ(fp_new.getNbMaterials(), 2);
    EXPECT_GT(fp_new.getFracPres(0, 0), fp_new.getFracPres(1, 0));
    EXPECT_GT(fp_new.getFracPres(0, 1), fp_new.getFracPres(1, 1));

    std::vector<std::string> namesbis;
    namesbis.push_back(std::string("copper"));

    elg3d::FracPres fp_newbis;
    elg3d::FracPresEnforcement_maintain_xD(&cellIDs, &fp, &fp_newbis, namesbis);

    EXPECT_EQ(fp_newbis.getNbMaterials(), 2);
    EXPECT_LT(fp_newbis.getFracPres(0, 0), fp_newbis.getFracPres(1, 0));
    EXPECT_LT(fp_newbis.getFracPres(0, 1), fp_newbis.getFracPres(1, 1));
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresEnforcementTest, fuse)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);


    elg3d::FracPres fp_new;

    kmds::GrowingView <kmds::TCellID> cellIDs("CELLIDS", 2);
//    cellIDs.resize(2);
//    cellIDs.set(0, 0);
//    cellIDs.set(1, 1);
    cellIDs.push_back(0);
    cellIDs.push_back(1);
    std::vector<std::string> names;
    names.push_back(std::string("iron"));
    names.push_back(std::string("copper"));

    elg3d::FracPresEnforcement_fuse_xD(&cellIDs, &fp, &fp_new, names, "fusedMat", false);

    EXPECT_EQ(fp_new.getNbMaterials(), 1);
    EXPECT_DOUBLE_EQ(fp_new.getFracPres(0, 0), 1.);
    EXPECT_DOUBLE_EQ(fp_new.getFracPres(0, 1), 0.8);


    elg3d::FracPres fp_newbis;
    elg3d::FracPresEnforcement_fuse_xD(&cellIDs, &fp, &fp_newbis, names, "fusedMat", true);

    EXPECT_EQ(fp_newbis.getNbMaterials(), 3);
    EXPECT_DOUBLE_EQ(fp_newbis.getFracPres(0, 0), 0.);
    EXPECT_DOUBLE_EQ(fp_newbis.getFracPres(1, 0), 0.);
    EXPECT_DOUBLE_EQ(fp_newbis.getFracPres(0, 1), 0.);
    EXPECT_DOUBLE_EQ(fp_newbis.getFracPres(1, 1), 0.);
    EXPECT_DOUBLE_EQ(fp_newbis.getFracPres(2, 0), 1.);
    EXPECT_DOUBLE_EQ(fp_newbis.getFracPres(2, 1), 0.8);
}
/*----------------------------------------------------------------------------*/