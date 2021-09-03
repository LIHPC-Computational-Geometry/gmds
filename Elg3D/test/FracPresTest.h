/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class FracPresTest : public ::testing::Test
{
 protected:
        FracPresTest()
        {
                ;
        }
        virtual ~FracPresTest()
        {
                ;
        }
};
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, createMat)
{
        elg3d::FracPres fp;

        fp.createMaterial("iron");

        EXPECT_EQ(fp.getNbMaterials(), 1);
        EXPECT_EQ(fp.getMaterialID("iron"), 0);
        EXPECT_EQ(fp.getMaterialName(0), "iron");

        fp.createMaterial("copper");
        EXPECT_EQ(fp.getNbMaterials(), 2);
        EXPECT_EQ(fp.getMaterialID("copper"), 1);
        EXPECT_EQ(fp.getMaterialName(1), "copper");

        fp.clear();
        EXPECT_EQ(fp.getNbMaterials(), 0);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, values)
{
        elg3d::FracPres fp;

        fp.createMaterial("iron");
        fp.createMaterial("copper");

        fp.setFracPres(0, 0, 0.8);
        fp.setFracPres(0, 1, 0.7);
        EXPECT_EQ(fp.getFracPres(0,1), 0.7);
        EXPECT_EQ(fp.getMatFracPres(0).size(), 2);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, maxValues)
{
        elg3d::FracPres fp;

        fp.createMaterial("iron");
        fp.createMaterial("copper");

        fp.setFracPres(0, 0, 0.2);
        fp.setFracPres(1, 0, 0.8);
        fp.setFracPres(0, 1, 0.7);
        fp.setFracPres(1, 1, 0.1);
        EXPECT_EQ(fp.getMaxMatFracPresIndex(0), 1);
        EXPECT_EQ(fp.getMaxMatFracPresIndex(1), 0);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, getter)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);

    std::map<kmds::TCellID , double> fpmat0 = fp.getMatFracPres(0);
    std::map<kmds::TCellID , double> fpmat1 = fp.getMatFracPres(1);

    EXPECT_EQ(fpmat0.size(), 2);
    EXPECT_EQ(fpmat1.size(), 2);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, setter)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);

    std::map<kmds::TCellID , double> fpmat0 = fp.getMatFracPres(0);
    std::map<kmds::TCellID , double> fpmat1 = fp.getMatFracPres(1);

    elg3d::FracPres fp_new;
    fp_new.setMaterialList(fp.getMaterialList());
    fp_new.setMatFracPres(0, fp.getMatFracPres(0));
    fp_new.setMatFracPres(1, fp.getMatFracPres(1));

    EXPECT_EQ(fp_new.getFracPres(0,1), fp.getFracPres(0,1));
    EXPECT_EQ(fp_new.getMatFracPres(0).size(), fp.getMatFracPres(0).size());
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, normalization_0)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);
    EXPECT_FALSE(fp.checkNormalized());

    fp.normalize(1);
    EXPECT_TRUE(fp.checkNormalized());
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, normalization_1)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);
    EXPECT_FALSE(fp.checkNormalized());
    EXPECT_TRUE(fp.normalize());
    EXPECT_TRUE(fp.checkNormalized());
    EXPECT_FALSE(fp.normalize());
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, traceAmounts)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 10e-10);
    fp.setFracPres(1, 1, 1.);
    EXPECT_TRUE(fp.isMixedCell(1));

    fp.removeTraceAmounts();
    EXPECT_FALSE(fp.isMixedCell(1));
    EXPECT_EQ(0, fp.getFracPres(0,1));
    EXPECT_EQ(1, fp.getFracPres(1,1));
}
/*----------------------------------------------------------------------------*/