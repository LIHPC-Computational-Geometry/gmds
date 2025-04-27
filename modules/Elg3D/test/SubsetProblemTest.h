/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/SubsetProblem.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class SubsetProblemTest : public ::testing::Test
{
 protected:
        SubsetProblemTest()
        {
                ;
        }
        virtual ~SubsetProblemTest()
        {
                ;
        }

    static void
    SetUpTestCase()
    {
    }

    static void
    TearDownTestCase()
    {
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(SubsetProblemTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
//TEST_F(SubsetProblemTest, majorityCriteria_3x3_2D)
//{
//        kmds::Mesh mesh;
//        elg3d::FracPres fp;
//        elg3d::initData_3x3_2D(&mesh, &fp);
//
//        elg3d::MaterialAssignment ma(mesh.getFaceCapacity());
//
//        elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);
//
//        EXPECT_EQ(ma.getMaterial(0), 0);
//        EXPECT_EQ(ma.getMaterial(1), 0);
//        EXPECT_EQ(ma.getMaterial(3), 0);
//        EXPECT_EQ(ma.getMaterial(4), 0);
//        EXPECT_EQ(ma.getMaterial(2), 1);
//        EXPECT_EQ(ma.getMaterial(5), 1);
//        EXPECT_EQ(ma.getMaterial(6), 2);
//        EXPECT_EQ(ma.getMaterial(7), 2);
//        EXPECT_EQ(ma.getMaterial(8), 2);
//}
/*----------------------------------------------------------------------------*/
//TEST_F(SubsetProblemTest, majorityCriteria_2x2_2D_nonmanifold)
//{
//    kmds::Mesh mesh;
//    elg3d::FracPres fp;
//    elg3d::initData_2x2_2D_nonmanifold(&mesh, &fp);
//
//    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());
//
//    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);
//
//    EXPECT_EQ(ma.getMaterial(0), 0);
//    EXPECT_EQ(ma.getMaterial(1), 1);
//    EXPECT_EQ(ma.getMaterial(2), 1);
//    EXPECT_EQ(ma.getMaterial(3), 0);
//}
/*----------------------------------------------------------------------------*/
TEST_F(SubsetProblemTest, 3x3x3_3D_0)
{
    kmds::Mesh mesh_ref;
    kmds::Mesh mesh_new;
    elg3d::FracPres fp_ref;
    elg3d::FracPres fp_new;
    elg3d::initData_3x3x3_3D(&mesh_ref, &fp_ref);

    kmds::GrowingView <kmds::TCellID> cellIDs("CELLIDS", 2);
    cellIDs.push_back(0);
    cellIDs.push_back(1);

    elg3d::SubsetProblem_extract_3D(&cellIDs,
                                    &mesh_ref,
                                    &fp_ref,
                                    &mesh_new,
                                    &fp_new,
                                    true
    );

    EXPECT_EQ(2, mesh_new.getNbRegions());
    EXPECT_EQ(12, mesh_new.getNbNodes());

    EXPECT_EQ(3, fp_new.getNbMaterials());
    EXPECT_EQ(fp_ref.getFracPres(0, 0), fp_new.getFracPres(0, 0));
    EXPECT_EQ(fp_ref.getFracPres(1, 0), fp_new.getFracPres(1, 0));
    EXPECT_EQ(fp_ref.getFracPres(2, 0), fp_new.getFracPres(2, 0));
    EXPECT_EQ(fp_ref.getFracPres(0, 1), fp_new.getFracPres(0, 1));
    EXPECT_EQ(fp_ref.getFracPres(1, 1), fp_new.getFracPres(1, 1));
    EXPECT_EQ(fp_ref.getFracPres(2, 1), fp_new.getFracPres(2, 1));
}
/*----------------------------------------------------------------------------*/
TEST_F(SubsetProblemTest, 3x3x3_3D_1)
{
    kmds::Mesh mesh_ref;
    kmds::Mesh mesh_new;
    elg3d::FracPres fp_ref;
    elg3d::FracPres fp_new;
    elg3d::initData_3x3x3_3D(&mesh_ref, &fp_ref);

    kmds::GrowingView <kmds::TCellID> cellIDs("CELLIDS", 2);
    cellIDs.push_back(0);
    cellIDs.push_back(1);

    elg3d::SubsetProblem_extract_3D(&cellIDs,
                                    &mesh_ref,
                                    &fp_ref,
                                    &mesh_new,
                                    &fp_new,
                                    false
    );

    EXPECT_EQ(2, mesh_new.getNbRegions());
    EXPECT_EQ(12, mesh_new.getNbNodes());

    EXPECT_EQ(2, fp_new.getNbMaterials());
    EXPECT_EQ(fp_ref.getFracPres(0, 0), fp_new.getFracPres(0, 0));
    EXPECT_EQ(fp_ref.getFracPres(1, 0), fp_new.getFracPres(1, 0));
    EXPECT_EQ(fp_ref.getFracPres(0, 1), fp_new.getFracPres(0, 1));
    EXPECT_EQ(fp_ref.getFracPres(1, 1), fp_new.getFracPres(1, 1));
}
/*----------------------------------------------------------------------------*/
TEST_F(SubsetProblemTest, 3x3x3_3D_2)
{
    kmds::Mesh mesh_ref;
    kmds::Mesh mesh_new;
    elg3d::FracPres fp_ref;
    elg3d::FracPres fp_new;
    elg3d::initData_3x3x3_3D(&mesh_ref, &fp_ref);

    kmds::GrowingView <kmds::TCellID> cellIDs("CELLIDS", 3);
    cellIDs.push_back(13);
    cellIDs.push_back(26);
    cellIDs.push_back(12);

    elg3d::SubsetProblem_extract_3D(&cellIDs,
                                    &mesh_ref,
                                    &fp_ref,
                                    &mesh_new,
                                    &fp_new,
                                    false
    );

    EXPECT_EQ(3, mesh_new.getNbRegions());
    EXPECT_EQ(19, mesh_new.getNbNodes());

    EXPECT_EQ(3, fp_new.getNbMaterials());
    EXPECT_EQ(fp_ref.getFracPres(0, cellIDs.get(0)), fp_new.getFracPres(0, 0));
    EXPECT_EQ(fp_ref.getFracPres(1, cellIDs.get(0)), fp_new.getFracPres(1, 0));
    EXPECT_EQ(fp_ref.getFracPres(2, cellIDs.get(0)), fp_new.getFracPres(2, 0));
    EXPECT_EQ(fp_ref.getFracPres(0, cellIDs.get(1)), fp_new.getFracPres(0, 1));
    EXPECT_EQ(fp_ref.getFracPres(1, cellIDs.get(1)), fp_new.getFracPres(1, 1));
    EXPECT_EQ(fp_ref.getFracPres(2, cellIDs.get(1)), fp_new.getFracPres(2, 1));
    EXPECT_EQ(fp_ref.getFracPres(0, cellIDs.get(2)), fp_new.getFracPres(0, 2));
    EXPECT_EQ(fp_ref.getFracPres(1, cellIDs.get(2)), fp_new.getFracPres(1, 2));
    EXPECT_EQ(fp_ref.getFracPres(2, cellIDs.get(2)), fp_new.getFracPres(2, 2));
}
/*----------------------------------------------------------------------------*/
TEST_F(SubsetProblemTest, big_3D)
{
    kmds::Mesh mesh_ref;
    kmds::Mesh mesh_new;
    elg3d::FracPres fp_ref;
    elg3d::FracPres fp_new;

    kmds::TCoord xyzmin[3] = {0., 0., 0.};
    kmds::TCoord xyzmax[3] = {100., 20., 20.};

    kmds::InitTools_createGrid_3D(&mesh_ref, xyzmin, xyzmax, 100, 20, 20);

    kmds::GrowingView <kmds::TCellID> cellIDs("CELLIDS", mesh_ref.getNbRegions());

    for(int i=10; i<67; i++) {
        for(int j=5; j<17; j++) {
            for(int k=9; k<14; k++) {
                const kmds::TCellID id = i*20*20 + j*20 + k;
                cellIDs.push_back(id);
            }
        }
    }

    elg3d::SubsetProblem_extract_3D(&cellIDs,
                                    &mesh_ref,
                                    &fp_ref,
                                    &mesh_new,
                                    &fp_new,
                                    false
    );

    EXPECT_EQ(3420, mesh_new.getNbRegions());
    EXPECT_EQ(4524, mesh_new.getNbNodes());

    EXPECT_EQ(0, fp_new.getNbMaterials());
}
/*----------------------------------------------------------------------------*/
TEST_F(SubsetProblemTest, assignCells_IxJxI_3D_nonmanifold)
{
//    const int nb_i = 30;
//    const int nb_j = 10;
//
//    kmds::Mesh mesh;
//    elg3d::FracPres fp;
//    elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, nb_i, nb_j);
//
//    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());
//
//    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);
//
////    EXPECT_EQ(ma.getMaterial(0), 0);
////    EXPECT_EQ(ma.getMaterial(1), 1);
////    EXPECT_EQ(ma.getMaterial(2), 1);
////    EXPECT_EQ(ma.getMaterial(3), 0);
//
//
//    elg3d::Tools_write_3D(&mesh, &fp, &ma, "SubsetProblemTest_assignCells_IxJxI_3D_majority");



}
/*----------------------------------------------------------------------------*/