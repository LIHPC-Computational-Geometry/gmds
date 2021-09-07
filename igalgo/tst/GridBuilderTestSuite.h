/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/GridBuilder.h>
/*----------------------------------------------------------------------------*/
TEST(GridBuildOpClass, test2D)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::F|gmds::F2N));

    gmds::GridBuilder gb(&m,2);

    ASSERT_TRUE(gb.isValid());

    gb.execute(3,1.0, 4, 1.0);

    ASSERT_EQ(m.getNbNodes(),12);
    ASSERT_EQ(m.getNbFaces(),6);

}
/*----------------------------------------------------------------------------*/
TEST(GridBuildOpClass, test3D)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

    gmds::GridBuilder gb(&m,3);

    ASSERT_TRUE(gb.isValid());

    gb.execute(3,1.0, 4, 1.0, 3, 2.0);

    ASSERT_EQ(m.getNbNodes(),36);
    ASSERT_EQ(m.getNbRegions(),12);

}
/*----------------------------------------------------------------------------*/
