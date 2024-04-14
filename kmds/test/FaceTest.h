/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class FaceTest : public ::testing::Test
{
protected:
    FaceTest()
    {
        ;
    }
    virtual ~FaceTest()
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
TEST_F(FaceTest, FaceCreation_quad) {
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;

    m.updateNodeCapacity(nx*ny);
    m.updateFaceCapacity(ni*nj);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {

            m.newQuad(i * ny + j,
                      (i + 1) * ny + j,
                      (i + 1) * ny + (j + 1),
                      i * ny + (j + 1)
            );

        }
    }

    EXPECT_EQ(nx * ny, m.getNbNodes());
    EXPECT_EQ(ni * nj, m.getNbFaces());
    EXPECT_EQ(m.getNbFaces(), m.getNbQuads());

    kmds::Face f = m.getFace(0);
    Kokkos::View<kmds::TCellID *> ids;
    f.nodeIds(ids);

    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(0, ids(0));
    EXPECT_EQ(3, ids(1));
    EXPECT_EQ(4, ids(2));
    EXPECT_EQ(1, ids(3));

    f = m.getFace(1);
    f.nodeIds(ids);

    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(1, ids(0));
    EXPECT_EQ(4, ids(1));
    EXPECT_EQ(5, ids(2));
    EXPECT_EQ(2, ids(3));

    f = m.getFace(2);
    f.nodeIds(ids);

    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(3, ids(0));
    EXPECT_EQ(6, ids(1));
    EXPECT_EQ(7, ids(2));
    EXPECT_EQ(4, ids(3));

    f = m.getFace(3);
    f.nodeIds(ids);

    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(4, ids(0));
    EXPECT_EQ(7, ids(1));
    EXPECT_EQ(8, ids(2));
    EXPECT_EQ(5, ids(3));
}
/*----------------------------------------------------------------------------*/
TEST_F(FaceTest, FaceCreation_quadtriangle) {
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;

    m.updateNodeCapacity(nx*ny+1);
    m.updateFaceCapacity(ni*nj+2);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }
    int a = m.addNode();
    m.setNodeLocation(a, nx, 1., 0.);

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {

            m.newQuad(i * ny + j,
                      (i + 1) * ny + j,
                      (i + 1) * ny + (j + 1),
                      i * ny + (j + 1)
            );

        }
    }
    m.newTriangle((nx-1)*ny,
                  a,
                  (nx-1)*ny+1
    );
    m.newTriangle((nx-1)*ny+1,
                  a,
                  (nx-1)*ny+2
    );

    EXPECT_EQ(nx * ny + 1, m.getNbNodes());
    EXPECT_EQ(ni * nj + 2, m.getNbFaces());
    EXPECT_EQ(ni * nj, m.getNbQuads());
    EXPECT_EQ(2, m.getNbTriangles());

    kmds::Face f = m.getFace(0);
    Kokkos::View<kmds::TCellID *> ids;
    f.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_QUAD, f.computeType());
    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(0, ids(0));
    EXPECT_EQ(3, ids(1));
    EXPECT_EQ(4, ids(2));
    EXPECT_EQ(1, ids(3));

    f = m.getFace(1);
    f.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_QUAD, f.computeType());
    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(1, ids(0));
    EXPECT_EQ(4, ids(1));
    EXPECT_EQ(5, ids(2));
    EXPECT_EQ(2, ids(3));

    f = m.getFace(2);
    f.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_QUAD, f.computeType());
    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(3, ids(0));
    EXPECT_EQ(6, ids(1));
    EXPECT_EQ(7, ids(2));
    EXPECT_EQ(4, ids(3));

    f = m.getFace(3);
    f.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_QUAD, f.computeType());
    EXPECT_EQ(4, f.getNbNodes());
    EXPECT_EQ(4, ids.extent(0));
    EXPECT_EQ(4, ids.size());

    EXPECT_EQ(4, ids(0));
    EXPECT_EQ(7, ids(1));
    EXPECT_EQ(8, ids(2));
    EXPECT_EQ(5, ids(3));

    f = m.getFace(4);
    f.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_TRIANGLE, f.computeType());
    EXPECT_EQ(3, f.getNbNodes());
    EXPECT_EQ(3, ids.extent(0));
    EXPECT_EQ(3, ids.size());

    EXPECT_EQ(6, ids(0));
    EXPECT_EQ(9, ids(1));
    EXPECT_EQ(7, ids(2));

    f = m.getFace(5);
    f.nodeIds(ids);

    EXPECT_EQ(kmds::KMDS_TRIANGLE, f.computeType());
    EXPECT_EQ(3, f.getNbNodes());
    EXPECT_EQ(3, ids.extent(0));
    EXPECT_EQ(3, ids.size());

    EXPECT_EQ(7, ids(0));
    EXPECT_EQ(9, ids(1));
    EXPECT_EQ(8, ids(2));
}
/*----------------------------------------------------------------------------*/
TEST_F(FaceTest, edgeOrientation) {
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;

    m.updateNodeCapacity(nx * ny);
    m.updateFaceCapacity(ni * nj);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {

            m.newQuad(i * ny + j,
                      (i + 1) * ny + j,
                      (i + 1) * ny + (j + 1),
                      i * ny + (j + 1)
            );

        }
    }


    kmds::Face f = m.getFace(0);
    Kokkos::View<kmds::TCellID *> ids;
    f.nodeIds(ids);

    Kokkos::View<kmds::TCellID *> e("EDGE_NODES", 2);

    e(0) = ids(0);
    e(1) = ids(1);
    EXPECT_TRUE(f.isEdgeOrientedOutward(e));

    e(0) = ids(3);
    e(1) = ids(0);
    EXPECT_TRUE(f.isEdgeOrientedOutward(e));

    e(0) = ids(1);
    e(1) = ids(0);
    EXPECT_FALSE(f.isEdgeOrientedOutward(e));

    e(0) = ids(0);
    e(1) = ids(3);
    EXPECT_FALSE(f.isEdgeOrientedOutward(e));
}
/*----------------------------------------------------------------------------*/