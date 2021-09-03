/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
        Kokkos::InitArguments kargs;
        kargs.num_threads = 3;
//        int num_threads = 4;
//        int use_numa = 1;
//        int use_core = 1;
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