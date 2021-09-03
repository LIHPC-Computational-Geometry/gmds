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
class FakeTypesTest : public ::testing::Test
{
protected:
    FakeTypesTest()
    {
        ;
    }
    virtual ~FakeTypesTest()
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
TEST_F(FakeTypesTest, hexpyrtetprism)
{
    kmds::Mesh m;

    m.updateNodeCapacity(13);
    m.updateRegionCapacity(4);

    kmds::TCellID nodes[13];
    nodes[0] = m.newNode(0., 0., 0.);
    nodes[1] = m.newNode(1., 0., 0.);
    nodes[2] = m.newNode(1., 1., 0.);
    nodes[3] = m.newNode(0., 1., 0.);
    nodes[4] = m.newNode(0., 0., 1.);
    nodes[5] = m.newNode(1., 0., 1.);
    nodes[6] = m.newNode(1., 1., 1.);
    nodes[7] = m.newNode(0., 1., 1.);

    nodes[8] = m.newNode(0.5, 0.5, 2.);
    nodes[9] = m.newNode(1.5, 0.5, 1.5);

    nodes[10] = m.newNode(0., 1., 1.);
    nodes[11] = m.newNode(0., 1., 1.);
    nodes[12] = m.newNode(0., 1., 1.);


    kmds::TCellID hexId = m.newHexahedron(
            nodes[0],
            nodes[1],
            nodes[2],
            nodes[3],
            nodes[4],
            nodes[5],
            nodes[6],
            nodes[7]
    );

    kmds::TCellID pyrId = m.newPyramid(
            nodes[4],
            nodes[5],
            nodes[6],
            nodes[7],
            nodes[8]
    );

    kmds::TCellID tetId = m.newTetrahedron(
            nodes[5],
            nodes[6],
            nodes[8],
            nodes[9]
    );

    kmds::TCellID prismId = m.newPrism3(
            nodes[7],
            nodes[4],
            nodes[8],
            nodes[10],
            nodes[11],
            nodes[12]
    );

    kmds::Region r = m.getRegion(hexId);

    std::vector<kmds::FakeFace> fakeFaces = r.getFakeFaces();
    EXPECT_EQ(6, fakeFaces.size());

    for(auto ff: fakeFaces) {
        std::vector<kmds::FakeEdge> fakeEdges = ff.getFakeEdges();
        EXPECT_EQ(4, fakeEdges.size());
    }

    r = m.getRegion(pyrId);
    fakeFaces = r.getFakeFaces();
    EXPECT_EQ(5, fakeFaces.size());

    r = m.getRegion(tetId);
    fakeFaces = r.getFakeFaces();
    EXPECT_EQ(4, fakeFaces.size());

    r = m.getRegion(prismId);
    fakeFaces = r.getFakeFaces();
    EXPECT_EQ(5, fakeFaces.size());
}
/*----------------------------------------------------------------------------*/