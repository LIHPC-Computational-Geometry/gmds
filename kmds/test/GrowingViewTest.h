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