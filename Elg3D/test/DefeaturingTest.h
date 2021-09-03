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
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <KM/Utils/InitTools.h>

#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class DefeaturingTest : public ::testing::Test
{
 protected:
        DefeaturingTest()
        {
                ;
        }
        virtual ~DefeaturingTest()
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
TEST_F(DefeaturingTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(DefeaturingTest, defeature_3x3x3_3D)
{
    kmds::Mesh mesh;

    kmds::TCoord minXYZ[3];
    minXYZ[0] = 0.;
    minXYZ[1] = 0.;
    minXYZ[2] = 0.;
    kmds::TCoord maxXYZ[3];
    maxXYZ[0] = 1.;
    maxXYZ[1] = 1.;
    maxXYZ[2] = 1.;
    kmds::InitTools_createGrid_3D(&mesh, minXYZ, maxXYZ, 3, 3, 3);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    ma.createMaterial("matA");
    ma.createMaterial("matB");

    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", mesh.getNbRegions());
    mesh.getRegionIDs_dummy(&cellIDs);

    const kmds::TCellID nbCells = cellIDs.getNbElems();

    for(auto i=0; i<nbCells; i++) {
        ma.setMaterial(0, cellIDs.get(i));
    }
    ma.setMaterial(1, 13);

    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();
    kmds::Connectivity *c_R2R_byN = mesh.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(c_N2R);

    EXPECT_EQ(ma.getMaterial(13), 1);

    elg3d::assignCells_defeaturing_3D(&mesh, c_R2R_byN, &ma);

    EXPECT_EQ(ma.getMaterial(13), 0);
}
/*----------------------------------------------------------------------------*/
