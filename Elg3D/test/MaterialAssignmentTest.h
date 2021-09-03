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