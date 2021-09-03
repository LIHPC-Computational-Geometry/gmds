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