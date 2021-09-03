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
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/ExtractGeomModel.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/Tools.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class ExtractGeomTest : public ::testing::Test
{
 protected:
        ExtractGeomTest()
        {
                ;
        }
        virtual ~ExtractGeomTest()
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
TEST_F(ExtractGeomTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(ExtractGeomTest, sample_3x3_2D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3_2D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

    elg3d::Tools_write_2D(&mesh, &fp, &ma, "ExtractGeomTest_sample_3x3_2D");


    kmds::Connectivity* c_E2F = mesh.createConnectivity(kmds::E2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildEandE2F_2D_variant_1();

    gmds::cad::FACManager AGeomModel;
    elg3d::extractGeomModel_extract_2D(&mesh, &ma, c_E2F, &AGeomModel);

//    EXPECT_EQ(4, AGeomModel.getNbPoints());
//    EXPECT_EQ(3, AGeomModel.getNbCurves());
}
/*----------------------------------------------------------------------------*/
TEST_F(ExtractGeomTest, sample_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    elg3d::Tools_write_3D(&mesh, &fp, &ma, "ExtractGeomTest_sample_3x3x3_3D");


    kmds::Connectivity* c_F2R = mesh.createConnectivity(kmds::F2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildFandF2R();

    kmds::Variable <std::uintptr_t> *geomassoc = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr),
                                                                                                   kmds::KMDS_NODE, "geomassoc");
    gmds::cad::FACManager geomModel;

    elg3d::extractGeomModel_extract_3D(&mesh, &ma, c_F2R, &geomModel, geomassoc);

    EXPECT_EQ(3, geomModel.getNbSurfaces());
    EXPECT_EQ(4, geomModel.getNbCurves());
    // this is a multipoint
    EXPECT_EQ(1, geomModel.getNbPoints());
}
/*----------------------------------------------------------------------------*/