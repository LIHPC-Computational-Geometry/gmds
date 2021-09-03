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
#include "ELG3D/ALGOCMPNT/OptimizationSmooth.h"
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class OptimizationSmoothTest : public ::testing::Test
{
protected:
    OptimizationSmoothTest()
    {
        ;
    }
    virtual ~OptimizationSmoothTest()
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
TEST_F(OptimizationSmoothTest, 3x3_0_2D) {
    kmds::Mesh mesh;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;

    double minXYZ[3] = {0.,0.,0.};
    double maxXYZ[3] = {2.,2.,0.};

    kmds::InitTools_createGrid_2D(&mesh, minXYZ, maxXYZ, ni, nj);

    mesh.setNodeLocation(4, 0.5, 0.5, 0.);

    kmds::Variable<std::uintptr_t>* varGeomAssoc = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_2D(&mesh, facGeomManager, varGeomAssoc);


    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();

    ch.buildN2N(kmds::F);
    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);

    kmds::GrowingView<kmds::TCellID> fixedNodes("FIXED_NODES", mesh.getNbNodes());

    elg3d::optimizationSmooth_2D(3,
                                 &mesh,
                                 c_N2F,
                                 c_N2N,
                                 varGeomAssoc,
                                 &fixedNodes);

    kmds::TCoord xyz[3];
    mesh.getNodeLocation(4, xyz[0], xyz[1], xyz[2]);
    std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;

    mesh.getNodeLocation(0, xyz[0], xyz[1], xyz[2]);
    std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("optimizationSmooth_test", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(OptimizationSmoothTest, 3x3x3_0_3D) {
    kmds::Mesh mesh;

    const int nx = 3;
    const int ny = 3;
    const int nz = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;
    const int nk = nz - 1;

    mesh.updateNodeCapacity(nx * ny * nz);
    mesh.updateRegionCapacity(ni * nj * nk);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int a = mesh.addNode();
                mesh.setNodeLocation(a, i, j, k);
            }
        }
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                mesh.newHexahedron(i * ny * nz + j * nz + k,
                                (i + 1) * ny * nz + j * nz + k,
                                (i + 1) * ny * nz + (j + 1) * nz + k,
                                i * ny * nz + (j + 1) * nz + k,
                                i * ny * nz + j * nz + k + 1,
                                (i + 1) * ny * nz + j * nz + k + 1,
                                (i + 1) * ny * nz + (j + 1) * nz + k + 1,
                                i * ny * nz + (j + 1) * nz + k + 1

                );
            }
        }
    }

    mesh.setNodeLocation(13, 0.5, 0.5, 0.5);

    kmds::Variable<std::uintptr_t>* varGeomAssoc = mesh.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&mesh, facGeomManager, varGeomAssoc);


    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    ch.buildN2N(kmds::R);
    kmds::Connectivity* c_N2N = mesh.getConnectivity(kmds::N2N);

    kmds::GrowingView<kmds::TCellID> fixedNodes("FIXED_NODES", mesh.getNbNodes());

    elg3d::optimizationSmooth_3D(3,
                                 &mesh,
                                 c_N2R,
                                 c_N2N,
                                 varGeomAssoc,
                                 &fixedNodes);

    kmds::TCoord xyz[3];
    mesh.getNodeLocation(13, xyz[0], xyz[1], xyz[2]);
    std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;

    mesh.getNodeLocation(0, xyz[0], xyz[1], xyz[2]);
    std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;

    kmds::VTKWriter<kmds::Mesh> w(mesh);
    w.write("optimizationSmooth_test", kmds::R);
}
/*----------------------------------------------------------------------------*/