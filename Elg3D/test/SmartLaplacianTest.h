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
#include "ELG3D/ALGOCMPNT/SmartLaplacian.h"
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class SmartLaplacianTest : public ::testing::Test
{
protected:
    SmartLaplacianTest()
    {
        ;
    }
    virtual ~SmartLaplacianTest()
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
TEST_F(SmartLaplacianTest, BuildFandF2R_grid) {
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


    elg3d::smartLaplacian_3D(10,
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
}
/*----------------------------------------------------------------------------*/