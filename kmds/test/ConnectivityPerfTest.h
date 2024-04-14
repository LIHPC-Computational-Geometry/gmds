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
class ConnectivityPerfTest : public ::testing::Test
{
protected:
    ConnectivityPerfTest()
    {
        ;
    }
    virtual ~ConnectivityPerfTest()
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
TEST_F(ConnectivityPerfTest, BuildN2R_grid_3D) {
    kmds::Mesh m;

    const int nx = 201;
    const int ny = 201;
    const int nz = 201;

    const int ni = nx - 1;
    const int nj = ny - 1;
    const int nk = nz - 1;

    m.updateNodeCapacity(nx * ny * nz);
    m.updateRegionCapacity(ni * nj * nk);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int a = m.addNode();
                m.setNodeLocation(a, i, j, k);
            }
        }
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                m.newHexahedron(i * ny * nz + j * nz + k,
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

    Kokkos::Timer timer;
    timer.reset();

    kmds::Connectivity *c_N2R = m.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2R_variant_0();

    std::cout<<"BuildN2R "<<timer.seconds()<<std::endl;

}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityPerfTest, BuildN2R_grid_3D_1) {
    kmds::Mesh m;

    const int nx = 201;
    const int ny = 201;
    const int nz = 201;

    const int ni = nx - 1;
    const int nj = ny - 1;
    const int nk = nz - 1;

    m.updateNodeCapacity(nx * ny * nz);
    m.updateRegionCapacity(ni * nj * nk);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int a = m.addNode();
                m.setNodeLocation(a, i, j, k);
            }
        }
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {
                m.newHexahedron(i * ny * nz + j * nz + k,
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

    Kokkos::Timer timer;
    timer.reset();

    kmds::Connectivity *c_N2R = m.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2R_variant_1();

    std::cout<<"BuildN2R "<<timer.seconds()<<std::endl;

}
/*----------------------------------------------------------------------------*/