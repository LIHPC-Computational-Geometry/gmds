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
#include <KM/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class WriterTest : public ::testing::Test
{
protected:
    WriterTest()
    {
        ;
    }
    virtual ~WriterTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
		  Kokkos::InitializationSettings kargs;
		  kargs.set_num_threads(3);
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
TEST_F(WriterTest, N2F_quadtri)
{
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;

    m.updateNodeCapacity(nx*ny+1);
    m.updateFaceCapacity(ni*nj+2);

    kmds::TCellID nodeIDs[nx*ny+1];
    int index = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            nodeIDs[index] = m.newNode(i, j, 0.);
            index++;
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

    kmds::VTKWriter<kmds::Mesh> w(m);
    w.write("vtkwritertest_2D", kmds::F);
}
/*----------------------------------------------------------------------------*/
TEST_F(WriterTest, N2R_hexpyrtetprism)
{
    kmds::Mesh m;

    m.updateNodeCapacity(13);
    m.updateRegionCapacity(4);

    kmds::TCellID nodeIDs[13];
    nodeIDs[0] = m.newNode(0., 0., 0.);
    nodeIDs[1] = m.newNode(1., 0., 0.);
    nodeIDs[2] = m.newNode(1., 1., 0.);
    nodeIDs[3] = m.newNode(0., 1., 0.);
    nodeIDs[4] = m.newNode(0., 0., 1.);
    nodeIDs[5] = m.newNode(1., 0., 1.);
    nodeIDs[6] = m.newNode(1., 1., 1.);
    nodeIDs[7] = m.newNode(0., 1., 1.);

    nodeIDs[8] = m.newNode(0.5, 0.5, 2.);
    nodeIDs[9] = m.newNode(1.5, 0.5, 1.5);

    nodeIDs[10] = m.newNode(0., 1., 1.);
    nodeIDs[11] = m.newNode(0., 1., 1.);
    nodeIDs[12] = m.newNode(0., 1., 1.);


    kmds::TCellID hexId = m.newHexahedron(
            nodeIDs[0],
            nodeIDs[1],
            nodeIDs[2],
            nodeIDs[3],
            nodeIDs[4],
            nodeIDs[5],
            nodeIDs[6],
            nodeIDs[7]
    );

    kmds::TCellID pyrId = m.newPyramid(
            nodeIDs[4],
            nodeIDs[5],
            nodeIDs[6],
            nodeIDs[7],
            nodeIDs[8]
    );

    kmds::TCellID tetId = m.newTetrahedron(
            nodeIDs[5],
            nodeIDs[6],
            nodeIDs[8],
            nodeIDs[9]
    );

    kmds::TCellID prismId = m.newPrism3(
            nodeIDs[7],
            nodeIDs[4],
            nodeIDs[8],
            nodeIDs[10],
            nodeIDs[11],
            nodeIDs[12]
    );

    kmds::Variable<int>* varint = m.createVariable<int>(-1, kmds::KMDS_REGION, "varint");
    (*varint)[hexId] = 17;
    (*varint)[pyrId] = 27;
    (*varint)[tetId] = 37;
    (*varint)[prismId] = 47;

    kmds::Variable<double>* vardouble = m.createVariable<double>(-1., kmds::KMDS_REGION, "vardouble");
    (*vardouble)[hexId] = 7.1;
    (*vardouble)[pyrId] = 7.2;
    (*vardouble)[tetId] = 7.3;
    (*vardouble)[prismId] = 7.4;

    kmds::Variable<gmds::math::Vector>* varvec = m.createVariable<gmds::math::Vector>(gmds::math::Vector({0., 0., 0.}), kmds::KMDS_REGION, "varvec");
    (*varvec)[hexId] = gmds::math::Vector ({1., 0., 0.});
    (*varvec)[pyrId] = gmds::math::Vector ({1., 1., 0.});
    (*varvec)[tetId] = gmds::math::Vector ({1., 1., 1.});
    (*varvec)[prismId] = gmds::math::Vector ({2., 0., 0.});

    kmds::VTKWriter<kmds::Mesh> w(m);
    w.write("vtkwritertest_3D", kmds::R);
}
/*----------------------------------------------------------------------------*/