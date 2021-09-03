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
class VariableTest : public ::testing::Test
{
protected:
    VariableTest()
    {
        ;
    }
    virtual ~VariableTest()
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
TEST_F(VariableTest, creation)
{
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx -1;
    const int nj = ny -1;

    m.updateNodeCapacity(nx*ny);
    m.updateRegionCapacity(ni*nj);

    kmds::TCellID nodes[nx*ny];
    int index = 0;

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            nodes[index] = m.newNode(i, j, 0.);
            index++;
        }
    }

    for(int i=0; i<ni; i++) {
        for(int j=0; j<nj; j++) {

            m.newQuad(nodes[i*ny+j],
                    nodes[(i+1)*ny+j],
                    nodes[(i+1)*ny+(j+1)],
                    nodes[i*ny+(j+1)]
            );

        }
    }

    struct ColorFace
    {
        kmds::GrowingView<kmds::TCellID>* cellIDs;
        kmds::Variable<int>* v;


        ColorFace(kmds::GrowingView<kmds::TCellID>* cellIDs_, kmds::Variable<int>* v_)
                : cellIDs(cellIDs_)
                , v(v_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(int i) const {
            if (cellIDs->get(i) % 2 == 0) {
                (*v)[cellIDs->get(i)] = 0;
            } else {
                (*v)[cellIDs->get(i)] = 1;
            }
        }
    };

    kmds::GrowingView<kmds::TCellID> cellIDs("CELLS", m.getNbFaces());
    m.getFaceIDs(&cellIDs);

    kmds::Variable<int>* v = m.createVariable<int>(-1, kmds::KMDS_FACE, "color");
    Kokkos::parallel_for(cellIDs.getNbElems(), ColorFace(&cellIDs, v));

    EXPECT_EQ(0, (*v)[0]);
    EXPECT_EQ(1, (*v)[1]);

    kmds::Variable<int>* v_bis = m.getVariable<int>(kmds::KMDS_FACE, "color");
    EXPECT_EQ(0, (*v_bis)[0]);
    EXPECT_EQ(1, (*v_bis)[1]);
}
/*----------------------------------------------------------------------------*/
