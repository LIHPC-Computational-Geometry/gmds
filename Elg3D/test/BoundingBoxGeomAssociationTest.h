/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
#include <gmds/cad/GeomManager.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class BoundingBoxGeomAssociationTest : public ::testing::Test
{
 protected:
    BoundingBoxGeomAssociationTest()
        {
                ;
        }
        virtual ~BoundingBoxGeomAssociationTest()
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
//                int num_threads = 4;
//                int use_numa = 1;
//                int use_core = 1;
//                Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
                Kokkos::initialize(kargs);
        }

        static void
        TearDownTestCase()
        {
                // Kokkos::Serial::finalize();
                // Kokkos::Threads::finalize();
//                Kokkos::OpenMP::finalize();
                Kokkos::finalize();
        }
};
/*----------------------------------------------------------------------------*/
TEST_F(BoundingBoxGeomAssociationTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(BoundingBoxGeomAssociationTest, createData_3x3_2D)
{
    const int nb_x = 4;
    const int nb_y = 4;
    const int nb_i = nb_x - 1;
    const int nb_j = nb_y - 1;

    kmds::Mesh m;

    m.updateNodeCapacity(nb_x*nb_y);
    m.updateFaceCapacity(nb_i*nb_j);

    for(int i=0; i<nb_x; i++) {
        for(int j=0; j<nb_y; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }

    for(int i=0; i<nb_i; i++) {
        for(int j=0; j<nb_j; j++) {

            m.newQuad(i*nb_y+j,
                      (i+1)*nb_y+j,
                      (i+1)*nb_y+(j+1),
                      i*nb_y+(j+1)
            );

        }
    }


    kmds::Variable<std::uintptr_t>* v = m.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_2D(&m, facGeomManager, v);

    kmds::TSize nbNodes = m.getNbNodes();
    kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", nbNodes);
    m.getNodeIDs_dummy(&nodeIDs);

    int nbNonAssociatedNodes = 0;
    int nbAssociatedNodesCurveBottom = 0;
    int nbAssociatedNodesCorner00 = 0;
    for(int i=0; i<nbNodes; i++) {
        kmds::TCellID nid = nodeIDs.get(i);

        if((*v)[nid] == reinterpret_cast<std::uintptr_t>(nullptr)) {
            nbNonAssociatedNodes++;
        } else {

            if ((reinterpret_cast<gmds::cad::GeomEntity *>((*v)[nid])->name() == std::string("curve_bottom"))) {
                nbAssociatedNodesCurveBottom++;
            }

            if ((reinterpret_cast<gmds::cad::GeomEntity *>((*v)[nid])->name() == std::string("corner_0"))) {
                nbAssociatedNodesCorner00++;
            }
        }
    }

    EXPECT_EQ(nbNonAssociatedNodes, (nb_x-2)*(nb_y-2));
    EXPECT_EQ(nbAssociatedNodesCurveBottom, (nb_x-2));
    EXPECT_EQ(nbAssociatedNodesCorner00, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(BoundingBoxGeomAssociationTest, createData_3x3x3_3D)
{
    const int nb_x = 4;
    const int nb_y = 4;
    const int nb_z = 4;
    const int nb_i = nb_x - 1;
    const int nb_j = nb_y - 1;
    const int nb_k = nb_z - 1;

    kmds::Mesh m;

    m.updateNodeCapacity(nb_x*nb_y*nb_z);
    m.updateRegionCapacity(nb_i*nb_j*nb_k);

    int indexNode = m.addNodes(nb_x*nb_y*nb_z);
    for(int i=0; i<nb_x; i++) {
        for(int j=0; j<nb_y; j++) {
            for(int k=0; k<nb_z; k++) {
                m.setNodeLocation(indexNode, i, j, k);
                indexNode++;
            }
        }
    }



    kmds::TCellID indexCell = m.addHexahedra(nb_i*nb_j*nb_k);
    for(int i=0; i<nb_i; i++) {
        for(int j=0; j<nb_j; j++) {
            for(int k=0; k<nb_k; k++) {
                kmds::Region h = m.getRegion(indexCell);
                kmds::TCellID v[8];
                v[0] = i * nb_y * nb_z + j * nb_z + k;
                v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                v[4] = i * nb_y * nb_z + j * nb_z + k +1;
                v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k +1;
                v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k +1;
                v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                h.setNodes(v, 8);
                indexCell++;
            }
        }
    }

    kmds::Variable<std::uintptr_t>* v = m.createVariable<std::uintptr_t>(reinterpret_cast<std::uintptr_t>(nullptr), kmds::KMDS_NODE, "geomEntity");
    gmds::cad::FACManager* facGeomManager = new gmds::cad::FACManager;

    elg3d::BoundingBoxGeomAssociation_init_3D(&m, facGeomManager, v);

    kmds::TSize nbNodes = m.getNbNodes();
    kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", nbNodes);
    m.getNodeIDs_dummy(&nodeIDs);

    int nbNonAssociatedNodes = 0;
    int nbAssociatedNodesBottom = 0;
    int nbAssociatedNodesCurve00 = 0;
    for(int i=0; i<nbNodes; i++) {
        kmds::TCellID nid = nodeIDs.get(i);

        if((*v)[nid] == reinterpret_cast<std::uintptr_t>(nullptr)) {
            nbNonAssociatedNodes++;
        } else {

            if ((reinterpret_cast<gmds::cad::GeomEntity *>((*v)[nid])->name() == std::string("surf_bottom"))) {
                nbAssociatedNodesBottom++;
            }

            if ((reinterpret_cast<gmds::cad::GeomEntity *>((*v)[nid])->name() == std::string("curve_0"))) {
                nbAssociatedNodesCurve00++;
            }
        }
    }

    EXPECT_EQ(nbNonAssociatedNodes, (nb_x-2)*(nb_y-2)*(nb_z-2));
    EXPECT_EQ(nbAssociatedNodesBottom, (nb_x-2)*(nb_y-2));
    EXPECT_EQ(nbAssociatedNodesCurve00, (nb_x-2));
}
/*----------------------------------------------------------------------------*/