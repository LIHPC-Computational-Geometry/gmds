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
class ConnectivityTest : public ::testing::Test
{
protected:
    ConnectivityTest()
    {
        ;
    }
    virtual ~ConnectivityTest()
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
TEST_F(ConnectivityTest, N2F_quad)
{
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx -1;
    const int nj = ny -1;

    m.updateNodeCapacity(nx*ny);
    m.updateFaceCapacity(ni*nj);

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }

    for(int i=0; i<ni; i++) {
        for(int j=0; j<nj; j++) {

            m.newQuad(i*ny+j,
                      (i+1)*ny+j,
                      (i+1)*ny+(j+1),
                      i*ny+(j+1)
            );

        }
    }

    kmds::Connectivity* N2F = m.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2F();

    kmds::Node nodes[nx*ny];
    for(int id=0; id<nx*ny; id++) {
        nodes[id] = m.getNode(id);
    }

    Kokkos::View<kmds::TCellID *> ids;
    nodes[0].faceIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(0, ids(0));

    nodes[1].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[2].faceIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(1, ids(0));

    nodes[3].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[4].faceIds(ids);
    EXPECT_EQ(4, ids.size());

    nodes[5].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[6].faceIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(2, ids(0));

    nodes[7].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[8].faceIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(3, ids(0));
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, N2F_quadtri)
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

    kmds::Connectivity* N2F = m.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2F();

    kmds::Node nodes[nx*ny+1];
    for(int id=0; id<nx*ny+1; id++) {
        nodes[id] = m.getNode(id);
    }

    Kokkos::View<kmds::TCellID *> ids;
    nodes[0].faceIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(0, ids(0));

    nodes[1].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[2].faceIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(1, ids(0));

    nodes[3].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[4].faceIds(ids);
    EXPECT_EQ(4, ids.size());

    nodes[5].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[6].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[7].faceIds(ids);
    EXPECT_EQ(4, ids.size());

    nodes[8].faceIds(ids);
    EXPECT_EQ(2, ids.size());

    nodes[9].faceIds(ids);
    EXPECT_EQ(2, ids.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, N2R_hexpyrtetprism)
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

    kmds::Connectivity* N2R = m.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2R();

    kmds::Node nodes[13];
    for (int i = 0; i < 13; i++) {
        nodes[i] = m.getNode(nodeIDs[i]);
    }

    Kokkos::View<kmds::TCellID *> ids;
    nodes[0].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(hexId, ids(0));
    nodes[1].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(hexId, ids(0));
    nodes[2].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(hexId, ids(0));
    nodes[3].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(hexId, ids(0));
    nodes[4].regionIds(ids);
    EXPECT_EQ(3, ids.size());
    nodes[5].regionIds(ids);
    EXPECT_EQ(3, ids.size());
    nodes[6].regionIds(ids);
    EXPECT_EQ(3, ids.size());
    nodes[7].regionIds(ids);
    EXPECT_EQ(3, ids.size());
    nodes[8].regionIds(ids);
    EXPECT_EQ(3, ids.size());
    nodes[9].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(tetId, ids(0));
    nodes[10].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(prismId, ids(0));
    nodes[11].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(prismId, ids(0));
    nodes[12].regionIds(ids);
    EXPECT_EQ(1, ids.size());
    EXPECT_EQ(prismId, ids(0));
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, R2R_byN_hexpyrtetprism)
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

    kmds::Connectivity* N2R = m.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2R();

    kmds::Connectivity* R2R_N = m.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(N2R);

    Kokkos::View<kmds::TCellID*> neighbors;
    R2R_N->get(hexId, neighbors);

    std::set<kmds::TCellID> neighborsSet;
    for(int i=0; i<neighbors.size(); i++) {
        neighborsSet.insert(neighbors(i));
    }

    EXPECT_EQ(neighborsSet.size(), 3);
    EXPECT_TRUE(neighborsSet.find(pyrId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(tetId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(prismId) != neighborsSet.end());

    R2R_N->get(pyrId, neighbors);
    neighborsSet.clear();
    for(int i=0; i<neighbors.size(); i++) {
        neighborsSet.insert(neighbors(i));
    }
    EXPECT_EQ(neighborsSet.size(), 3);
    EXPECT_TRUE(neighborsSet.find(hexId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(tetId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(prismId) != neighborsSet.end());

    R2R_N->get(prismId, neighbors);
    neighborsSet.clear();
    for(int i=0; i<neighbors.size(); i++) {
        neighborsSet.insert(neighbors(i));
    }
    EXPECT_EQ(neighborsSet.size(), 3);
    EXPECT_TRUE(neighborsSet.find(hexId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(tetId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(pyrId) != neighborsSet.end());

    R2R_N->get(tetId, neighbors);
    neighborsSet.clear();
    for(int i=0; i<neighbors.size(); i++) {
        neighborsSet.insert(neighbors(i));
    }
    EXPECT_EQ(neighborsSet.size(), 3);
    EXPECT_TRUE(neighborsSet.find(hexId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(pyrId) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(prismId) != neighborsSet.end());
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, R2R_byN_grid)
{
    kmds::Mesh m;

    const int nx = 4;
    const int ny = 4;
    const int nz = 4;

    const int ni = nx -1;
    const int nj = ny -1;
    const int nk = nz -1;

    m.updateNodeCapacity(nx*ny*nz);
    m.updateRegionCapacity(ni*nj*nk);

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                int a = m.addNode();
                m.setNodeLocation(a, i, j, k);
            }
        }
    }

    for(int i=0; i<ni; i++) {
        for(int j=0; j<nj; j++) {
            for(int k=0; k<nk; k++) {
                m.newHexahedron(i*ny*nz+j*nz+k,
                                (i+1)*ny*nz+j*nz+k,
                                (i+1)*ny*nz+(j+1)*nz+k,
                                i*ny*nz+(j+1)*nz+k,
                                i*ny*nz+j*nz+k+1,
                                (i+1)*ny*nz+j*nz+k+1,
                                (i+1)*ny*nz+(j+1)*nz+k+1,
                                i*ny*nz+(j+1)*nz+k+1

                );
            }
        }
    }

    kmds::Connectivity* N2R = m.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2R();

    kmds::Connectivity* R2R_N = m.createConnectivity(kmds::R2R_byN);
    ch.buildR2R_byN(N2R);

    Kokkos::View<kmds::TCellID*> neighbors;
    R2R_N->get(0, neighbors);

    std::set<kmds::TCellID> neighborsSet;
    for(int i=0; i<neighbors.size(); i++) {
        neighborsSet.insert(neighbors(i));
    }

    EXPECT_EQ(neighborsSet.size(), 7);
    EXPECT_TRUE(neighborsSet.find(0) == neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(1) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(3) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(4) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(9) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(10) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(12) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(13) != neighborsSet.end());

    R2R_N->get(1, neighbors);
    neighborsSet.clear();
    for(int i=0; i<neighbors.size(); i++) {
        neighborsSet.insert(neighbors(i));
    }
    EXPECT_EQ(neighborsSet.size(), 11);
    EXPECT_TRUE(neighborsSet.find(1) == neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(0) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(2) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(3) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(4) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(5) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(9) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(10) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(11) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(12) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(13) != neighborsSet.end());
    EXPECT_TRUE(neighborsSet.find(14) != neighborsSet.end());

    R2R_N->get(13, neighbors);
    neighborsSet.clear();
    for(int i=0; i<neighbors.size(); i++) {
        neighborsSet.insert(neighbors(i));
    }
    EXPECT_EQ(neighborsSet.size(), 26);
    EXPECT_TRUE(neighborsSet.find(13) == neighborsSet.end());
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildEandE2F_grid)
{
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx -1;
    const int nj = ny -1;

    m.updateNodeCapacity(nx*ny);
    m.updateFaceCapacity(ni*nj);

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }

    for(int i=0; i<ni; i++) {
        for(int j=0; j<nj; j++) {

            m.newQuad(i*ny+j,
                      (i+1)*ny+j,
                      (i+1)*ny+(j+1),
                      i*ny+(j+1)
            );

        }
    }

    kmds::Connectivity* c_E2F = m.createConnectivity(kmds::E2F);
    kmds::ConnectivityHelper ch(&m);
    ch.buildEandE2F();

    const kmds::TSize nbEdges = m.getNbEdges();
    kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", nbEdges);
    m.getEdgeIDs(&edgeIDs);

    EXPECT_EQ(nbEdges, 12);

    kmds::TCellID ea_id = kmds::NullID;
    std::vector<kmds::TCellID> ea = {0,1};
    kmds::TCellID eb_id = kmds::NullID;
    std::vector<kmds::TCellID> eb = {0,3};
    kmds::TCellID ec_id = kmds::NullID;
    std::vector<kmds::TCellID> ec = {4,7};
    kmds::TCellID ed_id = kmds::NullID;
    std::vector<kmds::TCellID> ed = {4,5};

    for(int i_e=0; i_e<nbEdges; i_e++) {
        kmds::TCellID eid =  edgeIDs.get(i_e);
        kmds::Edge e = m.getEdge(eid);

        kmds::TCellID nids[2];
        e.nodeIds(nids);

        if(((ea[0] == nids[0]) && (ea[1] == nids[1]))
           || ((ea[0] == nids[1]) && (ea[1] == nids[0]))) {
            ea_id = e.id;
        }

        if(((eb[0] == nids[0]) && (eb[1] == nids[1]))
           || ((eb[0] == nids[1]) && (eb[1] == nids[0]))) {
            eb_id = e.id;
        }

        if(((ec[0] == nids[0]) && (ec[1] == nids[1]))
           || ((ec[0] == nids[1]) && (ec[1] == nids[0]))) {
            ec_id = e.id;
        }

        if(((ed[0] == nids[0]) && (ed[1] == nids[1]))
           || ((ed[0] == nids[1]) && (ed[1] == nids[0]))) {
            ed_id = e.id;
        }
    }


    Kokkos::View<kmds::TCellID *> cellids;

    c_E2F->get(ea_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_E2F->get(eb_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_E2F->get(ec_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 2) && (cellids(1) == 3)) || ((cellids(0) == 3) && (cellids(1) == 2)));

    c_E2F->get(ed_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 1) && (cellids(1) == 3)) || ((cellids(0) == 3) && (cellids(1) == 1)));
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildEandE2F_2D_variant_0_grid)
{
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx -1;
    const int nj = ny -1;

    m.updateNodeCapacity(nx*ny);
    m.updateFaceCapacity(ni*nj);

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }

    for(int i=0; i<ni; i++) {
        for(int j=0; j<nj; j++) {

            m.newQuad(i*ny+j,
                      (i+1)*ny+j,
                      (i+1)*ny+(j+1),
                      i*ny+(j+1)
            );

        }
    }

    kmds::Connectivity* c_E2F = m.createConnectivity(kmds::E2F);
    kmds::ConnectivityHelper ch(&m);
    ch.buildEandE2F_2D_variant_0();

    const kmds::TSize nbEdges = m.getNbEdges();
    kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", nbEdges);
    m.getEdgeIDs(&edgeIDs);

    EXPECT_EQ(nbEdges, 12);

    kmds::TCellID ea_id = kmds::NullID;
    std::vector<kmds::TCellID> ea = {0,1};
    kmds::TCellID eb_id = kmds::NullID;
    std::vector<kmds::TCellID> eb = {0,3};
    kmds::TCellID ec_id = kmds::NullID;
    std::vector<kmds::TCellID> ec = {4,7};
    kmds::TCellID ed_id = kmds::NullID;
    std::vector<kmds::TCellID> ed = {4,5};

    for(int i_e=0; i_e<nbEdges; i_e++) {
        kmds::TCellID eid =  edgeIDs.get(i_e);
        kmds::Edge e = m.getEdge(eid);

        kmds::TCellID nids[2];
        e.nodeIds(nids);

        if(((ea[0] == nids[0]) && (ea[1] == nids[1]))
           || ((ea[0] == nids[1]) && (ea[1] == nids[0]))) {
            ea_id = e.id;
        }

        if(((eb[0] == nids[0]) && (eb[1] == nids[1]))
           || ((eb[0] == nids[1]) && (eb[1] == nids[0]))) {
            eb_id = e.id;
        }

        if(((ec[0] == nids[0]) && (ec[1] == nids[1]))
           || ((ec[0] == nids[1]) && (ec[1] == nids[0]))) {
            ec_id = e.id;
        }

        if(((ed[0] == nids[0]) && (ed[1] == nids[1]))
           || ((ed[0] == nids[1]) && (ed[1] == nids[0]))) {
            ed_id = e.id;
        }
    }


    Kokkos::View<kmds::TCellID *> cellids;

    c_E2F->get(ea_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_E2F->get(eb_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_E2F->get(ec_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 2) && (cellids(1) == 3)) || ((cellids(0) == 3) && (cellids(1) == 2)));

    c_E2F->get(ed_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 1) && (cellids(1) == 3)) || ((cellids(0) == 3) && (cellids(1) == 1)));
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildEandE2F_2D_variant_1_grid)
{
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx -1;
    const int nj = ny -1;

    m.updateNodeCapacity(nx*ny);
    m.updateFaceCapacity(ni*nj);

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            int a = m.addNode();
            m.setNodeLocation(a, i, j, 0.);
        }
    }

    for(int i=0; i<ni; i++) {
        for(int j=0; j<nj; j++) {

            m.newQuad(i*ny+j,
                      (i+1)*ny+j,
                      (i+1)*ny+(j+1),
                      i*ny+(j+1)
            );

        }
    }

    kmds::Connectivity* c_E2F = m.createConnectivity(kmds::E2F);
    kmds::ConnectivityHelper ch(&m);
    ch.buildEandE2F_2D_variant_1();

    const kmds::TSize nbEdges = m.getNbEdges();
    kmds::GrowingView<kmds::TCellID> edgeIDs("EDGES", nbEdges);
    m.getEdgeIDs(&edgeIDs);

    EXPECT_EQ(nbEdges, 12);

    kmds::TCellID ea_id = kmds::NullID;
    std::vector<kmds::TCellID> ea = {0,1};
    kmds::TCellID eb_id = kmds::NullID;
    std::vector<kmds::TCellID> eb = {0,3};
    kmds::TCellID ec_id = kmds::NullID;
    std::vector<kmds::TCellID> ec = {4,7};
    kmds::TCellID ed_id = kmds::NullID;
    std::vector<kmds::TCellID> ed = {4,5};

    for(int i_e=0; i_e<nbEdges; i_e++) {
        kmds::TCellID eid =  edgeIDs.get(i_e);
        kmds::Edge e = m.getEdge(eid);

        kmds::TCellID nids[2];
        e.nodeIds(nids);

        if(((ea[0] == nids[0]) && (ea[1] == nids[1]))
           || ((ea[0] == nids[1]) && (ea[1] == nids[0]))) {
            ea_id = e.id;
        }

        if(((eb[0] == nids[0]) && (eb[1] == nids[1]))
           || ((eb[0] == nids[1]) && (eb[1] == nids[0]))) {
            eb_id = e.id;
        }

        if(((ec[0] == nids[0]) && (ec[1] == nids[1]))
           || ((ec[0] == nids[1]) && (ec[1] == nids[0]))) {
            ec_id = e.id;
        }

        if(((ed[0] == nids[0]) && (ed[1] == nids[1]))
           || ((ed[0] == nids[1]) && (ed[1] == nids[0]))) {
            ed_id = e.id;
        }
    }


    Kokkos::View<kmds::TCellID *> cellids;

    c_E2F->get(ea_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_E2F->get(eb_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_E2F->get(ec_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 2) && (cellids(1) == 3)) || ((cellids(0) == 3) && (cellids(1) == 2)));

    c_E2F->get(ed_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 1) && (cellids(1) == 3)) || ((cellids(0) == 3) && (cellids(1) == 1)));
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildFandF2R_grid) {
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;
    const int nz = 3;

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

    kmds::Connectivity* c_F2R = m.createConnectivity(kmds::F2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildFandF2R();

    const kmds::TSize nbFaces = m.getNbFaces();
    kmds::GrowingView<kmds::TCellID> faceIDs("FACES", nbFaces);
    m.getFaceIDs(&faceIDs);

    EXPECT_EQ(nbFaces, 36);

    kmds::TCellID fa_id = kmds::NullID;
    std::vector<kmds::TCellID> fa = {0,3,4,1};
    kmds::TCellID fb_id = kmds::NullID;
    std::vector<kmds::TCellID> fb = {0,9,10,1};
    kmds::TCellID fc_id = kmds::NullID;
    std::vector<kmds::TCellID> fc = {3,12,13,4};
    kmds::TCellID fd_id = kmds::NullID;
    std::vector<kmds::TCellID> fd = {13,22,23,14};

    for(int i_f=0; i_f<nbFaces; i_f++) {
        kmds::TCellID fid =  faceIDs.get(i_f);
        kmds::Face f = m.getFace(fid);

        Kokkos::View<kmds::TCellID*> nids;
        f.nodeIds(nids);

        std::set<kmds::TCellID> nids_set;
        for(int i_n=0; i_n<nids.size(); i_n++) {

            nids_set.insert(nids(i_n));

            // find a specific face

            int found_fa_nb = 0;
            for(auto f_n: fa) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fa_nb++;
                }
            }
            if(found_fa_nb == 4) {
                fa_id = fid;
            }

            int found_fb_nb = 0;
            for(auto f_n: fb) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fb_nb++;
                }
            }
            if(found_fb_nb == 4) {
                fb_id = fid;
            }

            int found_fc_nb = 0;
            for(auto f_n: fc) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fc_nb++;
                }
            }
            if(found_fc_nb == 4) {
                fc_id = fid;
            }

            int found_fd_nb = 0;
            for(auto f_n: fd) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fd_nb++;
                }
            }
            if(found_fd_nb == 4) {
                fd_id = fid;
            }
        }
    }

    Kokkos::View<kmds::TCellID *> cellids;

    c_F2R->get(fa_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_F2R->get(fb_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_F2R->get(fc_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 0) && (cellids(1) == 2)) || ((cellids(0) == 2) && (cellids(1) == 0)));

    c_F2R->get(fd_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 5) && (cellids(1) == 7)) || ((cellids(0) == 7) && (cellids(1) == 5)));
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildFandF2R_variant_0_grid) {
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;
    const int nz = 3;

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

    kmds::Connectivity* c_F2R = m.createConnectivity(kmds::F2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildFandF2R_variant_0();

    const kmds::TSize nbFaces = m.getNbFaces();
    kmds::GrowingView<kmds::TCellID> faceIDs("FACES", nbFaces);
    m.getFaceIDs(&faceIDs);

    EXPECT_EQ(nbFaces, 36);

    kmds::TCellID fa_id = kmds::NullID;
    std::vector<kmds::TCellID> fa = {0,3,4,1};
    kmds::TCellID fb_id = kmds::NullID;
    std::vector<kmds::TCellID> fb = {0,9,10,1};
    kmds::TCellID fc_id = kmds::NullID;
    std::vector<kmds::TCellID> fc = {3,12,13,4};
    kmds::TCellID fd_id = kmds::NullID;
    std::vector<kmds::TCellID> fd = {13,22,23,14};

    for(int i_f=0; i_f<nbFaces; i_f++) {
        kmds::TCellID fid =  faceIDs.get(i_f);
        kmds::Face f = m.getFace(fid);

        Kokkos::View<kmds::TCellID*> nids;
        f.nodeIds(nids);

        std::set<kmds::TCellID> nids_set;
        for(int i_n=0; i_n<nids.size(); i_n++) {

            nids_set.insert(nids(i_n));

            // find a specific face

            int found_fa_nb = 0;
            for(auto f_n: fa) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fa_nb++;
                }
            }
            if(found_fa_nb == 4) {
                fa_id = fid;
            }

            int found_fb_nb = 0;
            for(auto f_n: fb) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fb_nb++;
                }
            }
            if(found_fb_nb == 4) {
                fb_id = fid;
            }

            int found_fc_nb = 0;
            for(auto f_n: fc) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fc_nb++;
                }
            }
            if(found_fc_nb == 4) {
                fc_id = fid;
            }

            int found_fd_nb = 0;
            for(auto f_n: fd) {
                if(nids_set.find(f_n) != nids_set.end()) {
                    found_fd_nb++;
                }
            }
            if(found_fd_nb == 4) {
                fd_id = fid;
            }
        }
    }

    Kokkos::View<kmds::TCellID *> cellids;

    c_F2R->get(fa_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_F2R->get(fb_id, cellids);
    EXPECT_EQ(cellids.size(), 1);
    EXPECT_EQ(cellids(0), 0);

    c_F2R->get(fc_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 0) && (cellids(1) == 2)) || ((cellids(0) == 2) && (cellids(1) == 0)));

    c_F2R->get(fd_id, cellids);
    EXPECT_EQ(cellids.size(), 2);
    EXPECT_TRUE(((cellids(0) == 5) && (cellids(1) == 7)) || ((cellids(0) == 7) && (cellids(1) == 5)));
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildN2F_grid) {
    kmds::Mesh m;

    const int nx = 4;
    const int ny = 4;
    const int nz = 4;

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

    kmds::Connectivity* c_F2R = m.createConnectivity(kmds::F2R);
    kmds::ConnectivityHelper ch(&m);
    ch.buildFandF2R_variant_0();
    kmds::Connectivity* c_N2F = m.createConnectivity(kmds::N2F);
    ch.buildN2F();

    std::vector<kmds::FakeFace> fakefacevec;
    std::map<kmds::FakeFace, int> fakefaceset;
    std::map<kmds::TCellID, std::set<kmds::FakeFace> > node2ff;

    int ffindex = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nk; k++) {

                kmds::TCellID n0 = i * ny * nz + j * nz + k;
                kmds::TCellID n1 = i * ny * nz + (j + 1) * nz + k;
                kmds::TCellID n2 = i * ny * nz + (j + 1) * nz + k + 1;
                kmds::TCellID n3 = i * ny * nz + j * nz + k + 1;

                std::vector<kmds::TCellID> nids = {n0,n1,n2,n3};
                kmds::FakeFace ff(nids);
                fakefacevec.push_back(ff);
                fakefaceset.emplace(ff, ffindex);
                ffindex++;
                node2ff[n0].insert(ff);
                node2ff[n1].insert(ff);
                node2ff[n2].insert(ff);
                node2ff[n3].insert(ff);
            }
        }
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nk; k++) {

                kmds::TCellID n0 = i * ny * nz + j * nz + k;
                kmds::TCellID n1 = (i + 1) * ny * nz + j * nz + k;
                kmds::TCellID n2 = (i + 1) * ny * nz + j * nz + k + 1;
                kmds::TCellID n3 = i * ny * nz + j * nz + k + 1;

                std::vector<kmds::TCellID> nids = {n0,n1,n2,n3};
                kmds::FakeFace ff(nids);
                fakefacevec.push_back(ff);
                fakefaceset.emplace(ff, ffindex);
                ffindex++;
                node2ff[n0].insert(ff);
                node2ff[n1].insert(ff);
                node2ff[n2].insert(ff);
                node2ff[n3].insert(ff);
            }
        }
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < nz; k++) {

                kmds::TCellID n0 = i * ny * nz + j * nz + k;
                kmds::TCellID n1 = (i + 1) * ny * nz + j * nz + k;
                kmds::TCellID n2 = (i + 1) * ny * nz + (j + 1) * nz + k;
                kmds::TCellID n3 = i * ny * nz + (j + 1) * nz + k;

                std::vector<kmds::TCellID> nids = {n0,n1,n2,n3};
                kmds::FakeFace ff(nids);
                fakefacevec.push_back(ff);
                fakefaceset.emplace(ff, ffindex);
                ffindex++;
                node2ff[n0].insert(ff);
                node2ff[n1].insert(ff);
                node2ff[n2].insert(ff);
                node2ff[n3].insert(ff);
            }
        }
    }

    EXPECT_EQ(fakefacevec.size(), 108);
    EXPECT_EQ(fakefaceset.size(), 108);


    const kmds::TSize nbFaces = m.getNbFaces();
    kmds::GrowingView<kmds::TCellID> faceIDs("FACES", nbFaces);
    m.getFaceIDs(&faceIDs);

    EXPECT_EQ(nbFaces, 108);

    for(int i_f=0; i_f<nbFaces; i_f++) {
        kmds::TCellID fid = faceIDs.get(i_f);
        kmds::Face f = m.getFace(fid);

        Kokkos::View<kmds::TCellID*> nids;
        f.nodeIds(nids);

        std::vector<kmds::TCellID> nids_tmp;
        for(int i_n=0; i_n<nids.size(); i_n++) {
            nids_tmp.push_back(nids(i_n));
        }

        kmds::FakeFace ff(nids_tmp);

        for(int i_n=0; i_n<nids.size(); i_n++) {

            kmds::TCellID nid = nids(i_n);
            EXPECT_FALSE(node2ff[nid].find(ff) == node2ff[nid].end());

        }
    }

    const kmds::TSize nbNodes = m.getNbNodes();
    kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", nbNodes);
    m.getNodeIDs(&nodeIDs);

    for(int i_n=0; i_n<nbNodes; i_n++) {
        kmds::TCellID nid = nodeIDs.get(i_n);

        Kokkos::View<kmds::TCellID *> fids;
        c_N2F->get(nid, fids);

        for(int i_f=0; i_f<fids.size(); i_f++) {

            kmds::Face f = m.getFace(fids(i_f));
            Kokkos::View<kmds::TCellID*> nids;
            f.nodeIds(nids);

            std::vector<kmds::TCellID> nids_tmp;
            for(int i_nbis=0; i_nbis<nids.size(); i_nbis++) {
                nids_tmp.push_back(nids(i_nbis));
            }

            kmds::FakeFace ff(nids_tmp);

            EXPECT_FALSE(node2ff[nid].find(ff) == node2ff[nid].end());
        }
    }

}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildN2N_quadtri_2D) {
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;

    const int ni = nx - 1;
    const int nj = ny - 1;

    m.updateNodeCapacity(nx * ny + 1);
    m.updateFaceCapacity(ni * nj + 2);

    kmds::TCellID nodeids[nx * ny + 1];
    int index = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            nodeids[index] = m.newNode(i, j, 0.);
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
    m.newTriangle((nx - 1) * ny,
                  a,
                  (nx - 1) * ny + 1
    );
    m.newTriangle((nx - 1) * ny + 1,
                  a,
                  (nx - 1) * ny + 2
    );

    kmds::Connectivity *c_N2N = m.createConnectivity(kmds::N2N);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2N(kmds::F);

    const kmds::TSize nbNodes = m.getNbNodes();
    kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", nbNodes);
    m.getNodeIDs(&nodeIDs);


    std::set<kmds::TCellID> neighbors_0 = {1,3};
    std::set<kmds::TCellID> neighbors_1 = {0,2,4};
    std::set<kmds::TCellID> neighbors_4 = {1,3,5,7};
    std::set<kmds::TCellID> neighbors_9 = {6,7,8};

    int nbneighbors_0 = 0;
    int nbneighbors_1 = 0;
    int nbneighbors_4 = 0;
    int nbneighbors_9 = 0;

    Kokkos::View<kmds::TCellID *> cellids;

    c_N2N->get(0, cellids);
    EXPECT_EQ(cellids.size(), neighbors_0.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_0.find(cellids(i)) != neighbors_0.end()) {
            nbneighbors_0++;
        }
    }
    EXPECT_EQ(nbneighbors_0, neighbors_0.size());

    c_N2N->get(1, cellids);
    EXPECT_EQ(cellids.size(), neighbors_1.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_1.find(cellids(i)) != neighbors_1.end()) {
            nbneighbors_1++;
        }
    }
    EXPECT_EQ(nbneighbors_1, neighbors_1.size());

    c_N2N->get(4, cellids);
    EXPECT_EQ(cellids.size(), neighbors_4.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_4.find(cellids(i)) != neighbors_4.end()) {
            nbneighbors_4++;
        }
    }
    EXPECT_EQ(nbneighbors_4, neighbors_4.size());

    c_N2N->get(9, cellids);
    EXPECT_EQ(cellids.size(), neighbors_9.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_9.find(cellids(i)) != neighbors_9.end()) {
            nbneighbors_9++;
        }
    }
    EXPECT_EQ(nbneighbors_9, neighbors_9.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(ConnectivityTest, BuildN2N_grid_3D) {
    kmds::Mesh m;

    const int nx = 3;
    const int ny = 3;
    const int nz = 3;

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

    kmds::Connectivity *c_N2N = m.createConnectivity(kmds::N2N);
    kmds::ConnectivityHelper ch(&m);
    ch.buildN2N(kmds::R);

    const kmds::TSize nbNodes = m.getNbNodes();
    kmds::GrowingView<kmds::TCellID> nodeIDs("NODES", nbNodes);
    m.getNodeIDs(&nodeIDs);


    std::set<kmds::TCellID> neighbors_0 = {1,3,9};
    std::set<kmds::TCellID> neighbors_1 = {0,2,4,10};
    std::set<kmds::TCellID> neighbors_4 = {1,3,5,7,13};
    std::set<kmds::TCellID> neighbors_13 = {4,10,12,14,16,22};

    int nbneighbors_0 = 0;
    int nbneighbors_1 = 0;
    int nbneighbors_4 = 0;
    int nbneighbors_13 = 0;

    Kokkos::View<kmds::TCellID *> cellids;

    c_N2N->get(0, cellids);
    EXPECT_EQ(cellids.size(), neighbors_0.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_0.find(cellids(i)) != neighbors_0.end()) {
            nbneighbors_0++;
        }
    }
    EXPECT_EQ(nbneighbors_0, neighbors_0.size());

    c_N2N->get(1, cellids);
    EXPECT_EQ(cellids.size(), neighbors_1.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_1.find(cellids(i)) != neighbors_1.end()) {
            nbneighbors_1++;
        }
    }
    EXPECT_EQ(nbneighbors_1, neighbors_1.size());

    c_N2N->get(4, cellids);
    EXPECT_EQ(cellids.size(), neighbors_4.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_4.find(cellids(i)) != neighbors_4.end()) {
            nbneighbors_4++;
        }
    }
    EXPECT_EQ(nbneighbors_4, neighbors_4.size());

    c_N2N->get(13, cellids);
    EXPECT_EQ(cellids.size(), neighbors_13.size());
    for(int i=0; i<cellids.size(); i++) {
        if(neighbors_13.find(cellids(i)) != neighbors_13.end()) {
            nbneighbors_13++;
        }
    }
    EXPECT_EQ(nbneighbors_13, neighbors_13.size());
}
/*----------------------------------------------------------------------------*/