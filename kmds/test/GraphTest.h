/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Google Test headers
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS headers
#include <KM/Utils/Graph.h>
/*----------------------------------------------------------------------------*/
// using namespace kmds;
/*----------------------------------------------------------------------------*/
class GraphTest : public ::testing::Test
{
protected:
    GraphTest()
    {
        ;
    }
    virtual ~GraphTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
		  Kokkos::InitializationSettings kargs;
		  kargs.set_num_threads(1);
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
TEST_F(GraphTest, creation)
{
    kmds::Graph gr("GRAPH_TEST", 1000, 20);

    kmds::GrowingView<kmds::TCellID> indset("INDEPENDENT_SET", gr.getNbVec());
    gr.getIndependentSet(&indset);

    EXPECT_EQ(indset.getNbElems(), 1000);
    EXPECT_TRUE(gr.checkIndependentSet(&indset));
}
/*----------------------------------------------------------------------------*/
TEST_F(GraphTest, populate)
{
    kmds::Graph gr("GRAPH_TEST", 1000, 20);

    gr.populateNonOriented(10.);

    kmds::GrowingView<kmds::TCellID> indset("INDEPENDENT_SET", gr.getNbVec());
    gr.getIndependentSet(&indset);

    EXPECT_TRUE(gr.checkIndependentSet(&indset));
}
/*----------------------------------------------------------------------------*/
TEST_F(GraphTest, coloring) {

    kmds::Graph gr("GRAPH_TEST", 20, 20);

    gr.populateNonOriented(10.);

    const int nbMaxColors = 20;
    gr.buildColoring(nbMaxColors);


    std::map<kmds::TCellID, int> vec2nb;

    for(int i=0; i<nbMaxColors; i++) {

        kmds::GrowingView<kmds::TCellID>* indset = gr.getColoring(i);
        EXPECT_TRUE(gr.checkIndependentSet(indset));

        std::cout<<"indset->getNbElems() "<<indset->getNbElems()<<std::endl;

        for(int j=0; j<indset->getNbElems(); j++) {

            kmds::TCellID vid = indset->get(j);
            if(vec2nb.find(vid) != vec2nb.end()) {
                vec2nb[vid] = vec2nb[vid] + 1;
            } else {
                vec2nb.emplace(vid, 1);
            }
        }
    }

    for(int i=0; i<gr.getNbVec(); i++) {

        EXPECT_TRUE(vec2nb.find(i) != vec2nb.end());
        EXPECT_TRUE(vec2nb[i] == 1);
    }
}
/*----------------------------------------------------------------------------*/
TEST_F(GraphTest, neighbors)
{
    kmds::Graph gr("GRAPH_TEST", 5, 20);
    EXPECT_EQ(gr.getNbVec(), 5);
    EXPECT_EQ(gr.getNbEdges(), 0);

    std::set<kmds::TCellID> nids0 = {1,4};
    gr.setNeighbors(0, nids0);
    EXPECT_EQ(gr.getNbEdges(), 2);
    EXPECT_TRUE(gr.isNeighbor(0,4));
    EXPECT_FALSE(gr.isNeighbor(0,3));

    std::vector<kmds::TCellID> nids3 = {0,1,4};
    gr.setNeighbors(3, nids3);
    EXPECT_EQ(gr.getNbEdges(), 5);
    EXPECT_TRUE(gr.isNeighbor(3,0));
    EXPECT_FALSE(gr.isNeighbor(0,3));

    Kokkos::View<kmds::TCellID*> nids2_get;
    gr.getNeighbors(2, nids2_get);
    Kokkos::View<kmds::TCellID*> nids4_get;
    gr.getNeighbors(4, nids4_get);
    EXPECT_EQ(nids2_get.extent(0), 0);
    EXPECT_EQ(nids4_get.extent(0), 0);

    Kokkos::View<kmds::TCellID*> nids3_get;
    gr.getNeighbors(3, nids3_get);
    EXPECT_EQ(nids3_get.extent(0), 3);

    std::set<kmds::TCellID> nids3_set;
    for(int i_n=0; i_n<nids3_get.extent(0); i_n++) {
        nids3_set.insert(nids3_get(i_n));
    }
    EXPECT_EQ(nids3_set.size(), 3);
    EXPECT_TRUE(nids3_set.find(0) != nids3_set.end());
    EXPECT_TRUE(nids3_set.find(1) != nids3_set.end());
    EXPECT_TRUE(nids3_set.find(4) != nids3_set.end());


    EXPECT_TRUE(gr.isNeighbor(3,1));
    std::vector<kmds::TCellID> nids3_bis = {0,4};
    gr.setNeighbors(3, nids3_bis);
    EXPECT_EQ(gr.getNbEdges(), 4);
    EXPECT_TRUE(gr.isNeighbor(3,0));
    EXPECT_FALSE(gr.isNeighbor(0,3));
    EXPECT_FALSE(gr.isNeighbor(3,1));
}
/*----------------------------------------------------------------------------*/