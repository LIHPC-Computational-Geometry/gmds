/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <GCoptimization.h>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

double gco_smoothCost_(int p, int q, int label_p, int label_q,void *extraData){
    if (label_q == label_p){
        return 0.;
    } else{
        return 1.;
    }
}

double gco_smoothCost_interface_(int p, int q, int label_p, int label_q,void *extraData){

    gmds::math::Vector n(1., -0.5);
    n.normalize();

    std::map<int, gmds::math::Point>* index2pts = reinterpret_cast<std::map<int, gmds::math::Point>* > (extraData);
    gmds::math::Point pt_p = (*index2pts)[p];
    gmds::math::Point pt_q = (*index2pts)[q];

    gmds::math::Vector v(pt_p, pt_q);
    v.normalize();

    double dotproduct = n.dot(v);

    if (label_q == label_p){
        return 0.;
    } else{
        return 1 - std::fabs(dotproduct);
    }
}

/*----------------------------------------------------------------------------*/
class GCOTest : public ::testing::Test
{
 protected:
        GCOTest()
        {
                ;
        }
        virtual ~GCOTest()
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
TEST_F(GCOTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(GCOTest, smallGraph)
{
    const int nbVert = 4;
    const int nbColor = 2;

    GCoptimizationGeneralGraph graphCut(nbVert, nbColor);

    for(int i=0; i<nbVert-1; i++) {
        graphCut.setNeighbors(i, i+1);
    }

    std::vector<double> dataCost(nbVert*nbColor, 1000000000.);
    dataCost[0] = 0;
    dataCost[1*nbColor] = 0;
    dataCost[2*nbColor+1] = 0;
    dataCost[3*nbColor+1] = 0;


    graphCut.setDataCost(dataCost.data());

    graphCut.setSmoothCost(&gco_smoothCost_, NULL);

    try { //graphCut.expansion(300);
        graphCut.swap();
    } catch (const GCException& e) {std::cout << e.message << '\n'; }

    EXPECT_EQ(0, graphCut.whatLabel(0));
    EXPECT_EQ(0, graphCut.whatLabel(1));
    EXPECT_EQ(1, graphCut.whatLabel(2));
    EXPECT_EQ(1, graphCut.whatLabel(3));
    EXPECT_NE(1, graphCut.whatLabel(0));
    EXPECT_NE(1, graphCut.whatLabel(1));
    EXPECT_NE(0, graphCut.whatLabel(2));
    EXPECT_NE(0, graphCut.whatLabel(3));
}
/*----------------------------------------------------------------------------*/
TEST_F(GCOTest, smallGraph_1)
{
    const int nbVert = 4;
    const int nbColor = 2;

    GCoptimizationGeneralGraph graphCut(nbVert, nbColor);

    for(int i=0; i<nbVert-1; i++) {
        graphCut.setNeighbors(i, i+1);
    }

    std::vector<double> dataCost(nbVert*nbColor, 1000000000.);
    dataCost[0] = 0;
    dataCost[1*nbColor] = 0.25;
    dataCost[1*nbColor+1] = 0.75;
    dataCost[2*nbColor] = 0.75;
    dataCost[2*nbColor+1] = 0.25;
    dataCost[3*nbColor+1] = 0;


    graphCut.setDataCost(dataCost.data());

    graphCut.setSmoothCost(&gco_smoothCost_, NULL);

    try { //graphCut.expansion(300);
        graphCut.swap();
    } catch (const GCException& e) {std::cout << e.message << '\n'; }

    EXPECT_EQ(0, graphCut.whatLabel(0));
    EXPECT_EQ(0, graphCut.whatLabel(1));
    EXPECT_EQ(1, graphCut.whatLabel(2));
    EXPECT_EQ(1, graphCut.whatLabel(3));
    EXPECT_NE(1, graphCut.whatLabel(0));
    EXPECT_NE(1, graphCut.whatLabel(1));
    EXPECT_NE(0, graphCut.whatLabel(2));
    EXPECT_NE(0, graphCut.whatLabel(3));

    EXPECT_DOUBLE_EQ(1.5, graphCut.compute_energy());
    EXPECT_DOUBLE_EQ(0.5, graphCut.giveDataEnergy());
    EXPECT_DOUBLE_EQ(1., graphCut.giveSmoothEnergy());
}
/*----------------------------------------------------------------------------*/
TEST_F(GCOTest, smallGraph_2) {

    const int nsub = 10;

    const int nbVert = nsub * nsub + 3 * nsub;
    const int nbColor = 2;

    GCoptimizationGeneralGraph graphCut(nbVert, nbColor);

    // init (i,j) -> index
    int indices[nsub+2][nsub+2];

    for(int i=0; i<nsub+2; i++) {
        for(int j=0; j<nsub+2; j++) {
            indices[i][j] = -1;
        }
    }
    int index = 0;
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            indices[i][j] = index;
            index++;
        }
    }
    // left
    for(int j=1; j<=nsub; j++) {
        indices[0][j] = index;
        index++;
    }
    //right
    for(int j=1; j<=nsub; j++) {
        indices[nsub+1][j] = index;
        index++;
    }
    //top
    for(int i=1; i<=nsub; i++) {
        indices[i][nsub+1] = index;
        index++;
    }

    const double fpA = 0.7;
    const double fpB = 0.3;

    std::vector<double> dataCost(nbVert*nbColor, 1000000.);
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            dataCost[indices[i][j] * nbColor] = 1. - fpA;
            dataCost[indices[i][j] * nbColor + 1] = 1. - fpB;
        }
    }
    // left A
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[0][j] * nbColor] = 0.;
    }
    //right B
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[nsub+1][j] * nbColor + 1] = 0.;
    }
    //top A
    for(int i=1; i<=nsub; i++) {
        dataCost[indices[i][nsub+1] * nbColor] = 0.;
    }


    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            const int current_index = indices[i][j];
            const int right_index = indices[i+1][j];
            const int top_index = indices[i][j+1];

            graphCut.setNeighbors(current_index, right_index);
            graphCut.setNeighbors(current_index, top_index);
        }
    }
    for(int j=1; j<=nsub; j++) {
        const int current_index = indices[0][j];
        const int right_index = indices[1][j];

        graphCut.setNeighbors(current_index, right_index);
    }

    graphCut.setDataCost(dataCost.data());

    graphCut.setSmoothCost(&gco_smoothCost_, NULL);

    try { //graphCut.expansion(300);
        graphCut.swap();
    } catch (const GCException& e) {std::cout << e.message << '\n'; }


    std::cout<<"=============="<<std::endl;

    for(int j=nsub+1; j>=0; j--) {
        for(int i=0; i<nsub+2; i++) {

            if(indices[i][j] == -1) {
                std::cout<<" ";
            } else {
                std::cout<<graphCut.whatLabel(indices[i][j]);
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<"=============="<<std::endl;

    // cost fp-wise
    int nblabel_A = 0;
    int nblabel_B = 0;

    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            if(graphCut.whatLabel(indices[i][j]) == 0) {
                nblabel_A++;
            }
            if(graphCut.whatLabel(indices[i][j]) == 1) {
                nblabel_B++;
            }
        }
    }
    std::cout<<"fpA "<<(double)nblabel_A/(nsub*nsub)<<" fpB "<<(double)nblabel_B/(nsub*nsub)<<std::endl;

    std::cout<<"energy tot    "<<graphCut.compute_energy()<<std::endl;
    std::cout<<"energy data   "<<graphCut.giveDataEnergy()<<std::endl;
    std::cout<<"energy smooth "<<graphCut.giveSmoothEnergy()<<std::endl;
}
/*----------------------------------------------------------------------------*/
TEST_F(GCOTest, smallGraph_3) {

    const int nsub = 10;

    const int nbVert = nsub * nsub + 3 * nsub;
    const int nbColor = 2;

    GCoptimizationGeneralGraph graphCut(nbVert, nbColor);

    std::map<int, gmds::math::Point> pts_A;
    std::map<int, gmds::math::Point> pts_B;

    // init (i,j) -> index
    int indices[nsub+2][nsub+2];

    for(int i=0; i<nsub+2; i++) {
        for(int j=0; j<nsub+2; j++) {
            indices[i][j] = -1;
        }
    }
    int index = 0;
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            indices[i][j] = index;
            index++;
        }
    }
    // left
    for(int j=1; j<=nsub; j++) {
        indices[0][j] = index;
        pts_A.emplace(index, gmds::math::Point(0.,j,0.));
        index++;
    }
    //right
    for(int j=1; j<=nsub; j++) {
        indices[nsub+1][j] = index;
        pts_B.emplace(index, gmds::math::Point(nsub+1,j,0.));
        index++;
    }
    //top
    for(int i=1; i<=nsub; i++) {
        indices[i][nsub+1] = index;
        pts_A.emplace(index, gmds::math::Point(i,nsub+1,0.));
        index++;
    }

    const double fpA = 0.7;
    const double fpB = 0.3;

    std::vector<double> dataCost(nbVert*nbColor, 1000000.);
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            const gmds::math::Point pt(i,j,0.);
            double minlength_A = HUGE_VALF;
            int index_a = -1;
            for(auto pa: pts_A) {
                if(pt.distance(pa.second) < minlength_A) {
                    minlength_A = pt.distance(pa.second);
                    index_a = pa.first;
                }
            }
            double minlength_B = HUGE_VALF;
            int index_b = -1;
            for(auto pb: pts_B) {
                if(pt.distance(pb.second) < minlength_B) {
                    minlength_B = pt.distance(pb.second);
                    index_b = pb.first;
                }
            }

            dataCost[indices[i][j] * nbColor] = (1. - fpA)*minlength_A;
            dataCost[indices[i][j] * nbColor + 1] = (1. - fpB)*minlength_B;
        }
    }
    // left A
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[0][j] * nbColor] = 0.;
    }
    //right B
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[nsub+1][j] * nbColor + 1] = 0.;
    }
    //top A
    for(int i=1; i<=nsub; i++) {
        dataCost[indices[i][nsub+1] * nbColor] = 0.;
    }


    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            const int current_index = indices[i][j];
            const int right_index = indices[i+1][j];
            const int top_index = indices[i][j+1];

            graphCut.setNeighbors(current_index, right_index);
            graphCut.setNeighbors(current_index, top_index);
        }
    }
    for(int j=1; j<=nsub; j++) {
        const int current_index = indices[0][j];
        const int right_index = indices[1][j];

        graphCut.setNeighbors(current_index, right_index);
    }

    graphCut.setDataCost(dataCost.data());

    graphCut.setSmoothCost(&gco_smoothCost_, NULL);

    try { //graphCut.expansion(300);
        graphCut.swap();
    } catch (const GCException& e) {std::cout << e.message << '\n'; }


    std::cout<<"=============="<<std::endl;

    for(int j=nsub+1; j>=0; j--) {
        for(int i=0; i<nsub+2; i++) {

            if(indices[i][j] == -1) {
                std::cout<<" ";
            } else {
                std::cout<<graphCut.whatLabel(indices[i][j]);
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<"=============="<<std::endl;

    // cost fp-wise
    int nblabel_A = 0;
    int nblabel_B = 0;

    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            if(graphCut.whatLabel(indices[i][j]) == 0) {
                nblabel_A++;
            }
            if(graphCut.whatLabel(indices[i][j]) == 1) {
                nblabel_B++;
            }
        }
    }
    std::cout<<"fpA "<<(double)nblabel_A/(nsub*nsub)<<" fpB "<<(double)nblabel_B/(nsub*nsub)<<std::endl;

    std::cout<<"energy tot    "<<graphCut.compute_energy()<<std::endl;
    std::cout<<"energy data   "<<graphCut.giveDataEnergy()<<std::endl;
    std::cout<<"energy smooth "<<graphCut.giveSmoothEnergy()<<std::endl;
}
/*----------------------------------------------------------------------------*/
TEST_F(GCOTest, smallGraph_3_bis) {

    const int nsub = 10;

    const int nbVert = nsub * nsub + 2 * nsub;
    const int nbColor = 2;

    GCoptimizationGeneralGraph graphCut(nbVert, nbColor);

    std::map<int, gmds::math::Point> pts_A;
    std::map<int, gmds::math::Point> pts_B;

    // init (i,j) -> index
    int indices[nsub+2][nsub+2];

    for(int i=0; i<nsub+2; i++) {
        for(int j=0; j<nsub+2; j++) {
            indices[i][j] = -1;
        }
    }
    int index = 0;
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            indices[i][j] = index;
            index++;
        }
    }
    // left
    for(int j=1; j<=nsub; j++) {
        indices[0][j] = index;
        pts_A.emplace(index, gmds::math::Point(0.,j,0.));
        index++;
    }
    //right
    for(int j=1; j<=nsub; j++) {
        indices[nsub+1][j] = index;
        pts_B.emplace(index, gmds::math::Point(nsub+1,j,0.));
        index++;
    }

    const double fpA = 0.7;
    const double fpB = 0.3;

    std::vector<double> dataCost(nbVert*nbColor, 1000000.);
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            const gmds::math::Point pt(i,j,0.);
            double minlength_A = HUGE_VALF;
            int index_a = -1;
            for(auto pa: pts_A) {
                if(pt.distance(pa.second) < minlength_A) {
                    minlength_A = pt.distance(pa.second);
                    index_a = pa.first;
                }
            }
            double minlength_B = HUGE_VALF;
            int index_b = -1;
            for(auto pb: pts_B) {
                if(pt.distance(pb.second) < minlength_B) {
                    minlength_B = pt.distance(pb.second);
                    index_b = pb.first;
                }
            }

            dataCost[indices[i][j] * nbColor] = (1. - fpA)*minlength_A;
            dataCost[indices[i][j] * nbColor + 1] = (1. - fpB)*minlength_B;

            if(i==8) {
                std::cout<<"dcostA "<<dataCost[indices[i][j] * nbColor]<<std::endl;
                std::cout<<"dcostB "<<dataCost[indices[i][j] * nbColor + 1]<<std::endl;
            }
        }
    }
    // left A
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[0][j] * nbColor] = 0.;
    }
    //right B
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[nsub+1][j] * nbColor + 1] = 0.;
    }


    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<nsub; j++) {
            const int current_index = indices[i][j];
            const int right_index = indices[i+1][j];
            const int top_index = indices[i][j+1];

            graphCut.setNeighbors(current_index, right_index);
            graphCut.setNeighbors(current_index, top_index);
        }
    }
    for(int j=1; j<=nsub; j++) {
        const int current_index = indices[0][j];
        const int right_index = indices[1][j];

        graphCut.setNeighbors(current_index, right_index);
    }

    graphCut.setDataCost(dataCost.data());

    graphCut.setSmoothCost(&gco_smoothCost_, NULL);

    try { //graphCut.expansion(300);
        graphCut.swap();
    } catch (const GCException& e) {std::cout << e.message << '\n'; }


    std::cout<<"=============="<<std::endl;

    for(int j=nsub+1; j>=0; j--) {
        for(int i=0; i<nsub+2; i++) {

            if(indices[i][j] == -1) {
                std::cout<<" ";
            } else {
                std::cout<<graphCut.whatLabel(indices[i][j]);
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<"=============="<<std::endl;

    // cost fp-wise
    int nblabel_A = 0;
    int nblabel_B = 0;

    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            if(graphCut.whatLabel(indices[i][j]) == 0) {
                nblabel_A++;
            }
            if(graphCut.whatLabel(indices[i][j]) == 1) {
                nblabel_B++;
            }
        }
    }
    std::cout<<"fpA "<<(double)nblabel_A/(nsub*nsub)<<" fpB "<<(double)nblabel_B/(nsub*nsub)<<std::endl;

    std::cout<<"energy tot    "<<graphCut.compute_energy()<<std::endl;
    std::cout<<"energy data   "<<graphCut.giveDataEnergy()<<std::endl;
    std::cout<<"energy smooth "<<graphCut.giveSmoothEnergy()<<std::endl;
}
/*----------------------------------------------------------------------------*/
TEST_F(GCOTest, smallGraph_4) {

    const int nsub = 10;

    const int nbVert = nsub * nsub + 3 * nsub;
    const int nbColor = 2;

    GCoptimizationGeneralGraph graphCut(nbVert, nbColor);

    std::map<int, gmds::math::Point> pts_A;
    std::map<int, gmds::math::Point> pts_B;

    std::map<int, gmds::math::Point> pts_all;

    // init (i,j) -> index
    int indices[nsub+2][nsub+2];

    for(int i=0; i<nsub+2; i++) {
        for(int j=0; j<nsub+2; j++) {
            indices[i][j] = -1;
        }
    }
    int index = 0;
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            indices[i][j] = index;
            pts_all.emplace(index, gmds::math::Point(i,j,0.));
            index++;
        }
    }
    // left
    for(int j=1; j<=nsub; j++) {
        indices[0][j] = index;
        pts_A.emplace(index, gmds::math::Point(0.,j,0.));
        pts_all.emplace(index, gmds::math::Point(0.,j,0.));
        index++;
    }
    //right
    for(int j=1; j<=nsub; j++) {
        indices[nsub+1][j] = index;
        pts_B.emplace(index, gmds::math::Point(nsub+1,j,0.));
        pts_all.emplace(index, gmds::math::Point(nsub+1,j,0.));
        index++;
    }
    //top
    for(int i=1; i<=nsub; i++) {
        indices[i][nsub+1] = index;
        pts_A.emplace(index, gmds::math::Point(i,nsub+1,0.));
        pts_all.emplace(index, gmds::math::Point(i,nsub+1,0.));
        index++;
    }

    const double fpA = 0.7;
    const double fpB = 0.3;

    std::vector<double> dataCost(nbVert*nbColor, 1000000.);
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            const gmds::math::Point pt(i,j,0.);
            double minlength_A = HUGE_VALF;
            int index_a = -1;
            for(auto pa: pts_A) {
                if(pt.distance(pa.second) < minlength_A) {
                    minlength_A = pt.distance(pa.second);
                    index_a = pa.first;
                }
            }
            double minlength_B = HUGE_VALF;
            int index_b = -1;
            for(auto pb: pts_B) {
                if(pt.distance(pb.second) < minlength_B) {
                    minlength_B = pt.distance(pb.second);
                    index_b = pb.first;
                }
            }

            dataCost[indices[i][j] * nbColor] = (1. - fpA)*minlength_A;
            dataCost[indices[i][j] * nbColor + 1] = (1. - fpB)*minlength_B;
        }
    }
    // left A
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[0][j] * nbColor] = 0.;
    }
    //right B
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[nsub+1][j] * nbColor + 1] = 0.;
    }
    //top A
    for(int i=1; i<=nsub; i++) {
        dataCost[indices[i][nsub+1] * nbColor] = 0.;
    }


    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            const int current_index = indices[i][j];
            const int right_index = indices[i+1][j];
            const int top_index = indices[i][j+1];

            graphCut.setNeighbors(current_index, right_index);
            graphCut.setNeighbors(current_index, top_index);
        }
    }
    for(int j=1; j<=nsub; j++) {
        const int current_index = indices[0][j];
        const int right_index = indices[1][j];

        graphCut.setNeighbors(current_index, right_index);
    }

    graphCut.setDataCost(dataCost.data());

    graphCut.setSmoothCost(&gco_smoothCost_interface_, &pts_all);

    try { //graphCut.expansion(300);
        graphCut.swap();
    } catch (const GCException& e) {std::cout << e.message << '\n'; }


    std::cout<<"=============="<<std::endl;

    for(int j=nsub+1; j>=0; j--) {
        for(int i=0; i<nsub+2; i++) {

            if(indices[i][j] == -1) {
                std::cout<<" ";
            } else {
                std::cout<<graphCut.whatLabel(indices[i][j]);
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<"=============="<<std::endl;

    // cost fp-wise
    int nblabel_A = 0;
    int nblabel_B = 0;

    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            if(graphCut.whatLabel(indices[i][j]) == 0) {
                nblabel_A++;
            }
            if(graphCut.whatLabel(indices[i][j]) == 1) {
                nblabel_B++;
            }
        }
    }
    std::cout<<"fpA "<<(double)nblabel_A/(nsub*nsub)<<" fpB "<<(double)nblabel_B/(nsub*nsub)<<std::endl;

    std::cout<<"energy tot    "<<graphCut.compute_energy()<<std::endl;
    std::cout<<"energy data   "<<graphCut.giveDataEnergy()<<std::endl;
    std::cout<<"energy smooth "<<graphCut.giveSmoothEnergy()<<std::endl;
}
/*----------------------------------------------------------------------------*/
TEST_F(GCOTest, smallGraph_4_bis) {

    const int nsub = 10;

    const int nbVert = nsub * nsub + 2 * nsub;
    const int nbColor = 2;

    GCoptimizationGeneralGraph graphCut(nbVert, nbColor);

    std::map<int, gmds::math::Point> pts_A;
    std::map<int, gmds::math::Point> pts_B;

    std::map<int, gmds::math::Point> pts_all;

    // init (i,j) -> index
    int indices[nsub+2][nsub+2];

    for(int i=0; i<nsub+2; i++) {
        for(int j=0; j<nsub+2; j++) {
            indices[i][j] = -1;
        }
    }
    int index = 0;
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            indices[i][j] = index;
            pts_all.emplace(index, gmds::math::Point(i,j,0.));
            index++;
        }
    }
    // left
    for(int j=1; j<=nsub; j++) {
        indices[0][j] = index;
        pts_A.emplace(index, gmds::math::Point(0.,j,0.));
        pts_all.emplace(index, gmds::math::Point(0.,j,0.));
        index++;
    }
    //right
    for(int j=1; j<=nsub; j++) {
        indices[nsub+1][j] = index;
        pts_B.emplace(index, gmds::math::Point(nsub+1,j,0.));
        pts_all.emplace(index, gmds::math::Point(nsub+1,j,0.));
        index++;
    }

    const double fpA = 0.7;
    const double fpB = 0.3;

    std::vector<double> dataCost(nbVert*nbColor, 1000000.);
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            const gmds::math::Point pt(i,j,0.);
            double minlength_A = HUGE_VALF;
            int index_a = -1;
            for(auto pa: pts_A) {
                if(pt.distance(pa.second) < minlength_A) {
                    minlength_A = pt.distance(pa.second);
                    index_a = pa.first;
                }
            }
            double minlength_B = HUGE_VALF;
            int index_b = -1;
            for(auto pb: pts_B) {
                if(pt.distance(pb.second) < minlength_B) {
                    minlength_B = pt.distance(pb.second);
                    index_b = pb.first;
                }
            }

            dataCost[indices[i][j] * nbColor] = (1. - fpA)*minlength_A;
            dataCost[indices[i][j] * nbColor + 1] = (1. - fpB)*minlength_B;
        }
    }
    // left A
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[0][j] * nbColor] = 0.;
    }
    //right B
    for(int j=1; j<=nsub; j++) {
        dataCost[indices[nsub+1][j] * nbColor + 1] = 0.;
    }


    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<nsub; j++) {
            const int current_index = indices[i][j];
            const int right_index = indices[i+1][j];
            const int top_index = indices[i][j+1];

            graphCut.setNeighbors(current_index, right_index);
            graphCut.setNeighbors(current_index, top_index);
        }
    }
    for(int j=1; j<=nsub; j++) {
        const int current_index = indices[0][j];
        const int right_index = indices[1][j];

        graphCut.setNeighbors(current_index, right_index);
    }

    graphCut.setDataCost(dataCost.data());

    graphCut.setSmoothCost(&gco_smoothCost_interface_, &pts_all);

    try { //graphCut.expansion(300);
        graphCut.swap();
    } catch (const GCException& e) {std::cout << e.message << '\n'; }



    std::cout<<"=============="<<std::endl;

    for(int j=nsub+1; j>=0; j--) {
        for(int i=0; i<nsub+2; i++) {

            if(indices[i][j] == -1) {
                std::cout<<" ";
            } else {
                std::cout<<graphCut.whatLabel(indices[i][j]);
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<"=============="<<std::endl;

    // cost fp-wise
    int nblabel_A = 0;
    int nblabel_B = 0;

    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            if(graphCut.whatLabel(indices[i][j]) == 0) {
                nblabel_A++;
            }
            if(graphCut.whatLabel(indices[i][j]) == 1) {
                nblabel_B++;
            }
        }
    }
    std::cout<<"fpA "<<(double)nblabel_A/(nsub*nsub)<<" fpB "<<(double)nblabel_B/(nsub*nsub)<<std::endl;

    std::cout<<"energy tot    "<<graphCut.compute_energy()<<std::endl;
    std::cout<<"energy data   "<<graphCut.giveDataEnergy()<<std::endl;
    std::cout<<"energy smooth "<<graphCut.giveSmoothEnergy()<<std::endl;
}
/*----------------------------------------------------------------------------*/