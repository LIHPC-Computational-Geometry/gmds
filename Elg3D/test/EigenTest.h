/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <Eigen/Eigen>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class EigenTest : public ::testing::Test
{
 protected:
    EigenTest()
        {
                ;
        }
        virtual ~EigenTest()
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
TEST_F(EigenTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(EigenTest, sparseMatProductVect)
{

    int n = 4;
    int m = 5;
    int nnz = 10;

    Eigen::SparseMatrix<double> spMat(n,m);
    spMat.reserve(nnz);

    std::vector<Eigen::Triplet<double> > tpList(8);
    tpList[0] = (Eigen::Triplet<double> (0, 0, 1.));
    tpList[1] = (Eigen::Triplet<double> (0, 4, 1.));
    tpList[2] = (Eigen::Triplet<double> (1, 0, 1.));
    tpList[3] = (Eigen::Triplet<double> (1, 4, 2.));
    tpList[4] = (Eigen::Triplet<double> (2, 1, 1.));
    tpList[5] = (Eigen::Triplet<double> (2, 2, 1.));
    tpList[6] = (Eigen::Triplet<double> (3, 2, 2.));
    tpList[7] = (Eigen::Triplet<double> (3, 3, 1.));
    spMat.setFromTriplets(tpList.begin(), tpList.end());

    EXPECT_EQ(8, spMat.nonZeros());

    std::vector<double> coefs;
    Eigen::VectorXd x(5);
    x(0) = 1.;
    x(1) = 2.;
    x(2) = 3.;
    x(3) = 4.;
    x(4) = 5.;

    Eigen::VectorXd b(4);

    b = spMat * x;

    EXPECT_DOUBLE_EQ(6., b(0));
    EXPECT_DOUBLE_EQ(11., b(1));
    EXPECT_DOUBLE_EQ(5., b(2));
    EXPECT_DOUBLE_EQ(10., b(3));
}
/*----------------------------------------------------------------------------*/