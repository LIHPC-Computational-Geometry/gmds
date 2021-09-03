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