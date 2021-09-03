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
#include <set>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/AssignCells.h"
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/ALGOCMPNT/MaterialInterfaces.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class MaterialInterfacesTest : public ::testing::Test
{
protected:
    MaterialInterfacesTest()
    {
        ;
    }
    virtual ~MaterialInterfacesTest()
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
//        int num_threads = 3;
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
TEST_F(MaterialInterfacesTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialInterfacesTest, majorityCriteria_3x3_2D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3_2D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2F, &ma, &nodesInterfaces);

    EXPECT_EQ(6, nodesInterfaces.getNbElems());

    std::set<kmds::TCellID> nodeset;
    for(int i=0; i<nodesInterfaces.getNbElems(); i++) {
        nodeset.insert(nodesInterfaces.get(i));
    }

    EXPECT_TRUE(nodeset.find(2) != nodeset.end());
    EXPECT_TRUE(nodeset.find(6) != nodeset.end());
    EXPECT_TRUE(nodeset.find(8) != nodeset.end());
    EXPECT_TRUE(nodeset.find(9) != nodeset.end());
    EXPECT_TRUE(nodeset.find(10) != nodeset.end());
    EXPECT_TRUE(nodeset.find(11) != nodeset.end());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matA("NODES_ON_INTERFACES_MATA", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matA"), &mesh, c_N2F, &ma, &nodesInterfaces_matA);

    EXPECT_EQ(5, nodesInterfaces_matA.getNbElems());

    std::set<kmds::TCellID> nodeset_matA;
    for(int i=0; i<nodesInterfaces_matA.getNbElems(); i++) {
        nodeset_matA.insert(nodesInterfaces_matA.get(i));
    }

    EXPECT_TRUE(nodeset_matA.find(2) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(6) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(8) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(9) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(10) != nodeset_matA.end());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matB("NODES_ON_INTERFACES_MATB", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matB"), &mesh, c_N2F, &ma, &nodesInterfaces_matB);

    EXPECT_EQ(4, nodesInterfaces_matB.getNbElems());

    std::set<kmds::TCellID> nodeset_matB;
    for(int i=0; i<nodesInterfaces_matB.getNbElems(); i++) {
        nodeset_matA.insert(nodesInterfaces_matB.get(i));
    }

    EXPECT_TRUE(nodeset_matA.find(2) != nodeset_matB.end());
    EXPECT_TRUE(nodeset_matA.find(6) != nodeset_matB.end());
    EXPECT_TRUE(nodeset_matA.find(10) != nodeset_matB.end());
    EXPECT_TRUE(nodeset_matA.find(11) != nodeset_matB.end());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matC("NODES_ON_INTERFACES_MATC", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matC"), &mesh, c_N2F, &ma, &nodesInterfaces_matC);

    EXPECT_EQ(4, nodesInterfaces_matC.getNbElems());

    std::set<kmds::TCellID> nodeset_matC;
    for(int i=0; i<nodesInterfaces_matC.getNbElems(); i++) {
        nodeset_matA.insert(nodesInterfaces_matC.get(i));
    }

    EXPECT_TRUE(nodeset_matA.find(8) != nodeset_matC.end());
    EXPECT_TRUE(nodeset_matA.find(9) != nodeset_matC.end());
    EXPECT_TRUE(nodeset_matA.find(10) != nodeset_matC.end());
    EXPECT_TRUE(nodeset_matA.find(11) != nodeset_matC.end());
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialInterfacesTest, majorityCriteria_2x2_2D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x2_2D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getFaceCapacity());

    elg3d::assignCellsMajorityCriteria_2D(&mesh, &fp, &ma);

    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2F();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2F, &ma, &nodesInterfaces);

    EXPECT_EQ(5, nodesInterfaces.getNbElems());

    std::set<kmds::TCellID> nodeset;
    for(int i=0; i<nodesInterfaces.getNbElems(); i++) {
        nodeset.insert(nodesInterfaces.get(i));
    }

    EXPECT_TRUE(nodeset.find(1) != nodeset.end());
    EXPECT_TRUE(nodeset.find(3) != nodeset.end());
    EXPECT_TRUE(nodeset.find(4) != nodeset.end());
    EXPECT_TRUE(nodeset.find(5) != nodeset.end());
    EXPECT_TRUE(nodeset.find(7) != nodeset.end());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matA("NODES_ON_INTERFACES_MATA", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matA"), &mesh, c_N2F, &ma, &nodesInterfaces_matA);

    EXPECT_EQ(5, nodesInterfaces_matA.getNbElems());

    std::set<kmds::TCellID> nodeset_matA;
    for(int i=0; i<nodesInterfaces_matA.getNbElems(); i++) {
        nodeset_matA.insert(nodesInterfaces_matA.get(i));
    }

    EXPECT_TRUE(nodeset_matA.find(1) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(3) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(4) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(5) != nodeset_matA.end());
    EXPECT_TRUE(nodeset_matA.find(7) != nodeset_matA.end());

    kmds::GrowingView<kmds::TCellID> nodesInterfaces_matB("NODES_ON_INTERFACES_MATB", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(ma.getMaterialID("matB"), &mesh, c_N2F, &ma, &nodesInterfaces_matB);

    EXPECT_EQ(5, nodesInterfaces_matB.getNbElems());

    std::set<kmds::TCellID> nodeset_matB;
    for(int i=0; i<nodesInterfaces_matB.getNbElems(); i++) {
        nodeset_matB.insert(nodesInterfaces_matB.get(i));
    }

    EXPECT_TRUE(nodeset_matB.find(1) != nodeset_matB.end());
    EXPECT_TRUE(nodeset_matB.find(3) != nodeset_matB.end());
    EXPECT_TRUE(nodeset_matB.find(4) != nodeset_matB.end());
    EXPECT_TRUE(nodeset_matB.find(5) != nodeset_matB.end());
    EXPECT_TRUE(nodeset_matB.find(7) != nodeset_matB.end());
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialInterfacesTest, majorityCriteria_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    EXPECT_EQ(24, nodesInterfaces.getNbElems());

    std::set<kmds::TCellID> nodeset;
    for(int i=0; i<nodesInterfaces.getNbElems(); i++) {
        nodeset.insert(nodesInterfaces.get(i));
    }

    EXPECT_TRUE(nodeset.find(2) != nodeset.end());
    EXPECT_TRUE(nodeset.find(6) != nodeset.end());
    EXPECT_TRUE(nodeset.find(10) != nodeset.end());
    EXPECT_TRUE(nodeset.find(14) != nodeset.end());

    EXPECT_TRUE(nodeset.find(18) != nodeset.end());
    EXPECT_TRUE(nodeset.find(22) != nodeset.end());
    EXPECT_TRUE(nodeset.find(26) != nodeset.end());
    EXPECT_TRUE(nodeset.find(30) != nodeset.end());

    EXPECT_TRUE(nodeset.find(32) != nodeset.end());
    EXPECT_TRUE(nodeset.find(36) != nodeset.end());
    EXPECT_TRUE(nodeset.find(40) != nodeset.end());
    EXPECT_TRUE(nodeset.find(44) != nodeset.end());

    EXPECT_TRUE(nodeset.find(33) != nodeset.end());
    EXPECT_TRUE(nodeset.find(37) != nodeset.end());
    EXPECT_TRUE(nodeset.find(41) != nodeset.end());
    EXPECT_TRUE(nodeset.find(45) != nodeset.end());

    EXPECT_TRUE(nodeset.find(34) != nodeset.end());
    EXPECT_TRUE(nodeset.find(38) != nodeset.end());
    EXPECT_TRUE(nodeset.find(42) != nodeset.end());
    EXPECT_TRUE(nodeset.find(46) != nodeset.end());

    EXPECT_TRUE(nodeset.find(35) != nodeset.end());
    EXPECT_TRUE(nodeset.find(39) != nodeset.end());
    EXPECT_TRUE(nodeset.find(43) != nodeset.end());
    EXPECT_TRUE(nodeset.find(47) != nodeset.end());
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialInterfacesTest, majorityCriteria_face_3x3x3_3D)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_3x3x3_3D(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    kmds::Connectivity* c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    EXPECT_EQ(24, nodesInterfaces.getNbElems());

    std::set<kmds::TCellID> nodeset = {2,6,10,14,18,22,26,30,32,36,40,44,33,37,41,45,34,38,42,46,35,39,43,47};
    int nbnodesfound = 0;
    for(int i=0; i<nodesInterfaces.getNbElems(); i++) {
        if(nodeset.find(nodesInterfaces.get(i)) != nodeset.end()) {
            nbnodesfound++;
        }
    }

    EXPECT_EQ(nodeset.size(), nbnodesfound);


    ch.buildFandF2R();
    kmds::Connectivity* c_F2R = mesh.getConnectivity(kmds::F2R);
//    kmds::Connectivity* c_N2F = mesh.createConnectivity(kmds::N2F);
//    ch.buildN2F();
//
//    kmds::TSize fcapacity = mesh.getFaceCapacity();
//    Kokkos::UnorderedMap<kmds::TCellID, void> interfaceFacesSet(fcapacity);
//    Kokkos::parallel_for(nodesInterfaces.getNbElems(), KOKKOS_LAMBDA(const int i) {
//        kmds::TCellID nid = nodesInterfaces.get(i);
//        Kokkos::View<kmds::TCellID *> fids;
//        c_N2F->get(nid, fids);
//        for(int i_f=0; i_f<fids.size(); i_f++) {
//
//            Kokkos::View<kmds::TCellID *> cids;
//            c_F2R->get(fids(i_f), cids);
//
//            // check if it is an interface face
//            // we ignore boundary faces
//            if(cids.size() == 2) {
//
//                if(ma.getMaterial(cids(0)) != ma.getMaterial(cids(1))) {
//                    Kokkos::UnorderedMapInsertResult res = interfaceFacesSet.insert(fids(i_f));
//                    bool success = res.success();
//                    bool exist = res.existing();
//                    bool fail = res.failed();
//                    assert(success || exist);
//                    assert(!fail);
//                }
//            }
//        }
//    });

    kmds::GrowingView<kmds::TCellID> facesInterfaces("FACES_ON_INTERFACES", mesh.getNbFaces());
    elg3d::MaterialInterfaces_getFaceOnInterfaces(&mesh, c_F2R, &ma, &facesInterfaces);

    EXPECT_EQ(facesInterfaces.getNbElems(), 15);

//    std::cout<<"nbinterfacenodes "<<nodesInterfaces.getNbElems()<<std::endl;
//    std::cout<<"facescapacity "<<fcapacity<<std::endl;
//    std::cout<<"nbfaces "<<mesh.getNbFaces()<<std::endl;
//    std::cout<<"nbinterfacesfaces "<<interfaceFacesSet.size()<<std::endl;
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialInterfacesTest, majorityCriteria_2x1x2_3D_nonmanifold)
{
    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_2x1x2_3D_nonmanifold(&mesh, &fp);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    EXPECT_EQ(10, nodesInterfaces.getNbElems());

    std::set<kmds::TCellID> nodeset;
    for (int i = 0; i < nodesInterfaces.getNbElems(); i++) {
        nodeset.insert(nodesInterfaces.get(i));
    }

    EXPECT_TRUE(nodeset.find(1) != nodeset.end());
    EXPECT_TRUE(nodeset.find(4) != nodeset.end());
    EXPECT_TRUE(nodeset.find(6) != nodeset.end());
    EXPECT_TRUE(nodeset.find(7) != nodeset.end());
    EXPECT_TRUE(nodeset.find(8) != nodeset.end());
    EXPECT_TRUE(nodeset.find(9) != nodeset.end());
    EXPECT_TRUE(nodeset.find(10) != nodeset.end());
    EXPECT_TRUE(nodeset.find(11) != nodeset.end());
    EXPECT_TRUE(nodeset.find(13) != nodeset.end());
    EXPECT_TRUE(nodeset.find(16) != nodeset.end());
}
/*----------------------------------------------------------------------------*/
TEST_F(MaterialInterfacesTest, majorityCriteria_IxJxI_3D_nonmanifold)
{
    const int nb_i = 30;
    const int nb_j = 30;

    kmds::Mesh mesh;
    elg3d::FracPres fp;
    elg3d::initData_IxJxI_3D_nonmanifold(&mesh, &fp, nb_i, nb_j);

    elg3d::MaterialAssignment ma(mesh.getRegionCapacity());

    elg3d::assignCellsMajorityCriteria_3D(&mesh, &fp, &ma);

    kmds::Connectivity *c_N2R = mesh.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh);
    ch.buildN2R();

    kmds::GrowingView<kmds::TCellID> nodesInterfaces("NODES_ON_INTERFACES", mesh.getNbNodes());
    elg3d::MaterialInterfaces_getNodeOnInterfaces(&mesh, c_N2R, &ma, &nodesInterfaces);

    EXPECT_EQ((nb_i-1)*3*(nb_j+1)+2*(nb_j+1), nodesInterfaces.getNbElems());
}
/*----------------------------------------------------------------------------*/