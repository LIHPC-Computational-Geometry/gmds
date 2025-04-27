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
#include <KM/DS/Mesh.h>

#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InitData.h"
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class ToolsTest : public ::testing::Test
{
 protected:
        ToolsTest()
        {
                ;
        }
        virtual ~ToolsTest()
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
TEST_F(ToolsTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(ToolsTest, computeMidPoint_2D)
{
    kmds::Mesh mesh;

    const double xyz_min[3] = {0,0,0};
    const double xyz_max[3] = {1,1,1};
    const int ni = 7;
    const int nj = 3;

    kmds::InitTools_createGrid_2D(&mesh, xyz_min, xyz_max, ni, nj);

    kmds::Variable<gmds::math::Point>* varMidpoints =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_FACE, "varMidPoint");

    const kmds::TSize nbCells = mesh.getNbFaces();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getFaceIDs(&cellIDs);

    elg3d::Tools_computeMidPointVar_2D(&cellIDs, &mesh, varMidpoints);

    EXPECT_EQ((*varMidpoints)[cellIDs.get(0)], mesh.getFace(cellIDs.get(0)).midpoint());

    bool foundInvalid = false;
    for(int i=0; i<nbCells; i++) {
        foundInvalid = (gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF) == (*varMidpoints)[cellIDs.get(i)]);
        if(foundInvalid) {
            break;
        }
    }
    EXPECT_FALSE(foundInvalid);

    // partially compute the midpoints
    kmds::Variable<gmds::math::Point>* varMidpoints_trunc =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_FACE, "varMidPoint_trunc");
    kmds::GrowingView <kmds::TCellID> cellIDs_trunc("CELLS", nbCells);
    cellIDs_trunc.push_back(1);
    cellIDs_trunc.push_back(7);
    cellIDs_trunc.push_back(17);

    elg3d::Tools_computeMidPointVar_2D(&cellIDs_trunc, &mesh, varMidpoints_trunc);
    EXPECT_EQ((*varMidpoints_trunc)[cellIDs_trunc.get(0)], mesh.getFace(cellIDs_trunc.get(0)).midpoint());

    int nbInvalid = 0;
    for(int i=0; i<nbCells; i++) {
        if(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF) == (*varMidpoints_trunc)[cellIDs_trunc.get(i)]) {
            nbInvalid++;
        }
    }
    EXPECT_EQ(nbInvalid, nbCells - 3);
}
/*----------------------------------------------------------------------------*/
TEST_F(ToolsTest, computeMidPoint_3D)
{
    kmds::Mesh mesh;

    const double xyz_min[3] = {0,0,0};
    const double xyz_max[3] = {1,1,1};
    const int ni = 7;
    const int nj = 3;
    const int nk = 3;

    kmds::InitTools_createGrid_3D(&mesh, xyz_min, xyz_max, ni, nj, nk);

    kmds::Variable<gmds::math::Point>* varMidpoints =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "varMidPoint");

    const kmds::TSize nbCells = mesh.getNbRegions();
    kmds::GrowingView <kmds::TCellID> cellIDs("CELLS", nbCells);
    mesh.getRegionIDs(&cellIDs);

    elg3d::Tools_computeMidPointVar_3D(&cellIDs, &mesh, varMidpoints);

    EXPECT_EQ((*varMidpoints)[cellIDs.get(0)], mesh.getRegion(cellIDs.get(0)).midpoint());

    bool foundInvalid = false;
    for(int i=0; i<nbCells; i++) {
        foundInvalid = (gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF) == (*varMidpoints)[cellIDs.get(i)]);
        if(foundInvalid) {
            break;
        }
    }
    EXPECT_FALSE(foundInvalid);

    // partially compute the midpoints
    kmds::Variable<gmds::math::Point>* varMidpoints_trunc =
            mesh.createVariable<gmds::math::Point>(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF), kmds::KMDS_REGION, "varMidPoint_trunc");
    kmds::GrowingView <kmds::TCellID> cellIDs_trunc("CELLS", nbCells);
    cellIDs_trunc.push_back(1);
    cellIDs_trunc.push_back(7);
    cellIDs_trunc.push_back(17);

    elg3d::Tools_computeMidPointVar_3D(&cellIDs_trunc, &mesh, varMidpoints_trunc);
    EXPECT_EQ((*varMidpoints_trunc)[cellIDs_trunc.get(0)], mesh.getRegion(cellIDs_trunc.get(0)).midpoint());

    int nbInvalid = 0;
    for(int i=0; i<nbCells; i++) {
        if(gmds::math::Point(-HUGE_VALF, -HUGE_VALF, -HUGE_VALF) == (*varMidpoints_trunc)[cellIDs_trunc.get(i)]) {
            nbInvalid++;
        }
    }
    EXPECT_EQ(nbInvalid, nbCells - 3);
}
/*----------------------------------------------------------------------------*/
TEST_F(ToolsTest, loadstoreNodePos)
{
    kmds::Mesh mesh;


}
/*----------------------------------------------------------------------------*/