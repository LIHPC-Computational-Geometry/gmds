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
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
class FracPresTest : public ::testing::Test
{
 protected:
        FracPresTest()
        {
                ;
        }
        virtual ~FracPresTest()
        {
                ;
        }
};
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, dummy)
{
        EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, createMat)
{
        elg3d::FracPres fp;

        fp.createMaterial("iron");

        EXPECT_EQ(fp.getNbMaterials(), 1);
        EXPECT_EQ(fp.getMaterialID("iron"), 0);
        EXPECT_EQ(fp.getMaterialName(0), "iron");

        fp.createMaterial("copper");
        EXPECT_EQ(fp.getNbMaterials(), 2);
        EXPECT_EQ(fp.getMaterialID("copper"), 1);
        EXPECT_EQ(fp.getMaterialName(1), "copper");

        fp.clear();
        EXPECT_EQ(fp.getNbMaterials(), 0);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, values)
{
        elg3d::FracPres fp;

        fp.createMaterial("iron");
        fp.createMaterial("copper");

        fp.setFracPres(0, 0, 0.8);
        fp.setFracPres(0, 1, 0.7);
        EXPECT_EQ(fp.getFracPres(0,1), 0.7);
        EXPECT_EQ(fp.getMatFracPres(0).size(), 2);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, maxValues)
{
        elg3d::FracPres fp;

        fp.createMaterial("iron");
        fp.createMaterial("copper");

        fp.setFracPres(0, 0, 0.2);
        fp.setFracPres(1, 0, 0.8);
        fp.setFracPres(0, 1, 0.7);
        fp.setFracPres(1, 1, 0.1);
        EXPECT_EQ(fp.getMaxMatFracPresIndex(0), 1);
        EXPECT_EQ(fp.getMaxMatFracPresIndex(1), 0);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, getter)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);

    std::map<kmds::TCellID , double> fpmat0 = fp.getMatFracPres(0);
    std::map<kmds::TCellID , double> fpmat1 = fp.getMatFracPres(1);

    EXPECT_EQ(fpmat0.size(), 2);
    EXPECT_EQ(fpmat1.size(), 2);
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, setter)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);

    std::map<kmds::TCellID , double> fpmat0 = fp.getMatFracPres(0);
    std::map<kmds::TCellID , double> fpmat1 = fp.getMatFracPres(1);

    elg3d::FracPres fp_new;
    fp_new.setMaterialList(fp.getMaterialList());
    fp_new.setMatFracPres(0, fp.getMatFracPres(0));
    fp_new.setMatFracPres(1, fp.getMatFracPres(1));

    EXPECT_EQ(fp_new.getFracPres(0,1), fp.getFracPres(0,1));
    EXPECT_EQ(fp_new.getMatFracPres(0).size(), fp.getMatFracPres(0).size());
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, normalization_0)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);
    EXPECT_FALSE(fp.checkNormalized());

    fp.normalize(1);
    EXPECT_TRUE(fp.checkNormalized());
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, normalization_1)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 0.7);
    fp.setFracPres(1, 1, 0.1);
    EXPECT_FALSE(fp.checkNormalized());
    EXPECT_TRUE(fp.normalize());
    EXPECT_TRUE(fp.checkNormalized());
    EXPECT_FALSE(fp.normalize());
}
/*----------------------------------------------------------------------------*/
TEST_F(FracPresTest, traceAmounts)
{
    elg3d::FracPres fp;

    fp.createMaterial("iron");
    fp.createMaterial("copper");

    fp.setFracPres(0, 0, 0.2);
    fp.setFracPres(1, 0, 0.8);
    fp.setFracPres(0, 1, 10e-10);
    fp.setFracPres(1, 1, 1.);
    EXPECT_TRUE(fp.isMixedCell(1));

    fp.removeTraceAmounts();
    EXPECT_FALSE(fp.isMixedCell(1));
    EXPECT_EQ(0, fp.getFracPres(0,1));
    EXPECT_EQ(1, fp.getFracPres(1,1));
}
/*----------------------------------------------------------------------------*/