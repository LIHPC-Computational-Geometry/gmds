//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/padding/SelectivePadding.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(PaddingTestClass, testSBP1)
{
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | E2R| N2R | N2F | N2E));

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);

    //==================================================================
    // MESH PREPARATION
    //==================================================================
    MeshDoctor doctor(&m);
    doctor.buildFacesAndR2F();
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();

    //====================================================================
    // Init hard constraints
    Variable<double>* quality_var = m.newVariable<double, GMDS_REGION>("GMDS_Quality");
    Variable<int>* hard = m.newVariable<int, GMDS_FACE>("GMDS_HARD");

    hard->set(0, 1); //only face 0 is hard constrained
    //==================================================================
    // SELECTIVE PADDING ALGORITHM
    //==================================================================
    SelectivePadding algo(&m);
    Variable<int>* xfi = m.newVariable<int, GMDS_FACE>("GMDS_XFI");

    algo.setHardFaces(hard);
    algo.setPaddingFaces(xfi);
    algo.execute(SelectivePadding::Option::SBP);

    int nb_padded_faces=0;
    for(auto i:m.faces()){
        if(xfi->value(i)==1){
            nb_padded_faces++;
        }
    }
    ASSERT_EQ(6, nb_padded_faces);
}
