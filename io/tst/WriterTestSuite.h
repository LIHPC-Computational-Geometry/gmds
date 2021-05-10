//
// Created by ledouxf on 1/22/19.
//

/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/math/Point.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;

/*----------------------------------------------------------------------------*/
TEST(MeshClass, testWriterVTK_NF)
{
    // WE WRITE
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));

    gmds::math::Vector3d normal(0,0,1);
    gmds::Variable<int>* v_n = m.newVariable<int,gmds::GMDS_NODE>("val_int");
    gmds::Variable<double>* v_f = m.newVariable<double,gmds::GMDS_FACE>("val_double");


    gmds::Variable<gmds::math::Vector3d>* v_nv = m.newVariable<gmds::math::Vector3d,gmds::GMDS_NODE>("val_vec");
    gmds::Variable<gmds::math::Vector3d>* v_fv = m.newVariable<gmds::math::Vector3d,gmds::GMDS_FACE>("val_vec");

    gmds::Node n0 = m.newNode(0,0,0);
    gmds::Node n1 = m.newNode(1,0,0);
    gmds::Node n2 = m.newNode(1,1,0);
    gmds::Node n3 = m.newNode(0,1,0);
    gmds::Node n4 = m.newNode(2,0,0);

    gmds::Face f0 = m.newTriangle(n1,n2,n4);

    gmds::Face f1 = m.newQuad(n0,n1,n2,n3);

    v_f->set(f0.id(),2.5);
    v_f->set(f1.id(),4.2);

    v_fv->set(f1.id(),normal);
    for(auto i:m.nodes()){
        v_n->set(i,10*i);
        v_nv->set(i,normal);
    }

    gmds::CellGroup<gmds::Face>* cg = m.newGroup<gmds::Face>("surface1");
    cg->add(f0);

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("toto");

    // WE READ
    gmds::Mesh m2(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));
    gmds::IGMeshIOService ioService2(&m2);
    gmds::VTKReader vtkReader(&ioService2);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.setDataOptions(gmds::N|gmds::F);
    vtkReader.read("toto");

    ASSERT_EQ(m2.getNbNodes(),5);
    ASSERT_EQ(m2.getNbFaces(),2);


    gmds::Variable<int>* v2_n = m2.getVariable<int,gmds::GMDS_NODE>("val_int");

    gmds::Variable<gmds::math::Vector3d>* v2_nv = m.getVariable<gmds::math::Vector3d,gmds::GMDS_NODE>("val_vec");
    gmds::Variable<gmds::math::Vector3d>* v2_fv = m.getVariable<gmds::math::Vector3d,gmds::GMDS_FACE>("val_vec");


    ASSERT_EQ((*v2_n)[n1.id()],10);
    ASSERT_EQ((*v2_n)[n3.id()],30);

    for(auto i:m2.nodes()){
        ASSERT_EQ((*v2_nv)[i],normal);
    }
    ASSERT_EQ(m2.getAllVariables(gmds::GMDS_NODE).size(),3);

    ASSERT_EQ((*v2_fv)[0],gmds::math::Vector3d(0,0,0));
    ASSERT_EQ((*v2_fv)[1],normal);

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testWriterVTK_NR)
{
    // WE WRITE
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::N|gmds::R2N));

    gmds::math::Vector3d normal(0,0,1);
    gmds::Variable<int>* v_n = m.newVariable<int,gmds::GMDS_NODE>("val_int");
    gmds::Variable<double>* v_f = m.newVariable<double,gmds::GMDS_REGION>("val_double");


    gmds::Variable<gmds::math::Vector3d>* v_nv = m.newVariable<gmds::math::Vector3d,gmds::GMDS_NODE>("val_vec");
    gmds::Variable<gmds::math::Vector3d>* v_fv = m.newVariable<gmds::math::Vector3d,gmds::GMDS_REGION>("val_vec");
    int nb_nodes = 5;
    gmds::TCellID n[nb_nodes][nb_nodes][nb_nodes];
    for(auto i=0;i<nb_nodes;i++) {
        for (auto j = 0; j < nb_nodes; j++) {
            for (auto k = 0; k < nb_nodes; k++) {
                n[i][j][k] = m.newNode(i, j, k).id();
                v_n->set(n[i][j][k], 10 * n[i][j][k]);
                v_nv->set(n[i][j][k], normal);
            }
        }
    }
    for(auto i=0;i<nb_nodes-1;i++) {
        for (auto j = 0; j < nb_nodes-1; j++) {
            for (auto k = 0; k < nb_nodes-1; k++) {
                gmds::Region rijk= m.newHex(n[i][j][k],n[i+1][j][k],
                                      n[i+1][j+1][k],n[i][j+1][k],
                                      n[i][j][k+1],n[i+1][j][k+1],
                                      n[i+1][j+1][k+1],n[i][j+1][k+1]);
                v_fv->set(rijk.id(), normal);

            }
        }
    }


    gmds::IGMeshIOService ioService(&m);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::R);
    vtkWriter.write("toto_3d");

    // WE READ
    gmds::Mesh m2(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::N|gmds::R2N));
    gmds::IGMeshIOService ioService2(&m2);
    gmds::VTKReader vtkReader(&ioService2);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read("toto_3d");

    ASSERT_EQ(m2.getNbNodes(),125);
    ASSERT_EQ(m2.getNbRegions(),64);


    gmds::Variable<int>* v2_n = m2.getVariable<int,gmds::GMDS_NODE>("val_int");

    gmds::Variable<gmds::math::Vector3d>* v2_nv = m2.getVariable<gmds::math::Vector3d,gmds::GMDS_NODE>("val_vec");
    gmds::Variable<gmds::math::Vector3d>* v2_fv = m2.getVariable<gmds::math::Vector3d,gmds::GMDS_REGION>("val_vec");


    ASSERT_EQ((*v2_n)[1],10);
    ASSERT_EQ((*v2_n)[2],20);

    for(auto i:m2.nodes()){
        ASSERT_EQ((*v2_nv)[i],normal);
    }
    ASSERT_EQ(m2.getAllVariables(gmds::GMDS_NODE).size(),3);

    ASSERT_EQ((*v2_fv)[0],gmds::math::Vector3d(0,0,1));
    ASSERT_EQ((*v2_fv)[1],normal);

}
