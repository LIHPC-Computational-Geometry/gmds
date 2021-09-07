/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/25/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/BoundaryExtractor3D.h>
#include <gmds/igalgo/BoundaryExtractor2D.h>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(BoundaryExClass, testCubeBoundary)
{
    // WE WRITE
    gmds::Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                                     R2N|R2F|R2E|
                                     F2N|F2R|F2E|
                                     E2F|E2N|N2E));
    gmds::Mesh m_surf(gmds::MeshModel(DIM3|F|E|N|F2N|F2E|E2F|N2E|E2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/tet_in_box.vtk";
    gmds::IGMeshIOService ioService(&m_vol);

    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    gmds::BoundaryExtractor3D boundary_extractor(&m_vol,&m_surf);
    ASSERT_TRUE(boundary_extractor.isValid());


    Variable<int>* node_on_pnt = m_surf.newVariable<int, GMDS_NODE>("on_point"  );
    Variable<int>* node_on_crv = m_surf.newVariable<int, GMDS_NODE>("on_curve"  );
    Variable<int>* node_on_srf = m_surf.newVariable<int, GMDS_NODE>("on_surface");
    Variable<int>* edge_on_crv = m_surf.newVariable<int, GMDS_EDGE>("on_curve"  );
    Variable<int>* face_on_srf = m_surf.newVariable<int, GMDS_FACE>("on_surface");
   boundary_extractor.setColorOption(node_on_pnt, edge_on_crv, face_on_srf, node_on_srf,node_on_crv);

    boundary_extractor.execute();

    std::map<int, int> nb_faces_per_surface;
    for(auto f_id:m_surf.faces()){
        nb_faces_per_surface[face_on_srf->value(f_id)]++;
    }
    ASSERT_EQ(nb_faces_per_surface.size(),6);
    ASSERT_EQ(nb_faces_per_surface[1],8);
    ASSERT_EQ(nb_faces_per_surface[2],8);
    ASSERT_EQ(nb_faces_per_surface[3],8);
    ASSERT_EQ(nb_faces_per_surface[4],8);
    ASSERT_EQ(nb_faces_per_surface[5],8);
    ASSERT_EQ(nb_faces_per_surface[6],8);


}
/*----------------------------------------------------------------------------*/
TEST(BoundaryExClass, testHalfDisk)
{
    // WE WRITE
    gmds::Mesh m_surf(gmds::MeshModel(DIM3|F|E|N|F2N|N2F|E2N));
    gmds::Mesh m_bnd (gmds::MeshModel(DIM3|E|N|N2E|E2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/half_disk.vtk";
    gmds::IGMeshIOService ioService(&m_surf);

    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m_surf);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();

    gmds::BoundaryExtractor2D boundary_extractor(&m_surf,&m_bnd);
    ASSERT_TRUE(boundary_extractor.isValid());


    Variable<int>* node_on_pnt = m_bnd.newVariable<int, GMDS_NODE>("on_point"  );
    Variable<int>* node_on_crv = m_bnd.newVariable<int, GMDS_NODE>("on_curve"  );
    Variable<int>* edge_on_crv = m_bnd.newVariable<int, GMDS_EDGE>("on_curve"  );
    boundary_extractor.setColorOption(node_on_pnt, edge_on_crv, node_on_crv);

    boundary_extractor.execute();

    std::set<int> crv_ids, pnt_ids;
    for(auto n_id:m_bnd.nodes()){
        crv_ids.insert(node_on_crv->value(n_id));
        pnt_ids.insert(node_on_pnt->value(n_id));
    }
    ASSERT_EQ(pnt_ids.size(),3);
    ASSERT_EQ(crv_ids.size(),3);


}/*----------------------------------------------------------------------------*/
TEST(BoundaryExClass, testHolesInSquare)
{
    // WE WRITE
    gmds::Mesh m_surf(gmds::MeshModel(DIM3|F|E|N|F2N|N2F|E2N));
    gmds::Mesh m_bnd (gmds::MeshModel(DIM3|E|N|N2E|E2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/HolesInSquare0.vtk";
    gmds::IGMeshIOService ioService(&m_surf);

    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m_surf);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();

    gmds::BoundaryExtractor2D boundary_extractor(&m_surf,&m_bnd);
    ASSERT_TRUE(boundary_extractor.isValid());


    Variable<int>* node_on_pnt = m_bnd.newVariable<int, GMDS_NODE>("on_point"  );
    Variable<int>* node_on_crv = m_bnd.newVariable<int, GMDS_NODE>("on_curve"  );
    Variable<int>* edge_on_crv = m_bnd.newVariable<int, GMDS_EDGE>("on_curve"  );
    boundary_extractor.setColorOption(node_on_pnt, edge_on_crv, node_on_crv);

    boundary_extractor.execute();

    std::set<int> crv_ids, pnt_ids;
    for(auto n_id:m_bnd.nodes()){
        crv_ids.insert(node_on_crv->value(n_id));
        pnt_ids.insert(node_on_pnt->value(n_id));
    }
    ASSERT_EQ(pnt_ids.size(),6);
    ASSERT_EQ(crv_ids.size(),7);


}