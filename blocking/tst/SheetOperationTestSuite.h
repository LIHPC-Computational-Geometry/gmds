/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>
#include "util_functions.h"
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_validate_pillow)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    gmds::cad::GeomVolume* v = geom_model.getVolumes()[0];

    auto [x_min,y_min,z_min,x_max,y_max,z_max] = v->BBox();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);
    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;
    //We pick all the boundary faces
    for(auto f:all_faces){
        if(f->info().geom_dim==2){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.validate_pillowing_surface(surf));

    surf.clear();
    //We pick a full boundary surface on coord X=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.validate_pillowing_surface(surf));
    surf.pop_back();
    ASSERT_FALSE(bl.validate_pillowing_surface(surf));
    surf.pop_back();
    ASSERT_FALSE(bl.validate_pillowing_surface(surf));
    surf.pop_back();
    ASSERT_FALSE(bl.validate_pillowing_surface(surf));
}

/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    ASSERT_EQ(bl.get_nb_cells<3>(),1);
    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(1);
    //export_vtk(bl,gmds::N | gmds::F, "pillow_1_surf.vtk");
    ASSERT_EQ(bl.get_nb_cells<3>(),2);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 4);
    ASSERT_EQ(nb_nodes_on_surface, 0);
    ASSERT_EQ(nb_nodes_in_volume, 0);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 16);
    ASSERT_EQ(nb_edges_on_surface, 4);
    ASSERT_EQ(nb_edges_in_volume, 0);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 10);
    ASSERT_EQ(nb_faces_in_volume, 1);
}
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_2)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }
    ASSERT_EQ(bl.get_nb_cells<3>(),8);
    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;
    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_2_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),16);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 16);
    ASSERT_EQ(nb_nodes_on_surface, 14);
    ASSERT_EQ(nb_nodes_in_volume, 4);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 28);
    ASSERT_EQ(nb_edges_on_surface, 44);
    ASSERT_EQ(nb_edges_in_volume, 19);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 36);
    ASSERT_EQ(nb_faces_in_volume, 30);
}

/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_3)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_3_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),3);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 4);
    ASSERT_EQ(nb_nodes_on_surface, 2);
    ASSERT_EQ(nb_nodes_in_volume, 0);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 16);
    ASSERT_EQ(nb_edges_on_surface, 8);
    ASSERT_EQ(nb_edges_in_volume, 1);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 12);
    ASSERT_EQ(nb_faces_in_volume, 3);}

/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_4)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);
    //export_vtk(bl,gmds::N | gmds::F, "pillow_4_surf.vtk");
    ASSERT_EQ(bl.get_nb_cells<3>(),4);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 3);
    ASSERT_EQ(nb_nodes_on_surface, 3);
    ASSERT_EQ(nb_nodes_in_volume, 1);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 15);
    ASSERT_EQ(nb_edges_on_surface, 9);
    ASSERT_EQ(nb_edges_in_volume, 4);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 12);
    ASSERT_EQ(nb_faces_in_volume, 6);
}

/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_5)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()+5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);
    //export_vtk(bl,gmds::N | gmds::F, "pillow_5_surf.vtk");
    ASSERT_EQ(bl.get_nb_cells<3>(),5);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 2);
    ASSERT_EQ(nb_nodes_on_surface, 4);
    ASSERT_EQ(nb_nodes_in_volume, 2);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 14);
    ASSERT_EQ(nb_edges_on_surface, 10);
    ASSERT_EQ(nb_edges_in_volume, 7);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 12);
    ASSERT_EQ(nb_faces_in_volume, 9);
}
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_6)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()+5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()+5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_6_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),6);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 0);
    ASSERT_EQ(nb_nodes_on_surface, 4);
    ASSERT_EQ(nb_nodes_in_volume, 4);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 12);
    ASSERT_EQ(nb_edges_on_surface, 8);
    ASSERT_EQ(nb_edges_in_volume, 12);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 10);
    ASSERT_EQ(nb_faces_in_volume, 13);
}
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_7)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();

    ASSERT_TRUE(bl.pillow(all_faces));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_7_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),7);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 0);
    ASSERT_EQ(nb_nodes_on_surface, 0);
    ASSERT_EQ(nb_nodes_in_volume, 8);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 12);
    ASSERT_EQ(nb_edges_on_surface, 0);
    ASSERT_EQ(nb_edges_in_volume, 20);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 6);
    ASSERT_EQ(nb_faces_in_volume, 18);
}
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_8)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0
    for(auto f:all_faces){
        auto nb_adj = bl.get_blocks_of_face(f);
        if(nb_adj.size()==1)
            surf.push_back(f);
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_8_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),32);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 12);
    ASSERT_EQ(nb_nodes_on_surface, 6);
    ASSERT_EQ(nb_nodes_in_volume, 27);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 24);
    ASSERT_EQ(nb_edges_on_surface, 24);
    ASSERT_EQ(nb_edges_in_volume, 80);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 24);
    ASSERT_EQ(nb_faces_in_volume, 84);
}
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_9)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick a full boundary surface on coord X=5.0 and Y=5.0
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Y()-5)<0.1){
            surf.push_back(f);
        }
        else if(fabs(ci.Z()-5)<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_9_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),20);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 15);
    ASSERT_EQ(nb_nodes_on_surface, 15);
    ASSERT_EQ(nb_nodes_in_volume, 8);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 27);
    ASSERT_EQ(nb_edges_on_surface, 45);
    ASSERT_EQ(nb_edges_in_volume, 31);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 36);
    ASSERT_EQ(nb_faces_in_volume, 42);
}
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_10)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick an inner surface
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X())<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_EQ(surf.size(),4);

    //export_vtk(bl,gmds::N | gmds::F, "pillow_10_surf_before.vtk");
    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::E, "pillow_10_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),12);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 16);
    ASSERT_EQ(nb_nodes_on_surface, 10);
    ASSERT_EQ(nb_nodes_in_volume, 2);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 28);
    ASSERT_EQ(nb_edges_on_surface, 36);
    ASSERT_EQ(nb_edges_in_volume, 11);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 32);
    ASSERT_EQ(nb_faces_in_volume, 20);
}

/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_11)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick an inner surface
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X())<0.1 && ci.Y()<0 && ci.Z()<0){
            surf.push_back(f);
        }
        else if(ci.X() <0 && fabs(ci.Y())<0.1 && ci.Z()<0){
            surf.push_back(f);
        }
        else if(ci.X()<0 && ci.Y()<0 && fabs(ci.Z())<0.1){
            surf.push_back(f);
        }
    }
    ASSERT_EQ(surf.size(),3);

    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

	 //export_vtk(bl,gmds::N | gmds::E, "pillow_11_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),11);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 15);
    ASSERT_EQ(nb_nodes_on_surface, 9);
    ASSERT_EQ(nb_nodes_in_volume, 2);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 27);
    ASSERT_EQ(nb_edges_on_surface, 33);
    ASSERT_EQ(nb_edges_in_volume, 10);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 30);
    ASSERT_EQ(nb_faces_in_volume, 18);
}

/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_pillow_12)
{
    gmds::cad::FACManager geom_model;
    set_up(&geom_model,"/tet_in_box.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model, true);
    gmds::blocking::CurvedBlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gmds::blocking::CurvedBlocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for(auto par_edges_i:edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    ASSERT_EQ(bl.get_nb_cells<3>(),8);

    std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
    std::vector<gmds::blocking::CurvedBlocking::Face> surf;

    surf.clear();
    //We pick an inner surface
    for(auto f:all_faces){
        gmds::math::Point ci = bl.get_center_of_face(f);
        if(fabs(ci.X())<0.1 ){
            surf.push_back(f);
        }
        else if(ci.X() <0 && fabs(ci.Y()-5)<0.1 ){
            surf.push_back(f);
        }
    }
    ASSERT_EQ(surf.size(),6);

    ASSERT_TRUE(bl.pillow(surf));

    bl.smooth(10);

    //export_vtk(bl,gmds::N | gmds::E, "pillow_12_surf.vtk");

    ASSERT_EQ(bl.get_nb_cells<3>(),14);
    auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
    ASSERT_EQ(nb_nodes_on_vertex, 8);
    ASSERT_EQ(nb_nodes_on_curve, 16);
    ASSERT_EQ(nb_nodes_on_surface, 12);
    ASSERT_EQ(nb_nodes_in_volume, 3);
    auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
    ASSERT_EQ(nb_edges_on_curve, 28);
    ASSERT_EQ(nb_edges_on_surface, 40);
    ASSERT_EQ(nb_edges_in_volume, 15);
    auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
    ASSERT_EQ(nb_faces_on_surface, 34);
    ASSERT_EQ(nb_faces_in_volume, 25);
}
/*----------------------------------------------------------------------------*/
TEST(SheetOperationTestSuite, test_notch)
{
	 gmds::cad::FACManager geom_model;
	 set_up(&geom_model,"/Notch/notch_tet.vtk");
	 gmds::blocking::CurvedBlocking bl(&geom_model, false);

	 gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

	 gmds::IGMeshIOService ioService(&m);
	 gmds::VTKReader vtkReader(&ioService);
	 vtkReader.setCellOptions(gmds::N | gmds::R);
	 std::string dir(TEST_SAMPLES_DIR);
	 vtkReader.read(dir+"/Notch/notch_blocks.vtk");

	 bl.init_from_mesh(m);
	 gmds::blocking::CurvedBlockingClassifier cl(&bl);
	 cl.classify();

	 //the classification does not consider volume right now. We enforce unclassified cells to be on volume 1 so.
	 std::vector<gmds::blocking::CurvedBlocking::Edge> all_edges = bl.get_all_edges();
	 std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = bl.get_all_faces();
	 std::vector<gmds::blocking::CurvedBlocking::Block> all_blocks = bl.get_all_blocks();
	 for(auto &c:all_edges){
		  if(c->info().geom_dim==4){
			   c->info().geom_dim=3;
			   c->info().geom_id=1;
		  }
	 }
	 for(auto &c:all_faces){
		  if(c->info().geom_dim==4){
			   c->info().geom_dim=3;
			   c->info().geom_id=1;
		  }
	 }	 for(auto &c:all_blocks){
		  if(c->info().geom_dim==4){
			   c->info().geom_dim=3;
			   c->info().geom_id=1;
		  }
	 }
	 ASSERT_EQ(bl.get_nb_cells<3>(),7);
	 auto [nb_nodes_on_vertex,nb_nodes_on_curve,nb_nodes_on_surface,nb_nodes_in_volume] = get_node_statistics(bl);
	 ASSERT_EQ(nb_nodes_on_vertex, 12);
	 ASSERT_EQ(nb_nodes_on_curve, 11);
	 ASSERT_EQ(nb_nodes_on_surface, 3);
	 ASSERT_EQ(nb_nodes_in_volume, 0);
	 auto [nb_edges_on_curve,nb_edges_on_surface,nb_edges_in_volume] = get_edge_statistics(bl);
	 ASSERT_EQ(nb_edges_on_curve, 29);
	 ASSERT_EQ(nb_edges_on_surface, 19);
	 ASSERT_EQ(nb_edges_in_volume, 3);
	 auto [nb_faces_on_surface,nb_faces_in_volume] = get_face_statistics(bl);
	 ASSERT_EQ(nb_faces_on_surface, 24);
	 ASSERT_EQ(nb_faces_in_volume, 9);

	 std::vector<gmds::blocking::CurvedBlocking::Face> surf;
	 //We pick the set of faces we want to use to pillow
	 gmds::math::Point ref[4] = {gmds::math::Point(3.29331, 3, 1.70669),
		                          gmds::math::Point(2.29331, 3, 0.706695),
		                          gmds::math::Point(3.19709, 1, 1.79022),
		                          gmds::math::Point(2.24882, 1, 0.824707)};
	 for(auto f:all_faces){
		  gmds::math::Point ci = bl.get_center_of_face(f);
		  if(ci.distance(ref[0])<0.1 ||
		      ci.distance(ref[1])<0.1 ||
		      ci.distance(ref[2])<0.1 ||
		      ci.distance(ref[3])<0.1 ){
			   surf.push_back(f);
		  }
	 }
	 ASSERT_EQ(surf.size(),4);
	 ASSERT_TRUE(bl.pillow(surf));
	 bl.smooth(10);

	 auto [nb_nodes_on_vertex2,nb_nodes_on_curve2,nb_nodes_on_surface2,nb_nodes_in_volume2] = get_node_statistics(bl);
	 ASSERT_EQ(nb_nodes_on_vertex2, 12);
	 ASSERT_EQ(nb_nodes_on_curve2, 15);
	 ASSERT_EQ(nb_nodes_on_surface2, 7);
	 ASSERT_EQ(nb_nodes_in_volume2, 1);
	 auto [nb_edges_on_curve2,nb_edges_on_surface2,nb_edges_in_volume2] = get_edge_statistics(bl);
	 ASSERT_EQ(nb_edges_on_curve2, 33);
	 ASSERT_EQ(nb_edges_on_surface2, 31);
	 ASSERT_EQ(nb_edges_in_volume2, 8);
	 auto [nb_faces_on_surface2,nb_faces_in_volume2] = get_face_statistics(bl);
	 ASSERT_EQ(nb_faces_on_surface2, 32);
	 ASSERT_EQ(nb_faces_in_volume2, 17);
	 //export_vtk(bl,gmds::N | gmds::E, "pillow_notch.vtk");
}
