//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/padding/SelectivePadding.h>
#include <gtest/gtest.h>
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
    ASSERT_EQ(4, nb_padded_faces);
}
/*----------------------------------------------------------------------------*/
TEST(PaddingTestClass, test_cylinder1)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | E2R| N2R | N2F | N2E));

	std::string vtk_file_geom = "/home/calderans/dev/gmds/padding/cylinder_1.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F|gmds::R);
	vtkReader.read(vtk_file_geom);

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
	Variable<double>* quality_edges = m.newVariable<double, GMDS_EDGE>("GMDS_Quality_edges");
	Variable<double>* quality_faces = m.newVariable<double, GMDS_FACE>("GMDS_Quality_faces");

	Variable<int>* hard = m.newVariable<int, GMDS_FACE>("GMDS_HARD");

	for(auto e : m.edges()){
		Edge edge = m.get<Edge>(e);
		if(edge.nbRegions() == 1){
			std::vector<Face> e_faces = edge.get<Face>();
			if(e_faces.size() != 2){
				throw GMDSException("Erreur : une arête voisine à un seul hex n'a pas deux faces voisines.");
			}
			math::Vector3d norm0 = e_faces[0].normal();
			math::Vector3d norm1 = e_faces[1].normal();

			double dot = fabs(norm0.dot(norm1));

			if(dot>=0.5){
				quality_edges->set(e,dot);
			}
		}else if(edge.nbRegions() == 3){
			std::vector<Face> faces;
			for(auto const &f : edge.get<Face>()){
				if(f.nbRegions() == 2) faces.push_back(f);
			}

			math::Vector3d norm0 = faces[0].normal();
			math::Vector3d norm1 = faces[1].normal();

			double dot = norm0.dot(norm1);
			double abs_dot = fabs(norm0.dot(norm1));

			if(dot == -abs_dot){
				dot = norm0.dot(-norm1);
			}

			if(dot<=-0.5){
				quality_edges->set(e,dot);
			}
		}
	}

	for(auto f : m.faces()){
		Face face = m.get<Face>(f);
		std::vector<TCellID> edges = face.getIDs<Edge>();
		double max_quality = 0;
		for(auto e : edges){
			if(quality_edges->value(e) > max_quality) max_quality = quality_edges->value(e);
		}
		quality_faces->set(f, max_quality);
	}

	//manually hard constrained facets
	hard->set(343, 1); hard->set(377, 1); hard->set(43, 1); hard->set(2, 1);
	hard->set(344, 1); hard->set(376, 1); hard->set(44, 1); hard->set(3, 1);
	hard->set(641, 1); hard->set(667, 1); hard->set(414, 1); hard->set(381, 1);
	hard->set(642, 1); hard->set(666, 1); hard->set(415, 1); hard->set(382, 1);
	hard->set(930, 1); hard->set(955, 1); hard->set(703, 1); hard->set(670, 1);
	hard->set(931, 1); hard->set(956, 1); hard->set(704, 1); hard->set(671, 1);
	hard->set(1219, 1); hard->set(1244, 1); hard->set(992, 1); hard->set(959, 1);
	hard->set(1220, 1); hard->set(1245, 1); hard->set(993, 1); hard->set(960, 1);
	hard->set(1508, 1); hard->set(1533, 1); hard->set(1281, 1); hard->set(1248, 1);
	hard->set(1509, 1); hard->set(1534, 1); hard->set(1282, 1); hard->set(1249, 1);
	hard->set(1797, 1); hard->set(1822, 1); hard->set(1570, 1); hard->set(1537, 1);
	hard->set(1798, 1); hard->set(1823, 1); hard->set(1571, 1); hard->set(1538, 1);
	hard->set(2086, 1); hard->set(2111, 1); hard->set(1859, 1); hard->set(1826, 1);
	hard->set(2087, 1); hard->set(2112, 1); hard->set(1860, 1); hard->set(1827, 1);

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

	IGMeshIOService ios(&m);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("/home/calderans/dev/test_cylinder.vtk");

	//ASSERT_EQ(4, nb_padded_faces);
}
