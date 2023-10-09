//
// Created by rochec on 26/09/23.
//

#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/claire/MFEMMeshWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(MFEMMeshWriterTestClass, MFEM_Test)
{
	// Mesh
	Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));
	std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir+"/Aero/2D/APOLLO_2D_toFit.vtk";
	//std::string vtk_file = dir+"/Aero/2D/APOLLO_2D_20k_MESH.vtk";
	std::string vtk_file = dir+"/Aero/2D/Test2.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	//MFEMMeshWriter MeshWriter(&m, "APOLLO_2D_toFit");
	//MFEMMeshWriter MeshWriter(&m, "APOLLO_2D_20k_MESH");
	MFEMMeshWriter MeshWriter(&m, "Test2");
	MFEMMeshWriter::STATUS res = MeshWriter.execute();

	ASSERT_EQ(res, MeshWriter.SUCCESS);

}


TEST(MFEMMeshWriterTestClass, MFEM_BG_and_LS)
{
	Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Aero/2D/Stardust_2D_inCircle_Int.vtk";

	std::string output_file = "Stardust_2D" ;

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.setDataOptions(gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	Variable<int> *var_material = m.getVariable<int, GMDS_FACE>("CellEntityIds");
	Variable<double> *var_dist = m.newVariable<double, GMDS_NODE>("GMDS_Distance");
	std::vector<TCellID> bnd;
	std::vector<TCellID> bnd_edges;
	TInt markTreatedNodes = m.newMark<gmds::Node>();

	for (auto n_id:m.nodes())
	{
		Node n = m.get<Node>(n_id);
		//std::cout << n.get<Face>().size() << std::endl;
		int compteur_mat1(0);
		int compteur_mat2(0);
		for (auto const f:n.get<Face>())
		{
			//std::cout << "Face 1, mat " << var_material->value(f.id()) << std::endl;
			if (var_material->value(f.id()) == 1)
			{
				compteur_mat1++;
			}
			else if (var_material->value(f.id())==-1)
			{
				compteur_mat2++;
			}
		}
		//std::cout << "Hello" << std::endl;
		if (compteur_mat1 != 0 && compteur_mat2 != 0)
		{
			bnd.push_back(n_id);
			var_dist->set(n_id, 0);
			m.mark(n, markTreatedNodes);
		}
	}

	// Compute boundary edges
	for (auto e_id:m.edges())
	{
		Edge e = m.get<Edge>(e_id);
		std::vector<Face> e_faces = e.get<Face>();
		if ( e_faces.size() == 2
		    && (var_material->value(e_faces[0].id()) != var_material->value(e_faces[1].id())) )
		{
			bnd_edges.push_back(e.id());
		}
	}

	// Compute dist for other nodes
	for (auto n_id:m.nodes())
	{
		Node n = m.get<Node>(n_id);
		if (!m.isMarked(n, markTreatedNodes))
		{
			double min_dist = std::numeric_limits<double>::max() ;
			// Check the dist to the edges
			for (auto e_id:bnd_edges)
			{
				math::Segment e_seg = m.get<Edge>(e_id).segment();
				math::Point p_proj = n.point();
				p_proj = e_seg.project(p_proj);
				if ((n.point()-p_proj).norm() < min_dist)
				{
					min_dist = (n.point()-p_proj).norm() ;
				}
			}
			var_dist->set(n_id, var_material->value(n.get<Face>()[0].id())*min_dist);
			m.mark(n, markTreatedNodes);
		}
	}


	m.unmarkAll<Node>(markTreatedNodes);
	m.freeMark<Node>(markTreatedNodes);


	// .mesh writing of background mesh for MFEM
	MFEMMeshWriter mfemWriter = MFEMMeshWriter(&m, output_file+"_BG");
	MFEMMeshWriter::STATUS res_writing = mfemWriter.execute();
	ASSERT_EQ(res_writing, MFEMMeshWriter::SUCCESS);


	// .gf writing of the LS for MFEM
	{
		std::ofstream stream = std::ofstream(output_file+"_LS.gf", std::ios::out);
		stream.precision(15);

		// Header
		stream << "FiniteElementSpace\n";
		stream << "FiniteElementCollection: H1_2D_P1\n";
		stream << "VDim: 1\n";
		stream << "Ordering: 0\n";
		stream << "\n";

		for (auto n : m.nodes()) {
			stream << var_dist->value(n) << "\n";
		}
	}

}

TEST(MFEMMeshWriterTestClass, MFEM_BG_and_LS_3D)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Aero/3D/Apollo_3D_inSphere_Int.vtk";

	std::string output_file = "Apollo_3D" ;

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.setDataOptions(gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	Variable<int> *var_material = m.getVariable<int, GMDS_REGION>("CellEntityIds");
	Variable<double> *var_dist = m.newVariable<double, GMDS_NODE>("GMDS_Distance");
	std::vector<TCellID> bnd;
	std::vector<TCellID> bnd_faces;
	TInt markTreatedNodes = m.newMark<gmds::Node>();

	for (auto n_id:m.nodes())
	{
		Node n = m.get<Node>(n_id);
		//std::cout << n.get<Face>().size() << std::endl;
		int compteur_mat1(0);
		int compteur_mat2(0);
		for (auto const f:n.get<Face>())
		{
			//std::cout << "Face 1, mat " << var_material->value(f.id()) << std::endl;
			if (var_material->value(f.id()) == 1)
			{
				compteur_mat1++;
			}
			else if (var_material->value(f.id())==-1)
			{
				compteur_mat2++;
			}
		}
		//std::cout << "Hello" << std::endl;
		if (compteur_mat1 != 0 && compteur_mat2 != 0)
		{
			bnd.push_back(n_id);
			var_dist->set(n_id, 0);
			m.mark(n, markTreatedNodes);
		}
	}

	// Compute boundary faces
	for (auto f_id:m.faces())
	{
		Face f = m.get<Face>(f_id);
		std::vector<Region> f_regions = f.get<Region>();
		if ( f_regions.size() == 2
		    && (var_material->value(f_regions[0].id()) != var_material->value(f_regions[1].id())) )
		{
			bnd_faces.push_back(f.id());
		}
	}

	// Compute dist for other nodes
	for (auto n_id:m.nodes())
	{
		Node n = m.get<Node>(n_id);
		if (!m.isMarked(n, markTreatedNodes))
		{
			double min_dist = std::numeric_limits<double>::max() ;
			// Check the dist to the edges
			for (auto f_id:bnd_faces)
			{
				Face f = m.get<Face>(f_id);
				math::Point p_proj = n.point();
				p_proj = f.project(p_proj);
				if ((n.point()-p_proj).norm() < min_dist)
				{
					min_dist = (n.point()-p_proj).norm() ;
				}
			}
			var_dist->set(n_id, var_material->value(n.get<Region>()[0].id())*min_dist);
			m.mark(n, markTreatedNodes);
		}
	}


	m.unmarkAll<Node>(markTreatedNodes);
	m.freeMark<Node>(markTreatedNodes);


	// .mesh writing of background mesh for MFEM
	MFEMMeshWriter mfemWriter = MFEMMeshWriter(&m, output_file+"_BG");
	MFEMMeshWriter::STATUS res_writing = mfemWriter.execute();
	ASSERT_EQ(res_writing, MFEMMeshWriter::SUCCESS);


	// .gf writing of the LS for MFEM
	{
		std::ofstream stream = std::ofstream(output_file+"_LS.gf", std::ios::out);
		stream.precision(15);

		// Header
		stream << "FiniteElementSpace\n";
		stream << "FiniteElementCollection: H1_3D_P1\n";
		stream << "VDim: 1\n";
		stream << "Ordering: 0\n";
		stream << "\n";

		for (auto n : m.nodes()) {
			stream << var_dist->value(n) << "\n";
		}
	}

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write(output_file+"_BG.vtk");

}