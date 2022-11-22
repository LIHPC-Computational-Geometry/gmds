//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroPipeline_3D.h>
#include <gmds/claire/AeroBoundaries_3D.h>
#include <gmds/claire/AeroExtrusion_3D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <iostream>
#include <chrono>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroPipeline_3D::AeroPipeline_3D(ParamsAero Aparams) :
	AbstractAeroPipeline(Aparams)
{
	m_meshTet = new Mesh(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                  F2E | E2F | R2E | N2R | N2F | N2E));
	m_meshHex = new Mesh(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                     F2E | E2F | R2E | N2R | N2F | N2E));
	m_couche_id = m_meshHex->newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	m_meshHex->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	m_Bnd = new AeroBoundaries_3D(m_meshTet) ;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::STATUS
AeroPipeline_3D::execute(){

	clock_t t_start, t_end;

	LectureMaillage();

	t_start = clock();
	m_Bnd->execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	//m_manager.initAndLinkFrom3DMesh(&m_mTetra,&m_linker_TG);

	// Calcul du level set
	std::cout << "-> Calcul des Level Sets" << std::endl;
	t_start = clock();
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(m_meshTet, m_Bnd->getMarkParoi(), m_Bnd->getMarkAmont(),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	lsCombined.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;


	// Calcul du gradient du champ de Level Set
	std::cout << "-> Calcul du gradient du champ des Level Sets" << std::endl;
	t_start = clock();
	m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad3D(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                                       m_meshTet->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	grad3D.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;

	// Generate the blocking of the geometry surface.
	GeometrySurfaceBlockingGeneration();

	// Extrusion
	std::cout << "-> Extrusion" << std::endl;
	t_start = clock();
	AeroExtrusion_3D aero_extrusion(m_meshTet, m_meshHex, m_params,
	                                m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                                m_meshTet->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	aero_extrusion.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	// Write the final mesh.
	EcritureMaillage();

	return AbstractAeroPipeline::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::LectureMaillage(){

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_meshTet);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(m_params.input_file);

	gmds::MeshDoctor doctor(m_meshTet);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	// Ecriture du maillage généré
	gmds::IGMeshIOService ioService(m_meshHex);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(m_params.output_file);

	// Ecriture du maillage initial (tetra)
	ioService = m_meshTet;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline3D_Tetra.vtk");

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::GeometrySurfaceBlockingGeneration()
{
	//-------------------------------------//
	// Special case of the C2_3D geometry  //
	//	Surface Blocking 1						//
	//-------------------------------------//
	{
		Node n0 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n1 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n2 = m_meshHex->newNode({0.5, 0.5, -0.5});
		Node n3 = m_meshHex->newNode({0.5, -0.5, -0.5});

		Node n4 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n5 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n6 = m_meshHex->newNode({0.5, 0.5, 0.5});
		Node n7 = m_meshHex->newNode({0.5, -0.5, 0.5});

		Face f0 = m_meshHex->newQuad(n0, n1, n2, n3);     // F->N (x4)
		n0.add<Face>(f0);
		n1.add<Face>(f0);
		n2.add<Face>(f0);
		n3.add<Face>(f0);

		Face f1 = m_meshHex->newQuad(n4, n5, n6, n7);
		n4.add<Face>(f1);
		n5.add<Face>(f1);
		n6.add<Face>(f1);
		n7.add<Face>(f1);

		Face f2 = m_meshHex->newQuad(n0, n1, n5, n4);
		n0.add<Face>(f2);
		n1.add<Face>(f2);
		n5.add<Face>(f2);
		n4.add<Face>(f2);

		Face f3 = m_meshHex->newQuad(n0, n3, n7, n4);
		n0.add<Face>(f3);
		n3.add<Face>(f3);
		n7.add<Face>(f3);
		n4.add<Face>(f3);

		Face f4 = m_meshHex->newQuad(n3, n2, n6, n7);
		n3.add<Face>(f4);
		n2.add<Face>(f4);
		n6.add<Face>(f4);
		n7.add<Face>(f4);

		Face f5 = m_meshHex->newQuad(n2, n1, n5, n6);
		n2.add<Face>(f5);
		n1.add<Face>(f5);
		n5.add<Face>(f5);
		n6.add<Face>(f5);

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}
	}


	//-------------------------------------//
	// Special case of the C2_3D geometry  //
	//	Surface Blocking 2						//
	//-------------------------------------//
	/*
	{
		Node n1 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n2 = m_meshHex->newNode({0.0, -0.5, -0.5});
		Node n3 = m_meshHex->newNode({0.5, -0.5, -0.5});

		Node n4 = m_meshHex->newNode({-0.5, 0, -0.5});
		Node n5 = m_meshHex->newNode({0.0, 0, -0.5});
		Node n6 = m_meshHex->newNode({0.5, 0, -0.5});

		Node n7 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n8 = m_meshHex->newNode({0.0, 0.5, -0.5});
		Node n9 = m_meshHex->newNode({0.5, 0.5, -0.5});

		Node n10 = m_meshHex->newNode({-0.5, -0.5, 0});
		Node n11 = m_meshHex->newNode({0.0, -0.5, 0});
		Node n12 = m_meshHex->newNode({0.5, -0.5, 0});

		Node n13 = m_meshHex->newNode({-0.5, 0, 0});
		Node n14 = m_meshHex->newNode({0.5, 0, 0});

		Node n15 = m_meshHex->newNode({-0.5, 0.5, 0});
		Node n16 = m_meshHex->newNode({0.0, 0.5, 0});
		Node n17 = m_meshHex->newNode({0.5, 0.5, 0});

		Node n18 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n19 = m_meshHex->newNode({0.0, -0.5, 0.5});
		Node n20 = m_meshHex->newNode({0.5, -0.5, 0.5});

		Node n21 = m_meshHex->newNode({-0.5, 0, 0.5});
		Node n22 = m_meshHex->newNode({0.0, 0, 0.5});
		Node n23 = m_meshHex->newNode({0.5, 0, 0.5});

		Node n24 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n25 = m_meshHex->newNode({0.0, 0.5, 0.5});
		Node n26 = m_meshHex->newNode({0.5, 0.5, 0.5});



		Face f1 = m_meshHex->newQuad(n1, n2, n5, n4);     // F->N (x4)
		n1.add<Face>(f1);
		n2.add<Face>(f1);
		n5.add<Face>(f1);
		n4.add<Face>(f1);

		Face f2 = m_meshHex->newQuad(n2, n3, n6, n5);
		n2.add<Face>(f2);
		n3.add<Face>(f2);
		n6.add<Face>(f2);
		n5.add<Face>(f2);

		Face f3 = m_meshHex->newQuad(n4, n5, n8, n7);
		n4.add<Face>(f3);
		n5.add<Face>(f3);
		n8.add<Face>(f3);
		n7.add<Face>(f3);

		Face f4 = m_meshHex->newQuad(n5, n6, n9, n8);
		n5.add<Face>(f4);
		n6.add<Face>(f4);
		n9.add<Face>(f4);
		n8.add<Face>(f4);

		Face f4 = m_meshHex->newQuad(n3, n2, n6, n7);
		n3.add<Face>(f4);
		n2.add<Face>(f4);
		n6.add<Face>(f4);
		n7.add<Face>(f4);

		Face f5 = m_meshHex->newQuad(n2, n1, n5, n6);
		n2.add<Face>(f5);
		n1.add<Face>(f5);
		n5.add<Face>(f5);
		n6.add<Face>(f5);

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}
	}
	*/


}

