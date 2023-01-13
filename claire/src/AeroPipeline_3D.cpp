//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
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
	                                  F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	m_meshHex = new Mesh(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                     F2E | E2F | R2E | E2R | N2R | N2F | N2E));
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
	LeastSquaresGradientComputation grad3D(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
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

	// Erase node connected to nothing
	math::Utils::MeshCleaner(m_meshTet);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	// Ecriture du maillage généré
	gmds::IGMeshIOService ioService(m_meshHex);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
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
	if (m_params.block_surface_3D==1)
	{
		Node n0 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n1 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n2 = m_meshHex->newNode({0.5, 0.5, -0.5});
		Node n3 = m_meshHex->newNode({0.5, -0.5, -0.5});

		Node n4 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n5 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n6 = m_meshHex->newNode({0.5, 0.5, 0.5});
		Node n7 = m_meshHex->newNode({0.5, -0.5, 0.5});

		// Creates the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n1.id(), n2.id(), n3.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n5.id(), n6.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n1.id(), n5.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n3.id(), n7.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n2.id(), n6.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n1.id(), n5.id(), n6.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}
	}

	//-------------------------------------//
	// Special case of the C2_3D geometry  //
	//	Surface Blocking 2						//
	//-------------------------------------//

	else if (m_params.block_surface_3D==2)
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

		// Create the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n5.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n3.id(), n6.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n5.id(), n8.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n6.id(), n9.id(), n8.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n19.id(), n22.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n19.id(), n20.id(), n23.id(), n22.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n22.id(), n25.id(), n24.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n22.id(), n23.id(), n26.id(), n25.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n11.id(), n10.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n3.id(), n12.id(), n11.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n11.id(), n19.id(), n18.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n12.id(), n20.id(), n19.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n6.id(), n14.id(), n12.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n9.id(), n17.id(), n14.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n14.id(), n23.id(), n20.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n14.id(), n17.id(), n26.id(), n23.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n9.id(), n8.id(), n16.id(), n17.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n7.id(), n15.id(), n16.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n17.id(), n16.id(), n25.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n16.id(), n15.id(), n24.id(), n25.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n1.id(), n10.id(), n13.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n4.id(), n13.id(), n15.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n10.id(), n18.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n15.id(), n13.id(), n21.id(), n24.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}
	}


	//-------------------------------------//
	// Special case of the C3_3D geometry  //
	//	Surface Blocking 1						//
	//-------------------------------------//

	else if (m_params.block_surface_3D==3)
	{
		Node n0 = m_meshHex->newNode({0, 0, 0});

		Node n1 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n2 = m_meshHex->newNode({0.5, -0.5, 0.5});
		Node n3 = m_meshHex->newNode({0.0, 0.0, 0.5});
		Node n4 = m_meshHex->newNode({0.5, 0.0, 0.5});
		Node n5 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n6 = m_meshHex->newNode({0.0, 0.5, 0.5});

		Node n7 = m_meshHex->newNode({0.5, 0.0, 0.0});
		Node n8 = m_meshHex->newNode({0.0, 0.5, 0.0});
		Node n9 = m_meshHex->newNode({0.5, 0.5, 0.0});

		Node n10 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n11 = m_meshHex->newNode({0.5, -0.5, -0.5});
		Node n12 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n13 = m_meshHex->newNode({0.5, 0.5, -0.5});

		// Create the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n11.id(), n10.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n10.id(), n12.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n11.id(), n13.id(), n12.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n4.id(), n7.id(), n0.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n7.id(), n9.id(), n8.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n0.id(), n8.id(), n6.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n3.id(), n6.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n2.id(), n4.id(), n3.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n6.id(), n8.id(), n12.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n9.id(), n13.id(), n12.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n4.id(), n7.id(), n11.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n9.id(), n13.id(), n11.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

		gmds::IGMeshIOService ioService(m_meshHex);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N|gmds::F);
		vtkWriter.setDataOptions(gmds::N|gmds::F);
		vtkWriter.write("test.vtk");

	}


	//-------------------------------------//
	// Special case of the C3_3D geometry  //
	//	Surface Blocking 2						//
	//-------------------------------------//

	else if (m_params.block_surface_3D==4)
	{
		Node n0 = m_meshHex->newNode({0, 0, 0});

		Node n1 = m_meshHex->newNode({-0.5, -0.5, 0.5});
		Node n2 = m_meshHex->newNode({0.5, -0.5, 0.5});
		Node n3 = m_meshHex->newNode({0.0, 0.0, 0.5});
		Node n4 = m_meshHex->newNode({0.5, 0.0, 0.5});
		Node n5 = m_meshHex->newNode({-0.5, 0.5, 0.5});
		Node n6 = m_meshHex->newNode({0.0, 0.5, 0.5});

		Node n7 = m_meshHex->newNode({0.5, 0.0, 0.0});
		Node n8 = m_meshHex->newNode({0.0, 0.5, 0.0});
		Node n9 = m_meshHex->newNode({0.5, 0.5, 0.0});

		Node n10 = m_meshHex->newNode({-0.5, -0.5, -0.5});
		Node n11 = m_meshHex->newNode({0.5, -0.5, -0.5});
		Node n12 = m_meshHex->newNode({-0.5, 0.5, -0.5});
		Node n13 = m_meshHex->newNode({0.5, 0.5, -0.5});

		Node n20 = m_meshHex->newNode({0, -0.5, 0.5});
		Node n21 = m_meshHex->newNode({-0.5, -0.5, 0});
		Node n22 = m_meshHex->newNode({0.0, -0.5, 0});
		Node n23 = m_meshHex->newNode({0.5, -0.5, 0});
		Node n24 = m_meshHex->newNode({0.0, -0.5, -0.5});
		Node n25 = m_meshHex->newNode({0.25, 0.0, 0.5});
		Node n26 = m_meshHex->newNode({0.0, 0.0, 0.25});
		Node n27 = m_meshHex->newNode({0.25, 0.0, 0.25});
		Node n28 = m_meshHex->newNode({0.5, 0.0, 0.25});
		Node n29 = m_meshHex->newNode({0.25, 0.0, 0.0});
		Node n30 = m_meshHex->newNode({-0.5, 0.0, 0.5});
		Node n31 = m_meshHex->newNode({0.0, 0.25, 0.5});
		Node n32 = m_meshHex->newNode({-0.5, 0.0, 0.0});
		Node n33 = m_meshHex->newNode({0.0, 0.25, 0.25});
		Node n34 = m_meshHex->newNode({0.0, 0.25, 0.0});
		Node n35 = m_meshHex->newNode({0.25, 0.25, 0.0});
		Node n36 = m_meshHex->newNode({0.5, 0.25, 0.0});
		Node n37 = m_meshHex->newNode({0.5, 0.0, -0.5});
		Node n38 = m_meshHex->newNode({-0.5, 0.0, -0.5});
		Node n39 = m_meshHex->newNode({0.0, 0.0, -0.5});
		Node n40 = m_meshHex->newNode({0.0, 0.5, 0.25});
		Node n41 = m_meshHex->newNode({-0.5, 0.5, 0.0});
		Node n42 = m_meshHex->newNode({0.25, 0.5, 0.0});
		Node n43 = m_meshHex->newNode({0.0, 0.5, -0.5});

		// Create the faces
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n3.id(), n31.id(), n30.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n30.id(), n31.id(), n6.id(), n5.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n20.id(), n25.id(), n3.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n2.id(), n4.id(), n25.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n29.id(), n35.id(), n34.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n29.id(), n7.id(), n36.id(), n35.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n34.id(), n35.id(), n42.id(), n8.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n36.id(), n9.id(), n42.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n24.id(), n39.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n24.id(), n11.id(), n37.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n43.id(), n39.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n39.id(), n37.id(), n13.id(), n43.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n4.id(), n28.id(), n23.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n7.id(), n28.id(), n23.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n7.id(), n36.id(), n37.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n9.id(), n36.id(), n37.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n31.id(), n33.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n34.id(), n33.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n40.id(), n33.id(), n34.id(), n8.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n40.id(), n33.id(), n31.id(), n6.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n30.id(), n32.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n38.id(), n32.id(), n21.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n38.id(), n32.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n30.id(), n32.id(), n41.id());

		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n6.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n8.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n8.id(), n42.id(), n43.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n9.id(), n42.id(), n43.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n29.id(), n27.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n25.id(), n27.id(), n26.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n28.id(), n27.id(), n25.id(), n4.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n28.id(), n27.id(), n29.id(), n7.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n24.id(), n22.id(), n21.id(), n10.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n24.id(), n22.id(), n23.id(), n11.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n22.id(), n23.id(), n2.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n22.id(), n21.id(), n1.id());

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

		gmds::IGMeshIOService ioService(m_meshHex);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N|gmds::F);
		vtkWriter.setDataOptions(gmds::N|gmds::F);
		vtkWriter.write("Surface_3D.vtk");

	}

	else if (m_params.block_surface_3D==5)
	{
		Node n0 = m_meshHex->newNode({-1, -1, 1});
		Node n1 = m_meshHex->newNode({0, -1, 1});
		Node n2 = m_meshHex->newNode({1, -1, 1});
		Node n3 = m_meshHex->newNode({-1, 0, 1});
		Node n4 = m_meshHex->newNode({0, 0, 1});
		Node n5 = m_meshHex->newNode({1, 0, 1});
		Node n6 = m_meshHex->newNode({-1, -1, 0});
		Node n7 = m_meshHex->newNode({0, -1, 0});
		Node n8 = m_meshHex->newNode({1, -1, 0});
		Node n9 = m_meshHex->newNode({-1, 0, 0});
		Node n10 = m_meshHex->newNode({0, 0, 0});
		Node n11 = m_meshHex->newNode({1, 0, 0});
		Node n12 = m_meshHex->newNode({-1, 1, 0});
		Node n13 = m_meshHex->newNode({0, 1, 0});
		Node n14 = m_meshHex->newNode({-1, -1, -1});
		Node n15 = m_meshHex->newNode({0, -1, -1});
		Node n16 = m_meshHex->newNode({1, -1, -1});
		Node n17 = m_meshHex->newNode({-1, 0, -1});
		Node n18 = m_meshHex->newNode({0, 0, -1});
		Node n19 = m_meshHex->newNode({1, 0, -1});
		Node n20 = m_meshHex->newNode({-1, 1, -1});
		Node n21 = m_meshHex->newNode({0, 1, -1});
		Node n30 = m_meshHex->newNode({-0.5, -1, 1});
		Node n31 = m_meshHex->newNode({0.5, -1, 1});
		Node n32 = m_meshHex->newNode({-1, -1, 0.5});
		Node n33 = m_meshHex->newNode({-0.5, -1, 0.5});
		Node n34 = m_meshHex->newNode({0, -1, 0.5});
		Node n35 = m_meshHex->newNode({0.5, -1, 0.5});
		Node n36 = m_meshHex->newNode({1, -1, 0.5});
		Node n37 = m_meshHex->newNode({-0.5, -1, 0});
		Node n38 = m_meshHex->newNode({0.5, -1, 0});
		Node n39 = m_meshHex->newNode({-1, -1, -0.5});
		Node n40 = m_meshHex->newNode({-0.5, -1, -0.5});
		Node n41 = m_meshHex->newNode({0, -1, -0.5});
		Node n42 = m_meshHex->newNode({0.5, -1, -0.5});
		Node n43 = m_meshHex->newNode({1, -1, -0.5});
		Node n44 = m_meshHex->newNode({-0.5, -1, -1});
		Node n45 = m_meshHex->newNode({0.5, -1, -1});
		Node n46 = m_meshHex->newNode({-0.5, 0, 1});
		Node n47 = m_meshHex->newNode({0.5, 0, 1});
		Node n48 = m_meshHex->newNode({-1, 0, 0.5});
		Node n49 = m_meshHex->newNode({-0.5, 0, 0.5});
		Node n50 = m_meshHex->newNode({0, 0, 0.5});
		Node n51 = m_meshHex->newNode({0.5, 0, 0.5});
		Node n52 = m_meshHex->newNode({1, 0, 0.5});
		Node n53 = m_meshHex->newNode({-0.5, 0, 0});
		Node n54 = m_meshHex->newNode({0.5, 0, 0});
		Node n55 = m_meshHex->newNode({-1, 0, -0.5});
		Node n56 = m_meshHex->newNode({0, 0, -0.5});
		Node n57 = m_meshHex->newNode({0.5, 0, -0.5});
		Node n58 = m_meshHex->newNode({1, 0, -0.5});
		Node n59 = m_meshHex->newNode({-0.5, 0, -1});
		Node n60 = m_meshHex->newNode({0.5, 0, -1});
		Node n61 = m_meshHex->newNode({-1, 0.5, 0});
		Node n62 = m_meshHex->newNode({-0.5, 0.5, 0});
		Node n63 = m_meshHex->newNode({0, 0.5, 0});
		Node n64 = m_meshHex->newNode({-1, 0.5, -0.5});
		Node n65 = m_meshHex->newNode({0, 0.5, -0.5});
		Node n66 = m_meshHex->newNode({-1, 0.5, -1});
		Node n67 = m_meshHex->newNode({-0.5, 0.5, -1});
		Node n68 = m_meshHex->newNode({0, 0.5, -1});
		Node n69 = m_meshHex->newNode({-0.5, 1, 0});
		Node n70 = m_meshHex->newNode({-1, 1, -0.5});
		Node n71 = m_meshHex->newNode({-0.5, 1, -0.5});
		Node n72 = m_meshHex->newNode({0, 1, -0.5});
		Node n73 = m_meshHex->newNode({-0.5, 1, -1});

		// Creates the faces
		// x=-1
		TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n32.id(), n48.id(), n3.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n48.id(), n32.id(), n6.id(), n9.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n9.id(), n55.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n14.id(), n17.id(), n55.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n55.id(), n17.id(), n66.id(), n64.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n61.id(), n9.id(), n55.id(), n64.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n61.id(), n12.id(), n70.id(), n64.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n66.id(), n20.id(), n70.id(), n64.id());
		// x=1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n5.id(), n52.id(), n36.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n11.id(), n52.id(), n36.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n8.id(), n11.id(), n58.id(), n43.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n16.id(), n19.id(), n58.id(), n43.id());
		// x=0
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n63.id(), n65.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n68.id(), n65.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n68.id(), n65.id(), n72.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n63.id(), n65.id(), n72.id());

		// z=-1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n17.id(), n14.id(), n44.id(), n59.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n15.id(), n44.id(), n59.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n15.id(), n45.id(), n60.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n19.id(), n16.id(), n45.id(), n60.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n17.id(), n59.id(), n67.id(), n66.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n73.id(), n67.id(), n66.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n73.id(), n67.id(), n68.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n59.id(), n67.id(), n68.id());
		//z=1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n3.id(), n46.id(), n30.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n4.id(), n46.id(), n30.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n4.id(), n47.id(), n31.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n2.id(), n5.id(), n47.id(), n31.id());
		// z=0
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n63.id(), n62.id(), n53.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n9.id(), n61.id(), n62.id(), n53.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n61.id(), n62.id(), n69.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n63.id(), n62.id(), n69.id());

		// y=-1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n14.id(), n44.id(), n40.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n37.id(), n40.id(), n39.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n6.id(), n37.id(), n33.id(), n32.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n0.id(), n30.id(), n33.id(), n32.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n1.id(), n30.id(), n33.id(), n34.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n37.id(), n33.id(), n34.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n7.id(), n37.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n15.id(), n44.id(), n40.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n45.id(), n15.id(), n41.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n41.id(), n7.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n34.id(), n7.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n34.id(), n1.id(), n31.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n36.id(), n2.id(), n31.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n35.id(), n36.id(), n8.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n43.id(), n8.id(), n38.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n42.id(), n43.id(), n16.id(), n45.id());
		// y = 0
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n3.id(), n46.id(), n49.id(), n48.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n9.id(), n53.id(), n49.id(), n48.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n53.id(), n49.id(), n50.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n46.id(), n49.id(), n50.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n4.id(), n47.id(), n51.id(), n50.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n5.id(), n47.id(), n51.id(), n52.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n54.id(), n51.id(), n52.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n11.id(), n54.id(), n57.id(), n58.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n19.id(), n60.id(), n57.id(), n58.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n18.id(), n60.id(), n57.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n54.id(), n57.id(), n56.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n10.id(), n54.id(), n51.id(), n50.id());
		// y = 1
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n12.id(), n69.id(), n71.id(), n70.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n20.id(), n73.id(), n71.id(), n70.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n21.id(), n73.id(), n71.id(), n72.id());
		f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshHex, n13.id(), n69.id(), n71.id(), n72.id());


		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

		gmds::IGMeshIOService ioService(m_meshHex);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N|gmds::F);
		vtkWriter.setDataOptions(gmds::N|gmds::F);
		vtkWriter.write("Surface_3D.vtk");

	}

}
/*------------------------------------------------------------------------*/

