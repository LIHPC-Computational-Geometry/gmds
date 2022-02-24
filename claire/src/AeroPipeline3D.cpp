//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroPipeline3D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <unit_test_config.h>
#include <iostream>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroPipeline3D::AeroPipeline3D(ParamsAero Aparams) :
	AbstractAeroPipeline(Aparams),
  m_mTetra(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
                                        F2E | E2F | R2E | N2R | N2F | N2E)),
  m_mHexa( gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
                         F2E | E2F | R2E | N2R | N2F | N2E) )
{
	m_couche_id = m_mHexa.newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	m_mesh = &m_mTetra;
	m_meshGen = &m_mHexa;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::execute(){
	LectureMaillage();
	EcritureMaillage();
	m_isOver = true;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::LectureMaillage(){

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(m_params.input_file);

	gmds::MeshDoctor doctor(m_mesh);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	// Ecriture du maillage généré
	gmds::IGMeshIOService ioService(m_meshGen);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	std::string dir(TEST_SAMPLES_DIR);
	vtkWriter.write(m_params.output_file);

	// Ecriture du maillage initial (tetra)
	ioService = m_mesh;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline3D_Test1_Tetra.vtk");

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::InitialisationFronts(){

	std::cout << "-> Initialisation des fronts" << std::endl;

	// Marques sur les deux fronts d'intérêt pour l'aéro :
	// Paroi -> Noeuds sur la paroi
	// Ext -> Noeuds sur la frontière extérieure
	m_markFrontNodesParoi = m_mesh->newMark<gmds::Node>();
	m_markFrontNodesExt = m_mesh->newMark<gmds::Node>();

}
/*------------------------------------------------------------------------*/