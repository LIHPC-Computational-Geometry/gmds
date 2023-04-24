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
#include <gmds/claire/IntervalAssignment_3D.h>
#include <gmds/claire/MFEMMeshWriter.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/smoothy/LaplacianSmoother.h>
#include <iostream>
#include <chrono>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroPipeline_3D::AeroPipeline_3D(std::string Aparams, std::string &Aworking_dir) :
	AbstractAeroPipeline(Aparams, Aworking_dir)
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
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	//m_manager.initAndLinkFrom3DMesh(&m_mTetra,&m_linker_TG);
	PreTraitementMeshTet();

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
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;


	// Calcul du gradient du champ de Level Set
	std::cout << "-> Calcul du gradient du champ des Level Sets" << std::endl;
	t_start = clock();
	ComputeVectorFieldForExtrusion();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;

	GeometrySurfaceBlockingGeneration();	// Generate the blocking of the geometry surface.
	SurfaceBlockingClassification();			// Link the surface blocking to the geometry.

	// Write the surface blocking to check the classification
	if (m_params.with_debug_files)
	{
		gmds::IGMeshIOService ioService(m_meshHex);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F);
		vtkWriter.setDataOptions(gmds::N | gmds::F);
		vtkWriter.write("Surface_3D_GEOM_CLASSIFICATION.vtk");
	}

	// Extrusion
	std::cout << "-> Extrusion" << std::endl;
	Variable<math::Vector3d>* var_VectorsForExtrusion = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");
	t_start = clock();
	AeroExtrusion_3D aero_extrusion(m_meshTet, m_meshHex, m_params,
	                                m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                                var_VectorsForExtrusion);
	aero_extrusion.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	BlockingGeometricClassification();	// Geometric classification of the final blocking, according to the geometric classification of the input tet mesh
	//smoothy::LaplacianSmoother smoother(m_linker_HG);
	//smoother.smoothCurves();
	//smoother.smoothSurfaces();
	//smoother.smoothVolumes(10);

	// Interval Assignment
	std::cout << "-> Interval Assignment" << std::endl;
	t_start = clock();
	IntervalAssignment_3D intAss(m_meshHex, m_params,
	                                m_meshHex->newVariable<int,GMDS_EDGE>("GMDS_EdgeDiscretization"));
	intAss.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*double(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
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
void
AeroPipeline_3D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	// Ecriture du maillage généré
	gmds::IGMeshIOService ioService(m_meshHex);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write(m_params.output_file);

	gmds::VTKWriter vtkWriter_edges(&ioService);
	vtkWriter_edges.setCellOptions(gmds::N|gmds::E);
	vtkWriter_edges.setDataOptions(gmds::N|gmds::E);
	vtkWriter_edges.write("EdgesDiscretization.vtk");

	// Ecriture du maillage initial (tetra)
	ioService = IGMeshIOService(m_meshTet);
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline3D_Tetra.vtk");

	// TEST FOR MFEM
	/*
	if (m_params.with_debug_files)
	{
		MeshDoctor doc(m_meshHex);
		doc.buildFacesAndR2F();
		doc.updateUpwardConnectivity();
		doc.orient2DFaces();

		std::cout << "Nbr elements avant: " << m_meshHex->getNbRegions() << std::endl;
		for (auto r_id:m_meshHex->regions())
		{
			math::Utils::orientRegion(m_meshHex, m_meshHex->get<Region>(r_id));
		}
		std::cout << "Nbr elements après: " << m_meshHex->getNbRegions() << std::endl;

		gmds::IGMeshIOService ioService_MFEM(m_meshHex);
		gmds::VTKWriter vtkWriter_MFEM(&ioService_MFEM);
		vtkWriter_MFEM.setCellOptions(gmds::N | gmds::R);
		vtkWriter_MFEM.setDataOptions(gmds::N | gmds::R);
		std::string dir(".");
		vtkWriter_MFEM.write("Apollo_3D_Blocks.vtk");
	}
	 */

	{
		std::cout << "MFEM writing..." << std::endl;
		FastLocalize fl = FastLocalize(m_meshTet);
		Variable<int>* var_couche = m_meshHex->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
		Variable<math::Vector3d>* var_VectorsForExtrusion = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");
		std::map<TCellID,TCellID> new_nodes;
		for (auto n_id:m_meshHex->nodes())
		{
			if (var_couche->value(n_id)==0)
			{
				Cell::Data data = fl.find(m_meshHex->get<Node>(n_id).point());
				//Node n_closest = m_meshHex->get<Node>(data.id);
				math::Vector3d v = var_VectorsForExtrusion->value(data.id).normalize();
				math::Point p_new = m_meshHex->get<Node>(n_id).point() + -1*(m_params.delta_cl/2.0)*v;
				Node n_new = m_meshHex->newNode(p_new);
				new_nodes[n_id] = n_new.id();
				var_couche->set(n_new.id(), -1);
			}
		}
		for (auto f_id:m_meshHex->faces())
		{
			Face f = m_meshHex->get<Face>(f_id);
			std::vector<Node> f_nodes = f.get<Node>();
			if (var_couche->value(f_nodes[0].id())==0
			    && var_couche->value(f_nodes[1].id())==0
			    && var_couche->value(f_nodes[2].id())==0
			    && var_couche->value(f_nodes[3].id())==0)
			{
				Region r = m_meshHex->newHex(f_nodes[0].id(), f_nodes[1].id(), f_nodes[2].id(), f_nodes[3].id(),
				                  new_nodes[f_nodes[0].id()], new_nodes[f_nodes[1].id()], new_nodes[f_nodes[2].id()], new_nodes[f_nodes[3].id()]);
				var_couche->set(r.id(), -1);
			}
		}

		/*
		math::Point p;
		int compteur(0);
		for (auto n_id:m_meshHex->nodes())
		{
			if (var_couche->value(n_id)==0)
			{
				Node n = m_meshHex->get<Node>(n_id);
				p = p + n.point();
				compteur++;
			}
		}
		if (compteur != 0)
		{
			p.setX(p.X() / compteur);
			p.setY(p.Y() / compteur);
		}
		Node n_new = m_meshHex->newNode(p);
		for (auto f_id:m_meshHex->faces())
		{
			Face f = m_meshHex->get<Face>(f_id);
			std::vector<Node> f_nodes = f.get<Node>();
			if (var_couche->value(f_nodes[0].id())==0
			    && var_couche->value(f_nodes[1].id())==0
			    && var_couche->value(f_nodes[2].id())==0
			    && var_couche->value(f_nodes[3].id())==0)
			{
				m_meshHex->newPyramid(f_nodes[0], f_nodes[1], f_nodes[2], f_nodes[3], n_new);
			}
		}
		*/

		if (m_params.with_debug_files)
		{
			MeshDoctor doc(m_meshHex);
			doc.buildFacesAndR2F();
			doc.updateUpwardConnectivity();
			doc.orient2DFaces();

			std::cout << "Nbr elements avant: " << m_meshHex->getNbRegions() << std::endl;
			for (auto r_id:m_meshHex->regions())
			{
				math::Utils::orientRegion(m_meshHex, m_meshHex->get<Region>(r_id));
			}
			std::cout << "Nbr elements après: " << m_meshHex->getNbRegions() << std::endl;

			gmds::IGMeshIOService ioService_MFEM(m_meshHex);
			gmds::VTKWriter vtkWriter_MFEM(&ioService_MFEM);
			vtkWriter_MFEM.setCellOptions(gmds::N | gmds::R);
			vtkWriter_MFEM.setDataOptions(gmds::N | gmds::R);
			std::string dir(".");
			vtkWriter_MFEM.write("AeroPipeline3D_Blocking_withPyramids.vtk");
		}

		MFEMMeshWriter mfemwriter = MFEMMeshWriter(m_meshHex, "AeroPipeline3D_Blocking_toFit");
		mfemwriter.execute();
	}

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::GeometrySurfaceBlockingGeneration()
{
	if (m_params.block_surface_3D==0)
	{
		// Lecture du maillage
		std::cout << "-> Lecture du maillage de surface ..." << std::endl;

		gmds::IGMeshIOService ioService(m_meshHex);
		gmds::VTKReader vtkReader(&ioService);
		vtkReader.setCellOptions(gmds::N|gmds::F);
		vtkReader.read(m_params.input_file_3D_surface);

		for (auto n_id : m_meshHex->nodes()) {
			m_couche_id->set(n_id, 0);
		}

		// Build the edges and the connectivities
		math::Utils::buildEfromFandConnectivies(m_meshHex);

	}

	//-------------------------------------//
	// Special case of the C2_3D geometry  //
	//	Surface Blocking 1						//
	//-------------------------------------//
	else if (m_params.block_surface_3D==1)
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

	}

	// Build the edges and the connectivities
	math::Utils::buildEfromFandConnectivies(m_meshHex);

	// Write the surface block structure
	gmds::IGMeshIOService ioService(m_meshHex);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("Surface_3D.vtk");

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::PreTraitementMeshTet()
{
	for (auto r_id:m_meshTet->regions())
	{
		Region r = m_meshTet->get<Region>(r_id);
		std::vector<Node> r_nodes = r.get<Node>();
		if ( m_Bnd->getNodeColor(r_nodes[0].id()) != 0
		    && (m_Bnd->getNodeColor(r_nodes[0].id()) == m_Bnd->getNodeColor(r_nodes[1].id()))
		    && (m_Bnd->getNodeColor(r_nodes[0].id()) == m_Bnd->getNodeColor(r_nodes[2].id()))
		    && (m_Bnd->getNodeColor(r_nodes[0].id()) == m_Bnd->getNodeColor(r_nodes[3].id())) )
		{
			//std::cout << "Tetra " << r_id << " pré traité" << std::endl;

			math::Point P_bar = r.center();
			m_meshTet->deleteRegion(r);
			Node n_new = m_meshTet->newNode(P_bar);

			r_nodes[0].remove<Region>(r);
			r_nodes[1].remove<Region>(r);
			r_nodes[2].remove<Region>(r);
			r_nodes[3].remove<Region>(r);

			// Get the 4 faces of the tet
			TCellID f1_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[0].id(), r_nodes[1].id(), r_nodes[2].id());
			TCellID f2_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[1].id(), r_nodes[2].id(), r_nodes[3].id());
			TCellID f3_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[0].id(), r_nodes[2].id(), r_nodes[3].id());
			TCellID f4_id = math::Utils::CommonFace3Nodes(m_meshTet, r_nodes[0].id(), r_nodes[1].id(), r_nodes[2].id());

			// New tets
			Region r_1 = m_meshTet->newTet(n_new, r_nodes[0], r_nodes[1], r_nodes[2]);
			Region r_2 = m_meshTet->newTet(n_new, r_nodes[1], r_nodes[2], r_nodes[3]);
			Region r_3 = m_meshTet->newTet(n_new, r_nodes[0], r_nodes[2], r_nodes[3]);
			Region r_4 = m_meshTet->newTet(n_new, r_nodes[0], r_nodes[1], r_nodes[3]);

			// The 3 new faces
			Face f_new_01 = m_meshTet->newTriangle(n_new, r_nodes[0], r_nodes[1]);
			Face f_new_02 = m_meshTet->newTriangle(n_new, r_nodes[0], r_nodes[2]);
			Face f_new_03 = m_meshTet->newTriangle(n_new, r_nodes[0], r_nodes[3]);
			Face f_new_12 = m_meshTet->newTriangle(n_new, r_nodes[1], r_nodes[2]);
			Face f_new_13 = m_meshTet->newTriangle(n_new, r_nodes[1], r_nodes[3]);
			Face f_new_23 = m_meshTet->newTriangle(n_new, r_nodes[2], r_nodes[3]);

			// The 4 new edges
			Edge e_new_0 = m_meshTet->newEdge(n_new, r_nodes[0]);
			Edge e_new_1 = m_meshTet->newEdge(n_new, r_nodes[1]);
			Edge e_new_2 = m_meshTet->newEdge(n_new, r_nodes[2]);
			Edge e_new_3 = m_meshTet->newEdge(n_new, r_nodes[3]);

			// Get the old tet edges
			TCellID e_01_id = math::Utils::CommonEdge(m_meshTet, r_nodes[0].id(), r_nodes[1].id()) ;
			TCellID e_02_id = math::Utils::CommonEdge(m_meshTet, r_nodes[0].id(), r_nodes[2].id()) ;
			TCellID e_03_id = math::Utils::CommonEdge(m_meshTet, r_nodes[0].id(), r_nodes[2].id()) ;
			TCellID e_12_id = math::Utils::CommonEdge(m_meshTet, r_nodes[1].id(), r_nodes[2].id()) ;
			TCellID e_13_id = math::Utils::CommonEdge(m_meshTet, r_nodes[1].id(), r_nodes[3].id()) ;
			TCellID e_23_id = math::Utils::CommonEdge(m_meshTet, r_nodes[2].id(), r_nodes[3].id()) ;

			// Update the connectivities N->R
			r_nodes[0].add<Region>(r_1);
			r_nodes[0].add<Region>(r_3);
			r_nodes[0].add<Region>(r_4);

			r_nodes[1].add<Region>(r_1);
			r_nodes[1].add<Region>(r_2);
			r_nodes[1].add<Region>(r_4);

			r_nodes[2].add<Region>(r_1);
			r_nodes[2].add<Region>(r_2);
			r_nodes[2].add<Region>(r_3);

			r_nodes[3].add<Region>(r_2);
			r_nodes[3].add<Region>(r_3);
			r_nodes[3].add<Region>(r_4);

			// Update the connectivities N->R (x4) for the new node
			n_new.add<Region>(r_1);
			n_new.add<Region>(r_2);
			n_new.add<Region>(r_3);
			n_new.add<Region>(r_4);

			// Update the connectivities, N->F (x6) and N->E (x6)
			n_new.add<Face>(f_new_01);
			n_new.add<Face>(f_new_02);
			n_new.add<Face>(f_new_03);
			n_new.add<Face>(f_new_12);
			n_new.add<Face>(f_new_13);
			n_new.add<Face>(f_new_23);

			n_new.add<Edge>(e_01_id);
			n_new.add<Edge>(e_02_id);
			n_new.add<Edge>(e_03_id);
			n_new.add<Edge>(e_12_id);
			n_new.add<Edge>(e_13_id);
			n_new.add<Edge>(e_23_id);

			// Update the connectivities for the nodes of the old tet
			r_nodes[0].add<Edge>(e_new_0);
			r_nodes[0].add<Face>(f_new_01);
			r_nodes[0].add<Face>(f_new_02);
			r_nodes[0].add<Face>(f_new_03);

			r_nodes[1].add<Edge>(e_new_1);
			r_nodes[1].add<Face>(f_new_01);
			r_nodes[1].add<Face>(f_new_12);
			r_nodes[1].add<Face>(f_new_13);

			r_nodes[2].add<Edge>(e_new_2);
			r_nodes[2].add<Face>(f_new_02);
			r_nodes[2].add<Face>(f_new_12);
			r_nodes[2].add<Face>(f_new_23);

			r_nodes[3].add<Edge>(e_new_3);
			r_nodes[3].add<Face>(f_new_03);
			r_nodes[3].add<Face>(f_new_13);
			r_nodes[3].add<Face>(f_new_23);

			// Update the connectivities R->F (x3) and R->E (x3)
			r_1.add<Face>(f1_id);
			r_1.add<Face>(f_new_01.id());
			r_1.add<Face>(f_new_12.id());
			r_1.add<Edge>(e_01_id);
			r_1.add<Edge>(e_12_id);
			r_1.add<Edge>(e_02_id);

			r_2.add<Face>(f2_id);
			r_2.add<Face>(f_new_12.id());
			r_2.add<Face>(f_new_23.id());
			r_2.add<Edge>(e_12_id);
			r_2.add<Edge>(e_02_id);
			r_2.add<Edge>(e_01_id);

			r_3.add<Face>(f3_id);
			r_3.add<Face>(f_new_02.id());
			r_3.add<Face>(f_new_23.id());
			r_3.add<Edge>(e_02_id);
			r_3.add<Edge>(e_23_id);
			r_3.add<Edge>(e_03_id);

			r_4.add<Face>(f4_id);
			r_4.add<Face>(f_new_01.id());
			r_4.add<Face>(f_new_13.id());
			r_4.add<Edge>(e_01_id);
			r_4.add<Edge>(e_13_id);
			r_4.add<Edge>(e_03_id);

		}
	}

	// Ecriture du maillage initial (tetra)
	gmds::IGMeshIOService ioService(m_meshTet);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroPipeline3D_Tetra_PreTraite.vtk");
}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::SurfaceBlockingClassification()
{
	// Init the geometry manager and the linker
	m_manager->initAndLinkFrom3DMesh(m_meshTet, m_linker_TG);

	// Init the linker for the Blocking
	m_linker_HG->setGeometry(m_manager);
	m_linker_HG->setMesh(m_meshHex);

	//------------------------------//
	//	BLOCK CORNERS CLASSIFICATION //
	//------------------------------//
	for (auto n_id:m_meshHex->nodes())
	{
		std::vector<cad::GeomPoint*> points;
		m_manager->getPoints(points);

		std::vector<cad::GeomCurve*> curves;
		m_manager->getCurves(curves);

		std::vector<cad::GeomSurface*> surfaces;
		m_manager->getSurfaces(surfaces);

		Node n = m_meshHex->get<Node>(n_id);
		math::Point p = n.point();
		math::Point p_proj = n.point();

		double min_dist(std::numeric_limits<double>::max());
		int geom_dim = surfaces[0]->dim();
		int geom_id = 3;

		// Check on the geometrical points first
		for (auto point:points)
		{
			p_proj = n.point();
			point->project(p_proj);
			if ( (n.point()-p_proj).norm() < min_dist )
			{
				min_dist = (n.point()-p_proj).norm();
				geom_dim = 1;
				geom_id = point->id();
			}
		}

		// Check on the geometrical curves
		for (auto curve:curves)
		{
			p_proj = n.point();
			curve->project(p_proj);
			if ( (n.point()-p_proj).norm() < min_dist )
			{
				min_dist = (n.point()-p_proj).norm();
				geom_dim = 2;
				geom_id = curve->id();
			}
		}

		// Check on the geometrical surfaces
		for (auto surface:surfaces)
		{
			p_proj = n.point();
			surface->project(p_proj);
			if ( (n.point()-p_proj).norm() < min_dist
			    && min_dist > pow(10,-20))		// Need to add a tolerence here.
			{
				min_dist = (n.point()-p_proj).norm();
				geom_dim = 3;
				geom_id = surface->id();
			}
		}

		// Link the node of the Blocking to the geometrical entity that
		// minimize the distance
		if(geom_dim == 1){
			m_linker_HG->linkNodeToPoint(n_id, geom_id);
		}
		else if(geom_dim==2){
			m_linker_HG->linkNodeToCurve(n_id, geom_id);
		}
		else if(geom_dim==3){
			m_linker_HG->linkNodeToSurface(n_id, geom_id);
		}

	}

	//------------------------------//
	//	 BLOCK EDGES CLASSIFICATION  //
	//------------------------------//
	for (auto e_id:m_meshHex->edges())
	{
		std::vector<Node> e_nodes = m_meshHex->get<Edge>(e_id).get<Node>();
		int geom_id_n0 = m_linker_HG->getGeomId<Node>(e_nodes[0].id()) ;
		int geom_id_n1 = m_linker_HG->getGeomId<Node>(e_nodes[1].id()) ;
		int geom_dim_n0 = m_linker_HG->getGeomDim<Node>(e_nodes[0].id()) ;
		int geom_dim_n1 = m_linker_HG->getGeomDim<Node>(e_nodes[1].id());

		int geom_dim;
		int geom_id;

		if (geom_dim_n0 > geom_dim_n1)
		{
			geom_dim = geom_dim_n0;
			geom_id = geom_id_n0;
		}
		else if (geom_dim_n1 > geom_dim_n0)
		{
			geom_dim = geom_dim_n1;
			geom_id = geom_id_n1;
		}
		else if (geom_dim_n0 == geom_dim_n1
		         && geom_id_n0 == geom_id_n1)
		{
			geom_dim = geom_dim_n0;
			geom_id = geom_id_n0;
		}
		else if (geom_dim_n0 == geom_dim_n1
		         && geom_dim_n1 == 1)	// Find the common curve between the two geometrical points
		{
			cad::GeomPoint* geom_p0 = m_manager->getPoint(geom_id_n0);
			cad::GeomPoint* geom_p1 = m_manager->getPoint(geom_id_n1);
			geom_id = m_manager->getCommonCurve(geom_p0, geom_p1 );
			geom_dim = 2;
		}

		// Link the edge of the Blocking to the geometrical entity
		if(geom_dim==2){
			m_linker_HG->linkEdgeToCurve(e_id, geom_id);
		}
		else if(geom_dim==3){
			m_linker_HG->linkEdgeToSurface(e_id, geom_id);
		}

	}

	//------------------------------//
	//	 BLOCK FACES CLASSIFICATION  //
	//------------------------------//
	for (auto f_id:m_meshHex->faces())
	{
		std::vector<Edge> f_edges = m_meshHex->get<Face>(f_id).get<Edge>();
		int geom_id_e0 = m_linker_HG->getGeomId<Edge>(f_edges[0].id()) ;
		int geom_id_e1 = m_linker_HG->getGeomId<Edge>(f_edges[1].id()) ;
		int geom_id_e2 = m_linker_HG->getGeomId<Edge>(f_edges[2].id()) ;
		int geom_id_e3 = m_linker_HG->getGeomId<Edge>(f_edges[3].id()) ;
		int geom_dim_e0 = m_linker_HG->getGeomDim<Edge>(f_edges[0].id()) ;
		int geom_dim_e1 = m_linker_HG->getGeomDim<Edge>(f_edges[1].id()) ;
		int geom_dim_e2 = m_linker_HG->getGeomDim<Edge>(f_edges[2].id()) ;
		int geom_dim_e3 = m_linker_HG->getGeomDim<Edge>(f_edges[3].id()) ;

		int geom_dim(3);
		int geom_id;

		if (geom_dim_e0 == 3)
		{
			geom_id = geom_id_e0;
		}
		else if (geom_dim_e1 == 3)
		{
			geom_id = geom_id_e1;
		}
		else if (geom_dim_e2 == 3)
		{
			geom_id = geom_id_e2;
		}
		else if (geom_dim_e3 == 3)
		{
			geom_id = geom_id_e3;
		}
		else
		{
			cad::GeomCurve* curve_0 = m_manager->getCurve(geom_id_e0);
			cad::GeomCurve* curve_1 = m_manager->getCurve(geom_id_e1);
			geom_id = m_manager->getCommonSurface(curve_0, curve_1);
		}
		// Link the face to the surface
		m_linker_HG->linkFaceToSurface(f_id, geom_id);
	}

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::ComputeVectorFieldForExtrusion()
{
	Variable<math::Vector3d>* var_VectorsForExtrusion = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");

	if (m_params.vectors_field <= 0 || m_params.vectors_field > 4)
	{
		//m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
		LeastSquaresGradientComputation grad3D(m_meshTet,
		                                       m_meshTet->getVariable<double, GMDS_NODE>("GMDS_Distance"),
		   												var_VectorsForExtrusion);
		grad3D.execute();
	}

	else if (m_params.vectors_field == 1)
	{
		//m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
		LeastSquaresGradientComputation grad3D(m_meshTet,
		                                       m_meshTet->getVariable<double, GMDS_NODE>("GMDS_Distance_Int"),
		   									var_VectorsForExtrusion);
		grad3D.execute();
	}

	else if (m_params.vectors_field == 2)
	{
		Variable<math::Vector3d>* var_VectorField_1 = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_1");
		Variable<math::Vector3d>* var_VectorField_2 = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_2");

		LeastSquaresGradientComputation grad2D_1(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
		                                         var_VectorField_1);
		grad2D_1.execute();
		LeastSquaresGradientComputation grad2D_2(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
		                                         var_VectorField_2);
		grad2D_2.execute();

		for (auto n_id:m_meshTet->nodes())
		{
			math::Vector3d vec_1 = var_VectorField_1->value(n_id);
			math::Vector3d vec_2 = var_VectorField_2->value(n_id);

			vec_1.normalize();
			vec_2.normalize();

			Node n = m_meshTet->get<Node>(n_id);
			math::Point p = n.point();

			if (p.X() <= m_params.x_VectorField_Z1)
			{
				var_VectorsForExtrusion->set(n_id, vec_1);
			}
			else if (p.X() >= m_params.x_VectorField_Z2)
			{
				var_VectorsForExtrusion->set(n_id, vec_2);
			}
			else
			{
				// Compute the transition field
				double alpha = (p.X() - m_params.x_VectorField_Z1)/(m_params.x_VectorField_Z2-m_params.x_VectorField_Z1) ;
				math::Vector3d v_transit = alpha*vec_2 + (1.0-alpha)*vec_1 ;
				v_transit.normalize();
				var_VectorsForExtrusion->set(n_id, v_transit);
			}

		}

		m_meshTet->deleteVariable(GMDS_NODE, var_VectorField_1);
		m_meshTet->deleteVariable(GMDS_NODE, var_VectorField_2);

	}

}
/*------------------------------------------------------------------------*/
void
AeroPipeline_3D::BlockingGeometricClassification()
{
	std::vector<cad::GeomSurface*> surfaces;
	m_manager->getSurfaces(surfaces);

	// Get the exterior surface
	double xyz_min[3];
	double xyz_max[3];
	cad::GeomSurface* surf_ext = surfaces[0];
	surf_ext->computeBoundingBox(xyz_min, xyz_max);
	for (auto const surf:surfaces)
	{
		double xyz_min_surf[3];
		double xyz_max_surf[3];
		surf->computeBoundingBox(xyz_min_surf, xyz_max_surf);
		if (xyz_min_surf[0] <= xyz_min[0]
		    && xyz_max[0] <= xyz_max_surf[0]
		    && xyz_min_surf[1] <= xyz_min[1]
		    && xyz_max[1] <= xyz_max_surf[1]
		    && xyz_min_surf[2] <= xyz_min[2]
		    && xyz_max[2] <= xyz_max_surf[2])
		{
			surf_ext = surf;
			xyz_min[0] = xyz_min_surf[0];
			xyz_min[1] = xyz_min_surf[1];
			xyz_min[2] = xyz_min_surf[2];
			xyz_max[0] = xyz_max_surf[0];
			xyz_max[1] = xyz_max_surf[1];
			xyz_max[2] = xyz_max_surf[2];
		}
	}

	//--------------------------------//
	//	  BLOCK NODES CLASSIFICATION   //
	//--------------------------------//
	for (auto const n_id:m_meshHex->nodes())
	{
		if (m_couche_id->value(n_id) == m_params.nbr_couches)
		{
			m_linker_HG->linkNodeToSurface(n_id, surf_ext->id());
		}
	}

	//--------------------------------//
	//	  BLOCK EDGES CLASSIFICATION   //
	//--------------------------------//
	for (auto const e_id:m_meshHex->edges())
	{
		std::vector<Node> e_nodes = m_meshHex->get<Edge>(e_id).get<Node>();
		if (m_couche_id->value(e_nodes[0].id()) == m_params.nbr_couches
		    && m_couche_id->value(e_nodes[1].id()) == m_params.nbr_couches)
		{
			m_linker_HG->linkEdgeToSurface(e_id, surf_ext->id());
		}
	}

	//--------------------------------//
	//	  BLOCK FACES CLASSIFICATION   //
	//--------------------------------//
	for (auto const f_id:m_meshHex->faces())
	{
		std::vector<Node> f_nodes = m_meshHex->get<Face>(f_id).get<Node>();
		if (m_couche_id->value(f_nodes[0].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[1].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[2].id()) == m_params.nbr_couches
		    && m_couche_id->value(f_nodes[3].id()) == m_params.nbr_couches)
		{
			m_linker_HG->linkFaceToSurface(f_id, surf_ext->id());
		}
	}

}
/*------------------------------------------------------------------------*/