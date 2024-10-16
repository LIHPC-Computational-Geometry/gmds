//
// Created by rochec on 27/08/24.
//

#include <gmds/aero/ComputeBezierDegree_3D.h>
#include <gmds/aero/ComputeBezierCtrlPtstoInterpolateSurface.h>
#include <gmds/aero/Blocking3D.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gtest/gtest.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
TEST(AeroTestClass, Test_ComputeBezierCtrlPtstoInterpolateSurface)
{
	gmds::Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R | F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	Node n1 = m.newNode(0, 0, 0);
	Node n2 = m.newNode(1, 0, 0);
	Node n3 = m.newNode(1, 1, 0);
	Node n4 = m.newNode(0, 1, 0);
	Node n5 = m.newNode(0, 0, 1);
	Node n6 = m.newNode(1, 0, 1);
	Node n7 = m.newNode(1, 1, 1);
	Node n8 = m.newNode(0, 1, 1);

	Region h1 = m.newHex(n1, n2, n3, n4, n5, n6, n7, n8);

	MeshDoctor doc_mesh(&m);
	doc_mesh.buildFacesAndR2F();
	doc_mesh.buildEdgesAndX2E();
	doc_mesh.updateUpwardConnectivity();

	gmds::Mesh m_surface = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                      F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	double amp(0.1);
	int nx(100);
	int ny(50);
	double dx(float(2)/float(nx-1));
	double dy(float(1)/float(ny-1));
	Array2D<math::Point> pts(nx,ny);
	for (int i=0;i<nx;i++)
	{
		for (int j=0;j<ny;j++)
		{
			pts(i,j).setX(i*dx);
			pts(i,j).setY(j*dy);
			pts(i,j).setZ(amp*cos(i*dx*4.0));
		}
	}

	Array2D<TCellID> n_ids(nx,ny);
	Array2D<TCellID> n_ids_2(nx,ny);
	for (int i=0;i<nx;i++)
	{
		for (int j=0;j<ny;j++)
		{
			Node n = m_surface.newNode(pts(i,j));
			n_ids(i,j) = n.id();
			n = m_surface.newNode(pts(i,j));
			n.setZ(pts(i,j).Z()+1.0);
			n_ids_2(i,j) = n.id();
		}
	}

	for (auto i=0;i<nx-1;i++)
	{
		for (int j=0;j<ny-1;j++)
		{
			//math::Utils::GetOrCreateQuadAndConnectivities(&m_surface, n_ids(i,j), n_ids(i+1,j), n_ids(i+1,j+1), n_ids(i,j+1));
			m_surface.newHex(n_ids(i,j), n_ids(i+1,j), n_ids(i+1,j+1), n_ids(i,j+1),
			                 n_ids_2(i,j), n_ids_2(i+1,j), n_ids_2(i+1,j+1), n_ids_2(i,j+1));
		}
	}

	MeshDoctor doc(&m_surface);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	cad::FACManager manager;
	manager.initFrom3DMesh(&m_surface);

	// Write the mesh that represents the surface
	IGMeshIOService ioService_surf(&m_surface);
	VTKWriter writer_surf(&ioService_surf);
	writer_surf.setCellOptions(N|F);
	writer_surf.setDataOptions(N|F);
	writer_surf.write("TEST_ComputeBezierCtrlPtstoInterpolateSurface_Surface.vtk");

	std::cout << "Nbr surfaces " << manager.getNbSurfaces() << std::endl;
	std::cout << "Nbr curves " << manager.getNbCurves() << std::endl;

	cad::GeomMeshLinker linker;
	linker.setGeometry(&manager);
	linker.setMesh(&m);

	ASSERT_EQ(manager.getNbSurfaces(), 6);
	ASSERT_EQ(manager.getNbCurves(), 12);

	TCellID f_1234 = math::Utils::CommonFace(&m, n1.id(), n2.id(), n3.id(), n4.id());
	std::vector<TCellID> faces_to_fit;
	faces_to_fit.push_back(f_1234);

	linker.linkFaceToSurface(f_1234, 1);

	cad::GeomSurface* surf = manager.getSurface(1);
	math::Point p = n1.point();
	surf->project(p);
	n1.setPoint(p);
	p = n2.point();
	surf->project(p);
	n2.setPoint(p);
	p = n3.point();
	surf->project(p);
	n3.setPoint(p);
	p = n4.point();
	surf->project(p);
	n4.setPoint(p);

	gmds::Mesh m_linearBlock = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                      F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	gmds::Node n1_linearBlock = m_linearBlock.newNode(n1.point());
	gmds::Node n2_linearBlock = m_linearBlock.newNode(n2.point());
	gmds::Node n3_linearBlock = m_linearBlock.newNode(n3.point());
	gmds::Node n4_linearBlock = m_linearBlock.newNode(n4.point());
	m_linearBlock.newQuad(n1_linearBlock, n2_linearBlock, n3_linearBlock, n4_linearBlock);

	// Write the linear block face
	IGMeshIOService ioService(&m_linearBlock);
	VTKWriter writer(&ioService);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("TEST_ComputeBezierCtrlPtstoInterpolateSurface_LinearBlock.vtk");

	Array2D<math::Point> InitCtrlPts(4,4);
	InitCtrlPts(0,0) = n1_linearBlock.point() ;
	InitCtrlPts(3,0) = n2_linearBlock.point() ;
	InitCtrlPts(3,3) = n3_linearBlock.point() ;
	InitCtrlPts(0,3) = n4_linearBlock.point() ;
	ComputeBezierCtrlPtstoInterpolateSurface algo = ComputeBezierCtrlPtstoInterpolateSurface(surf, &InitCtrlPts);
	ComputeBezierCtrlPtstoInterpolateSurface::STATUS res = algo.execute();

	ASSERT_EQ(res, ComputeBezierCtrlPtstoInterpolateSurface::SUCCESS);

	Array2D<math::Point> CtrlPts = algo.getCtrlPts();

	Blocking2D blocks = Blocking2D();
	gmds::Node bc_1 = blocks.newBlockCorner(CtrlPts(0,0));
	gmds::Node bc_2 = blocks.newBlockCorner(CtrlPts(3,0));
	gmds::Node bc_3 = blocks.newBlockCorner(CtrlPts(3,3));
	gmds::Node bc_4 = blocks.newBlockCorner(CtrlPts(0,3));

	Blocking2D::Block b = blocks.newBlock(bc_1,bc_2,bc_3,bc_4);
	b.setNbDiscretizationI(4);
	b.setNbDiscretizationJ(4);
	blocks.initializeGridPoints();

	b = blocks.block(0);

	// Write the projected ctrl pts positions on the surface
	gmds::Mesh m_initctrlpts = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                     F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	// Write the curved block face ctrl pts
	IGMeshIOService ioService_initctrlpts(&blocks);
	VTKWriter writer_initctrlpts(&ioService_initctrlpts);
	writer_initctrlpts.setCellOptions(N|F);
	writer_initctrlpts.setDataOptions(N|F);
	writer_initctrlpts.write("TEST_ComputeBezierCtrlPtstoInterpolateSurface_InitCtrlPts.vtk");

	for (auto i=0;i<b.getNbDiscretizationI();i++)
	{
		for (auto j=0;j<b.getNbDiscretizationJ();j++)
		{
			math::Point p = b(i,j).point();
			surf->project(p);
			b(i,j).setPoint(p) ;
		}
	}

	// Write the curved block face init ctrl pts projected onto geometric surface
	IGMeshIOService ioService_initctrlptsproj(&blocks);
	VTKWriter writer_initctrlptsproj(&ioService_initctrlptsproj);
	writer_initctrlptsproj.setCellOptions(N|F);
	writer_initctrlptsproj.setDataOptions(N|F);
	writer_initctrlptsproj.write("TEST_ComputeBezierCtrlPtstoInterpolateSurface_InitCtrlPtsProjected.vtk");


	for (auto i=0;i<b.getNbDiscretizationI();i++)
	{
		for (auto j=0;j<b.getNbDiscretizationJ();j++)
		{
			b(i,j).setPoint(CtrlPts(i,j)) ;
		}
	}

	// Write the curved block face ctrl pts
	IGMeshIOService ioService_ctrlpts(&blocks);
	VTKWriter writer_ctrlpts(&ioService_ctrlpts);
	writer_ctrlpts.setCellOptions(N|F);
	writer_ctrlpts.setDataOptions(N|F);
	writer_ctrlpts.write("TEST_ComputeBezierCtrlPtstoInterpolateSurface_CtrlPts.vtk");


	gmds::Mesh m_curved = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                          F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	// Create hex cells on the curved blocks
	int sample(20);
	math::BezierSurface bs(algo.getCtrlPts());
	for (int i=0;i<sample;i++)
	{
		for (int j=0;j<sample;j++)
		{
			math::Point p0 = bs(float(i)/float(sample), float(j)/float(sample));
			math::Point p1 = bs(float(i+1)/float(sample), float(j)/float(sample));
			math::Point p2 = bs(float(i+1)/float(sample), float(j+1)/float(sample));
			math::Point p3 = bs(float(i)/float(sample), float(j+1)/float(sample));
			Node n0 = m_curved.newNode(p0);
			Node n1 = m_curved.newNode(p1);
			Node n2 = m_curved.newNode(p2);
			Node n3 = m_curved.newNode(p3);
			m_curved.newQuad(n0,n1,n2,n3);
		}
	}

	// Write the curved block face ctrl pts
	IGMeshIOService ioService_curved(&m_curved);
	VTKWriter writer_curved(&ioService_curved);
	writer_curved.setCellOptions(N|F);
	writer_curved.setDataOptions(N|F);
	writer_curved.write("TEST_ComputeBezierCtrlPtstoInterpolateSurface_FinalMesh.vtk");


}
/*----------------------------------------------------------------------------*/
