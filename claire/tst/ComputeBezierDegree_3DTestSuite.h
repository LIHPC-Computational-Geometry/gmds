//
// Created by rochec on 17/01/24.
//
#include <gmds/claire/ComputeBezierDegree_3D.h>
#include <gmds/claire/ComputeBezierCtrlPtstoInterpolateSurface.h>
#include <gmds/claire/Blocking3D.h>
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
TEST(ClaireTestClass, Test_computeBezierDegree3D)
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

	Node n9 = m.newNode(2, 0, 0);
	Node n10 = m.newNode(2, 1, 0);
	Node n11 = m.newNode(2, 0, 1);
	Node n12 = m.newNode(2, 1, 1);

	Region h1 = m.newHex(n1, n2, n3, n4, n5, n6, n7, n8);
	Region h2 = m.newHex(n2, n3, n7, n6, n9, n10, n12, n11);

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
	writer_surf.write("TEST_ComputeBezierDegree3D_Surface.vtk");

	std::cout << "Nbr surfaces " << manager.getNbSurfaces() << std::endl;
	std::cout << "Nbr curves " << manager.getNbCurves() << std::endl;

	cad::GeomMeshLinker linker;
	linker.setGeometry(&manager);
	linker.setMesh(&m);

	ASSERT_EQ(manager.getNbSurfaces(), 6);
	ASSERT_EQ(manager.getNbCurves(), 12);

	TCellID f_1234 = math::Utils::CommonFace(&m, n1.id(), n2.id(), n3.id(), n4.id());
	TCellID f_23109 = math::Utils::CommonFace(&m, n2.id(), n3.id(), n10.id(), n9.id());
	std::vector<TCellID> faces_to_fit;
	faces_to_fit.push_back(f_1234);
	faces_to_fit.push_back(f_23109);

	linker.linkFaceToSurface(f_1234, 1);
	linker.linkFaceToSurface(f_23109, 1);

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
	p = n9.point();
	surf->project(p);
	n9.setPoint(p);
	p = n10.point();
	surf->project(p);
	n10.setPoint(p);

	// Write the mesh results
	IGMeshIOService ioService(&m);
	VTKWriter writer(&ioService);
	writer.setCellOptions(N|R);
	writer.setDataOptions(N|R);
	writer.write("TEST_ComputeBezierDegree3D_LinearBlocks.vtk");

	ComputeBezierDegree_3D algo = ComputeBezierDegree_3D(&m,&manager,&linker,&faces_to_fit,0.02,10);
	ComputeBezierDegree_3D::STATUS res = algo.execute();

	std::cout << "Final degree: " << algo.getDegree() << std::endl;

	ASSERT_EQ(res, ComputeBezierDegree_3D::SUCCESS);

	// Create the blocking
	Blocking3D blocking = Blocking3D();
	std::map<TCellID,TCellID> map_corresponding_ids;
	for (auto n_id:m.nodes())
	{
		Node bc = blocking.newBlockCorner(m.get<Node>(n_id).point());
		map_corresponding_ids[n_id] = bc.id();
	}
	for (auto r_id:m.regions())
	{
		std::vector<Node> r_nodes = m.get<Region>(r_id).get<Node>();
		Blocking3D::Block b = blocking.newBlock(map_corresponding_ids[r_nodes[0].id()],
		                                        map_corresponding_ids[r_nodes[1].id()],
		                                        map_corresponding_ids[r_nodes[2].id()],
		                                        map_corresponding_ids[r_nodes[3].id()],
		                                        map_corresponding_ids[r_nodes[4].id()],
		                                        map_corresponding_ids[r_nodes[5].id()],
		                                        map_corresponding_ids[r_nodes[6].id()],
		                                        map_corresponding_ids[r_nodes[7].id()]) ;
		b.setNbDiscretizationI(algo.getDegree()+1);
		b.setNbDiscretizationJ(algo.getDegree()+1);
		b.setNbDiscretizationK(algo.getDegree()+1);
	}
	blocking.initializeGridPoints();

	/*
	for (auto b:blocking.allBlocks())
	{
		for (auto i=0;i<b.getNbDiscretizationI();i++)
		{
			for (auto j=0;j<b.getNbDiscretizationJ();j++)
			{
				for (auto k=0;k<b.getNbDiscretizationK();k++)
				{
					math::Point p = b(i, j, k).point();
					surf->project(p);
					b(i,j,k).setPoint(p);
				}
			}
		}
	}
	 */

	gmds::Mesh m_curved = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                     F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	for (auto f_id:faces_to_fit)
	{
		Face f = m.get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();
		TCellID bf_id = math::Utils::CommonFace(&blocking, f_nodes[0].id(), f_nodes[1].id(), f_nodes[2].id(), f_nodes[3].id());
		Blocking3D::BlockFace bf = blocking.blockFace(bf_id);

		Array2D<math::Point> ctrl_pts(bf.getNbDiscretizationI(),bf.getNbDiscretizationJ()) ;
		for (int i=0;i<bf.getNbDiscretizationI();i++)
		{
			for (int j=0;j<bf.getNbDiscretizationJ();j++)
			{
				ctrl_pts(i,j) = bf(i,j).point();
			}
		}

		ComputeBezierCtrlPtstoInterpolateSurface algo_bezier = ComputeBezierCtrlPtstoInterpolateSurface(surf,
																																		&ctrl_pts);
		ComputeBezierCtrlPtstoInterpolateSurface::STATUS res_bezier = algo_bezier.execute();

		for (int i=0;i<bf.getNbDiscretizationI();i++)
		{
			for (int j=0;j<bf.getNbDiscretizationJ();j++)
			{
				ctrl_pts(i,j) = bf(i,j).point();
			}
		}

		ASSERT_EQ(res_bezier, ComputeBezierCtrlPtstoInterpolateSurface::SUCCESS);

		math::BezierSurface f_curved = math::BezierSurface(algo_bezier.getCtrlPts());

		for (int i=0;i<=algo.getDegree();i++)
		{
			for (int j=0;j<=algo.getDegree();j++)
			{
				bf(i,j).setPoint(algo_bezier.getCtrlPts()(i,j));
			}
		}

		int sample(50);
		for (int i=0;i<sample;i++)
		{
			for (int j=0;j<sample;j++)
			{
				p = f_curved(float(i)/float(sample-1), float(j)/float(sample-1));
				m_curved.newNode(p);
			}
		}

	}

	// Re compute positions of the control points on faces
	TInt mark_bndfaces = blocking.newMark<gmds::Face>();
	for (auto f_id:faces_to_fit)
	{
		blocking.mark(blocking.get<Face>(f_id), mark_bndfaces);
	}
	for (auto f_id:blocking.faces())
	{
		if (!blocking.isMarked(blocking.get<Face>(f_id), mark_bndfaces))
		{
			/*
			std::cout << "------------------------" << std::endl;
			std::cout << "Face " << f_id << std::endl;
			Blocking3D::BlockFace bf = blocking.blockFace(f_id);
			Array2D<math::Point> points(bf.getNbDiscretizationI(), bf.getNbDiscretizationJ());
			std::cout << "Nx " << bf.getNbDiscretizationI() << std::endl;
			std::cout << "Ny " << bf.getNbDiscretizationJ() << std::endl;
			for (int i=0;i<bf.getNbDiscretizationI();i++)
			{
				points(i,0) = bf(i,0).point() ;
				points(i,bf.getNbDiscretizationJ()-1) = bf(i,bf.getNbDiscretizationJ()-1).point() ;
			}
			for (int j=0;j<bf.getNbDiscretizationJ();j++)
			{
				points(0,j) = bf(0,j).point() ;
				points(bf.getNbDiscretizationI()-1,j) = bf(bf.getNbDiscretizationI()-1,j).point() ;
			}
			math::TransfiniteInterpolation::computeQuad(points);
			for (int i=1;i<bf.getNbDiscretizationI()-1;i++)
			{
				for (int j=1;j<bf.getNbDiscretizationJ()-1;j++)
				{
					bf(i,j).setPoint(points(i,j));
				}
			}
			*/
		}
	}

	for (auto b:blocking.allBlocks())
	{
		b.computeInnerBlockNodesPoints();
	}

	blocking.unmarkAll<Face>(mark_bndfaces);
	blocking.freeMark<Face>(mark_bndfaces);

	/*
	srand (time(NULL));
	for (auto b:blocking.allBlocks())
	{
		for (auto i=0;i<b.getNbDiscretizationI();i++)
		{
			for (auto j=0;j<b.getNbDiscretizationJ();j++)
			{
				for (auto k=0;k<b.getNbDiscretizationK();k++)
				{
					double dx = 0.2*((double) rand() / (RAND_MAX)) -0.1;
					double dy = 0.2*((double) rand() / (RAND_MAX)) -0.1;
					double dz = 0.2*((double) rand() / (RAND_MAX)) -0.1;
					//std::cout << "dx: " << dx << ", dy: " << dy << ", dz: " << dz << std::endl;
					b(i,j,k).setX(b(i,j,k).X()+dx);
					b(i,j,k).setY(b(i,j,k).Y()+dx);
					b(i,j,k).setZ(b(i,j,k).Z()+dx);
				}
			}
		}
	}
	gmds::Mesh m_tmp = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                     F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	                                     */

	// Create hex cells on the curved blocks
	int sample(20);
	for (auto b:blocking.allBlocks())
	{
		Array3D<math::Point> ctrl_pts(b.getNbDiscretizationI(), b.getNbDiscretizationJ(), b.getNbDiscretizationK());
		for (int i=0;i<b.getNbDiscretizationI();i++)
		{
			for (int j=0;j<b.getNbDiscretizationJ();j++)
			{
				for (int k=0;k<b.getNbDiscretizationK();k++)
				{
					ctrl_pts(i,j,k) = b(i,j,k).point();
				}
			}
		}
		math::BezierHex bh(ctrl_pts);
		for (int i=0;i<sample;i++)
		{
			for (int j=0;j<sample;j++)
			{
				for (int k=0;k<sample;k++)
				{
					math::Point p0 = bh(float(i)/float(sample), float(j)/float(sample), float(k)/float(sample));
					math::Point p1 = bh(float(i+1)/float(sample), float(j)/float(sample), float(k)/float(sample));
					math::Point p2 = bh(float(i+1)/float(sample), float(j+1)/float(sample), float(k)/float(sample));
					math::Point p3 = bh(float(i)/float(sample), float(j+1)/float(sample), float(k)/float(sample));
					math::Point p4 = bh(float(i)/float(sample), float(j)/float(sample), float(k+1)/float(sample));
					math::Point p5 = bh(float(i+1)/float(sample), float(j)/float(sample), float(k+1)/float(sample));
					math::Point p6 = bh(float(i+1)/float(sample), float(j+1)/float(sample), float(k+1)/float(sample));
					math::Point p7 = bh(float(i)/float(sample), float(j+1)/float(sample), float(k+1)/float(sample));
					Node n0 = m_curved.newNode(p0);
					Node n1 = m_curved.newNode(p1);
					Node n2 = m_curved.newNode(p2);
					Node n3 = m_curved.newNode(p3);
					Node n4 = m_curved.newNode(p4);
					Node n5 = m_curved.newNode(p5);
					Node n6 = m_curved.newNode(p6);
					Node n7 = m_curved.newNode(p7);
					m_curved.newHex(n0,n1,n2,n3,n4,n5,n6,n7);

					/*
					math::Point p0_tmp(float(i)/float(sample), float(j)/float(sample), float(k)/float(sample));
					math::Point p1_tmp(float(i+1)/float(sample), float(j)/float(sample), float(k)/float(sample));
					math::Point p2_tmp(float(i+1)/float(sample), float(j+1)/float(sample), float(k)/float(sample));
					math::Point p3_tmp(float(i)/float(sample), float(j+1)/float(sample), float(k)/float(sample));
					math::Point p4_tmp(float(i)/float(sample), float(j)/float(sample), float(k+1)/float(sample));
					math::Point p5_tmp(float(i+1)/float(sample), float(j)/float(sample), float(k+1)/float(sample));
					math::Point p6_tmp(float(i+1)/float(sample), float(j+1)/float(sample), float(k+1)/float(sample));
					math::Point p7_tmp(float(i)/float(sample), float(j+1)/float(sample), float(k+1)/float(sample));
					Node n0_tmp = m_tmp.newNode(p0_tmp);
					Node n1_tmp = m_tmp.newNode(p1_tmp);
					Node n2_tmp = m_tmp.newNode(p2_tmp);
					Node n3_tmp = m_tmp.newNode(p3_tmp);
					Node n4_tmp = m_tmp.newNode(p4_tmp);
					Node n5_tmp = m_tmp.newNode(p5_tmp);
					Node n6_tmp = m_tmp.newNode(p6_tmp);
					Node n7_tmp = m_tmp.newNode(p7_tmp);
					m_tmp.newHex(n0_tmp,n1_tmp,n2_tmp,n3_tmp,n4_tmp,n5_tmp,n6_tmp,n7_tmp);
					 */
				}
			}
		}
	}


	// Write the mesh results
	IGMeshIOService ioService_curved(&m_curved);
	//IGMeshIOService ioService_curved(&m_curved);
	VTKWriter writer_curved(&ioService_curved);
	writer_curved.setCellOptions(N|R);
	writer_curved.setDataOptions(N|R);
	writer_curved.write("TEST_ComputeBezierDegree3D_FinalMesh.vtk");

	// Write the control points of the blocking
	IGMeshIOService ioService_blocking_ctrlpts(&blocking);
	VTKWriter writer_blocking_ctrlpts(&ioService_blocking_ctrlpts);
	writer_blocking_ctrlpts.setCellOptions(N|R);
	writer_blocking_ctrlpts.setDataOptions(N|R);
	writer_blocking_ctrlpts.write("TEST_ComputeBezierDegree3D_CtrlPts.vtk");


}
/*----------------------------------------------------------------------------*/
