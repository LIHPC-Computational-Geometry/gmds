//
// Created by rochec on 13/02/24.
//

#include <gmds/aero/ComputeBezierCurveCtrlPtstoInterpolateCurve.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/math/BezierCurve.h>
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
TEST(AeroTestClass, Test_ComputeBezierCurveCtrlPtstoInterpolateCurve)
{
	//---> Build the Curve 1
	int Nx(30);
	int Nx2(60);
	double rayon(0.15);
	double dr(0.1);
	double ymax(0.15);
	double a = (rayon-ymax)/(rayon-1.0);
	double b = ymax-a;

	std::vector<math::Point> curve_discrete(Nx+Nx2);
	std::vector<math::Point> curve_discrete_opp(Nx+Nx2);
	for (int i=0;i<Nx;i++)
	{
		double theta = (float(i)/float(Nx-1.0))*M_PI/2.0 ;
		curve_discrete[i].setX(rayon-rayon*cos(theta));
		curve_discrete[i].setY(rayon*sin(theta));
		curve_discrete[i].setZ(0.0);

		curve_discrete_opp[i].setX(rayon-(dr+rayon)*cos(theta));
		curve_discrete_opp[i].setY((dr+rayon)*sin(theta));
		curve_discrete_opp[i].setZ(0.0);
	}
	for (int i=Nx;i<curve_discrete.size();i++)
	{
		curve_discrete[i].setX(rayon+(1.0-rayon)*float(i-Nx)/float(Nx2-1.0));
		curve_discrete[i].setY(a*curve_discrete[i].X()+b);
		curve_discrete[i].setZ(0);

		curve_discrete_opp[i].setX(rayon+(1.0-rayon)*float(i-Nx)/float(Nx2-1.0));
		curve_discrete_opp[i].setY(dr+a*curve_discrete[i].X()+b);
		curve_discrete_opp[i].setZ(0);
	}

	gmds::Mesh m_curve = Mesh(MeshModel(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E));
	std::map<int,TCellID> curve_nodes;
	std::map<int,TCellID> curve_nodes_opp;
	for (int i=0;i<curve_discrete.size();i++)
	{
		Node n = m_curve.newNode(curve_discrete[i]);
		curve_nodes[i] = n.id() ;
		n = m_curve.newNode(curve_discrete_opp[i]);
		curve_nodes_opp[i] = n.id() ;
	}
	for (int i=0;i<curve_discrete.size()-1;i++)
	{
		Face f = m_curve.newQuad(curve_nodes[i], curve_nodes[i+1], curve_nodes_opp[i+1], curve_nodes_opp[i]);
	}
	// <--- Curve 1 built

	MeshDoctor doc(&m_curve);
	doc.buildEfromF();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	cad::FACManager manager;
	cad::GeomMeshLinker linker;
	manager.initAndLinkFrom2DMesh(&m_curve, &linker);

	std::cout << "Nbr curves " << manager.getNbCurves() << std::endl;

	// Write the mesh that represents the curve
	IGMeshIOService ioService_curve(&m_curve);
	VTKWriter writer_curve(&ioService_curve);
	writer_curve.setCellOptions(N|E);
	writer_curve.setDataOptions(N|E);
	writer_curve.write("TEST_ComputeBezierCurveCtrlPtstoInterpolateCurve_Curve.vtk");

	// Init the control points of the curve
	/*
	int degree(1);
	std::vector<math::Point> ctrl_pts(degree+1);
	ctrl_pts[0].setX(-dr);
	ctrl_pts[0].setY(0.0);
	ctrl_pts[0].setZ(0.0);
	ctrl_pts[degree].setX(1.0);
	ctrl_pts[degree].setY(ymax+dr);
	ctrl_pts[degree].setZ(0.0);
	for (int i=1;i<degree;i++)
	{
		ctrl_pts[i].setX(ctrl_pts[0].X() + float(i)*(ctrl_pts[degree].X()-ctrl_pts[0].X())/float(degree));
		ctrl_pts[i].setY(ctrl_pts[0].Y() + float(i)*(ctrl_pts[degree].Y()-ctrl_pts[0].Y())/float(degree));
		ctrl_pts[i].setZ(ctrl_pts[0].Z() + float(i)*(ctrl_pts[degree].Z()-ctrl_pts[0].Z())/float(degree));
	}

	cad::GeomCurve* geom_curve = manager.getCurve(2);

	ComputeBezierCurveCtrlPtstoInterpolateCurve algo = ComputeBezierCurveCtrlPtstoInterpolateCurve(geom_curve, &ctrl_pts);
	ComputeBezierCurveCtrlPtstoInterpolateCurve::STATUS res_algo = algo.execute();
	ASSERT_EQ(res_algo, ComputeBezierCurveCtrlPtstoInterpolateCurve::SUCCESS);
	*/
	cad::GeomCurve* geom_curve = manager.getCurve(2);
	int degree(6);
	//int degree_loc(1);
	std::vector<math::Point> ctrl_pts(2);
	ctrl_pts[0].setX(-dr);
	ctrl_pts[0].setY(0.0);
	ctrl_pts[0].setZ(0.0);
	ctrl_pts[1].setX(1.0);
	ctrl_pts[1].setY(ymax+dr);
	ctrl_pts[1].setZ(0.0);
	std::map<int,std::vector<math::Point>> iterativ_ctrl_pts;
	iterativ_ctrl_pts[1] = ctrl_pts;
	std::vector<math::Point> ctrl_pts_2(3);
	ctrl_pts_2[0].setX(-dr);
	ctrl_pts_2[0].setY(0.0);
	ctrl_pts_2[0].setZ(0.0);
	ctrl_pts_2[1].setX(-0.1);
	ctrl_pts_2[1].setY(0.25);
	ctrl_pts_2[1].setZ(0.0);
	ctrl_pts_2[2].setX(1.0);
	ctrl_pts_2[2].setY(ymax+dr);
	ctrl_pts_2[2].setZ(0.0);
	iterativ_ctrl_pts[2] = ctrl_pts_2;

	//-----------//
	// METHODE 4 //
	//-----------//
	/*
	std::vector<math::Point> ctrl_pts_p2(3);
	math::BezierCurve loc_bc = math::BezierCurve(iterativ_ctrl_pts[2]);
	ctrl_pts_p2[0] = iterativ_ctrl_pts[2][0] ;
	ctrl_pts_p2[1] = loc_bc(0.5);
	ctrl_pts_p2[2] = iterativ_ctrl_pts[2][2] ;
	ComputeBezierCurveCtrlPtstoInterpolateCurve loc_algo = ComputeBezierCurveCtrlPtstoInterpolateCurve(geom_curve, &ctrl_pts_p2);
	ComputeBezierCurveCtrlPtstoInterpolateCurve::STATUS res_loc_algo = loc_algo.execute();
	iterativ_ctrl_pts[2] = loc_algo.getCtrlPts();
	std::cout << "New pos: " << loc_algo.getCtrlPts()[1] << std::endl;
	 */

	for (int loc_degree=3;loc_degree<=degree;loc_degree++)
	{
		math::BezierCurve loc_bc = math::BezierCurve(iterativ_ctrl_pts[loc_degree-1]);
		std::cout << "Bezier Curve Degree " << loc_degree-1 << ", Max Error: "
		          << math::Utils::maxErrorBtwBezierCurveandGeomCurve(&loc_bc, geom_curve, 100) << std::endl;
		std::vector<math::Point> loc_ctrl_pts(loc_degree+1);
		for (int i=0;i<loc_degree+1;i++)
		{
			loc_ctrl_pts[i] = loc_bc(float(i)/float(loc_degree)) ;
		}
		ComputeBezierCurveCtrlPtstoInterpolateCurve loc_algo = ComputeBezierCurveCtrlPtstoInterpolateCurve(geom_curve, &loc_ctrl_pts);
		ComputeBezierCurveCtrlPtstoInterpolateCurve::STATUS res_loc_algo = loc_algo.execute();
		ASSERT_EQ(res_loc_algo, ComputeBezierCurveCtrlPtstoInterpolateCurve::SUCCESS);
		loc_ctrl_pts = loc_algo.getCtrlPts();
		iterativ_ctrl_pts[loc_degree] = loc_ctrl_pts;
	}


	// Plot the ctrl pts
	gmds::Mesh m_ctrl_pts = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R | F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	for (int i=0;i<=degree;i++)
	{
		//m_ctrl_pts.newNode(algo.getCtrlPts()[i]);
		m_ctrl_pts.newNode(iterativ_ctrl_pts[degree][i]);
	}

	// Write the control points
	IGMeshIOService ioService_ctrlpts(&m_ctrl_pts);
	VTKWriter writer_ctrlpts(&ioService_ctrlpts);
	writer_ctrlpts.setCellOptions(N|E);
	writer_ctrlpts.setDataOptions(N|E);
	writer_ctrlpts.write("TEST_ComputeBezierCurveCtrlPtstoInterpolateCurve_CtrlPts.vtk");

	// Build the Bezier Curve
	//math::BezierCurve bc = math::BezierCurve(algo.getCtrlPts());
	math::BezierCurve bc = math::BezierCurve(iterativ_ctrl_pts[degree]);

	std::cout << "Bezier Curve Degree " << degree << ", Max Error: "
	          << math::Utils::maxErrorBtwBezierCurveandGeomCurve(&bc, geom_curve, 100) << std::endl;

	int sample(100);
	std::vector<math::Point> bezier_curve_nodes(sample);
	std::map<int,TCellID> bc_nodes;
	gmds::Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R | F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	for (int i=0;i<sample;i++)
	{
		double u = float(i)/float(sample-1);
		bezier_curve_nodes[i] = bc(u) ;
		Node n = m.newNode(bezier_curve_nodes[i]);
		bc_nodes[i] = n.id();
	}
	for (int i=0;i<sample-1;i++)
	{
		m.newEdge(bc_nodes[i], bc_nodes[i+1]);
	}
	// Write the control points
	IGMeshIOService ioService_m(&m);
	VTKWriter writer_m(&ioService_m);
	writer_m.setCellOptions(N|E);
	writer_m.setDataOptions(N|E);
	writer_m.write("TEST_ComputeBezierCurveCtrlPtstoInterpolateCurve_FinalMesh.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(AeroTestClass, Test_ComputeBezierCurveCtrlPtstoInterpolateCurve_2)
{
	//---> Build the CURVE 2
	int Nx(60);
	//double ymax(0.4); // Curves 2 et 3
	double ymax(0.2); // Curves 4

	std::vector<math::Point> curve_discrete(Nx);
	std::vector<math::Point> curve_discrete_opp(Nx);
	for (int i=0;i<Nx;i++)
	{
		//double theta = M_PI / 2.0;	// Curve 2
		//double theta = M_PI; // Curve 3
		double theta = 2.0*M_PI; // Curve 4
		curve_discrete[i].setX(float(i)/float(Nx-1));
		curve_discrete[i].setY(ymax*sin(theta*float(i)/float(Nx-1)));
		curve_discrete[i].setZ(0.0);

		curve_discrete_opp[i].setX(float(i)/float(Nx-1));
		curve_discrete_opp[i].setY(ymax+ymax*sin(theta*float(i)/float(Nx-1)));
		curve_discrete_opp[i].setZ(0.0);
	}

	gmds::Mesh m_curve = Mesh(MeshModel(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E));
	std::map<int,TCellID> curve_nodes;
	std::map<int,TCellID> curve_nodes_opp;
	for (int i=0;i<curve_discrete.size();i++)
	{
		Node n = m_curve.newNode(curve_discrete[i]);
		curve_nodes[i] = n.id() ;
		n = m_curve.newNode(curve_discrete_opp[i]);
		curve_nodes_opp[i] = n.id() ;
	}
	for (int i=0;i<curve_discrete.size()-1;i++)
	{
		Face f = m_curve.newQuad(curve_nodes[i], curve_nodes[i+1], curve_nodes_opp[i+1], curve_nodes_opp[i]);
	}
	// <--- Curve 2 built

	MeshDoctor doc(&m_curve);
	doc.buildEfromF();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	cad::FACManager manager;
	cad::GeomMeshLinker linker;
	manager.initAndLinkFrom2DMesh(&m_curve, &linker);

	std::cout << "Nbr curves " << manager.getNbCurves() << std::endl;

	// Write the mesh that represents the curve
	IGMeshIOService ioService_curve(&m_curve);
	VTKWriter writer_curve(&ioService_curve);
	writer_curve.setCellOptions(N|E);
	writer_curve.setDataOptions(N|E);
	writer_curve.write("TEST_ComputeBezierCurveCtrlPtstoInterpolateCurve_Curve.vtk");

	// Init the control points of the curve
	//----------//
	// METHOD 1 //
	//----------//
	/*
	int degree(5);
	std::vector<math::Point> ctrl_pts(degree+1);
	ctrl_pts[0].setX(0);
	ctrl_pts[0].setY(0.0);
	ctrl_pts[0].setZ(0.0);
	ctrl_pts[degree].setX(1.0);
	ctrl_pts[degree].setY(0.0);
	ctrl_pts[degree].setZ(0.0);
	for (int i=1;i<degree;i++)
	{
	   ctrl_pts[i].setX(ctrl_pts[0].X() + float(i)*(ctrl_pts[degree].X()-ctrl_pts[0].X())/float(degree));
	   ctrl_pts[i].setY(ctrl_pts[0].Y() + float(i)*(ctrl_pts[degree].Y()-ctrl_pts[0].Y())/float(degree));
	   ctrl_pts[i].setZ(ctrl_pts[0].Z() + float(i)*(ctrl_pts[degree].Z()-ctrl_pts[0].Z())/float(degree));
	}

	cad::GeomCurve* geom_curve = manager.getCurve(1);

	ComputeBezierCurveCtrlPtstoInterpolateCurve algo = ComputeBezierCurveCtrlPtstoInterpolateCurve(geom_curve, &ctrl_pts);
	ComputeBezierCurveCtrlPtstoInterpolateCurve::STATUS res_algo = algo.execute();
	ASSERT_EQ(res_algo, ComputeBezierCurveCtrlPtstoInterpolateCurve::SUCCESS);
	*/

	//----------//
	// METHOD 2 //
	//----------//
	cad::GeomCurve* geom_curve = manager.getCurve(1);

	math::Vector3d t0 = geom_curve->computeTangent(0);
	math::Vector3d t1 = geom_curve->computeTangent(1);

	math::Point p0(0.0,0.0,0.0) ;
	math::Point p1(1.0,0.0,0.0) ;
	math::Line l0 = math::Line(p0, t0);
	math::Line l1 = math::Line(p1, t1);

	math::Point pi ;
	double param;
	std::cout << l0.intersect2D(l1, pi, param) << std::endl;

	int degree(5);
	//int degree_loc(1);
	std::vector<math::Point> ctrl_pts(2);
	ctrl_pts[0].setX(0.0);
	ctrl_pts[0].setY(0.0);
	ctrl_pts[0].setZ(0.0);
	ctrl_pts[1].setX(1.0);
	ctrl_pts[1].setY(0.0);
	ctrl_pts[1].setZ(0.0);
	std::map<int,std::vector<math::Point>> iterativ_ctrl_pts;
	iterativ_ctrl_pts[1] = ctrl_pts;
	std::vector<math::Point> ctrl_pts_2(3);
	ctrl_pts_2[0].setX(ctrl_pts[0].X());
	ctrl_pts_2[0].setY(ctrl_pts[0].Y());
	ctrl_pts_2[0].setZ(ctrl_pts[0].Z());
	ctrl_pts_2[1].setX(pi.X());
	ctrl_pts_2[1].setY(pi.Y());
	ctrl_pts_2[1].setZ(0.0);
	//std::cout << ctrl_pts_2[1] << std::endl;
	//ctrl_pts_2[1].setX(0.5);
	//ctrl_pts_2[1].setY(0.6);
	//ctrl_pts_2[1].setZ(0.0);
	ctrl_pts_2[2].setX(ctrl_pts[1].X());
	ctrl_pts_2[2].setY(ctrl_pts[1].Y());
	ctrl_pts_2[2].setZ(ctrl_pts[1].Z());
	iterativ_ctrl_pts[2] = ctrl_pts_2;

	//-----------//
	// METHODE 4 //
	//-----------//
	std::vector<math::Point> ctrl_pts_p2(3);
	math::BezierCurve loc_bc = math::BezierCurve(iterativ_ctrl_pts[2]);
	ctrl_pts_p2[0] = iterativ_ctrl_pts[2][0] ;
	ctrl_pts_p2[1] = loc_bc(0.5);
	ctrl_pts_p2[2] = iterativ_ctrl_pts[2][2] ;
	ComputeBezierCurveCtrlPtstoInterpolateCurve loc_algo = ComputeBezierCurveCtrlPtstoInterpolateCurve(geom_curve, &ctrl_pts_p2);
	ComputeBezierCurveCtrlPtstoInterpolateCurve::STATUS res_loc_algo = loc_algo.execute();
	iterativ_ctrl_pts[2] = loc_algo.getCtrlPts();
	std::cout << "New pos: " << loc_algo.getCtrlPts()[1] << std::endl;

	for (int loc_degree=3;loc_degree<=degree;loc_degree++)
	{
		math::BezierCurve loc_bc = math::BezierCurve(iterativ_ctrl_pts[loc_degree-1]);
		std::cout << "Bezier Curve Degree " << loc_degree-1 << ", Max Error: "
		          << math::Utils::maxErrorBtwBezierCurveandGeomCurve(&loc_bc, geom_curve, 100) << std::endl;
		std::vector<math::Point> loc_ctrl_pts(loc_degree+1);
		for (int i=0;i<loc_degree+1;i++)
		{
			loc_ctrl_pts[i] = loc_bc(float(i)/float(loc_degree)) ;
		}
		ComputeBezierCurveCtrlPtstoInterpolateCurve loc_algo = ComputeBezierCurveCtrlPtstoInterpolateCurve(geom_curve, &loc_ctrl_pts);
		ComputeBezierCurveCtrlPtstoInterpolateCurve::STATUS res_loc_algo = loc_algo.execute();
		ASSERT_EQ(res_loc_algo, ComputeBezierCurveCtrlPtstoInterpolateCurve::SUCCESS);
		loc_ctrl_pts = loc_algo.getCtrlPts();
		iterativ_ctrl_pts[loc_degree] = loc_ctrl_pts;
	}




	// Plot the ctrl pts
	gmds::Mesh m_ctrl_pts = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R | F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	for (int i=0;i<=degree;i++)
	{
		//m_ctrl_pts.newNode(algo.getCtrlPts()[i]);
		m_ctrl_pts.newNode(iterativ_ctrl_pts[degree][i]);
	}

	// Write the control points
	IGMeshIOService ioService_ctrlpts(&m_ctrl_pts);
	VTKWriter writer_ctrlpts(&ioService_ctrlpts);
	writer_ctrlpts.setCellOptions(N|E);
	writer_ctrlpts.setDataOptions(N|E);
	writer_ctrlpts.write("TEST_ComputeBezierCurveCtrlPtstoInterpolateCurve_CtrlPts.vtk");

	// Build the Bezier Curve
	//math::BezierCurve bc = math::BezierCurve(algo.getCtrlPts());
	math::BezierCurve bc = math::BezierCurve(iterativ_ctrl_pts[degree]);

	std::cout << "Bezier Curve Degree " << degree << ", Max Error: "
	          << math::Utils::maxErrorBtwBezierCurveandGeomCurve(&bc, geom_curve, 100) << std::endl;

	int sample(100);
	std::vector<math::Point> bezier_curve_nodes(sample);
	std::map<int,TCellID> bc_nodes;
	gmds::Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R | F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	for (int i=0;i<sample;i++)
	{
		double u = float(i)/float(sample-1);
		bezier_curve_nodes[i] = bc(u) ;
		Node n = m.newNode(bezier_curve_nodes[i]);
		bc_nodes[i] = n.id();
	}
	for (int i=0;i<sample-1;i++)
	{
		m.newEdge(bc_nodes[i], bc_nodes[i+1]);
	}
	// Write the control points
	IGMeshIOService ioService_m(&m);
	VTKWriter writer_m(&ioService_m);
	writer_m.setCellOptions(N|E);
	writer_m.setDataOptions(N|E);
	writer_m.write("TEST_ComputeBezierCurveCtrlPtstoInterpolateCurve_FinalMesh.vtk");

}
/*----------------------------------------------------------------------------*/