//
// Created by calderans on 06/07/22.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blockMesher/BlockMesher2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/cad/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include "gmds/blockMesher/BlockClassificator.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(BlockMesher2DTestSuite, test_1_quart_disk)
{
	Mesh triMesh(MeshModel(gmds::DIM2|gmds::F|gmds::E|gmds::N|gmds::N2E|gmds::E2N|gmds::E2F|gmds::N2F|gmds::F2N));

	std::string vtk_file = "/home/calderans/dev/v5.vtk";

	gmds::IGMeshIOService ioService(&triMesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&triMesh);
	doc.updateUpwardConnectivity();

	Mesh blockMesh(MeshModel(DIM2|F|E|N|N2E|E2N|E2F|N2F|F2N));

	cad::FACManager manager;
	cad::GeomMeshLinker linkerTri;
	manager.initAndLinkFrom2DMesh(&triMesh, &linkerTri);

	Node n0 = blockMesh.newNode(0,0,0);
	Node n1 = blockMesh.newNode(1,0,0);
	Node n2 = blockMesh.newNode(1,1,0);
	Node n3 = blockMesh.newNode(0,1,0);

	Node n4 = blockMesh.newNode(0.5,0,0);
	Node n5 = blockMesh.newNode(0,0.5,0);
	Node n6 = blockMesh.newNode(0.5,0.5,0);

	Face f0 = blockMesh.newQuad(n0,n4,n6,n5);
	Face f1 = blockMesh.newQuad(n4,n1,n2,n6);
	Face f2 = blockMesh.newQuad(n5,n6,n2,n3);

	Edge e0 = blockMesh.newEdge(n0,n4);
	e0.add<Face>(f0);
	Edge e1 = blockMesh.newEdge(n4,n1);
	e1.add<Face>(f1);
	Edge e2 = blockMesh.newEdge(n1,n2);
	e2.add<Face>(f1);
	Edge e3 = blockMesh.newEdge(n2,n3);
	e3.add<Face>(f2);
	Edge e4 = blockMesh.newEdge(n3,n5);
	e4.add<Face>(f2);
	Edge e5 = blockMesh.newEdge(n5,n0);
	e5.add<Face>(f0);
	Edge e6 = blockMesh.newEdge(n4,n6);
	e6.add<Face>(f0);
	e6.add<Face>(f1);
	Edge e7 = blockMesh.newEdge(n5,n6);
	e7.add<Face>(f0);
	e7.add<Face>(f2);
	Edge e8 = blockMesh.newEdge(n2,n6);
	e8.add<Face>(f1);
	e8.add<Face>(f2);

	gmds::MeshDoctor doc2(&blockMesh);
	doc2.updateUpwardConnectivity();

	cad::GeomMeshLinker linker;
	linker.setGeometry(&manager);
	linker.setMesh(&blockMesh);

	linker.linkNodeToPoint(n0.id(),1);
	linker.linkNodeToPoint(n1.id(),2);
	linker.linkNodeToPoint(n3.id(),3);
	linker.linkNodeToCurve(n2.id(),3);

	BlockMesher2D bm(&blockMesh,&linker,&manager);

	bm.updateClassification();

	bm.projectNodes();

	bm.executeMeshing();

	std::cout<<"nb faces "<<n6.nbFaces()<<std::endl;
	std::cout<<"nb edges "<<n6.nbEdges()<<std::endl;

	std::cout<<"nb edges mesh "<<blockMesh.getNbEdges()<<std::endl;

	std::cout<<"nb points "<<manager.getNbPoints()<<std::endl;
	std::cout<<"nb curves "<<manager.getNbCurves()<<std::endl;

	IGMeshIOService ios(&triMesh);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("/home/calderans/dev/test_mesh.vtk");

	IGMeshIOService iosb(&blockMesh);
	VTKWriter writerb(&iosb);
	writerb.setCellOptions(N|F);
	writerb.setDataOptions(N|F);
	writerb.write("/home/calderans/dev/test_blocks.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(BlockMesher2DTestSuite, test_2_anchor)
{
	Mesh triMesh(MeshModel(gmds::DIM2|gmds::F|gmds::E|gmds::N|gmds::N2E|gmds::E2N|gmds::E2F|gmds::N2F|gmds::F2N|F2E));

	std::string vtk_file_geom = "/home/calderans/dev/v3_old.vtk";

	gmds::IGMeshIOService ioService(&triMesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReader.read(vtk_file_geom);

	gmds::MeshDoctor doc(&triMesh);
	doc.updateUpwardConnectivity();

	Mesh blockMesh(MeshModel(DIM2|F|E|N|N2E|E2N|E2F|N2F|F2N|F2E));

	cad::FACManager manager;
	cad::GeomMeshLinker linkerTri;
	manager.initAndLinkFrom2DMesh(&triMesh, &linkerTri);


	std::string vtk_file_blocks = "/home/calderans/dev/exec025_action17_on_1_b.vtk";

	gmds::IGMeshIOService ioServiceb(&blockMesh);
	gmds::VTKReader vtkReaderb(&ioServiceb);
	vtkReaderb.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReaderb.read(vtk_file_blocks);

	gmds::MeshDoctor doc2(&blockMesh);
	doc2.buildE();
	doc2.buildF2E(MeshModel(DIM2|F|E|N|N2E|E2N|E2F|N2F|F2N|F2E));
	doc2.updateUpwardConnectivity();


	for(auto e : blockMesh.edges()){
		Edge edge = blockMesh.get<Edge>(e);
		Node n0 = edge.get<Node>()[0];
		Node n1 = edge.get<Node>()[1];

		std::vector<TCellID> n0_faces = n0.getIDs<Face>();
		std::vector<TCellID> n1_faces = n1.getIDs<Face>();

		std::set<TCellID> common_faces;

		for(auto f0 : n0_faces){
			for(auto f1 : n1_faces){
				if(f0 == f1){
					common_faces.insert(f0);
				}
			}
		}
		for(auto f : common_faces){
			//edge.add<Face>(f);
			Face face = blockMesh.get<Face>(f);
			face.add<Edge>(e);
		}
	}

	std::cout<<"nb edges mesh "<<blockMesh.getNbEdges()<<std::endl;
	std::cout<<"nb faces mesh "<<blockMesh.getNbFaces()<<std::endl;



	/*for(auto e : blockMesh.edges()){
		Edge edge = blockMesh.get<Edge>(e);
		Node n0 = edge.get<Node>()[0];
		Node n1 = edge.get<Node>()[1];
		n0.add<Edge>(e);
		n1.add<Edge>(e);
	}*/

	cad::GeomMeshLinker linker;
	linker.setGeometry(&manager);
	linker.setMesh(&blockMesh);

	linker.linkNodeToPoint(2,4);
	linker.linkNodeToPoint(3,3);
	linker.linkNodeToPoint(5,1);
	linker.linkNodeToPoint(7,2);
	for(auto i = 8; i<=17;i++){
		linker.linkNodeToCurve(i,1);
	}
	linker.linkNodeToCurve(0,3);
	linker.linkNodeToCurve(1,3);
	linker.linkNodeToCurve(4,3);
	linker.linkNodeToCurve(6,3);

	BlockMesher2D bm(&blockMesh,&linker,&manager);

	bm.updateClassification();

	bm.projectNodes();

	//std::cout<<"nb faces "<<n6.nbFaces()<<std::endl;
	//std::cout<<"nb edges "<<n6.nbEdges()<<std::endl;


	std::cout<<"nb points "<<manager.getNbPoints()<<std::endl;
	std::cout<<"nb curves "<<manager.getNbCurves()<<std::endl;
	//std::cout<<"nb surfaces "<<manager.getNbSurfaces()<<std::endl;

	for (int i = 0; i < 1; ++i) {
		for (auto n : blockMesh.nodes()) {/*
			if (linker.getGeomDim<Node>(n) == 2) {
				Node node = blockMesh.get<Node>(n);
				std::vector<math::Point> voisins;
				voisins.push_back(node.point());
				std::cout<<"nb edges "<<node.get<Edge>().size()<<std::endl;
				for (auto e_n : node.get<Edge>()) {
					Node voisin = e_n.get<Node>()[0] != node ? e_n.get<Node>()[0] : e_n.get<Node>()[1];
					if (linker.getGeomDim<Node>(voisin.id()) > 0) {
						voisins.push_back(voisin.point());
					}
				}
				std::cout<<"point avant "<<node.point()<<std::endl;
				math::Point center(math::Point::massCenter(voisins));
				cad::GeomCurve *c = manager.getCurve(linker.getGeomId<Node>(n));
				std::cout<<"point avant proj "<<center<<std::endl;

				c->project(center);
				std::cout<<"point après proj "<<center<<std::endl;
				node.setPoint(center);
				std::cout<<"point du noeud apres "<<node.point()<<std::endl;

			}*/
			/*else */if (linker.getGeomDim<Node>(n) == 0) {
				Node node = blockMesh.get<Node>(n);
				std::vector<math::Point> voisins;
				voisins.push_back(node.point());
				for (auto const &e_n : node.get<Edge>()) {
					Node voisin = e_n.get<Node>()[0] != node ? e_n.get<Node>()[0] : e_n.get<Node>()[1];
						voisins.push_back(voisin.point());
				}
				math::Point center(math::Point::massCenter(voisins));
				node.setPoint(center);
			}
		}
	}
	//smoother.smoothSurfaces(10);

	bm.executeMeshing();

	IGMeshIOService ios(&triMesh);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("/home/calderans/dev/test_mesh.vtk");

	IGMeshIOService iosb(&blockMesh);
	VTKWriter writerb(&iosb);
	writerb.setCellOptions(N|F);
	writerb.setDataOptions(N|F);
	writerb.write("/home/calderans/dev/test_blocks.vtk");

Mesh* quad(bm.getMesh());

	for(auto n : quad->nodes()){
		//quad->newTriangle(n,n,n);
	}

	IGMeshIOService iosQ(quad);
	VTKWriter writerQ(&iosQ);
	writerQ.setCellOptions(N|F);
	writerQ.setDataOptions(N|F);
	writerQ.write("/home/calderans/dev/test_quad.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(BlockMesher2DTestSuite, test_3_anchor2)
{
	Mesh triMesh(MeshModel(gmds::DIM2|gmds::F|gmds::E|gmds::N|gmds::N2E|gmds::E2N|gmds::E2F|gmds::N2F|gmds::F2N|F2E));

	std::string vtk_file_geom = "/home/calderans/dev/v6.vtk";

	gmds::IGMeshIOService ioService(&triMesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReader.read(vtk_file_geom);

	gmds::MeshDoctor doc(&triMesh);
	doc.updateUpwardConnectivity();

	Mesh blockMesh(MeshModel(DIM2|F|E|N|N2E|E2N|E2F|N2F|F2N|F2E));

	cad::FACManager manager;
	cad::GeomMeshLinker linkerTri;
	manager.initAndLinkFrom2DMesh(&triMesh, &linkerTri);


	std::string vtk_file_blocks = "/home/calderans/dev/exec021_action17_on_1_b.vtk";

	gmds::IGMeshIOService ioServiceb(&blockMesh);
	gmds::VTKReader vtkReaderb(&ioServiceb);
	vtkReaderb.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReaderb.read(vtk_file_blocks);

	gmds::MeshDoctor doc2(&blockMesh);
	doc2.buildE();
	doc2.buildF2E(MeshModel(DIM2|F|E|N|N2E|E2N|E2F|N2F|F2N|F2E));
	doc2.updateUpwardConnectivity();


	for(auto e : blockMesh.edges()){
		Edge edge = blockMesh.get<Edge>(e);
		Node n0 = edge.get<Node>()[0];
		Node n1 = edge.get<Node>()[1];

		std::vector<TCellID> n0_faces = n0.getIDs<Face>();
		std::vector<TCellID> n1_faces = n1.getIDs<Face>();

		std::set<TCellID> common_faces;

		for(auto f0 : n0_faces){
			for(auto f1 : n1_faces){
				if(f0 == f1){
					common_faces.insert(f0);
				}
			}
		}
		for(auto f : common_faces){
			//edge.add<Face>(f);
			Face face = blockMesh.get<Face>(f);
			face.add<Edge>(e);
		}
	}

	std::cout<<"nb edges mesh "<<blockMesh.getNbEdges()<<std::endl;
	std::cout<<"nb faces mesh "<<blockMesh.getNbFaces()<<std::endl;



	/*for(auto e : blockMesh.edges()){
	   Edge edge = blockMesh.get<Edge>(e);
	   Node n0 = edge.get<Node>()[0];
	   Node n1 = edge.get<Node>()[1];
	   n0.add<Edge>(e);
	   n1.add<Edge>(e);
	}*/

	cad::GeomMeshLinker linker;
	linker.setGeometry(&manager);
	linker.setMesh(&blockMesh);

	linker.linkNodeToPoint(2,6);
	linker.linkNodeToPoint(3,5);
	linker.linkNodeToPoint(5,3);
	linker.linkNodeToPoint(7,4);
	linker.linkNodeToPoint(11,1);
	linker.linkNodeToPoint(13,2);

	linker.linkNodeToCurve(9,1);
	linker.linkNodeToCurve(10,1);
	linker.linkNodeToCurve(12,1);

	linker.linkNodeToCurve(16,2);
	linker.linkNodeToCurve(17,2);

	linker.linkNodeToCurve(8,3);
	linker.linkNodeToCurve(14,3);
	linker.linkNodeToCurve(15,3);

	linker.linkNodeToCurve(0,5);
	linker.linkNodeToCurve(1,5);
	linker.linkNodeToCurve(4,5);
	linker.linkNodeToCurve(6,5);

	BlockMesher2D bm(&blockMesh,&linker,&manager);

	bm.updateClassification();

	bm.projectNodes();

	//std::cout<<"nb faces "<<n6.nbFaces()<<std::endl;
	//std::cout<<"nb edges "<<n6.nbEdges()<<std::endl;


	std::cout<<"nb points "<<manager.getNbPoints()<<std::endl;
	std::cout<<"nb curves "<<manager.getNbCurves()<<std::endl;
	//std::cout<<"nb surfaces "<<manager.getNbSurfaces()<<std::endl;

	for (int i = 0; i < 5; ++i) {
		for (auto n : blockMesh.nodes()) {/*
			 if (linker.getGeomDim<Node>(n) == 2) {
			    Node node = blockMesh.get<Node>(n);
			    std::vector<math::Point> voisins;
			    voisins.push_back(node.point());
			    std::cout<<"nb edges "<<node.get<Edge>().size()<<std::endl;
			    for (auto e_n : node.get<Edge>()) {
			       Node voisin = e_n.get<Node>()[0] != node ? e_n.get<Node>()[0] : e_n.get<Node>()[1];
			       if (linker.getGeomDim<Node>(voisin.id()) > 0) {
			          voisins.push_back(voisin.point());
			       }
			    }
			    std::cout<<"point avant "<<node.point()<<std::endl;
			    math::Point center(math::Point::massCenter(voisins));
			    cad::GeomCurve *c = manager.getCurve(linker.getGeomId<Node>(n));
			    std::cout<<"point avant proj "<<center<<std::endl;

			    c->project(center);
			    std::cout<<"point après proj "<<center<<std::endl;
			    node.setPoint(center);
			    std::cout<<"point du noeud apres "<<node.point()<<std::endl;

			 }*/
			/*else */if (linker.getGeomDim<Node>(n) == 0) {
				Node node = blockMesh.get<Node>(n);
				std::vector<math::Point> voisins;
				voisins.push_back(node.point());
				for (auto const &e_n : node.get<Edge>()) {
					Node voisin = e_n.get<Node>()[0] != node ? e_n.get<Node>()[0] : e_n.get<Node>()[1];
					voisins.push_back(voisin.point());
				}
				math::Point center(math::Point::massCenter(voisins));
				node.setPoint(center);
			}
		}
	}
	//smoother.smoothSurfaces(10);

	bm.executeMeshing();

	IGMeshIOService ios(&triMesh);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("/home/calderans/dev/test_mesh.vtk");

	IGMeshIOService iosb(&blockMesh);
	VTKWriter writerb(&iosb);
	writerb.setCellOptions(N|F);
	writerb.setDataOptions(N|F);
	writerb.write("/home/calderans/dev/test_blocks.vtk");

	Mesh* quad(bm.getMesh());

	for(auto n : quad->nodes()){
		//quad->newTriangle(n,n,n);
	}

	IGMeshIOService iosQ(quad);
	VTKWriter writerQ(&iosQ);
	writerQ.setCellOptions(N|F);
	writerQ.setDataOptions(N|F);
	writerQ.write("/home/calderans/dev/test_quad.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(BlockMesher2DTestSuite, test_blockCreation){
	Mesh triMesh(MeshModel(gmds::DIM3|R|gmds::F|gmds::E|gmds::N|gmds::N2E|gmds::E2N|gmds::E2F|gmds::N2F|gmds::F2N|F2E|N2R|R2N|E2R|R2E|R2F|F2R));

	std::string vtk_file_geom = "/home/calderans/dev/test_autoblocks.vtk";

	gmds::IGMeshIOService ioService(&triMesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(DIM3|gmds::N|gmds::E|gmds::F|R);
	vtkReader.read(vtk_file_geom);

	gmds::MeshDoctor doc(&triMesh);
	doc.buildEdgesAndX2E();
	doc.buildFacesAndR2F();
	doc.updateUpwardConnectivity();

	Mesh blockMesh(MeshModel(DIM3|R|F|E|N|N2E|E2N|E2F|N2F|F2N|F2E|N2R|R2N|E2R|R2E|F2R|R2F));

	cad::FACManager manager;
	cad::GeomMeshLinker linkerTri;
	manager.initAndLinkFrom3DMesh(&triMesh,&linkerTri);


	IGMeshIOService ios(&triMesh);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("/home/calderans/dev/test_mesh.vtk");

	IGMeshIOService iosQ(&manager.getMeshView());
	VTKWriter writerQ(&iosQ);
	writerQ.setCellOptions(N|F);
	writerQ.setDataOptions(N|F);
	writerQ.write("/home/calderans/dev/test_geom.vtk");

	cad::GeomMeshLinker linker;
	linker.setGeometry(&manager);
	linker.setMesh(&blockMesh);

	BlockClassificator classificator(&blockMesh,&linker,&manager);

	int matrix[5][5][5];
	for(int k = 0; k<5; k++) {
		for (int j = 0; j < 5; j++) {
			for (int i = 0; i < 5; i++) {
				matrix[i][j][k] = 0;
			}
		}
	}

	matrix[0][0][0] = 1;
	matrix[1][0][0] = 1;
	matrix[2][0][0] = 1;
	matrix[0][1][0] = 1;

	/*for(int i = 0; i<3; i++){
		for(int j = 0; j<3; j++){
			matrix[i][j][0] = 1;
		}
	}
	matrix[1][1][1] = 1;
	for(int i = 0; i<3; i++){
		for(int j = 0; j<3; j++){
			matrix[i][j][2] = 1;
		}
	}
	matrix[1][1][3] = 1;*/



	classificator.blockCreation(matrix);

	gmds::MeshDoctor docB(&blockMesh);
	docB.buildEdgesAndX2E();
	docB.buildFacesAndR2F();
	docB.updateUpwardConnectivity();



	std::cout<<"Nodes "<<blockMesh.getNbNodes()<<std::endl;
	std::cout<<"Edges "<<blockMesh.getNbEdges()<<std::endl;
	std::cout<<"Faces "<<blockMesh.getNbFaces()<<std::endl;
	std::cout<<"Regions "<<blockMesh.getNbRegions()<<std::endl;

	classificator.blockClassification();
	classificator.getBlocks();

	IGMeshIOService iosb(&blockMesh);
	VTKWriter writerb(&iosb);
	writerb.setCellOptions(N|F);
	writerb.setDataOptions(N|F);
	writerb.write("/home/calderans/dev/test_blocks_auto.vtk");
}