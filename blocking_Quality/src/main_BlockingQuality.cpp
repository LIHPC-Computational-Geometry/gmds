//
// Created by bourmaudp on 02/12/22.
//
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <string>


#include <gmds/ig/Mesh.h>
#include <gmds/ig/Node.h>
#include <gmds/quality/QuadQuality.h>
#include <gmds/quality/HexQuality.h>


#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/ig/MeshDoctor.h>

//================================================================================
using namespace gmds;
int main(){
	std::cout<<"TEST"<<std::endl;

	/*
	Mesh mImprint (	MeshModel(DIM3|F|E|N|F2N|N2F|E2N));
	IGMeshIOService ioServiceRead(&mImprint);
	VTKReader vtkReader(&ioServiceRead);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	//vtkReader.read("/home/bourmaudp/Documents/GMDS/gmds/test_samples/HolesInSquare0.vtk");
	vtkReader.read("/home/bourmaudp/Documents/GMDS/gmds/test_samples/Cube.vtk");
	 */
	Mesh simpleCube(MeshModel(DIM3|F|E|N|F2N|N2F));
	Node n0 = simpleCube.newNode(math::Point(0,0,0));
	Node n1 = simpleCube.newNode(math::Point(1,0,0));
	Node n2 = simpleCube.newNode(math::Point(1,0,-1));
	Node n3 = simpleCube.newNode(math::Point(0,0,-1));
	Node n4 = simpleCube.newNode(math::Point(0,1,0));
	Node n5 = simpleCube.newNode(math::Point(1,1,0));
	Node n6 = simpleCube.newNode(math::Point(1,1,-1));
	Node n7 = simpleCube.newNode(math::Point(0,1,-1));

	simpleCube.newQuad(n0,n1,n2,n3);
	simpleCube.newQuad(n4,n5,n6,n7);
	simpleCube.newQuad(n1,n2,n6,n5);
	simpleCube.newQuad(n0,n3,n7,n4);
	simpleCube.newQuad(n0,n1,n5,n4);
	simpleCube.newQuad(n3,n2,n6,n7);

	Mesh simpleRectangle(MeshModel(DIM3|F|E|N|F2N|N2F));
	Node nr0 = simpleRectangle.newNode(math::Point(0,0,0));
	Node nr1 = simpleRectangle.newNode(math::Point(3,0,0));
	Node nr2 = simpleRectangle.newNode(math::Point(3,0,-1));
	Node nr3 = simpleRectangle.newNode(math::Point(0,0,-1));
	Node nr4 = simpleRectangle.newNode(math::Point(0,1,0));
	Node nr5 = simpleRectangle.newNode(math::Point(3,1,0));
	Node nr6 = simpleRectangle.newNode(math::Point(3,1,-1));
	Node nr7 = simpleRectangle.newNode(math::Point(0,1,-1));

	simpleRectangle.newQuad(nr0,nr1,nr2,nr3);
	simpleRectangle.newQuad(nr4,nr5,nr6,nr7);
	simpleRectangle.newQuad(nr1,nr2,nr6,nr5);
	simpleRectangle.newQuad(nr0,nr3,nr7,nr4);
	simpleRectangle.newQuad(nr0,nr1,nr5,nr4);
	simpleRectangle.newQuad(nr3,nr2,nr6,nr7);

	Mesh doubleRectangle(MeshModel(DIM3|F|E|N|F2N|N2F));
	Node dnr0 = doubleRectangle.newNode(math::Point(0,0,0));
	Node dnr1 = doubleRectangle.newNode(math::Point(3,0,0));
	Node dnr2 = doubleRectangle.newNode(math::Point(3,0,-1));
	Node dnr3 = doubleRectangle.newNode(math::Point(0,0,-1));
	Node dnr4 = doubleRectangle.newNode(math::Point(0,1,0));
	Node dnr5 = doubleRectangle.newNode(math::Point(3,1,0));
	Node dnr6 = doubleRectangle.newNode(math::Point(3,1,-1));
	Node dnr7 = doubleRectangle.newNode(math::Point(0,1,-1));
	Node dnr8 = doubleRectangle.newNode(math::Point(6,0,0));
	Node dnr9 = doubleRectangle.newNode(math::Point(6,0,-1));
	Node dnr10 = doubleRectangle.newNode(math::Point(6,1,0));
	Node dnr11 = doubleRectangle.newNode(math::Point(6,1,-1));

	doubleRectangle.newQuad(dnr0,dnr1,dnr2,dnr3);
	doubleRectangle.newQuad(dnr4,dnr5,dnr6,dnr7);
	doubleRectangle.newQuad(dnr1,dnr2,dnr6,dnr5);
	doubleRectangle.newQuad(dnr0,dnr3,dnr7,dnr4);
	doubleRectangle.newQuad(dnr0,dnr1,dnr5,dnr4);
	doubleRectangle.newQuad(dnr3,dnr2,dnr6,dnr7);

	doubleRectangle.newQuad(dnr1,dnr8,dnr9,dnr3);
	doubleRectangle.newQuad(dnr1,dnr8,dnr10,dnr5);
	doubleRectangle.newQuad(dnr5,dnr10,dnr11,dnr6);
	doubleRectangle.newQuad(dnr3,dnr9,dnr11,dnr6);
	doubleRectangle.newQuad(dnr8,dnr9,dnr11,dnr10);




	double qualityValue = 0;
	/*
	std::cout<<"List faces : "<<std::endl;
	for (auto f : simpleCube.faces()){
		std::cout<<"Face : "<<f<<std::endl;
		std::vector<Node> listNodes = simpleCube.get<Face>(f).get<Node>();
		std::cout<<"List nodes : "<<std::endl;
		for (auto n : listNodes){
			std::cout<<"Node : "<<std::endl;

		}
	}*/

	//HexQuality simpleCube
	//quality::HexQuality he = quality::HexQuality::build(n0.point(),n1.point(),n2.point(),n3.point(),n4.point(),n5.point(),n6.point(),n7.point());

	//HexQuality simpleRectangle
	//quality::HexQuality he = quality::HexQuality::build(nr0.point(),nr1.point(),nr2.point(),nr3.point(),nr4.point(),nr5.point(),nr6.point(),nr7.point());

	//HexQuality bigRectangle
	quality::HexQuality he = quality::HexQuality::build(dnr0.point(),dnr1.point(),dnr2.point(),dnr3.point(),dnr4.point(),dnr5.point(),dnr6.point(),dnr7.point());

	//HexQuality 2 rectangles -> bigRectangle
	//quality::HexQuality he = quality::HexQuality::build(dnr0.point(),dnr1.point(),dnr2.point(),dnr3.point(),dnr4.point(),dnr5.point(),dnr6.point(),dnr7.point());
	//quality::HexQuality he2 = quality::HexQuality::build(dnr1.point(),dnr8.point(),dnr9.point(),dnr3.point(),dnr5.point(),dnr10.point(),dnr11.point(),dnr6.point());


	/*
	std::cout<<"Quality bloc 1"<<std::endl;
	std::cout<<"Le volume "<<he.volume()<<std::endl;
	std::cout<<"La diagonal "<<he.diagonal()<<std::endl;
	std::cout<<"Le edgeRatio "<<he.edgeRatio()<<std::endl;
	std::cout<<"Le maxi edgeRatio "<<he.maximumEdgeRatio()<<std::endl;
	std::cout<<"Le Jacobien "<<he.jacobian()<<std::endl;
	std::cout<<"Le Scaled Jacobien "<<he.scaledJacobian()<<std::endl;
	std::cout<<"Le shape "<<he.shape()<<std::endl;
	std::cout<<"Le maxiAF "<<he.maximumAspectFrobenius()<<std::endl;
	std::cout<<"Le meanAF "<<he.meanAspectFrobenius()<<std::endl;
	std::cout<<"Le oddy "<<he.oddy()<<std::endl;
	std::cout<<"Le shear "<<he.shear()<<std::endl;
	std::cout<<"Le skew "<<he.skew()<<std::endl;
	std::cout<<"Le stretch "<<he.stretch()<<std::endl;
	std::cout<<"Le taper "<<he.taper()<<std::endl;
	 */
	std::cout<<"TEST EDGE RATIO MODIFIED "<<he.edgeRatioModified()<<std::endl;

	/*
	std::cout<<"Quality bloc 2"<<std::endl;
	std::cout<<"Le volume "<<he2.volume()<<std::endl;
	std::cout<<"La diagonal "<<he2.diagonal()<<std::endl;
	std::cout<<"Le edgeRatio "<<he2.edgeRatio()<<std::endl;
	std::cout<<"Le maxi edgeRatio "<<he2.maximumEdgeRatio()<<std::endl;
	std::cout<<"Le Jacobien "<<he2.jacobian()<<std::endl;
	std::cout<<"Le Scaled Jacobien "<<he2.scaledJacobian()<<std::endl;
	std::cout<<"Le shape "<<he2.shape()<<std::endl;
	std::cout<<"Le maxiAF "<<he2.maximumAspectFrobenius()<<std::endl;
	std::cout<<"Le meanAF "<<he2.meanAspectFrobenius()<<std::endl;
	std::cout<<"Le oddy "<<he2.oddy()<<std::endl;
	std::cout<<"Le shear "<<he2.shear()<<std::endl;
	std::cout<<"Le skew "<<he2.skew()<<std::endl;
	std::cout<<"Le stretch "<<he2.stretch()<<std::endl;
	std::cout<<"Le taper "<<he2.taper()<<std::endl;
	 */
	/*
	gmds::MeshDoctor doc_imprint(&simpleCube);
	doc_imprint.updateUpwardConnectivity();
	doc_imprint.buildBoundaryCells();

	std::cout<<"Le volume apres doctor "<<he.volume()<<std::endl;
	*/

	/*
	gmds::IGMeshIOService ioService_write(&simpleCube);
	gmds::VTKWriter vtkWriterGba(&ioService_write);
	vtkWriterGba.setCellOptions(gmds::N|gmds::F);
	vtkWriterGba.setDataOptions(gmds::N|gmds::F);
	vtkWriterGba.write("cubeSimple.vtk");
	*/
}