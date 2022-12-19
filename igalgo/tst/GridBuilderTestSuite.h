/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Blocking2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
TEST(GridBuildOpClass, test2D)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::F|gmds::F2N));

    gmds::GridBuilder gb(&m,2);

    ASSERT_TRUE(gb.isValid());

    gb.execute(3,1.0, 4, 1.0);

    ASSERT_EQ(m.getNbNodes(),12);
    ASSERT_EQ(m.getNbFaces(),6);

}
/*----------------------------------------------------------------------------*/
TEST(GridBuildOpClass, test3D)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

    gmds::GridBuilder gb(&m,3);

    ASSERT_TRUE(gb.isValid());

    gb.execute(3,1.0, 4, 1.0, 3, 2.0);

    ASSERT_EQ(m.getNbNodes(),36);
    ASSERT_EQ(m.getNbRegions(),12);

}
/*----------------------------------------------------------------------------*/
TEST(GridTestSuite, test_blocking2D_output)
{
    Blocking2D m;
    Node n1 = m.newBlockCorner(0,0);
    Node n2 = m.newBlockCorner(1,0);
    Node n3 = m.newBlockCorner(1,1);
    Node n4=  m.newBlockCorner(0,1);

    Blocking2D::Block b1 = m.newBlock(n1,n2,n3,n4);

    Node n5 = m.newBlockCorner(2,0,0);
    Node n6 = m.newBlockCorner(2,1.5,0);
    Blocking2D::Block b2 = m.newBlock(n2,n5,n6,n3);
    b1.setNbDiscretizationI(10);
    b1.setNbDiscretizationJ(10);
    b2.setNbDiscretizationI(10);
    b2.setNbDiscretizationJ(10);

    m.initializeGridPoints();

    IGMeshIOService ios(&m);
    VTKWriter writer(&ios);
    writer.setCellOptions(N|F);
    writer.setDataOptions(N|F);
    writer.write("blocking2D_sample.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(GridTestSuite, testWriterVTK_Val){
	std::cout << "==== Main Valentin Postat ====" << std::endl;
	std::cout << "limite min : " << std::numeric_limits<int64_t>::min() << std::endl;
	std::cout << "limite max : " << std::numeric_limits<int32_t>::max() << std::endl;

	// Generate mesh with connectivity
	gmds::Mesh m((gmds::MeshModel(DIM3|R|F|E|N|R2N|R2F|R2E|F2N|F2R|F2E|E2F|E2N|N2E)));
	// Call grid builder with mesh
	gmds::GridBuilder gb(&m,3);
	// Number of bloc for each dim
	TInt x_n=4,y_n=4,z_n=4;

	gb.execute(x_n,1.0, y_n, 1.0, z_n, 1.0);

	std::cout << "Nodes " << m.getNbNodes() << std::endl;
	std::cout << "Regions " << m.getNbRegions() << std::endl ;
	std::cout << "Edges " << m.getNbEdges() << std::endl;
	std::cout << "Faces " << m.getNbFaces() << std::endl;
	// Correct numberof bloc
	x_n=x_n-1,y_n=y_n-1,z_n=z_n-1;
	// Mesh doctor to generate connectivity
	gmds::MeshDoctor doc(&m);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	std::cout << "Nodes " << m.getNbNodes() << std::endl;
	std::cout << "Edges " << m.getNbEdges() << std::endl;
	std::cout << "Faces " << m.getNbFaces() << std::endl;
	std::cout << "Regions " << m.getNbRegions() << std::endl; // Well region looks to be Hexa

	// Assign boolean value to blocs (exists or not) 0 no 1 yes
	Variable<int>* bloc_exist = m.newVariable<int, gmds::GMDS_REGION>("bloc_exist");
	// Assign boolean value to face (Constrained Strong or not) 0 no 1 yes
	Variable<int>* face_constraint = m.newVariable<int, gmds::GMDS_FACE>("face_constraint");
	// Assign boolean value to face (In solution or not) 0 no 1 yes
	Variable<int>* face_in_solution = m.newVariable<int, gmds::GMDS_FACE>("face_in_solution");

	// Writing
	std::string fname("test0.vtk");
	IGMeshIOService ios(&m);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F|R);
	writer.setDataOptions(N|F|R);
	writer.write(fname); // paraview ctrl + space (Shrink) pour avoir uniquement les hex (enlever les F dans l'écriture)

	// Generate our problem boolean value
	int exist_tens[x_n][y_n][z_n];
	int constraints_tens[x_n + 1][y_n + 1][z_n + 1];

	// Given a tensor (n_xn_yn_z) labelise blocs that exists
	for (auto i = 0; i < x_n; i++) {
		for (auto j = 0; j < y_n; j++) {
			for (auto k = 0; k < z_n; k++) {
				exist_tens[i][j][k] = 1;
			}
		}
	}

	// Init variables
	face_constraint->setValuesTo(0);
	bloc_exist->setValuesTo(0);
	face_in_solution->setValuesTo(0);


	std::cout << "Faces" << std::endl;
	for(auto i:m.faces()){
		std::cout << i << std::endl;
		face_constraint->set(i,i);
	}
	std::cout << "Faces" << std::endl;
	for(auto i:m.faces()){
		std::cout << i << std::endl;
		face_in_solution->set(i,i);
	}
	std::cout << "Regions" << std::endl;
	for(auto i:m.regions()){
		std::cout << i << std::endl;
		bloc_exist->set(i,i);
	}
	face_constraint->compact();
	face_in_solution->compact();
	bloc_exist->compact();

	std::cout << face_constraint->getNbValues() << " " << face_in_solution->getNbValues() << " " << bloc_exist->getNbValues() << std::endl;

	// En x de 1 et y de 1
	// Given a tensor (n_x+1n_y+1*n_z+1) labelize faces that are constraints SC
}