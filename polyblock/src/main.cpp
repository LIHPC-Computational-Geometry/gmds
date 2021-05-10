/*----------------------------------------------------------------------------*/
#include <iostream>
#include <string>
/*----------------------------------------------------------------------------*/
#include <Eigen/Eigen>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds//ig/MeshDoctor.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/igalgo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
#include <gmds/polyblock/Gregson2011.h>
#include <gmds/polyblock/PolycubeToolbox.h>
/*----------------------------------------------------------------------------*/
// make -j4 && echo && ./Polycube/Polycube ../MeshData/FilConducteur/cylindre.mesh
/*----------------------------------------------------------------------------*/
using namespace std;
using namespace Eigen;
using namespace gmds;
typedef MatrixXd matrice;
typedef VectorXd vecteur;
/*----------------------------------------------------------------------------*/
void ex(Mesh &mesh){
	for(auto e_id : mesh.edges()){
		Edge e = mesh.get<Edge>(e_id);
		std::cout << e.id() << " : " << e.nbFaces() << '\n';
	}
	std::cout << mesh.getNbEdges() << '\n';
}
/*----------------------------------------------------------------------------*/
Mesh meshloader(string AFileName){
    Mesh mesh(MeshModel(DIM3|R|F|E|N|
                                R2N|R2F|F2R|R2E|F2N|E2N|E2F|F2E|E2R|N2F|N2E|N2R));
	std::cout << "read mesh file: " << AFileName << std::endl;
  	IGMeshIOService ioService(&mesh);
    MeditReader reader(&ioService);
    reader.setCellOptions(N|R);
    reader.setDataOptions(N|R);
  	reader.read(AFileName);

  	std::cout<<"Updating and completing mesh info ...";
    std::cout.flush();

  	MeshDoctor doc(&mesh);
  	doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
  	doc.updateUpwardConnectivity();

  	std::cout<<" Loading done."<<std::endl;
  	return mesh;
}
/*----------------------------------------------------------------------------*/
Mesh meshloader2D(string AFileName){
  	Mesh mesh(MeshModel(DIM3|F|E|N|
  						  F2N|E2N|E2F|F2E|N2F|N2E));
	std::cout << "read mesh file: " << AFileName << std::endl;
  	IGMeshIOService ioService(&mesh);
  	MeditReader reader(&ioService);
    reader.setCellOptions(N|F);
    reader.setDataOptions(N|F);
    reader.read(AFileName);
  	std::cout<<"Updating and completing mesh info ...";
    std::cout.flush();

  	MeshDoctor doc(&mesh);
    doc.buildEdgesAndX2E();
  	doc.updateUpwardConnectivity();
  	std::cout<<" Loading done."<<std::endl;
  	return mesh;
}
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[]){
    cout<<"===== Polycube executable ====="<<endl;

    std::string file_name(argv[1]);

    if (argc > 2 ){
        std::string dimension(argv[2]);
        if (dimension == "2D" || dimension == "2"){
            Mesh mesh2D = meshloader2D(file_name);
            IGMeshIOService ioService(&mesh2D);
            VTKWriter w_2D(&ioService);
            w_2D.setCellOptions(N|E|F);
            w_2D.setDataOptions(N|E|F);
            w_2D.write("result0_2D.vtk");

            Gregson2011_2D *algo2D = new Gregson2011_2D(&mesh2D);
            algo2D->transformMeshToPolycube();

            w_2D.write("result2D_final");
        }
    }
    else {
        Mesh mesh = meshloader(file_name);
        IGMeshIOService ioService(&mesh);
        VTKWriter w(&ioService);
        w.setCellOptions(N|F);
        w.setDataOptions(N|F);
        w.write("result0_3D.vtk");

        std::cout<<"Mesh info :"<<std::endl;
        std::cout<<"-Nb regions= "<<mesh.getNbRegions()<<std::endl;
        std::cout<<"-Nb faces  = "<<mesh.getNbFaces()<<std::endl;
        std::cout<<"-Nb edges  = "<<mesh.getNbEdges()<<std::endl;
        std::cout<<"-Nb nodes  = "<<mesh.getNbNodes()<<std::endl;


        PolycubeToolbox *algo = new PolycubeToolbox(&mesh);

        //algo->graphCutRun();
         algo->gregson2011Run();

        std::cout << "Writing results ..." << '\n';
        w.write("resultat_final.vtk");

        /* HOW TO EXTRACT CORNERS */
        Mesh mesh_test = meshloader(file_name);
        std::vector<std::vector<TCellID>> corners = algo->get_faces_of_block();

        mesh_test.deleteVariable(GMDS_NODE,"corner");
        gmds::Variable<int>* 	corner = mesh_test.newVariable<int,GMDS_NODE>("corner");
        std::set<TCellID> corner_ids;
        for (auto c : corners /* this give use a vector of corner for each chart*/){
            for (auto p : c /* this give use the TCellID of the point */){
                (*corner)[p] = 1;
                corner_ids.insert(p);
            }
        }
        std::vector<TCellID> all_point_corners;
        std::map<TCellID,int> vec_position;
        for (auto c : corner_ids /* this give use a vector of corner for each chart*/) {
            all_point_corners.push_back(c);
            vec_position[c]=all_point_corners.size();
        }
        //Nodes are now written as they are stored in the vec
        for(auto c:all_point_corners) {
            math::Point pc = mesh_test.get<Node>(c).getPoint();
            std::cout << "create node "<<vec_position[c]<< " at location " << pc.X()<<" "<<pc.Y()<<" "<<pc.Z() << std::endl;
        }
        //and now each face
        for (auto c : corners /* this give use a vector of corner for each chart*/){
            std::cout << "create face node ";
            std::set<TCellID> only1;
            for (auto p : c /* this give use the TCellID of the point */){
                only1.insert(p);
            }
            for (auto p : only1 /* this give use the TCellID of the point */){
                std::cout<<vec_position[p]<<" ";
            }
            std::cout<<std::endl;
        }
        IGMeshIOService ioService2(&mesh_test);
        VTKWriter w2(&ioService2);
        w2.setCellOptions(N|F);
        w2.setDataOptions(N|F);
        w2.write("DISPLAY_OF_CORNERS");


        Mesh bnd_mesh(MeshModel(DIM3|F|N|F2N));
        for(auto f_id: mesh.faces()){
            Face f = mesh.get<Face>(f_id);
            if(f.get<Region>().size()==1){
                std::vector<Node> nf = f.get<Node>();
                Node n0 = bnd_mesh.newNode(nf[0].getPoint());
                Node n1 = bnd_mesh.newNode(nf[1].getPoint());
                Node n2 = bnd_mesh.newNode(nf[2].getPoint());
                bnd_mesh.newTriangle(n0,n1,n2);
            }
        }

        IGMeshIOService ioServiceBnd(&bnd_mesh);
        VTKWriter bw(&ioServiceBnd);
        bw.setCellOptions(N|F);
        bw.setDataOptions(N|F);
        bw.write("bnd_final");
    }

    return 0;
}

