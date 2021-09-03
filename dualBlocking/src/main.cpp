//
// Created by calderans on 25/03/19.
//


/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include "gmds/dualBlocking/DualBlockingSession.h"

using namespace gmds;
using namespace db;
/*----------------------------------------------------------------------------*/
int main(int argc, char** argv)
{
    std::cout << "Dual Blocking" << std::endl;
    if (argc != 2 )
    {
        std::cout << "Usage: toto *.vtk" << std::endl;
        exit(0);
    }
    std::string param_file(argv[1]);

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();

    Mesh imesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doci(&imesh);
    doci.buildFacesAndR2F();
    doci.buildEdgesAndX2E();
    doci.updateUpwardConnectivity();

    std::cout<<"Nb cells: "<<mesh.getNbRegions()<<std::endl;
    std::cout<<"Nb nodes: "<<mesh.getNbNodes()<<std::endl;

    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    for(auto n_id: mesh.nodes()){
        Node n = mesh.get<Node>(n_id);
        if(n.getIDs<Edge>().size()==0){
            std::cout<<"HERE MAIN 2) One node "<<n_id<<" has no incident edges!!"<<std::endl;
            std::cout<<n.getPoint()<<std::endl;
            //exit(0);
        }
    }


    DualBlockingSession session(&mesh,&hmesh);
    //session.init();

    //auto_intersect
    //session.createSurface(79574,x,1);

    //session.createBoundary(26321,2);
    //session.createBoundary(46091,3);

    //S3
    //session.createBoundary(289320,17);
    //session.createBoundary(716969,18);

    //S7
    //session.createBoundary(130115,8);
    //session.createBoundary(130918,9);

    //B8
    //session.createBoundary(69930,7);

    //session.createSurface(22080,x,1);
    //session.createSurface(22080,y,2);
    //session.createSurface(4559,y,3);
    //session.createSurface(21122,x,4);
    //session.createSurface(57194,z,5);

    //B44
    //session.createSurface(113345,x,1);
    //session.createSurface(113345,z,2);
    //session.classifyZone(0);
    //session.createSurface(119769,x,3);
    //session.createSurface(119769,z,4);
    //session.createSurface(118924,y,5);
    //session.createSurface(94003,y,6);

    //S34
    //session.createBoundary(104068,10);

    //S37
    //session.createBoundary(583212,7);

    //S40
    //session.createBoundary(376044,21);
    //session.createBoundary(378419,22);
    //session.createBoundary(384685,23);

    //S41
    //session.createBoundary(353647,27);
    //session.createBoundary(311745,28);
    //session.createBoundary(345176,29);
    //session.createBoundary(367717,30);

    //B28
    //session.createBoundary(300157,21);
    //session.createBoundary(76463,22);
    //session.createBoundary(318754,23);

    //B10
    //session.createBoundary(229690,10);

    //B31_2
    //session.createBoundary(329808,11);
    //session.createBoundary(364274,12);
    //session.createBoundary(282801,13);

    //EDF
    /*session.createSurface(33419,x,1);

    session.createSurface(33419,z,2);
    session.createSurface(5835,z,3);
    session.createSurface(24236,z,4);
    session.createSurface(30304,x,5);
    session.createSurface(35016,x,6);
    session.createBoundary(15091,7);
    session.createBoundary(29100,8);*/

    //B41
    //session.createSurface(18011,x,1);
    //session.createSurface(18011,y,2);
    //session.createSurface(33511,z,3);

    //session.createSurface(12769,x,1);
    //session.createSurface(7353,y,2);
    //session.createSurface(7353,z,3);
    //session.createSurface(18526,z,4);

    //B19
    //session.createSurface(2313,z,1);

    //S7
    //session.createSurface(535729,x,1);
    //session.createSurface(77888,x,2);
    //session.createSurface(77888,z,3);
    //session.createSurface(277616,y,4);
    /*session.createSurface(493887,z,5);
    session.createSurface(295049,z,6);
    session.createSurface(256395,x,7);
    session.createSurface(292856,x,8);
    session.createSurface(332623,z,9);
    session.createSurface(19490,z,10);
    session.createSurface(323089,z,11);
    session.createSurface(162930,y,12);
    session.createBoundary(369673,13);
    session.createBoundary(126610,14);*/

    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::R);
    //vtkWriter.write("/ccc/temp/cont001/ocre/calderans/Results_debug/block.vtk");

    //session.classifyZone(0);

    //session.createSurface(12769,x,1);
    //session.createSurface(7353,y,2);
    //session.createSurface(7353,z,3);
    //session.createSurface(18526,z,4);


    //B31_1
    //session.createBoundary(284754,13);
    //session.createBoundary(28928,14);

    //B31_2
    //session.createBoundary(364274,12);
    //session.createBoundary(282801,13);

    //S38
    //session.createBoundary(171924,17);
    //session.createBoundary(165129,18);
    //session.createBoundary(197073,19);


    //session.createBoundary(45591,1);
    //session.createBoundary(20372,1);

    //session.refineSurfaceSheet();

    //doc.buildFacesAndR2F();
    //doc.buildEdgesAndX2E();
    //doc.updateUpwardConnectivity();

    //session.getSheetVisualisation(21);
    //session.getSheetVisualisation(22);
    //session.getSheetVisualisation(23);



    /*if(!session.colorDual()){
        vtkWriter.write("/ccc/temp/cont001/ocre/calderans/Results_debug/block.vtk");
        std::cout<<"Color failed"<<std::endl;
	}*/
    //session.createBlock();



    std::cout<<"End of program"<<std::endl;

}
