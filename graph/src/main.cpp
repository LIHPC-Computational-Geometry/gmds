//
// Created by calderans on 10/03/20.
//
#include <iostream>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/igalgo/GridBuilder.h>
#include "gmds/graph/TopologyGraph.h"

#include "gmds/graph/MinCut.h"

#include <gmds/smoothy/LaplacianSmoother.h>

using namespace gmds;
using namespace graph;

int main(int argc, char** argv)
{

    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));


    IGMeshIOService ioService(&m_vol);
    gmds::VTKReader vtkReader(&ioService);

    if (argc == 1 )
    {
        GridBuilder gb(&m_vol, 3);
        gb.execute(3,1.,3,1.,3,1.);
    }else if(argc == 2){
        std::cout<<argv[1]<<std::endl;
        std::string param_file(argv[1]);

        vtkReader.setCellOptions(gmds::N|gmds::R);
        vtkReader.setDataOptions(gmds::N|gmds::R);
        vtkReader.read(param_file);
    }else if(argc == 3){
        if(strcmp(argv[2],"s") == 0){
            std::string param_file(argv[1]);

            vtkReader.setCellOptions(gmds::N|gmds::R);
            //vtkReader.setDataOptions(gmds::N|gmds::R);
            vtkReader.read(param_file);

            MeshDoctor doc(&m_vol);
            doc.buildFacesAndR2F();
            doc.buildEdgesAndX2E();
            doc.updateUpwardConnectivity();

            cad::FACManager manager;
            cad::GeomMeshLinker linker;

            manager.initAndLinkFrom3DMesh(&m_vol,&linker);

            std::cout<<"Smoothing"<<std::endl;

            cad::GeomSmoother smoother(&linker);
            smoother.smoothVolumes(100);

            VTKWriter vtkWriter(&ioService);
            vtkWriter.setCellOptions(gmds::N|gmds::R);
            vtkWriter.setDataOptions(gmds::N|gmds::R);

            std::size_t pos = param_file.find(".vtk");
            std::string result_name(param_file.substr(0,pos)+"_smoothed.vtk");
            std::cout<<result_name<<std::endl;
            vtkWriter.write(result_name);

            return 0;
        }
    }


    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker,0.9);

    TopologyGraph graph(&manager, &linker, &m_vol);
    //graph.cutFace(419,445);
    //graph.colorSurfaces();
    //graph.writeGML();
    //graph.findIso();


    /*graph.pickSurf(1921);
    graph.pickSurf(1945);
    graph.pickSurf(3282);
    graph.pickSurf(3409);
    graph.pickSurf(2116);*/

    /*graph.pickSurf(110545);
    graph.pickSurf(48432);
    graph.pickSurf(135927);
    graph.pickSurf(137144);
    graph.pickSurf(95977);*/

    /*graph.pickSurf(718683);
    graph.pickSurf(574498);
    graph.pickSurf(671756);
    graph.pickSurf(906813);
    graph.pickSurf(719583);*/


    //raidisseur
    //source 12
    //graph.pickSurf(27667);
    //graph.pickSurf(24946);
    //graph.pickSurf(25809);

    //raidisseur2
    //source
    /*graph.pickSurf(8223);
    graph.pickSurf(11426);
    graph.pickSurf(15201);
    graph.pickSurf(21126);*/

    //raidisseur3
    //source
    //graph.pickSurf(28272);
    //graph.pickSurf(29264);
    //graph.pickSurf(13407);

    //graph.pickSurf(44291);
    //graph.pickSurf(39420);

    //graph.generateCut();

    /*DualGraph dual(&m_vol);
    dual.buildDualGraph();
*/

    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("/home/simon/Data/Results_debug/cut_face_res.vtk");

    //graph.cutFace(734,762);
    //vtkWriter.write("/home/simon/Data/Results_debug/cut_face_res2.vtk");


    /*std::string param_file(argv[1]);
    std::size_t pos = param_file.find(".vtk");
    std::string result_name(param_file.substr(0,pos)+"_cut.vtk");
    vtkWriter.write(result_name);*/

    return 0;
}