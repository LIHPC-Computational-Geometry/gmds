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

    if (argc != 2 )
    {
        GridBuilder gb(&m_vol, 3);
        gb.execute(3,1.,3,1.,3,1.);
    }else{
        std::string param_file(argv[1]);

        vtkReader.setCellOptions(gmds::N|gmds::R);
        vtkReader.setDataOptions(gmds::N|gmds::R);
        vtkReader.read(param_file);
    }


    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

    TopologyGraph graph(&manager,&m_vol);
    graph.buildGraph();
    graph.printGraph();
    graph.colorSurfaces();
    graph.writeGML();
    graph.findIso();

    VTKWriter vtkWriter_i(&ioService);
    vtkWriter_i.setCellOptions(gmds::N|gmds::F);
    vtkWriter_i.setDataOptions(gmds::N|gmds::F);
    vtkWriter_i.write("/home/simon/Dev/Results/graph_test.vtk");

    return 0;
}