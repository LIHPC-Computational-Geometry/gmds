/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 7/27/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/frame/CrossFieldGeneration2D.h>

/*----------------------------------------------------------------------------*/
#include <iostream>
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "==== 2D Cross Field Generator====" << std::endl;
    Mesh m(MeshModel(DIM3|F|N|E|N2E|N2F|F2N|E2N|F2E|E2F));


    //==================================================================
    // MESH READING
    //==================================================================
    std::cout << "Reading " << std::endl;
    std::string fIn, fOut;
    fOut = fIn;

    if (argc < 3)
        throw gmds::GMDSException("Wrong parameters: need input file + output file (.vtk)");

    if (argc ==3) {
        fIn = std::string(argv[1]);
        fOut = std::string(argv[2]);
    }


    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(N|F);
    vtkReader.read(fIn);


    gmds::MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();
    doc.orient2DFaces();

    CrossFieldGeneration2D field_generator(&m);
    //directory where the debug files are going to be written
    field_generator.setDebugPrefix("cross2D-debug-");
    field_generator.execute(CrossFieldGeneration2D::laplace_solve);


    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write(fOut);

}