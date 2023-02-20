/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/smoothy/LaplacianSmoother.h>
#include <gmds/igalgo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
#include <iostream>

/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "============== Morph 2D ================" << std::endl;

    //==================================================================
    // PARAMETERS' PARSING
    //==================================================================
    std::string file_input;
    if (argc != 2) {
        std::cout << "Require two parameters : \n";
        std::cout << "  - [IN ] input mesh (.vtk) that we must deform, \n";
        exit(0);
    }

    file_input = std::string(argv[1]);
    std::cout << "Parameters " << std::endl;
    std::cout << "  - Input file : " << file_input << std::endl;
    std::cout << "=======================================" << std::endl;

    //==================================================================
    // GEOMETRY READING
    //==================================================================
    std::cout<<"> Start geometry reading"<<std::endl;
    Mesh mesh_in(MeshModel(DIM3 | F | E | N |
                            F2N | F2E | E2N | N2E | N2R | N2F));

    IGMeshIOService ioService(&mesh_in);
    VTKReader vtkReader(&ioService);
	 vtkReader.setCellOptions(N| F | E);
	// vtkReader.setDataOptions(N| F | E);
    vtkReader.read(file_input);
	 std::cout<<"Stats: ("<<mesh_in.getNbFaces()<<", "<<mesh_in.getNbEdges()<<")"<<std::endl;
    MeshDoctor doc(&mesh_in);
//    doc.buildEdgesAndX2E();
 //   doc.updateUpwardConnectivity();


	std::string file_output ="morph2d_result.vtk";
    std::cout<<"> Write output mesh in: "<<file_output<<std::endl;
    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(N|F|E);
    vtkWriter.write(file_output);
    std::cout << "======== Task done =========" << std::endl;
}