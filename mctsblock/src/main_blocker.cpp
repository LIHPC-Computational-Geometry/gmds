/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/mctsblock/Blocking.h>

using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "============== MCTS locker ================" << std::endl;

//==================================================================
// PARAMETERS' PARSING
//==================================================================
    std::string file_geom, file_block_out;
    if (argc != 3) {
        std::cout << "Require 2 paramaters : \n";
        std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
        std::cout << "  - [OUT] the final blocks (.vtk),\n"<< std::endl;
        throw gmds::GMDSException("Wrong number of parameters");
    }

    file_geom = std::string(argv[1]);
    file_block_out = std::string(argv[7]);
    std::cout << "Parameters " << std::endl;
    std::cout << "  - Geometry file: " << file_geom << std::endl;
    std::cout << "  - Output block file  : " << file_block_out << std::endl;
    std::cout << "=======================================" << std::endl;

//==================================================================
// GEOMETRY READING
//==================================================================
    std::cout<<"> Start geometry reading"<<std::endl;
    Mesh geometry(MeshModel(DIM3 | R | F | E | N |
                            R2N | R2F | R2E |
                            F2N | F2R | F2E |
                            E2F | E2N | N2E | N2R | N2F));

    IGMeshIOService ioService(&geometry);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(N| R);
    vtkReader.read(file_geom);
    MeshDoctor doc(&geometry);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager geom_manager;
    geom_manager.initFrom3DMesh(&geometry);

	 mctsblock::Blocking blocking(&geom_manager,true);

    std::cout << "======== Task done by blocker =========" << std::endl;
}