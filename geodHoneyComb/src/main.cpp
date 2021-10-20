/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/geodHoneyComb/GeodHexMesher.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/ig/MeshDoctor.h>

using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{

    std::cout << "==== 3D TYPE DISPLAYER====" << std::endl;
    Mesh m(MeshModel(DIM3 | R | N | R2N));

    //==================================================================
    // MESH READING
    //==================================================================
    std::cout << "Reading " << std::endl;
    std::string fIn, fOut;

    fIn = std::string(argv[1]);
    fOut = std::string(argv[2]);

    std::cout << "INPUT FILE: " << fIn << std::endl;
    std::cout << "OUTPUT FILE: " << fOut << std::endl;


    IGMeshIOService ioService(&m);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(fIn);



    Variable<int>* var_type =   m.getOrCreateVariable<int,GMDS_REGION>("GMDS_TYPE");

    for(auto r_id:m.regions()){
        Region r = m.get<Region>(r_id);
        var_type->value(r_id)=r.type();
    }
    std::cout<<"\t nb hexes= "<<m.getNbHexahedra()<<std::endl;
    std::cout<<"\t nb tets = "<<m.getNbTetrahedra()<<std::endl;
    VTKWriter writer2(&ioService);
    writer2.setCellOptions(gmds::N|gmds::R);
    writer2.setDataOptions(gmds::N|gmds::R);
    writer2.write("GeneratedPoints.vtk");


}
