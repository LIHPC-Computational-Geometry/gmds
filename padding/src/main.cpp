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
#include <gmds/padding/SelectivePadding.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "=== SMART PADDING ====" << std::endl;
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | F2N | E2N | R2F | F2R |
                     F2E | E2F | R2E | E2R| N2R | N2F | N2E));

    //==================================================================
    // MESH READING
    //==================================================================
    std::cout << "Reading " << std::endl;
    std::string fIn, fOut;

    if (argc !=3)
        throw gmds::GMDSException("[Wrong parameters] usage should be PaddingExe in.vtk out.vtk");

    fIn = std::string(argv[1]);
    fOut = std::string(argv[2]);


    IGMeshIOService ioService(&m);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(fIn);



    //==================================================================
    // MESH PREPARATION
    //==================================================================
    MeshDoctor doctor(&m);
    doctor.buildFacesAndR2F();
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();

    //====================================================================
    // Init hard constraints
    double min_qual = 0.3;
    Variable<double>* quality_var = m.newVariable<double, GMDS_REGION>("GMDS_Quality");
    Variable<int>* hard = m.newVariable<int, GMDS_FACE>("GMDS_HARD");

    for(auto cell_id:m.regions()) {
        Region r = m.get<Region>(cell_id);
        double qi = r.computeScaledJacobian();
        (*quality_var)[cell_id]= qi;
        if(qi<min_qual) {
            //We get boundary faces of r and constrained them
            std::vector<Face> faces_of_r = r.get<Face>();
            for(auto f:faces_of_r){
                if(f.getIDs<Region>().size()==1){
                    hard->set(f.id(), 1);
                }
            }
        }
    }
    //==================================================================
    // SELECTIVE PADDING ALGORITHM
    //==================================================================
    SelectivePadding algo(&m);

    Variable<int>* xfi = m.newVariable<int, GMDS_FACE>("GMDS_XFI");

    algo.setHardFaces(hard);
    algo.setPaddingFaces(xfi);
    algo.execute();
    //==================================================================
    // MESH WRITING
    //==================================================================
    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N | gmds::F);
    vtkWriter.setDataOptions(gmds::N | gmds::F);
    vtkWriter.write(fOut);

}
