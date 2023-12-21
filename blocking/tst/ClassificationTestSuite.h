/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
/**@brief setup function that initialize a geometric model using the faceted
 * representation and an input vtk file name. The vtk file must contain a
 * tetrahedral mesh
 *
 * @param AGeomModel geometric model we initialize
 * @param AFileName vtk filename
 */
void set_up_Geom(gmds::cad::FACManager* AGeomModel, const std::string AFileName)
{
    gmds::Mesh vol_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
                                        | gmds::E2F | gmds::E2N | gmds::N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir +"/"+ AFileName;
    gmds::IGMeshIOService ioService(&vol_mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(vtk_file);
    gmds::MeshDoctor doc(&vol_mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();
    AGeomModel->initFrom3DMesh(&vol_mesh);

}
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingClassifierTestSuite, testClass) {

    gmds::cad::FACManager geom_model;
    set_up_Geom(&geom_model,"Cube.vtk");
    gmds::blocking::CurvedBlocking bl(&geom_model,true);
    gmds::blocking::CurvedBlockingClassifier classifier(&bl);

}