/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <fstream>
#include <bitset>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/MeshTransformation.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_mesh_transformation)
{
  TInt border = std::numeric_limits<TInt>::min();
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/simpleCubeV2.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();

  //adding metric to the mesh
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);
  //
  const gmds::BitVector & nodesIds = simplexMesh.getBitVectorNodes();
  double epsilon = 1E-1;
  Variable<Eigen::Matrix3d>* metricNode = simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  for(unsigned int nodeId = 0 ; nodeId < nodesIds.capacity() ; nodeId++)
  {
    if(nodesIds[nodeId] == 1)
    {
      double x = epsilon + SimplicesNode(&simplexMesh, nodeId).getCoords()[0];
      (*metricNode)[nodeId](0,0) = x;
      (*metricNode)[nodeId](1,1) = x;
      (*metricNode)[nodeId](2,2) = x;
    }
  }
  MeshTransformation mT(&simplexMesh);
  mT.transformation();
}
