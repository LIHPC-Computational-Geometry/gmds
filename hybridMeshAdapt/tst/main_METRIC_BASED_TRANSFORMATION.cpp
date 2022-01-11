/*----------------------------------------------------------------------------*/
#include <fstream>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/MeshTransformation.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {

  //==================================================================
  // METRIC BASED TRANSFORMATION MESH READING
  //==================================================================
  std::cout << "==== Metric Based Transformation Running ====" << std::endl;
  if(argc != 2)
  {
      throw gmds::GMDSException("NO INPUT FILE");
  }
  std::string fOut;
  const std::string fIn = std::string(argv[1]);
  const std::string extansion = ".vtk";
  if (fIn.find(extansion) == std::string::npos) {
    throw gmds::GMDSException("NOT A .vtk FILE");
  }
  std::cout << "INPUT FILE: " << fIn << std::endl;
  std::size_t position = fIn.find(extansion);
  fOut = fIn.substr(0,position) + "_TRANSFORMED.vtk";

  SimplexMesh simplexMesh = SimplexMesh();

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(fIn);
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

  std::cout << "OUTPUT FILE: " << fOut << std::endl;
  gmds::ISimplexMeshIOService ioService2(&simplexMesh);
  VTKWriter writer2(&ioService2);
  writer2.setCellOptions(gmds::N|gmds::R);
  writer2.setDataOptions(gmds::N);
  writer2.write(fOut);
  return 0;
}
/*----------------------------------------------------------------------------*/
