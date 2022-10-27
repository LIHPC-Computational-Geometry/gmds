/*----------------------------------------------------------------------------*/
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/DelaunayPointInsertion.h>
#include <gmds/hybridMeshAdapt/Octree.h>
#include <gmds/hybridMeshAdapt/MetricFFPointgeneration.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{

    std::cout << "==== METRIC BASED ADAPTATION  ====" << std::endl;
    std::string fIn, fOut;
    if(argc != 2)
    {
        throw gmds::GMDSException("MISSING PARAMETERS : <mesh_in> ");
    }
    fIn = std::string(argv[1]);
    if (fIn.find('.vtk') == std::string::npos) {
      throw gmds::GMDSException("<mesh_in> NOT A .vtk FILE");
    }

    //==================================================================
    // MESH READING
    //==================================================================
    std::cout << "Reading " << std::endl;
    std::string extansion(".vtk");
    std::size_t position = fIn.find(extansion);
    fOut = fIn.substr(0,position) + "_Points_Generated.vtk";
    std::cout << "INPUT FILE: " << fIn << std::endl;
    std::cout << "OUTPUT FILE: " << fOut << std::endl;

    SimplexMesh simplexMesh = SimplexMesh();
    gmds::ISimplexMeshIOService ioService(&simplexMesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::R|gmds::N);
    vtkReader.setDataOptions(gmds::N);
    vtkReader.read(fIn);
    simplexMesh.buildAdjInfoGlobal();
    simplexMesh.initializeEdgeStructure();
    simplexMesh.buildSimplexHull();


    Octree oc(&simplexMesh, 50);
    simplexMesh.setOctree(&oc);
    //Variable<int>* BND_METRIC = simplexMesh.getVariable<Eigen::Matrix3d,SimplicesNode>("metric");

    //==================================================================
    // MODIFICATION OF THE INPUT MESH'S METRIC
    //==================================================================
    Variable<Eigen::Matrix3d>* metricNode = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
    metricNode->setValuesTo(m);
    const gmds::BitVector& meshNode = simplexMesh.getBitVectorNodes();
    for(unsigned int nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
    {
      if(meshNode[nodeId] == 1)
      {
        simplexMesh.setAnalyticMetric(nodeId);
      }
    }
    //////////////////////////////////////////////////////////////////////////////
    MetricFFPointgeneration p(&simplexMesh);
    p.execute();


    //std::cout << "MESH VALIDITY CHECK" << std::endl;
    //simplexMesh.checkMesh();

    gmds::VTKWriter vtkWriterMA(&ioService);
    vtkWriterMA.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMA.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMA.write(fOut);

    return 0;
}

/*----------------------------------------------------------------------------*/
