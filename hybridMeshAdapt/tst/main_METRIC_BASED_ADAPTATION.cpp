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
#include <gmds/hybridMeshAdapt/MetricAdaptation.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
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
    fOut = fIn.substr(0,position) + "_ADAPTED.vtk";
    std::cout << "INPUT FILE: " << fIn << std::endl;
    std::cout << "OUTPUT FILE: " << fOut << std::endl;


    //////////////////////////////////////////////
    /*SimplexMesh regularTet = SimplexMesh();
    TInt nodeA = regularTet.addNode(0.0, 0.0, 0.0);
    TInt nodeB = regularTet.addNode(0.0, 1.0 / sqrt(2), 1.0 / sqrt(2));
    TInt nodeC = regularTet.addNode(1.0 / sqrt(2), 0.0, 1.0 / sqrt(2));
    TInt nodeD = regularTet.addNode(1.0 / sqrt(2), 1.0 / sqrt(2), 0.0);
    TSimplexID tet = regularTet.addTetraedre(nodeA, nodeB, nodeC, nodeD);
    Variable<Eigen::Matrix3d>* metric = regularTet.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    Eigen::Matrix3d metricUnit =  Eigen::MatrixXd::Identity(3, 3);
    //m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    metric->setValuesTo(metricUnit);
    std::cout << "tet quality -> " << regularTet.computeQualityElement(tet) << std::endl;
    gmds::ISimplexMeshIOService ioServiceTet(&regularTet);
    gmds::VTKWriter vtkWriterTET(&ioServiceTet);
    vtkWriterTET.setCellOptions(gmds::N|gmds::R);
    vtkWriterTET.write("REGULAR_TET.vtk");
    return 0 ;*/
    //////////////////////////////////////////////


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
    m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
    metricNode->setValuesTo(m);
    const gmds::BitVector& meshNode = simplexMesh.getBitVectorNodes();
    //double metricMin = 0.5;
    //double metricMax = 0.5;
    double meanEdge = 0.0; double minEdge = 0.0; double maxEdge = 0.0;
    simplexMesh.getEdgeSizeInfo(meanEdge, minEdge, maxEdge);
    std::cout << "meanEdge START -> " << meanEdge << std::endl;
    std::cout << "minEdge  START -> " << minEdge << std::endl;
    std::cout << "maxEdge  START -> " << maxEdge << std::endl;

    for(unsigned int nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
    {
      if(meshNode[nodeId] == 1)
      {
        simplexMesh.setAnalyticMetric(nodeId, simplexMesh.getOctree());
      }
    }
    //////////////////////////////////////////////////////////////////////////////
    std::cout << "INITIAL NODE SIZE IN MESH --> " << meshNode.capacity() << std::endl;
    MetricAdaptation mA(&simplexMesh);
    //mA.execute();
    mA.executeCustomMethod();
    std::cout << "NODE SIZE IN MESH AFTER METRIC ADAPTATION--> " << meshNode.capacity() << std::endl;


    std::cout << "MESH VALIDITY CHECK" << std::endl;
    simplexMesh.checkMesh();

    //////
    const gmds::BitVector& tetIDS = simplexMesh.getBitVectorTet();
    std::vector<TSimplexID> v{};
    for(unsigned int tet = 0 ; tet < tetIDS.capacity() ; tet++)
    {
      if(tetIDS[tet] != 0)
      {
        const SimplicesCell cell(&simplexMesh, tet);
        if(cell.isSliver())
        {
          v.push_back(tet);
        }
      }
    }
    simplexMesh.deleteAllSimplicesBut(v);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::R);
    vtkWriter.write("SILVER.vtk");
    //////


    gmds::VTKWriter vtkWriterMA(&ioService);
    vtkWriterMA.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMA.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMA.write(fOut);
}

/*----------------------------------------------------------------------------*/
