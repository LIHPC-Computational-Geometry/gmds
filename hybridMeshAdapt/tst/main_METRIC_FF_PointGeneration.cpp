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
#include <gmds/math/Hexahedron.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace math;
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
    std::cout << "Reading " << std::endl;
    std::string extansion(".vtk");
    std::cout << "INPUT FILE: " << fIn << std::endl;


    //vector of lambda that capture metric and ffield
    std::vector<std::function<std::vector<double>()>> metricXYZ_functors{};
    metricXYZ_functors.push_back([] { return std::vector<double>{0.1, 0.1, 0.1}; });
    metricXYZ_functors.push_back([] { return std::vector<double>{0.2, 0.2, 0.2}; });

    std::vector<std::function<std::vector<math::Vector3d>()>> frameXYZ_functor{};
    frameXYZ_functor.push_back([] { return std::vector<math::Vector3d>{ math::Vector3d({1.0, 0.0, 0.0}),  math::Vector3d({0.0, 1.0, 0.0}),  math::Vector3d({0.0, 0.0, 1.0})}; });
    frameXYZ_functor.push_back([] { return std::vector<math::Vector3d>{ math::Vector3d({sqrt(2.0)/2, 0.0, -sqrt(2.0)/2}) , math::Vector3d({0.0, 1.0, 0.0}) , math::Vector3d({sqrt(2.0)/2, 0.0, sqrt(2.0)/2})}; });
    //==================================================================
    // MODIFICATION OF THE INPUT MESH'S METRIC
    //==================================================================
    std::vector<std::vector<double>> metricXYZ{};
    std::vector<std::vector<math::Vector3d>> frameXYZ{};
    //for(auto const ff : frameXYZ)
    {
      //for(auto const metric : metricXYZ)
      {
        //==================================================================
        // MESH READING
        //==================================================================
        SimplexMesh simplexMesh = SimplexMesh();
        gmds::ISimplexMeshIOService ioService(&simplexMesh);
        gmds::VTKReader vtkReader(&ioService);
        vtkReader.setCellOptions(gmds::R|gmds::N);
        vtkReader.setDataOptions(gmds::N);
        vtkReader.read(fIn);
        simplexMesh.buildAdjInfoGlobal();
        simplexMesh.initializeEdgeStructure();
        simplexMesh.buildSimplexHull();
        simplexMesh.setSurfacesAndCurvesIndx();

        Octree oc(&simplexMesh, 10);
        simplexMesh.setOctree(&oc);


        Variable<Eigen::Matrix3d>* metricNode = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
        Eigen::Matrix3d m = Eigen::MatrixXd::Identity(3, 3);
        metricNode->setValuesTo(m);
        const gmds::BitVector& meshNode = simplexMesh.getBitVectorNodes();

        for(unsigned int nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
        {
          if(meshNode[nodeId] == 1)
          {
            simplexMesh.setAnalyticMetric(nodeId, simplexMesh.getOctree());
          }
        }
        //////////////////////////////////////////////////////////////////////////////
        std::cout << "FRONTAL ALGO STARTING ..." << std::endl;
        MetricFFPointgeneration p(&simplexMesh);
        p.execute();

      }
    }

    //std::cout << "MESH VALIDITY CHECK" << std::endl;
    //simplexMesh.checkMesh();
    return 0;
}

/*----------------------------------------------------------------------------*/
