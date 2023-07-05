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
    if(!(argc == 3 || argc == 4))
    {
        throw gmds::GMDSException("MISSING PARAMETERS : <mesh_in> <interpolation_value> <size>");
    }
    fIn = std::string(argv[1]);
    if (fIn.find('.vtk') == std::string::npos) {
      throw gmds::GMDSException("<mesh_in> NOT A .vtk FILE");
    }
    if(std::atof(argv[2]) > 1.0 || std::atof(argv[2]) < 0.0) {
      throw gmds::GMDSException("<interpolation_value> should be higher that 0 and lower than 1 !");
    }
    std::cout << "Reading " << std::endl;
    std::string extansion(".vtk");
    std::cout << "INPUT FILE: " << fIn << std::endl;

    SimplexMesh simplexMesh = SimplexMesh();
    gmds::ISimplexMeshIOService ioService(&simplexMesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::R|gmds::N);
    vtkReader.setDataOptions(gmds::N);
    vtkReader.read(fIn);
    simplexMesh.buildAdjInfoGlobal();
    std::cout << "Adjacent info built" << std::endl;
    simplexMesh.initializeEdgeStructure();
    std::cout << "Edge structure built" << std::endl;
    simplexMesh.buildSimplexHull();
    std::cout << "Triangle hull built" << std::endl;
    simplexMesh.setSurfacesAndCurvesIndx();
    std::cout << "Surface and Curves indices set" << std::endl;

    //double minSizeEdge = simplexMesh.findMinSizeEdgeSurface();
    double meanSizeEdge = simplexMesh.findMeanSizeEdgeSurface();

    const std::vector<double> min_maxXYZ = simplexMesh.getMinMaxCoord();
    std::cout << "minimum size Edge on surfaces is : " << meanSizeEdge << std::endl;
    //vector of lambda that capture metric and ffield
    std::vector<std::function<std::vector<double>(double x_)>> metricXYZ_functors{};
    if(argc == 4){
      double size = std::stof(argv[3]);
      metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{size, size, size}; });
      //metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{0.1, 0.1, 0.1};});
    }
    else{
      std::cout << "size 0.25 * meanSizeEdge-> " << 0.25 * meanSizeEdge << std::endl;
      std::cout << "size 0.5 * meanSizeEdge -> " << 0.5 * meanSizeEdge << std::endl;
      std::cout << "size 0.75 * meanSizeEdge -> " << 0.75 * meanSizeEdge << std::endl;
      std::cout << "size meanSizeEdge -> " <<  meanSizeEdge << std::endl;
      metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{0.5*meanSizeEdge, 0.5*meanSizeEdge, 0.5*meanSizeEdge}; });
      //metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{0.1*x_ + 0.0125*(1.0 - x_), 0.1, 0.1};});
      //metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{0.05, 0.05, 0.05};});
      //metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{0.5*meanSizeEdge, 0.5*meanSizeEdge, 0.5*meanSizeEdge}; });
      //metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{0.5*meanSizeEdge, 0.25*meanSizeEdge, 0.5*meanSizeEdge}; });
      //metricXYZ_functors.push_back([&](double x_) { return std::vector<double>{0.5*meanSizeEdge, 0.35*meanSizeEdge, 0.25*meanSizeEdge}; });
    }


    std::vector<std::function<std::vector<math::Vector3d>()>> frameXYZ_functor{};
    frameXYZ_functor.push_back([] { return std::vector<math::Vector3d>{}; });
    const double alpha0 = 0.0;
    const double alpha1 = M_PI / 8.0;
    const double alpha2 = M_PI / 4.0;
    //frameXYZ_functor.push_back([&] { return std::vector<math::Vector3d>{ math::Vector3d({1.0, 0.0, 0.0}),  math::Vector3d({0.0, cos(alpha0), sin(alpha0)}),  math::Vector3d({0.0, -sin(alpha0), cos(alpha0)})}; });
    //frameXYZ_functor.push_back([&] { return std::vector<math::Vector3d>{ math::Vector3d({1.0, 0.0, 0.0}),  math::Vector3d({0.0, cos(alpha1), sin(alpha1)}),  math::Vector3d({0.0, -sin(alpha1), cos(alpha1)})}; });
    //frameXYZ_functor.push_back([&] { return std::vector<math::Vector3d>{ math::Vector3d({1.0, 0.0, 0.0}),  math::Vector3d({0.0, cos(alpha2), sin(alpha2)}),  math::Vector3d({0.0, -sin(alpha2), cos(alpha2)})}; });


    //==================================================================
    // MODIFICATION OF THE INPUT MESH'S METRIC
    //==================================================================
    unsigned int cpt = 0;
    const double d = std::atof(argv[2]);
    std::cout << "interpolation factor -> " << d << std::endl;
    for(auto const ff : frameXYZ_functor)
    {
      for(auto const metric : metricXYZ_functors)
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
          if(meshNode[nodeId] != 0)
          {
            //simplexMesh.setAnalyticMetric(nodeId, simplexMesh.getOctree());
            Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
            double x = SimplicesNode(&simplexMesh, nodeId).getCoords().X();
            m(0, 0) = 1.0/(metric(x)[0]*metric(x)[0]); m(1, 1) = 1.0/(metric(x)[1]*metric(x)[1]) ; m(2, 2) = 1.0/(metric(x)[2]*metric(x)[2]);
            simplexMesh.setAnalyticMetric(nodeId, m);
            simplexMesh.setFrames(nodeId, ff());
          }
        }

        //////////////////////////////////////////////////////////////////////////////
        std::cout << "FRONTAL ALGO STARTING ..." << std::endl;
        std::string name =+ "HEX_" + std::to_string(cpt);
        std::cout << " name -> " << name << std::endl;

        MetricFFPointgeneration p(&simplexMesh, name, d);
        p.execute();
        ++cpt;
        std::cout << std::endl;
        std::cout << std::endl;
      }
    }

    //std::cout << "MESH VALIDITY CHECK" << std::endl;
    //simplexMesh.checkMesh();
    return 0;
}

/*----------------------------------------------------------------------------*/
