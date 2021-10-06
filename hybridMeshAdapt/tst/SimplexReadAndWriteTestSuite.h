/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <fstream>
#include <bitset>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/PointSmoothing.h>
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/EdgeCollapse.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/EdgeInsertion.h>
#include <gmds/hybridMeshAdapt/DelaunayPointInsertion.h>
#include <gmds/hybridMeshAdapt/Octree.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
//#include <chrono>
//#include "sys/types.h"
//#include "sys/sysinfo.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
TEST(SimplexReadAndWriteTestClass, test_write_simple_triangle)
{
  gmds::hybrid::SimplexMesh simplexMesh;

  simplexMesh.addNode(math::Point(0.0, 0.0, 0.0)); //0
  simplexMesh.addNode(math::Point(1.0, 0.0, 0.0)); //1
  simplexMesh.addNode(math::Point(0.0, 0.0, 1.0)); //2

  simplexMesh.addTriangle(0, 1, 2);

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F);
  vtkWriter.write("simple_triangle.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexReadAndWriteTestClass, test)
{
  gmds::hybrid::SimplexMesh simplexMesh;

  simplexMesh.addNode(math::Point(-10.6112, -5, -10.0744)); //0


  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F);
  vtkWriter.write("simple_.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexReadAndWriteTestClass, test_write_simple_triangles)
{
  gmds::hybrid::SimplexMesh simplexMesh;

  simplexMesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
  simplexMesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
  simplexMesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
  simplexMesh.addNode(math::Point(0.0, 1.0, 0.0)); // Node 3
  simplexMesh.addNode(math::Point(0.0, -1.0, 0.0)); // Node 4
  simplexMesh.addNode(math::Point(0.0, 0.0, -1.0)); // Node 5
  simplexMesh.addNode(math::Point(1.0, 0.0, 1.0)); // Node 6

  simplexMesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
  simplexMesh.addTetraedre(0, 1 , 2, 4); // Tetra 1
  simplexMesh.addTetraedre(0, 1 , 4, 5); // Tetra 2
  simplexMesh.addTetraedre(0, 1 , 3, 5); // Tetra 3

  // the triangle must be oriented with a face of tetraedron
  simplexMesh.addTriangle(0, 1, 2); //-1
  simplexMesh.addTriangle(0, 1, 5); //-2
  simplexMesh.addTriangle(0, 1, 3); //-
  simplexMesh.addTriangle(1, 4, 2); //-4
  simplexMesh.addTriangle(1, 2, 3); //-5

  simplexMesh.addTetraedre(1, 3 , 2, 6); // Tetra 4*/
  Octree oc(&simplexMesh, 1);
  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriter.write("simple_triangles.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, test_volume)
{
  gmds::hybrid::SimplexMesh simplexMesh;
  simplexMesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
  simplexMesh.addNode(math::Point(2.0, 0.0, 0.0)); // Node 1
  simplexMesh.addNode(math::Point(0.0, 0.0, 2.0)); // Node 2
  simplexMesh.addNode(math::Point(0.0, 2.0, 0.0)); // Node 3

  simplexMesh.addTetraedre(0, 1 , 2, 3); // Tetra 0

  SimplicesCell cell = SimplicesCell(&simplexMesh, 0);
  double volumeCell = cell.getVolumeOfCell();
}

TEST(SimplexMeshTestClass, test_hexa_generation_on_modelCAD5)
{
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  //std::string vtk_file = dir+"/halfCylinder/HalfCylinder.vtk";
  //std::string vtk_file = dir+"/Dome/Dome.vtk";
  //std::string vtk_file = dir+"/ModeleEnS/ModeleEnS.vtk";
  //std::string vtk_file = dir+"/CylindreTordu/CylindreTordu.vtk";
  //std::string vtk_file = dir+"/ModeleCAD1/ModeleCAD1.vtk";
  //std::string vtk_file = dir+"/ModeleCAD2/ModeleCAD2.vtk";
  //std::string vtk_file = dir+"/ModeleCAD3/ModeleCAD3.vtk";
  //std::string vtk_mesh = dir+"/ModeleCAD6/ModeleCAD6.vtk";
  //std::string vtk_mesh = dir+"/ModeleCAD7/ModeleCAD7.vtk";
  std::string vtk_mesh = dir+"/ModeleCAD5/ModeleCAD5.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();

  //return;
  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);
  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);
  //////////////////////////////////////////////////////////////////////////////


  //std::string vtk_node = dir+"/ModeleCAD5/ModeleCAD5_generatedPoints_2585.vtk";
  //std::string vtk_node = dir+"/ModeleCAD6/ModeleCAD6_generatedPoints_4601.vtk";
  //std::string vtk_node = dir+"/ModeleCAD7/ModeleCAD7_generatedPoints_2128.vtk";
  //std::string vtk_node = dir+"/ModeleCAD5/ModeleCAD5_generatedPoints_8396.vtk";
  //std::string vtk_node = dir+"/ModeleCAD5/ModeleCAD5_edge_generatedPoints_8396.vtk";
  std::string vtk_node = dir+"/ModeleCAD5/ModeleCAD5_edge_generatedPoints_15701.vtk";

  SimplexMesh simplexNodes = SimplexMesh();
  gmds::ISimplexMeshIOService ioServiceNodes(&simplexNodes);
  gmds::VTKReader vtkReaderNodes(&ioServiceNodes);
  vtkReaderNodes.setCellOptions(gmds::R|gmds::N|gmds::F);
  vtkReaderNodes.setDataOptions(gmds::N);
  vtkReaderNodes.read(vtk_node);
  Variable<int>* BND_CURVE_COLOR_NODES   = simplexNodes.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR_NODES = simplexNodes.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  const gmds::BitVector& nodesToAddIds = simplexNodes.getBitVectorNodes();

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());

  std::vector<std::vector<TInt>> edges{};
  const gmds::BitVector triIdx = simplexNodes.getBitVectorTri();
  for(unsigned int tri = 1 ; tri < triIdx.capacity() ; tri++)
  {
    if(triIdx[tri] != 0)
    {
      std::vector<TInt> nodes = SimplicesTriangle(&simplexNodes,tri).getNodes();
      std::vector<TInt> edge{nodes[0], nodes[1]};
      edges.push_back(edge);
    }
  }


  for(unsigned int idx = 0 ; idx < nodesToAddIds.capacity() ; idx++)
  {
    if(nodesToAddIds[idx] != 0)
    {
      const gmds::BitVector & nodesIds = simplexMesh.getBitVectorNodes();
      math::Point point = SimplicesNode(&simplexNodes, idx).getCoords();

      bool alreadyAdd = false;
      std::vector<TSimplexID> tetraContenaingPt{};
      TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);

      if(!alreadyAdd)
      {
        for(auto & edge : edges)
        {
          if(edge[0] == idx){edge[0] = node;}
          if(edge[1] == idx){edge[1] = node;}
        }
        if((*BND_CURVE_COLOR_NODES)[idx] != 0) {BND_CURVE_COLOR->set(node, (*BND_CURVE_COLOR_NODES)[idx]);}
        else if((*BND_SURFACE_COLOR_NODES)[idx] != 0) {BND_SURFACE_COLOR->set(node, (*BND_SURFACE_COLOR_NODES)[idx]);}

        //forcage de node a cause d'une mauvaise labelisation dans fram3D pour le model M5 et les 8396 points 
        /*if(node == 10506)
        {
          BND_SURFACE_COLOR->set(node, 9);
        }
        else if(node == 8448)
        {
          BND_SURFACE_COLOR->set(node, 5);
        }*/

        simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
        bool status = false;
        std::vector<TSimplexID> deletedSimplex{};


        DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt, status, nodesAdded, deletedSimplex);

        if(status)
        {
          if(nodesAdded.capacity() != nodesIds.capacity())
          {
            nodesAdded.resize(nodesIds.capacity());
          }
          nodesAdded.assign(node);
        }
      }
      else
      {
        if(nodesAdded.capacity() != nodesIds.capacity())
        {
          nodesAdded.resize(nodesIds.capacity());
        }
        nodesAdded.assign(node);
      }
    }
  }


  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  //vtkWriterDI.write("ModeleCAD5_DI_8396.vtk");
  vtkWriterDI.write("ModeleCAD5_DI_15701.vtk");

  simplexMesh.edgesRemove(nodesAdded);

  gmds::VTKWriter vtkWriterER(&ioService);
  vtkWriterER.setCellOptions(gmds::N|gmds::R);
  vtkWriterER.setDataOptions(gmds::N|gmds::R);
  //vtkWriterER.write("Dome_ER_2922.vtk");
  //vtkWriterER.write("HalfCylinder_ER_6683.vtk");
  //vtkWriterER.write("ModeleEnS_ER_4466.vtk");
  //vtkWriterER.write("CylindreTordu_ER_5686.vtk");
  //vtkWriterER.write("ModeleCAD_ER_5686.vtk");
  //vtkWriterER.write("ModeleCAD2_ER_5199.vtk");
  //vtkWriterER.write("ModeleCAD3_ER_5853.vtk");
  //vtkWriterER.write("ModeleCAD5_ER_2585.vtk");
  //vtkWriterER.write("ModeleCAD7_ER_2128_SCV.vtk");
  //vtkWriterER.write("ModeleCAD5_ER_8396_SCV.vtk");
  vtkWriterER.write("ModeleCAD5_ER_15701_SCV.vtk");
  std::cout << "vtkWriterER done !" << std::endl;

  std::cout << "buildEdges start !" << std::endl;
  simplexMesh.buildEdges(edges, nodesAdded);
  std::cout << "buildEdges done !" << std::endl;

  gmds::VTKWriter vtkWriterEB(&ioService);
  vtkWriterEB.setCellOptions(gmds::N|gmds::R);
  vtkWriterEB.setDataOptions(gmds::N|gmds::R);
  //vtkWriterEB.write("ModeleCAD5_EB_8396_SCV.vtk");
  vtkWriterEB.write("ModeleCAD5_EB_15701_SCV.vtk");

}
