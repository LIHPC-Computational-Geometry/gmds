#include "gmds/hybridMeshAdapt/DelaunayPointInsertion.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
/******************************************************************/
#include <chrono>
/******************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesCell;
using namespace simplicesTriangle;
using namespace simplicesNode;
using namespace math;
/******************************************************************/
DelaunayPointInsertion::DelaunayPointInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNode, const CriterionRAIS& criterion,
                                              std::vector<TSimplexID>& initialCavity, bool& status, const gmds::BitVector& markedNodes,std::vector<TSimplexID>& deletedSimplex,
                                              const std::multimap<TInt, std::pair<TInt, TInt>>& facesAlreadyBuilt, std::vector<TSimplexID> markedSimplex)
{
  if(simplexMesh != nullptr)
  {
    /*Si simpliceNode n'est pas a linterrieur de simplexMeshon ne fait rien*/
    std::vector<TSimplexID> initCavity;
    int border = std::numeric_limits<int>::min();
    gmds::BitVector cellAlreadyAdd(simplexMesh->tetCapacity());
    gmds::BitVector triangleAlreadyAdd(simplexMesh->triCapacity());

    if(initialCavity.size() == 0)
    {
      SimplicesNode::nodeNeighborInfo nodeInfo;
      if(!simplexMesh->checkSimplicesContenaing(simpliceNode.getCoords(), initCavity))
      {
        std::cout << "SIMPLICE NODE : " << simpliceNode << "IS NOT ON ANY TETRAEDRE" << std::endl;
        return;
      }
    }
    else
    {
      initCavity = initialCavity;
    }
    bool flag = true;
    std::vector<TSimplexID> cavity = initCavity;
    CavityReduction cavReduction;
    cavReduction.cavInitial = initCavity;


    if(flag)
    {
      for(auto const & simplex : cavity)
      {
        if(simplex >= 0){cellAlreadyAdd.assign(simplex);}
        else{triangleAlreadyAdd.assign(-simplex);}
      }
      for(unsigned int iter = 0; iter < cavity.size() ; iter++) // je suis obligé d'utiliser ce type de loop car je push_back dans le vector cavity dans le loop
      {
        TSimplexID tet = cavity[iter];
        if(!(tet < 0)) // is not triangle
        {
          const std::vector<TSimplexID>&& adjaSimplex = SimplicesCell(simplexMesh, tet).adjacentTetra();
          for(auto const & simplex : adjaSimplex)
          {
            if(simplex != border)
            {
              if(simplex < 0)
              {
                if(triangleAlreadyAdd[-simplex] == 0)
                {
                  cavity.push_back(simplex);
                  triangleAlreadyAdd.assign(-simplex);
                }
              }
              else
              {
                if(cellAlreadyAdd[simplex] == 0)
                {
                  if(isNodeInCircumSphere(simplexMesh, simpliceNode, simplex))
                  {
                    cavity.push_back(simplex);
                  }
                  cellAlreadyAdd.assign(simplex);
                }
              }
            }
          }
        }
      }
      std::vector<TSimplexID> cellsCreated{};
      PointInsertion pi(simplexMesh, simpliceNode, criterion, status, cavity, markedNodes, deletedSimplex, facesAlreadyBuilt, cellsCreated, markedSimplex/*, cavReduction*/);
    }
  }
  else
  {
    //exception TODO
    std::cout << "simplexMesh == nullptr in delaunaypointInsertion.cpp" << std::endl;
  }
}


bool DelaunayPointInsertion::isNodeInCircumSphere(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& node, const TSimplexID& simplex)
{
  bool flag = false;
  std::vector<TInt>&& nodes = SimplicesCell(simplexMesh, simplex).getNodes();

  if(nodes.size() == 4)
  {
    SimplicesNode S0 = SimplicesNode(simplexMesh, nodes[0]);
    SimplicesNode S1 = SimplicesNode(simplexMesh, nodes[1]);
    SimplicesNode S2 = SimplicesNode(simplexMesh, nodes[2]);
    SimplicesNode S3 = SimplicesNode(simplexMesh, nodes[3]);

    Vector3d v1 = S0.getCoords(); Vector3d v2 = S1.getCoords(); Vector3d v3 = S2.getCoords(); Vector3d v4 = S3.getCoords();
    double x1 = v1.X(); double y1 = v1.Y(); double z1 = v1.Z();
    double x2 = v2.X(); double y2 = v2.Y(); double z2 = v2.Z();
    double x3 = v3.X(); double y3 = v3.Y(); double z3 = v3.Z();
    double x4 = v4.X(); double y4 = v4.Y(); double z4 = v4.Z();

    double Ox = 0.0;
    double Oy = 0.0;
    double Oz = 0.0;
    double det = 0.0;
    double d12, d13, d14 = 0.0;
    double x12, x13, x14 = 0.0;
    double y12, y13, y14 = 0.0;
    double z12, z13, z14 = 0.0;

    Eigen::Matrix3d m0;
    gmds::Variable<Eigen::Matrix3d>* var =  simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("metric");
    if(var != nullptr)
    {
      m0 = var->value(nodes[0]);
      double m00 = m0(0, 0); double m01 = m0(0, 1); double m02 = m0(0, 2);
      double m10 = m0(1, 0); double m11 = m0(1, 1); double m12 = m0(1, 2);
      double m20 = m0(2, 0); double m21 = m0(2, 1); double m22 = m0(2, 2);


      x12 = 2.0*(m00*(x2 - x1) + m01*(y2 - y1) + m02*(z2 - z1));
      y12 = 2.0*(m11*(y2 - y1) + m01*(x2 - x1) + m12*(z2 - z1));
      z12 = 2.0*(m22*(z2 - z1) + m02*(x2 - x1) + m12*(y2 - y1));

      x13 = 2.0*(m00*(x3 - x1) + m01*(y3 - y1) + m02*(z3 - z1));
      y13 = 2.0*(m11*(y3 - y1) + m01*(x3 - x1) + m12*(z3 - z1));
      z13 = 2.0*(m22*(z3 - z1) + m02*(x3 - x1) + m12*(y3 - y1));

      x14 = 2.0*(m00*(x4 - x1) + m01*(y4 - y1) + m02*(z4 - z1));
      y14 = 2.0*(m11*(y4 - y1) + m01*(x4 - x1) + m12*(z4 - z1));
      z14 = 2.0*(m22*(z4 - z1) + m02*(x4 - x1) + m12*(y4 - y1));

      d12 = m00*x2*x2 + m11*y2*y2 + m22*z2*z2 + 2.0*m01*x2*y2 + 2.0*m02*x2*z2 + 2.0*m12*y2*z2 - (m00*x1*x1 + m11*y1*y1 + m22*z1*z1 + 2.0*m01*x1*y1 + 2.0*m02*x1*z1 + 2.0*m12*y1*z1);
      d13 = m00*x3*x3 + m11*y3*y3 + m22*z3*z3 + 2.0*m01*x3*y3 + 2.0*m02*x3*z3 + 2.0*m12*y3*z3 - (m00*x1*x1 + m11*y1*y1 + m22*z1*z1 + 2.0*m01*x1*y1 + 2.0*m02*x1*z1 + 2.0*m12*y1*z1);
      d14 = m00*x4*x4 + m11*y4*y4 + m22*z4*z4 + 2.0*m01*x4*y4 + 2.0*m02*x4*z4 + 2.0*m12*y4*z4 - (m00*x1*x1 + m11*y1*y1 + m22*z1*z1 + 2.0*m01*x1*y1 + 2.0*m02*x1*z1 + 2.0*m12*y1*z1);


    }
    else
    {
      x12 = 2.0*(x2 - x1); y12 = 2.0*(y2 - y1); z12 = 2.0*(z2 - z1);
      x13 = 2.0*(x3 - x1); y13 = 2.0*(y3 - y1); z13 = 2.0*(z3 - z1);
      x14 = 2.0*(x4 - x1); y14 = 2.0*(y4 - y1); z14 = 2.0*(z4 - z1);

      d12 = (x2*x2 + y2*y2 + z2*z2) - (x1*x1 + y1*y1 + z1*z1);
      d13 = (x3*x3 + y3*y3 + z3*z3) - (x1*x1 + y1*y1 + z1*z1);
      d14 = (x4*x4 + y4*y4 + z4*z4) - (x1*x1 + y1*y1 + z1*z1);
    }

    Ox = d12 *(y13*z14 - y14*z13) + d13*(z12*y14 - y12*z14) + d14*(y12*z13 - z12 * y13);
    Oy = d12 *(z13*x14 - z14*x13) + d13*(x12*z14 - z12*x14) + d14*(z12*x13 - x12 * z13);
    Oz = d12 *(x13*y14 - x14*y13) + d13*(y12*x14 - x12*y14) + d14*(x12*y13 - y12 * x13);

    det = x12 * (y13*z14  - y14*z13)  - x13 * (y12*z14 - y14*z12) + x14 * (y12*z13-y13*z12);

    if(det != 0)
    {
        Ox /= det; Oy /= det; Oz /= det;
        Metric<Eigen::Matrix3d> metric0 =  Metric<Eigen::Matrix3d>(m0);
        const math::Vector3d O = math::Vector3d(Ox, Oy, Oz);
        double dist = metric0.metricDist(node.getCoords(), O, metric0);
        double alpha = dist / metric0.metricDist(S0.getCoords(), O, metric0);
        flag = (alpha <= 1.0)? true : false;
    }
    else
    {
      flag = true;
      //std::cout << "adjTet volume is nul " << std::endl;
      //TODO
      //std::cout << "det == 0" << std::endl;
    }

    return flag;
  }
  else
  {
    //TODO
    std::cout << "nodes.size() == 4" << std::endl;
  }

  return flag;
}
