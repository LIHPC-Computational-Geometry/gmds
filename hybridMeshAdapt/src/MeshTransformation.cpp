#include "gmds/hybridMeshAdapt/MeshTransformation.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace simplicesCell;
using namespace simplicesTriangle;
using namespace operators;
/******************************************************************************/
MeshTransformation::MeshTransformation(SimplexMesh* simplexMesh):
m_simplexMesh(simplexMesh)
{

}
/******************************************************************************/
void MeshTransformation::transformation()
{

  /*std::vector<double> res{};
  Variable<Eigen::Matrix3d>* NODE_METRIC  = m_simplexMesh->getVariable<Eigen::Matrix3d,SimplicesNode>("NODE_METRIC");

  //extraction des coordonnÃ©es X,Y,Z grÃ¢ce aux systeme de ressort pondÃ©rÃ©
  const BitVector& nodeBitVector = m_simplexMesh->getBitVectorNodes();
  const unsigned int nodes_size  = nodeBitVector.size();
  Eigen::VectorXd constraintsU = Eigen::VectorXd::Zero(nodes_size);
  Eigen::VectorXd constraintsV = Eigen::VectorXd::Zero(nodes_size);
  Eigen::VectorXd constraintsW = Eigen::VectorXd::Zero(nodes_size);

  math::Vector3d x = math::Vector3d(1.0, 0.0, 0.0);
  math::Vector3d y = math::Vector3d(0.0, 1.0, 0.0);
  math::Vector3d z = math::Vector3d(0.0, 0.0, 1.0);

  for(unsigned int nodeIter = 0 ; nodeIter < nodeBitVector.capacity() ; nodeIter++)
  {
    if(nodeBitVector[nodeIter] != 0)
    {
      const SimplicesNode node = SimplicesNode(m_simplexMesh, nodeIter);
      std::vector<TInt> neighboorNodes = SimplicesNode(m_simplexMesh, nodeIter).getNeighboorNodes();

      double scalarX = 0.0;
      double scalarY = 0.0;
      double scalarZ = 0.0;


      double mi00 = (*NODE_METRIC)[nodeIter](0,0);
      double mi11 = (*NODE_METRIC)[nodeIter](1,1);
      double mi22 = (*NODE_METRIC)[nodeIter](2,2);

      for(auto const nNode : neighboorNodes)
      {
        double m00 = (*NODE_METRIC)[nNode](0,0);
        double m11 = (*NODE_METRIC)[nNode](1,1);
        double m22 = (*NODE_METRIC)[nNode](2,2);

        const math::Point neighboorNodeCoord = SimplicesNode(m_simplexMesh, nNode).getCoords();
        //math::Vector3d vec = neighboorNodeCoord - SimplicesNode(m_simplexMesh, nodeIter).getCoords();

        Eigen::Matrix3d Mi = (*NODE_METRIC)[nodeIter];
        Eigen::Matrix3d Mj = (*NODE_METRIC)[nNode];
        math::Vector3d vec = neighboorNodeCoord - SimplicesNode(m_simplexMesh, nodeIter).getCoords();

        double alpha_ijX = vec.dot(x) * (0.5*m00 + 0.5*mi00);
        double alpha_ijY = vec.dot(y) * (0.5*m11 + 0.5*mi11);
        double alpha_ijZ = vec.dot(z) * (0.5*m22 + 0.5*mi22);


        std::cout << "alpha_ijX ---> " << alpha_ijX << std::endl;
        std::cout << "alpha_ijY ---> " << alpha_ijY << std::endl;
        std::cout << "alpha_ijZ ---> " << alpha_ijZ << std::endl;

        scalarX += (alpha_ijX) / neighboorNodes.size();
        scalarY += (alpha_ijY) / neighboorNodes.size();
        scalarZ += (alpha_ijZ) / neighboorNodes.size();
      }
      constraintsU[nodeIter] = scalarX;
      constraintsV[nodeIter] = scalarY;
      constraintsW[nodeIter] = scalarZ;
    }
  }

  Eigen::VectorXd u = resolveSystem(constraintsU, 0);
  Eigen::VectorXd v = resolveSystem(constraintsV, 1);
  Eigen::VectorXd w = resolveSystem(constraintsW, 2);



  //apply transformation to each node
  unsigned int cptNode = 0;
  for(unsigned int nodeIter = 0 ; nodeIter < nodeBitVector.capacity() ; nodeIter++)
  {
    if(nodeBitVector[nodeIter] == 1)
    {
      math::Point newPosition = math::Point(u[cptNode], v[cptNode], w[cptNode]);
      m_simplexMesh->changeNodeCoord(nodeIter, newPosition);
      cptNode++;
    }
  }

  const math::Vector3d node105Coord = SimplicesNode(m_simplexMesh, 105).getCoords();
  const math::Vector3d node117Coord = SimplicesNode(m_simplexMesh, 117).getCoords();
  const math::Vector3d vec = node105Coord - node117Coord;
  std::cout << "lenght --> " << vec.norm() << std::endl;*/
  resolveSystemByIterativeGradient();
}
/******************************************************************************/
Eigen::VectorXd MeshTransformation::resolveSystem(const Eigen::VectorXd& constraints, const TInt wichDimension)
{
  //Resolution AX=b (b = constraints)
  const BitVector& nodeBitVector = m_simplexMesh->getBitVectorNodes();
  const unsigned int nodes_size  = nodeBitVector.size();
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  //==================================================================
  // STEP 1 / Filling A
  //==================================================================
  for(unsigned int nodeIter = 0 ; nodeIter < nodeBitVector.capacity() ; nodeIter++)
  {
    if(nodeBitVector[nodeIter] != 0)
    {
      std::vector<TInt> neighboorNodes;// = SimplicesNode(m_simplexMesh, nodeIter).getNeighboorNodes();
      //filling the diagonal matrix with 1.0;
      tripletList.push_back(T(nodeIter,nodeIter,1.0));
      const TInt sizeNeighboorNode = neighboorNodes.size();
      double den = sizeNeighboorNode/*0.0*/;
      /*for(auto const nNode : neighboorNodes)
      {
        const math::Point neighboorNodeCoord = SimplicesNode(m_simplexMesh, nNode).getCoords();
        math::Vector3d vec = neighboorNodeCoord - SimplicesNode(m_simplexMesh, nodeIter).getCoords();
        den      +=  vec.norm();
      }*/
      for(auto const nNode : neighboorNodes)
      {
        math::Point neighboorNodeCoord = SimplicesNode(m_simplexMesh, nNode).getCoords();
        {
          //lenght weight coefficient
          //math::Vector3d vec = neighboorNodeCoord - SimplicesNode(m_simplexMesh, nodeIter).getCoords();
          //double v_ij        =  - vec.norm() / den;
          //uniform coefficient
          double v_ij = - 1.0 / sizeNeighboorNode;
          tripletList.push_back(T(nodeIter,nNode,v_ij));
        }
      }
    }
  }

  Eigen::SparseMatrix<double> A(nodes_size,nodes_size);
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  //std::cout << Eigen::MatrixXd(A) << std::endl;
  A.makeCompressed();
  //==================================================================
  // STEP 2 / System Solving
  //==================================================================
  //Cholesky methode (wrong because A is not orthonormal)
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > cholA(A);
  Eigen::VectorXd X = cholA.solve(constraints);

  //LU methode
  /*Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  X = solver.solve(constraints);*/

  return X;
}
/******************************************************************************/
void MeshTransformation::resolveSystemByIterativeGradient()
{
  Variable<Eigen::Matrix3d>* NODE_METRIC   = m_simplexMesh->getVariable<Eigen::Matrix3d,SimplicesNode>("NODE_METRIC");
  Variable<int>* BND_VERTEX_COLOR          = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR           = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR         = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  double epsilon = 10E-4;
  double step    = 10E-3;
  unsigned int iterationMax = 3000;
  const BitVector& nodeBitVector = m_simplexMesh->getBitVectorNodes();
  std::vector<math::Point> newsPositions(nodeBitVector.capacity());
  std::vector<double> dE_vecs(nodeBitVector.capacity());

  for(unsigned int cpt = 0 ; cpt < iterationMax ; cpt++)
  {
    std::cout << "cpt --> " << cpt << std::endl;
    for(unsigned int nodeIter = 0 ; nodeIter < nodeBitVector.capacity() ; nodeIter++)
    {
      if(nodeBitVector[nodeIter] != 0)
      {
        std::vector<TInt> neighboorNodes;// = SimplicesNode(m_simplexMesh, nodeIter).getNeighboorNodes();
        Eigen::Matrix3d Mi = (*NODE_METRIC)[nodeIter];
        Eigen::Vector3d dE = Eigen::Vector3d(0.0, 0.0, 0.0);

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int dimNode = ((*BND_VERTEX_COLOR)[nodeIter] != 0)? SimplexMesh::topo::CORNER : ((*BND_CURVE_COLOR)[nodeIter] != 0)? SimplexMesh::topo::RIDGE : ((*BND_SURFACE_COLOR)[nodeIter] != 0)? SimplexMesh::topo::SURFACE : SimplexMesh::topo::VOLUME;
        unsigned int indexNode;

        if     (dimNode == SimplexMesh::topo::CORNER ){  indexNode =  (*BND_VERTEX_COLOR )[nodeIter] ; }
        else if(dimNode == SimplexMesh::topo::RIDGE  ){  indexNode =  (*BND_CURVE_COLOR  )[nodeIter] ; }
        else if(dimNode == SimplexMesh::topo::SURFACE){  indexNode =  (*BND_SURFACE_COLOR)[nodeIter] ; }
        else                                          {  indexNode =  0;}

        //sorting the node by nature
        std::vector<TInt> sortNeighboorNodes;
        for(auto const nNode : neighboorNodes)
        {
          if(dimNode == SimplexMesh::topo::SURFACE && ((*BND_VERTEX_COLOR)[nNode] != 0 || (*BND_CURVE_COLOR)[nNode] != 0 || (*BND_SURFACE_COLOR)[nNode] == indexNode))
          {
            sortNeighboorNodes.push_back(nNode);
          }
        }
        if(sortNeighboorNodes.size() == 0)
        {
          std::copy(neighboorNodes.begin(), neighboorNodes.end(), std::back_inserter(sortNeighboorNodes));
        }

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        double minLenght = 0.1;
        double mi00 = std::fabs(Mi(0,0));
        double mi11 = std::fabs(Mi(1,1));
        double mi22 = std::fabs(Mi(2,2));

        for(auto const nNode : sortNeighboorNodes)
        {
          const math::Point neighboorNodeCoord = SimplicesNode(m_simplexMesh, nNode).getCoords();
          const math::Point currentNodeCoord  = SimplicesNode(m_simplexMesh, nodeIter).getCoords();

          Eigen::Matrix3d Mj = (*NODE_METRIC)[nNode];
          double mj00 = std::fabs(Mj(0,0));
          double mj11 = std::fabs(Mj(1,1));
          double mj22 = std::fabs(Mj(2,2));


          Eigen::Vector3d vec = Eigen::Vector3d(neighboorNodeCoord.X(), neighboorNodeCoord.Y(), neighboorNodeCoord.Z()) - Eigen::Vector3d(currentNodeCoord.X(), currentNodeCoord.Y(), currentNodeCoord.Z());
          Eigen::Vector3d lij = Eigen::Vector3d(vec.x() / (0.5 * mi00 + 0.5*mj00), vec.y() / (0.5 * mi11 + 0.5*mj11), vec.z() / (0.5 * mi22 + 0.5*mj22));
          const double  lij_Norme = lij.norm();
          double alpha_ij = 1.0;//0.5 * sqrt(vec.dot(Mi * vec)) + 0.5 * sqrt(vec.dot(Mj * vec));

          dE = vec  * (1.0 - lij_Norme / vec.norm()) + dE;
        }

        math::Point vec = math::Point(dE.x(), dE.y(), dE.z());
        math::Point newPosition = SimplicesNode(m_simplexMesh, nodeIter).getCoords() - step * vec;
        const std::vector<TSimplexID>&& ball = SimplicesNode(m_simplexMesh, nodeIter).ballOf();
        bool insideACell = false;
        double epsilonVol = 0.0;

        if(dimNode == SimplexMesh::topo::VOLUME)
        {
          for(auto const simplex : ball)
          {
            if(simplex >= 0)
            {
              double u = SimplicesCell(m_simplexMesh, simplex).signedBarycentric(0,newPosition);
              double v = SimplicesCell(m_simplexMesh, simplex).signedBarycentric(1,newPosition);
              double w = SimplicesCell(m_simplexMesh, simplex).signedBarycentric(2,newPosition);
              double t = SimplicesCell(m_simplexMesh, simplex).signedBarycentric(3,newPosition);

              if(u >= epsilonVol && v >= epsilonVol && w >= epsilonVol && t >= epsilonVol)
              {
                insideACell = true;
                break;
              }
            }
          }

          if(insideACell == false)
          {
            dE = Eigen::Vector3d(0.0, 0.0, 0.0);
          }
        }

        newsPositions[nodeIter] = newPosition;
        dE_vecs[nodeIter] = dE.norm();
        if(dE.norm() > epsilon)
        {
            //m_simplexMesh->changeNodeCoord(nodeIter, newPosition);
        }
      }
    }

    gmds::ISimplexMeshIOService ioService(m_simplexMesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    //std::string name = "TRANSFORMATION_CYLINDRETORDU_FINAL_" + std::to_string(cpt) + ".vtk";
    //std::string name = "TRANSFORMATION_CUBE_FINAL_" + std::to_string(cpt) + ".vtk";
    std::string name = "TRANSFORMATION_MODELECAD1_FINAL_" + std::to_string(cpt) + ".vtk";
    vtkWriter.write(name);
  }
}
/******************************************************************************/
Eigen::Matrix3d MeshTransformation::interpolateMetric(const math::Point& currentNode, const math::Point& nodeStart, const math::Point& nodeEnd, const Eigen::Matrix3d& Mi, const Eigen::Matrix3d& Mj)
{
  const math::Vector3d vec = nodeStart - nodeEnd;
  const math::Vector3d vec0 = nodeStart - currentNode;
  const math::Vector3d vec1 = nodeEnd - currentNode;
  const double lenght = vec.norm();
  const double  bar0  = vec0.norm() / lenght;
  const double  bar1  = vec1.norm() / lenght;

  Eigen::Matrix3d M0;
  Eigen::Matrix3d M1;
  M0 << exp(bar0 * log(Mi(0,0))), 0.0, 0.0, 0.0, exp(bar0 * log(Mi(1,1))), 0.0, 0.0, 0.0, exp(bar0 * log(Mi(2,2)));
  M1 << exp(bar1 * log(Mj(0,0))), 0.0, 0.0, 0.0, exp(bar1 * log(Mj(1,1))), 0.0, 0.0, 0.0, exp(bar1 * log(Mj(2,2)));

  Eigen::Matrix3d interpolateMetric = M0 * M1;
  return interpolateMetric;
}
/*
void MeshTransformation::resolveSystemByIterativeGradient()
{
  Variable<Eigen::Matrix3d>* NODE_METRIC   = m_simplexMesh->getVariable<Eigen::Matrix3d,SimplicesNode>("NODE_METRIC");
  Variable<int>* BND_VERTEX_COLOR          = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR           = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR         = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  double epsilon = 10E-3;
  double step    = 10E-4;
  unsigned int iterationMax = 200;
  const BitVector& nodeBitVector = m_simplexMesh->getBitVectorNodes();

  for(unsigned int cpt = 0 ; cpt < iterationMax ; cpt++)
  {
    for(unsigned int nodeIter = 0 ; nodeIter < nodeBitVector.capacity() ; nodeIter++)
    {
      if(nodeBitVector[nodeIter] != 0)
      {
        const math::Point currentNodeCoord = SimplicesNode(m_simplexMesh, nodeIter).getCoords();
        std::vector<TInt> neighboorNodes = SimplicesNode(m_simplexMesh, nodeIter).getNeighboorNodes();
        Eigen::Matrix3d Mi = (*NODE_METRIC)[nodeIter];
        Eigen::Vector3d dE = Eigen::Vector3d(0.0, 0.0, 0.0);

        const math::Vector3d x = math::Vector3d(1.0, 0.0, 0.0);
        const math::Vector3d y = math::Vector3d(0.0, 1.0, 0.0);
        const math::Vector3d z = math::Vector3d(0.0, 0.0, 1.0);

        const double a = Mi(0,0);
        const double b = Mi(1,1);
        const double c = Mi(2,2);



        unsigned int dimNode = ((*BND_VERTEX_COLOR)[nodeIter] != 0)? SimplexMesh::topo::CORNER : ((*BND_CURVE_COLOR)[nodeIter] != 0)? SimplexMesh::topo::RIDGE : ((*BND_SURFACE_COLOR)[nodeIter] != 0)? SimplexMesh::topo::SURFACE : SimplexMesh::topo::VOLUME;
        unsigned int indexNode;

        if     (dimNode == SimplexMesh::topo::CORNER ){  indexNode =  (*BND_VERTEX_COLOR )[nodeIter] ; }
        else if(dimNode == SimplexMesh::topo::RIDGE  ){  indexNode =  (*BND_CURVE_COLOR  )[nodeIter] ; }
        else if(dimNode == SimplexMesh::topo::SURFACE){  indexNode =  (*BND_SURFACE_COLOR)[nodeIter] ; }
        else                                          {  indexNode =  0;}

        //sorting the node by nature
        std::vector<TInt> sortNeighboorNodes;
        if(nodeIter == 84){std::cout << "dimNode nodeIter -> " << dimNode << std::endl;}
        if(nodeIter == 84){std::cout << "indexNode nodeIter -> " << indexNode << std::endl;}
        for(auto const nNode : neighboorNodes)
        {
          if(dimNode == SimplexMesh::topo::SURFACE && ((*BND_VERTEX_COLOR)[nNode] != 0 || (*BND_CURVE_COLOR)[nNode] != 0 || (*BND_SURFACE_COLOR)[nNode] == indexNode))
          {
            sortNeighboorNodes.push_back(nNode);
          }
        }

        if(sortNeighboorNodes.size() == 0)
        {
          std::copy(neighboorNodes.begin(), neighboorNodes.end(), std::back_inserter(sortNeighboorNodes));
        }
        for(auto const nNode : sortNeighboorNodes)
        {
          const math::Point neighboorNodeCoord = SimplicesNode(m_simplexMesh, nNode).getCoords();
          const math::Vector3d vec = neighboorNodeCoord - currentNodeCoord;
          const double lenghtVec = 1.0;//vec.norm();
          //projection dans la base propre de la metric (ici i,j,k pour les test)
          const double x0_norm = (std::fabs(vec.X()) == 0.0) ? 1.0 : std::fabs(vec.X());
          const double y0_norm = (std::fabs(vec.Y()) == 0.0) ? 1.0 : std::fabs(vec.Y());
          const double z0_norm = (std::fabs(vec.Z()) == 0.0) ? 1.0 : std::fabs(vec.Z());

          const double alpha = (1.0 - a > 0.0)? 1.0 : -1.0;
          const double beta  = (1.0 - b > 0.0)? 1.0 : -1.0;
          const double gamma = (1.0 - c > 0.0)? 1.0 : -1.0;



          dE.x() = (std::fabs(1.0 - a) * vec.X() * alpha / x0_norm)  + dE.x();
          dE.y() = (std::fabs(1.0 - b) * vec.Y() * beta  / y0_norm)  + dE.y();
          dE.z() = (std::fabs(1.0 - c) * vec.Z() * gamma / z0_norm)  + dE.z();

          if(nodeIter == 84){
            std::cout << "nNode        -> " << nNode<< std::endl;
            std::cout << "vec        -> " << vec<< std::endl;
            std::cout << "dE.x()        -> " << dE.x()<< std::endl;
            std::cout << "dE.y()        -> " << dE.y()<< std::endl;
            std::cout << "dE.z()        -> " << dE.z()<< std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
          }
        }

        if(nodeIter == 84){std::cout << "dE.norm() -> " << dE.norm() << std::endl;}
        if(nodeIter == 84){std::cout << "dE        -> " << dE << std::endl;}
        if(nodeIter == 84){std::cout << "epsilon   -> " << epsilon << std::endl;}

        if(dE.norm() > epsilon)
        {
          math::Point vec = math::Point(dE.x(), dE.y(), dE.z());
          math::Point previousPosition = SimplicesNode(m_simplexMesh, nodeIter).getCoords();
          math::Point newPosition      = previousPosition - step * vec;
          if(nodeIter == 84){std::cout << "previousPosition -> " << previousPosition << std::endl;}
          if(nodeIter == 84){std::cout << "- step * vec -> " << - step * vec << std::endl;}
          if(nodeIter == 84){std::cout << "newPosition -> " << newPosition << std::endl;}

          m_simplexMesh->changeNodeCoord(nodeIter, newPosition);

          //changement de la mÃ©tric
          (*NODE_METRIC)[nodeIter](0,0) = (*NODE_METRIC)[nodeIter](0,0) + step * dE.x();
          (*NODE_METRIC)[nodeIter](1,1) = (*NODE_METRIC)[nodeIter](1,1) + step * dE.y();
          (*NODE_METRIC)[nodeIter](2,2) = (*NODE_METRIC)[nodeIter](2,2) + step * dE.z();
        }
      }
    }

    gmds::ISimplexMeshIOService ioService(m_simplexMesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    std::string name = "transformedCube_" + std::to_string(cpt) + ".vtk";
    vtkWriter.write(name);

    const math::Point neighboorNodeCoord = SimplicesNode(m_simplexMesh, 107).getCoords();
    const math::Point currentNodeCoord = SimplicesNode(m_simplexMesh, 112).getCoords();
    const math::Vector3d vecBis = neighboorNodeCoord - currentNodeCoord;

    const math::Point neighboorNodeCoord0 = SimplicesNode(m_simplexMesh, 81).getCoords();
    const math::Point currentNodeCoord0 = SimplicesNode(m_simplexMesh, 82).getCoords();
    const math::Vector3d vecBis0 = neighboorNodeCoord0 - currentNodeCoord0;

    std::cout << "VEC_FRONT_FACE.norm() --> " << vecBis.norm() << std::endl;
    std::cout << "VEC_BACK_FACE.norm() --> " << vecBis0.norm() << std::endl;
    std::cout << std::endl;

  }
}*/
