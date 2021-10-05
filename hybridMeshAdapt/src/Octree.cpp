#include "gmds/hybridMeshAdapt/Octree.h"

////////////////////////////////////////////////////////////////////////////////
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;

Octree::Octree(SimplexMesh* simplexMesh,
       const unsigned int numbersMaxSimplices) :
       m_simplexMesh(simplexMesh),
       m_numbersMaxSimplices(numbersMaxSimplices)
{
  initialize();
}
/******************************************************************************/
Octree::Octree(
      double xmin, double xmax,
      double ymin, double ymax,
      double zmin, double zmax,
      SimplexMesh* simplexMesh,
      const unsigned int numbersMaxSimplices) :
      m_xmin(xmin), m_xmax(xmax),
      m_ymin(ymin), m_ymax(ymax),
      m_zmin(zmin), m_zmax(zmax),
      m_simplexMesh(simplexMesh),
      m_numbersMaxSimplices(numbersMaxSimplices)
{
  preprocess();
}
/******************************************************************************/
Octree::~Octree()
{
  for(unsigned int i = 0 ; i < m_ocs.size() ; i++)
  {
    Octree* oc = m_ocs[i];
    delete(oc);
  }
}
/******************************************************************************/
void Octree::initialize()
{
  const BitVector& nodesVector = m_simplexMesh->getBitVectorNodes();
  double epsilon = 10E-2;
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  std::vector<double> nodes{max, min, max, min, max, min}; // --> xmin, xmax, ymin, ymax, zmin, zmax
  m_nodes.reserve(nodesVector.size());

  for(unsigned int nodeIds = 0; nodeIds < nodesVector.capacity(); nodeIds++)
  {
    if(nodesVector[nodeIds] != 0)
    {
      m_nodes.push_back(nodeIds);
      math::Point pt = SimplicesNode(m_simplexMesh, nodeIds).getCoords();
      double x = pt.X();
      double y = pt.Y();
      double z = pt.Z();

      nodes[0] = (x < nodes[0])? x : nodes[0];
      nodes[1] = (x > nodes[1])? x : nodes[1];

      nodes[2] = (y < nodes[2])? y : nodes[2];
      nodes[3] = (y > nodes[3])? y : nodes[3];

      nodes[4] = (z < nodes[4])? z : nodes[4];
      nodes[5] = (z > nodes[5])? z : nodes[5];
    }
  }

  m_xmin = nodes[0] - epsilon; m_xmax = nodes[1] + epsilon;
  m_ymin = nodes[2] - epsilon; m_ymax = nodes[3] + epsilon;
  m_zmin = nodes[4] - epsilon; m_zmax = nodes[5] + epsilon;

  if(m_nodes.size() > m_numbersMaxSimplices)
  {

    Octree * oc0 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc1 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc2 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc3 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);

    Octree * oc4 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc5 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc6 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc7 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);

    m_ocs[0] = oc0 ; m_ocs[1] = oc1 ; m_ocs[2] = oc2 ; m_ocs[3] = oc3 ; m_ocs[4] = oc4 ; m_ocs[5] = oc5 ; m_ocs[6] = oc6 ; m_ocs[7] = oc7;
  }
}
/******************************************************************************/
void Octree::preprocess()
{
  const BitVector& nodesVector = m_simplexMesh->getBitVectorNodes();
  std::vector<TInt> nodesTmp{};

  for(unsigned int nodeIds = 0; nodeIds < nodesVector.capacity(); nodeIds++)
  {
    if(nodesVector[nodeIds] != 0)
    {
      math::Point pt = SimplicesNode(m_simplexMesh, nodeIds).getCoords();
      if(belongToOc(pt))
      {
        nodesTmp.push_back(nodeIds);
      }
    }
  }
  std::copy(nodesTmp.begin(), nodesTmp.end(), std::back_inserter(m_nodes));
  if(m_nodes.size() > m_numbersMaxSimplices)
  {
    Octree * oc0 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc1 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc2 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc3 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);

    Octree * oc4 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc5 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) /2.0, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc6 = new Octree(m_xmin, (m_xmax + m_xmin) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);
    Octree * oc7 = new Octree((m_xmax + m_xmin) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) /2.0, m_zmax, m_simplexMesh, m_numbersMaxSimplices);

    m_ocs[0] = oc0 ; m_ocs[1] = oc1 ; m_ocs[2] = oc2 ; m_ocs[3] = oc3 ; m_ocs[4] = oc4 ; m_ocs[5] = oc5 ; m_ocs[6] = oc6 ; m_ocs[7] = oc7;
  }
  else
  {
    return;
  }
}
/******************************************************************************/
bool Octree::belongToOc(const math::Point& pt)
{
  double x = pt.X();
  double y = pt.Y();
  double z = pt.Z();

  if(x >= m_xmin && x <= m_xmax)
  {
    if(y >= m_ymin && y <= m_ymax)
    {
      if(z >= m_zmin && z <= m_zmax)
      {
          return true;
      }
    }
  }
  return false;
}
/******************************************************************************/
TSimplexID Octree::findSimplexNextTo(const math::Point& pt)
{
  TSimplexID simplex = std::numeric_limits<int>::min();
  for(auto const & oc : m_ocs)
  {
    if(oc != nullptr)
    {
      if(oc->belongToOc(pt))
      {
        simplex = oc->findSimplexNextTo(pt);
      }
    }
    else
    {
      const BitVector& nodesVector = m_simplexMesh->getBitVectorNodes();
      for(auto const node : m_nodes)
      {
        if(nodesVector[node] != 0)
        {
          return m_simplexMesh->getSimplexFromBase(node);;
        }
      }
    }
  }

  return simplex;
}
/******************************************************************************/
std::vector<TInt> Octree::findNodesNextTo(const math::Point& pt)
{
  std::vector<TSimplexID> nodes{};
  for(auto const & oc : m_ocs)
  {
    if(oc != nullptr)
    {
      if(oc->belongToOc(pt))
      {
        nodes = oc->findNodesNextTo(pt);
      }
    }
    /*else
    {
      return m_nodes;
    }*/
  }

  return m_nodes;
}
/******************************************************************************/
