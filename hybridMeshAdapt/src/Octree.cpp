#include "gmds/hybridMeshAdapt/Octree.h"
////////////////////////////////////////////////////////////////////////////////
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace math;
////////////////////////////////////////////////////////////////////////////////
Octree::Octree(SimplexMesh* simplexMesh,
       const unsigned int numbersMaxSimplices) :
       m_simplexMesh(simplexMesh),
       m_numbersMaxSimplices(numbersMaxSimplices),
       m_rootOc(nullptr),
       m_parentOc(nullptr)
{
  std::vector<std::vector<Node>> nodes{};
  initialize();

  Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));



  writeOctree(m, nodes);
  for(auto const & node : nodes)
  {
    m.newHex(node[0], node[1], node[2], node[3], node[4], node[5], node[6], node[7]);
  }
  gmds::IGMeshIOService ioService(&m);
  gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("Octree.vtk");
}
/******************************************************************************/
Octree::Octree(SimplexMesh* simplexMesh,
      unsigned int subdivision,
      double xmin, double xmax,
      double ymin, double ymax,
      double zmin, double zmax):
    m_xmin(xmin), m_xmax(xmax),
    m_ymin(ymin), m_ymax(ymax),
    m_zmin(zmin), m_zmax(zmax),
    m_simplexMesh(simplexMesh)
{
  --subdivision;
  if(subdivision < 1)
    return;

  Octree * oc0 = new Octree(simplexMesh, subdivision, m_xmin, (m_xmax + m_xmin) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) /2.0);
  Octree * oc1 = new Octree(simplexMesh, subdivision, (m_xmax + m_xmin) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) /2.0);
  Octree * oc2 = new Octree(simplexMesh, subdivision, m_xmin, (m_xmax + m_xmin) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) /2.0, m_zmax);
  Octree * oc3 = new Octree(simplexMesh, subdivision, (m_xmax + m_xmin) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) /2.0, m_zmax);

  Octree * oc4 = new Octree(simplexMesh, subdivision, m_xmin, (m_xmax + m_xmin) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) /2.0);
  Octree * oc5 = new Octree(simplexMesh, subdivision, (m_xmax + m_xmin) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) /2.0);
  Octree * oc6 = new Octree(simplexMesh, subdivision, m_xmin, (m_xmax + m_xmin) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) /2.0, m_zmax);
  Octree * oc7 = new Octree(simplexMesh, subdivision, (m_xmax + m_xmin) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) /2.0, m_zmax);

  oc0->setParentOctree(this); oc1->setParentOctree(this); oc2->setParentOctree(this); oc3->setParentOctree(this);
  oc4->setParentOctree(this); oc5->setParentOctree(this); oc6->setParentOctree(this); oc7->setParentOctree(this);

  m_ocs[0] = oc0 ; m_ocs[1] = oc1 ; m_ocs[2] = oc2 ; m_ocs[3] = oc3 ; m_ocs[4] = oc4 ; m_ocs[5] = oc5 ; m_ocs[6] = oc6 ; m_ocs[7] = oc7;
}
/******************************************************************************/
void Octree::addNode(TInt node, std::unordered_set<TInt>& seen)
{
  const math::Point pt = SimplicesNode(m_simplexMesh, node).getCoords();

  for(Octree* oc : m_ocs)
  {
    if(oc != nullptr)
    {
      if(oc->belongToOc(pt))
      {
        oc->addNode(node, seen);
      }

    }
    else
    {
      if(seen.find(node) == seen.end())
      {
        seen.insert(node);
        m_nodes.push_back(node);
      }
      return;
    }
  }
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
void Octree::writeOctree(Mesh& mesh, std::vector<std::vector<Node>>& nodes) const
{
  const Point p0(m_xmin, m_ymin, m_zmin); const Point p1(m_xmin, m_ymax, m_zmin);
  const Point p2(m_xmin, m_ymax, m_zmax); const Point p3(m_xmin, m_ymin, m_zmax);

  const Point p4(m_xmax, m_ymin, m_zmin); const Point p5(m_xmax, m_ymax, m_zmin);
  const Point p6(m_xmax, m_ymax, m_zmax); const Point p7(m_xmax, m_ymin, m_zmax);

  Node n0 = mesh.newNode(p0); Node n1 = mesh.newNode(p1); Node n2 = mesh.newNode(p2); Node n3 = mesh.newNode(p3);
  Node n4 = mesh.newNode(p4); Node n5 = mesh.newNode(p5); Node n6 = mesh.newNode(p6); Node n7 = mesh.newNode(p7);

  nodes.push_back(std::vector<Node>{n0, n1, n2, n3, n4, n5, n6, n7});
  for(auto const & oc : m_ocs)
  {
    if(oc != nullptr)
    {
      oc->writeOctree(mesh, nodes);
    }
  }
}
/******************************************************************************/
void Octree::initialize()
{
  const BitVector& nodesVector = m_simplexMesh->getBitVectorNodes();
  double epsilon = 10E-4;
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();
  std::set<TSimplexID> s{};
  std::vector<double> nodes{max, min, max, min, max, min}; // --> xmin, xmax, ymin, ymax, zmin, zmax

  for(unsigned int nodeIds = 0; nodeIds < nodesVector.capacity(); nodeIds++)
  {
    if(nodesVector[nodeIds] != 0)
    {
      m_nodes.push_back(nodeIds);
      math::Point pt = SimplicesNode(m_simplexMesh, nodeIds).getCoords();
      std::vector<TSimplexID> ball = SimplicesNode(m_simplexMesh, nodeIds).ballOf();
      for(auto const simplex : ball)
      {
        if(s.find(simplex) == s.end())
        {
          s.insert(simplex);
          m_simplices.push_back(simplex);
        }
      }

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

    oc0->setRootOctree(this); oc1->setRootOctree(this); oc2->setRootOctree(this); oc3->setRootOctree(this);
    oc4->setRootOctree(this); oc5->setRootOctree(this); oc6->setRootOctree(this); oc7->setRootOctree(this);

    oc0->setParentOctree(this); oc1->setParentOctree(this); oc2->setParentOctree(this); oc3->setParentOctree(this);
    oc4->setParentOctree(this); oc5->setParentOctree(this); oc6->setParentOctree(this); oc7->setParentOctree(this);
    //root and parent is the same for direct child of root

    m_ocs[0] = oc0 ; m_ocs[1] = oc1 ; m_ocs[2] = oc2 ; m_ocs[3] = oc3 ; m_ocs[4] = oc4 ; m_ocs[5] = oc5 ; m_ocs[6] = oc6 ; m_ocs[7] = oc7;
  }
}
/******************************************************************************/
void Octree::preprocess()
{
  const BitVector& nodesVector = m_simplexMesh->getBitVectorNodes();
  std::vector<TInt> nodesTmp{};
  std::set<TSimplexID> s{};

  for(unsigned int nodeIds = 0; nodeIds < nodesVector.capacity(); nodeIds++)
  {
    if(nodesVector[nodeIds] != 0)
    {
      math::Point pt = SimplicesNode(m_simplexMesh, nodeIds).getCoords();
      if(belongToOc(pt))
      {
        nodesTmp.push_back(nodeIds);
        std::vector<TSimplexID> ball = SimplicesNode(m_simplexMesh, nodeIds).ballOf();
        for(auto const simplex : ball)
        {
          if(s.find(simplex) == s.end())
          {
            s.insert(simplex);
            m_simplices.push_back(simplex);
          }
        }
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

    oc0->setRootOctree(this->m_rootOc); oc1->setRootOctree(this->m_rootOc); oc2->setRootOctree(this->m_rootOc); oc3->setRootOctree(this->m_rootOc);
    oc4->setRootOctree(this->m_rootOc); oc5->setRootOctree(this->m_rootOc); oc6->setRootOctree(this->m_rootOc); oc7->setRootOctree(this->m_rootOc);

    oc0->setParentOctree(this); oc1->setParentOctree(this); oc2->setParentOctree(this); oc3->setParentOctree(this);
    oc4->setParentOctree(this); oc5->setParentOctree(this); oc6->setParentOctree(this); oc7->setParentOctree(this);

    m_ocs[0] = oc0 ; m_ocs[1] = oc1 ; m_ocs[2] = oc2 ; m_ocs[3] = oc3 ; m_ocs[4] = oc4 ; m_ocs[5] = oc5 ; m_ocs[6] = oc6 ; m_ocs[7] = oc7;
  }
  else
  {
    return;
  }
}
/******************************************************************************/
bool Octree::belongToOc(const math::Point& pt) const
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
          return m_simplexMesh->getSimplexFromBase(node);
        }
      }
    }
  }
  return simplex;
}
/******************************************************************************/
std::vector<TSimplexID> Octree::findSimplicesInOc(const math::Point& pt)
{
  std::vector<TSimplexID> ans {};
  std::unordered_set<TSimplexID> us{};
  for(auto const & oc : m_ocs)
  {
    if(oc != nullptr)
    {
      if(oc->belongToOc(pt))
      {
        return oc->findSimplicesInOc(pt);
      }
    }
    else
    {
      const BitVector& nodesVector = m_simplexMesh->getBitVectorNodes();
      ans = m_parentOc->m_simplices;
      /*for(auto const node : m_parentOc->m_nodes)
      {
        std::vector<TSimplexID> ball = SimplicesNode(m_simplexMesh, node).ballOf();
        for(auto const simplex : ball)
        {
          if(us.find(simplex) == us.end())
          {
            us.insert(simplex);
            ans.push_back(simplex);
          }
        }
      }*/
    }
  }
  return ans;
}
/******************************************************************************/
std::vector<TInt> Octree::findNodesNextTo(const math::Point& pt)
{
  for(auto const & oc : m_ocs)
  {
    if(oc != nullptr)
    {
      if(oc->belongToOc(pt))
      {
        return oc->findNodesNextTo(pt);
      }
    }
  }

  std::vector<TInt> ans{};
  for(auto const oc : m_parentOc->m_parentOc->m_ocs)
  {
    for(auto const o : oc->m_ocs)
    {
      double xmin = o->getBorderOctree()[0] ; double xmax = o->getBorderOctree()[1] ;
      double ymin = o->getBorderOctree()[2] ; double ymax = o->getBorderOctree()[3] ;
      double zmin = o->getBorderOctree()[4] ; double zmax = o->getBorderOctree()[5] ;
      if((xmax == m_xmin || xmin == m_xmin || xmin == m_xmax) &&
      (ymax == m_ymin || ymin == m_ymin || ymin == m_ymax) &&
      (zmax == m_zmin || zmin == m_zmin || zmin == m_zmax))
      {
        for(auto const n : o->m_nodes)
        {
          ans.push_back(n);
        }
      }
    }
  }

  return ans;

}
/******************************************************************************/
std::vector<double> Octree::getBorderOctree() const
{
  std::vector<double> ans{m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax};
  return ans;
}
