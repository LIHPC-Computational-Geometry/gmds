/******************************************************************************/
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace simplicesCell;
using namespace simplicesTriangle;

ISimplexMeshIOService::ISimplexMeshIOService(SimplexMesh* simplexMesh):
m_simplex_mesh(simplexMesh)
{
}


void ISimplexMeshIOService::createNodes(const std::vector<double>& AX,
                        const std::vector<double>& AY,
                        const std::vector<double>& AZ)
{
  for(size_t i = 0; i < AX.size(); i++)
  {
    m_simplex_mesh->addNode(AX[i],AY[i],AZ[i]);
  }
}


void ISimplexMeshIOService::createTriangle(const TCellID& AID1,
                           const TCellID& AID2,
                           const TCellID& AID3)
{
  m_simplex_mesh->addTriangle(AID1, AID2, AID3);
}


void ISimplexMeshIOService::createTet(const TCellID& AID1,
                      const TCellID& AID2,
                      const TCellID& AID3,
                      const TCellID& AID4)
{
  m_simplex_mesh->addTetraedre(AID1, AID2, AID3, AID4);
}


void ISimplexMeshIOService::getNodes(std::vector<IMeshIOService::NodeInfo>& AInfo)
{
  AInfo.clear();
  AInfo.reserve(m_simplex_mesh->getNbNodes());

  for(TInt nodeIndx = m_simplex_mesh->getFirstNode() ; nodeIndx <= m_simplex_mesh->nodeCapacity() ; nodeIndx = m_simplex_mesh->nextNode())
  {
    IMeshIOService::NodeInfo info;
    info.id    = nodeIndx;
    info.point = SimplicesNode(m_simplex_mesh,nodeIndx).getCoords();
    AInfo.push_back(info);
  }
}


void ISimplexMeshIOService::getEdges(std::vector<IMeshIOService::EdgeInfo>& AInfo){}


void ISimplexMeshIOService::getFaces(std::vector<IMeshIOService::CellInfo>& AInfo)
{
  AInfo.clear();
  AInfo.reserve(m_simplex_mesh->getNbTriangle());

  for(TInt triIndx = m_simplex_mesh->getFirstTri() ; triIndx <= m_simplex_mesh->triCapacity() ; triIndx = m_simplex_mesh->nextTri())
  {
    IMeshIOService::CellInfo info;;
    info.id       = triIndx;
    info.type     = GMDS_TRIANGLE;
    info.node_ids = SimplicesTriangle(m_simplex_mesh,triIndx).nodes();;
    AInfo.push_back(info);
  }
}


void ISimplexMeshIOService::getRegions(std::vector<IMeshIOService::CellInfo>& AInfo)
{
  AInfo.clear();
  AInfo.reserve(m_simplex_mesh->getNbTetra());

  for (TInt tetIndx = m_simplex_mesh->getFirstTet() ; tetIndx <= m_simplex_mesh->tetCapacity() ; tetIndx = m_simplex_mesh->nextTet())
  {

      IMeshIOService::CellInfo info;
      info.id = tetIndx;
      info.type = GMDS_TETRA;
      info.node_ids = SimplicesCell(m_simplex_mesh, tetIndx).nodes();
      AInfo.push_back(info);
  }
}


void ISimplexMeshIOService::getDataNodes(DataID& ADataID,
             std::vector<DataInt> &ADataInt,
             std::vector<DataReal> &ADataReal,
             std::vector<DataVector>& ADataVec)
{
  const gmds::BitVector& bitVecNode = m_simplex_mesh->getBitVectorNodes();
  unsigned int ANbCells = bitVecNode.size();
  ADataInt.clear();
  ADataInt.reserve(ANbCells);

  ADataReal.clear();
  ADataReal.reserve(ANbCells);

  ADataVec.clear();
  ADataVec.reserve(ANbCells);

  ADataID.values.clear();
  ADataID.values.resize(ANbCells);


  for(auto i_node : bitVecNode) {
      ADataID.values[i_node] = i_node;
  }


  std::vector<VariableItf*> vars = m_simplex_mesh->getAllVariables(GMDS_NODE);
  for(auto current_var : vars){
      switch (IMeshIOService::getType(current_var)){
          case(IMeshIOService::var_double_vec) :
          {
              Variable<math::Vector3d>* v_vec = dynamic_cast<Variable<math::Vector3d>*> (current_var);
              IMeshIOService::DataVector data;
              data.name = v_vec->getName();
              data.values.resize(ANbCells);
              for(auto i_node : bitVecNode) {
                  data.values[i_node] = (*v_vec)[i_node];
              }
              ADataVec.push_back(data);
          }
              break;
          default:
              break;
      }
  }
}

void ISimplexMeshIOService::getDataEdges(DataID& ADataID,
                         std::vector<DataInt>& ADataInt,
                         std::vector<DataReal>& ADataReal,
                         std::vector<DataVector>& ADataVec){}
void ISimplexMeshIOService::getDataFaces(DataID& ADataID,
                         std::vector<DataInt>& ADataInt,
                         std::vector<DataReal>& ADataReal,
                         std::vector<DataVector>& ADataVec){}
void ISimplexMeshIOService::getDataRegions(DataID& ADataID,
                           std::vector<DataInt>& ADataInt,
                           std::vector<DataReal>& ADataReal,
                           std::vector<DataVector>& ADataVec){}

void ISimplexMeshIOService::addDataIntNodes(DataInt& AData) {}
void ISimplexMeshIOService::addDataRealNodes(DataReal& AData) {}
void ISimplexMeshIOService::addDataVectorNodes(DataVector& AData) {}

void ISimplexMeshIOService::addDataIntEdges(DataInt& AData) {}
void ISimplexMeshIOService::addDataRealEdges(DataReal& AData) {}
void ISimplexMeshIOService::addDataVectorEdges(DataVector& AData) {}

void ISimplexMeshIOService::addDataIntFaces(DataInt& AData) {}
void ISimplexMeshIOService::addDataRealFaces(DataReal& AData) {}
void ISimplexMeshIOService::addDataVectorFaces(DataVector& AData) {}

void ISimplexMeshIOService::addDataIntRegions(DataInt& AData) {}
void ISimplexMeshIOService::addDataRealRegions(DataReal& AData) {}
void ISimplexMeshIOService::addDataVectorRegions(DataVector& AData) {}
