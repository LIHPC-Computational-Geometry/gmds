/******************************************************************************/
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace simplicesCell;
using namespace simplicesNode;
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
  m_simplex_mesh->addTetraedre(AID1, AID2, AID3, AID4, false);
}


void ISimplexMeshIOService::getNodes(std::vector<IMeshIOService::NodeInfo>& AInfo)
{
  const gmds::BitVector& nodeVector = m_simplex_mesh->getBitVectorNodes();
  AInfo.clear();
  AInfo.reserve(nodeVector.size());
  for(unsigned int node = 0 ; node < nodeVector.capacity() ; node++)
  {
    if(nodeVector[node] != 0)
    {
      IMeshIOService::NodeInfo info;
      info.id    = node;
      info.point = SimplicesNode(m_simplex_mesh, node).getCoords();
      AInfo.push_back(info);
    }
  }
}


void ISimplexMeshIOService::getEdges(std::vector<IMeshIOService::EdgeInfo>& AInfo){}


void ISimplexMeshIOService::getFaces(std::vector<IMeshIOService::CellInfo>& AInfo)
{
  const gmds::BitVector& triVector = m_simplex_mesh->getBitVectorTri();
  const int triVectorSize = triVector.size() - 1;
  AInfo.clear();
  AInfo.reserve(triVectorSize);

  /*we start at 1 and not at 0 because 0 is reserved for the tetraedron*/
  for(unsigned int tri = 1 ; tri < triVector.capacity() ; tri ++)
  {
    if(triVector[tri] != 0)
    {
      IMeshIOService::CellInfo info;
      info.id       = tri;
      info.type     = GMDS_TRIANGLE;
      info.node_ids = SimplicesTriangle(m_simplex_mesh,tri).nodes();
      AInfo.push_back(info);
    }
  }
}


void ISimplexMeshIOService::getRegions(std::vector<IMeshIOService::CellInfo>& AInfo)
{
  const gmds::BitVector& tetVector = m_simplex_mesh->getBitVectorTet();
  AInfo.clear();
  AInfo.reserve(tetVector.size());
  for(unsigned int tet = 0 ; tet < tetVector.capacity() ; tet ++)
  {
    if(tetVector[tet] != 0)
    {
      IMeshIOService::CellInfo info;
      info.id = tet;
      info.type = GMDS_TETRA;
      info.node_ids = SimplicesCell(m_simplex_mesh, tet).nodes();
      AInfo.push_back(info);
    }
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
  unsigned int cpt = 0;

  for(unsigned int nodeId = 0; nodeId < bitVecNode.capacity() ; nodeId++)
  {
    if(bitVecNode[nodeId] != 0)
    {
      ADataID.values[cpt] = nodeId;
      cpt++;
    }
  }

  std::vector<VariableItf*> vars = m_simplex_mesh->getAllVariables(GMDS_NODE);
  for(auto const &current_var : vars){
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
              break;
          }
          case(IMeshIOService::var_double) :
          {
            Variable<double>* v_vec = dynamic_cast<Variable<double>*> (current_var);
            IMeshIOService::DataReal data;
            data.name = v_vec->getName();
            data.values.resize(ANbCells);
            for(auto i_node : bitVecNode) {
                data.values[i_node] = (*v_vec)[i_node];
            }
            ADataReal.push_back(data);
            break;
          }
          case(IMeshIOService::var_int) :
          {
            Variable<int>* v_vec = dynamic_cast<Variable<int>*> (current_var);
            IMeshIOService::DataInt data;
            data.name = v_vec->getName();
            data.values.resize(ANbCells);
            for(auto i_node : bitVecNode) {
                data.values[i_node] = (*v_vec)[i_node];
            }
            ADataInt.push_back(data);
            break;
          }
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
                         std::vector<DataVector>& ADataVec)
{
  const gmds::BitVector& bitVecTri = m_simplex_mesh->getBitVectorTri();
  unsigned int ANbFaces = bitVecTri.size() - 1;

  ADataInt.clear();
  ADataInt.reserve(ANbFaces);

  ADataReal.clear();
  ADataReal.reserve(ANbFaces);

  ADataVec.clear();
  ADataVec.reserve(ANbFaces);

  ADataID.values.clear();
  ADataID.values.resize(ANbFaces);
  unsigned int cpt = 0;

  for(unsigned int triId = 1; triId < bitVecTri.capacity() ; triId++)
  {
    if(bitVecTri[triId] == 1)
    {
      ADataID.values[cpt] = triId;
      cpt++;
    }
  }
/*
  std::vector<VariableItf*> vars = m_simplex_mesh->getAllVariables(GMDS_TRIANGLE);
  for(auto current_var : vars){
      switch (IMeshIOService::getType(current_var)){
          case(IMeshIOService::var_double_vec) :
          {
              Variable<math::Vector3d>* v_vec = dynamic_cast<Variable<math::Vector3d>*> (current_var);
              IMeshIOService::DataVector data;
              data.name = v_vec->getName();
              data.values.resize(ANbFaces);
              for(auto i_tri : bitVecTri) {
                  data.values[i_tri] = (*v_vec)[i_tri];
              }
              ADataVec.push_back(data);
          }
              break;
          case(IMeshIOService::var_double) :
          {
            Variable<double>* v_vec = dynamic_cast<Variable<double>*> (current_var);
            IMeshIOService::DataReal data;
            data.name = v_vec->getName();
            data.values.resize(ANbFaces);
            for(auto i_tri : bitVecTri) {
                data.values[i_tri] = (*v_vec)[i_tri];
            }

            ADataReal.push_back(data);
          }
          case(IMeshIOService::var_int) :
          {
            std::cout << "var_int" << std::endl;
            Variable<int>* v_vec = dynamic_cast<Variable<int>*> (current_var);
            IMeshIOService::DataInt data;
            data.name = v_vec->name();
            data.values.resize(ANbFaces);
            for(unsigned int tri = 1 ; tri < bitVecTri.capacity() ; tri++)
            {
              if(bitVecTri[tri] == 1)
              {
                data.values[tri - 1] = (*v_vec)[tri];
              }
            }
            std::cout << "var_int end" << std::endl;
            ADataInt.push_back(data);
            std::cout << "push_back" << std::endl;
            break;
          }
          default:
              break;
      }
  }
*/
}

void ISimplexMeshIOService::getDataRegions(DataID& ADataID,
                           std::vector<DataInt>& ADataInt,
                           std::vector<DataReal>& ADataReal,
                           std::vector<DataVector>& ADataVec)
{
  const gmds::BitVector& bitVecTet = m_simplex_mesh->getBitVectorTet();
  unsigned int ANbCells = bitVecTet.size();
  ADataInt.clear();
  ADataInt.reserve(ANbCells);

  ADataReal.clear();
  ADataReal.reserve(ANbCells);

  ADataVec.clear();
  ADataVec.reserve(ANbCells);

  ADataID.values.clear();
  ADataID.values.resize(ANbCells);

  unsigned int cpt = 0;
  for(unsigned int tetId = 0; tetId < bitVecTet.capacity() ; tetId++)
  {
    if(bitVecTet[tetId] != 0)
    {
      ADataID.values[cpt] = tetId;
      cpt++;
    }
  }

  /*std::vector<VariableItf*> vars = m_simplex_mesh->getAllVariables(GMDS_NODE);
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
          case(IMeshIOService::var_double) :
          {
            Variable<double>* v_vec = dynamic_cast<Variable<double>*> (current_var);
            IMeshIOService::DataReal data;
            data.name = v_vec->name();
            data.values.resize(ANbCells);
            for(auto i_node : bitVecNode) {
                data.values[i_node] = (*v_vec)[i_node];
            }
            ADataReal.push_back(data);
          }
          default:
              break;
      }
  }
  */
}
/*----------------------------------------------------------------------------*/
//penser a creer un autre type pour ma structure de donnée plutot que d'utiliser ECellType et faire des etude de cas if pour savoir switch
//cela correspont a SimplicesNode,....
template<typename TBasicType, ECellType TCellType, typename TDataType>
struct addDataPolicy{
    void execute(gmds::hybrid::SimplexMesh* simplexMesh, TDataType& AData){
        Variable<TBasicType> *var;
        if(TCellType == GMDS_NODE)
        {
          var = simplexMesh->newVariable<TBasicType, SimplicesNode>(AData.name);
        }
        int nb_values = AData.values.size();
        var->setDomain(nb_values);
        for (auto i = 0; i < nb_values; i++) {
            (*var)[i] = AData.values[i];
        }
    }
};
/*----------------------------------------------------------------------------*/
void ISimplexMeshIOService::addDataIntNodes(DataInt& AData)
{
  addDataPolicy<int, GMDS_NODE,IMeshIOService::DataInt> policy;
  policy.execute(m_simplex_mesh,AData);
}
/*----------------------------------------------------------------------------*/
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
