/*----------------------------------------------------------------------------*/
//
// Created by totoro on 2019-01-12.
//
/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
#include <gmds/io/IGMeshIOService.h>
#include <gmds/ig/Mesh.h>
#include <gmds/utils/Variable.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
IGMeshIOService::IGMeshIOService(Mesh *AMesh)
        :m_mesh(AMesh){}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::createNodes(
        const std::vector<double> &AX,
        const std::vector<double> &AY,
        const std::vector<double> &AZ)
{
    unsigned int nb_nodes = AX.size();
    for(unsigned int i = 0; i<nb_nodes; i++) {
        m_mesh->newNode(AX[i], AY[i], AZ[i]);
    }
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::createEdge(const TCellID& AID1,
                                 const TCellID& AID2) {
    m_mesh->newEdge(AID1, AID2);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::createTriangle(const TCellID& AID1,
                                     const TCellID& AID2,
                                     const TCellID& AID3) {
    m_mesh->newTriangle(AID1, AID2, AID3);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::createQuad(const TCellID& AID1,
                const TCellID& AID2,
                const TCellID& AID3,
                const TCellID& AID4){
    m_mesh->newQuad(AID1, AID2, AID3, AID4);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::createPolygon(const std::vector<gmds::TCellID> &ANodes)
{
    m_mesh->newPolygon(ANodes);
}

/*----------------------------------------------------------------------------*/
void IGMeshIOService:: createTet(const TCellID& AID1,
                                 const TCellID& AID2,
                                 const TCellID& AID3,
                                 const TCellID& AID4)
{
    m_mesh->newTet(AID1,AID2,AID3,AID4);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::createHex(const TCellID& AID1,
                                const TCellID& AID2,
                                const TCellID& AID3,
                                const TCellID& AID4,
                                const TCellID& AID5,
                                const TCellID& AID6,
                                const TCellID& AID7,
                                const TCellID& AID8)
{
    m_mesh->newHex(AID1,AID2,AID3,AID4,AID5,AID6,AID7,AID8);
}

/*----------------------------------------------------------------------------*/
void IGMeshIOService::createPyramid(const TCellID& AID1,
                                    const TCellID& AID2,
                                    const TCellID& AID3,
                                    const TCellID& AID4,
                                    const TCellID& AID5)
{
    m_mesh->newPyramid(AID1,AID2,AID3,AID4,AID5);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getNodes(std::vector<IMeshIOService::NodeInfo> &AInfo) {
    AInfo.clear();

    AInfo.reserve(m_mesh->getNbNodes());

    for (auto node_id:m_mesh->nodes()) {
        Node n = m_mesh->get<Node>(node_id);
        IMeshIOService::NodeInfo info;
        info.id = node_id;
        info.point = n.point();
        AInfo.push_back(info);
    }
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getEdges(std::vector<gmds::IMeshIOService::EdgeInfo> &AInfo) {
    AInfo.clear();
    AInfo.reserve(m_mesh->getNbEdges());

    for (auto edge_id:m_mesh->edges()) {
        Edge e = m_mesh->get<Edge>(edge_id);
        IMeshIOService::EdgeInfo info;
        info.id = edge_id;
        std::vector<TCellID> n_ids = e.getIDs<Node>();
        info.node_ids[0]=n_ids[0];
        info.node_ids[1]=n_ids[1];
        AInfo.push_back(info);
    }
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getFaces(std::vector<gmds::IMeshIOService::CellInfo> &AInfo) {
    AInfo.clear();
    AInfo.reserve(m_mesh->getNbFaces());

    for (auto face_id:m_mesh->faces()) {
        Face f = m_mesh->get<Face>(face_id);
        IMeshIOService::CellInfo info;
        info.id = face_id;
        info.type = f.type();
        info.node_ids = f.getIDs<Node>();
        AInfo.push_back(info);
    }
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getRegions(std::vector<gmds::IMeshIOService::CellInfo> &AInfo) {
    AInfo.clear();
    AInfo.reserve(m_mesh->getNbRegions());

    for (auto r_id:m_mesh->regions()) {
        Region r = m_mesh->get<Region>(r_id);
        IMeshIOService::CellInfo info;
        info.id = r_id;
        info.type = r.type();
        info.node_ids = r.getIDs<Node>();
        AInfo.push_back(info);
    }
}
/*----------------------------------------------------------------------------*/
template <typename TCellContainer, ECellType TCellType> struct GetDataPolicy{
    void execute(Mesh* AMesh,
                 int ANbCells,
                 TCellContainer& ACells,
                 IGMeshIOService::DataID& ADataID,
                 std::vector<IMeshIOService::DataInt> &ADataInt,
                 std::vector<IMeshIOService::DataReal> &ADataReal,
                 std::vector<IGMeshIOService::DataVector>& ADataVec)
    {
        ADataInt.clear();
        ADataInt.reserve(ANbCells);

        ADataReal.clear();
        ADataReal.reserve(ANbCells);

        ADataVec.clear();
        ADataVec.reserve(ANbCells);

        ADataID.values.clear();

        for( auto i_node : ACells) {
            ADataID.values[i_node] = i_node;
        }

        std::vector<VariableItf*> vars = AMesh->getAllVariables(TCellType);
        for(auto current_var : vars){
            switch (IMeshIOService::getType(current_var)){
                case(IMeshIOService::var_int) :
                {
                    auto* v_int = dynamic_cast<Variable<int>*> (current_var);
                    IMeshIOService::DataInt data;
                    data.name = v_int->getName();
                    for(auto i_node : ACells) {
                        data.values[i_node] = (*v_int)[i_node];
                    }
                    ADataInt.push_back(data);

                }
                    break;
                case(IMeshIOService::var_double) :
                {
                    auto* v_double = dynamic_cast<Variable<double>*> (current_var);
                    IMeshIOService::DataReal data;
                    data.name = v_double->getName();
                    for(auto i_node : ACells) {
                        data.values[i_node] = (*v_double)[i_node];
                    }
                    ADataReal.push_back(data);
                }
                    break;
                case(IMeshIOService::var_double_vec) :
                {
                    auto* v_vec = dynamic_cast<Variable<math::Vector3d>*> (current_var);
                    IMeshIOService::DataVector data;
                    data.name = v_vec->getName();
                    for(auto i_node : ACells) {
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
};
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getDataNodes(DataID& ADataID,
             std::vector<DataInt> &ADataInt,
             std::vector<DataReal> &ADataReal,
             std::vector<DataVector>& ADataVec)
{
    GetDataPolicy<NodeContainer, GMDS_NODE> policy;
    policy.execute(m_mesh, m_mesh->getNbNodes(), m_mesh->nodes(),
                   ADataID, ADataInt, ADataReal, ADataVec);

}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getDataEdges(DataID& ADataID,
             std::vector<gmds::IMeshIOService::DataInt> &ADataInt,
             std::vector<gmds::IMeshIOService::DataReal> &ADataReal,
             std::vector<DataVector>& ADataVec)
{
    GetDataPolicy<EdgeContainer, GMDS_EDGE> policy;
    policy.execute(m_mesh, m_mesh->getNbEdges(), m_mesh->edges(),
                   ADataID, ADataInt, ADataReal, ADataVec);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getDataFaces(DataID& ADataID,
             std::vector<gmds::IMeshIOService::DataInt> &ADataInt,
             std::vector<gmds::IMeshIOService::DataReal> &ADataReal,
             std::vector<DataVector>& ADataVec)
{
    GetDataPolicy<FaceContainer, GMDS_FACE> policy;
    policy.execute(m_mesh, m_mesh->getNbFaces(), m_mesh->faces(),
                   ADataID, ADataInt, ADataReal, ADataVec);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::
getDataRegions(gmds::IMeshIOService::DataID &ADataID,
               std::vector<gmds::IMeshIOService::DataInt> &ADataInt,
               std::vector<gmds::IMeshIOService::DataReal> &ADataReal,
               std::vector<gmds::IMeshIOService::DataVector> &ADataVec)
{
    GetDataPolicy<RegionContainer, GMDS_REGION> policy;
    policy.execute(m_mesh, m_mesh->getNbRegions(), m_mesh->regions(),
                   ADataID, ADataInt, ADataReal, ADataVec);
}

/*----------------------------------------------------------------------------*/
template<typename TBasicType, ECellType TCellType, typename TDataType>
struct addDataPolicy{
    void execute(Mesh* AMesh, TDataType& AData){
        Variable<TBasicType> *var = AMesh->newVariable<TBasicType, TCellType>(AData.name);
        int nb_values = AData.values.size();
        var->setDomain(nb_values);
        for (auto i = 0; i < nb_values; i++) {
            (*var)[i] = AData.values[i];
        }
    }
};
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataIntNodes(gmds::IMeshIOService::DataInt &AData) {
   addDataPolicy<int, GMDS_NODE,IMeshIOService::DataInt> policy;
   policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataRealNodes(gmds::IMeshIOService::DataReal &AData)
{
    addDataPolicy<double, GMDS_NODE,IMeshIOService::DataReal> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataVectorNodes(DataVector& AData){
    addDataPolicy<math::Vector3d, GMDS_NODE,IMeshIOService::DataVector> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataIntEdges(gmds::IMeshIOService::DataInt &AData) {
    addDataPolicy<int, GMDS_EDGE,IMeshIOService::DataInt> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataRealEdges(gmds::IMeshIOService::DataReal &AData)
{
    addDataPolicy<double, GMDS_EDGE, IMeshIOService::DataReal> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataVectorEdges(DataVector& AData){
    addDataPolicy<math::Vector3d, GMDS_EDGE,IMeshIOService::DataVector> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataIntFaces(gmds::IMeshIOService::DataInt &AData) {
    addDataPolicy<int, GMDS_FACE,IMeshIOService::DataInt> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataRealFaces(gmds::IMeshIOService::DataReal &AData)
{
    addDataPolicy<double, GMDS_FACE, IMeshIOService::DataReal> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataVectorFaces(DataVector& AData){
    addDataPolicy<math::Vector3d, GMDS_FACE, IMeshIOService::DataVector> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataIntRegions(gmds::IMeshIOService::DataInt &AData) {
    addDataPolicy<int, GMDS_REGION,IMeshIOService::DataInt> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataRealRegions(gmds::IMeshIOService::DataReal &AData)
{
    addDataPolicy<double, GMDS_REGION, IMeshIOService::DataReal> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
void IGMeshIOService::addDataVectorRegions(DataVector& AData){
    addDataPolicy<math::Vector3d, GMDS_REGION, IMeshIOService::DataVector> policy;
    policy.execute(m_mesh,AData);
}
/*----------------------------------------------------------------------------*/
