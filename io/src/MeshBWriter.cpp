/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/io/MeshBWriter.h>
#include <fstream>
#include <libmeshb7.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
MeshBWriter::MeshBWriter(IMeshIOService *AMeshService)
        :IWriter(AMeshService)
{}
/*----------------------------------------------------------------------------*/
MeshBWriter::~MeshBWriter()
{}
/*----------------------------------------------------------------------------*/
void MeshBWriter::initialize(const std::string &AFileName) {

    m_index = GmfOpenMesh(AFileName.c_str(), GmfWrite, 1, 3 );
    if(!m_index){
        std::string s ="[LibMeshb error]Impossible to create a meshb.7 file  "+AFileName;
        throw GMDSException(s);
    }

}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeNodes() {
    std::vector<IMeshIOService::NodeInfo> nodes_info;
    m_mesh_service->getNodes(nodes_info);

    GmfSetKwd(m_index, GmfVertices, nodes_info.size() );

    auto meshb_node_id=0;
    for (auto info : nodes_info) {
        math::Point p = info.point;
        GmfSetLin( m_index, GmfVertices, p.X(), p.Y(), p.Z(), 1);
        m_node_ids_mapping[info.id] = meshb_node_id++;
    }
}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeEdges() {
    m_mesh_service->getEdges(m_edges_info);
}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeFaces()
{
    m_mesh_service->getFaces(m_faces_info);
    std::vector<int64_t> tri_idx;
    std::vector<int64_t> quad_idx;
    for(auto i=0; i<m_faces_info.size();i++){
        if (m_faces_info[i].type==GMDS_TRIANGLE)
            tri_idx.push_back(i);
        if (m_faces_info[i].type==GMDS_QUAD)
            quad_idx.push_back(i);
    }
    if(!tri_idx.empty()) {
        GmfSetKwd(m_index, GmfTriangles, tri_idx.size());
        for(auto current_idx: tri_idx) {
            IMeshIOService::CellInfo tri = m_faces_info[current_idx];
            GmfSetLin(m_index, GmfTriangles,
                      m_node_ids_mapping[tri.node_ids[0]]+1,
                      m_node_ids_mapping[tri.node_ids[1]]+1,
                      m_node_ids_mapping[tri.node_ids[2]]+1,
                      1);
        }
    }
    if(!quad_idx.empty()) {
        GmfSetKwd(m_index, GmfQuadrilaterals, quad_idx.size());
        for(auto current_idx: quad_idx) {
            IMeshIOService::CellInfo tri = m_faces_info[current_idx];
            GmfSetLin(m_index, GmfQuadrilaterals,
                      m_node_ids_mapping[tri.node_ids[0]]+1,
                      m_node_ids_mapping[tri.node_ids[1]]+1,
                      m_node_ids_mapping[tri.node_ids[2]]+1,
                      m_node_ids_mapping[tri.node_ids[3]]+1,
                      1);
        }
    }
}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeRegions()
{
    m_mesh_service->getRegions(m_regions_info);
    std::vector<int64_t> tet_idx;
    for(auto i=0; i<m_regions_info.size();i++){
        if (m_regions_info[i].type==GMDS_TETRA)
            tet_idx.push_back(i);
    }
    if(!tet_idx.empty()) {
        GmfSetKwd(m_index, GmfTetrahedra, tet_idx.size());
        for(auto current_idx: tet_idx) {
            IMeshIOService::CellInfo tet = m_regions_info[current_idx];
            GmfSetLin(m_index, GmfTetrahedra,
                      m_node_ids_mapping[tet.node_ids[0]]+1,
                      m_node_ids_mapping[tet.node_ids[1]]+1,
                      m_node_ids_mapping[tet.node_ids[2]]+1,
                      m_node_ids_mapping[tet.node_ids[3]]+1,
                      1);
        }
    }
}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeDataNodes() {

    IMeshIOService::DataID data_id;
    std::vector<IMeshIOService::DataInt> data_int;
    std::vector<IMeshIOService::DataReal> data_real;
    std::vector<IMeshIOService::DataVector> data_vec;
    m_mesh_service->getDataNodes(data_id,data_int, data_real, data_vec);
}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeDataEdges()  {

    m_mesh_service->getDataEdges(m_cell_id[0],
                                 m_cell_var_int[0],
                                 m_cell_var_real[0],
                                 m_cell_var_vec[0]);

}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeDataFaces()
{
    m_mesh_service->getDataFaces(m_cell_id[1],
                                 m_cell_var_int[1],
                                 m_cell_var_real[1],
                                 m_cell_var_vec[1]);
}
/*----------------------------------------------------------------------------*/
void MeshBWriter::writeDataRegions()
{
    m_mesh_service->getDataRegions(m_cell_id[2],
                                   m_cell_var_int[2],
                                   m_cell_var_real[2],
                                   m_cell_var_vec[2]);
}
/*----------------------------------------------------------------------------*/
void MeshBWriter::finalize() {
    GmfCloseMesh(m_index);
}