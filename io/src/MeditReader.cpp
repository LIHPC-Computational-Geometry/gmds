/*----------------------------------------------------------------------------*/
/** \file    MeditReader.cpp
 *  \author  F. LEDOUX
 *  \date    09/11/2008
 */
/*----------------------------------------------------------------------------*/
#include <gmds/io/MeditReader.h>
#include <map>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
MeditReader::
MeditReader(gmds::IMeshIOService *AMeshService, const bool AUseMeditLabels)
        : IReader(AMeshService), m_with_labels(AUseMeditLabels)
{}
/*----------------------------------------------------------------------------*/
MeditReader::~MeditReader()
= default;
/*----------------------------------------------------------------------------*/
void MeditReader::setOptionWithLabel(const bool AWithLabels)
{
    m_with_labels=AWithLabels;
}
/*----------------------------------------------------------------------------*/
bool MeditReader::preCheckFormat() {
    bool find_dim = moveStreamOntoFirst("Dimension");

    if(!find_dim){
        return false;
    }
    *m_stream >> m_mesh_dimension;

    return true;
}
/*----------------------------------------------------------------------------*/
void MeditReader::readNodes() {
    bool find_vertices = moveStreamOntoFirst("Vertices");

    if (!find_vertices) {
        throw GMDSException("Error - no vertices field in a medit file");
    }
    TInt nb_nodes;
    (*m_stream) >> nb_nodes;

    std::vector<double> x, y, z;
    x.resize(nb_nodes);
    y.resize(nb_nodes);
    z.resize(nb_nodes);

    std::map<int, std::string> labels;
    int ref=0;
    if (m_mesh_dimension == 2) {
        for (auto i = 0; i < nb_nodes; i++) {
            (*m_stream) >> x[i] >> y[i] >> ref;
            z[i]=0;
        }
    }
    else if (m_mesh_dimension == 3) {
        for (auto i = 0; i < nb_nodes; i++) {
            (*m_stream) >> x[i] >> y[i] >> z[i] >> ref;
        }
    }

    m_mesh_service->createNodes(x, y, z);

}
/*----------------------------------------------------------------------------*/
void MeditReader::readEdges() {


    bool find_edges= moveStreamOntoFirst("Edges");
    if (!find_edges) {
        throw GMDSException("Error - no edges field in a medit file");
    }

    TInt nb_edges;
    *m_stream>>nb_edges;

    for(auto i = 0; i < nb_edges; i++) {
        TCellID x1, x2;
        int ref;
        *m_stream >> x1 >> x2 >> ref;
        m_mesh_service->createEdge(x1 - 1, x2 - 1);
    }
}
/*----------------------------------------------------------------------------*/
void MeditReader::readFaces() {

    bool find_triangles = moveStreamOntoFirst("Triangles");
    if (find_triangles) {
        readTriangles();
    }
    bool find_quads= moveStreamOntoFirst("Quadrilaterals");
    if (find_quads) {
        readQuadrilaterals();
    }

    if (!find_triangles && !find_quads) {
        throw GMDSException("Error - no triangles and quads field in a medit file");
    }
}
/*----------------------------------------------------------------------------*/
void MeditReader::readRegions()
{
    bool find_tet= moveStreamOntoFirst("Tetrahedra");
    if (find_tet) {
        readTetrahedra();
    }
    bool find_hex= moveStreamOntoFirst("Hexahedra");
    if (find_hex) {
        readHexahedra();
    }

    if (!find_tet && !find_hex) {
        throw GMDSException("Error - no tet and hex field in a medit file");
    }

}
/*----------------------------------------------------------------------------*/
void MeditReader::readQuadrilaterals() {
    TInt nb_quads;
    *m_stream>>nb_quads;
    for(auto i = 0; i < nb_quads; i++) {
        TCellID x1, x2, x3, x4;
        int ref;
        *m_stream>>x1>>x2>>x3>>x4>>ref;
        m_mesh_service->createQuad(x1-1,x2-1,x3-1,x4-1);
    }

}
/*----------------------------------------------------------------------------*/
void MeditReader::readTriangles() {
    TInt nb_triangles;
    *m_stream>>nb_triangles;
    for(auto i = 0; i < nb_triangles; i++) {
        TCellID x1, x2, x3;
        int ref;
        *m_stream>>x1>>x2>>x3>>ref;
        m_mesh_service->createTriangle(x1-1, x2-1, x3-1);
    }
}
/*----------------------------------------------------------------------------*/
void MeditReader::readHexahedra() {
    TInt nb_hex;
    *m_stream>>nb_hex;
    for(auto i = 0; i < nb_hex; i++)
    {
        TCellID x1, x2, x3, x4, x5, x6, x7, x8;
        int ref;

        *m_stream>>x1>>x2>>x3>>x4>>x5>>x6>>x7>>x8>>ref;
        m_mesh_service->createHex(x1-1,x2-1,x3-1,x4-1,x5-1,x6-1,x7-1,x8-1);
    }

}
/*----------------------------------------------------------------------------*/
void MeditReader::readTetrahedra() {
    TInt nb_tets;
    *m_stream>>nb_tets;
    for(auto i = 0; i < nb_tets; i++)
    {
        TCellID x1, x2, x3, x4;
        int ref;
        *m_stream>>x1>>x2>>x3>>x4>>ref;
        m_mesh_service->createTet(x1-1,x2-1,x3-1,x4-1);
    }

}
/*----------------------------------------------------------------------------*/
