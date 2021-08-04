/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryExtractor3D.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/math/Orientation.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;

/*----------------------------------------------------------------------------*/
BoundaryExtractor3D::
BoundaryExtractor3D(gmds::Mesh *AFromMesh, gmds::Mesh *AToMesh)
        : m_from_mesh(AFromMesh), m_to_mesh(AToMesh) {
    m_with_color = false;
    m_with_mapping = false;

    m_color_node_on_pnt = NULL;
    m_color_node_on_curv = NULL;
    m_color_node_on_surf = NULL;
    m_color_edge_on_curv = NULL;
    m_color_face_on_surf = NULL;

    m_node_map = NULL;
    m_edge_map = NULL;
    m_face_map = NULL;

    m_node_map_inv = NULL;
    m_edge_map_inv = NULL;
    m_face_map_inv = NULL;
}

/*----------------------------------------------------------------------------*/
BoundaryExtractor3D::~BoundaryExtractor3D() {
}

/*----------------------------------------------------------------------------*/
bool BoundaryExtractor3D::isValid() const {
    //check the input first
    bool valid_input = true;
    MeshModel model = m_from_mesh->getModel();
    if (!model.has(R) || !model.has(F2R) || !model.has(F2E) || !model.has(E2F))
        valid_input = false;

    if (!valid_input)
        return false;

    //and now, the output
    model = m_to_mesh->getModel();
    if (model.has(R))
        return false;
    if (!model.has(F) || !model.has(E) || !model.has(F2E) || !model.has(E2F))
        return false;


    return true;
}

/*----------------------------------------------------------------------------*/
void BoundaryExtractor3D::setColorOption(Variable<int> *ANodeColor,
                                         Variable<int> *AEdgeColor,
                                         Variable<int> *AFaceColor,
                                         Variable<int> *ANodeColorSurf,
                                         Variable<int> *ANodeColorCurv) {

    if (!m_to_mesh->hasVariable(GMDS_NODE, ANodeColor->getName())) {
        std::string mess = "The node variable " + ANodeColor->getName();
        mess += " is unknown by the target surface mesh";
        throw GMDSException(mess);
    }
    if (!m_to_mesh->hasVariable(GMDS_NODE, ANodeColorCurv->getName())) {
        std::string mess = "The node variable " + ANodeColorCurv->getName();
        mess += " is unknown by the target surface mesh";
        throw GMDSException(mess);
    }
    if (!m_to_mesh->hasVariable(GMDS_NODE, ANodeColorSurf->getName())) {
        std::string mess = "The node variable " + ANodeColorSurf->getName();
        mess += " is unknown by the target surface mesh";
        throw GMDSException(mess);
    }
    if (!m_to_mesh->hasVariable(GMDS_EDGE, AEdgeColor->getName())) {
        std::string mess = "The edge variable " + AEdgeColor->getName();
        mess += " is unknown by the target surface mesh";
        throw GMDSException(mess);
    }
    if (!m_to_mesh->hasVariable(GMDS_FACE, AFaceColor->getName())) {
        std::string mess = "The face variable " + AFaceColor->getName();
        mess += " is unknown by the target surface mesh";
        throw GMDSException(mess);
    }

    m_color_node_on_pnt = ANodeColor;
    m_color_edge_on_curv = AEdgeColor;
    m_color_face_on_surf = AFaceColor;
    m_color_node_on_surf = ANodeColorSurf;
    m_color_node_on_curv = ANodeColorCurv;
    m_with_color = true;
}

/*----------------------------------------------------------------------------*/
void BoundaryExtractor3D::setMappings(std::map<TCellID, TCellID> *ANodeMap,
                                      std::map<TCellID, TCellID> *AEdgeMap,
                                      std::map<TCellID, TCellID> *AFaceMap,
                                      std::map<TCellID, TCellID> *ANodeMapInv,
                                      std::map<TCellID, TCellID> *AEdgeMapInv,
                                      std::map<TCellID, TCellID> *AFaceMapInv) {
    m_node_map = ANodeMap;
    m_edge_map = AEdgeMap;
    m_face_map = AFaceMap;
    m_node_map_inv = ANodeMapInv;
    m_edge_map_inv = AEdgeMapInv;
    m_face_map_inv = AFaceMapInv;
    m_with_mapping = true;
}

/*----------------------------------------------------------------------------*/
void BoundaryExtractor3D::execute(const double AAngle) {
    int mark_node_on_surf = m_from_mesh->newMark<Node>();
    int mark_node_on_curv = m_from_mesh->newMark<Node>();
    int mark_node_on_pnt = m_from_mesh->newMark<Node>();
    int mark_node_isolated = m_from_mesh->newMark<Node>();

    int mark_node_frame = m_from_mesh->newMark<Node>();
    int mark_node_hard = m_from_mesh->newMark<Node>();

    int mark_edge_on_surf = m_from_mesh->newMark<Edge>();
    int mark_edge_on_curv = m_from_mesh->newMark<Edge>();

    int mark_face_on_surf = m_from_mesh->newMark<Face>();

    BoundaryOperator boundaryOp(m_from_mesh, AAngle);

    if (!boundaryOp.isValid()) {
        std::cout << "Invalid model for boundary operations" << std::endl;
        throw GMDSException("Invalid model for boundary operations");
    }

    //==================================================================
    // Mark boundary cells
    boundaryOp.markCellOnGeometry(mark_face_on_surf,
                                  mark_edge_on_surf,
                                  mark_node_on_surf,
                                  mark_edge_on_curv,
                                  mark_node_on_curv,
                                  mark_node_on_pnt,
                                  mark_node_isolated);

    //==================================================================
    Variable<int> *from_face_color =
            m_from_mesh->newVariable<int, GMDS_FACE>("BoundaryExtractor3D::exec");
    Variable<int> *from_edge_color =
            m_from_mesh->newVariable<int, GMDS_EDGE>("BoundaryExtractor3D::exec");
    Variable<int> *from_node_color =
            m_from_mesh->newVariable<int, GMDS_NODE>("BoundaryExtractor3D::exec");

    if (m_with_color) {
        // color faces on surfaces, edges on curves, nodes on points
        boundaryOp.colorFaces(mark_face_on_surf, mark_edge_on_curv,
                              from_face_color);
        boundaryOp.colorEdges(mark_edge_on_curv, mark_node_on_pnt,
                              from_edge_color);
        boundaryOp.colorNodes(mark_node_on_pnt, from_node_color);
    }

    //We have now to build the skin mesh and to add fill colors and mappings.
    std::map<TCellID, TCellID> node_map;


    //first we build nodes, node mapping and color of nodes on pnts
    for (auto n_id: m_from_mesh->nodes()) {
        if (m_from_mesh->isMarked<Node>(n_id, mark_node_on_surf)) {
            Node from_n = m_from_mesh->get<Node>(n_id);
            Node to_n = m_to_mesh->newNode(from_n.getPoint());
            if (m_with_mapping) {
                (*m_node_map)[n_id] = to_n.id();
                (*m_node_map_inv)[to_n.id()] = n_id;
            } else {
                node_map[n_id] = to_n.id();
            }
            if (m_with_color) {
                (*m_color_node_on_pnt)[to_n.id()] = (*from_node_color)[n_id];
            }
        }
    }
    //then we build edges, edge mapping and color of edges and nodes on curve
    for (auto e_id: m_from_mesh->edges()) {
        if (m_from_mesh->isMarked<Edge>(e_id, mark_edge_on_curv)) {
            Edge from_e = m_from_mesh->get<Edge>(e_id);
            std::vector<TCellID> from_ns = from_e.getIDs<Node>();
            Edge to_e;
            if (m_with_mapping) {
                to_e = m_to_mesh->newEdge((*m_node_map)[from_ns[0]],
                                          (*m_node_map)[from_ns[1]]);

                (*m_edge_map)[e_id] = to_e.id();
                (*m_edge_map_inv)[to_e.id()] = e_id;
            } else {
                to_e = m_to_mesh->newEdge(node_map[from_ns[0]],
                                          node_map[from_ns[1]]);
            }

            if (m_with_color) {
                (*m_color_edge_on_curv)[to_e.id()] = (*from_edge_color)[e_id];
                if (m_with_mapping) {
                    for (auto i = 0; i < 2; i++) {
                        if (m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_curv) &&
                            !m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_pnt)) {
                            (*m_color_node_on_curv)[(*m_node_map)[from_ns[i]]] =
                                    (*from_edge_color)[e_id];
                        }
                    }
                } else {
                    for (auto i = 0; i < 2; i++) {
                        if (m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_curv) &&
                            !m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_pnt)) {
                            (*m_color_node_on_curv)[node_map[from_ns[i]]] =
                                    (*from_edge_color)[e_id];
                        }
                    }
                }
            }
        }
    }

    //and now the faces
    for (auto f_id: m_from_mesh->faces()) {
        if (m_from_mesh->isMarked<Face>(f_id, mark_face_on_surf)) {
            Face from_f = m_from_mesh->get<Face>(f_id);
            std::vector<TCellID> from_ns = from_f.getIDs<Node>();
            std::vector<TCellID> to_ns;
            to_ns.reserve(from_ns.size());

            //we invert the created face if necessary to get outward normals
            //from_f should be adjacent to only one tet
            Region adj_tet = from_f.get<Region>()[0];
            std::vector<TCellID> adj_nodes = adj_tet.getIDs<Node>();
            //we get the opposite node in this tet
            TCellID opposite_id = NullID;
            for (auto adj_n:adj_nodes) {
                bool found = false;
                for (auto face_nid:from_ns) {
                    if (face_nid == adj_n)
                        found = true;
                }
                if (!found)
                    opposite_id = adj_n;
            }
            math::Point opp_point = m_from_mesh->get<Node>(opposite_id).getPoint();
            math::Point p0 = m_from_mesh->get<Node>(from_ns[0]).getPoint();
            math::Point p1 = m_from_mesh->get<Node>(from_ns[1]).getPoint();
            math::Point p2 = m_from_mesh->get<Node>(from_ns[2]).getPoint();
            if (math::Orientation::orient3d(opp_point, p0, p1, p2) == math::Orientation::NEGATIVE) {
                std::vector<TCellID> temp_ns = from_ns;
                from_ns[1] = temp_ns[temp_ns.size() - 1];
                from_ns[from_ns.size() - 1] = temp_ns[1];

            }
            Face to_f;
            if (m_with_mapping) {
                for (auto i:from_ns)
                    to_ns.push_back((*m_node_map)[i]);

                to_f = m_to_mesh->newFace(to_ns);

                (*m_face_map)[f_id] = to_f.id();
                (*m_face_map_inv)[to_f.id()] = f_id;
            } else {
                for (auto i:from_ns)
                    to_ns.push_back(node_map[i]);

                to_f = m_to_mesh->newFace(to_ns);

                //NO NODE->FACE CONNECTIVITY

            }
            if (m_with_color) {
                (*m_color_face_on_surf)[to_f.id()] = (*from_face_color)[f_id];
                if (m_with_mapping) {
                    for (auto i:from_ns) {
                        if (m_from_mesh->isMarked<Node>(i, mark_node_on_surf) &&
                            !m_from_mesh->isMarked<Node>(i, mark_node_on_curv) &&
                            !m_from_mesh->isMarked<Node>(i, mark_node_on_pnt)) {
                            (*m_color_node_on_surf)[(*m_node_map)[i]] =
                                    (*from_face_color)[f_id];
                        }
                    }
                } else {
                    for (auto i:from_ns) {
                        if (m_from_mesh->isMarked<Node>(i, mark_node_on_surf) &&
                            !m_from_mesh->isMarked<Node>(i, mark_node_on_curv) &&
                            !m_from_mesh->isMarked<Node>(i, mark_node_on_pnt)) {
                            (*m_color_node_on_surf)[node_map[i]] =
                                    (*from_face_color)[f_id];
                        }
                    }
                }
            }
        }
    }


    m_from_mesh->unmarkAll<Node>(mark_node_on_surf);
    m_from_mesh->unmarkAll<Node>(mark_node_on_curv);
    m_from_mesh->unmarkAll<Node>(mark_node_on_pnt);
    m_from_mesh->unmarkAll<Node>(mark_node_isolated);
    m_from_mesh->unmarkAll<Node>(mark_node_frame);
    m_from_mesh->unmarkAll<Node>(mark_node_hard);
    m_from_mesh->unmarkAll<Edge>(mark_edge_on_surf);
    m_from_mesh->unmarkAll<Edge>(mark_edge_on_curv);
    m_from_mesh->unmarkAll<Face>(mark_face_on_surf);

    m_from_mesh->freeMark<Node>(mark_node_on_surf);
    m_from_mesh->freeMark<Node>(mark_node_on_curv);
    m_from_mesh->freeMark<Node>(mark_node_on_pnt);
    m_from_mesh->freeMark<Node>(mark_node_isolated);
    m_from_mesh->freeMark<Node>(mark_node_frame);
    m_from_mesh->freeMark<Node>(mark_node_hard);
    m_from_mesh->freeMark<Edge>(mark_edge_on_surf);
    m_from_mesh->freeMark<Edge>(mark_edge_on_curv);
    m_from_mesh->freeMark<Face>(mark_face_on_surf);

    m_from_mesh->deleteVariable(GMDS_NODE, from_node_color);
    m_from_mesh->deleteVariable(GMDS_EDGE, from_edge_color);
    m_from_mesh->deleteVariable(GMDS_FACE, from_face_color);

}
