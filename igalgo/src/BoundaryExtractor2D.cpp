/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryExtractor2D.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/ig/MeshDoctor.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
#include <gmds/math/Orientation.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
BoundaryExtractor2D::
BoundaryExtractor2D(gmds::Mesh *AFromMesh, gmds::Mesh *AToMesh)
: m_from_mesh(AFromMesh), m_to_mesh(AToMesh)
{
    m_with_color = false;
    m_with_mapping=false;

    m_color_node_on_pnt  = nullptr ;
    m_color_node_on_curv = nullptr;
    m_color_edge_on_curv = nullptr;

    m_node_map = nullptr;
    m_edge_map = nullptr;

    m_node_map_inv = nullptr;
    m_edge_map_inv = nullptr;
}
/*----------------------------------------------------------------------------*/
BoundaryExtractor2D::~BoundaryExtractor2D()
= default;
/*----------------------------------------------------------------------------*/
bool BoundaryExtractor2D::isValid() const
{
    //check the input first
    bool valid_input=true;
    MeshModel model = m_from_mesh->getModel();
    if (!model.has(F) || !model.has(N2F) || !model.has(F2N))
        valid_input= false;

    if(!valid_input)
        return false;

    //and now, the output
    model = m_to_mesh->getModel();
    if (model.has(R))
        return false;
    if (!model.has(E) || !model.has(E2N) || !model.has(N2E))
        return false;


    return true;
}

/*----------------------------------------------------------------------------*/
void BoundaryExtractor2D::setColorOption(Variable<int> *ANodeColor,
                                         Variable<int> *AEdgeColor,
                                         Variable<int>* ANodeColorCurv)
{

    if(!m_to_mesh->hasVariable(GMDS_NODE, ANodeColor->getName())){
        std::string mess = "The node variable "+ANodeColor->getName();
        mess +=" is unknown by the target surface mesh";
        throw GMDSException(mess);
    }
    if(!m_to_mesh->hasVariable(GMDS_NODE, ANodeColorCurv->getName())){
        std::string mess = "The node variable "+ANodeColorCurv->getName();
        mess +=" is unknown by the target surface mesh";
        throw GMDSException(mess);
    }
    if(!m_to_mesh->hasVariable(GMDS_EDGE, AEdgeColor->getName())){
        std::string mess = "The edge variable "+AEdgeColor->getName();
        mess +=" is unknown by the target surface mesh";
        throw GMDSException(mess);
    }

    m_color_node_on_pnt = ANodeColor;
    m_color_edge_on_curv = AEdgeColor;
    m_color_node_on_curv = ANodeColorCurv;
    m_with_color = true;
}
/*----------------------------------------------------------------------------*/
void BoundaryExtractor2D::setMappings(std::map<TCellID, TCellID> *ANodeMap,
                                      std::map<TCellID, TCellID> *AEdgeMap,
                                      std::map<TCellID, TCellID> *ANodeMapInv,
                                      std::map<TCellID, TCellID> *AEdgeMapInv)
{
    m_node_map      = ANodeMap;
    m_edge_map      = AEdgeMap;
    m_node_map_inv  = ANodeMapInv;
    m_edge_map_inv  = AEdgeMapInv;
    m_with_mapping  = true;
}

/*----------------------------------------------------------------------------*/
void BoundaryExtractor2D::execute()
{
    //the call to the boundary operator is generic in 2D and 3D for many
    //treatments. So we use her more marks that we could expect first
    TInt mark_node_on_curv   = m_from_mesh->newMark<Node>();
	 TInt mark_node_on_pnt    = m_from_mesh->newMark<Node>();
	 TInt mark_node_isolated  = m_from_mesh->newMark<Node>();

	 TInt mark_edge_on_curv   = m_from_mesh->newMark<Edge>();


    BoundaryOperator2D boundaryOp(m_from_mesh);

    if (!boundaryOp.isValid()) {
        std::cout << "Invalid model for boundary operations" << std::endl;
        throw GMDSException("Invalid model for boundary operations");
    }

    //==================================================================
    // Mark boundary cells
    boundaryOp.markCellOnGeometry(mark_edge_on_curv,
                                  mark_node_on_curv,
                                  mark_node_on_pnt,
                                  mark_node_isolated);

    //==================================================================
    Variable<int>* from_edge_color =
            m_from_mesh->newVariable<int, GMDS_EDGE>("BoundaryExtractor3D::exec");
    Variable<int>* from_node_color =
            m_from_mesh->newVariable<int, GMDS_NODE>("BoundaryExtractor3D::exec");

    if (m_with_color) {
        // color  edges on curves, nodes on points
        boundaryOp.colorEdges(mark_edge_on_curv, mark_node_on_pnt,
                              from_edge_color);
        boundaryOp.colorNodes(mark_node_on_pnt, from_node_color);
    }

    //We have now to build the skin mesh and to add fill colors and mappings.
    std::map<TCellID , TCellID > node_map;

    //first we build nodes, node mapping and color of nodes on pnts
    for(auto n_id: m_from_mesh->nodes()){
        if(m_from_mesh->isMarked<Node>(n_id,mark_node_on_curv)||
           m_from_mesh->isMarked<Node>(n_id,mark_node_on_pnt)){
            Node from_n = m_from_mesh->get<Node>(n_id);
            Node to_n   = m_to_mesh->newNode(from_n.point());
            if(m_with_mapping){
                (*m_node_map)[n_id] = to_n.id();
                (*m_node_map_inv)[to_n.id()] = n_id;
            }
            else{
                node_map[n_id] = to_n.id();
            }
            if(m_with_color){
                (*m_color_node_on_pnt)[to_n.id()] = (*from_node_color)[n_id];
            }
        }
    }
    //then we build edges, edge mapping and color of edges and nodes on curve
    for(auto e_id: m_from_mesh->edges()){
        if(m_from_mesh->isMarked<Edge>(e_id,mark_edge_on_curv)){
            Edge from_e = m_from_mesh->get<Edge>(e_id);
            std::vector<TCellID> from_ns = from_e.getIDs<Node>();
            Edge to_e;
            if(m_with_mapping){
                to_e = m_to_mesh->newEdge((*m_node_map)[from_ns[0]],
                                          (*m_node_map)[from_ns[1]]);

                (*m_edge_map)[e_id] = to_e.id();
                (*m_edge_map_inv)[to_e.id()] = e_id;
            }
            else{
                to_e = m_to_mesh->newEdge(node_map[from_ns[0]],
                                          node_map[from_ns[1]]);
            }

            if(m_with_color){
                (*m_color_edge_on_curv)[to_e.id()] = (*from_edge_color)[e_id];
                if(m_with_mapping){
                    for(auto i=0;i<2;i++) {
                        if (m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_curv) &&
                            !m_from_mesh->isMarked<Node>(from_ns[i], mark_node_on_pnt)) {
                            (*m_color_node_on_curv)[(*m_node_map)[from_ns[i]]] =
                                    (*from_edge_color)[e_id];
                        }
                    }
                }
                else{
                    for(auto i=0;i<2;i++) {
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

    m_from_mesh->unmarkAll<Node>(mark_node_on_curv);
    m_from_mesh->unmarkAll<Node>(mark_node_on_pnt);
    m_from_mesh->unmarkAll<Node>(mark_node_isolated);
    m_from_mesh->unmarkAll<Edge>(mark_edge_on_curv);

    m_from_mesh->freeMark<Node>(mark_node_on_curv);
    m_from_mesh->freeMark<Node>(mark_node_on_pnt);
    m_from_mesh->freeMark<Node>(mark_node_isolated);
    m_from_mesh->freeMark<Edge>(mark_edge_on_curv);

    m_from_mesh->deleteVariable(GMDS_NODE, from_node_color);
    m_from_mesh->deleteVariable(GMDS_EDGE, from_edge_color);

}
