/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryOperator2D.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
BoundaryOperator2D::BoundaryOperator2D(Mesh* AMesh)
        :m_mesh(AMesh)
{}
/*----------------------------------------------------------------------------*/
BoundaryOperator2D::~BoundaryOperator2D()
= default;
/*----------------------------------------------------------------------------*/
bool BoundaryOperator2D::isValid() const
{
    MeshModel model = m_mesh->getModel();
    if(model.has(R)) {
            return false;
    }
    if (!model.has(F) ||
        !model.has(N) ||
        !model.has(F2N) ||
        !model.has(N2F)  ){
        return false;
    }

    return true;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator2D::getBoundaryNodes(std::vector<TCellID>& ANodeIDs)
{
    //We only use N2F and F2N
    ANodeIDs.clear();
    for (auto n_id : m_mesh->nodes())  {
        Node n = m_mesh->get<Node>(n_id);
        std::vector<Face> faces_n = n.get<Face>();
        bool is_boundary=false;
        for(auto f_i =0; f_i<faces_n.size() && !is_boundary;f_i++){
            Node adj_n1, adj_n2;
            faces_n[f_i].getAdjacentNodes(n,adj_n1,adj_n2);
            std::vector<TCellID> adj_faces=m_mesh->getCommonFaces(n,adj_n1);
            if(adj_faces.size()==1){
                ANodeIDs.push_back(n.id());
                is_boundary=true;
            }
            else{
                adj_faces=m_mesh->getCommonFaces(n,adj_n2);
                if(adj_faces.size()==1){
                    ANodeIDs.push_back(n.id());
                    is_boundary= true;
                }
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator2D::
markCellOnGeometry(int AMarkEOnCurve,
                   int AMarkNOnCurve,
                   int AMarkNOnPnt,
                   int AMarkAN)
{
    markCellsOnCurves(AMarkEOnCurve, AMarkNOnCurve);

    markNodesOnPoint(AMarkEOnCurve, AMarkNOnCurve, AMarkNOnPnt);

    markAloneNodes(AMarkAN);
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator2D::markAloneNodes(int AMarkAlone)
{
    int cpt = 0;
    for (auto n_id:m_mesh->nodes())
    {
        Node n = m_mesh->get<Node>(n_id);
        if (n.get<Face>().empty()) {
            m_mesh->mark(n, AMarkAlone);
            cpt++;
        }

    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator2D::
markCellsOnCurves(int AMarkCE, //mark for edges on curves //OUT
                  int AMarkCN) //mark for nodes on curves //OUT
{
    if(m_mesh->getModel().has(E)) {
        for (auto e_id: m_mesh->edges()) {
            Edge e = m_mesh->get<Edge>(e_id);
            std::vector<Node> ne = e.get<Node>();
            std::vector<TCellID> adj_faces = m_mesh->getCommonFaces(ne[0], ne[1]);
            if (adj_faces.size() == 1) {
                //2D boundary edge
                m_mesh->mark(e, AMarkCE);
                std::vector<TCellID> e_nodes = e.getIDs<Node>();
                m_mesh->mark<Node>(e_nodes[0], AMarkCN);
                m_mesh->mark<Node>(e_nodes[1], AMarkCN);
            }
        }
    }
    else{
        //we only have nodes and no edges
        std::vector<TCellID> v;
        getBoundaryNodes(v);
        for(auto v_id:v){
            m_mesh->mark<Node>(v_id, AMarkCN);
        }
    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator2D::markNodesOnPoint(int AMarkCE,// edge on curve IN
                                        int AMarkCN,// node on curve IN
                                        int AMarkPN)// node on vertex OUT
{
    for (auto n_id:m_mesh->nodes()) {
        if (m_mesh->isMarked<Node>(n_id, AMarkCN)){
            Node n = m_mesh->get<Node>(n_id);
            //We get its incident boundary nodes
            std::vector<Face> adj_faces = n.get<Face>();
            std::vector<Node> adj_nodes;
            int nb_adj=0;
            for (auto i_f=0; i_f<adj_faces.size() && nb_adj<2; i_f++) {
                Node adj_n1, adj_n2;
                adj_faces[i_f].getAdjacentNodes(n, adj_n1, adj_n2);
                std::vector<TCellID> common_f = m_mesh->getCommonFaces(n, adj_n1);
                if (common_f.size() == 1) {
                    adj_nodes.push_back(adj_n1);
                    nb_adj++;
                }
                common_f = m_mesh->getCommonFaces(n, adj_n2);
                if (common_f.size() == 1) {
                    adj_nodes.push_back(adj_n2);
                    nb_adj++;
                }
            }
            //So we have boundary node and its two boundary neighbours
            //We check if we have a brutal normal change

            Node n0 = adj_nodes[0];
            Node n1 = adj_nodes[1];

            math::Vector3d v0=n0.point()- n.point();
            math::Vector3d v1=n1.point()- n.point();
            v0.normalize();
            v1.normalize();
            double dotProduct = v0.dot(v1);
            // if we have a brutal normal change
            if (dotProduct > (-sqrt(2.0) / 2.0)){
                m_mesh->mark(n, AMarkPN);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator2D::
colorEdges(int AMarkEOnCurv, int AMarkNOnPnt, Variable<int>* AColor)
{
    Variable<int>* var_color = AColor;
    if(var_color==nullptr) {
        try {
            var_color = m_mesh->newVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");
        } catch (GMDSException &e) {
            var_color = m_mesh->getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");
        }
    }
    int color = 0; //Default value is 0
    TInt markDone = m_mesh->newMark<Edge>();

    //as we do not have N2E, we build it temporarily
    std::map<TCellID , std::vector<TCellID> > local_N2E;

    for (auto e_id:m_mesh->edges())
    {
        Edge ei = m_mesh->get<Edge>(e_id);
        std::vector<TCellID> ei_nodes = ei.getIDs<Node>();
        for(auto nj:ei_nodes)
            local_N2E[nj].push_back(ei.id());
    }

    std::vector<Edge> done_edges;
    for (auto e_id:m_mesh->edges())
    {
        // We only go through edges classified on curves and that have not been
        // yet handled
        if ( m_mesh->isMarked<Edge>(e_id, AMarkEOnCurv) &&
             !m_mesh->isMarked<Edge>(e_id, markDone)) {
            Edge e = m_mesh->get<Edge>(e_id);
            //new curve
            color++; // so new color
            m_mesh->mark(e, markDone);
            (*var_color)[e_id] = color;

            //propagation to curve edges sharing a point with e
            std::vector<Edge> next;
            next.push_back(e);

            while (!next.empty()){
                Edge current = next.back();
                next.pop_back();
                //We get the ajacent edges that are on a curve but not yet done
                std::vector<Node> current_nodes = current.get<Node>();

                for (const auto& ni: current_nodes) {

                    if (!m_mesh->isMarked(ni, AMarkNOnPnt)){
                        //If it is not a node classified on a point, we can found
                        // a next edge on this curve

                        std::vector<TCellID> n_edges = local_N2E[ni.id()];
                        for (auto ej:n_edges){
                            if (m_mesh->isMarked<Edge>(ej, AMarkEOnCurv) &&
                                !m_mesh->isMarked<Edge>(ej, markDone)){
                                m_mesh->mark<Edge>(ej, markDone);
                                (*var_color)[ej] = color;
                                next.push_back(m_mesh->get<Edge>(ej));

                            }
                        }
                    }
                }
            }
        }
    }
    m_mesh->unmarkAll<Edge>(markDone);
    m_mesh->freeMark<Edge>(markDone);
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator2D::
colorNodes(int AMarkNOnPnt, Variable<int>* AColor)
{
    Variable<int>* v_color=AColor;
    if(v_color==nullptr){
        try {
            v_color = m_mesh->newVariable<int, GMDS_NODE>("BND_VERTEX_COLOR");
        } catch (GMDSException &e) {
            v_color = m_mesh->getVariable<int, GMDS_NODE>("BND_VERTEX_COLOR");
        }
    }

    int color = 0; //Default value is 0
    for (auto n_id : m_mesh->nodes() ){
        if (m_mesh->isMarked<Node>(n_id, AMarkNOnPnt)) {
            //new point
            color++; // so new color
            (*v_color)[n_id] = color;
        }
    }
}